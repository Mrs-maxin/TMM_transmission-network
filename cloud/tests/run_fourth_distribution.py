"""IEEE118 OP2 complete propagation and frozen-GMM fourth-moment comparison.

Retains every standardized quadrature point and weight in small compressed blocks.
Focused ECDF files retain every sorted value and cumulative probability.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import argparse,json,sys,time,copy
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
import pandas as pd
from tmm_cloud import mixture as mx,montecarlo as m,network as nw,matpower_case as mc,state_model as sm,propagation as pg,run_record as rr,violation_events as ve
from tmm_cloud.gmm_fourth_moment_tmm import fit_gmm_fourth_moment_tmm
from tmm_cloud.fourth_comparison_metrics import PROBS,TAGS,distribution,inverse

O=rr.CLOUD_ROOT/'results/ieee118_fourth_moment_op2/20261030T120000'
D=rr.CLOUD_ROOT/'results/phase1_rts118/20261030T120000/stage_d_montecarlo/OP2'
H=rr.CLOUD_ROOT/'results/phase1b_tmm_op2/20261022T125454/stage_h_fit'
B=rr.CLOUD_ROOT/'results/phase1_rts118/20261021T201240/stage_b_network'
N=1<<20
ENTROPY=2026102201
KEYS=['Vm','Va','Pf','Qf','Pt','Qt','loading','lam','Pg']
TOL={'Vm':1e-6,'Va':1e-5,'Pf':1e-3,'Qf':1e-3,'Pt':1e-3,'Qt':1e-3,'loading':1e-6,'lam':1e-6,'Pg':1e-6}
UNITS={'Vm':'p.u.','Va':'deg','Pf':'MW','Qf':'Mvar','Pt':'MW','Qt':'Mvar','loading':'fraction of rateA','lam':'MW','Pg':'MW'}

def covariance(model):
    mu=model.weights@model.means
    delta=model.means-mu
    factors=np.where(np.isinf(model.df),1.,model.df/(model.df-2))
    factors[np.isinf(model.df)]=1.
    cov=np.einsum('k,kij,k->ij',model.weights,model.scatter,factors)+np.einsum('k,ki,kj->ij',model.weights,delta,delta)
    return mu,cov


def setup():
    net=nw.from_matpower(mc.PGLIB118);op=m.load_operating_point(B/'operating_point_OP2.npz');outmap=sm.OutputMap(net,op)
    z=np.load(H/'coordinates_OP2.npz');coords=mx.Coordinates(z['center'],z['basis'],z['scale'],z['eigenvalues'],1e-10,0.)
    assert coords.d==157 and len(coords.center)==235
    g=mx.from_npz(np.load(H/'model_GMM_K4_OP2.npz'));old=mx.from_npz(np.load(H/'model_TMM_K8_OP2.npz'))
    paths=sorted((O/'modeling_blocks').glob('modeling_chunk*.npz'))
    x,dm=sm.scenario_state(paths,net);assert x.shape==(5256,235)
    frozen=copy.deepcopy(g);new=fit_gmm_fourth_moment_tmm(coords.to_z(x),g)
    for a in ('weights','means','scatter','df'):assert np.array_equal(getattr(g,a),getattr(frozen,a))
    mg,cg=covariance(g);mt,ct=covariance(new)
    assert np.array_equal(mg,mt) and np.allclose(cg,ct,rtol=1e-12,atol=1e-14)
    rr.write_json(O/'calibration.json',new.audit)
    pd.DataFrame(new.audit['components']).to_csv(O/'component_moments.csv',index=False)
    np.savez_compressed(O/'model_fourth_TMM_K4_OP2.npz',**mx.to_npz_dict(new))
    nr=np.setdiff1d(np.arange(net.n),[net.ref])
    branches=[f'branch_index0_{j}_index1_{j+1}_{net.bus_id[f]}-{net.bus_id[t]}' for j,(f,t) in enumerate(zip(net.f,net.t))]
    names={'Vm':[f'bus{net.bus_id[i]}' for i in nr],'Va':[f'bus{net.bus_id[i]}' for i in nr],**{k:branches for k in ['Pf','Qf','Pt','Qt','loading']},'lam':['system'],'Pg':[f'unit_bus{net.bus_id[b]}' for b in net.gen_bus[net.Pmax>0]]}
    focus={('Vm',j) for j,i in enumerate(nr) if net.bus_id[i] in [9,17,83,37]}
    focus|={(key,j) for j,(f,t) in enumerate(zip(net.f,net.t)) if (net.bus_id[f],net.bus_id[t]) in [(100,103),(77,82)] for key in ['Pf','Qf','Pt','Qt','loading']}
    for key,j in sorted(focus):print('focus',key,j,names[key][j],flush=True)
    vset=np.full(net.n,np.nan);vset[net.gen_bus]=op.Vset
    rr.write_json(O/'control_checks.json',dict(source_GMM_unchanged=True,coordinate_refit=False,dimension=coords.d,source_state_dimension=len(coords.center),mean_max_abs_change=float(np.max(abs(mg-mt))),covariance_max_abs_change=float(np.max(abs(cg-ct))),component_covariance_preserved=True,posterior_update=False,source_commit='32bc9ae18e6ce8712e5f8a9c66d4833ab08a1d54',physical_definitions='source OutputMap and Stage D unmodified; max both-end apparent power/rateA; bus83 event Vmax+1e-6',modeling_rows=len(x)))
    return net,op,outmap,coords,{'GMM_K4':g,'fourth_TMM_K4':new,'MLE_TMM_K8':old},names,focus,nr,vset


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--method',choices=['GMM_K4','fourth_TMM_K4','MLE_TMM_K8'],required=True);ap.add_argument('--skip-identity',action='store_true');args=ap.parse_args()
    started=time.perf_counter();method=args.method
    if not args.skip_identity:assert json.loads((O/'identity.json').read_text())['identical']
    net,op,outmap,coords,models,names,focus,nr,vset=setup()
    model=models[method];seed=ENTROPY+1 if method=='MLE_TMM_K8' else ENTROPY
    dest=O/method;dest.mkdir(exist_ok=True);points=dest/'points';points.mkdir(exist_ok=True);cdf=dest/'cdf';cdf.mkdir(exist_ok=True)
    reference=outmap.from_stage_d(m.load_set(sorted(D.glob('check_heavy_chunk*.npz'))));assert len(reference['Vm'])==100000
    print(method,'generating',N,'Sobol points',flush=True)
    t=time.perf_counter();Z,w,alloc=pg.model_points(model,N,seed);w=w/w.sum();sampling=time.perf_counter()-t
    if method=='fourth_TMM_K4' and np.isinf(model.df).all():
        zg,wg,_=pg.model_points(models['GMM_K4'],N,seed);assert np.array_equal(Z,zg) and np.array_equal(w,wg);del zg,wg
    mm={};shape=outmap.evaluate(coords.to_x(Z[:1]));scratch=dest/'temporary_arrays';scratch.mkdir(exist_ok=True)
    for k in KEYS:mm[k]=np.lib.format.open_memmap(scratch/(k+'.npy'),mode='w+',dtype=np.float64,shape=(N,shape[k].shape[1]))
    thr=ve.thresholds(net,op);evsum=np.zeros(len(ve.event_names((1.,1.05,1.10))));event_names=ve.event_names((1.,1.05,1.10))
    t=time.perf_counter()
    for s in range(0,N,4096):
        e=min(s+4096,N)
        np.savez_compressed(points/f'points_{s//4096:04d}.npz',Z=Z[s:e],weight=w[s:e],index_start=np.array(s))
        q=outmap.evaluate(coords.to_x(Z[s:e]))
        for k in KEYS:mm[k][s:e]=q[k]
        ev,_=ve.indicators(q,thr,(1.,1.05,1.10));evsum+=w[s:e]@ev
        if s%(4096*16)==0:print(method,'propagation',e,'/',N,flush=True)
    for a in mm.values():a.flush()
    propagation=time.perf_counter()-t;del Z
    erows=[];mrows=[];fixed=[];events=[];allq=[];plotrows=[];histrows=[]
    ref_ev,_=ve.indicators(reference,thr,(1.,1.05,1.10))
    for j,name in enumerate(event_names):events.append(dict(method=method,event=name,probability=float(evsum[j]),MC_probability=float(ref_ev[:,j].mean()),MC_count=int(ref_ev[:,j].sum()),error=float(evsum[j]-ref_ev[:,j].mean())))
    t=time.perf_counter()
    for k in KEYS:
        Y=mm.pop(k);R=reference[k]
        for j in range(Y.shape[1]):
            y=np.asarray(Y[:,j]);r=R[:,j];order=np.argsort(y,kind='stable');ys=y[order];cw=np.cumsum(w[order]);cw[-1]=1.
            mean=float(w@y);var=float(w@((y-mean)**2));rm=float(r.mean());rv=float(r.var());sd=float(np.sqrt(rv))
            meta=dict(method=method,quantity=k,object=names[k][j],column=j,unit=UNITS[k],focus=(k,j) in focus,boundary=(k=='Vm' and names[k][j]=='bus37'))
            row={**meta,'mean':mean,'variance':var,'MC_mean':rm,'MC_variance':rv,'mean_error':mean-rm,'variance_error':var-rv,'variance_relative_error':var/rv-1 if rv>0 else np.nan}
            mrows.append(row)
            qs=inverse(ys,cw,PROBS);rq=np.quantile(r,PROBS,method='inverted_cdf')
            allq.append({**meta,**{tag:float(q) for tag,q in zip(TAGS,qs)},**{tag+'_MC':float(q) for tag,q in zip(TAGS,rq)}})
            if sd<TOL[k] or (k=='Vm' and net.bus_id[nr[j]] in net.bus_id[net.gen_bus]):
                fixed.append({**meta,'MC_std':sd,'near_fixed':sd<TOL[k],'MC_share_at_voltage_setpoint':float((abs(r-vset[nr[j]])<1e-6).mean()) if k=='Vm' and np.isfinite(vset[nr[j]]) else np.nan})
            if (k,j) in focus:
                erows.append({**meta,**distribution(ys,cw,r)})
                np.savez_compressed(cdf/f'{k}_{j:03d}.npz',sorted_values=ys,cumulative_probability=cw)
                sr=np.sort(r.astype(float));np.savez_compressed(cdf/f'MC_{k}_{j:03d}.npz',sorted_values=sr,cumulative_probability=np.arange(1,len(sr)+1)/len(sr))
                grid=np.unique(np.r_[np.linspace(0,1,1001),np.geomspace(1/N,.05,500),1-np.geomspace(1/N,.05,500)])
                xx=inverse(ys,cw,grid);rrq=sr[np.minimum(np.maximum(np.ceil(grid*len(sr)).astype(int)-1,0),len(sr)-1)]
                plotrows.extend(dict(method=method,quantity=k,object=names[k][j],probability=float(p),value=float(a),MC_value=float(b)) for p,a,b in zip(grid,xx,rrq))
                low=min(float(y.min()),float(r.min()));high=max(float(y.max()),float(r.max()))
                if high==low:low-=TOL[k];high+=TOL[k]
                bins=np.linspace(low,high,201);hm,_=np.histogram(y,bins,weights=w);hr,_=np.histogram(r,bins)
                histrows.extend(dict(method=method,quantity=k,object=names[k][j],bin_left=float(bins[i]),bin_right=float(bins[i+1]),density=float(hm[i]/(bins[i+1]-bins[i])),MC_density=float(hr[i]/len(r)/(bins[i+1]-bins[i]))) for i in range(len(hm)))
        del Y
        (scratch/(k+'.npy')).unlink()
        print(method,'completed',k,'objects',R.shape[1],flush=True)
        pd.DataFrame(erows).to_csv(dest/'focus_distribution_errors.csv',index=False)
        pd.DataFrame(mrows).to_csv(dest/'all_network_mean_variance.csv',index=False)
        pd.DataFrame(allq).to_csv(dest/'all_network_quantiles.csv',index=False)
    scratch.rmdir()
    pd.DataFrame(fixed).to_csv(dest/'fixed_and_concentrated.csv',index=False)
    pd.DataFrame(events).to_csv(dest/'event_probabilities.csv',index=False)
    pd.DataFrame(plotrows).to_csv(dest/'cdf_plot_data.csv',index=False)
    pd.DataFrame(histrows).to_csv(dest/'pdf_histogram_plot_data.csv',index=False)
    summary=dict(method=method,K=model.K,points=N,component_allocation=alloc.tolist(),entropy=seed,probability_weights='pi_k/N_k, normalized by floating sum',sobol='original propagation.model_points: scrambled scipy Sobol, per-component PCG64 SeedSequence(entropy).spawn(K); dimensions d for Gaussian, d+1 for finite t; original largest-remainder allocation',precision='float64 points and electrical quantities; original complex64/float32 AC reference',quantiles='generalized inverse; MC actual ECDF at tied values',tail_w1='exact union of weighted-model and MC probability breaks, divided by 0.049',sampling_seconds=sampling,propagation_and_point_save_seconds=propagation,metrics_seconds=time.perf_counter()-t,total_seconds=time.perf_counter()-started)
    rr.write_json(dest/'summary.json',summary)
    files=[mc.PGLIB118,B/'operating_point_OP2.npz',H/'coordinates_OP2.npz',H/'model_GMM_K4_OP2.npz',H/'model_TMM_K8_OP2.npz',O/'model_fourth_TMM_K4_OP2.npz',*sorted((O/'modeling_blocks').glob('*.npz'))]
    rr.write_manifest(dest,command=f'OPENBLAS_NUM_THREADS=1 python cloud/tests/run_fourth_distribution.py --method {method}',data_files=files,parameters=summary,random_settings={'entropy':seed,'count':N})
    (dest/'README.md').write_text('完整点集与权重保存在points；重点对象完整排序值和累计概率保存在cdf。完整传播使用原OutputMap。比较表是数值证据，诊断绘图数据不作正式论文图。\n')
    print(json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':main()
