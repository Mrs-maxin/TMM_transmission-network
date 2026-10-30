"""Verify new blocks against the committed original Stage G physical evidence."""
import json,sys,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
import pandas as pd
from tmm_cloud import run_record as rr,montecarlo as m,network as nw,matpower_case as mc,tail_stats as ts
O=rr.CLOUD_ROOT/'results/ieee118_fourth_moment_op2/20261030T120000'
D=rr.CLOUD_ROOT/'results/phase1_rts118/20261030T120000/stage_d_montecarlo/OP2'
OLD=rr.CLOUD_ROOT/'results/phase1_rts118/20261021T202200/stage_d_montecarlo/OP2'

def main():
    started=time.perf_counter()
    old=json.loads((OLD/'stage_d_summary_OP2.json').read_text())
    rows=[]
    for name,info in old['sets'].items():
        src=O/'modeling_blocks' if name=='modeling' else D
        paths=sorted(src.glob(name+'_chunk*.npz'))
        d=m.load_set(paths,['converged','fail_code','iterations','rounds'])
        rk,rc=np.unique(d['rounds'],return_counts=True)
        new={'scenarios':len(d['converged']),'converged':int(d['converged'].sum()),'failed_scenario_index':np.flatnonzero(~d['converged']).tolist(),'iterations_mean':float(d['iterations'].mean()),'iterations_max':int(d['iterations'].max()),'rounds_distribution':{str(int(a)):int(b) for a,b in zip(rk,rc)}}
        for k,v in new.items():
            same=abs(v-info[k])<1e-12 if k=='iterations_mean' else v==info[k]
            rows.append(dict(set=name,quantity=k,original=info[k],current=v,identical=bool(same)))
    net=nw.from_matpower(mc.PGLIB118)
    H=m.load_set(sorted(D.glob('check_heavy_chunk*.npz')),['V','Pf_MW','Qf_Mvar','Pt_MW','Qt_Mvar','imbalance_MW','converged'])
    assert H['converged'].all() and len(H['converged'])==100000
    L=np.maximum(np.hypot(H['Pf_MW'],H['Qf_Mvar']),np.hypot(H['Pt_MW'],H['Qt_Mvar']))/net.rateA[None]
    counts={'any_branch_loading_gt_1.05':int((L>1.05).any(axis=1).sum()),'imbalance_below_-325.75MW':int((H['imbalance_MW'] < -325.75).sum())}
    nr=np.setdiff1d(np.arange(net.n),[net.ref])
    table=pd.read_csv(rr.CLOUD_ROOT/'results/phase1_rts118/20261021T203828/stage_e_tail/tail_statistics_OP2.csv')
    diffs={}
    for k,X in {'Vm':np.abs(H['V'][:,nr]),'Pf':H['Pf_MW'],'Qf':H['Qf_Mvar'],'Pt':H['Pt_MW'],'Qt':H['Qt_Mvar']}.items():
        s=ts.column_statistics(X.astype(float));ref=table[table.quantity==k]
        diffs[k]={f:float(np.max(np.abs(s[f]-ref[f+'_heavy'].values))) for f in ('mean','std','excess_kurtosis','q005','q995')}
    good=all(r['identical'] for r in rows) and counts=={'any_branch_loading_gt_1.05':153,'imbalance_below_-325.75MW':228} and max(v for a in diffs.values() for v in a.values())<1e-9
    result=dict(identical=good,comparison=rows,event_counts=counts,statistic_max_abs_differences=diffs,seconds=time.perf_counter()-started)
    rr.write_json(O/'identity.json',result)
    print(json.dumps(result,indent=2),flush=True)
    assert good,'Regeneration differs from original physical evidence'
if __name__=='__main__':main()
