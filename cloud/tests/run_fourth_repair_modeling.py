"""Recreate only the unreadable modeling block with unchanged Stage D draws."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from tmm_cloud import error_model as em,montecarlo as m,run_record as rr
b=rr.latest_stage_dir('phase1_rts118','stage_b_network')
c=rr.latest_stage_dir('phase1_rts118','stage_c_error_model')
op_path=b/'operating_point_OP2.npz'
op=m.load_operating_point(op_path)
model=em.from_npz(np.load(c/'error_model_OP2.npz'))
rng=np.random.default_rng(np.random.SeedSequence(2026102112).spawn(4)[0])
e,w=em.sample_errors(model,5256,rng,'heavy')
Pren=(op.Pren[None]*100+e)/100
kload=1+0.01*rng.standard_normal((5256,len(op.Pd)))
out=rr.CLOUD_ROOT/'results/ieee118_fourth_moment_op2/20261030T120000/modeling_blocks'
out.mkdir(parents=True,exist_ok=True)
info=m.run_set(op_path,Pren,kload,4,1000,out,'modeling',{'strength_w':w.astype(np.float32)},single=False)
rr.write_json(out/'summary.json',info)
