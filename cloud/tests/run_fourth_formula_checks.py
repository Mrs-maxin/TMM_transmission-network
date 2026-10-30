"""Deterministic moment, ECDF, quantile and exact W1 checks, not an AC study."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from tmm_cloud.gmm_fourth_moment_tmm import df_from_radial_fourth,covariance_to_scatter_factor
from tmm_cloud.fourth_comparison_metrics import tail_w1,distribution,ecdf,inverse
for d in (1,2,157):
    for nu in (4.5,5.,8.,20.,200.):
        beta=d*(d+2)*(nu-2)/(nu-4)
        got,_=df_from_radial_fourth(beta,d)
        assert np.isclose(got,nu,rtol=1e-12)
        assert np.isclose(covariance_to_scatter_factor(got)*nu/(nu-2),1)
    assert df_from_radial_fourth(d*(d+2),d)[0]==np.inf
    assert df_from_radial_fourth(.9*d*(d+2),d)[0]==np.inf
assert df_from_radial_fourth(9.,1)[0]==5.
y=np.array([0.,0.,2.,5.]);c=np.arange(1,5)/4
assert ecdf(y,c,0.)==.5
assert inverse(y,c,.5)==0.
assert tail_w1(y,c,y,.001,.05)==0.
assert np.isclose(tail_w1(y+3,c,y,.95,.999),3.)
for a,b in ((.001,.05),(.95,.999)):
    r=np.array([-2.,0.,0.,1.,4.]);w=np.array([.1,.4,.5]);v=np.array([-1.,2.,6.]);cw=w.cumsum()
    grid=(np.arange(200000)+.5)/200000*(b-a)+a
    approx=np.mean(np.abs(inverse(v,cw,grid)-r[np.ceil(grid*len(r)).astype(int)-1]))
    assert abs(tail_w1(v,cw,r,a,b)-approx)<1e-4
print('PASS: 15 known t moments, Gaussian boundary, univariate kurtosis, ties, generalized inverse, exact tail W1')
