"""Full weighted ECDF and exact tail W1 for frozen IEEE118 comparisons."""
import numpy as np

PROBS=np.array([.005,.010,.020,.030,.040,.050,.950,.960,.970,.980,.990,.995])
TAGS=['q005','q010','q020','q030','q040','q050','q950','q960','q970','q980','q990','q995']


def inverse(y,c,p):
    return y[np.minimum(np.searchsorted(c,p,side='left'),len(y)-1)]


def ecdf(y,c,x):
    j=np.searchsorted(y,x,side='right')
    return np.r_[0.,c][j]


def tail_w1(y,c,r,a,b):
    """Exact step-function integral on the union of model and MC probability breaks."""
    n=len(r)
    cm=c[(c>a)&(c<b)]
    cr=np.arange(int(np.floor(a*n))+1,int(np.ceil(b*n)),dtype=float)/n
    cr=cr[(cr>a)&(cr<b)]
    breaks=np.unique(np.r_[a,cm,cr,b])
    mid=(breaks[:-1]+breaks[1:])/2
    qm=inverse(y,c,mid)
    qr=r[np.minimum(np.ceil(mid*n).astype(int)-1,n-1)]
    return float(np.diff(breaks)@np.abs(qm-qr)/(b-a))


def distribution(y,c,r):
    r=np.sort(np.asarray(r,dtype=float))
    out={}
    for side,(a,b) in {'left':(.001,.05),'right':(.95,.999)}.items():
        p=np.linspace(a,b,101)
        x=r[np.minimum(np.ceil(p*len(r)).astype(int)-1,len(r)-1)]
        fr=np.searchsorted(r,x,side='right')/len(r)
        fm=ecdf(y,c,x)
        out[side+'_cdf_rmse']=float(np.sqrt(np.mean((fm-fr)**2)))
        out[side+'_w1']=tail_w1(y,c,r,a,b)
    qm=inverse(y,c,PROBS)
    qr=r[np.minimum(np.ceil(PROBS*len(r)).astype(int)-1,len(r)-1)]
    for tag,v,ref in zip(TAGS,qm,qr):
        out[tag]=float(v)
        out[tag+'_MC']=float(ref)
        out[tag+'_error']=float(v-ref)
    return out
