import collections,sys   # time-mean (t>=1) leaf balance of a long run: usage rad_balance_series.py DIR...
def series(d):
    S=collections.OrderedDict()
    for l in open(d+'/run.log',errors='replace'):
        if l.startswith('RADLEAF'):
            f=l.split(); t=float(f[1]); s=S.setdefault(t,[0,0,0,0]); s[0]-=float(f[4]); s[1]+=float(f[13]); s[2]+=float(f[15]); s[3]+=1
    return S
import numpy as np
for d in sys.argv[1:]:
    S=series(d); ks=[t for t in S if t>=1.0]; a=np.array([S[t][:3] for t in ks])
    print(d,'t>=1: mean loss %.2f walls %.2f open %.2f  resid(mean) %.2f  |resid| rms %.2f  (resid/loss %.3f)'%(a[:,0].mean(),a[:,1].mean(),a[:,2].mean(),(a[:,0]-a[:,1]-a[:,2]).mean(),np.sqrt(((a[:,0]-a[:,1]-a[:,2])**2).mean()),(a[:,0]-a[:,1]-a[:,2]).mean()/a[:,0].mean()))
