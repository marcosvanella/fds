import collections,sys   # time-weighted (trapezoid) mean leaf balance over [T0,T1] with block SE: rad_balance_series.py [--t0 1 --t1 2 --nb 4] DIR...
import numpy as np       # loss = -sum SQ, walls = sum SNA, open = sum SNO, resid = loss-walls-open (kW); same mean as val_compare.py --mode mean
def series(d):
    S=collections.OrderedDict()
    for l in open(d+'/run.log',errors='replace'):
        if l.startswith('RADLEAF'):
            f=l.split()
            if len(f)<16: continue
            t=float(f[1]); s=S.setdefault(t,[0.,0.,0.]); s[0]-=float(f[4]); s[1]+=float(f[13]); s[2]+=float(f[15])
    return S
def tmean(t,v): return np.trapezoid(v,t)/(t[-1]-t[0])
a=sys.argv[1:]; o={'--t0':1.0,'--t1':2.0,'--nb':4}
while a and a[0] in o: o[a[0]]=float(a[1]); a=a[2:]
for d in [x for x in a if __import__("os").path.isfile(x+"/run.log")]:
    S=series(d); t=np.array([x for x in S if o['--t0']<=x<=o['--t1']])
    if len(t)<2: print(d,'no RADLEAF samples in [%g,%g]'%(o['--t0'],o['--t1'])); continue
    v=np.array([S[x] for x in t]); r=v[:,0]-v[:,1]-v[:,2]; m=[tmean(t,v[:,i]) for i in range(3)]; mr=tmean(t,r)
    e=np.linspace(t[0],t[-1],int(o['--nb'])+1); bm=[]
    for i in range(len(e)-1):
        k=(t>=e[i])&(t<=e[i+1])
        if k.sum()>=2: bm.append(tmean(t[k],r[k]))
    se=np.std(bm,ddof=1)/np.sqrt(len(bm)) if len(bm)>1 else float('nan')
    print(d,'t=[%.3f,%.3f] n=%d time-weighted: loss %.3f walls %.3f open %.3f resid %.3f +- %.3f (block SE, %d blocks) resid/loss %.4f'
          %(t[0],t[-1],len(t),m[0],m[1],m[2],mr,se,len(bm),mr/m[0]))
