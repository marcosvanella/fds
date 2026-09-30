import struct, sys, numpy as np
def read(fn, load=False):
    f = open(fn,'rb'); assert f.read(8)==b'FDSDUMP1'
    nrec = struct.unpack('i', f.read(4))[0]; recs=[]
    r=-1
    while True:
        r+=1
        if nrec>0 and r>=nrec: break
        h = f.read(16+4*4+8*4+8*4+8)
        if len(h)<104: break  # name,icyc,pred,first,kvar (16+16) ; T,DT,RMIN,RMAX ; X(4) ; nb,nk
        name=h[:16].decode().strip(); icyc,pred,first,kvar=struct.unpack('4i',h[16:32])
        T,DT,rmin,rmax=struct.unpack('4d',h[32:64]); X=struct.unpack('4d',h[64:96]); nb,nk=struct.unpack('2i',h[96:104])
        rec=dict(name=name,icyc=icyc,pred=pred,first=first,kvar=kvar,T=T,DT=DT,X=X,bef={},aft={})
        for lst,n in (('bef',nb),('aft',nk)):
            for a in range(n):
                nm=f.read(16).decode().strip(); rank,=struct.unpack('i',f.read(4)); lb=struct.unpack('4i',f.read(16)); ub=struct.unpack('4i',f.read(16))
                shp=[ub[i]-lb[i]+1 for i in range(rank)]; n_=int(np.prod(shp))
                d=np.frombuffer(f.read(8*n_),dtype='<f8')
                rec[lst][nm]=(lb,ub,rank,d.reshape(shp,order='F') if load else None)
        recs.append(rec)
    return recs
if __name__=='__main__':
    for r in read(sys.argv[1]):
        print(r['name'],r['icyc'],r['pred'],r['first'],'T=%.10g DT=%.10g'%(r['T'],r['DT']),r['X'],len(r['bef']),list(r['aft'].keys()))
