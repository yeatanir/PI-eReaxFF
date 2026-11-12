import numpy as np, pathlib
BOHR=0.529177210903; KB=8.617333262145179e-5; HBAR=0.6582119569509067; AMU=103.64269667161
P=8; T=300.0; K=1.0; k=0.26; m=1.00794; kT=KB*T
def frames(path):
    L=pathlib.Path(path).read_text().splitlines()
    return np.array([[[float(x) for x in l.split()[1:4]] for l in L[i+2:i+4]] for i in range(0,len(L),4)])
w=2*P*kT/HBAR*np.sin(np.pi*np.arange(1,P)/P)
ref={"X2":3*kT/K,"centroid_to_X2":3*kT/k,"Rg2":3*kT*np.sum(1/(m*AMU*w**2+k))}
vals={x:[] for x in ref}
for s in (1,2,3,4):
    d=pathlib.Path("seed%d"%s)
    pos=np.array([frames(d/("sim.pos_%d.xyz"%b)) for b in range(P)])*BOHR   # (P,frames,2,3) in A
    n=pos.shape[1]; pos=pos[:,n//5:]
    X=pos[0,:,0,:]; x=pos[:,:,1,:]; xc=x.mean(axis=0)
    vals["X2"].append(np.mean(np.sum(X**2,axis=1)))
    vals["centroid_to_X2"].append(np.mean(np.sum((xc-X)**2,axis=1)))
    vals["Rg2"].append(np.mean(np.mean(np.sum((x-xc)**2,axis=2),axis=0)))
    print("seed",s,"frames",n)
for key in ref:
    v=np.array(vals[key]); mm=v.mean(); se=v.std(ddof=1)/np.sqrt(len(v))
    print("i-PI P=8 %-15s mean %.6g +- %.2g  exact %.6g  dev %+.2f %% (%+.1f SE)"%(key,mm,se,ref[key],(mm-ref[key])/ref[key]*100,(mm-ref[key])/se))
