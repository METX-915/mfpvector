import sys, re, numpy as np
R={'C':0.170,'F':0.147,'H':0.120,'Cl':0.175}; PROBE=0.14
def read_itp(p):
    txt=open(p).read(); sec=txt.split('[ atoms ]')[1].split('[')[0]
    q={}
    for l in sec.splitlines():
        l=l.split(';')[0].split()
        if len(l)>=7: q[l[4]]=float(l[6])
    return q
def read_gro(p):
    L=open(p).read().splitlines(); n=int(L[1]); rows=L[2:2+n]
    names=[r[10:15].strip() for r in rows]; xyz=np.array([[float(r[20:28]),float(r[28:36]),float(r[36:44])] for r in rows])
    box=np.array([float(x) for x in L[2+n].split()[:3]]); return names,xyz,box
def sphere(n=240):
    i=np.arange(n)+0.5; phi=np.arccos(1-2*i/n); th=np.pi*(1+5**0.5)*i
    return np.c_[np.cos(th)*np.sin(phi),np.sin(th)*np.sin(phi),np.cos(phi)]
def run(gro,itp,label):
    names,x,box=read_gro(gro); q=read_itp(itp)
    el=np.array(['Cl' if n.upper().startswith('CL') else n[0] for n in names]); chg=np.array([q[n] for n in names])
    rad=np.array([R[e] for e in el])+PROBE; S=sphere(); acc=np.zeros(len(x))
    from scipy.spatial import cKDTree
    # periodic in x,y only: replicate neighbours
    shifts=[(i,j) for i in (-1,0,1) for j in (-1,0,1)]
    X=np.vstack([x+[i*box[0],j*box[1],0] for i,j in shifts]); RR=np.tile(rad,9)
    T=cKDTree(X); rmax=RR.max()
    for k in range(len(x)):
        pts=x[k]+rad[k]*S
        nb=[m for m in T.query_ball_point(x[k],rad[k]+rmax) if m!=k]
        if not nb: acc[k]=4*np.pi*rad[k]**2; continue
        d=np.linalg.norm(pts[:,None,:]-X[nb][None,:,:],axis=2)
        free=~(d<RR[nb][None,:]).any(1); acc[k]=4*np.pi*rad[k]**2*free.mean()
    area=2*box[0]*box[1]  # two faces
    tot=acc.sum()
    out={'label':label,'faces_nm2':area,'SASA_total_nm2':tot}
    for e in sorted(set(el)):
        m=el==e; out[f'{e}_area_frac']=acc[m].sum()/tot; out[f'{e}_exposed_per_nm2']=(acc[m]>0.01).sum()/area
    out['area_weighted_mean_charge_e']=(acc*chg).sum()/tot
    out['positive_area_frac']=acc[chg>0].sum()/tot
    out['area_weighted_mean_charge_pos_e']=(acc[chg>0]*chg[chg>0]).sum()/max(acc[chg>0].sum(),1e-9)
    out['area_weighted_mean_charge_neg_e']=(acc[chg<0]*chg[chg<0]).sum()/acc[chg<0].sum()
    return out

if __name__ == "__main__":
    # usage: python code/sasa.py <slab.gro> <polymer.itp> <label>
    o = run(sys.argv[1], sys.argv[2], sys.argv[3])
    print({k: (round(float(v), 4) if not isinstance(v, str) else v) for k, v in o.items()})
