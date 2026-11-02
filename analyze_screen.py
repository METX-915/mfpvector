#!/usr/bin/env python3
"""Recompute all screen statistics for the JHM revision.
Usage: python stats.py raw.csv outdir
Writes experiment_summary.json, results_table.csv, and LaTeX table snippets.
"""
import sys, os, json, numpy as np, pandas as pd
from scipy import stats
src, out = sys.argv[1], sys.argv[2]; os.makedirs(out, exist_ok=True)
df = pd.read_csv(src).dropna(subset=["totalIE_kJmol"])
POLY = ["PTFE","PVDF","ECTFE","PE"]; PF = ["PFOA","PFBS","GenX","FTOH62"]
LAB = {"PFOA":"PFOA","PFBS":"PFBS","GenX":"GenX","FTOH62":"6:2 FTOH"}
rng = np.random.default_rng(7)
def boot(x, n=10000):
    x = np.asarray(x); bs = [rng.choice(x, len(x), replace=True).mean() for _ in range(n)]
    return float(x.mean()), float(np.percentile(bs,2.5)), float(np.percentile(bs,97.5))
def cohd(a,b):
    a,b=np.asarray(a),np.asarray(b); sp=np.sqrt(((len(a)-1)*a.var(ddof=1)+(len(b)-1)*b.var(ddof=1))/(len(a)+len(b)-2))
    return float((a.mean()-b.mean())/sp)
# ---- per-cell
rows=[]
for p in POLY:
    for f in PF:
        s=df[(df.polymer==p)&(df.pfas==f)]
        m,lo,hi=boot(s.totalIE_kJmol.values)
        rows.append(dict(polymer=p,pfas=f,n=len(s),meanIE=m,ci_lo=lo,ci_hi=hi,sd=float(s.totalIE_kJmol.std(ddof=1)),
                         meanCoul=float(s.coulSR_kJmol.mean()),sdCoul=float(s.coulSR_kJmol.std(ddof=1)),
                         meanLJ=float(s.ljSR_kJmol.mean()),sdLJ=float(s.ljSR_kJmol.std(ddof=1)),
                         n_desorbed=int((s.totalIE_kJmol.abs()<1.0).sum())))
summ=pd.DataFrame(rows); summ.to_csv(f"{out}/results_table.csv",index=False)
ie=summ.meanIE.values; summ["A_ads_norm"]=(ie.max()-ie)/(ie.max()-ie.min())
summ["A_ads_norm"]=0.05+0.95*summ["A_ads_norm"]  # 0.05 floor as in manuscript
# ---- polymer level (mean over 4 compound means), and over 20 replicates
poly={}
for p in POLY:
    c=summ[summ.polymer==p]; r=df[df.polymer==p]
    poly[p]=dict(mean_of_cells=float(c.meanIE.mean()),sd_of_cells=float(c.meanIE.std(ddof=1)),
                 mean_reps=float(r.totalIE_kJmol.mean()),sd_reps=float(r.totalIE_kJmol.std(ddof=1)),
                 sem_reps=float(r.totalIE_kJmol.std(ddof=1)/np.sqrt(len(r))),n=len(r),
                 coul=float(r.coulSR_kJmol.mean()),coul_sd=float(r.coulSR_kJmol.std(ddof=1)),
                 lj=float(r.ljSR_kJmol.mean()),lj_sd=float(r.ljSR_kJmol.std(ddof=1)),
                 coul_cells=float(c.meanCoul.mean()),lj_cells=float(c.meanLJ.mean()))
# ---- pairwise contrasts, compound-blocked paired t (df=3) and replicate-level Welch
def blocked(p1,p2,col):
    a=summ[summ.polymer==p1].set_index("pfas").loc[PF,col].values
    b=summ[summ.polymer==p2].set_index("pfas").loc[PF,col].values
    d=a-b; t,pv=stats.ttest_rel(a,b); return dict(delta=float(d.mean()),t=float(t),p=float(pv),dz=float(d.mean()/d.std(ddof=1)))
def replev(p1,p2,col):
    a=df[df.polymer==p1][col].values; b=df[df.polymer==p2][col].values
    return dict(p=float(stats.ttest_ind(a,b,equal_var=False)[1]),d=cohd(a,b),
                p_mw=float(stats.mannwhitneyu(a,b,alternative="two-sided")[1]))
pairs=[("PVDF","PE"),("ECTFE","PE"),("PTFE","PE"),("PVDF","ECTFE"),("PVDF","PTFE"),("ECTFE","PTFE")]
contr=[]
for p1,p2 in pairs:
    for col,name in (("meanIE","totalIE_kJmol"),("meanCoul","coulSR_kJmol"),("meanLJ","ljSR_kJmol")):
        b=blocked(p1,p2,col); r=replev(p1,p2,name)
        contr.append(dict(contrast=f"{p1} vs {p2}",term=name.split("_")[0],**b,p_rep=r["p"],cohen_d=r["d"],p_mw=r["p_mw"]))
contr=pd.DataFrame(contr)
# PTFE vs mean of other three (blocked over compounds)
ptfe_vs_rest=[]
for col,name in (("meanIE","total"),("meanCoul","coul"),("meanLJ","lj")):
    a=summ[summ.polymer=="PTFE"].set_index("pfas").loc[PF,col].values
    b=np.mean([summ[summ.polymer==q].set_index("pfas").loc[PF,col].values for q in ["PVDF","ECTFE","PE"]],axis=0)
    d=a-b; t,pv=stats.ttest_rel(a,b)
    ptfe_vs_rest.append(dict(term=name,delta=float(d.mean()),t=float(t),p=float(pv),dz=float(d.mean()/d.std(ddof=1))))
# Friedman across 4 polymers with compound as block (IE, Coul)
fried={}
for col,name in (("meanIE","total"),("meanCoul","coul"),("meanLJ","lj")):
    arr=[summ[summ.polymer==p].set_index("pfas").loc[PF,col].values for p in POLY]
    fried[name]=float(stats.friedmanchisquare(*arr)[1])
# Fluoro vs PE per compound (Welch), as in original analyze_rq2
eff=[]
for p in ["PTFE","PVDF","ECTFE"]:
    for f in PF:
        a=df[(df.polymer==p)&(df.pfas==f)].totalIE_kJmol.values; b=df[(df.polymer=="PE")&(df.pfas==f)].totalIE_kJmol.values
        eff.append(dict(polymer=p,pfas=f,cohens_d=cohd(a,b),p_welch=float(stats.ttest_ind(a,b,equal_var=False)[1])))
# fluorine content correlation (wt% F): PTFE 76.0, PVDF 59.4, ECTFE 39.3 (C4H4F4Cl... ETFE-like), PE 0
wtF={"PTFE":76.0,"PVDF":59.4,"ECTFE":39.3,"PE":0.0}
x=[wtF[p] for p in POLY]; y=[poly[p]["mean_of_cells"] for p in POLY]
rho,prho=stats.spearmanr(x,y); rp,pp=stats.pearsonr(x,y)
summary=dict(source=os.path.basename(src),n_trajectories=int(len(df)),n_desorbed=int((df.totalIE_kJmol.abs()<1.0).sum()),
             systems=summ.to_dict("records"),polymer_level=poly,contrasts=contr.to_dict("records"),
             ptfe_vs_rest=ptfe_vs_rest,friedman_p=fried,effect_sizes_vs_PE=eff,
             fluorine_corr=dict(wtF=wtF,spearman_rho=float(rho),spearman_p=float(prho),pearson_r=float(rp),pearson_p=float(pp)))
json.dump(summary,open(f"{out}/experiment_summary.json","w"),indent=2)
# ---- console report
pd.set_option("display.width",200)
print(f"n={len(df)}  desorbed(|IE|<1)={summary['n_desorbed']}")
print("\nPer-cell:"); print(summ[["polymer","pfas","n","meanIE","ci_lo","ci_hi","sd","meanCoul","meanLJ"]].round(1).to_string(index=False))
print("\nPolymer level:")
for p in POLY: print(f"  {p:6s} cells {poly[p]['mean_of_cells']:6.1f}±{poly[p]['sd_of_cells']:4.1f} | reps {poly[p]['mean_reps']:6.1f}±{poly[p]['sd_reps']:4.1f} (sem {poly[p]['sem_reps']:.1f}) | Coul {poly[p]['coul']:6.2f}±{poly[p]['coul_sd']:.1f} LJ {poly[p]['lj']:6.1f}±{poly[p]['lj_sd']:.1f}")
print("\nContrasts (blocked paired t, df=3; replicate Welch):"); print(contr.round(4).to_string(index=False))
print("\nPTFE vs mean(other 3):"); print(pd.DataFrame(ptfe_vs_rest).round(4).to_string(index=False))
print("\nFriedman p:",{k:round(v,4) for k,v in fried.items()})
print(f"\nwt%F vs polymer mean IE: Spearman rho={rho:.2f} p={prho:.2f}; Pearson r={rp:.2f} p={pp:.2f}")
