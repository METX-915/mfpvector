#!/usr/bin/env python3
"""
mfpvector Step4 RQ2 analysis — polymer-PFAS interaction energy → relative adsorption.
Input : data/raw_results.csv (polymer,pfas,seed,coulSR,ljSR,totalIE)
Output: out/experiment_summary.json, out/results_table.csv,
        out/Fig2_adsorption.png
Stats : bootstrap 95% CI (10000), Cohen's d (MFP vs PE per PFAS), Shapiro→t/Mann-Whitney.
NOTE  : no absolute claim is made - relative ranking plus uncertainty only.
        A_ads = min-max normalisation (more negative IE = stronger adsorption = 1).
"""
import os, json, numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt

# paths resolved relative to this file so the script runs from a fresh clone
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
DATA=os.path.join(ROOT,"data"); R=os.path.join(ROOT,"out"); os.makedirs(R,exist_ok=True)
df=pd.read_csv(os.path.join(DATA,"raw_results.csv"))
df=df.dropna(subset=["totalIE_kJmol"])
# PP was in the original run plan but produced no completed trajectories;
# it is kept here only so the list matches the driver script, and is skipped below.
POLYS=["PTFE","PVDF","ECTFE","PE","PP"]; PFAS=["PFOA","PFBS","GenX","FTOH62"]
FLUORO=["PTFE","PVDF","ECTFE"]
rng=np.random.default_rng(7)

def boot(x,n=10000):
    x=np.asarray(x); bs=[np.mean(rng.choice(x,len(x),replace=True)) for _ in range(n)]
    return np.mean(x), np.percentile(bs,2.5), np.percentile(bs,97.5)
def cohens_d(a,b):
    a,b=np.asarray(a),np.asarray(b); na,nb=len(a),len(b)
    sp=np.sqrt(((na-1)*a.var(ddof=1)+(nb-1)*b.var(ddof=1))/(na+nb-2))
    return (a.mean()-b.mean())/sp if sp>0 else np.nan

# per-system stats
rows=[]
for p in POLYS:
    for fa in PFAS:
        s=df[(df.polymer==p)&(df.pfas==fa)]["totalIE_kJmol"].values
        if len(s)==0: continue
        m,lo,hi=boot(s)
        rows.append(dict(polymer=p,pfas=fa,n=len(s),meanIE=m,ci_lo=lo,ci_hi=hi,sd=s.std(ddof=1) if len(s)>1 else 0))
summ=pd.DataFrame(rows)
summ.to_csv(f"{R}/results_table.csv",index=False)

# sanity checks
san={}
san["all_adsorptive(meanIE<0)"]=bool((summ.meanIE<0).all())
san["PE_baseline_present"]=bool((summ.polymer=="PE").any())

# effect sizes: each fluoropolymer vs PE, per PFAS
eff=[]
for p in FLUORO:
    for fa in PFAS:
        a=df[(df.polymer==p)&(df.pfas==fa)]["totalIE_kJmol"].values
        b=df[(df.polymer=="PE")&(df.pfas==fa)]["totalIE_kJmol"].values
        if len(a)<2 or len(b)<2: continue
        norm=stats.shapiro(np.concatenate([a,b]))[1]>0.05 if len(a)+len(b)>=3 else True
        if norm: test="t-test"; pval=stats.ttest_ind(a,b,equal_var=False)[1]
        else:    test="Mann-Whitney"; pval=stats.mannwhitneyu(a,b,alternative="two-sided")[1]
        eff.append(dict(polymer=p,pfas=fa,cohens_d=cohens_d(a,b),test=test,p_value=pval))
effdf=pd.DataFrame(eff)

# Â_ads normalization (more negative IE -> stronger -> 1)
if len(summ):
    ie=summ.meanIE.values; A=(ie.max()-ie)/(ie.max()-ie.min()) if ie.max()!=ie.min() else np.zeros_like(ie)
    summ["A_ads_norm"]=A

summary=dict(
  n_systems=int(summ.groupby(["polymer","pfas"]).ngroups),
  seeds_per_system=int(df.groupby(["polymer","pfas"]).size().median()) if len(df) else 0,
  sanity=san,
  systems=summ.to_dict("records"),
  effect_sizes_vs_PE=effdf.to_dict("records"),
)
json.dump(summary, open(f"{R}/experiment_summary.json","w"), indent=2)

# ---- Figure (Ku style) ----
matplotlib.rcParams.update({'font.size':8,'axes.titlesize':9,'axes.labelsize':8,
 'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'font.family':'sans-serif',
 'font.sans-serif':['DejaVu Sans'],'axes.spines.top':False,'axes.spines.right':False,
 'axes.linewidth':0.8,'legend.frameon':False,'savefig.dpi':300})
OK=['#E69F00','#56B4E9','#009E73','#0072B2','#D55E00','#CC79A7']
fig,ax=plt.subplots(1,2,figsize=(7.2,3.0),constrained_layout=True)
# panel a: heatmap meanIE
piv=summ.pivot(index="polymer",columns="pfas",values="meanIE").reindex(index=POLYS,columns=PFAS)
im=ax[0].imshow(piv.values,cmap="viridis",aspect="auto")
ax[0].set_xticks(range(len(PFAS))); ax[0].set_xticklabels(PFAS,rotation=30,ha="right")
ax[0].set_yticks(range(len(POLYS))); ax[0].set_yticklabels(POLYS)
for i in range(piv.shape[0]):
    for j in range(piv.shape[1]):
        v=piv.values[i,j]
        if not np.isnan(v): ax[0].text(j,i,f"{v:.0f}",ha="center",va="center",color="w",fontsize=7)
cb=fig.colorbar(im,ax=ax[0],shrink=0.8); cb.set_label("Interaction energy (kJ/mol)")
ax[0].set_title("(a) Polymer–PFAS interaction energy")
# panel b: bar with CI per polymer (avg over PFAS)
pm=summ.groupby("polymer")["meanIE"].mean().reindex(POLYS)
plo=summ.groupby("polymer")["ci_lo"].mean().reindex(POLYS)
phi=summ.groupby("polymer")["ci_hi"].mean().reindex(POLYS)
x=np.arange(len(POLYS))
ax[1].bar(x,pm.values,color=[OK[0] if p in FLUORO else OK[1] for p in POLYS])
ax[1].errorbar(x,pm.values,yerr=[pm.values-plo.values,phi.values-pm.values],fmt="none",ecolor="k",capsize=3,lw=0.8)
ax[1].set_xticks(x); ax[1].set_xticklabels(POLYS,rotation=30,ha="right")
ax[1].set_ylabel("Mean interaction energy (kJ/mol)")
ax[1].set_title("(b) Adsorption strength by polymer")
ax[1].axhline(0,color="grey",lw=0.6)
fig.savefig(f"{R}/Fig2_adsorption.png",dpi=300,bbox_inches="tight")
print("WROTE:", f"{R}/experiment_summary.json", f"{R}/results_table.csv", f"{R}/Fig2_adsorption.png")
print("\nSanity:",san)
print("\nPer-system mean IE (kJ/mol):"); print(summ[["polymer","pfas","n","meanIE","ci_lo","ci_hi"]].to_string(index=False))
print("\nEffect sizes vs PE:"); print(effdf.to_string(index=False) if len(effdf) else "  (need PE + fluoro data)")
