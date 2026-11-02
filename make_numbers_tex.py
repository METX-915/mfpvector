"""Emit numbers.tex (LaTeX macros) from experiment_summary.json (+ optional pmf json, slab density file).
Usage: numbers.py outdir [pmf.json]"""
import sys,json,os
out=sys.argv[1]; pmfj=sys.argv[2] if len(sys.argv)>2 else None
S=json.load(open(f"{out}/experiment_summary.json")); T=json.load(open(f"{out}/tier.json"))
cell={(s["polymer"],s["pfas"]):s for s in S["systems"]}; PL=S["polymer_level"]
C={(c["contrast"],c["term"]):c for c in S["contrasts"]}; R={r["term"]:r for r in S["ptfe_vs_rest"]}
POLY=["PTFE","PVDF","ECTFE","PE"]; PF=["PFOA","PFBS","GenX","FTOH62"]
key={"PTFE":"ptfe","PVDF":"pvdf","ECTFE":"ectfe","PE":"pe","PFOA":"pfoa","PFBS":"pfbs","GenX":"genx","FTOH62":"ftoh"}
def pfmt(p): return f"{p:.3f}" if p>=0.001 else "<0.001"
def sgn(v,d=1): return f"{v:+.{d}f}"
L=[]
def m(name,val): L.append(f"\\newcommand{{\\{name}}}{{{val}}}")
m("nTraj",S["n_trajectories"])
for p in POLY:
    for f in PF:
        c=cell[(p,f)]; k=key[p]+key[f]
        m("IE"+k,f"{c['meanIE']:.0f}"); m("IEd"+k,f"{c['meanIE']:.1f}"); m("CIlo"+k,f"{c['ci_lo']:.1f}"); m("CIhi"+k,f"{c['ci_hi']:.1f}")
        m("SD"+k,f"{c['sd']:.1f}"); m("Coul"+k,sgn(c['meanCoul'])); m("LJ"+k,f"{c['meanLJ']:.1f}"); m("Ann"+k,f"{c['n']}")
        m("Aads"+k,f"{c['A_ads_norm']:.2f}")
    m("IE"+key[p],f"{PL[p]['mean_reps']:.0f}"); m("IEd"+key[p],f"{PL[p]['mean_reps']:.1f}"); m("IEsd"+key[p],f"{PL[p]['sd_reps']:.1f}")
    m("IEsem"+key[p],f"{PL[p]['sem_reps']:.1f}"); m("IEcellsd"+key[p],f"{PL[p]['sd_of_cells']:.1f}")
    m("Coul"+key[p],sgn(PL[p]['coul'])); m("Coulsd"+key[p],f"{PL[p]['coul_sd']:.1f}"); m("LJ"+key[p],f"{PL[p]['lj']:.1f}"); m("LJsd"+key[p],f"{PL[p]['lj_sd']:.1f}")
    m("Tier"+key[p],f"{T['tier_mean'][POLY.index(p)]:.2f}")
for (a,b) in [("PVDF","PE"),("ECTFE","PE"),("PTFE","PE"),("PVDF","ECTFE"),("PVDF","PTFE"),("ECTFE","PTFE")]:
    for term,tk in (("totalIE","IE"),("coulSR","Coul"),("ljSR","LJ")):
        c=C[(f"{a} vs {b}",term)]; k=tk+key[a]+key[b]
        m("p"+k,pfmt(c["p"])); m("dz"+k,f"{c['dz']:.2f}"); m("adz"+k,f"{abs(c['dz']):.2f}"); m("prep"+k,pfmt(c["p_rep"])); m("cohd"+k,f"{c['cohen_d']:.2f}"); m("delta"+k,sgn(c["delta"]))
for term,tk in (("total","IE"),("coul","Coul"),("lj","LJ")):
    m("pRest"+tk,pfmt(R[term]["p"])); m("dzRest"+tk,f"{R[term]['dz']:.2f}"); m("deltaRest"+tk,sgn(R[term]["delta"]))
    m("pFried"+tk,pfmt(S["friedman_p"][term]))
m("rhoF",f"{S['fluorine_corr']['spearman_rho']:.2f}"); m("pRhoF",f"{S['fluorine_corr']['spearman_p']:.2f}")
for i,p in enumerate(POLY):
    for j,f in enumerate(PF): m("Tier"+key[p]+key[f],f"{T['tier'][i][j]:.2f}")
# slab densities
dens={"ECTFE":(42,3.485,1.4488,1.55),"PE":(65,3.548,0.8038,0.855),"PTFE":(42,3.478,1.8828,1.98),"PVDF":(54,3.639,1.4825,1.68)}
for p,(n,h,r,ref) in dens.items():
    m("Nch"+key[p],n); m("Hslab"+key[p],f"{h:.2f}"); m("Rho"+key[p],f"{r:.3f}"); m("Rhoref"+key[p],f"{ref:.2f}"); m("Rhopct"+key[p],f"{100*r/ref:.0f}")
# PMF
if pmfj and os.path.exists(pmfj):
    P=json.load(open(pmfj))
    for p in POLY:
        m("dG"+key[p],f"{P[p]['dG']:.1f}"); m("dGerr"+key[p],f"{P[p]['err']:.1f}"); m("Ratio"+key[p],f"{abs(cell[(p,'PFOA')]['meanIE']/P[p]['dG']):.1f}")
        m("PMFns"+key[p],f"{P[p].get('ns',5)}")
    m("pmfPending","")
else:
    for p in POLY:
        m("dG"+key[p],"\\textcolor{red}{[dG]}"); m("dGerr"+key[p],"\\textcolor{red}{[err]}"); m("Ratio"+key[p],"\\textcolor{red}{[r]}"); m("PMFns"+key[p],"5")
    m("pmfPending","\\textcolor{red}{[PMF pending]}")
open(f"{out}/numbers.tex","w").write("\n".join(L)+"\n"); print(len(L),"macros")
