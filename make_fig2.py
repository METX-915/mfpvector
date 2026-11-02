"""Fig 2 (adsorption screen) for the JHM revision.
Usage: python make_fig2.py <outdir-with-experiment_summary.json> <raw.csv> <figdir> [pmf_summary.json]
(a) polymer x compound mean IE heatmap
(b) polymer-level Coulomb / Lennard-Jones decomposition (stacked) with total IE marker
(c) all replicates
(d) IE vs umbrella-sampling dG_ads for the four PFOA systems (placeholder if PMF json absent)
"""
import sys, os, json, csv, collections, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D; from matplotlib.patches import Patch
out, raw, figdir = sys.argv[1], sys.argv[2], sys.argv[3]; pmfj = sys.argv[4] if len(sys.argv)>4 else None
os.makedirs(figdir, exist_ok=True)
plt.rcParams.update({"font.size":13,"axes.spines.top":False,"axes.spines.right":False,"font.family":"DejaVu Sans",
                     "axes.titleweight":"bold","pdf.fonttype":42,"ps.fonttype":42})
BLUE,GREEN,ORANGE,GREY,RED="#0072B2","#009E73","#D55E00","#7f7f7f","#CC3311"
S=json.load(open(f"{out}/experiment_summary.json"))
POLY=["PTFE","PVDF","ECTFE","PE"]; comps=["PFOA","PFBS","GenX","6:2 FTOH"]; key={"PFOA":"PFOA","PFBS":"PFBS","GenX":"GenX","6:2 FTOH":"FTOH62"}
cell={(s["polymer"],s["pfas"]):s for s in S["systems"]}
M=np.array([[cell[(p,key[c])]["meanIE"] for c in comps] for p in POLY])
rep=collections.defaultdict(list)
for r in csv.DictReader(open(raw)): rep[(r["polymer"],r["pfas"])].append(float(r["totalIE_kJmol"]))
N=sum(len(v) for v in rep.values())
PL=S["polymer_level"]

fig,ax=plt.subplots(2,2,figsize=(13,10)); fig.subplots_adjust(hspace=0.42,wspace=0.40)
# (a)
e=np.arange(-0.5,4.5); im=ax[0,0].pcolormesh(e,e,M,cmap="viridis",shading="flat",edgecolors="face",linewidth=0)
ax[0,0].invert_yaxis(); ax[0,0].set_xticks(range(4)); ax[0,0].set_xticklabels(comps,rotation=25,ha="right")
ax[0,0].set_yticks(range(4)); ax[0,0].set_yticklabels(POLY)
for i in range(4):
    for j in range(4):
        v=M[i,j]; ax[0,0].text(j,i,f"{v:.0f}",ha="center",va="center",color="white" if v<-50 else "black",fontweight="bold")
cb=fig.colorbar(im,ax=ax[0,0],fraction=0.046,pad=0.04); cb.set_label("IE (kJ/mol)")
ax[0,0].set_title("(a) Polymer–PFAS interaction energy"); ax[0,0].spines[:].set_visible(False)
# (b) Coulomb / LJ decomposition
axb=ax[0,1]; x=np.arange(4)
lj=[PL[p]["lj"] for p in POLY]; co=[PL[p]["coul"] for p in POLY]; tot=[PL[p]["mean_reps"] for p in POLY]; sem=[PL[p]["sem_reps"] for p in POLY]
axb.bar(x,lj,0.6,color=BLUE,edgecolor="black",linewidth=0.6,label="Lennard-Jones")
axb.bar(x,co,0.6,color=RED,edgecolor="black",linewidth=0.6,label="Coulomb")
axb.errorbar(x,tot,yerr=sem,fmt="D",color="black",ms=7,capsize=4,lw=1.2,label="Total IE (mean ± s.e.m.)",zorder=5)
axb.axhline(0,color="grey",lw=0.8); axb.set_xticks(x); axb.set_xticklabels(POLY)
axb.set_ylabel("Energy (kJ/mol)"); axb.set_title("(b) Coulomb / Lennard-Jones decomposition"); axb.set_ylim(-80,30)
for i,p in enumerate(POLY):
    axb.text(i,co[i]+ (2.0 if co[i]>=0 else -2.0),f"{co[i]:+.1f}",ha="center",va="bottom" if co[i]>=0 else "top",fontsize=11,color=RED,fontweight="bold")
axb.legend(loc="lower right",fontsize=11,frameon=False); axb.set_ylim(-88,30)
# (c) replicates
MK={"PFOA":"o","PFBS":"s","GenX":"^","6:2 FTOH":"D"}; rng=np.random.default_rng(11); axc=ax[1,0]
for i,p in enumerate(POLY):
    for k,c in enumerate(comps):
        v=np.array(rep[(p,key[c])]); xx=i+(k-1.5)*0.16+rng.uniform(-0.03,0.03,len(v))
        axc.scatter(xx,v,marker=MK[c],s=42,facecolor=BLUE,edgecolor="black",linewidth=0.5,alpha=0.85,zorder=3)
    axc.plot([i-0.34,i+0.34],[tot[i]]*2,color="black",lw=2.4,zorder=5)
axc.axhline(0,color="grey",lw=0.8,ls="--"); axc.set_xticks(range(4)); axc.set_xticklabels(POLY)
axc.set_ylabel("Replicate IE (kJ/mol)"); axc.set_xlim(-0.6,3.6); axc.set_ylim(-105,22); axc.set_yticks(range(-100,1,20))
axc.set_title(f"(c) Per-replicate spread ({N} trajectories)")
axc.legend(handles=[Line2D([0],[0],marker=MK[c],ls="",mfc=BLUE,mec="black",ms=8,label=c) for c in comps]+
           [Line2D([0],[0],color="black",lw=2.4,label="polymer mean")],loc="upper center",fontsize=11,frameon=False,ncol=5,handletextpad=0.3,columnspacing=0.9,bbox_to_anchor=(0.5,1.0))
# (d) IE vs dG for PFOA on all four polymers
axd=ax[1,1]; ie=[cell[(p,"PFOA")]["meanIE"] for p in POLY]
P=json.load(open(pmfj)) if pmfj and os.path.exists(pmfj) else {}
dg=[P.get(p,{}).get("dG") for p in POLY]; dge=[P.get(p,{}).get("err") or 0 for p in POLY]
w=0.36
b1=axd.bar(x-w/2,ie,w,color=BLUE,edgecolor="black",linewidth=0.6,label="Short-range IE (PFOA)")
have=[i for i in range(4) if dg[i] is not None]; miss=[i for i in range(4) if dg[i] is None]
if have:
    axd.bar([x[i]+w/2 for i in have],[dg[i] for i in have],w,yerr=[dge[i] for i in have],capsize=4,color=GREEN,edgecolor="black",linewidth=0.6,label="PMF $\\Delta G_{\\mathrm{ads}}$")
    for i in have: axd.text(i+w/2,dg[i]-dge[i]-2,f"{abs(ie[i]/dg[i]):.1f}×",ha="center",va="top",fontsize=11.5,fontweight="bold",color=ORANGE)
if miss:
    axd.bar([x[i]+w/2 for i in miss],[-10]*len(miss),w,color="none",edgecolor=RED,linewidth=1.2,ls="--",label="PMF (re-run pending)" if not have else None)
    for i in miss: axd.text(i+w/2,-14,"re-run",ha="center",va="top",color=RED,fontsize=11)
for r in b1: axd.text(r.get_x()+r.get_width()/2,r.get_height()+2.4,f"{r.get_height():.0f}",ha="center",va="bottom",fontsize=11,fontweight="bold",color="white")
axd.axhline(0,color="grey",lw=0.8); axd.set_xticks(x); axd.set_xticklabels([f"{p}\n+ PFOA" for p in POLY],fontsize=12)
axd.set_ylabel("Energy (kJ/mol)"); axd.set_ylim(-85,28); axd.set_yticks(range(-80,1,20)); axd.set_title("(d) Interaction energy is not a free energy")
axd.legend(loc="upper center",fontsize=11,frameon=False,ncol=2,bbox_to_anchor=(0.5,1.0))
for ext,kw in (("pdf",{}),("png",{"dpi":300})):
    fig.savefig(f"{figdir}/Fig2.{ext}",bbox_inches="tight",**kw); print("wrote",f"{figdir}/Fig2.{ext}")
