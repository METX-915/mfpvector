"""Graphical abstract. Usage: make_ga.py outdir figdir"""
import sys,json,matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Ellipse
out,figdir=sys.argv[1],sys.argv[2]
S=json.load(open(f"{out}/experiment_summary.json")); PL=S["polymer_level"]
order=["PE","ECTFE","PVDF","PTFE"]; wf={"PE":0,"ECTFE":39,"PVDF":59,"PTFE":76}
col={"PE":"#7f8c8d","ECTFE":"#0072B2","PVDF":"#0072B2","PTFE":"#D55E00"}
fig=plt.figure(figsize=(13,5),dpi=300); fig.patch.set_facecolor("white")
fig.text(0.5,0.93,"Fluorine content does not set PFAS sorption to fluoropolymer microplastics",ha="center",va="center",fontsize=17,fontweight="bold",color="#2c3e50")
ax=fig.add_axes([0.04,0.12,0.44,0.7]); ax.axis("off"); ax.set_xlim(-0.6,4.9); ax.set_ylim(0,1.35)
ax.text(-0.5,1.28,"more fluorine →",fontsize=11,color="#7f8c8d")
for i,p in enumerate(order):
    ie=PL[p]["mean_reps"]; h=abs(ie)/70
    ax.add_patch(FancyBboxPatch((i-0.36,0.18),0.72,h,boxstyle="round,pad=0,rounding_size=0.04",fc=col[p],ec="none"))
    ax.text(i,0.18+h+0.04,f"{ie:.0f}",ha="center",fontsize=14,fontweight="bold",color=col[p])
    ax.text(i,0.09,p,ha="center",fontsize=13,fontweight="bold",color="#2c3e50"); ax.text(i,0.0,f"{wf[p]} wt% F",ha="center",fontsize=10,color="#7f8c8d")
co=PL["PTFE"]["coul"]
ax.text(3.45,0.62,"weakest\nsorbent",ha="left",va="center",fontsize=12,fontweight="bold",color="#D55E00")
ax.text(3.45,0.40,f"Coulomb term\n{co:+.0f} kJ mol⁻¹\n(repulsive)",ha="left",va="center",fontsize=10,color="#D55E00")
ax.text(-0.62,0.55,"interaction energy\n(kJ mol⁻¹)",rotation=90,ha="center",va="center",fontsize=10,color="#7f8c8d")
fig.add_artist(plt.Line2D([0.5,0.5],[0.15,0.85],color="#bdc3c7",lw=1))
bx=fig.add_axes([0.53,0.12,0.45,0.7]); bx.axis("off"); bx.set_xlim(0,10); bx.set_ylim(0,3)
bx.text(9.7,2.75,"total organic fluorine (TOF)",ha="right",fontsize=11,color="#7f8c8d")
bx.add_patch(FancyBboxPatch((0.3,1.7),1.0,0.9,boxstyle="round,pad=0,rounding_size=0.08",fc="#0072B2",ec="none"))
bx.add_patch(FancyBboxPatch((1.45,1.7),8.25,0.9,boxstyle="round,pad=0,rounding_size=0.08",fc="#2c3e50",ec="none"))
bx.text(1.8,2.15,"non-extractable",color="white",fontsize=13,fontweight="bold",va="center")
for k in range(4): bx.add_patch(Ellipse((6.7+k*0.85,2.15),0.6,0.3,fc="#009E73",ec="white",lw=1.5))
bx.text(0.8,1.45,"EOF",ha="center",fontsize=12,fontweight="bold",color="#0072B2")
bx.text(8.0,1.45,"fluoropolymer particles",ha="center",fontsize=12,fontweight="bold",color="#009E73")
bx.text(5,0.75,"MFP fluorine ≈ 1.7% of the EOF scale (median)",ha="center",fontsize=12,color="#2c3e50")
bx.text(5,0.3,"set by particle size and abundance, not chemistry",ha="center",fontsize=11,color="#7f8c8d")
fig.savefig(f"{figdir}/GA_jhm.tiff",dpi=300,pil_kwargs={"compression":"tiff_lzw"}); fig.savefig(f"{figdir}/GA_jhm.png",dpi=150); print("ok")
