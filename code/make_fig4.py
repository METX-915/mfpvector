"""Regenerate the evidence-tier figure (manuscript Fig 4, file Fig3.pdf) from data/.

  (a) PPARalpha docking scores          <- data/rq3_docking.json
  (b) evidence-tier matrix A_ads x H_receptor, including the PE baseline
  (c) mean tier score by polymer        <- data/mfp_vi.json

Usage:  python code/make_fig4.py
Output: figures/Fig3.pdf (vector) and figures/Fig3_preview.png
"""
import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec

import os
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
SRC=os.path.join(ROOT,"data"); FIG=os.path.join(ROOT,"figures"); os.makedirs(FIG,exist_ok=True)
dock=json.load(open(f"{SRC}/rq3_docking.json"))["best_affinity_kcal_per_mol"]
vi=json.load(open(f"{SRC}/mfp_vi.json"))["MFP_VI"]

PFAS=["PFOA","PFBS","GenX","FTOH62"]; LAB=["PFOA","PFBS","GenX","6:2 FTOH"]
POLY=["PTFE","PVDF","ECTFE","PE"]
C_DOCK=['#DD6E1E','#808080','#E8AC00','#1878B8']
C_POLY={'PTFE':'#808080','PVDF':'#1878B8','ECTFE':'#00A170','PE':'#B0B0B0'}

plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['DejaVu Sans'],
                     'mathtext.fontset':'dejavusans','axes.linewidth':1.1,
                     'xtick.major.width':1.1,'ytick.major.width':1.1})

fig=plt.figure(figsize=(9.0,5.47)); fig.patch.set_facecolor('white')
gs=gridspec.GridSpec(2,2,width_ratios=[1.0,1.28],height_ratios=[1,1],
                     wspace=0.30,hspace=0.42,left=0.085,right=0.90,top=0.93,bottom=0.10)

# ── (a) docking ──────────────────────────────────────────────
ax=fig.add_subplot(gs[0,0])
v=[dock[p] for p in PFAS]
ax.bar(range(4),v,color=C_DOCK,width=0.72,zorder=3)
for i,val in enumerate(v):
    ax.text(i,val+0.30,f"{val:.1f}",ha='center',va='bottom',fontsize=10.5,fontweight='bold',
            color='white' if i in (0,3) else 'black',zorder=4)
ax.set_xticks(range(4)); ax.set_xticklabels(LAB,rotation=28,ha='right',fontsize=10)
ax.set_ylim(-8.3,0); ax.set_ylabel('Docking (kcal/mol)\nmore negative = stronger',fontsize=9.5)
ax.set_title('(a) PPAR$\\alpha$ docking scores',fontsize=12,fontweight='bold',pad=8)
ax.spines[['top','right']].set_visible(False); ax.tick_params(labelsize=9.5)

# ── (c) mean tier ────────────────────────────────────────────
ax=fig.add_subplot(gs[1,0])
means=[np.mean([vi[p][f] for f in PFAS]) for p in POLY]
ax.bar(range(len(POLY)),means,color=[C_POLY[p] for p in POLY],width=0.62,zorder=3)
ax.set_xticks(range(len(POLY))); ax.set_xticklabels(POLY,fontsize=11)
ax.set_ylabel('Mean tier score',fontsize=10.5); ax.set_ylim(0,0.63)
ax.set_title('(c) Mean tier score by polymer',fontsize=12,fontweight='bold',pad=8)
ax.spines[['top','right']].set_visible(False); ax.tick_params(labelsize=10)

# ── (b) evidence-tier matrix ─────────────────────────────────
ax=fig.add_subplot(gs[:,1])
M=np.array([[vi[p][f] for f in PFAS] for p in POLY])
ex=np.arange(-0.5,4.5); ey=np.arange(-0.5,len(POLY)+0.5)
im=ax.pcolormesh(ex,ey,M,cmap='Reds',vmin=0,vmax=1.0,shading='flat',
                 linewidth=0,edgecolors='face',rasterized=False)
ax.set_aspect('auto'); ax.invert_yaxis()
for i in range(len(POLY)):
    for j in range(4):
        ax.text(j,i,f"{M[i,j]:.2f}",ha='center',va='center',fontsize=13,fontweight='bold',
                color='white' if M[i,j]>0.60 else 'black')
ax.set_xticks(range(4)); ax.set_xticklabels(LAB,rotation=28,ha='right',fontsize=11)
ax.set_yticks(range(len(POLY))); ax.set_yticklabels(POLY,fontsize=11)
ax.set_title('(b) Evidence-tier prioritization matrix',fontsize=12,fontweight='bold',pad=8)
cb=fig.colorbar(im,ax=ax,fraction=0.045,pad=0.04)
cb.set_label(r'$\hat{A}_{ads}\times\hat{H}_{receptor}$  (evidence tier)',fontsize=10)
cb.ax.tick_params(labelsize=9.5)

fig.savefig(os.path.join(FIG,'Fig3.pdf'),facecolor='white')
fig.savefig(os.path.join(FIG,'Fig3_preview.png'),dpi=200,facecolor='white')
print("mean tier:",[f"{m:.3f}" for m in means])
print("matrix rows:",[[f"{x:.2f}" for x in r] for r in M])
