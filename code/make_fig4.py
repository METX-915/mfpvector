"""Fig 3 of the manuscript (mass-balance scenario; file name kept from earlier numbering).
   Same logic and seed (20260710) as code/rq1_massbalance.py. Output: Fig4.pdf (vector) + Fig4_preview.png"""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from scipy.stats import spearmanr

# ── 원 스크립트와 동일 (seed 20260710, N=200,000) ──────────────
rng=np.random.default_rng(20260710); N=200_000
EOF=rng.uniform(180,780,N)
mp_per_kg=10**rng.uniform(np.log10(140),np.log10(27259),N)
mp_per_g=mp_per_kg/1000.0
mfp_frac=rng.uniform(0.02,0.08,N)
mfp_per_g=mp_per_g*mfp_frac
d_um=10**rng.uniform(np.log10(20),np.log10(300),N); d_cm=d_um*1e-4
vol_cm3=0.3*(np.pi/6.0)*d_cm**3
density=rng.uniform(1.78,2.20,N)
mass_particle_g=vol_cm3*density
fwt=rng.uniform(0.39,0.76,N)
mfp_F=mfp_per_g*mass_particle_g*fwt*1e9
ratio=mfp_F/EOF
med=np.median(ratio)
logout=np.log(mfp_F)
rho={'particle size':spearmanr(np.log(d_um),logout).statistic,
     'total MP abundance':spearmanr(np.log(mp_per_kg),logout).statistic,
     'MFP fraction':spearmanr(np.log(mfp_frac),logout).statistic,
     'F content':spearmanr(np.log(fwt),logout).statistic,
     'density':spearmanr(np.log(density),logout).statistic}

plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['DejaVu Sans'],
                     'mathtext.fontset':'dejavusans','axes.linewidth':1.1,
                     'xtick.major.width':1.1,'ytick.major.width':1.1})
BLUE='#5B9BD5'; BLUE_D='#1F6FB4'; ORANGE='#D2691E'; GREY='#808080'

fig,axes=plt.subplots(1,2,figsize=(10.0,4.66))
fig.subplots_adjust(left=0.075,right=0.975,top=0.90,bottom=0.155,wspace=0.28)

# ── (a) scenario distribution ────────────────────────────────
ax=axes[0]
bins=np.logspace(np.log10(ratio.min()),np.log10(ratio.max()),100)
ax.hist(ratio,bins=bins,color=BLUE,edgecolor='none')
ax.set_xscale('log')
ax.axvline(med,color=BLUE_D,ls='--',lw=2.4,zorder=5)
ax.text(0.60,0.975,f'central\n{med*100:.1f}%',transform=ax.transAxes,color=BLUE_D,
        fontsize=13,fontweight='bold',va='top')
ax.text(0.015,0.975,'order-of-magnitude\nscenario; MFP-F is\nnot part of EOF',
        transform=ax.transAxes,fontsize=9.5,style='italic',color='#7F7F7F',va='top')
ax.set_xlim(2e-6,5e1); ax.set_ylim(0,ax.get_ylim()[1]*1.22)
ax.xaxis.set_major_formatter(FuncFormatter(
    lambda x,_: f"{x*100:g}%" if x>=1e-4 else ""))
ax.set_ylabel('Scenario count',fontsize=12)
ax.set_xlabel('MFP-derived fluorine / EOF reference  (scale comparison)',fontsize=11.5)
ax.set_title('(a) MFP fluorine vs measured EOF scale',fontsize=13.5,fontweight='bold',pad=8)
ax.spines[['top','right']].set_visible(False); ax.tick_params(labelsize=11)

# ── (b) global sensitivity ───────────────────────────────────
ax=axes[1]
names=list(rho); vals=[abs(rho[k]) for k in names]
y=np.arange(len(names))[::-1]
ax.barh(y,vals,color=[ORANGE if k=='particle size' else GREY for k in names],
        height=0.62,edgecolor='black',lw=0.8,zorder=3)
for yy,v in zip(y,vals):
    ax.text(v+0.012,yy,f"{v:.2f}",va='center',fontsize=11.5)
ax.set_yticks(y); ax.set_yticklabels(names,fontsize=12.5)
ax.set_xlim(0,1.0); ax.set_xlabel(r'|Spearman $\rho$| with MFP fluorine',fontsize=12)
ax.set_title('(b) Sensitivity: what to measure first',fontsize=13.5,fontweight='bold',pad=8)
ax.text(0.985,0.30,'reflects prior widths,\nnot natural processes\n'
        r'$\rightarrow$ size-resolved inventory first',
        transform=ax.transAxes,fontsize=10,style='italic',color=ORANGE,ha='right',va='center')
ax.spines[['top','right']].set_visible(False); ax.tick_params(labelsize=11)

fig.savefig('Fig4.pdf',facecolor='white')
fig.savefig('Fig4_preview.png',dpi=200,facecolor='white')
print(f"median ratio = {med:.4f} ({med*100:.1f}%)")
print("rho:", {k:round(abs(v),3) for k,v in rho.items()})
