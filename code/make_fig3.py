"""Fig 3 (docking + evidence-tier matrix) from experiment_summary.json. Usage: make_fig3.py outdir figdir"""
import sys,json,os,numpy as np,matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
out,figdir=sys.argv[1],sys.argv[2]; os.makedirs(figdir,exist_ok=True)
plt.rcParams.update({"font.size":12,"axes.spines.top":False,"axes.spines.right":False,"font.family":"DejaVu Sans","axes.titleweight":"bold","pdf.fonttype":42,"ps.fonttype":42})
S=json.load(open(f"{out}/experiment_summary.json")); cell={(s["polymer"],s["pfas"]):s for s in S["systems"]}
POLY=["PTFE","PVDF","ECTFE","PE"]; comps=["PFOA","PFBS","GenX","6:2 FTOH"]; key={"PFOA":"PFOA","PFBS":"PFBS","GenX":"GenX","6:2 FTOH":"FTOH62"}
dock={"PFOA":-7.7,"PFBS":-6.5,"GenX":-7.1,"6:2 FTOH":-7.8}; H={"PFOA":0.92,"PFBS":0.20,"GenX":0.46,"6:2 FTOH":1.00}
A=np.array([[cell[(p,key[c])]["A_ads_norm"] for c in comps] for p in POLY]); T=A*np.array([H[c] for c in comps])
json.dump(dict(A_ads=A.round(3).tolist(),tier=T.round(3).tolist(),tier_mean=T.mean(1).round(3).tolist(),polymers=POLY,compounds=comps),open(f"{out}/tier.json","w"),indent=1)
fig=plt.figure(figsize=(9.5,5.6)); gs=fig.add_gridspec(2,2,width_ratios=[1,1.5],hspace=0.55,wspace=0.35)
a=fig.add_subplot(gs[0,0]); b=fig.add_subplot(gs[:,1]); c=fig.add_subplot(gs[1,0])
cols={"PFOA":"#E0781E","PFBS":"#8C8C8C","GenX":"#E8A33D","6:2 FTOH":"#0072B2"}
a.bar(range(4),[dock[x] for x in comps],color=[cols[x] for x in comps],edgecolor="black",lw=0.5)
for i,x in enumerate(comps): a.text(i,dock[x]+0.35,f"{dock[x]:.1f}",ha="center",color="white",fontweight="bold",fontsize=10)
a.set_xticks(range(4)); a.set_xticklabels(comps,rotation=25,ha="right",fontsize=10); a.set_ylabel("Docking (kcal/mol)\nmore negative = stronger",fontsize=9)
a.set_title("(a) PPARα docking scores",fontsize=11); a.axhline(0,color="grey",lw=0.6)
im=b.imshow(T,cmap="Reds",vmin=0,vmax=1,aspect="auto")
for i in range(4):
    for j in range(4): b.text(j,i,f"{T[i,j]:.2f}",ha="center",va="center",fontweight="bold",color="white" if T[i,j]>0.55 else "black")
b.set_xticks(range(4)); b.set_xticklabels(comps,rotation=25,ha="right"); b.set_yticks(range(4)); b.set_yticklabels(POLY); b.spines[:].set_visible(False)
b.set_title("(b) Evidence-tier prioritization matrix",fontsize=11)
cb=fig.colorbar(im,ax=b,fraction=0.046,pad=0.04); cb.set_label("$\\hat{A}_{\\mathrm{ads}}\\times\\hat{H}_{\\mathrm{receptor}}$ (evidence tier)",fontsize=9)
c.bar(range(4),T.mean(1),color=["#0072B2","#0072B2","#009E73","#8C8C8C"],edgecolor="black",lw=0.5)
c.set_xticks(range(4)); c.set_xticklabels(POLY,fontsize=10); c.set_ylabel("Mean tier score",fontsize=10); c.set_ylim(0,0.72); c.set_title("(c) Mean tier score by polymer",fontsize=11)
for ext,kw in (("pdf",{}),("png",{"dpi":300})): fig.savefig(f"{figdir}/Fig3.{ext}",bbox_inches="tight",**kw)
print(json.dumps(dict(A=A.round(2).tolist(),T=T.round(2).tolist(),mean=T.mean(1).round(2).tolist())))
