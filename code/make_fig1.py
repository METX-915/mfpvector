import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle, FancyArrowPatch
import numpy as np

plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['DejaVu Sans'],
                     'mathtext.fontset':'dejavusans'})

BLUE_D='#2E74B5'; BLUE_L='#BDD7EE'; ORANGE='#ED7D31'
GREEN_E='#548235'; GREEN_L='#E2EFDA'; ORG_L='#FFF2CC'; RED_E='#C00000'; RED_L='#FBE5E5'
BLUE_BOX='#DEEBF7'; GREY='#595959'

fig=plt.figure(figsize=(9.0,4.96)); fig.patch.set_facecolor('white')
ax=fig.add_axes([0,0,1,1]); ax.set_xlim(0,100); ax.set_ylim(0,55); ax.axis('off')

# ── panel titles ─────────────────────────────────────────────
ax.text(2.0,52.6,'(a) Fluorine accounting in environmental samples',
        fontsize=10.5,fontweight='bold',va='center',ha='left')
ax.text(51.5,52.6,'(b) Microfluoroplastics (MFPs) as a PFAS vector',
        fontsize=10.5,fontweight='bold',va='center',ha='left')
ax.plot([49.5,49.5],[3.5,51.0],color='#D9D9D9',lw=1.0)

# ── (a) TOF bar ──────────────────────────────────────────────
x0,y0,w,h=3.5,37.0,42.0,6.6
ax.text(x0+w/2,y0+h+1.5,'Total organic fluorine (TOF)',fontsize=11,
        fontweight='bold',ha='center',va='bottom')
eof_w=w*0.20
ax.add_patch(Rectangle((x0,y0),eof_w,h,facecolor=BLUE_D,edgecolor='black',lw=1.2,zorder=3))
ax.add_patch(Rectangle((x0+eof_w,y0),w-eof_w,h,facecolor=ORANGE,edgecolor='black',lw=1.2,zorder=3))
ax.text(x0+eof_w/2,y0+h/2+0.80,'EOF',fontsize=10,fontweight='bold',color='white',ha='center',va='center',zorder=4)
ax.text(x0+eof_w/2,y0+h/2-1.45,'(extracted)',fontsize=7.0,color='white',ha='center',va='center',zorder=4)
ax.text(x0+eof_w+(w-eof_w)/2,y0+h/2,'Non-extractable organofluorine',
        fontsize=10.8,fontweight='bold',color='white',ha='center',va='center',zorder=4)

ax.add_patch(FancyArrowPatch((x0+w*0.44,y0-6.6),(x0+w*0.40,y0-0.6),
        arrowstyle='-|>',mutation_scale=13,lw=1.1,color='black'))
ax.text(x0+w*0.44,y0-8.4,'microfluoroplastics expected here (insoluble,',
        fontsize=10,style='italic',ha='center',va='top')
ax.text(x0+w*0.44,y0-11.0,'not recovered by routine extraction);',fontsize=10,style='italic',ha='center',va='top')
ax.text(x0+w*0.44,y0-13.6,'partly outside TOF if combustion is incomplete',fontsize=10,style='italic',ha='center',va='top')

# (나) 적용: of total fluorine → of TF/TOF
ax.text(x0+w/2,15.0,'EOF = 0.04–15% of TF or TOF',fontsize=11.5,ha='center',va='center')
ax.text(x0+w/2,12.4,'(denominator differs between studies)',fontsize=9.5,color=GREY,ha='center',va='center')

# ── (b) particle ─────────────────────────────────────────────
cx,cy,r=63.0,42.0,6.4
ax.add_patch(Circle((cx,cy),r,facecolor=BLUE_L,edgecolor=BLUE_D,lw=1.2,alpha=0.95,zorder=2))
rng=np.random.default_rng(7)
n=34; ang=rng.uniform(0,2*np.pi,n); rad=r*np.sqrt(rng.uniform(0.02,0.82,n))
ax.scatter(cx+rad*np.cos(ang),cy+rad*np.sin(ang),s=15,color=BLUE_D,zorder=3)

for (px,py,lx,ly,lab) in [(cx+4.3,cy+3.4,78.5,45.6,'fluorine reservoir'),
                          (cx+5.9,cy+0.2,78.5,42.2,'PFAS carrier'),
                          (cx+4.3,cy-3.4,78.5,38.8,'leachate source')]:
    ax.plot([px,lx-0.8],[py,ly],color='black',lw=0.9,zorder=4)
    ax.scatter([px],[py],s=22,color='black',zorder=5)
    ax.text(lx,ly,lab,fontsize=10,va='center',ha='left')
ax.text(cx,cy-r-2.2,'Microfluoroplastic (PTFE, PVDF, ECTFE)',
        fontsize=10.5,fontweight='bold',ha='center',va='top')

# ── (b) 4-step pipeline ──────────────────────────────────────
steps=[('Adsorb',BLUE_BOX,BLUE_D),('Transport',GREEN_L,GREEN_E),
       ('Release\n(gut)',ORG_L,ORANGE),('Receptor\n(PPAR$\\alpha$)',RED_L,RED_E)]
bw,bh,gap,by=9.4,8.2,2.5,20.0
bx=52.0
for i,(lab,fc,ec) in enumerate(steps):
    x=bx+i*(bw+gap)
    ax.add_patch(FancyBboxPatch((x,by),bw,bh,boxstyle="round,pad=0.25,rounding_size=1.1",
                                facecolor=fc,edgecolor=ec,lw=1.4,zorder=3))
    ax.text(x+bw/2,by+2.0,lab,fontsize=10.5,fontweight='bold',ha='center',va='center',zorder=4)
    if i<3:
        ax.add_patch(FancyArrowPatch((x+bw+0.45,by+bh/2),(x+bw+gap-0.45,by+bh/2),
                     arrowstyle='-|>',mutation_scale=13,lw=1.3,color='black',zorder=4))
# step icons
ax.scatter([54.3,56.0,57.7,55.1,56.9],[26.6,27.2,26.6,25.6,25.6],s=13,color='#203864',zorder=4)
xs=np.linspace(64.4,70.6,60)
for off in (26.9,26.0,25.1):
    ax.plot(xs,off+0.42*np.sin((xs-64.4)*2.2),color=GREEN_E,lw=1.3,zorder=4)
ax.scatter([76.6,78.0,77.2,79.4],[25.4,26.2,26.9,25.8],s=11,color=ORANGE,zorder=4)
ax.scatter([81.6],[27.4],s=34,color=ORANGE,zorder=4)
ax.add_patch(FancyArrowPatch((78.9,25.9),(81.0,27.1),arrowstyle='-|>',
             mutation_scale=10,lw=1.1,color=ORANGE,zorder=4))
ax.add_patch(FancyBboxPatch((89.6,25.2),3.9,2.6,boxstyle="round,pad=0.1,rounding_size=0.5",
             facecolor='none',edgecolor=RED_E,lw=1.6,zorder=4))
ax.add_patch(Rectangle((91.0,25.9),1.1,1.2,facecolor='black',zorder=5))

# ── bottom formulas (네 항 모두 hat 표기) ─────────────────────
ax.text(50.0,8.6,r'Vector hazard $=\ \hat{A}_{\mathrm{ads}}\ \times\ \hat{D}_{\mathrm{desorp}}'
                 r'\ \times\ \hat{P}_{\mathrm{persist}}\ \times\ \hat{H}_{\mathrm{receptor}}$',
        fontsize=13.5,fontweight='bold',ha='center',va='center')
ax.text(50.0,4.4,'Environmental risk = particle abundance × vector hazard  (conceptual; not evaluated here)',
        fontsize=10.5,style='italic',color=GREY,ha='center',va='center')

fig.savefig('Fig1.pdf',dpi=600,facecolor='white',bbox_inches=None)
fig.savefig('Fig1.png',dpi=300,facecolor='white',bbox_inches=None)
print("OK")
