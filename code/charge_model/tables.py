"""Rows of SI Tables S17 (16 cells, both realisations) and S18 (charge models) as LaTeX (run from repository root)."""
import pandas as pd, numpy as np
from scipy import stats
o=pd.read_csv('data/raw_results.csv')
a=pd.read_csv('data/realization2/raw_results_realization2.csv')
b=pd.read_csv('data/realization2/results_realization2.csv')[['polymer','pfas','seed','coulSR_kJmol','ljSR_kJmol','totalIE_kJmol']]
r2=pd.concat([a,b])
for df in (o,r2): df.columns=['polymer','pfas','seed','coul','lj','ie']
P=['PTFE','PVDF','ECTFE','PE']; C=['PFOA','PFBS','GenX','FTOH62']
CN={'PFOA':'PFOA','PFBS':'PFBS','GenX':'GenX','FTOH62':'6:2 FTOH'}
f=lambda v:f"${v:+.1f}$" if round(v,1)>0 else (f"$0.0$" if round(v,1)==0 else f"${v:.1f}$")
g=lambda v:f"${v:.1f}$"
top=[];bot=[]
for c in C:
  for p in P:
    x=o[(o.polymer==p)&(o.pfas==c)]; y=r2[(r2.polymer==p)&(r2.pfas==c)].sort_values('seed')
    pv=stats.ttest_ind(y.ie,x.ie,equal_var=False).pvalue
    top.append(f"{p} & {CN[c]} & {g(x.ie.mean())} ({x.ie.std():.1f}) & {g(y.ie.mean())} ({y.ie.std():.1f}) & {pv:.2f} & {f(x.coul.mean())} / {f(y.coul.mean())} & {g(x.lj.mean())} / {g(y.lj.mean())} \\\\")
    bot.append(f"{p} & {CN[c]} & "+" & ".join(g(v) for v in y.ie)+" \\\\")
open('tab_r2_top.tex','w').write("\n".join(top)); open('tab_r2_bot.tex','w').write("\n".join(bot))
# charge table
d=pd.read_csv('data/charge_model/charge_sensitivity_3models.csv')
M=['am1bcc','gasteiger','mmff94']
rows=[];mx=0
for c in C:
  for p in P:
    s=d[(d.polymer==p)&(d.pfas==c)]
    if s.empty: rows.append(f"{p} & {CN[c]} & \\multicolumn{{6}}{{c}}{{trajectories not retained}} \\\\"); continue
    cm=[s[s.model==m].coul.mean() for m in M]; cs=[s[s.model==m].coul.sem() for m in M]; im=[s[s.model==m].ie.mean() for m in M]
    mx=max(mx,abs(cm[1]-cm[0]))
    rows.append(f"{p} & {CN[c]} & "+" & ".join(f"{f(m)} ({e:.1f})" for m,e in zip(cm,cs))+" & "+" & ".join(g(v) for v in im)+" \\\\")
open('tab_chg.tex','w').write("\n".join(rows)); print('max gast-am1 cell diff',mx)
print(open('tab_r2_top.tex').read()); print(open('tab_chg.tex').read())
