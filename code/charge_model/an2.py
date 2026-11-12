"""SI Section S11: compound-blocked statistics of both slab realisations (run from repository root)."""
import pandas as pd, numpy as np
from scipy import stats
o=pd.read_csv('data/raw_results.csv')
a=pd.read_csv('data/realization2/raw_results_realization2.csv')
b=pd.read_csv('data/realization2/results_realization2.csv')[['polymer','pfas','seed','coulSR_kJmol','ljSR_kJmol','totalIE_kJmol']]
r2=pd.concat([a,b]); 
for df in (o,r2): df.columns=['polymer','pfas','seed','coul','lj','ie']
P=['PTFE','PVDF','ECTFE','PE']; C=['PFOA','PFBS','GenX','FTOH62']
def block(df,excl=False):
    if excl: df=df[df.lj<=-35]
    cm=df.groupby(['pfas','polymer']).ie.mean().unstack()[P].loc[C]
    rest=cm[['PVDF','ECTFE','PE']].mean(axis=1); dd=cm.PTFE-rest
    t=stats.ttest_rel(cm.PTFE,rest)
    res=dict(n=len(df),means=cm.mean().round(1).to_dict(),delta=round(dd.mean(),1),p=round(t.pvalue,4),dz=round(dd.mean()/dd.std(ddof=1),2),
      fried=round(stats.friedmanchisquare(*[cm[p] for p in P]).pvalue,3))
    for q in ['PVDF','ECTFE','PE']: res['vs'+q]=(round((cm.PTFE-cm[q]).mean(),1),round(stats.ttest_rel(cm.PTFE,cm[q]).pvalue,4))
    for q,r in [('PVDF','PE'),('ECTFE','PE'),('PVDF','ECTFE')]: res[q+'-'+r]=round(stats.ttest_rel(cm[q],cm[r]).pvalue,3)
    return res
for nm,df in [('orig',o),('r2',r2)]:
    print(nm,block(df)); print(nm,'excl',block(df,True))
# per-cell excl
for p in P:
  for c in C:
    x=o[(o.polymer==p)&(o.pfas==c)]; y=r2[(r2.polymer==p)&(r2.pfas==c)]
    xe=x[x.lj<=-35]; ye=y[y.lj<=-35]
    if len(xe)<5 or len(ye)<5:
      print(p,c,'excl n',len(xe),len(ye),'diff',round(ye.ie.mean()-xe.ie.mean(),1),'p',round(stats.ttest_ind(ye.ie,xe.ie,equal_var=False).pvalue,2),'full diff',round(y.ie.mean()-x.ie.mean(),1))
# replicate SDs and seed table
for p in P:
  for c in C:
    y=r2[(r2.polymer==p)&(r2.pfas==c)].sort_values('seed'); print(p,c,[round(v,1) for v in y.ie], round(y.ie.std(),1), round(o[(o.polymer==p)&(o.pfas==c)].ie.std(),1))
# coul polymer level per realisation, replicate-level ranges
for nm,df in [('orig',o),('r2',r2)]:
  print(nm, df.groupby('polymer').coul.agg(['min','max']).round(1).to_dict())
# Spearman cell-mean across realisations
m1=o.groupby(['polymer','pfas']).ie.mean(); m2=r2.groupby(['polymer','pfas']).ie.mean()
print('spearman cells',stats.spearmanr(m1,m2), 'pearson',stats.pearsonr(m1,m2))
c1=o.groupby(['polymer','pfas']).coul.mean(); c2=r2.groupby(['polymer','pfas']).coul.mean()
print('coul spearman',stats.spearmanr(c1,c2))
# Welch PVDF-vs-PTFE contrast 4 compounds
for nm,df in [('orig',o),('r2',r2)]:
  cm=df.groupby(['pfas','polymer']).ie.mean().unstack()
  d=cm.PVDF-cm.PTFE; print(nm,'PVDF-PTFE',d.round(1).to_dict(),round(d.mean(),1),stats.ttest_rel(cm.PVDF,cm.PTFE).pvalue)
