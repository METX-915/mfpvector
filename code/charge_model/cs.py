"""SI Section S12: merges AM1-BCC (screen) with Gasteiger-Marsili / MMFF94 re-evaluations; writes data/charge_model/charge_sensitivity_3models.csv."""
import pandas as pd, numpy as np
from scipy import stats
cs=pd.read_csv('data/charge_model/charge_sensitivity.csv')[['polymer','pfas','seed','model','coulSR_kJmol','ljSR_kJmol']]
ab=pd.read_csv('data/realization2/raw_results_realization2.csv'); ab['model']='am1bcc'
ab=ab[['polymer','pfas','seed','model','coulSR_kJmol','ljSR_kJmol']]
d=pd.concat([cs,ab]); d.columns=['polymer','pfas','seed','model','coul','lj']; d['ie']=d.coul+d.lj
print(d.groupby('model').size())
lj=d.pivot_table(index=['polymer','pfas','seed'],columns='model',values='lj'); print('LJ maxdiff',(lj.max(axis=1)-lj.min(axis=1)).max())
P=['PTFE','PVDF','ECTFE','PE']; M=['am1bcc','gasteiger','mmff94']
cm=d.groupby(['polymer','pfas','model']).coul.mean().unstack()[M]
cs_=d.groupby(['polymer','pfas','model']).coul.sem().unstack()[M]
im=d.groupby(['polymer','pfas','model']).ie.mean().unstack()[M]
print(pd.concat([cm.round(1),cs_.round(1)],axis=1,keys=['mean','sem']))
print(im.round(1))
print('Spearman cell-mean coul'); print(cm.corr('spearman').round(2))
for m in M:
  s=d[d.model==m]; print(m,'positive frac by polymer',s.groupby('polymer').coul.apply(lambda x:f"{(x>0).sum()}/{len(x)}").to_dict(), 'PTFE range',s[s.polymer=='PTFE'].coul.min().round(1),s[s.polymer=='PTFE'].coul.max().round(1))
  pm=s.groupby('polymer')[['coul','ie']].mean().loc[P]; print(pm.round(1).T)
  # within compound: PTFE vs each other polymer coul (cells present)
  for c in ['PFBS','FTOH62','GenX']:
    x=s[(s.polymer=='PTFE')&(s.pfas==c)].coul
    out=[]
    for q in ['PVDF','ECTFE','PE']:
      y=s[(s.polymer==q)&(s.pfas==c)].coul
      if len(y): out.append(f"{q} d={x.mean()-y.mean():.1f} p={stats.ttest_ind(x,y,equal_var=False).pvalue:.1g}")
    print('  ',c,out)
  # compound-block (PFBS,FTOH62,GenX) polymer IE ranking and PTFE vs rest
  bl=s[s.pfas.isin(['PFBS','FTOH62','GenX'])].groupby(['pfas','polymer']).ie.mean().unstack()[P]
  rest=bl[['PVDF','ECTFE','PE']].mean(axis=1)
  print('  IE block means',bl.mean().round(1).to_dict(),'PTFE-rest',(bl.PTFE-rest).round(1).tolist())
d.to_csv('data/charge_model/charge_sensitivity_3models.csv',index=False)
