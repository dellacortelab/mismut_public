import sys, numpy as np, pandas as pd, torch
sys.path.insert(0,'src'); import iel_fm as I, protocols as P
D=P.Design(); pid={p:i for i,p in enumerate(D.pats.patient_id)}
T=I.merge(I.tcga()); cl=pd.read_csv('data/cptac/labels.csv',index_col=0)['label']; C={k:v for k,v in I.merge(I.cptac()).items() if k in cl.index}
sl=pd.read_csv('data/surgen/subset.csv',dtype={'case_id':str}).set_index('case_id'); S={k:v for k,v in I.merge(I.surgen()).items() if k in sl.index}
tr=sorted(k for k in T if pid[k] in set(D.train)); te=sorted(k for k in T if pid[k] in set(D.test))
ytr=np.array([D.y[pid[k]] for k in tr])
models=[I.train_abmil([T[k][1] for k in tr], ytr, s) for s in (0,1,2)]
co={'tcga_test':({k:T[k] for k in te},{k:D.y[pid[k]] for k in te}),'cptac':(C,{k:int(cl[k]) for k in C}),'surgen':(S,{k:int(sl.loc[k,'mmr_loss_binary']) for k in S})}
F=I.CF; rows=[]
for name,(bags,lab) in co.items():
    for k,(c,e) in bags.items():
        X=torch.from_numpy(e.astype(np.float32)).cuda()
        a=np.mean([I.abmil_attention(m,(X-mu)/sd,'cuda')[1] for m,mu,sd in models],0)
        top=np.argsort(-a)[:max(5,len(a)//20)]   # top 5% attention tiles
        r=dict(cohort=name,patient=k,y=lab[k])
        for j,f in enumerate(F): r[f'top_{f}']=c[top,j].mean(); r[f'all_{f}']=c[:,j].mean()
        rows.append(r)
d=pd.DataFrame(rows)
out=[]
for name,g in d.groupby('cohort',sort=False):
    for y,gg in g.groupby('y'):
        out.append({'cohort':name,'class':'MSI' if y else 'MSS',**{f: round(gg[f'top_{f}'].mean()/max(gg[f'all_{f}'].mean(),1e-6),2) for f in F}})
print('enrichment = mean count in top-5%-attention tiles / mean over all tiles')
print(pd.DataFrame(out).to_string(index=False))
d.to_csv('results/metrics/attention_enrichment.csv',index=False)
