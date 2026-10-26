"""Specificity controls for attention-gated counts + TCGA label audit. Writes results/metrics/gate_controls.csv,
results/metrics/gate_preds.csv, results/metrics/label_audit.csv."""
import sys, numpy as np, pandas as pd, torch
sys.path.insert(0, 'src')
import iel_fm as I, protocols as P
from external import auc_ci
D = P.Design(); pid = {p: i for i, p in enumerate(D.pats.patient_id)}
T = I.merge(I.tcga()); cl = pd.read_csv('data/cptac/labels.csv', index_col=0)['label']
C = {k: v for k, v in I.merge(I.cptac()).items() if k in cl.index}
sl = pd.read_csv('data/surgen/subset.csv', dtype={'case_id': str}).set_index('case_id')
S = {k: v for k, v in I.merge(I.surgen()).items() if k in sl.index}
tr = sorted(k for k in T if pid[k] in set(D.train)); te = sorted(k for k in T if pid[k] in set(D.test))
ytr = np.array([D.y[pid[k]] for k in tr])
rng = np.random.default_rng(7)
real = [I.train_abmil([T[k][1] for k in tr], ytr, s) for s in (0, 1, 2)]
perm = [I.train_abmil([T[k][1] for k in tr], rng.permutation(ytr), s) for s in (0, 1, 2)]
co = {'tcga_test': ({k: T[k] for k in te}, {k: D.y[pid[k]] for k in te}), 'cptac': (C, {k: int(cl[k]) for k in C}),
      'surgen': (S, {k: int(sl.loc[k, 'mmr_loss_binary']) for k in S})}
F = I.CF
rows, preds = [], []
for name, (bags, lab) in co.items():
    for k, (c, e) in bags.items():
        X = torch.from_numpy(e.astype(np.float32)).cuda()
        out = {}
        for tag, ms in (('real', real), ('perm', perm)):
            r = [I.abmil_attention(m, (X - mu) / sd, 'cuda') for m, mu, sd in ms]
            out[tag] = np.mean([a for _, a in r], 0)
            if tag == 'real':
                preds.append(dict(cohort=name, patient=k, y=lab[k], p=float(np.mean([p for p, _ in r]))))
        n = len(c)
        gates = {'all': np.arange(n)}
        for f in (1, 5, 10, 20):
            m_ = max(5, int(round(n * f / 100)))
            gates[f'att{f}'] = np.argsort(-out['real'])[:m_]
            gates[f'perm{f}'] = np.argsort(-out['perm'])[:m_]
            gates[f'epi{f}'] = np.argsort(-c[:, 1])[:m_]  # label-free: most epithelial tiles
        r = dict(cohort=name, patient=k, y=lab[k])
        for g, idx in gates.items():
            for j, f in enumerate(F):
                r[f'{g}:{f}'] = c[idx, j].mean()
        rows.append(r)
d = pd.DataFrame(rows); d.to_csv('results/metrics/gate_patient_counts.csv', index=False)
pd.DataFrame(preds).to_csv('results/metrics/gate_preds.csv', index=False)
out = []
for name, g in d.groupby('cohort', sort=False):
    y = g.y.values
    for col in d.columns[3:]:
        a, lo, hi, _ = auc_ci(y, g[col].values)
        gate, f = col.split(':')
        out.append(dict(cohort=name, gate=gate, feature=f, auc=a, lo=lo, hi=hi))
o = pd.DataFrame(out); o.to_csv('results/metrics/gate_controls.csv', index=False)
pv = o.pivot_table(index=['feature', 'gate'], columns='cohort', values='auc')[['tcga_test', 'cptac', 'surgen']]
print(pv.round(3).to_string())
pr = pd.DataFrame(preds)
for name, g in pr.groupby('cohort', sort=False):
    print(name, 'ABMIL AUC', round(auc_ci(g.y.values, g.p.values)[0], 3))
