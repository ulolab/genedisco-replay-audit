"""Post-hoc: (b) drop T-cell lines from Achilles feature; (a) remove near-threshold band of 'non-essential' genes."""
import sys, json, glob, numpy as np, pandas as pd
sys.path.insert(0, "/workspace/ulo-replay/replay")
from common import load_screen, hit_ratio, ROOT
A = pd.read_csv(f"{ROOT}/cache/achilles.csv", index_col=0); A.columns=[c.split(" (")[0] for c in A.columns]
T = [l.strip() for l in open(f"{ROOT}/investigation/tcell_in_achilles.txt")]
for scr in ["IL2","IFNG"]:
    g, hits, ess = load_screen(scr); pool = sorted(set(g.Gene)&set(A.columns))
    for lab, M in [("all 808 lines", A), (f"drop {len(T)} T-cell lines", A.drop(index=T)), ("T-cell lines ONLY", A.loc[T])]:
        s = -M[pool].mean(0).fillna(0); top = list(s.sort_values(ascending=False).index[:640])
        print(scr, f"{lab:24s}", "all=%.3f ne=%.3f" % hit_ratio(top, hits, ess))
    eff = A[pool].mean(0)  # mean gene effect (neg = essential)
    for band in [-1.0, -0.5, -0.3]:
        excl = set(eff[eff < band].index) | ess   # treat these as 'essential' too
        h2 = hits - excl
        row = []
        for m in ["prior_essentiality","ulo_pgp_noprior","gd_coreset","prior_pathway","random"]:
            v=[len((set(json.load(open(f))["rounds"][-1]) - excl) & h2)/len(h2) for f in glob.glob(f"{ROOT}/runs/LOCKED_{scr}/{m}__seed*.json")]
            row.append(f"{m}={np.mean(v):.3f}±{np.std(v):.3f}")
        print(scr, f"band: CEGv2 + effect<{band} (+{len(excl-ess)} genes, N/E hits {len(h2)}):", " ".join(row))
