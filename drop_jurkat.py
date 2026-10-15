import sys, numpy as np, pandas as pd
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from common import load_screen, hit_ratio, ROOT
A = pd.read_csv(f"{ROOT}/cache/achilles.csv", index_col=0); A.columns=[c.split(" (")[0] for c in A.columns]
for scr in ["IL2","IFNG"]:
    g, hits, ess = load_screen(scr); pool = sorted(set(g.Gene)&set(A.columns))
    for lab, M in [("all lines", A), ("drop ACH-000995", A.drop(index=["ACH-000995"]))]:
        s = -M[pool].mean(0).fillna(0); top = list(s.sort_values(ascending=False).index[:640])
        print(scr, lab, "all=%.3f ne=%.3f" % hit_ratio(top, hits, ess))
    j = -A.loc["ACH-000995", pool]; ne=list((hits-ess)&set(pool)); o=list(set(pool)-hits)
    print(scr, " ACH-000995 effect: N/E hits mean %.3f vs non-hits %.3f" % (-j[ne].mean(), -j[o].mean()))
