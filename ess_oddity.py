import sys, json, glob, numpy as np, pandas as pd
sys.path.insert(0, "/workspace/ulo-replay/replay")
from common import load_screen, hit_ratio, ROOT
for scr in ["IL2", "IFNG"]:
    g, hits, ess = load_screen(scr)
    sc = g.set_index("Gene")["Score"]
    from genedisco.datasets.features.achilles import Achilles
    A = Achilles.load_data(f"{ROOT}/cache"); names = A.get_row_names(); X = A.get_data()[0]
    idx = {n: i for i, n in enumerate(names)}
    pool = sorted(set(g.Gene) & set(names))
    essc = pd.Series(-X[[idx[n] for n in pool]].mean(1), index=pool)
    ne = hits - ess
    print(f"\n=== {scr}: genes {len(g)}, pool {len(pool)}, hits {len(hits)}, hits∩CEGv2 {len(hits&ess)}, N/E hits {len(ne)}, N/E hits in pool {len(ne & set(pool))}")
    print(" hit score sign: pos", sum(sc[h] > 0 for h in hits if h in sc), "neg", sum(sc[h] < 0 for h in hits if h in sc))
    print(" CEGv2 in pool", len(ess & set(pool)))
    # Achilles effect quantiles: hits vs non-hits
    pe = essc.rank(pct=True)
    print(" mean essentiality pct: N/E hits %.3f, other pool %.3f" % (pe[list(ne & set(pool))].mean(), pe[list(set(pool) - hits)].mean()))
    top = list(essc.sort_values(ascending=False).index[:640])
    a, n = hit_ratio(top, hits, ess)
    print(" ess top640 all %.3f ne %.3f; CEGv2 among top640: %d" % (a, n, len(set(top) & ess)))
    got = [x for x in top if x in ne]
    print(" N/E hits picked (%d):" % len(got), ", ".join(f"{x}({sc[x]:+.2f},eff={-essc[x]:.2f})" for x in got))
    # strong-essential non-CEG: Achilles effect < -0.5
    strong = set(essc[essc > 0.5].index)
    print(" N/E hits with Achilles effect<-0.5: %d / %d; picked ones that are strong: %d" % (len(ne & strong), len(ne & set(pool)), len(set(got) & strong)))
    # post-hoc: N/E excluding all Achilles-common-essential(effect<-0.5)
    for m in ["prior_essentiality", "ulo_pgp_noprior", "gd_coreset", "random", "prior_pathway"]:
        fs = sorted(glob.glob(f"{ROOT}/runs/LOCKED_{scr}/{m}__seed*.json"))
        if not fs: continue
        v = []
        for f in fs:
            p = set(json.load(open(f))["rounds"][-1])
            e2 = ess | strong
            h2 = hits - e2
            v.append((hit_ratio(p, hits, ess)[1], len((p - e2) & h2) / len(h2), len(p & ne & strong) / max(len(p & ne), 1)))
        v = np.array(v); print(f"  {m:20s} n={len(fs)} NE={v[:,0].mean():.3f} NE_strict(post-hoc)={v[:,1].mean():.3f} frac_of_NEhits_that_are_strong_ess={v[:,2].mean():.2f}")
