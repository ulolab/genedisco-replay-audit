import sys, numpy as np, pandas as pd
sys.path.insert(0, "/workspace/ulo-replay/replay")
from common import load_screen, ROOT
ceg = pd.read_csv(f"{ROOT}/BioDiscoveryAgent/CEGv2.txt", sep="\t")
h = pd.read_csv(f"{ROOT}/cache/hgnc_mapping.tsv", sep="\t", usecols=["hgnc_id","symbol"], dtype=str)
cur = set(h[h.hgnc_id.isin(ceg.HGNC_ID)].symbol)
for scr in ["IL2","IFNG"]:
    g, hits, ess = load_screen(scr)
    genes=set(g.Gene)
    print(scr, "CEGv2 symbols present in screen:", len(ess&genes), "/", len(ess), "; current-symbol CEGv2 in screen:", len(cur&genes))
    miss = sorted((cur - ess) & hits)
    print("  hits that are CEGv2 under renamed symbol (escape N/E filter):", len(miss), miss)
    print("  N/E denominator now", len(hits-ess), "-> with renamed symbols", len(hits-ess-cur))
