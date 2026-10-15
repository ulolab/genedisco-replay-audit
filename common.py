"""Shared: BDA-protocol data + scoring. Hit lists/ground truth come unmodified from BioDiscoveryAgent/datasets."""
import numpy as np, pandas as pd, os, json
ROOT = "/workspace/ulo-replay"
BDA = f"{ROOT}/BioDiscoveryAgent"
BATCH, ROUNDS = 128, 5

def load_screen(name):
    g = pd.read_csv(f"{BDA}/datasets/ground_truth_{name}.csv")
    g = g.groupby("Gene", as_index=False)["Score"].mean()
    hits = set(np.load(f"{BDA}/datasets/topmovers_{name}.npy", allow_pickle=True).tolist())
    ess = set(pd.read_csv(f"{BDA}/CEGv2.txt", sep="\t")["GENE"])
    return g, hits, ess

def hit_ratio(picked, hits, ess):
    """Same as BDA analyze.py: all = |pred∩hits|/|hits|; nonessential: remove CEGv2 from both."""
    p = set(picked)
    allr = len(p & hits) / len(hits)
    h2 = hits - ess
    ne = len((p - ess) & h2) / len(h2)
    return allr, ne

def check_not_degenerate(scores, name):
    """Hard stop on constant / NaN / degenerate acquisition scores (lesson #3)."""
    s = np.asarray(scores, dtype=float)
    if s.size == 0 or not np.all(np.isfinite(s)):
        raise RuntimeError(f"DEGENERATE[{name}]: empty or non-finite scores")
    if np.nanstd(s) < 1e-10 or len(np.unique(np.round(s, 12))) < 2:
        raise RuntimeError(f"DEGENERATE[{name}]: constant scores (std={np.nanstd(s):.3g})")
    return s

def save_run(out, method, seed, rounds_picks, meta):
    os.makedirs(out, exist_ok=True)
    json.dump({"method": method, "seed": seed, "rounds": rounds_picks, **meta},
              open(f"{out}/{method}__seed{seed}.json", "w"))
