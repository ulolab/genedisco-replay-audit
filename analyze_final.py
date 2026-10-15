"""Applies LOCKED_PLAN decision rules exactly. Usage: analyze_final.py RUNDIR SCREEN"""
import sys, os, json, glob, numpy as np, collections
sys.path.insert(0, os.path.dirname(__file__)); from common import *
d, screen = sys.argv[1], sys.argv[2]
g, hits, ess = load_screen(screen)
PUB = {"IFNG": {"Sonnet 3.5": (0.095, 0.107), "DiscoBAX": (0.088, 0.069)},
       "IL2": {"Sonnet 3.5": (0.104, 0.122), "DiscoBAX": (0.074, 0.057)}}  # (All, N/E), BDA Table 1
ULO = ["ulo_bandit", "ulo_pgp_noprior"]
R = collections.defaultdict(dict)
for f in glob.glob(f"{d}/*.json"):
    j = json.load(open(f)); R[j["method"]][j["seed"]] = j
seeds = list(range(int(sys.argv[3]) if len(sys.argv) > 3 else 10))
def hr(m, s, k): return hit_ratio(R[m][s]["rounds"][-1], hits, ess)[k]
def jac(a, b): a, b = set(a), set(b); return len(a & b) / len(a | b)
complete = {m for m in R if all(s in R[m] for s in seeds)}
missing = {m: [s for s in seeds if s not in R[m]] for m in R if m not in complete}
under = {m for m in R if any(R[m][s].get("under_budget") for s in R[m])}
out = []
P = out.append
P(f"Seeds required: 0-9. Incomplete methods: {missing or 'none'}. Under-budget (excluded from best-baseline): {sorted(under) or 'none'}\n")
P("| method | N/E mean ± SD | All mean ± SD | Jaccard vs prior_pathway | Jaccard vs prior_essentiality | mean sec/run |")
P("|---|---|---|---|---|---|")
stats = {}
for m in sorted(complete, key=lambda m: -np.mean([hr(m, s, 1) for s in seeds])):
    ne = np.array([hr(m, s, 1) for s in seeds]); al = np.array([hr(m, s, 0) for s in seeds])
    jp = np.mean([jac(R[m][s]["rounds"][-1], R["prior_pathway"][s]["rounds"][-1]) for s in seeds])
    je = np.mean([jac(R[m][s]["rounds"][-1], R["prior_essentiality"][s]["rounds"][-1]) for s in seeds])
    stats[m] = (ne, al, jp, je)
    P(f"| {m}{' (UNDER BUDGET)' if m in under else ''} | {ne.mean():.3f} ± {ne.std(ddof=1):.3f} | {al.mean():.3f} ± {al.std(ddof=1):.3f} | {jp:.2f} | {je:.2f} | {np.mean([R[m][s]['secs'] for s in seeds]):.0f} |")
ULO_SET = set(ULO) | {"ulo_pgp", "ulo_pgp_nofeedback", "ulo_bandit_nopgp"}
base = [m for m in complete if m not in ULO_SET and m not in under]
best = max([m for m in base if m.startswith("gd_")], key=lambda m: stats[m][0].mean())
rng = np.random.default_rng(12345)
def boot(x):
    bs = rng.choice(x, (10000, len(x)), replace=True).mean(1); return np.percentile(bs, [2.5, 97.5])
def verdict(lo, hi): return "WIN" if lo > 0 else ("LOSS" if hi < 0 else "TIE")
P(f"\nBest re-run GeneDisco baseline on N/E (not under budget): {best}\n")
P("| Ulo method | comparator | metric | mean diff | paired bootstrap 95% CI | verdict |")
P("|---|---|---|---|---|---|")
claims = {}
for u in ULO:
    ok = True
    for c in [best, "prior_essentiality", "prior_pathway"]:
        for k, name in [(0, "N/E"), (1, "All")]:
            diff = stats[u][k] - stats[c][k]  # paired by seed
            lo, hi = boot(diff); v = verdict(lo, hi)
            P(f"| {u} | {c} | {name} | {diff.mean():+.3f} | [{lo:+.3f}, {hi:+.3f}] | {v} |")
            if name == "N/E" and v != "WIN": ok = False
    collapsed = max(stats[u][2], stats[u][3]) >= 0.8
    claims[u] = ok and not collapsed
P("\n| Ulo method | published comparator (copied numbers) | metric | ours mean ± 2SE | published | label |")
P("|---|---|---|---|---|---|")
for u in ULO:
    for pm, (pa, pn) in PUB[screen].items():
        for k, name, pv in [(0, "N/E", pn), (1, "All", pa)]:
            x = stats[u][k]; se = x.std(ddof=1) / np.sqrt(len(x))
            lab = "above" if x.mean() - 2 * se > pv else ("below" if x.mean() + 2 * se < pv else "comparable")
            P(f"| {u} | {pm} | {name} | {x.mean():.3f} ± {2*se:.3f} | {pv:.3f} | {lab} (copied) |")
P("\nUlo 'win' claim (beats best re-run baseline AND both prior rankings on N/E, AND Jaccard<0.8 with simplest baselines):")
for u, v in claims.items(): P(f"- {u}: {'WIN CLAIM SUPPORTED' if v else 'NO WIN CLAIM'}")
print("\n".join(out))
