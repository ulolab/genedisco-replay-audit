"""Prior-ranking baselines + Ulo no-LLM methods under BDA protocol (128x5), candidate pool = Achilles ∩ screen
(identical to the GeneDisco baselines' pool). Feedback = BDA ground-truth Score of picked genes only.
Priors are pinned pre-2021: Reactome_2016 + KEGG_2019_Human (Enrichr libs), DepMap Achilles 20Q4 (Dec 2020).
Task pathways = the ones BDA's own human pathway baseline uses (baseline.py), chosen from the task text."""
import sys, os, time, argparse, numpy as np, pandas as pd, scipy.sparse as ss
sys.path.insert(0, os.path.dirname(__file__))
from common import *
from sklearn.linear_model import Ridge

TASK_PW = {"IFNG": ["R-HSA-877300", "R-HSA-1280215"], "IL2": ["R-HSA-451927", "R-HSA-1280215"]}
HIT_Z = 2.5          # proxy "hit" used for feedback: robust z of observed score (fixed a priori, not tuned)
MAX_PW = 300         # ignore giant pathways when building the graph

def read_gmt(p):
    out = {}
    for line in open(p):
        f = line.rstrip("\n").split("\t"); out[f[0]] = [x for x in f[2:] if x]
    return out

def setup(screen):
    g, hits, ess = load_screen(screen)
    from genedisco.datasets.features.achilles import Achilles
    A = Achilles.load_data(f"{ROOT}/cache")
    names = A.get_row_names(); X = A.get_data()[0]
    idx = {n: i for i, n in enumerate(names)}
    pool = sorted(set(g.Gene) & set(names))
    X = X[[idx[n] for n in pool]]
    y = g.set_index("Gene").loc[pool, "Score"].values
    gi = {n: i for i, n in enumerate(pool)}
    sets = {**read_gmt(f"{ROOT}/priors/Reactome_2016.gmt"), **read_gmt(f"{ROOT}/priors/KEGG_2019_Human.gmt")}
    rows, cols, vals = [], [], []
    for k, gs in sets.items():
        m = sorted({gi[x] for x in gs if x in gi})
        if len(m) < 2 or len(m) > MAX_PW: continue
        w = 1.0 / (len(m) - 1)
        for i in m:
            for j in m:
                if i != j: rows.append(i); cols.append(j); vals.append(w)
    W = ss.csr_matrix((vals, (rows, cols)), shape=(len(pool),) * 2); W.sum_duplicates()
    deg = np.asarray(W.sum(1)).ravel(); Dinv = ss.diags(1 / np.maximum(deg, 1e-12))
    P = (Dinv @ W).T.tocsr()  # column-stochastic-ish transition for propagation
    task = np.zeros(len(pool))
    for k, gs in sets.items():
        if any(k.endswith(pid) for pid in TASK_PW[screen]):
            for x in gs:
                if x in gi: task[gi[x]] += 1
    # Achilles 20Q4 mean gene effect: more negative = more essential
    ess_score = -X.mean(1)
    return dict(pool=pool, X=X, y=y, P=P, deg=deg, task=task, ess_score=ess_score, hits=hits, ess=ess)

def ppr(P, seed_vec, alpha=0.5, iters=30):
    s = seed_vec / max(seed_vec.sum(), 1e-12); r = s.copy()
    for _ in range(iters): r = alpha * s + (1 - alpha) * (P @ r)
    return r

def robust_z(v):
    med = np.median(v); mad = np.median(np.abs(v - med)) * 1.4826 + 1e-9
    return np.abs(v - med) / mad

def top_k(score, taken, k, rng, name):
    s = check_not_degenerate(score, name).copy()
    s = s + rng.uniform(0, 1e-9, len(s))  # seeded tie-break only
    s[list(taken)] = -np.inf
    return list(np.argsort(-s)[:k])

def prior_pathway_score(D):  # task-pathway membership, then graph connectivity to task genes
    return D["task"] * 10 + ppr(D["P"], (D["task"] > 0).astype(float)) * 1e3

def run(method, D, seed):
    rng = np.random.default_rng(seed); n = len(D["pool"]); taken = []
    y, P = D["y"], D["P"]
    rand = lambda: list(rng.choice([i for i in range(n) if i not in set(taken)], BATCH, replace=False))
    picks_by_round = []; arm_log = []
    bandit = {a: [1.0, 1.0] for a in ["pgp", "model", "prior", "random"]}
    for r in range(ROUNDS):
        T = set(taken)
        obs = np.array(taken, dtype=int)
        z = robust_z(y[obs]) if len(obs) else None
        def pgp_score(use_feedback=True, prior_w=0.2):
            seedv = np.zeros(n)
            if use_feedback and len(obs): seedv[obs] = np.maximum(z - 1.0, 0)
            pri = (D["task"] > 0).astype(float)
            s = ppr(P, seedv) if seedv.sum() > 0 else np.zeros(n)
            return s / max(s.max(), 1e-12) + prior_w * ppr(P, pri) / max(ppr(P, pri).max(), 1e-12)
        if method == "random": new = rand()
        elif method == "prior_essentiality": new = top_k(D["ess_score"], T, BATCH, rng, method)
        elif method == "prior_pathway": new = top_k(prior_pathway_score(D), T, BATCH, rng, method)
        elif method.startswith("ulo_pgp"):
            if r == 0:
                new = rand() if method == "ulo_pgp_noprior" else top_k(prior_pathway_score(D), T, BATCH, rng, method)
            elif method == "ulo_pgp_nofeedback":
                new = top_k(pgp_score(use_feedback=False, prior_w=1.0), T, BATCH, rng, method)
            else:
                new = top_k(pgp_score(prior_w=0.0 if method == "ulo_pgp_noprior" else 0.2), T, BATCH, rng, method)
        elif method.startswith("ulo_bandit"):
            if r == 0: new = top_k(prior_pathway_score(D), T, BATCH, rng, method); arm_of = {i: "prior" for i in new}
            else:
                # update arms with last round's outcomes
                for i in last_new:
                    a = arm_of[i]; hit = robust_z(y[obs])[list(obs).index(i)] > HIT_Z
                    bandit[a][0 if hit else 1] += 1
                m = Ridge(alpha=10.0).fit(D["X"][obs], z)
                ranked = {"model": m.predict(D["X"]), "pgp": pgp_score(), "prior": prior_pathway_score(D),
                          "random": rng.uniform(size=n)}
                arms = [a for a in ranked if not (method == "ulo_bandit_nopgp" and a == "pgp")]
                orders = {a: list(np.argsort(-(check_not_degenerate(ranked[a], f"{method}:{a}")))) for a in arms}
                new, arm_of = [], {}
                chosen = set(T)
                while len(new) < BATCH:
                    a = max(arms, key=lambda a: rng.beta(*bandit[a]))
                    while orders[a] and orders[a][0] in chosen: orders[a].pop(0)
                    i = orders[a].pop(0); chosen.add(i); new.append(i); arm_of[i] = a
                arm_log.append({a: sum(1 for i in new if arm_of[i] == a) for a in arms})
            last_new = new
        else: raise ValueError(method)
        assert len(set(new)) == BATCH and not (set(new) & T), "duplicate/repeat picks"
        taken += [int(i) for i in new]
        picks_by_round.append([D["pool"][i] for i in taken])
    return picks_by_round, arm_log

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--screen", default="IFNG"); ap.add_argument("--methods")
    ap.add_argument("--seeds"); ap.add_argument("--out"); a = ap.parse_args()
    D = setup(a.screen)
    for m in a.methods.split(","):
        for s in map(int, a.seeds.split(",")):
            f = f"{a.out}/{m}__seed{s}.json"
            if os.path.exists(f): continue  # resumable
            t0 = time.time(); rp, arms = run(m, D, s)
            save_run(a.out, m, s, rp, {"screen": a.screen, "secs": time.time() - t0, "arms": arms, "pool": "achilles∩screen"})
            print(m, s, hit_ratio(rp[-1], D["hits"], D["ess"]), round(time.time() - t0, 1), flush=True)
