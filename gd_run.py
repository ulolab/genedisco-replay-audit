"""Run a GeneDisco acquisition fn under BDA protocol: y = BDA ground-truth Score, test_ratio=0,
batch 128, 5 rounds (round 1 = GeneDisco's random initial batch, as in GeneDisco's loop).
Model/acq code is GeneDisco's, unmodified, except: degenerate-score guard on model predictions."""
import sys, os, pickle, time, json, argparse, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from common import *
import slingpy as sp
from slingpy.data_access.data_sources.hdf5_tools import HDF5Tools
from slingpy import HDF5DataSource
from genedisco.apps import single_cycle_application as sca
from genedisco.apps.active_learning_loop import ActiveLearningLoop
from genedisco.models import meta_models

ap = argparse.ArgumentParser(); ap.add_argument("--screen", default="IFNG"); ap.add_argument("--acq"); ap.add_argument("--seed", type=int)
ap.add_argument("--out"); a = ap.parse_args()
CACHE = f"{ROOT}/cache"
g, hits, ess = load_screen(a.screen)
h5 = f"{CACHE}/bda_{a.screen}.h5"
if not os.path.exists(h5):
    HDF5Tools.save_h5_file(h5, g[["Score"]].values.astype(np.float32), f"bda_{a.screen}",
                           column_names=["Score"], row_names=g.Gene.tolist())
sca.SingleCycleApplication.get_dataset_y = staticmethod(lambda n, c: HDF5DataSource(h5, duplicate_merge_strategy=sp.MeanMergeStrategy()))

# GeneDisco's own test-set evaluation is irrelevant to picking and crashes on an empty test set -> skip.
sca.SingleCycleApplication.evaluate_model = lambda self, *a, **k: {}
ActiveLearningLoop.prepare_hitratio_evaluation = lambda self, *a, **k: None
_orig = meta_models.PytorchMLPRegressorWithUncertainty.predict
def guarded(self, *args, **kw):
    r = _orig(self, *args, **kw)
    if kw.get("return_std_and_margin"):
        check_not_degenerate(np.ravel(r[0]), f"{a.acq}:pred_mean"); check_not_degenerate(np.ravel(r[1]), f"{a.acq}:pred_std")
    return r
meta_models.PytorchMLPRegressorWithUncertainty.predict = guarded

wd = f"{ROOT}/runs/gd_work/{a.screen}/{a.acq}/seed{a.seed}"; os.makedirs(wd, exist_ok=True)
t0 = time.time()
loop = ActiveLearningLoop(model_name="bayesian_mlp", acquisition_function_name=a.acq, acquisition_batch_size=BATCH,
    num_active_learning_cycles=ROUNDS, feature_set_name="achilles", dataset_name="schmidt_2021_ifng",
    cache_directory=CACHE, output_directory=wd, test_ratio=0.0, seed=a.seed)
loop.run()
rp = [pickle.load(open(f"{wd}/cycle_{k}/selected_indices.pickle", "rb")) for k in range(ROUNDS)]
sizes = [len(set(x)) for x in rp]
under = sizes != [BATCH*(k+1) for k in range(ROUNDS)]
if under: print("WARNING UNDER-BUDGET (GeneDisco acq returned duplicate genes; cumulative set deduped):", sizes)
save_run(a.out, f"gd_{a.acq}", a.seed, rp, {"screen": a.screen, "secs": time.time()-t0, "pool": "achilles∩screen", "unique_per_round": sizes, "under_budget": under})
print("done", a.acq, a.seed, hit_ratio(rp[-1], hits, ess), round(time.time()-t0))
