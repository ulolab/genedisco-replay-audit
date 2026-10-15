#!/bin/bash
# LOCKED confirmatory run. Resumable: finished (method,seed) JSONs are skipped. Max 2 concurrent GeneDisco jobs.
cd /workspace/ulo-replay && source .venv/bin/activate; mkdir -p logs/locked
SEEDS=0,1,2,3,4,5,6,7,8,9
for SC in IL2 IFNG; do
  python replay/ulo_run.py --screen $SC --methods random,prior_essentiality,prior_pathway,ulo_pgp,ulo_pgp_noprior,ulo_pgp_nofeedback,ulo_bandit,ulo_bandit_nopgp --seeds $SEEDS --out runs/LOCKED_$SC >> logs/locked/ulo_$SC.log 2>&1
  echo "ulo $SC exit=$?" >> logs/locked/status.log
done
for SC in IL2 IFNG; do for acq in random topuncertain softuncertain marginsample coreset badge kmeans_embedding kmeans_data; do for s in 0 1 2 3 4 5 6 7 8 9; do
  [ -f runs/LOCKED_$SC/gd_${acq}__seed$s.json ] || echo "$SC $acq $s"
done; done; done | xargs -P 2 -n 3 sh -c 'rm -rf runs/gd_work/$0/$1/seed$2; python replay/gd_run.py --screen $0 --acq $1 --seed $2 --out runs/LOCKED_$0 > logs/locked/gd_$0_$1_$2.log 2>&1; echo "$0 $1 $2 exit=$?" >> logs/locked/status.log'
echo LOCKEDDONE >> logs/locked/status.log
