# genedisco-replay-audit

An independent replay of the BioDiscoveryAgent benchmark on GeneDisco, by the Ulo Lab team.

- BioDiscoveryAgent commit: `19673c7`
- GeneDisco commit: `f2ee622`

## What's here
- `ESSENTIALITY_ODDITY.md`: what we observed about the essential-gene filter used in the N/E metric.
- Analysis scripts (`*.py`), the final results (`final_IL2.txt`, `final_IFNG.txt`, `ess_oddity_output.txt`).
- Per-seed outputs: `runs/LOCKED_IL2/` and `runs/LOCKED_IFNG/` (one JSON per method and seed, 16 methods x 10 seeds per screen).
- `investigation/tcell_in_achilles.txt`: the DepMap IDs of the T-cell lines in the Achilles feature (used by `tcell_and_band.py`).
- `IL2_essentiality_66_NE_hits.tsv`: the 66 IL-2 N/E hits from the essentiality ranking, with screen scores and DepMap effects.

## Setup
Run everything from the repo root with Python 3.9. Scripts find their inputs relative to the repo root (set `ULO_ROOT` to use another directory).

1. Clone the two upstream repos into the repo root at the pinned commits:
   ```
   git clone https://github.com/snap-stanford/BioDiscoveryAgent && git -C BioDiscoveryAgent checkout 19673c7
   git clone https://github.com/genedisco/genedisco && git -C genedisco checkout f2ee622
   ```
   `CEGv2.txt` and the IL-2/IFNG ground truth and hit lists (`datasets/`) come from BioDiscoveryAgent.
2. Install the pinned requirements (this also installs `./genedisco` in editable mode):
   ```
   python3.9 -m venv .venv && source .venv/bin/activate
   pip install -r requirements_freeze.txt --extra-index-url https://download.pytorch.org/whl/cpu
   ```
   The extra index is needed for `torch==2.8.0+cpu`.
3. Download the pathway priors used by `ulo_run.py` into `priors/` (Enrichr gene-set libraries; not shipped):
   ```
   mkdir -p priors
   for l in Reactome_2016 KEGG_2019_Human; do
     curl -fL "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=$l" -o priors/$l.gmt
   done
   ```
   Expected SHA-256: `Reactome_2016.gmt` 87b6cb14031230615d43bb377f7943db9dab4773c31811a657b6d73d1dcfb072, `KEGG_2019_Human.gmt` 331210d19babbfad3e0f66d405a3445279d3a42107d9c6babb5fb25831bab963.
4. Fill `cache/` (see Data below).

## Data: `cache/`
No DepMap or HGNC files are shipped; `cache/` is in `.gitignore`. Files the scripts read from `cache/`:
- `achilles.csv`: DepMap Achilles CRISPR gene effect, release **DepMap 20Q4 Public** (the release pinned in GeneDisco's `genedisco/datasets/features/achilles.py`: figshare doi:10.6084/m9.figshare.13237076.v1, file https://ndownloader.figshare.com/files/25494359, CC BY 4.0). Read directly by `drop_jurkat.py` and `tcell_and_band.py`. Official download page: https://depmap.org/portal/download/
- `hgnc_mapping.tsv`: HGNC protein-coding gene table, which GeneDisco downloads from http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt. Used by GeneDisco to update outdated gene symbols and read by `ceg_symbols.py`.
- `achilles.h5`: built by GeneDisco from the two files above; used by `ulo_run.py` and `ess_oddity.py`.
- `schmidt_2021_ifng.h5`, `bda_IL2.h5`, `bda_IFNG.h5`: written on the first `gd_run.py` run (only needed for the GeneDisco baselines).

GeneDisco fills the cache on its first run: `Achilles.load_data("cache")` downloads `achilles.csv` and `hgnc_mapping.tsv` if missing and then builds `achilles.h5`. To fill it without running a baseline:
```
mkdir -p cache
python -c "from genedisco.datasets.features.achilles import Achilles; Achilles.load_data('cache')"
```
The HGNC table is a live file, so a fresh download can differ from ours, which can change a few gene symbols. Our copies had SHA-256 `achilles.csv` 71f32e3515dfa24d3bff559ea99fd5a0fa102c264064a6377e13d8eecd54298f and `hgnc_mapping.tsv` 9a3136b623bdcaaa996e0d8b7dc309a8c26d6a52c135a91b0bacb20d4e770885.

## Reproduce
From the shipped per-seed JSONs (fast):
```
python analyze_final.py runs/LOCKED_IL2 IL2 > final_IL2.txt
python analyze_final.py runs/LOCKED_IFNG IFNG > final_IFNG.txt
python ess_oddity.py > ess_oddity_output.txt
python drop_jurkat.py
python tcell_and_band.py
python ceg_symbols.py
```
`analyze_final.py RUNDIR SCREEN [NSEEDS]` takes the run directory and the screen (`IL2` or `IFNG`); NSEEDS defaults to 10 (seeds 0-9).

Rerun the fast (no-LLM, no-GeneDisco-model) methods, about 1 minute per screen on a CPU, into a separate directory and compare with `runs/LOCKED_*`:
```
for SC in IL2 IFNG; do
  python ulo_run.py --screen $SC \
    --methods random,prior_essentiality,prior_pathway,ulo_pgp,ulo_pgp_noprior,ulo_pgp_nofeedback,ulo_bandit,ulo_bandit_nopgp \
    --seeds 0,1,2,3,4,5,6,7,8,9 --out runs/rerun_$SC
done
```
These runs are deterministic: apart from the `secs` timing field, every output JSON should match the shipped one exactly.

Full locked run, including the slow GeneDisco baselines (`gd_*`, minutes to hours per seed): `bash run_locked.sh`. It skips any (method, seed) JSON that already exists in `runs/LOCKED_*`, so move those aside first to recompute from scratch. One GeneDisco run is `python gd_run.py --screen IL2 --acq coreset --seed 0 --out runs/LOCKED_IL2`.

## Note
The sensitivity numbers here come from our own post-hoc analysis, not from the original authors. We're sharing them to ask a question, and we'd welcome corrections.

## Contact

Questions or corrections: research@ulolab.com
