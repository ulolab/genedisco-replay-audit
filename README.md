# genedisco-replay-audit

An independent replay of the BioDiscoveryAgent benchmark on GeneDisco, by the Ulo Lab team.

- BioDiscoveryAgent commit: `19673c7`
- GeneDisco commit: `f2ee622`

## What's here
- `ESSENTIALITY_ODDITY.md`: what we observed about the essential-gene filter used in the N/E metric.
- Analysis scripts (`*.py`), the final IL-2 results (`final_IL2.txt`), and per-seed outputs (`per_seed_IL2/`).

## Reproduce
Run everything from the repo root (Python 3.9; `pip install -r requirements_freeze.txt`).
1. Clone BioDiscoveryAgent (`19673c7`) and GeneDisco (`f2ee622`) into the repo root before running anything. `CEGv2.txt` and the IL-2/IFNG hit lists come from BioDiscoveryAgent.
2. `bash run_locked.sh` runs every method and seed with the exact settings used.
3. `python analyze_final.py` produces `final_IL2.txt` and `final_IFNG.txt`; `python ess_oddity.py` produces `ess_oddity_output.txt`.
4. `IL2_essentiality_66_NE_hits.tsv` lists the 66 IL-2 N/E hits from the essentiality ranking, with screen scores and DepMap effects.

## Data
No datasets are included. Get GeneDisco and DepMap from their original sources; DepMap is linked rather than redistributed because of its license.

## Note
The sensitivity numbers here come from our own post-hoc analysis, not from the original authors. We're sharing them to ask a question, and we'd welcome corrections.

Contact: research@ulolab.com
