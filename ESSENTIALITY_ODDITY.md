# Why does the essentiality ranking win IL-2 "non-essential" hits?

Post-hoc investigation, 2026-10-14 (PT). Read-only. Nothing under `replay/`, `runs/` or `logs/` was changed, and the IFNG run (pid 1171565) was left alone.
Scripts: `ess_oddity.py`, `ceg_symbols.py`, `drop_jurkat.py`, `tcell_and_band.py`. T-cell line lookup: `tcell_cellosaurus.tsv` and `tcell_in_achilles.txt`.

## Verdict
**This is not a bug in our adapter, and it is not leakage. The cause is the narrow "non-essential" definition, applied to a screen that is full of essential genes.**
"N/E" in BioDiscoveryAgent (BDA) means "not in CEGv2", and CEGv2 is a small list of 684 core-essential genes. Thousands of genes that are strongly essential in DepMap are not on it: ribosomal proteins (RPL/RPS), proteasome (PSMD), replication (POLE, GINS1, CDC6/7), and spliceosome genes. The IL-2 hit list is unusually full of these genes. So any ranking that sorts by DepMap essentiality collects "non-essential" hits that are really essential genes. Every one of the 66 N/E hits that prior_essentiality picks has a mean DepMap effect below -1.1. Once those genes are treated as essential, its N/E score falls to about 0.
The reported 0.137 is computed correctly under BDA's metric. It does not show that essentiality finds IL-2 biology.

## Evidence
1. **Our definition matches BDA's.** `replay/common.py:hit_ratio` removes CEGv2 (`BioDiscoveryAgent/CEGv2.txt`) from both the picks and the hits. BDA's `analyze.py` does the same with `--essential 0`. Both use the same file, the same set subtraction and the same denominator (|hits − CEGv2|), and the hit lists come from BDA unmodified.
   - Symbol drift (old CEGv2 symbols vs. current ones, e.g. AARS → AARS1) lets 0 IL-2 hits and 1 IFNG hit escape the filter. That is negligible.
2. **How much the IL-2 hits overlap with essential genes.**

   | | IL-2 | IFNG |
   |---|---|---|
   | Hits | 654 | 920 |
   | Hits in CEGv2 | 173 (26%) | 118 (13%) |
   | Mean DepMap essentiality percentile of N/E hits (pool baseline 0.49) | 0.70 | 0.60 |
   | N/E hits with mean effect < -0.5 | 167 of 437 in pool | 161 of 763 |

   IL-2 hits are heavily enriched for essential genes, inside CEGv2 and beyond it. The IL-2 readout is sorting on IL-2 production, and knocking out a housekeeping gene plausibly lowers the signal. That is a biological confound, not IL-2-specific regulation.
3. **What essentiality ranking actually picks.** It picks 66 N/E hits, and all 66 are core machinery: about 20 RPL/RPS genes, PSMD2/6/7, POLE, POLD3, RFC3, GINS1, CDC6/7/23, SNRPB/C, PRPF6, FBL, SMC3, ATP1A1, and so on. Their mean effects range from -1.1 to -2.1. None of them is a T-cell or IL-2 pathway gene. Full list: `ess_oddity.py` output.
4. **Leakage check.** DepMap Achilles 20Q4 does include Jurkat (ACH-000995) and 4 other T-cell lines: PF-382, SUP-T1, SMZ-1 and CTV-1. These were identified through Cellosaurus; the DepMap sample_info download returned 403.
   - The IL-2 and IFNG screens come from primary human T cells (Schmidt et al.), not from DepMap.
   - Removing the 5 T-cell lines does not change the ranking's result at all: IL-2 N/E 0.125 → 0.125, IFNG 0.057 → 0.057.
   - Using only the T-cell lines gives the same answer (IL-2 0.125).
   - So there is no meaningful leakage. The signal is pan-cancer essentiality.
   - Note: recomputing the ranking from `achilles.csv` gives 0.125 instead of the locked 0.137. The difference likely comes from NaN handling in GeneDisco's h5 vs. the raw CSV. It does not affect the comparison.
5. **Threshold-edge test (post-hoc).** I re-scored the locked picks after also treating genes with mean DepMap effect below a cutoff as essential.

   | IL-2 N/E (post-hoc) | locked | +effect<-1.0 | +effect<-0.5 | +effect<-0.3 |
   |---|---|---|---|---|
   | prior_essentiality | 0.137 | 0.003 | 0.000 | 0.000 |
   | ulo_pgp_noprior (random-start graph) | 0.116 | 0.087 ± 0.019 | 0.076 ± 0.045 | 0.079 ± 0.057 |
   | gd_coreset | 0.086 | 0.066 | 0.046 | 0.039 |
   | prior_pathway | 0.075 | 0.077 | 0.082 | 0.090 |
   | random | 0.035 | 0.033 | 0.032 | 0.032 |

   On IFNG, the -0.5 cutoff gives: prior_essentiality 0.003, ulo_pgp_noprior 0.069, gd_coreset 0.046, prior_pathway 0.074, random 0.036.
   Caveat: these cutoffs use the same DepMap feature that prior_essentiality ranks on. Dropping prior_essentiality to 0 is therefore partly true by construction. The meaningful result is the size of the drop for each of the other methods.

## What this means for each method's IL-2 claims
- **prior_essentiality:** its "beats every method" result (N/E 0.137, All 0.277) comes from the metric. It does not reflect IL-2 discovery. Report it, but footnote that its N/E hits are all strongly essential genes that are not on CEGv2.
- **ulo_pgp_noprior / random-start graph:** it is partly affected. About 54% of its IL-2 N/E hits are strongly essential genes (effect < -0.5), against 42% for random. Its 0.116 falls to about 0.076–0.087 under the stricter post-hoc definitions, with a large SD.
  - Even after the drop, it stays above gd_coreset (0.046–0.066) and random.
  - It is roughly tied with prior_pathway (0.077–0.090).
  - So its IL-2 "WIN vs gd_coreset" probably survives in direction, but the margin shrinks. Its "WIN vs prior_pathway" is not robust to this check.
  - On IFNG it is much less affected (0.078 → about 0.069).
- **gd_coreset (best GeneDisco baseline):** it is also inflated (0.086 → about 0.046 at -0.5; 66% of its N/E hits are strongly essential), because it learns from the Achilles features.
- **prior_pathway, ulo_pgp, ulo_bandit, and the pathway-seeded variants:** they are essentially unaffected; their scores rise slightly under the stricter definitions.
- **Locked results stand as computed.** All numbers here are exploratory and post-hoc, and must not replace the locked table. I recommend adding a footnote to the IL-2 write-up saying that the CEGv2-based N/E metric does not exclude broadly essential genes, and that on IL-2 this favours essentiality-correlated methods.
