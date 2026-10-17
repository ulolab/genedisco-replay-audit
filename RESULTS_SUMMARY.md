# Results summary (plain English)

**Bottom line: we make no win claim.** Our two no-LLM methods lost or tied on the clean test screen (IL-2) and won only narrowly on the screen we had used for practice (IFNG). The most useful thing we found is about the metric, not about our methods.

## What we tested
We replayed the GeneDisco benchmark as used in BioDiscoveryAgent (128 genes × 5 rounds, 10 seeds), with no language model. We compared two Ulo methods (`ulo_bandit` and `ulo_pgp_noprior`) against GeneDisco baselines that we re-ran ourselves and against two simple rankings that ignore the experiment entirely: DepMap essentiality and pathway membership. The plan, including the win rule, was fixed before the final runs (not publicly pre-registered). IFNG was used for practice and IL-2 was kept as the clean test.

**Win rule:** beat the best re-run GeneDisco baseline *and* both simple rankings on the non-essential (N/E) hit ratio, using paired bootstrap 95% confidence intervals.

## IL-2 (clean test): no win
- `ulo_bandit` (0.077 N/E) **lost** to Coreset (0.086).
- `ulo_pgp_noprior` (0.116) beat Coreset, but **lost** to the essentiality ranking (0.137).
- The essentiality ranking, which never looks at any experiment result, scored higher than every method we ran.

## IFNG (practice screen): narrow win, weaker evidence
- `ulo_bandit` (0.080) beat Coreset (0.074; difference +0.005 unrounded, CI [+0.002, +0.009]) and both simple rankings. That meets the win rule, but this screen was used during development, so it counts for much less than IL-2.
- `ulo_pgp_noprior` (0.078) tied with Coreset.
- Both methods are still below the published Sonnet 3.5 agent (0.107; copied from the paper, so not a head-to-head comparison).

## The finding that matters
On IL-2, the N/E metric removes only the 684 CEGv2 core-essential genes. Many IL-2 hits are other strongly essential "housekeeping" genes (ribosome, proteasome, DNA replication). These are genuine hits in the screen. The point is that the metric counts them as non-essential, so a ranking by essentiality alone scores well. Details are in `ESSENTIALITY_ODDITY.md`, and we raised it as a question with the authors in [BioDiscoveryAgent issue #8](https://github.com/snap-stanford/BioDiscoveryAgent/issues/8).

All numbers come from `final_IL2.txt` and `final_IFNG.txt`, and the README shows how to reproduce them.
