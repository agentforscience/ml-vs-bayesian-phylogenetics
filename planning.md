# Research Planning — Direction Budget

**Hypothesis.** ML and Bayesian phylogenetics agree on major clade structure in >85% of
cases but disagree on fine-scale (within-clade) branching order in 30–50% of outbreak
datasets, and these disagreements change the inferred transmission chain in ≥2 of 8
datasets in ways that would alter public health recommendations.

## Why the hypothesis is testable with the assembled resources

Eight curated viral outbreak alignments span an order of magnitude in *phylogenetic
information per taxon*, which is the variable that should control ML/Bayesian divergence:

| dataset | taxa | aln len | pars.-inf. sites | PIS/taxon | mean SNPs/pair | identical pairs |
|---|---|---|---|---|---|---|
| `sarscov2_usa_early2020`   | 120 | 29 921 |  57 | **0.47** | 10.7 | 0.008 |
| `fmdv_uk_2007`             |  45 |  8 176 |  23 | **0.51** |  7.8 | 0.041 |
| `mers_korea_2015`          |  34 | 30 119 |  20 | **0.59** |  4.2 | 0.071 |
| `ebov_sierraleone_2014`    | 110 | 18 959 |  78 | 0.71 |  9.7 | 0.015 |
| `h5n1_dairy_cattle_2024`   | 110 |  1 754 |  95 | 0.86 | 10.3 | 0.021 |
| `ebov_westafrica_makona`   | 150 | 18 979 | 205 | 1.37 | 19.0 | 0.002 |
| `fmdv_uk_2001`             |  43 |  8 196 |  93 | 2.16 | 23.9 | 0.000 |
| `zika_americas_2015_2016`  | 100 | 10 877 | 283 | **2.83** | 42.9 | 0.002 |

An unrooted tree of *n* taxa has *n*−3 internal branches. Five of eight datasets have
fewer than one parsimony-informative site per taxon, i.e. the topology is *chronically
underdetermined* — exactly the regime where the two inference paradigms should part
company. Datasets at the top of the table (Zika, FMDV 2001) act as positive controls
where both methods should converge.

## Pilot result (already obtained, `artifacts/`)

Running the full pipeline on two datasets confirms the machinery detects what the
hypothesis claims, and confirms the two-scale structure:

* **`mers_korea_2015`** — strongly supported splits: 2/2 ML splits recovered by Bayes,
  4/4 Bayesian splits recovered by ML, **0 conflicts**. But overall normalised
  RF = **0.548**, quartet distance = 0.134, and the Bayesian 95% credible set contains
  **14 179 distinct topologies** (max topology posterior = 0.000133). Major structure
  agrees; fine structure is essentially unidentified.
* **`fmdv_uk_2007`** — same pattern (0 strong-split conflicts, normalised RF = 0.548),
  **and the epidemiological conclusion changes**: temporal rooting of the ML tree puts
  the index premises at **IP2b**, the Bayesian tree at **IP1b**, reversing the direction
  of the IP1b/IP2b link. Cottam et al. (2008) identify IP1b as the source premises.

## Direction budget — enumeration, scoring, and pruning

Scored 1–5 on: **Ev** evidence from the literature that the effect exists;
**Rel** relevance to the stated hypothesis; **IG** expected information gain;
**Feas** implementation feasibility with the installed toolchain. Total = sum.

| # | Direction | Ev | Rel | IG | Feas | **Total** | Verdict |
|---|---|---|---|---|---|---|---|
| **D1** | **Two-scale topological concordance**: quantify strong-split agreement vs fine-scale RF/quartet disagreement across all 8 datasets, and locate the ML point estimate inside the Bayesian topological credible set | 5 | 5 | 5 | 5 | **20** | **KEEP** |
| **D2** | **Epidemiological consequence**: recompute introductions / transmission clusters / who-infected-whom / index case from the ML and Bayesian trees and count datasets where the action-relevant answer flips | 5 | 5 | 5 | 4 | **19** | **KEEP** |
| **D3** | **Uncertainty propagation & prior sensitivity**: repeat D1–D2 over whole tree *distributions* (UFBoot replicates vs posterior samples) rather than point estimates, and vary the branch-length prior | 5 | 4 | 5 | 4 | **18** | **KEEP** |
| D4 | Simulation study with a known true tree, to convert "disagreement" into "who is right" | 4 | 3 | 4 | 4 | 15 | pruned |
| D5 | Add BEAST2 time-tree/coalescent Bayesian arm alongside MrBayes | 3 | 3 | 3 | 3 | 12 | pruned |
| D6 | Substitution-model choice (JC/HKY/GTR, partitioning) as a separate factor | 3 | 2 | 3 | 4 | 12 | pruned |
| D7 | Full joint transmission-tree inference (TransPhylo / outbreaker2 / phybreak) | 4 | 3 | 3 | 1 | 11 | pruned |
| D8 | Sampling-density / subsampling sensitivity | 3 | 2 | 3 | 4 | 12 | pruned |
| D9 | Maximum parsimony and distance methods as extra arms | 2 | 2 | 2 | 5 | 11 | pruned |
| D10 | Pandemic-scale (10⁴–10⁶ taxa) comparison with UShER/matOptimize | 3 | 1 | 3 | 2 | 9 | pruned |
| D11 | Recombination / homoplasy masking as a confounder | 3 | 2 | 3 | 3 | 11 | pruned |
| D12 | Alignment-method sensitivity (MAFFT vs alternatives) | 2 | 1 | 2 | 4 | 9 | pruned |

### Reasons for pruning

* **D4** — the single most valuable *addition*, and the strongest candidate for promotion
  if D1–D3 return an ambiguous answer. Pruned only because the hypothesis is stated about
  *agreement between methods on real outbreak data*, not about accuracy against truth;
  answering it needs a simulation framework that is out of scope for the stated claim.
* **D5** — MrBayes already supplies the Bayesian arm and a full posterior topology
  distribution. BEAST2 would change *two* things at once (inference engine **and**
  tree prior), confounding the comparison. BEAST 2.7 + JDK 21 are nevertheless installed
  (`code/_src/beast`) so this can be added without new setup.
* **D6** — model choice is a genuine third factor but it multiplies the design; it is
  partly controlled already, since IQ-TREE selects the model by BIC and the MrBayes block
  uses GTR+I+G. Retained only as a fixed, documented control, not as a swept variable.
* **D7** — TransPhylo, outbreaker2 and phybreak are all R packages and **R is not
  installable in this environment (no root, and conda is disallowed)**. The repositories
  are cloned for reference, and `code/mlbayes/epi.py` implements the tree-derived
  epidemiological readouts natively in Python instead. This is the one real capability gap.
* **D8, D9, D11, D12** — real but second-order effects that would dilute the sample of 8
  datasets across too many cells to support the "≥2 of 8" claim.
* **D10** — pandemic-scale inference is a different computational problem; MCMC is not
  feasible there at all, so there is no Bayesian arm to compare against.

**Do not expand the search space** beyond D1–D3 unless new evidence invalidates this
ranking; if it does, update the ranking here and record the reason in `STATE.md`.

## Pre-registered analysis plan for D1–D3

**D1 — topological concordance.** For each dataset compute: normalised RF and quartet
distance between the ML tree and the Bayesian majority-rule consensus; the fraction of
UFBoot ≥ 95 splits present in the Bayesian tree and of PP ≥ 0.95 splits present in the ML
tree (→ the ">85% major clade agreement" claim); the fine-scale RF restricted to weakly
supported splits (→ the "30–50% fine-scale disagreement" claim); the size of the Bayesian
95% topological credible set and whether the ML topology is in it; and the PP-vs-UFBoot
relationship on shared splits.

**D2 — epidemiological consequence.** For each dataset, root by root-to-tip regression and
compute the readouts in `code/mlbayes/epi.py`: number of independent introductions per
region, transmission-cluster membership, directed group-to-group transitions, and the
index group. Count a dataset as "conclusion-changing" under a pre-declared rule —
a change in the *count* of introductions or clusters, a change in the index group, or a
reversal of a directed transmission edge. Target: ≥2 of 8.

**D3 — uncertainty propagation.** Repeat D2 over the 1000 UFBoot replicate trees and over
a thinned posterior sample, giving two *distributions* of each epidemiological quantity;
compare the distributions rather than the point estimates. Run the branch-length-prior arm
(`EXPONENTIAL_MODEL` vs `DEFAULT_MODEL` in `code/mlbayes/run_bayes.py`) to separate
genuine ML-vs-Bayes differences from prior artefacts.

### Confounder that must be controlled (measured, not hypothetical)

On `fmdv_uk_2007` the MrBayes **default independent exponential(10) branch-length prior
inflates consensus tree length 1070×** relative to ML (7.12 vs 0.0067 subs/site);
`exponential(1000)` gives 0.063 and the compound Dirichlet `gammadir(1,1000,1,1)` gives
0.0049. Patristic distances feed transmission-cluster thresholds directly, so the naive
default would manufacture an "epidemiological difference" that is a prior artefact rather
than an ML-vs-Bayes difference. `DEFAULT_MODEL` therefore uses the compound Dirichlet
prior and the exponential prior is retained as an explicit D3 sensitivity arm.

---

# Amendment 1 (experiment_runner phase, 2026-09-09) — a within-method null

**Status: added BEFORE any production result was inspected.** The amendment is
driven entirely by the pilot in `artifacts/pipeline_test/`, which is the evidence
`STATE.md` records at the end of the resource-finding phase.

## Why the pre-registered plan is not sufficient on its own

Two facts from the pilot make a bare ML-vs-Bayesian comparison uninterpretable:

1. **The posterior never samples the same topology twice.** On both pilot datasets
   the reported credible set contained 14 180 of 15 002 sampled trees and the
   maximum topology probability was 6.7 × 10⁻⁵ = 1/15 002. The "95 % credible set
   size" is therefore reading out the MCMC *sample size*, not the concentration of
   the posterior. Any statistic of that form is uninformative here.
2. **The same program disagrees with itself.** Re-running IQ-TREE on identical
   data changing only the random seed gives RF = 14–32 on `fmdv_uk_2007`
   (normalised 0.17–0.38), against a between-method normalised RF of 0.55–0.62.
   A raw between-method distance therefore has no scale: part of it is simply
   heuristic search noise and Monte-Carlo noise, not "model choice".

A finding that "ML and Bayes disagree" is only evidence about **model choice** if
the disagreement is larger than what the *same* method produces when merely re-run.
Without that comparison the headline claim cannot be attributed to the methods.

## What is added

Five null distributions, computed identically to the between-method statistic and
carried through *both* the topological analysis (D1) and the epidemiological
conclusion-change rule (D2):

| id | null | construction |
|---|---|---|
| **N1** | MCMC replicate | majority-rule consensus of MrBayes run 1 vs run 2 (two independent chains, same prior, same data) |
| **N2** | ML search replicate | all 10 pairs among 5 IQ-TREE runs differing only in `-seed` |
| **N3** | posterior spread | Bayesian consensus vs each posterior sample tree |
| **N4** | bootstrap spread | ML tree vs each UFBoot replicate |
| **P**  | prior arm | compound-Dirichlet vs exponential(10) branch-length prior (already planned as D3) |

## Amended decision rules

* **D1.** The between-method normalised RF is reported as a **percentile of N2,
  N3 and N4**, not as a bare number. The two-scale claim is tested by contrasting
  confident-split recovery (major structure) with the fine-scale RF, exactly as
  pre-registered — the nulls only add the missing scale.
* **D2.** The pre-declared conclusion-change rule (a change in the count of
  introductions or clusters, a change in the index group, or a reversal of a
  directed transmission edge) is applied unchanged to `(ML, Bayes)`, **and also**
  to every N1 and N2 pair. The claim "≥ 2 of 8 datasets are conclusion-changing"
  is then judged against the number of datasets in which re-seeding *one* method
  is already conclusion-changing.
* Sample-size-bound credible-set statistics are retained for comparability with the
  pilot but are **not** used as evidence; the number of distinct topologies in the
  sample and the maximum topology frequency are reported alongside them so the
  saturation is visible.

## Direction budget

This is an addition **inside** kept directions D1 and D3 (it is uncertainty
propagation, with the uncertainty of each method used as its own yardstick). No
pruned direction is revived and the search space is not widened: D4–D12 remain
pruned for the reasons given above.

## Other deviations from the resource-phase plan, and why

* **Chain lengths are per-dataset, not a uniform 2 000 000 generations.** A
  measured speed calibration (`artifacts/speedcal/speed.json`) showed a 21×
  spread in cost per generation across the eight alignments. Generations are set
  per dataset to give each chain a comparable slice of compute; the number of
  *retained samples* is held constant at 10 000/run instead, because credible-set
  size is a function of sample count and must not vary between datasets.
  Convergence is checked per dataset (ASDSF < 0.01, max PSRF < 1.02, min ESS > 200)
  and reported, not assumed.
* **`epi.root_by_tip_dates` is replaced by `src/fastroot.rtt_root`** for the
  distribution analyses. The reference implementation clones the tree once per
  candidate root edge, which is unusable for 1000-tree sets. The replacement gets
  the same TempEst criterion from a two-pass dynamic program and additionally
  optimises the root *position within* the chosen edge. `src/test_fastroot.py`
  validates it against the reference: same root edge on 2 of 3 test datasets and
  never a lower R².
* **A primary trait is designated per dataset** (`PRIMARY_TRAIT` in
  `src/analysis_epi.py`) — the grouping variable a public-health decision would
  actually turn on — so the conclusion-change rule is not diluted across
  incidental metadata columns. `mers_korea_2015` has no such variable (single
  country, no premises labels) and is analysed on transmission clusters only.
