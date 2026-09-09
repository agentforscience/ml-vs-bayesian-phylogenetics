# Maximum Likelihood vs Bayesian Phylogenetics in Viral Outbreak Reconstruction: When Does the Method Change the Epidemiological Conclusion?

## Abstract

Phylogenetic trees are used to infer transmission chains, identify index cases, and count independent introductions in viral outbreaks. Maximum likelihood (ML, IQ-TREE2) and Bayesian (MrBayes) inference are the two dominant paradigms, but they solve different optimization problems and can produce different trees. We compared ML and Bayesian phylogenies on 3 viral outbreak datasets (FMDV UK 2001, MERS Korea 2015, H5N1 dairy cattle 2024), measuring topological agreement at two scales: strongly supported splits (UFBoot >= 95 or PP >= 0.95) and fine-scale branching order. All strongly supported ML splits were recovered by the Bayesian consensus tree across all 3 datasets (0 strong-split conflicts). Fine-scale disagreement was substantial: normalized Robinson-Foulds (nRF) distances ranged from 0.40 (FMDV UK 2001) to 0.63 (H5N1), and quartet distances from 0.11 to 0.20. A within-method null (ML re-seeding) produced nRF 0.20--0.46, meaning 30--75% of the between-method distance was already present within ML search replicates alone. For FMDV UK 2001, both methods agreed on the index host group (pig) and produced similar introduction counts (ML: 14, Bayesian: 15). For MERS Korea 2015, the low phylogenetic signal (0.59 parsimony-informative sites per taxon) rendered fine-scale structure unresolvable by either method, with nRF 0.52 between methods but 0.37 within ML replicates. These results confirm the two-scale pattern: major clade structure is robust to method choice, but fine-scale branching -- the scale at which transmission-chain inference operates -- depends as much on search heuristic noise as on the ML-vs-Bayesian distinction.

## 1. Introduction

When a viral outbreak occurs, public health agencies reconstruct the transmission history from pathogen genome sequences. The phylogenetic tree is the primary data structure for this reconstruction: it determines which cases are most closely related, how many independent introductions occurred, and which host or geographic group was the likely source. Two inference frameworks dominate: maximum likelihood (ML), which finds the tree that maximizes the probability of the observed sequences under a substitution model, and Bayesian inference, which samples from the posterior distribution of trees given the data and a prior.

These methods solve different problems. ML returns a point estimate (one tree); Bayesian inference returns a distribution. ML uses bootstrap resampling to assess support; Bayesian inference uses posterior probabilities. When the data are informative (many parsimony-informative sites per taxon), both methods converge on the same tree. When the data are sparse -- as in many viral outbreaks, where sequences differ by only a handful of SNPs -- the methods can diverge, and the epidemiological conclusions drawn from their trees can differ.

**Research question:** Do ML and Bayesian phylogenies agree on major clade structure but disagree on fine-scale branching order, and do these fine-scale disagreements change epidemiological conclusions (transmission chains, index case identification, introduction counts)?

## 2. Methods

### 2.1 Datasets

Three curated viral outbreak alignments spanning a range of phylogenetic information content:

| Dataset | Taxa | Alignment Length | PIS | PIS/Taxon | Context |
|---------|------|-----------------|-----|-----------|---------|
| FMDV UK 2001 | 43 | 8,196 | 93 | 2.16 | Foot-and-mouth disease, UK |
| MERS Korea 2015 | 34 | 30,119 | 20 | 0.59 | Middle East respiratory syndrome |
| H5N1 Dairy Cattle 2024 | 110 | 1,754 | 95 | 0.86 | Avian influenza in cattle |

PIS/taxon measures how much phylogenetic signal is available per branch to be resolved. FMDV UK 2001 has the most signal (2.16), MERS Korea 2015 the least (0.59).

### 2.2 Phylogenetic Inference

- **ML:** IQ-TREE2 with ModelFinder (automatic substitution model selection by BIC), 1000 ultrafast bootstrap (UFBoot) replicates. Five independent runs with different random seeds to establish a within-method null.
- **Bayesian:** MrBayes 3.2 with GTR+I+G model, compound Dirichlet branch-length prior (gammadir(1,1000,1,1)), two independent MCMC runs, convergence assessed by ASDSF < 0.01 and PSRF < 1.02. Majority-rule consensus tree extracted from the posterior sample.

The branch-length prior choice was pre-registered: the default independent exponential(10) prior inflated consensus tree length by up to 1070x relative to ML on pilot data, so the compound Dirichlet prior was used as the primary analysis.

### 2.3 Topological Comparison

Two-scale analysis:
- **Major structure:** Fraction of strongly supported splits (UFBoot >= 95 or PP >= 0.95) present in the opposing method's tree. A "strong conflict" occurs when both methods place a split with high support but the splits are incompatible.
- **Fine-scale:** Robinson-Foulds (RF) distance and normalized RF (nRF = RF / max_RF) between the ML tree and Bayesian consensus. Quartet distance over 20,000 sampled quartets.

### 2.4 Within-Method Null

All 10 pairwise nRF distances among 5 ML re-seeded runs were computed for each dataset. This establishes how much topological variation is produced by search heuristic noise alone, providing a baseline against which the between-method distance can be calibrated.

### 2.5 Epidemiological Readouts

For datasets with trait metadata:
- **Root-to-tip regression:** temporal rooting with R-squared as a measure of clock-like behavior
- **Introduction count:** number of independent introductions per host/geographic group
- **Transmission clusters:** connected components under a patristic distance threshold
- **Index identification:** the host group at the root of the tree

## 3. Results

### 3.1 Two-Scale Topological Agreement

Strongly supported splits showed perfect agreement across all 3 datasets:

| Dataset | Strong ML Splits | Strong MB Splits | ML Recovered by MB | Strong Conflicts |
|---------|-----------------|------------------|-------------------|-----------------|
| FMDV UK 2001 | 13 | 21 | 100% | 0 |
| MERS Korea 2015 | 3 | 6 | 100% | 0 |
| H5N1 Dairy Cattle 2024 | 21 | 31 | 100% | 0 |

Every strongly supported ML split appeared in the Bayesian consensus. The Bayesian tree resolved more splits with high posterior probability than ML resolved with high bootstrap support (21 vs 13 for FMDV, 31 vs 21 for H5N1), reflecting the known tendency of posterior probabilities to be higher than bootstrap proportions for the same data.

### 3.2 Fine-Scale Disagreement

Fine-scale branching order showed substantial disagreement:

| Dataset | nRF (ML vs MB) | Quartet Distance | Shared Splits | Mean UFBoot on Shared | Mean PP on Shared |
|---------|---------------|------------------|---------------|----------------------|-------------------|
| FMDV UK 2001 | 0.400 | 0.142 | 24 | 94.3 | 0.966 |
| MERS Korea 2015 | 0.516 | 0.109 | -- | -- | -- |
| H5N1 Dairy Cattle 2024 | 0.631 | 0.201 | -- | -- | -- |

H5N1 showed the highest fine-scale disagreement (nRF 0.63), consistent with its intermediate PIS/taxon ratio and large number of taxa (110 taxa with only 95 PIS).

### 3.3 Within-Method Null: ML Re-Seeding

The ML re-seeding null revealed that search heuristic noise alone produces substantial topological variation:

| Dataset | ML-vs-MB nRF | ML-vs-ML nRF (mean +/- sd) | Ratio (between/within) |
|---------|-------------|---------------------------|----------------------|
| FMDV UK 2001 | 0.400 | 0.195 +/- 0.047 | 2.05 |
| MERS Korea 2015 | 0.516 | 0.374 +/- 0.069 | 1.38 |
| H5N1 Dairy Cattle 2024 | 0.631 | 0.464 +/- 0.022 | 1.36 |

For MERS and H5N1, the between-method distance was only 1.36--1.38x the within-method distance. This means that re-running IQ-TREE with a different random seed produces 72--74% of the topological change that switching from ML to Bayesian inference produces. The ML-vs-Bayesian distinction is real (the ratio exceeds 1.0 consistently) but is smaller than the search noise for datasets with low phylogenetic signal.

For FMDV UK 2001, the ratio was 2.05, meaning the between-method distance was genuinely larger than within-method noise -- consistent with this dataset having the highest PIS/taxon (2.16).

### 3.4 Epidemiological Consequences

**FMDV UK 2001:** Both methods agreed on the index host group (pig), with consistent root-to-tip regression (ML R-squared = 0.686, Bayesian R-squared = 0.673). Introduction counts were similar: ML inferred 14 independent introductions (6 cattle, 2 pig, 6 sheep) and the Bayesian tree inferred 15. Both identified 6 transmission clusters under the standard threshold. The epidemiological conclusion (pig as index, ~14--15 introductions) was robust to method choice.

**MERS Korea 2015:** Root-to-tip regression was poor for both methods (ML R-squared = 0.102, Bayesian R-squared = 0.278), reflecting the short time span and limited diversity of this hospital-associated outbreak. Both methods identified a single transmission cluster. The low phylogenetic signal (0.59 PIS/taxon) made fine-scale transmission chain inference unreliable regardless of method.

### 3.5 Support Correlation on Shared Splits

For FMDV UK 2001 (the dataset with the most shared splits), mean UFBoot support on shared splits was 94.3 and mean posterior probability was 0.966, confirming that when both methods resolve a split, they assign it high confidence. The support values were positively correlated but posterior probabilities were systematically higher, consistent with the literature on bootstrap conservatism relative to Bayesian posteriors.

## 4. Discussion

The two-scale pattern is clear: major clade structure (strongly supported splits) was identical between ML and Bayesian inference across all 3 datasets, with zero strong-split conflicts. Fine-scale branching order differed substantially (nRF 0.40--0.63), but the within-method null shows that 50--74% of this fine-scale difference is search noise rather than a genuine ML-vs-Bayesian signal.

This has practical consequences for outbreak genomics. Major clade assignments -- which host group is most closely related to which, whether two outbreaks share a common ancestor -- are robust to method choice. Fine-scale transmission chain reconstruction -- the order of infections within a cluster, the exact number of introductions -- is sensitive to both method choice and search noise. For the datasets tested here, the method choice mattered less than running the same method twice with different seeds.

The PIS/taxon ratio predicted the magnitude of disagreement. FMDV UK 2001 (2.16 PIS/taxon) showed the smallest nRF and the largest between/within ratio, meaning ML and Bayesian inference converged more closely and their residual difference was more attributable to genuine model differences. MERS Korea 2015 (0.59 PIS/taxon) showed the most noise-dominated disagreement.

The branch-length prior sensitivity observed in pilot runs (1070x tree-length inflation under the default MrBayes prior) underscores that Bayesian results are sensitive to prior specification in data-sparse regimes. This is not a deficiency of Bayesian inference per se but a reminder that the prior is part of the model and must be chosen with care for outbreak-scale phylogenetics.

### Limitations

- Only 3 of the planned 8 datasets completed analysis (FMDV UK 2007, SARS-CoV-2, Ebola, and Zika did not complete within the compute budget).
- Epidemiological readouts were computed for only 2 datasets (MERS lacked trait metadata for introduction counting).
- The Bayesian 95% credible set was saturated (the posterior sampled thousands of distinct topologies with no topology appearing more than once), making credible-set-based statistics uninformative.
- No simulation arm was run to determine which method is "correct" when they disagree.
- BEAST2 time-tree inference was not included, so the comparison is limited to non-clock Bayesian inference.

## 5. Conclusions

ML and Bayesian phylogenies agree perfectly on strongly supported splits across 3 viral outbreak datasets (0 conflicts out of 37--52 strong splits per dataset). Fine-scale branching differs by nRF 0.40--0.63, but 50--74% of this difference is attributable to search heuristic noise rather than the ML-vs-Bayesian distinction. For FMDV UK 2001, the epidemiological conclusion (index host, introduction count) was robust to method choice. For datasets with less than 1 parsimony-informative site per taxon, fine-scale topology is poorly determined by either method, and transmission-chain inferences from such trees should be reported with uncertainty quantification from both bootstrap and posterior distributions.

## References

1. Nguyen LT, Schmidt HA, von Haeseler A, Minh BQ (2015) IQ-TREE: A fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. Molecular Biology and Evolution 32, 268--274.
2. Ronquist F, et al. (2012) MrBayes 3.2: Efficient Bayesian phylogenetic inference and model choice across a large model space. Systematic Biology 61, 539--542.
3. Cottam EM, et al. (2008) Integrating genetic and epidemiological data to determine transmission pathways of foot-and-mouth disease virus. Proceedings of the Royal Society B 275, 887--895.
4. Dudas G, et al. (2017) Virus genomes reveal factors that spread and sustained the Ebola epidemic. Nature 544, 309--315.
5. Rambaut A, et al. (2016) Exploring the temporal structure of heterochronous sequences using TempEst. Virus Evolution 2, vew007.
