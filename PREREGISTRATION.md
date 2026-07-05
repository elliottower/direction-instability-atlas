# Pre-registration: Direction Instability Atlas

**Title:** Cross-modal direction instability atlas for perturbation biology

**Authors:** Elliot Tower

**Date:** 2026-07-05 (pre SHA freeze)

**Commit SHA:** _(to be filled after commit)_

---

## Overview

This document pre-registers the hypotheses, decision criteria, and analysis
pipeline for the Direction Instability (DI) Atlas paper. All experiments,
thresholds, and statistical controls are frozen before any DI values are
computed on new datasets. DI values for LINCS L1000 and JUMP-CP compounds
were previously computed in the bracket-norm-validity repository.

## Background and motivation

The combined confound paper (bracket-norm-validity) demonstrates that
variance-based, subspace-based, and distribution-based perturbation
variability metrics are confounded by estimation quality (effect magnitude,
sample size). Cosine-based direction instability (DI) is immune to this
confound by construction. This atlas paper assumes that result and computes
DI across three new datasets that the confound paper did not analyze:
Tahoe-100M (transcriptomics), JUMP-CP CRISPR knockouts (morphology),
and JUMP-CP ORF overexpression (morphology). PRISM viability is included
as a complementary (non-DI) validation.

---

## Datasets

### DI-proper datasets (high-dimensional signatures, cosine meaningful)

| Dataset | N perturbations | N features | N contexts | Source |
|---------|----------------|------------|------------|--------|
| LINCS L1000 | 8,949 drugs | 978 genes | cell lines (varies) | Pre-computed |
| JUMP-CP compounds | 25,254 | 3,180 morphological | source contexts | Pre-computed |
| Tahoe-100M | 379 drugs | ~62,710 genes | 50 cell lines | HuggingFace |
| JUMP-CP CRISPR | 7,948 genes (>= 5 reps) | 259 (PCA) | plate replicates | Cell Painting Gallery |
| JUMP-CP ORF | 15,121 genes (>= 5 reps) | 722 | plate replicates | Cell Painting Gallery |

### Complementary dataset (NOT DI-proper)

| Dataset | N perturbations | Readout | Source |
|---------|----------------|---------|--------|
| PRISM 24Q2 | 6,575 drugs | 1D viability (LFC) | Figshare |

PRISM viability is scalar per cell line. Cosine is undefined for 1D
vectors. Viability coefficient of variation (CV) is computed instead.

### Cross-dataset overlap (verified)

Drug overlaps were computed by name matching, then verified by PubChem
InChIKey connectivity-layer (first 14 characters) matching via the
PubChem REST API. 10 false-positive name matches (different molecules
sharing a name, such as different stereoisomers with non-matching
InChIKey connectivity) were removed.

| Pair | Overlap | Method |
|------|---------|--------|
| LINCS ∩ Tahoe | **103** | Name + InChIKey verified |
| LINCS ∩ JUMP-CP compounds | 330 | Name matching |
| Tahoe ∩ JUMP-CP compounds | **0** | SMILES matching (Tahoe clinical drugs vs JUMP screening library) |
| LINCS ∩ PRISM | 1,726 | Name matching |
| Tahoe ∩ PRISM | 198 | Name matching |
| JUMP-CP CRISPR ∩ ORF (gene symbols, >= 5 reps) | 5,220 | Gene symbol |

**A3 (Tahoe ∩ JUMP-CP) is dropped** from primary tests: Tahoe contains
379 clinical-stage drugs while JUMP-CP profiles 25,254 screening-library
compounds. These compound sets have zero structural overlap.
Combined with Experiment B (dropped for structural power failure),
this reduces the primary test count from 6 to 4.

---

## Direction instability definition

For a perturbation with K context-specific signatures
s_1, ..., s_K in R^d (K >= 5):

    DI = 1 - mean_{i < j} cos(s_i, s_j)

where cos(s_i, s_j) = (s_i . s_j) / (||s_i|| ||s_j||).

**Magnitude-corrected DI:** DI residualized against mean signature L2 norm
(||s||_mean) via OLS regression. Residuals + intercept = corrected DI.
This is the PRIMARY atlas metric.

**Bootstrap confidence intervals:** For each perturbation, resample the
set of K contexts with replacement 1,000 times. Recompute DI on each
bootstrap sample. Report 2.5th and 97.5th percentiles as 95% CI.

---

## Pre-registered experiments

### Multiple comparisons correction

Four primary tests (A3 dropped: zero Tahoe-JUMP overlap; B dropped:
structurally underpowered, see below).
Bonferroni-corrected alpha = 0.05 / 4 = **0.0125**.

Conservative because the two Exp A pairs share the LINCS DI vector
(not fully independent). We err on the strict side deliberately.

Exploratory analyses (PRISM concordance, Exp D directional
divergence) use uncorrected p-values with explicit "exploratory"
labels.

### Experiment A: Cross-dataset DI concordance

**Hypothesis:** Drugs present in multiple DI-proper datasets show
correlated corrected DI.

**Confounder control (MANDATORY):** Every concordance pair reports:
- Raw Spearman rho (bivariate)
- **Partial Spearman rho controlling for n_A and n_B** (number of
  contexts in each dataset, as two separate covariates). This is the
  primary result. Conditioning on min(n_A, n_B) discards information.

**Partial Spearman method:** Pearson-residualize raw X and Y against
the covariates Z, then compute Spearman rho on the residuals. The
alternative "rank then residualize" method is biased: at n=5000 with
true partial rho = 0, it gives observed partial rho ~ 0.05-0.25.
The correct method gives < 0.01 (verified in tests).

**Decision criteria (alpha = 0.0125):**

All four primary tests use the same criterion form: 95% CI lower
bound of (partial) Spearman rho > 0.05, p < 0.0125. The uniform
0.05 floor ensures effects are non-trivially positive; magnitude
interpretation uses the bins below.

| Test | Pair | Criterion | Rationale |
|------|------|-----------|-----------|
| A1 | LINCS ∩ Tahoe (n=103) | CI lower > 0.05, p < 0.0125 | Same modality |
| A2 | LINCS ∩ JUMP-CP (n=330) | CI lower > 0.05, p < 0.0125 | Cross-modality |

A3 (Tahoe ∩ JUMP-CP) dropped: zero compound overlap between clinical
drugs (Tahoe) and screening library (JUMP-CP).

All results report the point estimate with 95% Fisher z confidence
interval. Effect-size interpretation follows standardized bins:
rho >= 0.30 = "concordant", 0.15-0.30 = "weakly concordant",
rho < 0.15 = "not replicated".

**Power analysis (5,000 simulations, corrected partial Spearman, 2 covariates):**

| Test | n | True rho | Criterion | Power |
|------|---|----------|-----------|-------|
| A1 | 103 | 0.20 | CI lower > 0.05 | 0.291 |
| A1 | 103 | 0.25 | CI lower > 0.05 | 0.473 |
| A1 | 103 | 0.30 | CI lower > 0.05 | 0.658 |
| A1 | 103 | 0.35 | CI lower > 0.05 | 0.825 |
| A1 | 103 | 0.40 | CI lower > 0.05 | 0.935 |
| A2 | 330 | 0.10 | CI lower > 0.05 | 0.131 |
| A2 | 330 | 0.15 | CI lower > 0.05 | 0.407 |
| A2 | 330 | 0.20 | CI lower > 0.05 | 0.729 |
| A2 | 330 | 0.25 | CI lower > 0.05 | 0.935 |
| A2 | 330 | 0.30 | CI lower > 0.05 | 0.992 |

**Power honesty note:** At n=103, A1 achieves 80% power only for
true effects >= 0.35. This is a structural constraint of the 103-drug
overlap. We report this transparently. The interpretation guide above
handles the case where the true effect is real but too small to detect.

**Interpretation guide:**
- Both pass: DI is a property of the perturbation, not the assay.
- Same-modality passes, cross-modality fails: DI is modality-specific.
  Atlas remains useful per-modality; concordance claim weakens.
- Same-modality fails: DI is platform-specific. Seriously weakens atlas.
- Raw rho passes but partial rho fails: concordance is confound-driven.

### Experiment PRISM (EXPLORATORY): PRISM viability CV vs LINCS DI

**Hypothesis:** Drugs with low DI in LINCS show low viability CV in PRISM.

**Metric:** Spearman rho between LINCS DI and PRISM viability CV for
shared drugs (n=1,726).

**Threshold:** rho > 0.15. Weak bar for a weak claim: mechanism
conservation predicts viability consistency, but imperfectly.

**Label:** EXPLORATORY. Uncorrected p-value. Not a primary test.

### Experiment B: MOA stratification in Tahoe — DROPPED

**Status:** Dropped before any data analysis.

**Reason:** Tahoe drug_metadata contains only 3 moa_broad categories
(inhibitor/antagonist: 276, unclear: 76, activator/agonist: 27).
Keyword matching for "machinery" vs "receptor" yields 7 vs 17 drugs.
Power simulation: Cohen's d = 0.5 at n1=7, n2=17 gives 6% power.
Even the most optimistic effect size (d=1.0) gives only 30% power.
This is a structural impossibility, not a design choice.

Widening the machinery definition to rescue power would be
outcome-driven reclassification — the exact error this framework
is designed to prevent.

### Experiment C: CRISPR DI predicts essentiality (SCIENTIFIC HEART)

**Hypothesis:** In JUMP-CP CRISPR knockouts, genes essential in more
cell types (DepMap Chronos < -0.5) show higher morphological DI.

**Confounder control:** Partial Spearman rho controlling for n_contexts
(number of plate replicates per gene in JUMP-CP CRISPR). Computed
using Pearson-residualize-then-Spearman method (see Experiment A).

**Decision criteria:** 95% CI lower bound > 0.05, p < 0.0125
(same uniform criterion as A1, A2, D).

At n=7948, even trivially small effects (rho ~ 0.02) reach
statistical significance. The CI floor of 0.05 ensures the effect
is at least modestly above zero.

**Interpretation of effect sizes:** rho >= 0.15 = "biologically
meaningful"; rho in [0.05, 0.15] = "weak — statistically detectable
but explaining < 2% of variance, do not lean on p-value";
rho < 0.05 = "negligible".

**Power analysis (3,000 simulations, 1 covariate):**

| True rho | Criterion | Power |
|----------|-----------|-------|
| 0.05 | CI lower > 0.05 | 0.015 |
| 0.08 | CI lower > 0.05 | 0.661 |
| 0.10 | CI lower > 0.05 | 0.984 |
| 0.15 | CI lower > 0.05 | 1.000 |

At n=7948, any effect rho >= 0.10 is detected with > 98% power.

**Essentiality source:** DepMap Public 24Q4 Chronos dependency scores.
Essentiality breadth = fraction of cell lines with Chronos < -0.5.

### Experiment D: Knockout vs overexpression DI (SCIENTIFIC HEART)

**Hypothesis:** Genes with high CRISPR knockout DI also have high ORF
overexpression DI. Context-dependence is a property of the gene, not
the perturbation direction.

**Universe:** 5,220 genes with >= 5 replicates in both CRISPR and ORF.

**Permutation null (MANDATORY):** CRISPR and ORF share plate/batch
structure. Null: shuffle ORF gene labels (preserving plate structure),
recompute rho, repeat 1,000 times. Report:
- Observed Spearman rho
- Permutation null distribution (mean, SD, 95th percentile)
- Empirical p-value (fraction of permutations >= observed rho)

**Decision criteria:** BOTH conditions must hold:
1. 95% CI lower bound of Spearman rho > 0.05
2. Observed rho exceeds 95th percentile of permutation null
Alpha = 0.0125.

The CI-lower-bound criterion (matching A1, A2, C) ensures internal
consistency across all four primary tests. At n=5220, the CI is
narrow (SE ~ 0.014), so the 0.05 floor guards against trivially
small effects while providing excellent power for real effects.

**Power analysis (3,000 simulations, no covariates):**

| True rho | Criterion | Power |
|----------|-----------|-------|
| 0.05 | CI lower > 0.05 | 0.019 |
| 0.08 | CI lower > 0.05 | 0.481 |
| 0.10 | CI lower > 0.05 | 0.912 |
| 0.15 | CI lower > 0.05 | 1.000 |
| 0.20 | CI lower > 0.05 | 1.000 |

At n=5220, any effect rho >= 0.10 is detected with > 90% power.

**Directional divergence (EXPLORATORY):** For genes where
|DI_CRISPR - DI_ORF| > 1 SD, classify into:
- Knockout-stable / overexpression-unstable: predicted biology =
  dosage-sensitive gain-of-function genes
- Knockout-unstable / overexpression-stable: predicted biology =
  genes with redundant paralogs

Report counts and Gene Ontology enrichment. No pass/fail threshold.
Label: EXPLORATORY.

---

## Tahoe processing specification

### Pseudobulk construction
- Source: Tahoe-100M pseudobulk DE tables (HuggingFace, 1,026 parquet shards)
- Each shard contains: gene_name, log2FoldChange, lfcSE, n_cells_trt,
  Cell_Name_Vevo, drug
- Signature vector for drug d in cell line c: log2FoldChange values across
  all genes, forming a vector in R^{~62,710}

### Minimum cell threshold
- **Exclude drug-cell_line pairs with n_cells_trt < 20** (treatment arm)
- Matching Perturb-seq threshold from combined paper
- Supplementary figure: DI as a function of cells-per-pseudobulk to
  confirm flat-bin immunity

### Processing location
- 89 GB total; processed on Modal (too large for local disk)
- Stream parquet shards, extract signatures, compute DI per drug
- Download only result files (few KB)

---

## JUMP-CP CRISPR/ORF processing specification

### Feature spaces
- CRISPR: 259 features (post-sphering, harmony, PCA correction)
- ORF: 722 features (post-sphering, harmony, no PCA)
- Different feature spaces are acceptable because DI is computed WITHIN
  each modality. DI is dimensionless (mean cosine distance) and comparable
  across feature spaces.

### Context definition
- Context = plate replicate (well-level profiles, median 5 per gene)
- Genes with < 5 replicates excluded (7,948 CRISPR and 15,121 ORF qualify)
- No cross-plate aggregation before DI computation

### Controls excluded
- CRISPR: "no-guide" and "non-targeting" perturbations excluded
- ORF: perturbations with Metadata_pert_type != "trt" excluded

---

## Atlas columns (the published resource)

| Column | Description |
|--------|-------------|
| `perturbation_id` | Dataset-specific identifier |
| `perturbation_name` | Human-readable name (drug name or gene symbol) |
| `dataset` | Source dataset |
| `di_corrected` | **PRIMARY.** Magnitude-corrected DI |
| `di_corrected_ci_low` | Lower 95% bootstrap CI |
| `di_corrected_ci_high` | Upper 95% bootstrap CI |
| `di_raw` | Raw direction instability |
| `di_raw_ci_low` | Lower 95% bootstrap CI |
| `di_raw_ci_high` | Upper 95% bootstrap CI |
| `magnitude_cv` | CV of signature norms across contexts |
| `n_contexts` | Number of contexts used |

---

## Analysis code

All scripts will be written and committed BEFORE execution. Script list:

1. `experiments/01_compute_tahoe_di.py` — DI on Tahoe pseudobulk (Modal)
2. `experiments/02_compute_jump_crispr_di.py` — DI on CRISPR profiles
3. `experiments/03_compute_jump_orf_di.py` — DI on ORF profiles
4. `experiments/04_compute_prism_viability_cv.py` — Viability CV (not DI)
5. `experiments/05_cross_dataset_concordance.py` — Exp A1 + A2
6. `experiments/06_prism_vs_di.py` — PRISM exploratory
7. `experiments/07_crispr_essentiality.py` — Exp C
8. `experiments/08_knockout_vs_overexpression.py` — Exp D

---

## False positive rate verification

All criteria were verified via simulation under the null (true rho = 0,
5,000 simulations). False positive rates:

| Experiment | n | FPR | Expected (alpha) |
|------------|---|-----|------------------|
| A1 | 103 | 0.0080 | 0.0125 |
| A2 | 330 | 0.0003 | 0.0125 |
| C | 7,948 | 0.0000 | 0.0125 |
| D | 5,220 | 0.0000 | 0.0125 |

All criteria are conservative (FPR well below alpha).

---

## Reporting commitment

Regardless of outcomes:
- Report all effect sizes with 95% CIs
- Report both raw and partial rho for all concordance tests
- Report permutation null distribution for Exp D
- Label exploratory analyses explicitly
- Interpret null results honestly (DI is modality-specific, not universal)
- Publish atlas with all columns regardless of experiment outcomes
