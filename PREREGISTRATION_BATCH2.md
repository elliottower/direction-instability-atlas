# Pre-registration: DI Atlas Batch 2A — Mechanistic Correlates

**Title:** What explains direction instability variation? Paralog buffering, expression landscape, and drug selectivity

**Authors:** Elliot Tower

**Date:** 2026-07-06 (pre SHA freeze)

**Commit SHA:** (to be frozen after script commit)

**Relation to Batch 1:** Batch 1 (PREREGISTRATION.md, SHA `0d07a01`)
established that DI is cross-modally concordant (A2: rho=0.19,
n=378), inversely correlated with essentiality (C: rho=-0.33,
n=7,627), mechanism-specific rather than gene-intrinsic (D: rho=0.01,
n=5,220), and correlated with viability variance (PRISM: rho=0.36,
n=75). Batch 2A asks why DI varies: what biological properties of a
perturbation target predict high or low DI?

---

## Overview

Three pre-registered experiments test whether DI variation across genes
and drugs is explained by (E2) cross-cell-line expression variance of
the target gene, (E3) selective vs. broad drug killing, and (E5)
paralog-mediated genetic buffering. All three extend Batch 1 results
using data already computed (crispr_di.csv, orf_di.csv, prism_cv.csv,
LINCS DI) plus one new download each.

---

## Multiple comparisons correction

Three primary tests. Bonferroni-corrected alpha = 0.05 / 3 = **0.0167**.

---

## Data sources

### Already computed (Batch 1)

| Dataset | n | File |
|---------|---|------|
| CRISPR DI | 7,946 genes | results/crispr_di.csv |
| ORF DI | 12,590 genes | results/orf_di.csv |
| PRISM CV | 1,433 drugs | results/prism_cv.csv |
| LINCS DI | 8,949 drugs | drug-perturbation-geometry/zenodo_v1/drug_instability_8949.csv |

### New downloads required

| Dataset | Source | File | Size |
|---------|--------|------|------|
| DepMap expression (TPM) | DepMap Public 24Q4 | OmicsExpressionProteinCodingGenesTPMLogp1.csv | ~200 MB |
| Ensembl paralog pairs | BioMart REST API (GRCh38, Ensembl 113) | Programmatic download in script | ~5 MB |
| PRISM viability matrix | Already on disk | data/prism/Repurposing_Public_24Q2_LFC_COLLAPSED.csv | On disk |

---

## Experiment E2: Expression variance as a Waddington basin proxy

### Hypothesis

Genes whose expression varies more across cell lines (shallow attractor
basins, context-dependent regulation) produce higher morphological DI
when knocked out. Conversely, genes with stable cross-cell-line
expression (deep universal basins, housekeeping) produce low DI.

### Biological rationale

Batch 1 Experiment C showed that broadly essential genes have LOW DI
(rho=-0.33). Essential genes tend to be ubiquitously expressed. This
experiment tests the more general version: expression variance across
cell lines — not just essentiality — predicts DI, because expression
variance indexes how context-dependent a gene's regulation is.

### Metric

For each gene g:
- **Expression variance:** Variance of log1p(TPM) across DepMap cell
  lines (using OmicsExpressionProteinCodingGenesTPMLogp1.csv). This is
  the cross-cell-line variance of the gene's own expression, NOT the
  expression of its targets.
- **CRISPR DI:** di_corrected from crispr_di.csv

### Confounder control

Partial Spearman rho controlling for n_contexts (number of JUMP-CP
CRISPR plate replicates per gene). Method: Pearson-residualize X and
Y against covariates Z (with intercept column), then compute Spearman
rho on the residuals. NOT rank-then-residualize, which is biased
(gives false partial rho 0.05-0.25 at n=5000 under the null; see
Batch 1 PREREGISTRATION.md). One covariate.

### Decision criteria (alpha = 0.0167)

BOTH conditions must hold:
1. 95% CI lower bound of partial Spearman rho > 0.05
2. p < 0.0167

CI computed via Fisher z-transform with df adjustment for 1 covariate.
The CI-lower-bound gate is essential: at n ~ 7,000, trivially small
effects (rho ~ 0.03) reach statistical significance. The 0.05 floor
ensures the effect is non-trivially positive.

### Power analysis

At n ~ 7,000 (estimated overlap between CRISPR DI genes and DepMap
expression genes), SE ~ 0.012. Any effect rho >= 0.08 is detected
with > 90% power. Same power regime as Batch 1 Experiment C.

### Interpretation

| rho | Interpretation |
|-----|----------------|
| >= 0.15 | Strong support: expression landscape position predicts DI |
| 0.05–0.15 | Weak but detectable: expression variance explains some DI variation |
| < 0.05 | Null: DI is not primarily driven by expression-level context-dependence |

### Supplementary

Report scatter plot of expression variance vs. DI, colored by
essentiality breadth, to visualize whether expression variance and
essentiality provide independent or redundant information about DI.

---

## Experiment E3: Drug selectivity — does DI predict selective killing?

### Hypothesis

High-DI drugs (transcriptomically context-dependent) are selective
killers: they kill a subset of cell lines strongly while sparing
others. Low-DI drugs kill broadly and uniformly (or not at all).

### Biological rationale

Batch 1 PRISM exploratory showed that DI correlates with viability CV
(rho=0.36, n=75). CV conflates selectivity with potency variance. This
experiment decomposes the signal: is DI associated with the *pattern*
of killing (selective vs. broad), not just the *variance* of killing?

### Magnitude confound warning

Gini coefficient on raw kill depth (negative LFC) is a variance-like
statistic on 1D viability. Drugs with larger mean kill magnitude
mechanically have more dynamic range to be selective over, producing
higher Gini regardless of true selectivity. This is the same
magnitude→variance confound documented in Part I of the combined
paper (DepMap d=2.65 → sign-reversal after CV normalization). The
primary metric must therefore be **magnitude-corrected Gini**, using
the same OLS residualization pattern as `di_corrected`:

    gini_corrected = gini_raw - OLS_predict(gini_raw ~ mean_kill_magnitude) + intercept

This parallels the atlas convention: `di_corrected` residualizes DI
against mean signature L2 norm; `gini_corrected` residualizes Gini
against mean kill depth.

### Metrics

For each drug in LINCS ∩ PRISM:
- **DI:** di_corrected from LINCS (magnitude-corrected direction
  instability)
- **Gini raw:** Gini coefficient of the drug's negative-LFC
  distribution across cell lines (using only cell lines where
  LFC < 0, i.e., growth inhibition). Gini = 0 means uniform killing;
  Gini = 1 means perfectly selective.
- **Gini corrected (PRIMARY):** Gini raw residualized against mean
  kill magnitude (mean of negative LFC values) via OLS, plus
  intercept. This is the primary selectivity metric.
- **Kill fraction (sensitivity analysis):** Fraction of cell lines
  with LFC < -0.5 (viability < ~60%).

### Universe

LINCS ∩ PRISM overlap by drug name. Expected n ~ 1,726 (from
PREREGISTRATION.md cross-reference). Restrict to drugs with
|mean_lfc| >= 0.01 (exclude inert compounds), matching Batch 1 PRISM
filter.

### Confounder control

Partial Spearman rho (Pearson-residualize X and Y against covariates,
then Spearman on residuals — same method as Batch 1) controlling for
LINCS n_contexts. One covariate. Mean kill magnitude is already
removed from Gini via the correction above, so it is NOT a second
covariate (double-correcting would attenuate real signal).

### Decision criteria (alpha = 0.0167)

95% CI lower bound of partial Spearman rho > 0.05 AND p < 0.0167.
CI computed via Fisher z-transform with df adjustment for 1 covariate.

Primary metric: gini_corrected. Gini raw and kill fraction reported
as sensitivity analyses.

### Power analysis

At n ~ 1,000+ (conservative estimate after filtering), SE ~ 0.032.
Any effect rho >= 0.12 detected with > 80% power.

### Interpretation

| rho | Interpretation |
|-----|----------------|
| >= 0.20 | Strong: DI predicts drug selectivity (precision oncology hook) |
| 0.05–0.20 | Moderate: some selectivity signal, may be driven by specific drug classes |
| < 0.05 | Null: DI does not predict selective vs. broad killing |

### Supplementary

Report top 20 high-DI selective killers and top 20 low-DI broad
killers with known MOA annotations (from LINCS pert_info).

---

## Experiment E5: Paralog buffering predicts CRISPR DI

### Hypothesis

Genes with expressed paralogs in a subset of cell lines show higher
CRISPR knockout DI than genes without paralogs. Paralog compensation
buffers loss-of-function in some cellular contexts but not others,
producing context-dependent phenotypes (high DI).

### Biological rationale

Batch 1 Experiment D showed that KO-DI and OE-DI are uncorrelated
(rho=0.01). The exploratory analysis identified 289 genes that are
KO-unstable but OE-stable, consistent with paralog-buffered genes.
This experiment tests the paralog hypothesis directly.

### Metrics

For each gene g:
- **CRISPR DI:** di_corrected from crispr_di.csv
- **Paralog count:** Number of paralogs from Ensembl BioMart
  (hsapiens_gene_ensembl, attribute: hsapiens_paralog_ensembl_gene)
- **Has expressed paralog (binary):** Whether any paralog has mean
  TPM > 1 across DepMap cell lines (requires expression data from E2)
- **Paralog expression breadth:** For the most highly expressed
  paralog, the fraction of DepMap cell lines where it is expressed
  (TPM > 1). High breadth = universal compensation = expected LOW DI.
  Low breadth = selective compensation = expected HIGH DI.

### Tests

**Primary (E5a):** Partial Spearman rho between paralog count and
CRISPR DI, controlling for n_contexts and mean expression level of
the gene itself (mean TPM across DepMap cell lines). Two covariates.
Method: Pearson-residualize X and Y against covariates Z (with
intercept column), then Spearman on residuals. NOT
rank-then-residualize (biased; see Batch 1 PREREGISTRATION.md).

**Secondary (E5b):** Among genes WITH at least one expressed paralog
(TPM > 1 in any DepMap cell line), partial Spearman rho between
paralog expression breadth and CRISPR DI, controlling for n_contexts
and mean expression level. Same Pearson-residualize-then-Spearman
method. Expected sign: NEGATIVE (broader paralog expression = more
uniform compensation = lower DI).

### Decision criteria (alpha = 0.0167)

E5a is the primary test. BOTH conditions must hold:
1. 95% CI lower bound of partial Spearman rho > 0.05
2. p < 0.0167

CI via Fisher z-transform with df adjustment for 2 covariates. The
CI-lower-bound gate is essential at n ~ 7,000 (see E2 rationale).

E5b is secondary (same alpha, reported separately). Expected sign is
negative, so the criterion is:
1. CI upper bound < -0.05
2. p < 0.0167

### Power analysis

At n ~ 7,000 genes, same power as E2. Any effect rho >= 0.08 detected
with > 90% power.

### Interpretation

| E5a rho | E5b rho | Interpretation |
|---------|---------|----------------|
| > 0.05, significant | < -0.05, significant | Full support: paralogs buffer KO effects context-dependently |
| > 0.05, significant | n.s. | Partial: having paralogs matters, but buffering is not context-specific |
| n.s. | n.s. | Null: paralog buffering does not explain DI variation |

---

## Analysis code

Scripts will be written and committed BEFORE execution:

| Script | Experiment | Inputs |
|--------|-----------|--------|
| `experiments/09_expression_variance.py` | E2 | crispr_di.csv + DepMap expression |
| `experiments/10_prism_selectivity.py` | E3 | LINCS DI + PRISM LFC matrix |
| `experiments/11_paralog_buffering.py` | E5 | crispr_di.csv + BioMart paralogs + DepMap expression |

All scripts use `experiments/utils.py` (partial_spearman, ci_lower_bound,
fisher_z_pvalue, save_results) from Batch 1.

---

## Outputs

| File | Contents |
|------|----------|
| `results/expression_variance.json` | E2 result + scatter data |
| `results/prism_selectivity.json` | E3 result (Gini + kill fraction) |
| `results/paralog_buffering.json` | E5a + E5b results |

---

## Relationship to concurrent work

Raju (arXiv:2603.00678, arXiv:2604.16642) introduced Shesha, a
within-context coherence metric measuring cell-to-cell directional
agreement under perturbation. DI measures cross-context coherence.
These are complementary scales. Batch 2A does not test the DI-Shesha
relationship (that requires single-cell Perturb-seq data and is
deferred to Batch 2B). However, E2 (expression variance as basin
depth proxy) is motivated by the same Waddington landscape intuition
that Raju's papers invoke.

### Why E1 (Shesha vs DI) is excluded from this freeze

E1 (two-scale geometry) requires single-cell Perturb-seq data
(Replogle et al.), which is a different data modality requiring a
new compute pipeline. If E1 were folded into Batch 2A, a later
addition would retroactively change the multiple-comparisons family
(3 tests → 4 tests), invalidating the Bonferroni correction. E1
will receive its own pre-registration freeze when the Perturb-seq
pipeline exists (Batch 2B).

### Verified Raju references

- arXiv:2603.00678 — "From Syntax to Semantics: Geometric Stability
  as the Missing Axis of Perturbation Biology" (Feb 28, 2026)
- arXiv:2604.16642 — "Geometric coherence of single-cell CRISPR
  perturbations reveals regulatory architecture and predicts cellular
  stress" (Apr 16, 2026)
- Code: github.com/prashantcraju/geometric-stability-crispr
- All three verified via arXiv on 2026-07-06.

---

## Reporting commitment

Same as Batch 1: report all effect sizes with 95% CIs, report both
raw and partial rho, interpret null results honestly. All three
experiments are reported regardless of outcome.
