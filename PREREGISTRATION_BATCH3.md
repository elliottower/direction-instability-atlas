# Pre-registration: DI Atlas Batch 3 — Reviewer-Prompted Extensions

**Title:** MoA structure, cross-cell-line pilot, and convergent essentiality evidence

**Authors:** Elliot Tower

**Date:** 2026-07-09 (pre SHA freeze)

**Commit SHA:** `4c34699`

**Provenance:** These experiments were designed in response to peer
review of an earlier draft (desk-rejected from GigaScience for
scope/breadth, 2026-07-08). They are pre-registered before execution
but are *reviewer-prompted, post-hoc additions*, not part of the
original prospective plan. Timestamps can be verified against the
repository history. This distinction is disclosed in the manuscript.

**Relation to earlier batches:** Batch 1 (SHA `0d07a01`) established
cross-modal concordance, the essentiality sign-flip, mechanism
specificity, and PRISM viability concordance. Batch 2A (SHA `a17c125`)
tested expression variance, drug selectivity, and paralog buffering.
Batch 3 addresses three reviewer concerns: (F) demonstrated utility
of the atlas, (G) validation of the essentiality sign-flip on true
cross-cell-line morphological data, and (H) convergent evidence from
a different data modality.

---

## Multiple comparisons correction

Three primary tests. Bonferroni-corrected alpha = 0.05 / 3 = **0.0167**.

---

## Experiment F: MoA-stratified DI characterization

### Question

Does direction instability vary systematically across drug mechanism-
of-action classes, after controlling for profiling depth and signal
strength?

### Why this matters

Batch 1 and 2A showed that DI is not explained by gene-level
annotations (essentiality, expression variance, paralogs). Reviewer
concern: if DI is not explained by anything, is it just noise? This
experiment tests whether DI is *mechanism-structured* — whether drugs
sharing a pharmacological mechanism have similar DI — which would
demonstrate that DI captures biologically meaningful variation and
that the atlas is useful for stratifying drug classes.

### Framing and circularity disclosure

DI measures cross-cell-line signature consistency. Drugs sharing an
MoA may share similar targets and pathways, which could produce
similar cross-cell-line behavior and therefore similar DI. This means
"MoA explains DI variance" may partly reflect the structural fact
that drugs with similar mechanisms behave similarly. We frame the
result as **characterization** ("DI is mechanism-structured: core-
machinery drugs show low DI while receptor-signaling drugs show high
DI") rather than **prediction** ("MoA predicts DI"). The value is in
the interpretable biological pattern, not in the R² itself.

### Data

- LINCS DI: drug_instability_8949.csv (8,949 drugs)
- MoA annotations: `moa` column in LINCS DI CSV
- Restrict to MoA classes with >= 5 drugs (estimated ~80 classes,
  ~1,031 drugs based on pre-analysis of the MoA column)

### Confound controls (MANDATORY)

DI depends mechanically on the number of cell lines a drug was
profiled in, and MoA classes may differ systematically in profiling
depth. The primary analysis therefore uses **residualized DI**:

1. Compute `di_residualized` = residuals of OLS(di_corrected ~
   n_cell_lines + mean_norm) + intercept. This removes both the
   profiling-depth confound and any residual magnitude effect.
2. Compute one-way ANOVA of `di_residualized` on MoA class.
3. Report partial eta-squared (η²_partial) on residualized DI.

### Permutation null (MANDATORY)

One-way ANOVA η² is inflated by class-size imbalance (many small
groups with 5-6 drugs can diverge by chance). The permutation null
controls for this:

1. Shuffle MoA labels among the 1,031 drugs (preserving class sizes).
2. Recompute η² on the shuffled data.
3. Repeat 1,000 times.
4. Report: observed η², null distribution (mean, SD, 95th percentile),
   empirical p-value (fraction of null >= observed).

### Decision criteria (alpha = 0.0167)

BOTH conditions must hold:
1. Empirical permutation p < 0.0167
2. Observed partial η² on residualized DI exceeds the 99th percentile
   of the permutation null

The R² is NOT the primary result — the pattern of which MoA classes
are high vs. low DI is. The η² test establishes that the pattern is
non-random; the biological interpretation is the contribution.

### Reporting

- Table of top 10 and bottom 10 MoA classes by mean residualized DI
  (with SD and n per class)
- Raw η² vs. residualized η² (to quantify confound contribution)
- Biological interpretation: do low-DI MoAs target fundamental
  cellular machinery (expected: HDAC, topoisomerase, kinase) while
  high-DI MoAs target receptor-mediated signaling (expected: GPCR,
  nuclear receptor)?

### Interpretation guide

| Outcome | Interpretation |
|---------|----------------|
| η² significant, biologically interpretable | DI is mechanism-structured; atlas stratifies drug classes |
| η² significant, not interpretable | MoA structure exists but lacks clear biological pattern |
| η² not significant after confound control | MoA structure in raw ANOVA was confound-driven |

---

## Experiment G: cpg0000 pilot cross-cell-line CRISPR

### Question

Does the essentiality sign-flip (Batch 1 Experiment C: ρ = −0.33)
replicate in direction when DI is computed across cell lines rather
than plate replicates?

### Why this matters

Batch 1 Experiment C computed CRISPR DI over plate replicates within
a single cell line (U2OS, source_13 in cpg0016). Plate replicates
measure technical reproducibility, not biological context-dependence.
Reviewers asked whether the sign-flip reflects biology or a property
of plate-level noise. The cpg0000 pilot dataset contains CRISPR
profiles from TWO cell lines (A549 and U2OS), providing the only
available true cross-cell-line morphological CRISPR data in the JUMP
ecosystem.

### Data

- cpg0000 pilot: Cell Painting Gallery (AWS S3, cpg0000-jump-pilot)
- CRISPR profiles from A549 and U2OS (~160 genes with profiles in
  both cell lines)
- DepMap Chronos 24Q4: essentiality breadth per gene (already on disk)

### DI computation with K=2

With only 2 cell lines, DI degenerates to:

    DI = 1 - cos(s_A549, s_U2OS)

This is a single pairwise cosine distance, not an average over
multiple pairs. It is geometrically valid but extremely noisy per
gene (no averaging over contexts). The Batch 1 protocol requires
K >= 5; this experiment explicitly relaxes that requirement and
reports the limitation.

No magnitude correction is applied (OLS on n=~160 with 1 feature
would overfit). Raw DI is the primary metric.

No bootstrap CIs (K=2 admits only one bootstrap resample pattern).

### Statistical test

Spearman correlation between cross-cell-line DI and DepMap
essentiality breadth (fraction of cell lines with Chronos < -0.5).

No covariate control (n_contexts = 2 for all genes; no variation to
partial out).

### Decision criteria

This is a **directional convergent check**, not a confirmatory test.
At n ≈ 160, the study has ~80% power for |ρ| >= 0.22 and <50% power
for |ρ| <= 0.15. The Batch 1 result was ρ = −0.33.

Decision:
- If ρ is negative and p < 0.0167: directional replication
- If ρ is negative but p >= 0.0167: consistent direction, underpowered
- If ρ is positive or near zero: does not replicate; plate-replicate
  interpretation weakened

Report: Spearman ρ, 95% CI, p-value, n.

### Power analysis

| True rho | n | Power (one-sided, alpha=0.0167) |
|----------|---|------|
| -0.15 | 160 | 0.32 |
| -0.20 | 160 | 0.53 |
| -0.25 | 160 | 0.73 |
| -0.30 | 160 | 0.87 |
| -0.35 | 160 | 0.95 |

### Limitations (stated in advance)

1. K=2 cell lines makes per-gene DI extremely noisy
2. n ≈ 160 genes is underpowered vs. Batch 1 (n = 7,653)
3. Feature space may differ between cpg0000 and cpg0016 (pilot
   vs. production preprocessing)
4. Gene overlap between cpg0000 CRISPR and DepMap may be smaller
   than 160

These are disclosed pre-analysis. The experiment's value is as a
directional check, not a definitive test.

---

## Experiment H: DepMap essentiality variability as convergent evidence

### Question

Does the essentiality sign-flip generalize to a different data
modality? Specifically: do broadly essential genes show lower
cross-cell-line variability in dependency scores than selectively
essential genes?

### Why this matters

The Batch 1 sign-flip was computed on morphological DI (JUMP-CP
CRISPR plate replicates). The biological interpretation —
"functional constraint forces stereotyped disruption" — predicts
that the same pattern should appear in dependency scores: genes
essential in many cell types should show uniform dependency across
those cells (low variability), while selectively essential genes
should show variable dependency (high variability). Testing this on
DepMap Chronos provides convergent evidence from a fundamentally
different measurement.

### Estimand disclosure (CRITICAL)

DepMap Chronos dependency scores are a DIFFERENT geometric construct
from morphological DI:

- Morphological DI: mean pairwise cosine distance between
  high-dimensional Cell Painting signature vectors across contexts.
  Each context produces a vector in R^d.
- DepMap variability: coefficient of variation of scalar dependency
  scores across cell lines. Each cell line produces a single number.

These measure different things. A positive result is **convergent
evidence** for the shared biological hypothesis (functional constraint
suppresses context-dependence), NOT validation of morphological DI.
The manuscript must frame it this way.

### Data

- DepMap Chronos 24Q4: CRISPRGeneEffect.csv (1,178 cell lines,
  17,916 genes; already on disk)
- Essentiality breadth: fraction of cell lines with Chronos < -0.5
  (same definition as Batch 1 Experiment C)

### Metric

For each gene g:
- **Essentiality variability:** SD of Chronos dependency scores
  across cell lines (SD, not CV, because Chronos scores can be near
  zero or negative, making CV undefined)
- **Essentiality breadth:** fraction of cell lines with Chronos < -0.5

### Confound control

Partial Spearman rho between essentiality variability and essentiality
breadth, controlling for the number of non-missing Chronos values
per gene (some genes have missing data in some cell lines). Method:
Pearson-residualize then Spearman (same as all prior experiments).

### Decision criteria (alpha = 0.0167)

The biological prediction is NEGATIVE: broadly essential genes
(high breadth) should have LOW variability (uniform dependency).

BOTH conditions must hold:
1. 95% CI upper bound of partial Spearman rho < -0.05
2. p < 0.0167

### Interpretation

| Outcome | Interpretation |
|---------|----------------|
| Negative ρ, significant | Convergent: functional constraint suppresses context-dependence in dependency scores as in morphological profiles |
| Null or positive | Does not converge; morphological sign-flip may be modality-specific |

### Why this is not validation

Even a positive result here does not "validate" the morphological DI
sign-flip, because:
1. The geometric construct is different (scalar SD vs. vector cosine)
2. The data modality is different (dependency vs. morphology)
3. The contexts partially overlap (DepMap cell lines ⊃ JUMP-CP U2OS)

It tests the *biological hypothesis* (functional constraint suppresses
context-dependent phenotypes), not the *metric* (morphological DI).

---

## Analysis code

Scripts will be written and committed BEFORE execution:

| Script | Experiment | Inputs |
|--------|-----------|--------|
| `experiments/12_moa_stratification.py` | F | LINCS DI CSV |
| `experiments/13_cpg0000_pilot.py` | G | cpg0000 S3 data + DepMap |
| `experiments/14_depmap_convergent.py` | H | DepMap Chronos |

All scripts use `experiments/utils.py` from Batches 1-2.

---

## Outputs

| File | Contents |
|------|----------|
| `results/moa_stratification.json` | Exp F: η², permutation null, class table |
| `results/cpg0000_pilot.json` | Exp G: cross-cell-line essentiality ρ |
| `results/depmap_convergent.json` | Exp H: essentiality variability ρ |

---

## Reporting commitment

Same as Batches 1-2: report all effect sizes with 95% CIs, report
both raw and confound-controlled results, interpret null results
honestly. All three experiments are reported regardless of outcome.
Reviewer-prompted provenance is disclosed.
