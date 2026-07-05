# DI Atlas Paper Plan — V3 (HISTORICAL)

> **Superseded by PREREGISTRATION.md** for all decision criteria,
> alpha values, and power analyses. This document preserves the V3
> plan as reviewed by Perplexity. Key differences: PREREGISTRATION.md
> drops Experiment B (structurally underpowered), drops A3 (zero
> overlap), corrects alpha to 0.0125 (4 tests), and uses the
> corrected partial Spearman method (Pearson-residualize then
> Spearman, not rank-then-residualize).

V2 fixes (preserved):
1. PRISM DI is undefined (viability is 1D, "direction" is meaningless)
2. Experiment A threshold too weak (rho > 0.20 explains 4% of variance)
3. Raw DI as headline column propagates the confound Paper 1 warns about

V3 fixes (from Perplexity review of V2):
4. All concordance experiments need partial-correlation control for n_contexts
5. Exp C threshold too weak (rho > 0.10 is significant-but-trivial trap)
6. Exp D needs permutation null (shared plate structure inflates rho)
7. Multiple comparisons correction across pre-registered experiments
8. Atlas needs per-perturbation bootstrap CIs, not just point estimates
9. Tahoe pseudobulk needs minimum-cell threshold
10. Exp D needs directional divergence analysis (knockout-stable/OE-unstable vs reverse)

---

## Thesis (revised)

Direction instability is concordant across transcriptomic and morphological modalities
and generalizes from chemical to genetic perturbations. Magnitude-corrected DI
scores are published as a community resource for filtering perturbation experiments.

Three DI-proper modalities (LINCS transcriptomics, Tahoe transcriptomics, JUMP-CP
morphology), plus PRISM viability as a complementary (non-DI) validation.

This is NOT the combined confound paper (which proves DI is the right metric).
This paper assumes that result and says: "here's what DI looks like on everything."

## Repo

New repo: `direction-instability-atlas`.
Imports `geometry/bracket_norm.py` from bracket-norm-validity for DI computation.

---

## Datasets

### DI-proper datasets (high-dimensional signatures, cosine is meaningful)

1. **LINCS L1000** — 8,949 drugs, DI across cell lines (978-gene transcriptomics)
   - Already computed in bracket-norm-validity.

2. **JUMP-CP compounds** — 25,254 compounds, DI across source contexts (morphology)
   - Already computed in bracket-norm-validity.

3. **Tahoe-100M** — 1,100 drugs × 50 cancer cell lines (full-transcriptome single-cell)
   - Source: https://huggingface.co/datasets/tahoebio/Tahoe-100M
   - Pseudobulk differential expression per drug per cell line → DI across cell lines.
   - Directly comparable to LINCS (both transcriptomics, both drug perturbations).
   - Full transcriptome (~20k genes) vs LINCS (978 genes) — can also restrict to
     the 978 LINCS landmarks for a matched comparison.
   - **Key test:** LINCS ∩ Tahoe concordance on shared drugs. Same readout type,
     different platform, different cell line panels. Should be the highest
     cross-dataset correlation.

4. **JUMP-CP CRISPR knockouts** — 7,975 genes, morphological profiles
   - Same JUMP-CP dataset, genetic perturbation type.
   - DI across source contexts (plates/batches), same computation as compounds.

5. **JUMP-CP ORF overexpression** — 12,602 genes, morphological profiles
   - Gain-of-function complement to CRISPR loss-of-function.

### Complementary dataset (NOT DI-proper)

6. **PRISM Repurposing** — ~4,518 compounds × 578 cell lines (viability)
   - Viability is 1-dimensional per cell line. Cosine between scalars is undefined.
   - DI CANNOT be computed on PRISM. This is not a limitation of the analysis;
     it is a property of the data. Directional metrics require high-dimensional
     signatures.
   - **What we CAN compute:** profile variability (CV of viability across cell lines)
     or, if the secondary screen provides 8-dose curves, dose-response shape
     variability (CV of AUC, or Frechet variance of normalized dose curves).
   - **Role in atlas:** sanity check. Do drugs with low DI (in LINCS/JUMP-CP) also
     show low viability variability in PRISM? This tests whether mechanism
     conservation (DI) predicts phenotypic consistency (viability CV) — a weaker
     but distinct claim from cross-modal DI concordance.
   - **Explicitly not included in Experiment A concordance test.**

---

## Pre-registration experiments

### Experiment A: Cross-dataset DI concordance (revised)

**Hypothesis:** Drugs present in multiple DI-proper datasets will show correlated
DI across datasets.

Three concordance pairs, tiered by expected agreement:

| Pair | Modality match | Expected strength |
|------|---------------|-------------------|
| LINCS ∩ Tahoe | same (transcriptomics) | highest |
| LINCS ∩ JUMP-CP | cross (transcriptomics vs morphology) | moderate |
| Tahoe ∩ JUMP-CP | cross (transcriptomics vs morphology) | moderate |

**Confounder control (MANDATORY — applies to all concordance tests):**
Two datasets with correlated per-drug sample sizes (n_contexts) will produce
correlated DI even if DI is pure noise, because both estimates are noisy in the
same cell-count-dependent way. This is the exact confound Paper 1 characterizes.
For every concordance pair, report BOTH:
- Raw Spearman rho (bivariate)
- **Partial Spearman rho controlling for n_contexts in both datasets**
  (rank-based partial correlation, conditioning on BOTH n_A and n_B as two covariates;
  min() throws away information and is less defensible)

The partial rho is the primary result. Raw rho is reported for transparency.
If the partial rho drops below the threshold while raw rho passes, the
concordance is confound-driven, not biology-driven.

**Decision criteria (tiered, applied to PARTIAL rho):**
- **Same-modality (LINCS ∩ Tahoe):** Partial Spearman rho > 0.40.
  Both are transcriptomics; if DI is a property of the drug and not the cell
  line panel, these should correlate well after controlling for sample size.
  Partial rho < 0.40 would mean platform/cell-line-panel differences dominate
  over drug-intrinsic mechanism stability.
- **Cross-modality (LINCS ∩ JUMP-CP, Tahoe ∩ JUMP-CP):** Partial Spearman rho > 0.20.
  Transcriptomic and morphological instability need not agree perfectly —
  a drug could change gene expression consistently but morphology inconsistently,
  or vice versa. Partial rho > 0.20 supports the claim that mechanism stability
  has a modality-general component; partial rho < 0.20 means it's modality-specific.

**Interpretation guide:**
- All three partial rho pass → DI is a property of the drug. Strong atlas thesis.
- Same-modality passes, cross-modality fails → DI is modality-specific.
  Atlas is still useful (per-modality resource) but concordance claim weakens.
  This is framed positively: "DI is modality-specific, and the right assay
  depends on the question" is itself a useful resource conclusion.
- Same-modality fails → DI is platform-specific. Would seriously weaken the atlas.
- Raw rho passes but partial rho fails → concordance is confound-driven.
  Would require honest discussion of n_contexts as the dominant signal.

**PRISM excluded from this experiment.** PRISM viability CV is tested separately as
a complementary sanity check (not a DI concordance test).

### Experiment A2: PRISM viability CV vs DI (complementary, not primary)

**Hypothesis:** Drugs with low DI in LINCS will show low viability CV in PRISM.

**Decision criteria:** Spearman rho > 0.15 between LINCS DI and PRISM viability CV
for shared drugs. This is a weak bar because the claim is weak: mechanism
conservation predicts viability consistency, but imperfectly.

### Experiment B: Machinery/receptor stratification across new modalities

**Hypothesis:** The machinery < receptor separation (validated in LINCS and JUMP-CP
in the combined paper) extends to Tahoe.

**Decision criteria:** Cohen's d > 0.3 (machinery vs receptor) in Tahoe.

Note: PRISM is NOT included because we're computing viability CV, not DI.
The machinery/receptor distinction is about mechanism direction conservation,
which viability can't measure. We CAN do an exploratory check: do machinery
drugs have lower viability CV than receptor drugs? But this is labeled
EXPLORATORY, not pre-registered.

### Experiment C: JUMP-CP CRISPR DI predicts essentiality (SCIENTIFIC HEART)

**Hypothesis:** In JUMP-CP CRISPR knockouts, genes essential in more cell types
(per DepMap Chronos) will have higher morphological DI.

This is the cleanest "new results on Broad data" story. Nobody has asked whether
essential gene knockouts produce more context-dependent morphological changes.
The combined paper showed this in viability (DepMap) and transcriptomics
(Perturb-seq). Extending to morphology is novel.

**Confounder control:** Same as Experiment A — report partial Spearman rho
controlling for n_contexts (number of source plates/batches per gene in JUMP-CP).
Genes with more replicates produce better-estimated DI, and essentiality status
may correlate with experimental coverage.

**Decision criteria:** Partial Spearman rho > 0.20 between DepMap essentiality
breadth (fraction of cell lines with Chronos < -0.5) and JUMP-CP CRISPR DI,
controlling for n_contexts.

**Interpretation of weak effects:** If partial rho falls in [0.10, 0.20],
report as "weak/negligible — statistically significant due to large N but
explaining < 4% of variance." Do NOT lean on p-value to claim confirmation.
The combined paper exists to criticize exactly this reasoning (significant
effect that's trivially small). rho < 0.10 is reported as non-replication.

### Experiment D: Knockout vs overexpression DI correlation (SCIENTIFIC HEART)

**Hypothesis:** Genes with high DI in CRISPR knockout will also have high DI in ORF
overexpression (same gene, opposite perturbation direction, same morphological readout).

This is the most novel experiment in the paper. If knockout and overexpression
of the same gene produce similarly (in)consistent morphological responses across
contexts, that says something about the gene's role in cell biology — not about
the perturbation direction. Genes that are context-dependent are context-dependent
regardless of whether you remove or add them.

**Permutation null (MANDATORY):** CRISPR and ORF perturbations are profiled on
the same JUMP-CP plates and batches. Shared technical structure (well-position
effects, batch artifacts, cell-density confounds) will inflate rho regardless
of biology. Pre-registered null: shuffle gene labels in the ORF dataset
(breaking the gene-level pairing while preserving the plate structure) and
recompute rho 1,000 times. Report:
- Observed rho
- Permutation null distribution (mean, SD, 95th percentile)
- Empirical p-value (fraction of permutations >= observed rho)

The observed rho must exceed the 95th percentile of the null to pass.

**Decision criteria:** Spearman rho > 0.15 for shared genes AND observed rho
exceeds the 95th percentile of the permutation null. Both must hold.

**Directional divergence analysis (pre-registered exploratory):**
For genes where knockout DI and overexpression DI diverge substantially
(|DI_CRISPR - DI_ORF| > 1 SD), classify into:
- **Knockout-stable, overexpression-unstable:** gene removal has conserved
  morphological effects, but gene addition is context-dependent.
  Predicted biology: dosage-sensitive gain-of-function genes.
- **Knockout-unstable, overexpression-stable:** gene removal produces
  context-dependent compensatory responses, but gene addition has conserved
  effects. Predicted biology: genes with redundant paralogs.

Report the counts and enrichment of each class in Gene Ontology categories.
This analysis is labeled EXPLORATORY — no pass/fail threshold. The sign of
divergence is where the memorable biology lives.

---

## Statistical framework

### Multiple comparisons correction

Four pre-registered experiments (A, B, C, D) with the following primary criteria:

| Exp | Primary criterion | Test |
|-----|------------------|------|
| A | Partial rho (3 pairs, tiered thresholds) | 3 tests |
| B | Cohen's d > 0.3 in Tahoe | 1 test |
| C | Partial rho > 0.20 | 1 test |
| D | rho > 0.15 AND exceeds permutation null | 1 test |

Total primary tests: 6. Apply Bonferroni correction: alpha = 0.05/6 = 0.0083
for each individual test. Note: Bonferroni is conservative here because the
three Exp A pairs share the LINCS DI vector (not fully independent). We err
on the strict side deliberately. Pre-register that:
- All primary criteria use Bonferroni-corrected alpha = 0.0083
- Experiments A2 (PRISM complementary) and Exp D directional divergence
  are labeled EXPLORATORY and reported as hypothesis-generating only,
  with uncorrected p-values and explicit "exploratory" labels
- Effect sizes and CIs are reported for ALL tests regardless of significance

### Tahoe pseudobulk minimum-cell threshold

Given that the combined paper (Exp 1) demonstrated cell-count-dependent DI
estimation even for cosine-based metrics (rho = -0.16 with cell count),
pre-register:
- **Minimum cells per pseudobulk: >= 20 cells** (matching Perturb-seq threshold)
- Drug-cell-line combinations with < 20 cells excluded before DI computation
- Report DI as a function of cells-per-pseudobulk in a supplementary figure
  to confirm flat-bin immunity carries over from Perturb-seq to Tahoe

---

## Atlas columns (the resource)

For every perturbation in every dataset, the atlas publishes:

| Column | Description |
|--------|-------------|
| `di_corrected` | **PRIMARY.** Magnitude-corrected DI: DI residualized against mean signature L2 norm. This is the column users should default to. |
| `di_corrected_ci_low` | Lower bound of 95% bootstrap CI (resampling over contexts, 1000 iterations). |
| `di_corrected_ci_high` | Upper bound of 95% bootstrap CI. |
| `di_raw` | Raw direction instability (Eq. 1 from combined paper). Provided for transparency and reproducibility. |
| `di_raw_ci_low` | Lower bound of 95% bootstrap CI for raw DI. |
| `di_raw_ci_high` | Upper bound of 95% bootstrap CI for raw DI. |
| `magnitude_cv` | Coefficient of variation of signature norms across contexts. |
| `n_contexts` | Number of cell lines / source contexts used to compute DI. Users should filter by n_contexts >= 10 for reliable estimates. |

**Per-entry uncertainty rationale:** A drug with DI = 0.3 from 10 cell lines and
one with DI = 0.3 from 40 cell lines have very different reliability. If people
build on this atlas, they'll rank compounds by DI without knowing which rankings
are stable. Bootstrap CIs (resampling the set of contexts and recomputing DI)
let users filter by CI width or check whether two perturbations' CIs overlap.
No existing perturbation resource ships per-entry uncertainty — this is a
differentiator and the kind of rigor that elevates a resource paper.

**Corrected-DI-as-primary rationale:** The combined paper (Exp 2) showed that
magnitude correction drops DI's predictive AUROC from 0.945 to 0.754. That 0.19
gap is magnitude-driven signal. Publishing raw DI as the headline column would
hand users a metric that's ~20% magnitude-contaminated — exactly the confound
Paper 1 warns about. Leading with corrected DI is consistent with Paper 1's
message.

---

## Figures (revised)

1. **Panel: DI distributions across 5 DI-proper datasets + PRISM CV**
   Violin plots. DI datasets on same y-axis; PRISM CV on separate axis with
   explicit annotation "viability CV, not DI."

2. **Cross-dataset scatter: LINCS vs Tahoe DI** (same-modality, primary concordance)
   Annotate Spearman rho. Color by MOA class if available.

3. **Cross-dataset scatter: LINCS vs JUMP-CP DI** (cross-modality)
   Same format as Fig 2.

4. **Experiment C: CRISPR DI vs essentiality breadth**
   Scatter plot, JUMP-CP CRISPR DI vs fraction of DepMap cell lines where
   Chronos < -0.5. Color by pathway. This is the "new results on Broad data" figure.

5. **Experiment D: Knockout vs overexpression DI**
   Scatter plot, per-gene DI in CRISPR vs DI in ORF. Annotate rho.
   Highlight genes where knockout and overexpression DI diverge (interesting biology).

6. **MOA stratification panel** — machinery vs receptor in LINCS, JUMP-CP, Tahoe.
   Three columns. PRISM as exploratory fourth column with "viability CV" label.

7. **Heatmap: drugs × datasets** — for drugs present in 3+ datasets,
   heatmap of corrected DI. Rows ordered by mean corrected DI.

---

## Pre-commit work (before pre-registration freeze)

Same as V1 but with revised script list:

### Phase 0: Repo setup (30 min)
- [ ] Create `direction-instability-atlas` repo
- [ ] Copy `geometry/bracket_norm.py` (or pip-installable shared dependency)
- [ ] Set up `pyproject.toml` with uv
- [ ] Create `data/`, `results/`, `experiments/` structure

### Phase 1: Data download + metadata inspection (NO DI computation)
- [ ] Download PRISM 24Q2 (Figshare) — inspect format, confirm 1D viability
- [ ] Download Tahoe-100M pseudobulk tables (HuggingFace)
- [ ] Download JUMP-CP CRISPR + ORF profiles
- [ ] Inspect data formats, column names, dimensions
- [ ] Document in DATA_INVENTORY.md

### Phase 2: Cross-reference / matching (metadata only, no DI)
- [ ] Match drugs: LINCS ↔ Tahoe ↔ JUMP-CP ↔ PRISM by InChIKey
- [ ] Count overlap sizes (how many drugs in 2, 3, 4 datasets)
- [ ] Match JUMP-CP CRISPR genes ↔ ORF genes (shared gene set)
- [ ] Match JUMP-CP CRISPR genes ↔ DepMap essentiality labels
- [ ] Document all matching in CROSS_REFERENCE.md

### Phase 3: Write analysis scripts (frozen before execution)
- [ ] `experiments/01_compute_tahoe_di.py`
- [ ] `experiments/02_compute_jump_crispr_di.py`
- [ ] `experiments/03_compute_jump_orf_di.py`
- [ ] `experiments/04_compute_prism_viability_cv.py` (NOT DI)
- [ ] `experiments/05_cross_dataset_concordance.py` (Exp A, DI-proper only)
- [ ] `experiments/06_prism_vs_di.py` (Exp A2, complementary)
- [ ] `experiments/07_moa_stratification_tahoe.py` (Exp B)
- [ ] `experiments/08_crispr_essentiality.py` (Exp C)
- [ ] `experiments/09_knockout_vs_overexpression.py` (Exp D)

### Phase 4: Pre-registration freeze
- [ ] Write PREREGISTRATION.md
- [ ] Commit everything (scripts + prereg + data inventory + cross-reference)
- [ ] Record SHA
- [ ] THEN run experiments

---

## Scientific heart of the paper (what leads)

Per Perplexity's review: Experiments C and D are the strongest, most novel
contributions. The concordance stuff (Exp A) is breadth — it shows DI is
consistent. But C and D are *findings*:

- **Exp C:** Essential gene knockouts produce more context-dependent morphological
  responses. This is biology. Nobody has shown this in Cell Painting data.
- **Exp D:** Knockout and overexpression DI are correlated. Context-dependence is
  a property of the gene, not the perturbation direction. This is novel.

The paper should lead with C and D as the scientific contribution, with
concordance (A) and stratification (B) as supporting evidence that the atlas
is a reliable resource.

---

## Timeline estimate

Same as V1: ~2-3 weeks to submission-ready draft.

---

## What changed from V1 → V2

1. **PRISM demoted from DI-proper to complementary.** Viability is 1D; cosine is
   undefined. PRISM contributes viability CV, not DI. Excluded from Experiment A.
   Thesis revised from "four modalities" to "three modalities + complementary."

2. **Experiment A thresholds raised.** Same-modality (LINCS∩Tahoe) bar raised to
   rho > 0.40. Cross-modality kept at rho > 0.20 but framed as the weaker claim
   it is. PRISM split into separate Exp A2.

3. **Corrected DI is the primary atlas column.** Raw DI is secondary. Consistent
   with Paper 1's message that ~20% of DI's signal is magnitude-driven.

4. **C and D elevated to scientific heart.** Concordance is supporting breadth,
   not the headline. The novel biology (essentiality → morphological DI,
   knockout/overexpression correlation) leads the paper.

## What changed from V2 → V3

5. **Partial correlation control on ALL concordance experiments (A, C).**
   Every rho-based test now reports partial rho controlling for n_contexts.
   Without this, correlated sample sizes produce spurious concordance —
   the exact confound Paper 1 characterizes. Primary criteria use partial rho.

6. **Experiment C threshold raised from rho > 0.10 to rho > 0.20.**
   rho = 0.10 explains 1% of variance. Pre-registering it as "confirmed"
   would commit the significant-but-trivial sin the combined paper condemns.
   Effects in [0.10, 0.20] reported as "weak/negligible."

7. **Experiment D: permutation null added.**
   CRISPR and ORF share plates/batches. Shuffled-gene-label null (1000 perms)
   required to show observed rho exceeds technical-structure baseline.
   Both rho > 0.15 AND exceeding 95th percentile of null required.

8. **Experiment D: directional divergence analysis added.**
   Pre-registered exploratory analysis of genes where knockout and overexpression
   DI diverge — knockout-stable/OE-unstable vs reverse. This is where the
   memorable biology lives (dosage-sensitive GOF vs paralog redundancy).

9. **Multiple comparisons: Bonferroni across 6 primary tests.**
   Alpha = 0.0083 per test. Exploratory analyses (A2, Exp D divergence)
   reported with uncorrected p-values and explicit "exploratory" labels.

10. **Bootstrap CIs per perturbation in atlas columns.**
    `di_corrected_ci_low`, `di_corrected_ci_high`, `di_raw_ci_low`,
    `di_raw_ci_high` added. No existing perturbation resource ships
    per-entry uncertainty. Differentiator for Nature Methods vs Scientific Data.

11. **Tahoe pseudobulk minimum-cell threshold: >= 20 cells.**
    Matching Perturb-seq threshold. Supplementary figure confirms flat-bin
    immunity carries over. Pre-registered before data inspection.
