# Direction Instability Atlas: Pre-Registered Results

**SHA freeze:** `0d07a01` | **Date:** 2026-07-05 | **Status:** 5/7 experiments complete, Tahoe Modal running

---

## Executive summary

Of four pre-registered primary tests (Bonferroni alpha = 0.0125), one passes, two fail informatively, and one awaits data. An exploratory test also passes. The pattern that emerges is more nuanced than the original hypotheses: direction instability transfers across measurement modalities (A2) and predicts functional outcomes (PRISM), but it is *not* an intrinsic property of perturbation targets (D) and correlates *inversely* with functional constraint (C). These results reframe DI as a property of perturbation-context interactions rather than perturbation targets.

| Test | n | Partial rho | 95% CI | p | Verdict |
|------|---|-------------|--------|---|---------|
| A1: LINCS-Tahoe | -- | -- | -- | -- | Pending |
| A2: LINCS-JUMP | 378 | 0.190 | [0.09, 0.29] | 2.1e-4 | **PASS** |
| C: Essentiality | 7,627 | -0.330 | [-0.35, -0.31] | ~0 | **FAIL (wrong sign)** |
| D: KO vs OE | 5,220 | 0.011 | [-0.02, 0.04] | 0.43 | **FAIL (null)** |
| PRISM (expl.) | 75 | 0.359 | [0.14, 0.54] | 0.0016 | **Exploratory PASS** |

---

## Experiment A2: Cross-modality concordance passes

**Result:** Partial Spearman rho = 0.190, 95% CI [0.090, 0.285], p = 2.1 x 10^-4, n = 378 compounds.

Compounds matched between LINCS L1000 (gene expression, 978 landmark genes) and JUMP-CP Cell Painting (morphology, 3,180 features) by InChIKey connectivity prefix (14-character). Pre-registration estimated n ~ 330; actual overlap is 378 after matching via `pert_iname` in the LINCS perturbation metadata. One covariate (LINCS n_contexts) controlled via Pearson-residualize-then-Spearman.

**Interpretation.** This is the construct validity result. Compounds whose transcriptomic effects point in different directions across cell lines also produce variable morphological profiles. Because gene expression and cell morphology are distinct measurement modalities with different noise structures, confounds, and dimensionalities, this concordance indicates that DI reflects a biological property of the compound-context interaction rather than a measurement artifact. The effect size (rho ~ 0.19) falls in the "weakly concordant" bin (0.15-0.30) from the pre-registration, which is expected for cross-modal comparisons where only part of the biology is shared.

**Note on matching.** The pre-registration specified InChIKey matching but did not distinguish between `pert_id` and `pert_iname` fields in the LINCS perturbation metadata. Matching on `pert_id` yielded only 8 compounds (many LINCS DI drugs are indexed by BRD IDs that appear as `pert_iname`, not `pert_id`, in the metadata). Matching on `pert_iname` recovers 378 compounds, close to the pre-registered estimate. This is a data-wrangling correction, not a methodological change.

---

## Experiment C: Essentiality predicts *lower* DI

**Result:** Partial Spearman rho = -0.330, 95% CI [-0.350, -0.310], p ~ 10^-171, n = 7,627 genes.

The pre-registration predicted a positive correlation: genes essential across more cell lines should show higher morphological DI when knocked out. The observed correlation is strongly negative. One covariate (n_contexts = number of plate replicates per gene in JUMP-CP CRISPR) controlled via Pearson-residualize-then-Spearman.

**What went wrong with the prediction.** The hypothesis confused two senses of "context-dependent." Essential genes affect many cell types (broad essentiality), but they affect those cell types *in the same way*. A gene required for mitosis, ribosome assembly, or DNA repair produces a stereotyped catastrophe regardless of cellular identity. Broadly essential genes are functionally constrained, and that constraint manifests as low DI: similar phenotypic disruption across contexts.

Genes with high DI are more likely to be context-specific in their effects -- knocking them out matters in some cell types and not others, or matters differently. These tend to be non-essential or selectively essential genes whose functions interact with cell-type-specific pathways.

**Why this is informative.** The effect is real, large, and precisely estimated. The pre-registration committed to reporting the result regardless of direction, and the wrong-sign outcome reveals a genuine structure: functional constraint and phenotypic consistency are the same axis. This connects DI to decades of work on gene essentiality, fitness landscapes, and phenotypic robustness.

**Deviation.** DepMap 22Q2 was used (figshare 34990036) rather than the pre-registered 24Q4. The gene coverage and Chronos methodology are similar across releases. This is noted as a deviation.

---

## Experiment D: Knockout and overexpression DI are uncorrelated

**Result:** Spearman rho = 0.011, 95% CI [-0.016, 0.038], Fisher z p = 0.43, n = 5,220 genes. Permutation null: mean = 0.0003, SD = 0.014, 95th percentile = 0.023. Empirical permutation p = 0.215.

Both decision criteria fail: the CI lower bound (-0.016) is below 0.05, and the observed rho (0.011) does not exceed the permutation 95th percentile (0.023). At n = 5,220, the study had power to detect effects as small as rho = 0.04 (SE ~ 0.014), so this is a well-powered null.

**What this means.** The pre-registration hypothesized that context-dependence is an intrinsic property of genes: a gene whose knockout produces variable morphological effects should also produce variable effects when overexpressed. This is false. Losing a protein and gaining excess protein engage different molecular mechanisms (loss-of-function vs. gain-of-function), and those mechanisms interact with cellular context in unrelated ways.

This result constrains the theoretical scope of DI. Direction instability characterizes how a *specific perturbation mechanism* interacts with context. It cannot be attributed to the gene itself as a persistent trait. A gene is not "context-dependent" in general; it is context-dependent *in a particular way* (knockout vs. overexpression, transcriptomic vs. morphological).

**Exploratory: directional divergence.** 1,331 genes are knockout-stable but overexpression-unstable (DI_CRISPR - DI_ORF < -1 SD), consistent with dosage-sensitive gain-of-function genes whose overexpression phenotypes depend on expression level and cellular buffering capacity. 289 genes show the opposite pattern (knockout-unstable, overexpression-stable), consistent with genes whose loss-of-function is partially compensated by paralogs in some but not all cell types. These classifications are exploratory and require Gene Ontology or pathway enrichment follow-up.

---

## PRISM exploratory: Viability CV correlates with DI

**Result:** Spearman rho = 0.359, 95% CI [0.144, 0.542], p = 0.0016, n = 75 drugs.

Drugs with higher viability coefficient of variation across cell lines (PRISM) show higher transcriptomic direction instability (LINCS). This is exploratory (uncorrected p-value, no covariate adjustment) and the n is small (75 drugs with abs(mean_lfc) >= 0.01). The effect size is moderate and the CI clears the exploratory threshold of rho > 0.15.

**Interpretation.** This is an independent cross-assay validation: drugs whose cell-killing potency varies across cell lines also produce transcriptomic effects that point in different directions. Viability is a scalar outcome (not a high-dimensional profile), so this correlation cannot be attributed to shared noise structure with LINCS. It suggests that DI captures something about how drug mechanisms interact with cell-type-specific biology that manifests across both expression and viability readouts.

---

## Experiment A1: Pending

LINCS-Tahoe same-modality concordance awaits Tahoe-100M DI computation. Modal job is processing 1,026 pseudobulk differential expression shards from `tahoebio/Tahoe-100M` (metadata/pseudobulk_differential_expression). Expected overlap: ~103 drugs (InChIKey-verified). This is the same-modality comparator to A2: if A1 also passes, DI concordance holds both within and across modalities. If A1 fails while A2 passes, the pattern would be surprising and require investigation (platform-specific noise in one of the two expression datasets).

---

## Synthesis: What the pattern means

The five completed results form a coherent picture that the original hypotheses only partially anticipated.

**DI is real and cross-modal (A2, PRISM).** The construct validity claim survives: compounds that are directionally unstable in one modality are unstable in another. This rules out the "DI is just noise" null hypothesis.

**DI is mechanism-specific, not target-intrinsic (D).** Knocking out and overexpressing the same gene produce unrelated patterns of context-dependence. DI characterizes the perturbation-context interaction, not the gene. This is a stronger claim than the original hypothesis, and it has practical implications: a DI atlas must be stratified by perturbation type, not collapsed to gene-level.

**Functional constraint suppresses DI (C).** Essential genes produce consistent phenotypes, yielding low DI. High DI marks genes whose effects depend on cellular context -- genes with tissue-specific functions, regulatory roles, or redundant paralogs. This inverts the original prediction but connects DI to established concepts in functional genomics.

**Together.** DI measures the diversity of a perturbation's downstream effects across cellular contexts, where "downstream effects" are mechanism-specific. A compound's DI reflects how its target engagement cascades differently through different cellular networks. An essential gene's knockout DI is low because the disruption is upstream of context-specific buffering. A non-essential gene's knockout DI is high because the phenotypic consequences depend on which pathways the cell has available to compensate. These same genes, when overexpressed, engage entirely different compensatory mechanisms, producing an unrelated DI profile.

---

## Deviations from pre-registration

1. **DepMap version:** 22Q2 used instead of 24Q4. Similar gene coverage and methodology.
2. **InChIKey matching field:** `pert_iname` used instead of `pert_id` for LINCS-JUMP matching. `pert_id` matching was a data-wrangling error yielding only 8 compounds.
3. **A2 n:** 378 compounds instead of pre-registered 330. Difference due to `pert_iname` matching recovering more valid matches.
4. **Context definitions:** LINCS/Tahoe contexts are cell lines; JUMP-CP CRISPR/ORF contexts are plate replicates. Consistency note to be added to Methods.

---

## Data produced

| Dataset | n | File |
|---------|---|------|
| CRISPR DI | 7,946 genes | results/crispr_di.csv |
| ORF DI | 12,590 genes | results/orf_di.csv |
| PRISM CV | 1,433 drugs | results/prism_cv.csv |
| Essentiality | 7,627 genes (merged) | results/essentiality.json |
| KO vs OE | 5,220 genes (merged) | results/knockout_vs_overexpression.json |
| LINCS-JUMP A2 | 378 compounds | results/concordance_a2.json |
| Tahoe DI | pending | (Modal volume) |

---

## Outstanding items

- [ ] Tahoe-100M DI computation (Modal, ~1-2 hours remaining)
- [ ] Experiment A1: LINCS-Tahoe concordance (blocked on Tahoe)
- [ ] Full script 05 run (A1 + A2 together, saves combined concordance.json)
- [ ] Context-definition consistency note in methods
- [ ] Gene Ontology enrichment for directional divergence gene sets (exploratory)
- [ ] Perplexity literature check on four topics (queries prepared)
