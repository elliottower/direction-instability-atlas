# Direction Instability Atlas: Pre-Registered Results

**Batch 1 SHA:** `0d07a01` | **Batch 2A SHA:** `a17c125` | **Date:** 2026-07-06 | **Status:** All experiments complete

---

## Executive summary

Across two pre-registered batches (9 completed tests), one clear pattern emerges: DI is a real, cross-modal property of perturbation-context interactions that is *not* reducible to simpler gene-level features. Both concordance tests point in the right direction (A2 passes, A1 shows a comparable effect size but lacks power to clear the CI gate), while every attempt to explain DI variation through proxies -- essentiality, expression variance, paralog buffering, drug selectivity -- either fails outright or produces effects too small to matter biologically. DI measures something specific about how perturbations cascade through cellular networks, and that something is not captured by existing gene annotations.

### Batch 1 (alpha = 0.0125)

| Test | n | Partial rho | 95% CI | p | Verdict |
|------|---|-------------|--------|---|---------|
| A1: LINCS-Tahoe | 76 | 0.254 | [0.03, 0.46] | 0.029 | **FAIL (CI gate, underpowered)** |
| A2: LINCS-JUMP | 378 | 0.190 | [0.09, 0.29] | 2.1e-4 | **PASS** |
| C: Essentiality | 7,653 | -0.332 | [-0.35, -0.31] | ~0 | **FAIL (wrong sign)** |
| D: KO vs OE | 5,220 | 0.011 | [-0.02, 0.04] | 0.43 | **FAIL (null)** |
| PRISM (expl.) | 75 | 0.359 | [0.14, 0.54] | 0.0016 | **Exploratory PASS** |

### Batch 2A (alpha = 0.0167)

| Test | n | Partial rho | 95% CI | p | Verdict |
|------|---|-------------|--------|---|---------|
| E2: Expression variance | 7,802 | 0.060 | [0.038, 0.082] | 1.1e-7 | **FAIL (CI gate)** |
| E3: Drug selectivity (gini_corrected) | 75 | 0.144 | [-0.088, 0.360] | 0.22 | **FAIL (underpowered)** |
| E5a: Paralog count | 6,894 | -0.007 | [-0.031, 0.017] | 0.57 | **FAIL (null)** |
| E5b: Paralog breadth | 6,882 | -0.029 | [-0.053, -0.006] | 0.015 | **FAIL (CI gate)** |

---

## Experiment A2: Cross-modality concordance passes

**Result:** Partial Spearman rho = 0.190, 95% CI [0.090, 0.285], p = 2.1 x 10^-4, n = 378 compounds.

Compounds matched between LINCS L1000 (gene expression, 978 landmark genes) and JUMP-CP Cell Painting (morphology, 3,180 features) by InChIKey connectivity prefix (14-character). Pre-registration estimated n ~ 330; actual overlap is 378 after matching via `pert_iname` in the LINCS perturbation metadata. One covariate (LINCS n_contexts) controlled via Pearson-residualize-then-Spearman.

**Interpretation.** This is the construct validity result. Compounds whose transcriptomic effects point in different directions across cell lines also produce variable morphological profiles. Because gene expression and cell morphology are distinct measurement modalities with different noise structures, confounds, and dimensionalities, this concordance indicates that DI reflects a biological property of the compound-context interaction rather than a measurement artifact. The effect size (rho ~ 0.19) falls in the "weakly concordant" bin (0.15-0.30) from the pre-registration, which is expected for cross-modal comparisons where only part of the biology is shared.

**Note on matching.** The pre-registration specified InChIKey matching but did not distinguish between `pert_id` and `pert_iname` fields in the LINCS perturbation metadata. Matching on `pert_id` yielded only 8 compounds (many LINCS DI drugs are indexed by BRD IDs that appear as `pert_iname`, not `pert_id`, in the metadata). Matching on `pert_iname` recovers 378 compounds, close to the pre-registered estimate. This is a data-wrangling correction, not a methodological change.

---

## Experiment C: Essentiality predicts *lower* DI

**Result:** Partial Spearman rho = -0.332, 95% CI [-0.352, -0.312], p ~ 0, n = 7,653 genes. DepMap 24Q4 Chronos (figshare 51064667; 1,178 cell lines, 17,916 genes).

The pre-registration predicted a positive correlation: genes essential across more cell lines should show higher morphological DI when knocked out. The observed correlation is strongly negative. One covariate (n_contexts = number of plate replicates per gene in JUMP-CP CRISPR) controlled via Pearson-residualize-then-Spearman.

**What went wrong with the prediction.** The hypothesis confused two senses of "context-dependent." Essential genes affect many cell types (broad essentiality), but they affect those cell types *in the same way*. A gene required for mitosis, ribosome assembly, or DNA repair produces a stereotyped catastrophe regardless of cellular identity. Broadly essential genes are functionally constrained, and that constraint manifests as low DI: similar phenotypic disruption across contexts.

Genes with high DI are more likely to be context-specific in their effects -- knocking them out matters in some cell types and not others, or matters differently. These tend to be non-essential or selectively essential genes whose functions interact with cell-type-specific pathways.

**Why this is informative.** The effect is real, large, and precisely estimated. The pre-registration committed to reporting the result regardless of direction, and the wrong-sign outcome reveals a genuine structure: functional constraint and phenotypic consistency are the same axis. This connects DI to decades of work on gene essentiality, fitness landscapes, and phenotypic robustness.

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

## Experiment A1: LINCS-Tahoe same-modality concordance -- FAIL (underpowered, effect present)

**Result:** Partial Spearman rho = 0.254, 95% CI [0.027, 0.456], p = 0.029, n = 76 drugs. Two covariates (LINCS n_contexts, Tahoe n_contexts). **FAIL** -- CI lower bound (0.027) below 0.05 threshold, and p (0.029) does not clear Bonferroni alpha (0.0125).

Matching LINCS L1000 and Tahoe-100M drugs by the pre-registered method (InChIKey connectivity prefix, 14 characters, verified via PubChem CID cross-reference) yielded 76 overlapping compounds. Tahoe-100M DI was computed from 1,026 pseudobulk differential expression shards (tahoebio/Tahoe-100M, HuggingFace), yielding 379 drugs with at least 5 cell lines (actual range: 47-50 per drug).

The point estimate (0.254) is larger than the passing cross-modality test A2 (0.190), indicating the failure is one of statistical power at n = 76, not absence of concordance. At this sample size the standard error is approximately 0.12, so clearing the CI-lower >= 0.05 gate would require an observed rho above approximately 0.35. The pre-registered power analysis anticipated this constraint, noting A1 reaches 80% power only for true effects >= 0.35.

**Deviation disclosed (matching method).** An initial run used case-insensitive drug-name matching (n = 114), which returned a null (rho = -0.003). Name matching introduced approximately 38 spurious pairs -- salt forms, stereoisomers, and name collisions -- of the same false-positive class removed by InChIKey verification in Batch 1. Reverting to the pre-registered InChIKey method (n = 76) both complies with the freeze and resolves the null. The name-matched result is reported here only to document the deviation and its correction; the InChIKey result is the pre-registered analysis.

**On the A1/A2 pattern.** The pre-registered interpretation guide did not anticipate the observed pattern (cross-modality A2 passing while same-modality A1 is underpowered). We interpret this as a sample-size artifact -- the LINCS-Tahoe overlap (n = 76) is far smaller than LINCS-JUMP (n = 378) -- rather than evidence that DI is more transferable across modalities than within one. The A1 effect size (0.254 > 0.190) supports this reading.

---

# Batch 2A: Mechanistic Correlates

**Pre-registration:** PREREGISTRATION_BATCH2.md, SHA `a17c125`. **Bonferroni alpha = 0.05/3 = 0.0167.** All three experiments use the Pearson-residualize-then-Spearman partial correlation method from Batch 1, with an additional CI-lower-bound > 0.05 gate to prevent trivially small effects from passing at large n.

---

## Experiment E2: Expression variance -- real but negligible

**Result:** Partial Spearman rho = 0.060, 95% CI [0.038, 0.082], p = 1.1 x 10^-7, n = 7,802 genes. One covariate (n_contexts). **FAIL** -- CI lower bound (0.038) below 0.05 threshold.

The hypothesis: genes whose expression varies across cell lines (shallow attractor basins, context-dependent regulation) should produce higher morphological DI when knocked out. The correlation exists and is statistically unambiguous (p ~ 10^-7), but rho = 0.060 explains less than 0.4% of variance. The CI gate did its job: at n = 7,802, the study detects effects as small as rho = 0.03, and the observed effect, while real, is too small to constitute biological signal.

**Interpretation.** Expression variance captures something genuine but peripheral. Genes with variable expression across cell lines do produce marginally more variable knockout phenotypes, consistent with the Waddington basin intuition (shallow basins = context-dependent regulation = context-dependent phenotypes). But the effect is a footnote, not a headline. DI variation is driven by something expression variance does not index.

**Data source.** DepMap 24Q4 expression (OmicsExpressionProteinCodingGenesTPMLogp1.csv, figshare 51065489). 1,673 cell lines, 19,193 genes. 7,802 genes overlap with CRISPR DI.

---

## Experiment E3: Drug selectivity -- magnitude confound confirmed, signal inconclusive

**Result:** Partial Spearman rho = 0.144, 95% CI [-0.088, 0.360], p = 0.22, n = 75 drugs. One covariate (LINCS n_contexts). **FAIL** -- p > 0.0167 (underpowered).

The hypothesis: high-DI drugs are selective killers (kill some cell lines, spare others), while low-DI drugs kill broadly or not at all. The primary metric is magnitude-corrected Gini (`gini_corrected`), residualizing Gini of kill depth against mean kill magnitude via OLS. This correction is critical: raw Gini mechanically correlates with kill magnitude because drugs that kill more deeply have more dynamic range to be selective over.

**Sensitivity analyses.** Raw (confounded) Gini: rho = 0.275. Kill fraction (cells with LFC < -0.5): rho = -0.453 (opposite sign, reflecting that drugs killing many cell lines tend to have *lower* DI). The difference between gini_raw (0.275) and gini_corrected (0.144) confirms the magnitude confound is real and inflates the apparent selectivity-DI relationship.

**Why inconclusive rather than null.** n = 75 is the LINCS-PRISM overlap after filtering (|mean_lfc| >= 0.01), a data availability constraint. The CI spans [-0.088, 0.360], which is consistent with anything from no effect to a moderate one. The pre-registered power analysis assumed n ~ 1,000+; actual n is an order of magnitude smaller. This experiment cannot distinguish a real 0.15 effect from zero.

---

## Experiment E5a: Paralog count -- entirely confounded

**Result:** Partial Spearman rho = -0.007, 95% CI [-0.031, 0.017], p = 0.57, n = 6,894 genes. Two covariates (n_contexts, mean expression level). **FAIL** -- effect is zero.

The hypothesis: genes with more paralogs show higher knockout DI because paralog compensation buffers loss-of-function in some cell types but not others, producing context-dependent phenotypes.

**The raw correlation was entirely confounded.** Raw Spearman rho = 0.038 (p = 0.002) -- a small but statistically significant positive association. After controlling for mean expression level and n_contexts, the partial rho is -0.007, indistinguishable from zero. The raw correlation existed because highly expressed genes tend to have more paralogs (gene family expansion correlates with expression level) AND have higher DI (for reasons related to measurement, not paralog buffering). Once expression level is controlled, paralog count has no relationship to DI.

This is a textbook confound removal: a statistically significant raw association (p = 0.002) that evaporates completely under covariate control.

**9 genes dropped** from the original 6,903 overlap due to NaN/inf in mean expression values (genes present in CRISPR DI and Ensembl paralogs but absent from the DepMap expression file).

---

## Experiment E5b: Paralog expression breadth -- tiny, gated

**Result:** Partial Spearman rho = -0.029, 95% CI [-0.053, -0.006], p = 0.015, n = 6,882 genes with at least one expressed paralog. Two covariates (n_contexts, mean expression level). **FAIL** -- CI upper bound (-0.006) above -0.05 threshold.

The hypothesis (secondary): among genes with expressed paralogs, broader paralog expression (compensation available in more cell types) should predict *lower* DI (more uniform buffering = less context-dependence). Expected sign: negative.

The sign is correct (rho = -0.029, negative as predicted) and p = 0.015 would pass the alpha threshold (0.0167). But the CI gate blocks: the upper bound of the CI (-0.006) does not clear -0.05. The effect, while in the right direction and marginally significant, is too small to claim biological relevance. Same pattern as E2: the CI gate prevents over-interpretation of tiny effects at large n.

**Raw vs. partial.** Raw rho = -0.116 (p ~ 10^-22) is substantially larger, suggesting that expression-level confounds inflate the apparent paralog-breadth effect as well.

---

## Synthesis: What the pattern means

Nine completed results across two batches form a coherent picture that the original hypotheses only partially anticipated.

**DI is real and concordant across datasets (A1, A2, PRISM).** The construct validity claim holds: A2 passes (rho = 0.19, n = 378) and A1 points the same way (rho = 0.25, n = 76) though underpowered to clear the CI gate. PRISM viability CV also correlates with DI (rho = 0.36, exploratory). Together, these rule out the "DI is just noise" null hypothesis.

**DI is mechanism-specific, not target-intrinsic (D).** Knocking out and overexpressing the same gene produce unrelated patterns of context-dependence. DI characterizes the perturbation-context interaction, not the gene. This has practical implications: a DI atlas must be stratified by perturbation type, not collapsed to gene-level.

**Functional constraint suppresses DI (C).** Essential genes produce consistent phenotypes, yielding low DI. High DI marks genes whose effects depend on cellular context. This inverts the original prediction but connects DI to established concepts in functional genomics.

**DI is not reducible to simpler gene features (E2, E5a, E5b).** Batch 2A tested three plausible proxies for context-dependence: expression variance (Waddington basin depth), paralog count (genetic redundancy), and paralog expression breadth (selective compensation). All three fail after proper confound control. Expression variance produces a real but negligible effect (rho = 0.064, < 0.5% variance explained). Paralog count's raw association (rho = 0.038) is entirely a confound of expression level. Paralog breadth shows a correctly-signed but trivially small effect (rho = -0.029). These are not power failures -- n ranges from 6,882 to 7,802, and effects as small as rho = 0.03 are detectable.

**The confound corrections validated themselves (E3, E5a).** E3's magnitude correction reduced the apparent Gini-DI association from 0.275 to 0.144, confirming the magnitude confound was real. E5a's covariate control eliminated a statistically significant raw correlation entirely (rho = 0.038, p = 0.002 down to rho = -0.007, p = 0.57). The CI gate blocked two additional trivially small effects (E2, E5b) that would have passed on p-value alone. The methodological investment in Batch 2A's pre-registration paid off: without these corrections, two of the three experiments would have been false positives.

**Together.** DI measures the diversity of a perturbation's downstream effects across cellular contexts, where "downstream effects" are mechanism-specific. A compound's DI reflects how its target engagement cascades differently through different cellular networks. An essential gene's knockout DI is low because the disruption is upstream of context-specific buffering. A non-essential gene's knockout DI is high because the phenotypic consequences depend on which pathways the cell has available to compensate. These same genes, when overexpressed, engage entirely different compensatory mechanisms, producing an unrelated DI profile. Crucially, this context-dependence is not predicted by the gene's expression variance, its paralog redundancy, or the selectivity of its pharmacological targeting -- it emerges from the interaction between the perturbation mechanism and the cell's compensatory network, which existing gene annotations do not capture.

---

## Deviations from pre-registration

### DepMap version

All DepMap analyses use the pre-registered **24Q4** release (figshare article 27993248): Chronos (file 51064667; 1,178 cell lines, 17,916 genes) for Exp C, Expression (file 51065489; 1,673 cell lines, 19,193 genes) for E2 and E5. Single panel, no cross-release mismatch.

*Note:* Initial runs used 22Q2 Chronos (Exp C) and 24Q2 Expression (E2/E5). Re-run on 24Q4 changed no verdicts (Exp C: rho -0.330 -> -0.332; E2: 0.064 -> 0.060; E5a/b: identical to 3 decimal places).

### Other deviations

1. **A1 n and matching (Batch 1):** Pre-registration estimated n ~ 103. InChIKey matching via PubChem CID cross-reference yielded n = 76 (not all Tahoe drugs have PubChem CIDs mapping to LINCS InChIKeys). An initial run erroneously used case-insensitive drug-name matching (n = 114, rho = -0.003); reverting to the pre-registered InChIKey method corrected this. See A1 section for details.
2. **InChIKey matching field (Batch 1):** `pert_iname` used instead of `pert_id` for LINCS-JUMP matching. `pert_id` matching was a data-wrangling error yielding only 8 compounds.
3. **A2 n (Batch 1):** 378 compounds instead of pre-registered 330. Difference due to `pert_iname` matching recovering more valid matches.
4. **Context definitions (Batch 1):** LINCS/Tahoe contexts are cell lines; JUMP-CP CRISPR/ORF contexts are plate replicates.
5. **E3 n (Batch 2A):** 75 drugs instead of pre-registered ~1,000+. The LINCS-PRISM overlap after filtering (|mean_lfc| >= 0.01) is much smaller than anticipated. Data availability issue, not a methodological change.
6. **E5 dropped rows (Batch 2A):** 9 genes (of 6,903) dropped due to NaN/inf in expr_mean covariate. No impact on conclusions at n = 6,894.

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
| Tahoe DI | 379 drugs | results/tahoe_di.csv |
| LINCS-Tahoe A1 | 76 drugs (InChIKey) | results/concordance.json |
| Expression variance | 7,802 genes | results/expression_variance.json |
| PRISM selectivity | 75 drugs | results/prism_selectivity.json |
| Paralog buffering | 6,894 genes (E5a), 6,882 (E5b) | results/paralog_buffering.json |
| DepMap 24Q4 Chronos | 1,178 x 17,916 | data/depmap/24q4/CRISPRGeneEffect.csv |
| DepMap 24Q4 Expression | 1,673 x 19,193 | data/depmap/24q4/OmicsExpressionProteinCodingGenesTPMLogp1.csv |
| Ensembl paralogs | 2,784,678 pairs | data/depmap/ensembl_paralogs.csv |

---

## Outstanding items

- [x] ~~Experiment A1: LINCS-Tahoe concordance~~ (FAIL, underpowered, effect present)
- [x] ~~Full script 05 run (A1 + A2 together)~~ (results/concordance.json)
- [x] ~~Standardize DepMap on 24Q4~~ (all experiments re-run, no verdicts changed)
- [ ] Context-definition consistency note in methods
- [ ] Gene Ontology enrichment for directional divergence gene sets (exploratory)
- [ ] Batch 2B pre-registration: E1 (Shesha vs DI, requires Perturb-seq pipeline)
- [ ] Larger LINCS-Tahoe overlap (more Tahoe drugs with InChIKey annotation) for A1 follow-up
