# Zenodo Deposit: Direction Instability Atlas

## Upload metadata

**Title:** A Pre-Registered Direction Instability Atlas Across Three Perturbation Modalities

**Authors:** Tower, Elliot (ORCID: _fill in_)

**Upload type:** Dataset

**Publication date:** 2026-07-06

**DOI:** (auto-assigned by Zenodo)

**License:** CC-BY-4.0

**Keywords:**
direction instability, perturbation biology, LINCS, Tahoe-100M, JUMP Cell Painting, pre-registration, context-dependence, drug mechanism, gene essentiality, morphological profiling

**Description (paste into Zenodo abstract field):**

Cross-modal direction instability (DI) atlas for perturbation biology, spanning three measurement technologies and four perturbation classes. DI quantifies the directional consistency of a perturbation's effects across cellular contexts (DI = 1 - mean pairwise cosine similarity). This deposit provides magnitude-corrected DI values with bootstrap confidence intervals for 8,949 LINCS L1000 drugs, 379 Tahoe-100M drugs, 7,948 JUMP-CP CRISPR knockouts, and 12,590 JUMP-CP ORF overexpressions.

Nine pre-registered hypotheses (frozen at commit SHAs 0d07a01 and a17c125 before analysis) characterize what DI captures: cross-modal concordance between expression and morphology confirms that DI reflects biology (partial rho = 0.19, n = 378, p = 2.1e-4); essential genes show lower DI (rho = -0.33), revealing that functional constraint suppresses context-dependence; knockout and overexpression DI are uncorrelated (rho = 0.01, n = 5,220), establishing DI as mechanism-specific; and expression variance, paralog count, and paralog breadth each explain under 0.4% of DI variance after confound control.

This deposit contains the atlas data files (CSV), all experiment result files (JSON), pre-registration documents (frozen before analysis), the full analysis codebase (Python), and the manuscript.

**Related identifiers:**
- Companion paper: https://zenodo.org/records/21213700 (DOI: 10.5281/zenodo.21213700) (isPartOf)
- GitHub repository: https://github.com/_fill-in_/direction-instability-atlas (isSupplementedBy)

---

## File manifest

Upload these files/directories as a single .zip or as individual files:

### Atlas data (the resource)
| File | Description | Format |
|------|-------------|--------|
| `results/crispr_di.csv` | JUMP-CP CRISPR knockout DI, 7,948 genes | CSV |
| `results/orf_di.csv` | JUMP-CP ORF overexpression DI, 12,590 genes | CSV |
| `results/tahoe_di.csv` | Tahoe-100M drug DI, 379 drugs | CSV |
| `results/prism_cv.csv` | PRISM viability CV, 1,433 drugs | CSV |

Each CSV contains: perturbation ID, di_raw, di_corrected, bootstrap CI (low/high), magnitude_cv, mean_norm, n_contexts, dataset.

### Experiment results
| File | Description |
|------|-------------|
| `results/concordance.json` | A1 (LINCS-Tahoe) + A2 (LINCS-JUMP) concordance results |
| `results/essentiality.json` | Experiment C results |
| `results/knockout_vs_overexpression.json` | Experiment D results + permutation null |
| `results/expression_variance.json` | Experiment E2 results |
| `results/prism_selectivity.json` | Experiment E3 results |
| `results/paralog_buffering.json` | Experiment E5a/E5b results |
| `results/prism_concordance.json` | PRISM exploratory results |

### Pre-registration documents (frozen before analysis)
| File | SHA | Description |
|------|-----|-------------|
| `PREREGISTRATION.md` | `0d07a01` | Batch 1: concordance (A1, A2), essentiality (C), KO vs OE (D), PRISM |
| `PREREGISTRATION_BATCH2.md` | `a17c125` | Batch 2A: expression variance (E2), selectivity (E3), paralogs (E5) |

### Analysis code
| Directory/File | Description |
|----------------|-------------|
| `experiments/` | All 13 analysis scripts (01-11 + modal wrapper + utils) |
| `tests/` | Test suite (partial Spearman validation, audit checks, batch 2A) |
| `pyproject.toml` + `uv.lock` | Reproducible Python environment |

### Manuscript
| File | Description |
|------|-------------|
| `paper/atlas_paper_v2.tex` | Full manuscript (LaTeX) |

---

## Data-availability statement

Paste this into the manuscript (Section 2, end of Data availability subsection) and into the journal cover letter:

> All atlas data files, experiment result files, pre-registration
> documents, and analysis code are deposited at Zenodo
> (DOI: 10.5281/zenodo._atlas-DOI_) under a CC-BY-4.0 license.
> Pre-registration commits can be verified against the repository
> history at the SHA values reported in the text.
> LINCS L1000 DI values for 8,949 compounds were previously computed
> and are available from the companion study
> (DOI: 10.5281/zenodo.21213700).
> JUMP-CP compound DI values (25,254 compounds) are available from
> the same source.
> Source datasets are available from their original repositories:
> LINCS L1000 (GEO: GSE92742), Tahoe-100M (HuggingFace:
> tahoebio/Tahoe-100M), JUMP Cell Painting (Cell Painting Gallery,
> cpg0016), DepMap 24Q4 (figshare: 27993248), PRISM (figshare:
> Repurposing\_Public\_24Q2), Ensembl paralogs (BioMart, GRCh38).
