# Data Inventory

Generated 2026-07-05 during Phase 1-2 (pre-registration prep).
No DI values have been computed on new datasets. All numbers are metadata-only.

## Dataset sizes and formats

### 1. LINCS L1000 (already computed)
- **Source:** drug-perturbation-geometry/zenodo_v1/drug_instability_8949.csv
- **Format:** CSV with DI, magnitude_cv, frechet_variance, etc. per drug
- **Size:** 8,949 drugs, each with DI from >= 5 cell lines
- **Identifiers:** drug_name (BRD IDs + names), MOA, target
- **Status:** DI already computed in bracket-norm-validity

### 2. JUMP-CP compounds (already computed)
- **Source:** bracket-norm-validity/results/08_jump_cp/
- **Format:** npz (DI values), CSV (metadata)
- **Size:** 25,254 compounds, DI across source contexts
- **Identifiers:** JCP2022 IDs, InChIKey via compound_metadata.csv.gz
- **Status:** DI already computed

### 3. Tahoe-100M
- **Source:** HuggingFace tahoebio/Tahoe-100M
- **Format:** 1,026 parquet shards of pseudobulk differential expression
- **Size on disk:** ~89 GB (pseudobulk DE tables) — TOO LARGE FOR LOCAL
- **Drugs:** 379, median 15 cell lines per drug, ALL have >= 5
- **Cell lines:** 50 (102 in metadata, ~50 used in pseudobulk)
- **Genes per comparison:** ~62,710 (full transcriptome)
- **Identifiers:** drug name, PubChem CID (377/379), SMILES, MOA (broad + fine), targets
- **Key columns:** gene_name, log2FoldChange, lfcSE, n_cells_trt, Cell_Name_Vevo, drug
- **Cell lines have DepMap IDs** (Cell_ID_DepMap in cell_line_metadata)
- **Minimum cells per pseudobulk:** median 1,728 (well above 20-cell threshold)
- **Processing plan:** Stream on Modal (too large for laptop with 32GB free).
  Extract log2FoldChange vectors per drug×cell_line, compute DI.

### 4. PRISM Repurposing 24Q2
- **Source:** Figshare (article 25917643)
- **Format:** CSV
- **Size on disk:** ~1.5 GB total
- **Drugs:** 6,575 unique (across 3 screens: REP.300, REP.PRIMARY, REP.1M)
- **Cell lines:** 906 (859 passing QC)
- **Readout:** Single-dose viability (2.5 µM, 5-day), NOT dose-response curves
- **Identifiers:** Drug.Name, BRD IDs, MOA (6,115/6,790 with MOA), target (4,414/6,790)
- **No InChIKey in compound list** — will need to match by name or BRD ID
- **DI is NOT computable** (viability is 1D per cell line). Use viability CV instead.
- **LFC_COLLAPSED.csv (143 MB):** DOWNLOADED. Long format, 1.49M rows.
  Columns: row_id, broad_id, dose, compound_plate, screen, culture, LFC.
  Each row is one drug × cell_line × plate observation.

### 5. JUMP-CP CRISPR knockouts
- **Source:** Cell Painting Gallery, cpg0016-jump-assembled
- **Profile URL:** profiles_assembled/CRISPR/v1.0a/...sphering_harmony_PCA_corrected.parquet
- **Format:** Parquet, well-level profiles (post-sphering, harmony, PCA)
- **Size on disk:** 76.2 MB
- **Shape:** 51,185 wells × 259 features (PCA-reduced)
- **Perturbations:** 7,977 (7,948 with >= 5 reps)
- **Plates:** 148 (all source_13)
- **Gene overlap with ORF:** 5,220 genes with >= 5 reps in both
- **Identifiers:** JCP2022 IDs → gene symbols via metadata
- **Status:** DOWNLOADED

### 6. JUMP-CP ORF overexpression
- **Source:** Cell Painting Gallery, cpg0016-jump-assembled
- **Profile URL:** profiles_assembled/ORF/v1.0a/...sphering_harmony.parquet
- **Format:** Parquet, well-level profiles (post-sphering, harmony, no PCA)
- **Size on disk:** 345.1 MB
- **Shape:** 81,660 wells × 722 features
- **Perturbations:** 15,131 (15,121 with >= 5 reps, 12,597 treatment genes)
- **Plates:** 225 (all source_4)
- **Identifiers:** JCP2022 IDs → gene symbols via metadata
- **Note:** Different feature space from CRISPR (722 vs 259 features,
  no PCA correction). DI is computed within each modality so this is fine.
- **Status:** DOWNLOADED

## Cross-reference overlaps (by drug name matching)

| Pair | Overlap | Notes |
|------|---------|-------|
| LINCS ∩ Tahoe | **114** | Primary same-modality concordance test |
| LINCS ∩ PRISM | 1,726 | Large overlap; PRISM used for viability CV only |
| LINCS ∩ JUMP-CP | 330 | From combined paper MOA matching |
| Tahoe ∩ PRISM | 198 | |
| All three (LINCS ∩ Tahoe ∩ PRISM) | 106 | |

### Power implications
- **Exp A same-modality (n=114):** Partial rho > 0.40 with 2 covariates:
  at n=114, rho=0.40 gives p << 0.001. Well powered.
- **Exp A cross-modality (n=330):** Partial rho > 0.20: well powered.
- **Exp A2 (n=1726):** PRISM CV vs LINCS DI: massively powered.

## Disk space budget (32 GB free)

| Item | Size | Local? |
|------|------|--------|
| PRISM (full download) | ~1.5 GB | Yes |
| Tahoe metadata | ~5 MB | Yes (done) |
| Tahoe pseudobulk | ~89 GB | NO — use Modal |
| JUMP-CP CRISPR profiles | TBD (~5-10 GB?) | Maybe |
| JUMP-CP ORF profiles | TBD (~5-10 GB?) | Maybe |

**Tahoe must be processed on Modal.** Stream parquet shards, extract
log2FoldChange vectors, compute DI, download only the result (few KB).
