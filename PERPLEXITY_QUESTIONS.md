# Questions for Perplexity: DI Atlas Data Sourcing

Copy-paste these into Perplexity. Each section is a self-contained question.

---

## 1. PRISM data format and download

I'm computing "direction instability" (mean pairwise cosine distance across cell lines) for drugs in the PRISM Repurposing dataset from the Broad Institute DepMap project. I need to understand the exact data format before I can decide how to compute DI.

Questions:
- What is the exact data format of the PRISM Repurposing 24Q2 secondary screen? Is it a matrix of AUC values (one per drug × cell line), or raw dose-response curves (multiple doses per drug × cell line)?
- How many doses are in the secondary screen dose-response? Is it 8 doses as described in the 2020 paper?
- What identifier system does PRISM use for compounds? Is it InChIKey, broad_id, something else? I need to cross-reference with LINCS (which uses pert_id and InChIKey) and JUMP-CP (which uses InChIKey).
- Where exactly do I download the secondary screen data? The DepMap portal, Figshare, or somewhere else? Give me the exact URL and file names.
- Is there a metadata file that maps PRISM compounds to drug names, MOA classes, and target genes? I need MOA annotations to stratify by mechanism type (machinery vs receptor).
- Are there any known quality issues, batch effects, or gotchas with the 24Q2 release?

---

## 2. Tahoe-100M pseudobulk access

I want to compute direction instability on the Tahoe-100M dataset (1,100 drugs × 50 cell lines, single-cell transcriptomics). I do NOT need raw single-cell data — pseudobulk means per drug per cell line are sufficient.

Questions:
- Does the HuggingFace release (tahoebio/Tahoe-100M) include pre-computed pseudobulk expression profiles? The dataset card mentions "pseudobulk_differential_expression" — is this log-fold-change relative to DMSO control, or absolute expression?
- What is the gene space? All ~20k genes, or a subset? How does it compare to LINCS L1000's 978 landmark genes?
- How are drugs identified in Tahoe-100M? Drug names, InChIKeys, SMILES, or some internal ID? I need to cross-reference with LINCS and PRISM.
- How many drugs overlap between Tahoe-100M and LINCS L1000? Has anyone published this cross-reference?
- What cell lines are in Tahoe-100M? Are they DepMap cell lines with Chronos essentiality data available?
- What is the actual file size of the pseudobulk tables? Can I load them on a laptop with 16GB RAM or do I need cloud compute?
- Is there a Python API or do I download parquet/h5ad files directly?

---

## 3. JUMP-CP CRISPR and ORF profiles

I already computed direction instability for 25,254 JUMP-CP compound perturbations (morphological Cell Painting profiles). Now I want to extend to the genetic perturbations in the same dataset.

Questions:
- Where do I download the pre-computed morphological profiles (CellProfiler or deep-learning features) for JUMP-CP CRISPR knockouts (cpg0016)? Is it the Cell Painting Gallery on AWS, or is there a more convenient download?
- Same question for ORF overexpression profiles.
- How are the profiles formatted? Are they already normalized/aggregated to well-level or do I need to aggregate from single-cell?
- How many CRISPR knockouts have profiles in >= 5 source contexts (plates/batches)? I need >= 5 contexts to compute meaningful DI.
- How many ORF genes have profiles in >= 5 source contexts?
- What gene identifiers are used? Entrez, HUGO symbols, or something else? I need to cross-reference with DepMap (Chronos essentiality) and Perturb-seq (Replogle 2022).
- Has anyone published an analysis of cross-context consistency or direction instability on JUMP-CP genetic perturbations? I want to make sure I'm not duplicating existing work.
- The JUMP-CP dataset is 115 TB total. What's the size of just the pre-computed feature profiles for CRISPR + ORF?

---

## 4. Cross-referencing compounds across datasets

I need to match the same drugs/compounds across LINCS L1000, JUMP-CP, PRISM, and Tahoe-100M to compute cross-modality DI concordance.

Questions:
- Is there an existing cross-reference table mapping compounds across these four datasets? Has the Broad or anyone else published one?
- What is the best identifier to use for matching? InChIKey seems standard but do all four datasets provide it?
- The LINCS Repurposing Hub (clue.io) has MOA annotations. Does PRISM use the same annotation system? Does Tahoe-100M provide MOA labels?
- How many compounds are shared across all four datasets? Across any three? I need rough numbers to know if the cross-modality analysis is powered.
- Are there known issues with compound identity mismatches across these datasets (same InChIKey mapping to different compounds, stereoisomer confusion, salt form differences)?

---

## 5. Additional perturbation datasets I might be missing

I'm building a multi-modal direction instability atlas covering transcriptomics (LINCS, Tahoe-100M), morphology (JUMP-CP), and viability (PRISM). Are there other large public perturbation datasets I should include?

Questions:
- Is the L1000 Connectivity Map (CMap) data on GEO (GSE92742, GSE70138) the same data that's in the LINCS drug paper, or is there additional data I'm not using? What's the total compound count in the full CMap release vs what I have (8,949 drugs)?
- Are there other large-scale perturbation screens with drug/compound annotations that are publicly available? I'm thinking of:
  - GDSC (Genomics of Drug Sensitivity in Cancer) — is this redundant with PRISM?
  - CTRPv2 (Cancer Therapeutics Response Portal) — same question
  - PharmacoDB — does this aggregate the above?
  - NCI-60 — old but large
  - Any recent 2025-2026 releases I'm not aware of?
- For genetic perturbations specifically: besides Replogle 2022 Perturb-seq and JUMP-CP CRISPR, are there other large-scale genetic perturbation datasets with multi-context profiles?
- Has anyone published a "perturbation atlas" that tries to unify multiple datasets? I want to cite prior work and differentiate.

---

## 6. Confound controls data

The combined confound paper (separate from this atlas) identified that DI might just track target expression breadth or polypharmacology. I need data to control for these in the atlas paper too.

Questions:
- GTEx bulk RNA-seq: what's the easiest way to get tissue-breadth per gene (number of tissues where the gene is expressed above some threshold)? Is there a pre-computed table I can download rather than processing raw GTEx data?
- Human Protein Atlas: same question — is there a downloadable table of tissue-breadth or "detected in N tissues" per gene?
- ChEMBL: what API endpoint gives me the number of targets per compound (binding affinity < 1µM)? I need to cross-reference ~5,000 drugs. Is there a bulk download that's faster than per-compound API calls?
- STRING/BioGRID: for the "network centrality predicts low DI" hypothesis — what's the easiest way to get degree centrality or betweenness centrality per gene from the protein interaction network? Pre-computed table?

---

## 7. Compute requirements

Questions:
- Tahoe-100M pseudobulk tables: what's the memory footprint? Can I compute DI on a laptop (16GB RAM, Intel Mac, no GPU) or do I need cloud compute?
- JUMP-CP CRISPR/ORF profiles: same question. If profiles are pre-aggregated to well-level, what's the file size?
- PRISM secondary screen: same question.
- For any dataset that's too large for a laptop, what's the cheapest cloud option? I've used Modal and RunPod before.
