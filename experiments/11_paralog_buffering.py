"""Experiment E5: Paralog buffering predicts CRISPR DI.

Tests whether genes with expressed paralogs show higher knockout DI
(paralog compensation buffers LOF in some contexts but not others).

E5a (primary): partial Spearman(paralog_count, CRISPR DI) controlling
  for n_contexts and mean expression level. 2 covariates.
E5b (secondary): among genes with >= 1 expressed paralog, partial
  Spearman(paralog_expression_breadth, CRISPR DI). Expected: NEGATIVE
  (broader compensation = lower DI).

Method: Pearson-residualize then Spearman. NOT rank-then-residualize.
Decision criterion: CI lower bound > 0.05 (E5a) or CI upper < -0.05
(E5b), p < 0.0167.

Inputs:
  - CRISPR DI: results/crispr_di.csv (from script 02)
  - DepMap expression: data/depmap/OmicsExpressionProteinCodingGenesTPMLogp1.csv
  - Ensembl paralogs: BioMart REST API

Output: results/paralog_buffering.json
"""
import urllib.parse
import urllib.request
from datetime import datetime
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from experiments.utils import (
    ci_lower_bound,
    ci_upper_bound,
    fisher_z_pvalue,
    partial_spearman,
    save_results,
)

CRISPR_DI_PATH = Path("results/crispr_di.csv")
EXPRESSION_CACHE = Path("data/depmap/24q4/OmicsExpressionProteinCodingGenesTPMLogp1.csv")
PARALOG_CACHE = Path("data/depmap/ensembl_paralogs.csv")
RESULTS_PATH = Path("results/paralog_buffering.json")

ALPHA = 0.0167
CI_FLOOR = 0.05
EXPRESSED_TPM_THRESHOLD = 1.0


def fetch_ensembl_paralogs() -> pd.DataFrame:
    if PARALOG_CACHE.exists():
        return pd.read_csv(PARALOG_CACHE)

    PARALOG_CACHE.parent.mkdir(parents=True, exist_ok=True)
    print("Fetching paralogs from Ensembl BioMart...")

    query = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<!DOCTYPE Query>'
        '<Query virtualSchemaName="default" formatter="TSV" header="1"'
        ' uniqueRows="1" count="" datasetConfigVersion="0.6">'
        '<Dataset name="hsapiens_gene_ensembl" interface="default">'
        '<Attribute name="external_gene_name"/>'
        '<Attribute name="hsapiens_paralog_associated_gene_name"/>'
        '<Attribute name="hsapiens_paralog_perc_id"/>'
        '</Dataset></Query>'
    )

    encoded_query = urllib.parse.quote(query)
    url = f"https://www.ensembl.org/biomart/martservice?query={encoded_query}"
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req, timeout=300) as resp:
        text = resp.read().decode("utf-8")

    df = pd.read_csv(StringIO(text), sep="\t")
    df.columns = ["gene", "paralog_gene", "paralog_perc_id"]
    df = df.dropna(subset=["gene", "paralog_gene"])
    df = df[df["gene"] != df["paralog_gene"]]
    df.to_csv(PARALOG_CACHE, index=False)
    print(f"Saved {len(df)} paralog pairs to {PARALOG_CACHE}")
    return df


def compute_paralog_features(
    paralogs: pd.DataFrame, expr: pd.DataFrame
) -> pd.DataFrame:
    paralog_counts = paralogs.groupby("gene")["paralog_gene"].nunique().reset_index()
    paralog_counts.columns = ["gene", "paralog_count"]

    gene_mean_tpm = {}
    gene_breadth = {}
    for col in expr.columns:
        gene = col.split(" ")[0]
        values = expr[col].dropna().values
        gene_mean_tpm[gene] = float(np.mean(values))
        gene_breadth[gene] = float(np.mean(values > np.log1p(EXPRESSED_TPM_THRESHOLD)))

    paralog_expr_breadth = []
    for gene, group in paralogs.groupby("gene"):
        paralog_genes = group["paralog_gene"].unique()
        breadths = [gene_breadth.get(pg, 0.0) for pg in paralog_genes]
        expressed_paralogs = [b for b in breadths if b > 0]

        paralog_expr_breadth.append({
            "gene": gene,
            "n_expressed_paralogs": len(expressed_paralogs),
            "max_paralog_breadth": max(breadths) if breadths else 0.0,
            "mean_paralog_breadth": float(np.mean(breadths)) if breadths else 0.0,
        })

    breadth_df = pd.DataFrame(paralog_expr_breadth)

    result = paralog_counts.merge(breadth_df, on="gene", how="left")

    mean_tpm_df = pd.DataFrame([
        {"gene": g, "expr_mean": v} for g, v in gene_mean_tpm.items()
    ])
    result = result.merge(mean_tpm_df, on="gene", how="left")

    return result


def main():
    ts = lambda: datetime.now().strftime("%H:%M:%S")
    print(f"[{ts()}] Experiment E5: Paralog buffering vs CRISPR DI")

    crispr_di = pd.read_csv(CRISPR_DI_PATH)
    print(f"[{ts()}] CRISPR DI: {len(crispr_di)} genes")

    if not EXPRESSION_CACHE.exists():
        print(f"[{ts()}] ERROR: Expression data not found at {EXPRESSION_CACHE}")
        print("Run script 09 first to download expression data.")
        return

    expr = pd.read_csv(EXPRESSION_CACHE, index_col=0)
    print(f"[{ts()}] Expression: {expr.shape[0]} cell lines x {expr.shape[1]} genes")

    paralogs = fetch_ensembl_paralogs()
    print(f"[{ts()}] Paralogs: {len(paralogs)} pairs")

    paralog_features = compute_paralog_features(paralogs, expr)
    print(f"[{ts()}] Paralog features: {len(paralog_features)} genes")

    merged = crispr_di.merge(paralog_features, on="gene")
    print(f"[{ts()}] Overlap (genes with DI + paralog data): {len(merged)}")

    n_before = len(merged)
    merged = merged.dropna(subset=["di_corrected", "paralog_count", "n_contexts", "expr_mean"])
    merged = merged[np.isfinite(merged["expr_mean"].values) & np.isfinite(merged["di_corrected"].values)]
    if n_before != len(merged):
        print(f"[{ts()}] Dropped {n_before - len(merged)} rows with NaN/inf covariates")

    if len(merged) == 0:
        save_results({"experiment": "E5_paralog_buffering", "n": 0}, RESULTS_PATH)
        return

    # --- E5a: paralog count vs DI ---
    di_vals = merged["di_corrected"].values.astype(np.float64)
    pc_vals = merged["paralog_count"].values.astype(np.float64)
    n_ctx = merged["n_contexts"].values.astype(np.float64)
    expr_mean = merged["expr_mean"].values.astype(np.float64)
    Z = np.column_stack([n_ctx, expr_mean])
    n = len(merged)

    raw_rho_a, raw_p_a = stats.spearmanr(di_vals, pc_vals)
    partial_rho_a = partial_spearman(di_vals, pc_vals, Z)
    ci_low_a = ci_lower_bound(partial_rho_a, n, n_cov=2)
    ci_high_a = ci_upper_bound(partial_rho_a, n, n_cov=2)
    p_val_a = fisher_z_pvalue(partial_rho_a, n, n_cov=2)
    passes_a = ci_low_a > CI_FLOOR and p_val_a < ALPHA

    print(f"[{ts()}] E5a: paralog_count vs DI")
    print(f"  Raw rho = {raw_rho_a:.4f}, Partial rho = {partial_rho_a:.4f}")
    print(f"  CI = [{ci_low_a:.4f}, {ci_high_a:.4f}], p = {p_val_a:.2e}")
    print(f"  {'PASS' if passes_a else 'FAIL'}")

    # --- E5b: paralog expression breadth vs DI (subset with expressed paralogs) ---
    has_expressed = merged[merged["n_expressed_paralogs"] > 0]
    n_b = len(has_expressed)
    print(f"[{ts()}] E5b: {n_b} genes with >= 1 expressed paralog")

    if n_b > 30:
        di_b = has_expressed["di_corrected"].values.astype(np.float64)
        breadth_b = has_expressed["max_paralog_breadth"].values.astype(np.float64)
        n_ctx_b = has_expressed["n_contexts"].values.astype(np.float64)
        expr_b = has_expressed["expr_mean"].values.astype(np.float64)
        Z_b = np.column_stack([n_ctx_b, expr_b])

        raw_rho_b, raw_p_b = stats.spearmanr(di_b, breadth_b)
        partial_rho_b = partial_spearman(di_b, breadth_b, Z_b)
        ci_low_b = ci_lower_bound(partial_rho_b, n_b, n_cov=2)
        ci_high_b = ci_upper_bound(partial_rho_b, n_b, n_cov=2)
        p_val_b = fisher_z_pvalue(partial_rho_b, n_b, n_cov=2)
        passes_b = ci_high_b < -CI_FLOOR and p_val_b < ALPHA

        print(f"  Raw rho = {raw_rho_b:.4f}, Partial rho = {partial_rho_b:.4f}")
        print(f"  CI = [{ci_low_b:.4f}, {ci_high_b:.4f}], p = {p_val_b:.2e}")
        print(f"  {'PASS' if passes_b else 'FAIL'} (expected negative)")

        result_b = {
            "n": n_b,
            "raw_spearman_rho": round(raw_rho_b, 4),
            "raw_p_value": raw_p_b,
            "partial_spearman_rho": round(partial_rho_b, 4),
            "partial_ci_low": round(ci_low_b, 4),
            "partial_ci_high": round(ci_high_b, 4),
            "partial_p_value": p_val_b,
            "criterion": f"CI upper < -{CI_FLOOR} AND p < {ALPHA}",
            "passes": passes_b,
        }
    else:
        result_b = {"n": n_b, "note": "too few genes with expressed paralogs"}

    result = {
        "experiment": "E5_paralog_buffering",
        "method": "Pearson-residualize then Spearman, 2 covariates (n_contexts, expr_mean)",
        "alpha": ALPHA,
        "E5a_paralog_count": {
            "n": n,
            "raw_spearman_rho": round(raw_rho_a, 4),
            "raw_p_value": raw_p_a,
            "partial_spearman_rho": round(partial_rho_a, 4),
            "partial_ci_low": round(ci_low_a, 4),
            "partial_ci_high": round(ci_high_a, 4),
            "partial_p_value": p_val_a,
            "criterion": f"CI lower > {CI_FLOOR} AND p < {ALPHA}",
            "passes": passes_a,
            "label": "PRIMARY",
        },
        "E5b_paralog_breadth": result_b,
    }

    save_results(result, RESULTS_PATH)


if __name__ == "__main__":
    main()
