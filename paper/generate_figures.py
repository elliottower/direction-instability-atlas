"""Generate figures for the DI atlas paper (v4 Patterns).

Figures:
  fig_di_distributions.pdf      — DI distributions per dataset (3-panel)
  fig_forest_plot.pdf            — forest plot of all Batches 1-2A results
  fig_essentiality_scatter.pdf  — DI vs essentiality breadth (sign-flip)
  fig_ko_vs_oe_scatter.pdf      — knockout DI vs overexpression DI (null)
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

RESULTS = Path(__file__).parent.parent / "results"
FIGURES = Path(__file__).parent / "figures"
FIGURES.mkdir(exist_ok=True)

plt.rcParams.update({
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "figure.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.1,
})


def load_depmap_essentiality():
    depmap_path = RESULTS.parent / "data" / "depmap" / "CRISPRGeneEffect.csv"
    if not depmap_path.exists():
        depmap_path = RESULTS.parent / "data" / "depmap" / "24q4" / "CRISPRGeneEffect.csv"
    depmap = pd.read_csv(depmap_path, index_col=0)
    breadth = {}
    n_cl = depmap.shape[0]
    for col in depmap.columns:
        gene = col.split(" ")[0]
        breadth[gene] = (depmap[col] < -0.5).sum() / n_cl
    return pd.Series(breadth, name="essentiality_breadth")


def fig_distributions():
    print("Generating DI distributions...")
    crispr = pd.read_csv(RESULTS / "crispr_di.csv")
    orf = pd.read_csv(RESULTS / "orf_di.csv")
    tahoe = pd.read_csv(RESULTS / "tahoe_di.csv")

    datasets = [
        ("Tahoe-100M\n(379 drugs, cell lines)", tahoe["di_corrected"].dropna(), "#1b9e77"),
        ("JUMP-CP CRISPR\n(7,946 genes, plate reps)", crispr["di_corrected"].dropna(), "#d95f02"),
        ("JUMP-CP ORF\n(12,590 genes, plate reps)", orf["di_corrected"].dropna(), "#7570b3"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(11, 3.5), sharey=False)

    for ax, (label, data, color) in zip(axes, datasets):
        ax.hist(data.values, bins=50, color=color, alpha=0.75, edgecolor="white", linewidth=0.3)
        ax.set_xlabel("DI (corrected)")
        ax.set_title(label, fontsize=9, linespacing=1.3)
        ax.axvline(data.median(), color="black", linestyle="--", linewidth=0.8, alpha=0.6,
                   label=f"median = {data.median():.2f}")
        ax.axvline(1.0, color="gray", linestyle=":", linewidth=0.6, alpha=0.4)
        ax.text(0.97, 0.93, f"median = {data.median():.2f}\nn = {len(data):,}",
                transform=ax.transAxes, ha="right", va="top", fontsize=7.5)

    axes[0].set_ylabel("Count")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig_di_distributions.pdf")
    fig.savefig(FIGURES / "fig_di_distributions.png")
    plt.close(fig)
    print("  Saved.")


def fig_forest():
    print("Generating forest plot...")

    experiments = [
        {"label": "Cross-modal concordance", "rho": 0.190, "ci_lo": 0.090, "ci_hi": 0.285,
         "n": 378, "pass": True, "batch": "1"},
        {"label": "Same-modality concordance", "rho": 0.254, "ci_lo": 0.027, "ci_hi": 0.456,
         "n": 76, "pass": False, "batch": "1"},

        {"label": "Essentiality vs DI", "rho": -0.332, "ci_lo": -0.352, "ci_hi": -0.312,
         "n": 7653, "pass": False, "batch": "1"},
        {"label": "Knockout ≈ overexpression DI", "rho": 0.011, "ci_lo": -0.016, "ci_hi": 0.038,
         "n": 5220, "pass": False, "batch": "1"},

        {"label": "Expression variance", "rho": 0.060, "ci_lo": 0.038, "ci_hi": 0.082,
         "n": 7802, "pass": False, "batch": "2A"},
        {"label": "Paralog count", "rho": -0.007, "ci_lo": -0.031, "ci_hi": 0.017,
         "n": 6894, "pass": False, "batch": "2A"},
        {"label": "Paralog breadth", "rho": -0.029, "ci_lo": -0.053, "ci_hi": -0.006,
         "n": 6882, "pass": False, "batch": "2A"},
        {"label": "Drug selectivity", "rho": 0.144, "ci_lo": -0.088, "ci_hi": 0.360,
         "n": 75, "pass": False, "batch": "2A"},

        {"label": "Viability concordance (exploratory)", "rho": 0.359, "ci_lo": 0.144, "ci_hi": 0.542,
         "n": 75, "pass": None, "batch": "exp"},
    ]

    fig, ax = plt.subplots(figsize=(7, 5.5))

    y_positions = []
    y = 0
    prev_batch = None
    for exp in experiments:
        if prev_batch is not None and exp["batch"] != prev_batch:
            y -= 0.5
        y_positions.append(y)
        prev_batch = exp["batch"]
        y -= 1

    for i, (exp, yp) in enumerate(zip(experiments, y_positions)):
        if exp["pass"] is True:
            color = "#2ca02c"
            style = "-"
        elif exp["pass"] is None:
            color = "#666666"
            style = "--"
        else:
            color = "#1a1a1a"
            style = "-"

        ax.plot([exp["ci_lo"], exp["ci_hi"]], [yp, yp], color=color, linewidth=1.5,
                linestyle=style, solid_capstyle="butt")
        ax.plot(exp["rho"], yp, "o", color=color, markersize=5, zorder=5)

        ax.text(-0.52, yp, exp["label"], ha="right", va="center", fontsize=9)
        ax.text(0.72, yp, f"n = {exp['n']:,}", ha="left", va="center", fontsize=8,
                color="#555555")

    ax.axvspan(-0.05, 0.05, color="#e0e0e0", alpha=0.5, zorder=0)
    ax.axvline(0, color="gray", linewidth=0.5, alpha=0.5)

    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker="o", color="#2ca02c", markersize=5, linewidth=1.5,
               label="Passes criteria"),
        Line2D([0], [0], marker="o", color="#1a1a1a", markersize=5, linewidth=1.5,
               label="Fails criteria"),
        Line2D([0], [0], marker="o", color="#666666", markersize=5, linewidth=1.5,
               linestyle="--", label="Exploratory"),
        plt.Rectangle((0, 0), 1, 1, fc="#e0e0e0", alpha=0.5, label="±0.05 CI floor"),
    ]
    ax.legend(handles=legend_elements, loc="upper left", fontsize=8,
              framealpha=0.9, edgecolor="#cccccc")

    ax.set_xlabel("Partial Spearman ρ", fontsize=11)
    ax.set_xlim(-0.52, 0.72)
    ax.set_ylim(min(y_positions) - 0.5, max(y_positions) + 0.8)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.savefig(FIGURES / "fig_forest_plot.pdf")
    fig.savefig(FIGURES / "fig_forest_plot.png")
    plt.close(fig)
    print("  Saved.")


def fig_essentiality():
    print("Generating essentiality scatter...")
    crispr = pd.read_csv(RESULTS / "crispr_di.csv")
    ess = load_depmap_essentiality()

    merged = crispr.set_index("gene").join(ess, how="inner")
    merged = merged.dropna(subset=["di_corrected", "essentiality_breadth"])

    x = merged["essentiality_breadth"].values
    y = merged["di_corrected"].values

    rho, p = stats.spearmanr(x, y)

    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    ax.scatter(x, y, s=1.5, alpha=0.15, color="#2c7fb8", rasterized=True)

    z = np.polyfit(x, y, 1)
    xline = np.linspace(0, x.max(), 100)
    ax.plot(xline, np.polyval(z, xline), color="#d95f02", linewidth=1.5, zorder=5)

    ax.set_xlabel("Essentiality breadth\n(fraction of cell lines, Chronos < −0.5)", fontsize=10)
    ax.set_ylabel("Direction instability (corrected)")
    # Title shows partial Spearman (controlling for n_contexts) to match paper text
    with open(RESULTS / "essentiality.json") as f:
        ess_result = json.load(f)
    partial_rho = ess_result["partial_spearman_rho"]
    partial_n = ess_result["n"]
    ax.set_title(f"partial ρ = {partial_rho:.3f}, n = {partial_n:,}", fontsize=10, style="italic")
    ax.text(0.97, 0.97, "Pre-registered prediction: positive\nObserved: negative (sign-flip)",
            transform=ax.transAxes, ha="right", va="top", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.5))

    fig.savefig(FIGURES / "fig_essentiality_scatter.pdf")
    fig.savefig(FIGURES / "fig_essentiality_scatter.png")
    plt.close(fig)
    print(f"  Saved. rho={rho:.4f}, n={len(x)}")


def fig_ko_vs_oe():
    print("Generating KO vs OE scatter...")
    crispr = pd.read_csv(RESULTS / "crispr_di.csv").set_index("gene")
    orf = pd.read_csv(RESULTS / "orf_di.csv").set_index("gene")

    shared = crispr.index.intersection(orf.index)
    ko_di = crispr.loc[shared, "di_corrected"]
    oe_di = orf.loc[shared, "di_corrected"]

    both = pd.DataFrame({"ko": ko_di, "oe": oe_di}).dropna()
    rho, p = stats.spearmanr(both["ko"], both["oe"])

    clip_hi = both[["ko", "oe"]].quantile(0.995).max() * 1.05
    clip_lo = max(0, both[["ko", "oe"]].quantile(0.005).min() * 0.95)

    fig, ax = plt.subplots(figsize=(5, 4.5))
    ax.scatter(both["ko"], both["oe"], s=1.5, alpha=0.12, color="#7570b3", rasterized=True)

    z = np.polyfit(both["ko"], both["oe"], 1)
    xline = np.linspace(clip_lo, clip_hi, 100)
    ax.plot(xline, np.polyval(z, xline), color="#d95f02", linewidth=1.2, zorder=5,
            alpha=0.7)

    ax.set_xlim(clip_lo, clip_hi)
    ax.set_ylim(clip_lo, clip_hi)

    ax.set_xlabel("Knockout DI (JUMP-CP CRISPR, corrected)")
    ax.set_ylabel("Overexpression DI (JUMP-CP ORF, corrected)")
    ax.set_title(f"ρ = {rho:.3f}, n = {len(both):,}", fontsize=10, style="italic")
    ax.text(0.03, 0.97, "DI is mechanism-specific:\nKO and OE produce\nunrelated DI values",
            transform=ax.transAxes, ha="left", va="top", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#e0e0e0", alpha=0.5))

    fig.savefig(FIGURES / "fig_ko_vs_oe_scatter.pdf")
    fig.savefig(FIGURES / "fig_ko_vs_oe_scatter.png")
    plt.close(fig)
    print(f"  Saved. rho={rho:.4f}, n={len(both)}")


if __name__ == "__main__":
    fig_distributions()
    fig_forest()
    fig_essentiality()
    fig_ko_vs_oe()
    print("Done.")
