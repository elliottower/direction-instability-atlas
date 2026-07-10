"""Graphical abstract for the DI atlas paper (Patterns format).

Option C V3b: Two-panel biological diagram with shades of blue/orange per cell,
expression profile line plots showing overlap vs divergence.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np
from pathlib import Path

FIGURES = Path(__file__).parent / "figures"
FIGURES.mkdir(exist_ok=True)

C_LOW = "#2c7fb8"
C_HIGH = "#d95f02"
C_TEXT = "#1a1a1a"
C_SUB = "#555555"
C_LIGHT_BLUE = "#d0e4f0"
C_LIGHT_ORANGE = "#fde0c5"

BLUES = ["#08519c", "#3182bd", "#6baed6"]
ORANGES = ["#d94701", "#f16913", "#fdae6b"]


def draw_cell(ax, cx, cy, w, h, color, label, alpha=0.25):
    cell = FancyBboxPatch((cx - w/2, cy - h/2), w, h,
                           boxstyle="round,pad=0.02",
                           facecolor=color, edgecolor=mpatches.colors.to_rgba(color, 0.6),
                           alpha=alpha, linewidth=1.5, zorder=2)
    ax.add_patch(cell)
    ax.text(cx, cy, label, ha="center", va="center", fontsize=7.5,
            color=C_TEXT, fontweight="bold", zorder=3)


def main():
    fig = plt.figure(figsize=(10, 5.5), dpi=250)
    fig.patch.set_facecolor("white")

    fig.text(0.50, 0.96, "Which perturbations generalize across cell types?",
             ha="center", va="top", fontsize=15, fontweight="bold", color=C_TEXT)
    fig.text(0.50, 0.905, "~55,000 perturbations  ·  3 platforms  ·  11 pre-registered tests",
             ha="center", va="top", fontsize=7.5, color=C_SUB)

    # ── Left panel: GENERALIZES ──
    ax_l = fig.add_axes([0.03, 0.10, 0.44, 0.76])
    ax_l.set_xlim(0, 1)
    ax_l.set_ylim(0, 1)
    ax_l.axis("off")

    bg_l = FancyBboxPatch((0.02, 0.02), 0.96, 0.96,
                           boxstyle="round,pad=0.02",
                           facecolor=C_LIGHT_BLUE, edgecolor=C_LOW,
                           alpha=0.15, linewidth=1.5, zorder=0)
    ax_l.add_patch(bg_l)

    ax_l.text(0.50, 0.93, "GENERALIZES", ha="center", va="center",
              fontsize=13, fontweight="bold", color=C_LOW)
    ax_l.text(0.50, 0.86, "low Direction Instability", ha="center", va="center",
              fontsize=9, color=C_LOW, alpha=0.8)

    pill = FancyBboxPatch((0.25, 0.73), 0.50, 0.08,
                           boxstyle="round,pad=0.015",
                           facecolor=C_LOW, edgecolor="none", alpha=0.85, zorder=3)
    ax_l.add_patch(pill)
    ax_l.text(0.50, 0.77, "mTOR inhibitor", ha="center", va="center",
              fontsize=9, fontweight="bold", color="white", zorder=4)

    cell_xs = [0.25, 0.50, 0.75]
    arrow_origins = [0.30, 0.50, 0.70]
    for i, (cx, ox) in enumerate(zip(cell_xs, arrow_origins)):
        ax_l.annotate("", xy=(cx, 0.62), xytext=(ox, 0.73),
                      arrowprops=dict(arrowstyle="->,head_width=0.06", color=BLUES[i], lw=1.5))
        draw_cell(ax_l, cx, 0.56, 0.16, 0.09, BLUES[i], f"Cell {i+1}")

    genes = np.arange(7)
    base_profile = np.array([0.8, 0.35, 0.65, 0.9, 0.4, 0.75, 0.5])

    ax_expr_l = fig.add_axes([0.10, 0.22, 0.30, 0.18])
    for i, offset in enumerate([0.02, -0.01, 0.03]):
        profile = base_profile + offset
        ax_expr_l.plot(genes, profile, "-o", color=BLUES[i], alpha=0.75, markersize=4, lw=2.0,
                       label=f"Cell {i+1}")
    ax_expr_l.set_ylim(0, 1.1)
    ax_expr_l.set_xticks([])
    ax_expr_l.set_yticks([])
    ax_expr_l.set_xlabel("genes", fontsize=9, color=C_SUB)
    ax_expr_l.set_ylabel("expression", fontsize=9, color=C_SUB)
    ax_expr_l.spines["top"].set_visible(False)
    ax_expr_l.spines["right"].set_visible(False)
    ax_expr_l.spines["left"].set_color("#cccccc")
    ax_expr_l.spines["bottom"].set_color("#cccccc")
    ax_expr_l.set_title("Response profiles overlap", fontsize=7.5, color=C_LOW, fontweight="bold")
    ax_expr_l.legend(fontsize=6, loc="lower right", framealpha=0.7, edgecolor="none")
    ax_expr_l.patch.set_alpha(0)

    ax_l.text(0.50, 0.04, "mTOR inhibitors  ·  HDAC inhibitors\nEssential genes",
              ha="center", va="center", fontsize=7.5, color=C_LOW, linespacing=1.5,
              fontweight="bold")

    # ── Right panel: CONTEXT-DEPENDENT ──
    ax_r = fig.add_axes([0.53, 0.10, 0.44, 0.76])
    ax_r.set_xlim(0, 1)
    ax_r.set_ylim(0, 1)
    ax_r.axis("off")

    bg_r = FancyBboxPatch((0.02, 0.02), 0.96, 0.96,
                           boxstyle="round,pad=0.02",
                           facecolor=C_LIGHT_ORANGE, edgecolor=C_HIGH,
                           alpha=0.15, linewidth=1.5, zorder=0)
    ax_r.add_patch(bg_r)

    ax_r.text(0.50, 0.93, "CONTEXT-DEPENDENT", ha="center", va="center",
              fontsize=13, fontweight="bold", color=C_HIGH)
    ax_r.text(0.50, 0.86, "high Direction Instability", ha="center", va="center",
              fontsize=9, color=C_HIGH, alpha=0.8)

    pill_r = FancyBboxPatch((0.20, 0.73), 0.60, 0.08,
                             boxstyle="round,pad=0.015",
                             facecolor=C_HIGH, edgecolor="none", alpha=0.85, zorder=3)
    ax_r.add_patch(pill_r)
    ax_r.text(0.50, 0.77, "Serotonin agonist", ha="center", va="center",
              fontsize=9, fontweight="bold", color="white", zorder=4)

    for i, (cx, ox) in enumerate(zip(cell_xs, arrow_origins)):
        ax_r.annotate("", xy=(cx, 0.62), xytext=(ox, 0.73),
                      arrowprops=dict(arrowstyle="->,head_width=0.06", color=ORANGES[i], lw=1.5))
        draw_cell(ax_r, cx, 0.56, 0.16, 0.09, ORANGES[i], f"Cell {i+1}")

    profiles_diff = [
        base_profile,
        np.array([0.2, 0.85, 0.3, 0.15, 0.9, 0.2, 0.8]),
        np.array([0.5, 0.5, 0.9, 0.3, 0.6, 0.4, 0.1]),
    ]
    ax_expr_r = fig.add_axes([0.60, 0.22, 0.30, 0.18])
    for i, profile in enumerate(profiles_diff):
        ax_expr_r.plot(genes, profile, "-o", color=ORANGES[i], alpha=0.75, markersize=4, lw=2.0,
                       label=f"Cell {i+1}")
    ax_expr_r.set_ylim(0, 1.1)
    ax_expr_r.set_xticks([])
    ax_expr_r.set_yticks([])
    ax_expr_r.set_xlabel("genes", fontsize=9, color=C_SUB)
    ax_expr_r.set_ylabel("expression", fontsize=9, color=C_SUB)
    ax_expr_r.spines["top"].set_visible(False)
    ax_expr_r.spines["right"].set_visible(False)
    ax_expr_r.spines["left"].set_color("#cccccc")
    ax_expr_r.spines["bottom"].set_color("#cccccc")
    ax_expr_r.set_title("Response profiles diverge", fontsize=7.5, color=C_HIGH, fontweight="bold")
    ax_expr_r.legend(fontsize=6, loc="lower right", framealpha=0.7, edgecolor="none")
    ax_expr_r.patch.set_alpha(0)

    ax_r.text(0.50, 0.04, "Serotonin agonists  ·  Histamine agonists\nSelectively essential genes",
              ha="center", va="center", fontsize=7.5, color=C_HIGH, linespacing=1.5,
              fontweight="bold")

    # Bottom summary
    fig.text(0.50, 0.03,
             "Mechanism of action structures generalizability (η² = 0.20)  ·  "
             "Essential genes generalize (ρ = −0.33)  ·  "
             "Concordant across modalities (ρ = 0.19)",
             ha="center", va="center", fontsize=7.5, color=C_SUB)

    fig.savefig(FIGURES / "graphical_abstract.pdf", bbox_inches="tight", pad_inches=0.08)
    fig.savefig(FIGURES / "graphical_abstract.png", bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    print("Saved graphical_abstract.pdf / .png")


if __name__ == "__main__":
    main()
