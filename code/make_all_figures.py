#!/usr/bin/env python3
"""
Generate all manuscript figures for the A2A study.

Outputs (PDF for LaTeX, PNG for preview) into figures/:

  figure1_scatter.pdf              main text, full width  (figure*)
  figure2_prkn_decomposition.pdf   main text, one column  (figure)
  figureS1_threshold_heatmap.pdf   supplementary

Numbering note: the correlation forest plot has been dropped, and the
threshold heatmap moved to the supplement, so the main text carries two
figures. All values are computed from the results CSVs -- nothing is
hardcoded, so the figures cannot drift out of step with the tables.

Usage:
    python make_all_figures.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.stats import fisher_exact

# ------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------
PROJECT_DIR = Path("pathogenicity_abundance_project")
FIG_DIR = Path("figures")

PROTEINS = ["PTEN", "TPMT", "NUDT15", "ASPA", "PRKN", "VKORC1"]
ROLE = {"PTEN": "novel", "TPMT": "novel", "NUDT15": "novel",
        "ASPA": "control", "PRKN": "control", "VKORC1": "control"}
SCORE_COL = {p: "score" for p in PROTEINS}

PRED_T, EXP_T = 0.4, 0.8
PRED_GRID, EXP_GRID = [0.2, 0.3, 0.4], [0.7, 0.8, 0.9]

# Colourblind-safe: distinguished by hue *and* lightness, and legible in greyscale
GREY, ORANGE, NAVY, INK = "#9AA7B8", "#C1440E", "#1F4E8C", "#0B2545"
TEAL, PURPLE, BRICK, SLATE = "#1B7F5E", "#8C4A9E", "#B03A2E", "#6C7A89"

plt.rcParams["font.family"] = "DejaVu Sans"


def load(protein):
    """Load a protein's results, or return None with a warning if absent."""
    path = PROJECT_DIR / f"{protein}_final_results.csv"
    if not path.exists():
        print(f"  ! {path} not found - skipping {protein}")
        return None
    df = pd.read_csv(path)
    if "is_failure_case" not in df.columns:
        df["is_failure_case"] = (df.a2a_prediction < PRED_T) & (df[SCORE_COL[protein]] > EXP_T)
    return df


def fmt_p(p):
    if p == 0 or p < 1e-300:
        return r"$p < 10^{-300}$"
    e = int(np.floor(np.log10(p)))
    return rf"$p$ = {p/10**e:.1f} $\times$ 10$^{{{e}}}$"


def save(fig, stem):
    FIG_DIR.mkdir(exist_ok=True)
    for ext, kw in (("pdf", {}), ("png", {"dpi": 300})):
        fig.savefig(FIG_DIR / f"{stem}.{ext}", bbox_inches="tight",
                    facecolor="white", **kw)
    print(f"  saved figures/{stem}.pdf / .png")
    plt.close(fig)


# ==================================================================
# FIGURE 1 - scatter, six panels (main text, full width)
# ==================================================================
def make_figure1_scatter():
    # Sized for \textwidth in a two-column layout (~7in), so fonts render
    # near 1:1 rather than being shrunk to illegibility.
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.7))
    axes = axes.flatten()

    for ax, protein in zip(axes, PROTEINS):
        df = load(protein)
        if df is None:
            ax.set_visible(False)
            continue
        sc = SCORE_COL[protein]

        xmax = max(1.55, float(df[sc].max()) + 0.05)
        ax.add_patch(plt.Rectangle((EXP_T, 0), xmax - EXP_T, PRED_T,
                                   facecolor=ORANGE, alpha=0.09,
                                   edgecolor="none", zorder=0))

        other = df[~df.near_functional_site]
        site = df[df.near_functional_site]
        ax.scatter(other[sc], other.a2a_prediction, s=3, alpha=0.18,
                   color=GREY, linewidths=0, rasterized=True, zorder=2)
        ax.scatter(site[sc], site.a2a_prediction, s=9, alpha=0.85,
                   color=ORANGE, marker="D", linewidths=0, zorder=3)

        ax.axhline(PRED_T, color=NAVY, ls=":", lw=0.8, zorder=1)
        ax.axvline(EXP_T, color=NAVY, ls=":", lw=0.8, zorder=1)

        ax.set_title(f"{protein} ({ROLE[protein]})", fontsize=8.5,
                     fontweight="bold", color=INK, pad=4)
        ax.set_xlim(min(-0.3, float(df[sc].min()) - 0.05), xmax)
        ax.set_ylim(0, 1.0)
        ax.tick_params(labelsize=6.5, colors=INK)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)

    fig.supxlabel("Experimental abundance / stability score", fontsize=8.5, color=INK)
    fig.supylabel("A2A predicted abundance", fontsize=8.5, color=INK)

    fig.legend(handles=[
        plt.Line2D([], [], marker="o", ls="", color=GREY, ms=4, label="Other residues"),
        plt.Line2D([], [], marker="D", ls="", color=ORANGE, ms=4, label="Annotated functional site"),
        plt.Rectangle((0, 0), 1, 1, facecolor=ORANGE, alpha=0.15, label="Failure-case quadrant"),
    ], loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.06),
        frameon=False, fontsize=7.5)

    fig.tight_layout(rect=[0.02, 0.02, 1, 0.97])
    save(fig, "figure1_scatter")


# ==================================================================
# FIGURE 2 - PRKN decomposition (main text, single column)
# ==================================================================
def make_figure2_prkn():
    df = load("PRKN")
    if df is None:
        return

    CATALYTIC = [431]
    IBR = [332, 337, 352, 360, 365, 368, 373, 377]
    RING0 = [238, 241, 253, 257, 260, 263, 289, 293]
    RING2 = [418, 421, 436, 441, 446, 449, 457, 461]
    ALL_ZN = RING0 + IBR + RING2

    def stats(positions):
        site = df.position.isin(positions)
        ct = (pd.crosstab(df.is_failure_case, site)
                .reindex(index=[False, True], columns=[False, True], fill_value=0))
        table = ct.values.astype(float)
        orr, p = fisher_exact(ct.values)
        t = table + (0.5 if (table == 0).any() else 0.0)   # Haldane-Anscombe
        or_c = (t[1, 1] * t[0, 0]) / (t[1, 0] * t[0, 1])
        se = np.sqrt(sum(1.0 / x for x in t.flatten()))
        lo = 0.0 if (table == 0).any() else np.exp(np.log(or_c) - 1.96 * se)
        return dict(or_=orr, p=p, lo=lo, hi=np.exp(np.log(or_c) + 1.96 * se),
                    n=int(site.sum()),
                    n_fail=int(df.loc[site, "is_failure_case"].sum()))

    rows = [("Aggregate\n(25 pos.)", CATALYTIC + ALL_ZN, SLATE),
            ("Cys431\ncatalytic",    CATALYTIC,           TEAL),
            ("IBR Zn\n(8 pos.)",     IBR,                 PURPLE),
            ("All Zn\n(24 pos.)",    ALL_ZN,              BRICK)]
    res = [(lab, stats(pos), col) for lab, pos, col in rows]

    FLOOR = 0.03
    fig, ax = plt.subplots(figsize=(3.4, 3.6))
    ax.set_yscale("log")
    ax.set_ylim(FLOOR, 45)

    for i, (lab, s, col) in enumerate(res):
        if s["or_"] > 0:
            ax.bar(i, s["or_"], color=col, edgecolor="black", lw=0.7, width=0.6)
            ax.plot([i, i], [max(s["lo"], FLOOR), s["hi"]], color="black", lw=1.0, zorder=5)
            for y in (max(s["lo"], FLOOR), s["hi"]):
                ax.plot([i - 0.1, i + 0.1], [y, y], color="black", lw=1.0, zorder=5)
            ax.text(i, s["hi"] * 1.3, f"{s['or_']:.2f}", ha="center", va="bottom",
                    fontsize=8, fontweight="bold")
        else:
            # OR = 0 cannot be drawn on a log axis: show an open hatched bar
            # running to the upper 95% bound instead of a misleading floor value.
            ax.bar(i, s["hi"] - FLOOR, bottom=FLOOR, facecolor="white",
                   edgecolor=col, lw=1.2, hatch="///", width=0.6)
            ax.text(i, s["hi"] * 1.3, "0.00", ha="center", va="bottom",
                    fontsize=8, fontweight="bold", color=col)

    ax.axhline(1, color="gray", ls="--", lw=0.9)
    ax.text(len(res) - 0.55, 1.15, "OR = 1", fontsize=6.5, color="gray",
            ha="right", va="bottom")
    ax.set_xticks(range(len(res)))
    ax.set_xticklabels([r[0] for r in res], fontsize=7.5)
    ax.set_ylabel("Odds ratio (log scale)", fontsize=8.5)
    ax.tick_params(axis="y", labelsize=7.5)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    fig.tight_layout()
    save(fig, "figure2_prkn_decomposition")

    for lab, s, _ in res:
        print(f"    {lab.splitlines()[0]:<20} OR={s['or_']:.2f} "
              f"[{s['lo']:.2f}-{s['hi']:.2f}] p={s['p']:.2e} "
              f"n={s['n']} fail={s['n_fail']}")


# ==================================================================
# FIGURE S1 - threshold-sensitivity heatmaps (supplementary)
# ==================================================================
def make_figureS1_heatmap():
    grids, present = {}, []
    for protein in PROTEINS:
        path = PROJECT_DIR / f"{protein}_threshold_sensitivity.csv"
        if not path.exists():
            print(f"  ! {path} not found - skipping {protein}")
            continue
        sens = pd.read_csv(path)
        g = np.full((len(EXP_GRID), len(PRED_GRID)), np.nan)
        for i, e in enumerate(EXP_GRID):
            for j, pr in enumerate(PRED_GRID):
                row = sens[(sens.pred_threshold == pr) & (sens.exp_threshold == e)]
                if not row.empty:
                    g[i, j] = row.odds_ratio.values[0]
        grids[protein] = g
        present.append(protein)
    if not present:
        return

    allv = np.log2(np.concatenate([g.flatten() for g in grids.values()]))
    lim = np.nanmax(np.abs(allv))
    # Diverging colormap centred on log2(OR)=0 so that OR=1 renders neutral.
    # Without this, an asymmetric range (PRKN ~0.14, ASPA ~26) would colour
    # OR=1 as if it were already enriched.
    norm = TwoSlopeNorm(vmin=-lim, vcenter=0.0, vmax=lim)

    fig, axes = plt.subplots(2, 3, figsize=(7.0, 4.6))
    axes = axes.flatten()
    im = None
    for ax, protein in zip(axes, PROTEINS):
        if protein not in grids:
            ax.set_visible(False)
            continue
        g = grids[protein]
        im = ax.imshow(np.log2(g), cmap="RdBu_r", norm=norm, aspect="auto")
        for i in range(len(EXP_GRID)):
            for j in range(len(PRED_GRID)):
                v = g[i, j]
                shade = "white" if abs(np.log2(v)) > 0.62 * lim else "black"
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=7, color=shade)
        ax.set_xticks(range(len(PRED_GRID)))
        ax.set_xticklabels(PRED_GRID, fontsize=7)
        ax.set_yticks(range(len(EXP_GRID)))
        ax.set_yticklabels(EXP_GRID, fontsize=7)
        ax.set_title(f"{protein} ({ROLE[protein]})", fontsize=8.5,
                     fontweight="bold", color=INK, pad=4)

    fig.supxlabel("A2A prediction threshold", fontsize=8.5, color=INK)
    fig.supylabel("Experimental score threshold", fontsize=8.5, color=INK)
    cbar = fig.colorbar(im, ax=axes.tolist(), fraction=0.022, pad=0.02)
    cbar.set_label("log$_2$(odds ratio)", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    save(fig, "figureS1_threshold_heatmap")


if __name__ == "__main__":
    print("Figure 1 - scatter (six proteins, main text)")
    make_figure1_scatter()
    print("\nFigure 2 - PRKN decomposition (main text)")
    make_figure2_prkn()
    print("\nFigure S1 - threshold heatmaps (supplementary)")
    make_figureS1_heatmap()
    print("\nDone.")