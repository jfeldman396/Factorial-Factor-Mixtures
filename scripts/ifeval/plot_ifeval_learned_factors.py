#!/usr/bin/env python3
"""Visualize the learned IFEval mixture-factor model.

Creates:
  - pairwise factor scatter plots colored by MAP profile
  - interactive 3D factor scatter plot
  - marginal factor histograms with mixture group centers
  - factor-score and group heatmaps by LLM ordered by accuracy
"""

from __future__ import annotations

from pathlib import Path
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px


FIT_DIR = Path(os.environ.get(
    "FIT_DIR",
    Path(__file__).resolve().parents[1] / "results" / "reproduced_openeval_ifeval_H3_G3_interpretation",
))
OUT_DIR = Path(os.environ.get(
    "OUT_DIR",
    Path(__file__).resolve().parents[1] / "results" / "ifeval_3d_factor_visualizations",
))


def profile_palette(profiles: list[str]) -> dict[str, str]:
    colors = [
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
        "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
        "#4c78a8", "#f58518", "#54a24b", "#e45756", "#72b7b2",
        "#b279a2", "#ff9da6", "#9d755d", "#bab0ac", "#2f4b7c",
        "#a05195", "#f95d6a", "#ffa600",
    ]
    return {p: colors[i % len(colors)] for i, p in enumerate(sorted(profiles))}


def load_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    scores = pd.read_csv(FIT_DIR / "openeval_model_factor_scores_profiles.csv")
    groups = pd.read_csv(FIT_DIR / "openeval_factor_mixture_groups.csv")
    scores["profile_id"] = scores["profile_id"].astype(str)
    return scores, groups


def plot_pairwise(scores: pd.DataFrame, colors: dict[str, str]) -> None:
    pairs = [(1, 2), (1, 3), (2, 3)]
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.5), dpi=180)

    for ax, (a, b) in zip(axes, pairs):
        for profile, sub in scores.groupby("profile_id", sort=True):
            ax.scatter(
                sub[f"factor_{a}"],
                sub[f"factor_{b}"],
                s=18 + 80 * sub["accuracy"],
                color=colors[profile],
                alpha=0.82,
                edgecolor="white",
                linewidth=0.4,
                label=profile,
            )
        ax.axhline(0, color="0.78", lw=0.8, ls=":")
        ax.axvline(0, color="0.78", lw=0.8, ls=":")
        ax.set_xlabel(f"F{a}")
        ax.set_ylabel(f"F{b}")
        ax.set_title(f"F{a} vs F{b}")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, title="profile", loc="lower center", ncol=7, frameon=False, fontsize=7)
    fig.suptitle("IFEval learned mixture factors, colored by MAP profile", fontsize=14, fontweight="bold")
    fig.subplots_adjust(bottom=0.28, top=0.84, wspace=0.28)
    fig.savefig(OUT_DIR / "ifeval_pairwise_factor_scatter_profiles.png", bbox_inches="tight")
    plt.close(fig)


def plot_marginals(scores: pd.DataFrame, groups: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.7), dpi=180)
    for h, ax in enumerate(axes, start=1):
        x = scores[f"factor_{h}"]
        ax.hist(x, bins=22, color="#d8dee9", edgecolor="white")
        g = groups[groups["factor"].eq(f"F{h}")]
        if g.empty and "factor" in groups.columns:
            g = groups[groups["factor"].eq(h)]
        for _, row in g.iterrows():
            ax.axvline(row["mean"], color="#b23a48", lw=1.8, ls="--")
            ax.text(row["mean"], ax.get_ylim()[1] * 0.92, f"g{int(row['group'])}", rotation=90, va="top", ha="right", fontsize=7)
        ax.set_title(f"F{h} marginal")
        ax.set_xlabel("factor score")
        ax.set_ylabel("LLM count")
    fig.suptitle("IFEval marginal factor distributions with mixture centers", fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ifeval_factor_marginals_with_mixture_centers.png", bbox_inches="tight")
    plt.close(fig)


def plot_heatmaps(scores: pd.DataFrame) -> None:
    ordered = scores.sort_values("accuracy", ascending=False).reset_index(drop=True)
    factor_cols = [f"factor_{h}" for h in range(1, 4)]
    group_cols = [f"group_factor_{h}" for h in range(1, 4)]

    # Factor scores heatmap.
    mat = ordered[factor_cols].to_numpy().T
    vmax = np.nanmax(np.abs(mat))
    fig, ax = plt.subplots(figsize=(18, 4.8), dpi=180)
    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_yticks(range(3), labels=["F1", "F2", "F3"])
    ax.set_xticks(range(len(ordered)), labels=ordered["model_id"], rotation=90, fontsize=4)
    ax.set_title("IFEval factor scores by LLM, ordered by accuracy")
    fig.colorbar(im, ax=ax, fraction=0.018, pad=0.01, label="factor score")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ifeval_factor_score_heatmap_by_llm.png", bbox_inches="tight")
    plt.close(fig)

    # MAP group heatmap.
    gmat = ordered[group_cols].to_numpy().T
    fig, ax = plt.subplots(figsize=(18, 4.8), dpi=180)
    im = ax.imshow(gmat, aspect="auto", cmap="viridis", vmin=1, vmax=3)
    ax.set_yticks(range(3), labels=["F1 group", "F2 group", "F3 group"])
    ax.set_xticks(range(len(ordered)), labels=ordered["model_id"], rotation=90, fontsize=4)
    ax.set_title("IFEval MAP mixture groups by LLM, ordered by accuracy")
    cbar = fig.colorbar(im, ax=ax, fraction=0.018, pad=0.01)
    cbar.set_ticks([1, 2, 3])
    cbar.set_label("MAP group")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ifeval_map_group_heatmap_by_llm.png", bbox_inches="tight")
    plt.close(fig)


def plot_3d(scores: pd.DataFrame, colors: dict[str, str]) -> None:
    color_sequence = [colors[p] for p in sorted(colors)]
    hover = ["model_id", "accuracy", "profile_id", "group_factor_1", "group_factor_2", "group_factor_3"]
    fig = px.scatter_3d(
        scores,
        x="factor_1",
        y="factor_2",
        z="factor_3",
        color="profile_id",
        size="accuracy",
        hover_data=hover,
        color_discrete_sequence=color_sequence,
        title="IFEval learned mixture factors, colored by MAP profile",
    )
    fig.update_traces(marker=dict(opacity=0.84, line=dict(width=0.5, color="white")))
    fig.write_html(OUT_DIR / "ifeval_3d_factor_scatter_profiles_interactive.html")

    # Static 3D snapshot.
    fig_static = plt.figure(figsize=(9, 8), dpi=180)
    ax = fig_static.add_subplot(111, projection="3d")
    for profile, sub in scores.groupby("profile_id", sort=True):
        ax.scatter(
            sub["factor_1"],
            sub["factor_2"],
            sub["factor_3"],
            s=20 + 80 * sub["accuracy"],
            color=colors[profile],
            alpha=0.82,
            edgecolor="white",
            linewidth=0.4,
            label=profile,
        )
    ax.set_xlabel("F1")
    ax.set_ylabel("F2")
    ax.set_zlabel("F3")
    ax.set_title("IFEval learned mixture factors")
    ax.view_init(elev=22, azim=38)
    ax.legend(title="profile", loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=6, frameon=False)
    fig_static.savefig(OUT_DIR / "ifeval_3d_factor_scatter_profiles.png", bbox_inches="tight")
    plt.close(fig_static)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    scores, groups = load_data()
    colors = profile_palette(scores["profile_id"].unique().tolist())
    plot_pairwise(scores, colors)
    plot_marginals(scores, groups)
    plot_heatmaps(scores)
    plot_3d(scores, colors)

    profile_summary = (
        scores.groupby("profile_id")
        .agg(
            n_models=("model_id", "size"),
            mean_accuracy=("accuracy", "mean"),
            min_accuracy=("accuracy", "min"),
            max_accuracy=("accuracy", "max"),
        )
        .reset_index()
        .sort_values("mean_accuracy", ascending=False)
    )
    profile_summary.to_csv(OUT_DIR / "ifeval_visual_profile_summary.csv", index=False)
    print("Saved visualizations in:", OUT_DIR)
    print(profile_summary.to_string(index=False))


if __name__ == "__main__":
    main()
