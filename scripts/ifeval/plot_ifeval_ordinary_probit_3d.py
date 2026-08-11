#!/usr/bin/env python3
"""Visualize IFEval ordinary probit factors in 3D.

The ordinary probit coordinates are aligned to the selected IFEval mixture
factor coordinates.  Points are colored by the mixture model MAP profile so the
ordinary factor geometry can be compared with the mixture profile structure.
"""

from __future__ import annotations

import os
from itertools import permutations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px


BUNDLE_ROOT = Path(os.environ.get("IFEVAL_BUNDLE_ROOT", Path(__file__).resolve().parents[1]))
MIXTURE_DIR = Path(os.environ.get(
    "MIXTURE_DIR",
    BUNDLE_ROOT / "results" / "reproduced_openeval_ifeval_H3_G3_interpretation",
))
ORDINARY_DIR = Path(os.environ.get(
    "ORDINARY_DIR",
    BUNDLE_ROOT / "results" / "reproduced_openeval_ifeval_ordinary_probit_H3_visualization",
))
OUT_DIR = Path(os.environ.get(
    "OUT_DIR",
    BUNDLE_ROOT / "results" / "reproduced_ifeval_ordinary_vs_mixture_factor_visualization",
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


def best_alignment(ord_mat: np.ndarray, mix_mat: np.ndarray) -> tuple[list[int], np.ndarray, np.ndarray]:
    cor = np.corrcoef(ord_mat.T, mix_mat.T)[:3, 3:]
    best_perm = None
    best_score = -np.inf
    for perm in permutations(range(3)):
        score = sum(abs(cor[perm[j], j]) for j in range(3))
        if score > best_score:
            best_score = score
            best_perm = list(perm)
    signs = np.sign([cor[best_perm[j], j] for j in range(3)])
    signs[signs == 0] = 1
    aligned = ord_mat[:, best_perm] * signs
    return best_perm, signs, aligned


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    mix = pd.read_csv(MIXTURE_DIR / "openeval_model_factor_scores_profiles.csv")
    ordinary = pd.read_csv(ORDINARY_DIR / "ordinary_probit_factor_scores.csv").rename(
        columns={"model": "model_id", "1": "ordinary_raw_F1", "2": "ordinary_raw_F2", "3": "ordinary_raw_F3"}
    )
    dat = mix.merge(ordinary, on=["model_id", "accuracy"], how="inner")

    ord_mat = dat[["ordinary_raw_F1", "ordinary_raw_F2", "ordinary_raw_F3"]].to_numpy()
    mix_mat = dat[["factor_1", "factor_2", "factor_3"]].to_numpy()
    perm, signs, aligned = best_alignment(ord_mat, mix_mat)
    for h in range(3):
        dat[f"ordinary_aligned_F{h + 1}"] = aligned[:, h]

    cor = np.corrcoef(aligned.T, mix_mat.T)[:3, 3:]
    align_df = pd.DataFrame(
        {
            "mixture_factor": [f"F{i}" for i in range(1, 4)],
            "ordinary_factor_used": [f"ordinary_raw_F{p + 1}" for p in perm],
            "sign_applied": signs.astype(int),
            "score_correlation": [cor[i, i] for i in range(3)],
        }
    )
    align_df.to_csv(OUT_DIR / "ifeval_ordinary_to_mixture_factor_alignment.csv", index=False)
    dat.to_csv(OUT_DIR / "ifeval_ordinary_aligned_factor_coordinates.csv", index=False)

    colors = profile_palette(dat["profile_id"].astype(str).unique().tolist())
    color_sequence = [colors[p] for p in sorted(colors)]
    hover = ["model_id", "accuracy", "profile_id", "group_factor_1", "group_factor_2", "group_factor_3"]

    fig = px.scatter_3d(
        dat,
        x="ordinary_aligned_F1",
        y="ordinary_aligned_F2",
        z="ordinary_aligned_F3",
        color="profile_id",
        size="accuracy",
        hover_data=hover,
        color_discrete_sequence=color_sequence,
        title="IFEval ordinary probit factors, aligned and colored by mixture profile",
    )
    fig.update_traces(marker=dict(opacity=0.84, line=dict(width=0.5, color="white")))
    fig.write_html(OUT_DIR / "ifeval_ordinary_probit_3d_profiles_interactive.html")

    # Side-by-side static snapshot.
    fig_static = plt.figure(figsize=(16, 7.2), dpi=180)
    for idx, (cols, title) in enumerate(
        [
            (("factor_1", "factor_2", "factor_3"), "Mixture factors"),
            (("ordinary_aligned_F1", "ordinary_aligned_F2", "ordinary_aligned_F3"), "Ordinary probit factors, aligned"),
        ],
        start=1,
    ):
        ax = fig_static.add_subplot(1, 2, idx, projection="3d")
        for profile, sub in dat.groupby("profile_id", sort=True):
            ax.scatter(
                sub[cols[0]], sub[cols[1]], sub[cols[2]],
                s=18 + 80 * sub["accuracy"],
                color=colors[str(profile)],
                alpha=0.82,
                edgecolor="white",
                linewidth=0.4,
                label=profile,
            )
        ax.set_xlabel("F1")
        ax.set_ylabel("F2")
        ax.set_zlabel("F3")
        ax.set_title(title)
        ax.view_init(elev=22, azim=38)
    handles, labels = fig_static.axes[0].get_legend_handles_labels()
    fig_static.legend(handles, labels, title="mixture profile", loc="lower center", ncol=7, frameon=False, fontsize=7)
    fig_static.suptitle("IFEval factor spaces colored by mixture MAP profile", fontsize=14, fontweight="bold")
    fig_static.subplots_adjust(bottom=0.22, top=0.88, wspace=0.05)
    fig_static.savefig(OUT_DIR / "ifeval_mixture_vs_ordinary_probit_3d_profiles.png", bbox_inches="tight")
    plt.close(fig_static)

    print("Saved outputs in:", OUT_DIR)
    print(align_df.to_string(index=False))


if __name__ == "__main__":
    main()
