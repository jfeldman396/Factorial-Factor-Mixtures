#!/usr/bin/env zsh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$SCRIPT_DIR"

mkdir -p "$ROOT/results"

echo "Step 1: singular-value shelf and Bai-Ng ICp2 rank diagnostics"
MATRIX_PATH="$ROOT/data/ifeval/openeval_ifeval_only_binary_matrix.csv" \
OUT_DIR="$ROOT/results/reproduced_openeval_ifeval_rank_diagnostics" \
H_MAX=10 \
Rscript plot_ifeval_rank_diagnostics.R

echo "Step 2: ordinary binary probit held-out CV over H=1:10"
MATRIX_PATH="$ROOT/data/ifeval/openeval_ifeval_only_binary_matrix.csv" \
OUT_DIR="$ROOT/results/openeval_ifeval_cv_H1_10_ordinary" \
H_GRID=1:10 \
K_FOLDS=3 \
WORKERS=8 \
AUG_ITER=4 \
REFINE_ITER=5 \
LAMBDA_L1=2 \
Rscript openeval_ordinary_probit_cv_H_selection.R

echo "Step 3: summarize held-out log likelihood, held-out BIC, and training BIC rank selection"
ORDINARY_H_SUMMARY="$ROOT/results/openeval_ifeval_cv_H1_10_ordinary/ordinary_probit_H_summary.csv" \
ORDINARY_FOLD_SCORES="$ROOT/results/openeval_ifeval_cv_H1_10_ordinary/ordinary_probit_factor_fold_scores.csv" \
OUT_DIR="$ROOT/results/reproduced_openeval_ifeval_rank_selection_comparison" \
Rscript summarize_openeval_rank_selection_cv.R

echo "Step 4: tune entrywise Lambda sparsity penalty at H=3, G=3"
MATRIX_PATH="$ROOT/data/ifeval/openeval_ifeval_only_binary_matrix.csv" \
OUT_DIR="$ROOT/results/reproduced_openeval_ifeval_lambda_sparsity_tuning" \
H_FIXED=3 \
G_FIXED=3 \
WORKERS=8 \
Rscript tune_ifeval_lambda_sparsity.R

selected_lambda_penalty=2
if [[ -f "$ROOT/results/reproduced_openeval_ifeval_lambda_sparsity_tuning/selected_lambda_penalty.txt" ]]; then
  selected_lambda_penalty="$(< "$ROOT/results/reproduced_openeval_ifeval_lambda_sparsity_tuning/selected_lambda_penalty.txt")"
fi

echo "Step 5: fit and interpret selected independent-mixture probit model at H=3, G=3, lambda_l1=$selected_lambda_penalty"
MATRIX_PATH="$ROOT/data/ifeval/openeval_ifeval_only_binary_matrix.csv" \
ITEM_METADATA_PATH="$ROOT/data/ifeval/openeval_item_metadata.csv" \
OUT_DIR="$ROOT/results/reproduced_openeval_ifeval_H3_G3_interpretation" \
H_FIXED=3 \
G_FIXED=3 \
WORKERS=8 \
REFINEMENT_LAMBDA_L1_PENALTY="$selected_lambda_penalty" \
Rscript fit_interpret_ifeval_H3_G3.R

echo "Step 6: fit and visualize selected ordinary binary probit factor model at H=3"
MMLU_MATRIX="$ROOT/data/ifeval/openeval_ifeval_only_binary_matrix.csv" \
OUT_DIR="$ROOT/results/reproduced_openeval_ifeval_ordinary_probit_H3_visualization" \
H=3 \
WORKERS=8 \
AUG_ITER=4 \
REFINE_ITER=10 \
LAMBDA_L1=2 \
Rscript fit_visualize_ordinary_probit_factors.R

echo "Step 7: summarize ordinary-vs-mixture factor alignment"
ORDINARY_DIR="$ROOT/results/reproduced_openeval_ifeval_ordinary_probit_H3_visualization" \
MIXTURE_DIR="$ROOT/results/reproduced_openeval_ifeval_H3_G3_interpretation" \
OUT_DIR="$ROOT/results/reproduced_openeval_ordinary_vs_mixture_H3" \
Rscript summarize_openeval_ordinary_probit_factors.R

echo "Step 8: plot ordinary-vs-mixture factor comparison in 3D"
ORDINARY_DIR="$ROOT/results/reproduced_openeval_ifeval_ordinary_probit_H3_visualization" \
MIXTURE_DIR="$ROOT/results/reproduced_openeval_ifeval_H3_G3_interpretation" \
OUT_DIR="$ROOT/results/reproduced_ifeval_ordinary_vs_mixture_factor_visualization" \
Rscript plot_ifeval_ordinary_vs_mixture_baseR.R

echo "Step 9: plot mixture lambda heatmaps"
FIT_DIR="$ROOT/results/reproduced_openeval_ifeval_H3_G3_interpretation" \
OUT_DIR="$ROOT/results/reproduced_openeval_ifeval_H3_G3_interpretation" \
Rscript plot_openeval_mixture_lambda_heatmaps.R

echo "Step 10: plot side-by-side ordinary and mixture lambda heatmaps"
MATRIX_PATH="$ROOT/data/ifeval/openeval_ifeval_only_binary_matrix.csv" \
MIXTURE_LOADINGS_PATH="$ROOT/results/reproduced_openeval_ifeval_H3_G3_interpretation/openeval_item_intercepts_loadings_metadata.csv" \
ORDINARY_LAMBDA_PATH="$ROOT/results/reproduced_openeval_ifeval_ordinary_probit_H3_visualization/ordinary_probit_lambda.csv" \
ORDINARY_MIXTURE_COR_PATH="$ROOT/results/reproduced_openeval_ordinary_vs_mixture_H3/ordinary_mixture_factor_correlation_matrix.csv" \
OUT_DIR="$ROOT/results/reproduced_openeval_ordinary_vs_mixture_H3" \
Rscript plot_openeval_side_by_side_lambdas.R

echo "Step 11: summarize factor-only and cross-loading items"
LOADINGS_PATH="$ROOT/results/reproduced_openeval_ifeval_H3_G3_interpretation/openeval_item_intercepts_loadings_metadata.csv" \
OUT_DIR="$ROOT/results/loadings_crossloadings" \
Rscript summarize_ifeval_loadings_crossloadings.R

echo "Step 12: refresh learned-mixture 3D factor visualizations"
FIT_DIR="$ROOT/results/reproduced_openeval_ifeval_H3_G3_interpretation" \
OUT_DIR="$ROOT/results/ifeval_3d_factor_visualizations" \
Rscript plot_ifeval_learned_factors_baseR.R

if python3 -c "import matplotlib, plotly" >/dev/null 2>&1; then
  echo "Step 13: refresh optional Python/Plotly learned-mixture 3D factor visualizations"
  FIT_DIR="$ROOT/results/reproduced_openeval_ifeval_H3_G3_interpretation" \
  OUT_DIR="$ROOT/results/ifeval_3d_factor_visualizations" \
  python3 plot_ifeval_learned_factors.py

  echo "Step 14: refresh optional Python/Plotly ordinary-vs-mixture 3D visualization"
  MIXTURE_DIR="$ROOT/results/reproduced_openeval_ifeval_H3_G3_interpretation" \
  ORDINARY_DIR="$ROOT/results/reproduced_openeval_ifeval_ordinary_probit_H3_visualization" \
  OUT_DIR="$ROOT/results/reproduced_ifeval_ordinary_vs_mixture_factor_visualization" \
  python3 plot_ifeval_ordinary_probit_3d.py
else
  echo "Skipping optional Python/Plotly 3D refresh; existing 3D outputs are preserved."
fi

find "$ROOT" -type f | sort > "$ROOT/MANIFEST.txt"
