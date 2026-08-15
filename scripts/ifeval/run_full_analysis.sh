#!/usr/bin/env zsh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$SCRIPT_DIR"

MATRIX="$ROOT/data/ifeval/openeval_ifeval_only_binary_matrix.csv"
ITEMS="$ROOT/data/ifeval/openeval_item_metadata.csv"
OUT_ROOT="$ROOT/results/full/ifeval_reproduction"
CV_DIR="$OUT_ROOT/componentwise_cv"
SELECTED_DIR="$OUT_ROOT/selected_componentwise_H4_G3313_lambda4"
VIZ_DIR="$OUT_ROOT/selected_visualizations"

mkdir -p "$OUT_ROOT" "$CV_DIR" "$SELECTED_DIR" "$VIZ_DIR"

echo "Step 1: singular-value shelf diagnostic"
MATRIX_PATH="$MATRIX" \
OUT_DIR="$OUT_ROOT/rank_diagnostics" \
H_MAX=10 \
Rscript plot_ifeval_rank_diagnostics.R

echo "Step 2: held-out predictive likelihood CV for mixture models"
MATRIX_PATH="$MATRIX" \
ITEM_METADATA_PATH="$ITEMS" \
OUT_DIR="$CV_DIR" \
H_GRID=1:5 \
G_MODE=column_grid \
G_GRID=1,2,3 \
G_COMPONENT_VALUES=1,2,3 \
MAX_GAUSSIAN_COORDS=1 \
LAMBDA_L1_GRID=0,1,2,4,8,12 \
K_FOLDS=3 \
WORKERS=6 \
PRETRAIN_AUG_ITER=20 \
REFINE_ITER=20 \
MIXTURE_MAX_ITER=200 \
FIT_SELECTED_AFTER_CV=TRUE \
SAVE_FITS=FALSE \
RESUME_EXISTING=TRUE \
REFRESH_PLOTS=TRUE \
Rscript cv_ifeval_rank_lambda_models.R

echo "Step 3: direct refit of the interpretable component-wise model"
MATRIX_PATH="$MATRIX" \
ITEM_METADATA_PATH="$ITEMS" \
OUT_DIR="$SELECTED_DIR" \
H_FIXED=4 \
G_FIXED=3,3,1,3 \
WORKERS=8 \
PRETRAIN_AUG_ITER=20 \
REFINE_ITER=20 \
MIXTURE_MAX_ITER=200 \
REQUIRE_MIXTURE_CONVERGENCE=TRUE \
REFINEMENT_LAMBDA_L1_PENALTY=4 \
Rscript fit_interpret_ifeval_mixture.R

echo "Step 4: plot selected-model Lambda heatmaps"
FIT_DIR="$SELECTED_DIR" \
OUT_DIR="$VIZ_DIR" \
Rscript plot_openeval_mixture_lambda_heatmaps.R

echo "Step 5: plot selected-model marginal mixtures"
FIT_DIR="$SELECTED_DIR" \
OUT_DIR="$VIZ_DIR" \
Rscript plot_factor_marginal_mixtures.R

echo "Step 6: summarize selected-model cross-loading items"
LOADINGS_PATH="$SELECTED_DIR/openeval_item_intercepts_loadings_metadata.csv" \
OUT_DIR="$VIZ_DIR/loadings_crossloadings" \
Rscript summarize_ifeval_loadings_crossloadings.R

echo "Step 7: refresh selected-model factor visualizations"
FIT_DIR="$SELECTED_DIR" \
OUT_DIR="$VIZ_DIR/factor_visualizations" \
Rscript plot_ifeval_learned_factors_baseR.R

if python3 -c "import matplotlib, plotly" >/dev/null 2>&1; then
  echo "Step 8: refresh optional Python/Plotly selected-model visualizations"
  FIT_DIR="$SELECTED_DIR" \
  OUT_DIR="$VIZ_DIR/factor_visualizations" \
  python3 plot_ifeval_learned_factors.py
else
  echo "Skipping optional Python/Plotly 3D refresh; matplotlib/plotly unavailable."
fi

find "$OUT_ROOT" -type f | sort > "$OUT_ROOT/MANIFEST.txt"
echo "IFEval reproduction outputs saved in: $OUT_ROOT"
