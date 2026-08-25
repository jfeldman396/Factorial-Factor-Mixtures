#!/usr/bin/env zsh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$SCRIPT_DIR"

MATRIX="${MATRIX_PATH:-$ROOT/data/ifeval/openeval_ifeval_only_binary_matrix.csv}"
ITEMS="${ITEM_METADATA_PATH:-$ROOT/data/ifeval/openeval_item_metadata.csv}"
OUT_ROOT="${OUT_ROOT:-$ROOT/results/full/ifeval_reproduction}"
CV_DIR="${CV_DIR:-$OUT_ROOT/componentwise_cv}"
SELECTED_DIR="${SELECTED_DIR:-$OUT_ROOT/selected_componentwise_H4_G3313_lambda4}"
VIZ_DIR="${VIZ_DIR:-$OUT_ROOT/selected_visualizations}"

# To rebuild MATRIX and ITEMS from the full OpenEval model-by-item matrix, run
# scripts/data/format_openeval_binary_matrix.py followed by
# scripts/ifeval/create_ifeval_analysis_matrix.R.  See REPRODUCE.md.

mkdir -p "$OUT_ROOT" "$CV_DIR" "$SELECTED_DIR" "$VIZ_DIR"

echo "Step 1: held-out predictive likelihood CV for mixture models"
MATRIX_PATH="$MATRIX" \
ITEM_METADATA_PATH="$ITEMS" \
OUT_DIR="$CV_DIR" \
H_GRID="${H_GRID:-1:5}" \
G_MODE="${G_MODE:-column_grid}" \
G_GRID="${G_GRID:-1,2,3}" \
G_COMPONENT_VALUES="${G_COMPONENT_VALUES:-1,2,3}" \
MAX_GAUSSIAN_COORDS="${MAX_GAUSSIAN_COORDS:-1}" \
LAMBDA_L1_GRID="${LAMBDA_L1_GRID:-0,1,2,4,8,12}" \
K_FOLDS="${K_FOLDS:-3}" \
WORKERS="${WORKERS:-6}" \
PRETRAIN_AUG_ITER="${PRETRAIN_AUG_ITER:-20}" \
PRETRAIN_Z_UPDATE="${PRETRAIN_Z_UPDATE:-sample}" \
REFINE_ITER="${REFINE_ITER:-20}" \
MIXTURE_MAX_ITER="${MIXTURE_MAX_ITER:-200}" \
FIT_SELECTED_AFTER_CV="${FIT_SELECTED_AFTER_CV:-FALSE}" \
SAVE_FITS="${SAVE_FITS:-FALSE}" \
RESUME_EXISTING="${RESUME_EXISTING:-TRUE}" \
REFRESH_PLOTS="${REFRESH_PLOTS:-TRUE}" \
Rscript cv_ifeval_rank_lambda_models.R

if [[ "${USE_CV_SELECTED_FOR_VISUALS:-FALSE}" == "TRUE" ]]; then
  echo "Step 2: refit the CV-selected full-data mixture model"
  selected_info=("${(@f)$(Rscript -e 'args <- commandArgs(TRUE); s <- read.csv(args[1], stringsAsFactors = FALSE); g_label <- chartr(",", "-", s$G_config[1]); l_label <- chartr(".", "p", as.character(s$lambda_l1_penalty[1])); cat(s$H[1], "\n", s$G_config[1], "\n", s$lambda_l1_penalty[1], "\n", file.path(args[2], sprintf("selected_mixture_H%d_Gconfig%s_lambda%s", s$H[1], g_label, l_label)), "\n", sep = "")' "$CV_DIR/ifeval_rank_lambda_selected_by_heldout_ll.csv" "$CV_DIR")}")
  SELECTED_H_FROM_CV="${selected_info[1]}"
  SELECTED_G_FROM_CV="${selected_info[2]}"
  SELECTED_LAMBDA_FROM_CV="${selected_info[3]}"
  SELECTED_DIR="${selected_info[4]}"
  rds_files=("$SELECTED_DIR"/*_fit.rds(N))
  if [[ ${#rds_files[@]} -eq 0 ]]; then
    MATRIX_PATH="$MATRIX" \
    ITEM_METADATA_PATH="$ITEMS" \
    OUT_DIR="$SELECTED_DIR" \
    H_FIXED="$SELECTED_H_FROM_CV" \
    G_FIXED="$SELECTED_G_FROM_CV" \
    WORKERS="${SELECTED_WORKERS:-8}" \
    PRETRAIN_AUG_ITER="${SELECTED_PRETRAIN_AUG_ITER:-20}" \
    PRETRAIN_Z_UPDATE="${SELECTED_PRETRAIN_Z_UPDATE:-sample}" \
    REFINE_ITER="${SELECTED_REFINE_ITER:-20}" \
    MIXTURE_MAX_ITER="${SELECTED_MIXTURE_MAX_ITER:-200}" \
    REQUIRE_MIXTURE_CONVERGENCE="${SELECTED_REQUIRE_MIXTURE_CONVERGENCE:-TRUE}" \
    REFINEMENT_LAMBDA_L1_PENALTY="$SELECTED_LAMBDA_FROM_CV" \
    Rscript fit_interpret_ifeval_mixture.R
  else
    echo "CV-selected full-data fit already exists in: $SELECTED_DIR"
  fi
else
  echo "Step 2: direct refit of the interpretable component-wise model"
  MATRIX_PATH="$MATRIX" \
  ITEM_METADATA_PATH="$ITEMS" \
  OUT_DIR="$SELECTED_DIR" \
  H_FIXED="${SELECTED_H_FIXED:-4}" \
  G_FIXED="${SELECTED_G_FIXED:-3,3,1,3}" \
  WORKERS="${SELECTED_WORKERS:-8}" \
  PRETRAIN_AUG_ITER="${SELECTED_PRETRAIN_AUG_ITER:-20}" \
  PRETRAIN_Z_UPDATE="${SELECTED_PRETRAIN_Z_UPDATE:-sample}" \
  REFINE_ITER="${SELECTED_REFINE_ITER:-20}" \
  MIXTURE_MAX_ITER="${SELECTED_MIXTURE_MAX_ITER:-200}" \
  REQUIRE_MIXTURE_CONVERGENCE="${SELECTED_REQUIRE_MIXTURE_CONVERGENCE:-TRUE}" \
  REFINEMENT_LAMBDA_L1_PENALTY="${SELECTED_REFINEMENT_LAMBDA_L1_PENALTY:-4}" \
  Rscript fit_interpret_ifeval_mixture.R
fi

echo "Step 3: plot selected-model Lambda heatmaps"
FIT_DIR="$SELECTED_DIR" \
OUT_DIR="$VIZ_DIR" \
Rscript plot_openeval_mixture_lambda_heatmaps.R

echo "Step 4: plot selected-model marginal mixtures"
FIT_DIR="$SELECTED_DIR" \
OUT_DIR="$VIZ_DIR" \
Rscript plot_factor_marginal_mixtures.R

echo "Step 5: summarize selected-model cross-loading items"
LOADINGS_PATH="$SELECTED_DIR/openeval_item_intercepts_loadings_metadata.csv" \
OUT_DIR="$VIZ_DIR/loadings_crossloadings" \
Rscript summarize_ifeval_loadings_crossloadings.R

echo "Step 6: refresh selected-model factor visualizations"
FIT_DIR="$SELECTED_DIR" \
OUT_DIR="$VIZ_DIR/factor_visualizations" \
Rscript plot_ifeval_learned_factors_baseR.R

if python3 -c "import matplotlib, plotly" >/dev/null 2>&1; then
  echo "Step 7: refresh optional Python/Plotly selected-model visualizations"
  FIT_DIR="$SELECTED_DIR" \
  OUT_DIR="$VIZ_DIR/factor_visualizations" \
  python3 plot_ifeval_learned_factors.py
else
  echo "Skipping optional Python/Plotly 3D refresh; matplotlib/plotly unavailable."
fi

find "$OUT_ROOT" -type f | sort > "$OUT_ROOT/MANIFEST.txt"
echo "IFEval reproduction outputs saved in: $OUT_ROOT"
