#!/usr/bin/env zsh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Threshold label 0p67 uses the exact two-thirds cutoff.  A literal 0.67
# threshold would exclude scores equal to 2/3 and would coincide with 1.0.
THRESHOLD_LABELS=(${=THRESHOLD_LABELS:-0p5 0p67 1})
WORKERS="${WORKERS:-18}"
H_GRID="${H_GRID:-1:6}"
LAMBDA_L1_GRID="${LAMBDA_L1_GRID:-0,1,2,4,8,12}"
OUT_BASE="${OUT_BASE:-$ROOT/results/full/ifeval_threshold_sensitivity}"

for label in "${THRESHOLD_LABELS[@]}"; do
  matrix="$ROOT/data/ifeval_threshold_${label}/openeval_ifeval_only_binary_matrix.csv"
  items="$ROOT/data/ifeval_threshold_${label}/openeval_item_metadata.csv"
  out_root="$OUT_BASE/threshold_${label}"

  if [[ ! -f "$matrix" ]]; then
    echo "Missing matrix for threshold ${label}: $matrix" >&2
    exit 1
  fi

  echo "Running IFEval analysis for threshold label ${label}"
  MATRIX_PATH="$matrix" \
  ITEM_METADATA_PATH="$items" \
  OUT_ROOT="$out_root" \
  H_GRID="$H_GRID" \
  G_MODE=column_grid \
  G_GRID=1,2,3 \
  G_COMPONENT_VALUES=1,2,3 \
  MAX_GAUSSIAN_COORDS=1 \
  LAMBDA_L1_GRID="$LAMBDA_L1_GRID" \
  WORKERS="$WORKERS" \
  PRETRAIN_Z_UPDATE=sample \
  PRETRAIN_AUG_ITER="${PRETRAIN_AUG_ITER:-20}" \
  REFINE_ITER="${REFINE_ITER:-20}" \
  MIXTURE_MAX_ITER="${MIXTURE_MAX_ITER:-200}" \
  USE_CV_SELECTED_FOR_VISUALS=TRUE \
  RESUME_EXISTING=TRUE \
  zsh "$SCRIPT_DIR/run_full_analysis.sh"
done

echo "Threshold analyses saved under: $OUT_BASE"
