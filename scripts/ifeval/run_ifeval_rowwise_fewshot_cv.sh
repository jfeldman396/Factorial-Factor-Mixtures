#!/usr/bin/env zsh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

MATRIX_PATH="${MATRIX_PATH:-$ROOT/data/ifeval_threshold_1/openeval_ifeval_only_binary_matrix.csv}" \
OUT_DIR="${OUT_DIR:-$ROOT/results/full/ifeval_rowwise_fewshot_cv}" \
H_GRID="${H_GRID:-2:8}" \
G_VALUES="${G_VALUES:-1,2,3}" \
LAMBDA_L1_GRID="${LAMBDA_L1_GRID:-0,1,2,4,8,12}" \
K_FOLDS="${K_FOLDS:-5}" \
REVEAL_FRACTION="${REVEAL_FRACTION:-0.20}" \
N_REVEAL_REPEATS="${N_REVEAL_REPEATS:-3}" \
SEED="${SEED:-1}" \
WORKERS="${WORKERS:-18}" \
PARALLEL="${PARALLEL:-TRUE}" \
RESUME_EXISTING="${RESUME_EXISTING:-TRUE}" \
TASK_BATCH_SIZE="${TASK_BATCH_SIZE:-${WORKERS:-18}}" \
FIT_SELECTED_FULL_DATA="${FIT_SELECTED_FULL_DATA:-TRUE}" \
Rscript "$SCRIPT_DIR/run_ifeval_rowwise_fewshot_cv.R"
