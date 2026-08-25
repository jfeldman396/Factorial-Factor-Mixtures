#!/usr/bin/env zsh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

PYTHON_BIN="${PYTHON_BIN:-python3}"
LOCAL_SNAPSHOT_DIR="${LOCAL_SNAPSHOT_DIR:-$HOME/.cache/huggingface/hub/datasets--human-centered-eval--OpenEval/snapshots/91a9a3a78257a4b4c04e45f5a493796f8a2966b1}"
THRESHOLDS=(${=THRESHOLDS:-0.5 0.75 1})
MIN_ITEM_RESPONSE_PROP="${MIN_ITEM_RESPONSE_PROP:-0.25}"
MIN_MODEL_RESPONSE_PROP="${MIN_MODEL_RESPONSE_PROP:-0.25}"

cd "$ROOT"

for threshold in "${THRESHOLDS[@]}"; do
  label="${threshold//./p}"
  source_dir="$ROOT/data/openeval_ifeval_formatted_threshold_${label}"
  analysis_dir="$ROOT/data/ifeval_threshold_${label}"

  echo "Building OpenEval IFEval matrix for strict-accuracy threshold >= ${threshold}"
  "$PYTHON_BIN" scripts/data/format_openeval_binary_matrix.py \
    --local-snapshot-dir "$LOCAL_SNAPSHOT_DIR" \
    --benchmarks ifeval \
    --binary-threshold "$threshold" \
    --out-dir "$source_dir" \
    --min-item-response-prop "$MIN_ITEM_RESPONSE_PROP" \
    --min-model-response-prop "$MIN_MODEL_RESPONSE_PROP"

  OPENEVAL_FULL_MATRIX="$source_dir/openeval_binary_matrix_raw.csv" \
  OPENEVAL_ITEM_METADATA="$source_dir/openeval_item_metadata.csv" \
  OUT_DIR="$analysis_dir" \
  MIN_ITEM_RESPONSE_PROP="$MIN_ITEM_RESPONSE_PROP" \
  MIN_MODEL_RESPONSE_PROP="$MIN_MODEL_RESPONSE_PROP" \
  Rscript scripts/ifeval/create_ifeval_analysis_matrix.R
done

echo "Threshold-specific IFEval matrices written under data/ifeval_threshold_*"
