#!/bin/zsh

cd "/Users/joefeldman/Documents/Deep Factor Models" || exit 1

OUT_DIR="${OUT_DIR:-results/moderate_crossloading_joint_mfa_sample_size_p500_25reps_H3H4_G2G3_MAP}"

while true; do
  script_dir="$(cd "$(dirname "$0")" && pwd)"
  OUT_DIR="$OUT_DIR" /usr/local/bin/Rscript "$script_dir/plot_moderate_crossloading_checkpoint.R"
  if [ -f "$OUT_DIR/comparison_results.csv" ]; then
    exit 0
  fi
  sleep 600
done
