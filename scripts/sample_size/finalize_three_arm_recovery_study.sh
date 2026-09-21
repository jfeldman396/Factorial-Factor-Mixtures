#!/bin/zsh

# Assemble validated results and export the complete paper-ready figure set.

set -eu

export PATH="/usr/local/bin:/opt/homebrew/bin:/usr/bin:/bin:/usr/sbin:/sbin"

SCRIPT_DIR="${0:A:h}"
REPO_ROOT="${SCRIPT_DIR:h:h}"
RSCRIPT="${RSCRIPT:-$(command -v Rscript)}"
RUN_LABEL_PREFIX="${RUN_LABEL_PREFIX:-fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_npenalty3_5_full}"

cd "${REPO_ROOT}"
export RUN_LABEL_PREFIX
"${RSCRIPT}" scripts/sample_size/assemble_three_arm_recovery_results.R

plot_arm() {
  local arm="$1"
  local subtitle="$2"
  export RUN_LABEL="${RUN_LABEL_PREFIX}_${arm}"
  export RESULTS_FILE="${REPO_ROOT}/results/full/${RUN_LABEL}/comparison_results.csv"
  export PLOT_DIR="${REPO_ROOT}/results/selected_plots/sample_size/three_mixture_settings/${arm}"
  export TABLE_DIR="${REPO_ROOT}/results/selected_tables/sample_size/three_mixture_settings"
  export METHOD_FILTER="independent_marginal_mixture,viroli_laplace_gibbs"
  export G_COMPONENT_FILTER="2,3"
  export P_FILTER="500,1000,1500,2000"
  export OUTPUT_TAG="${arm}_product_map_vs_gibbs"
  export PLOT_SUBTITLE="${subtitle}"
  export PLOT_WIDTH=15.8
  export PLOT_HEIGHT=9.8
  export PLOT_DPI=300
  export WRITE_PDF=TRUE
  "${RSCRIPT}" scripts/sample_size/plot_fixed_ifeval_grouped_boxplot_panels.R
}

plot_arm separated \
  "Separated symmetric mixtures; Gibbs was run for p=500/1000"
plot_arm asymmetric_pi \
  "Asymmetric mixture probabilities; Gibbs has 25 reps at n=100 and 5 reps at n=200/400"
plot_arm moderate_overlap \
  "Moderately overlapping mixtures; Gibbs has 25 reps at n=100 and 5 reps at n=200/400"

"${RSCRIPT}" -e '
root <- "results/selected_plots/sample_size/three_mixture_settings"
arms <- c("separated", "asymmetric_pi", "moderate_overlap")
expected <- c("factor_score_rmse", "probability_rmse", "lambda_rmse", "alpha_rmse",
              "marginal_mu_rmse", "marginal_var_rmse", "marginal_weight_rmse", "seconds")
for (arm in arms) {
  files <- list.files(file.path(root, arm), full.names = TRUE)
  for (metric in expected) {
    for (ext in c("png", "pdf")) {
      hit <- files[grepl(paste0("_", metric, "\\\\.", ext, "$"), files)]
      if (length(hit) != 1L || file.info(hit)$size <= 0) {
        stop("Missing or empty ", arm, " ", metric, " ", ext, " export.")
      }
    }
  }
}
cat("Verified 48 non-empty paper-ready figure files.\n")
'

echo "Three-arm recovery finalization complete."
