# Sample-Size Simulation Scripts

This folder contains the runnable scripts for the fixed-DGP IFEval-like
simulation study.

## Canonical Workflow

1. Generate the true loading matrices and heatmaps:

   ```bash
   Rscript scripts/sample_size/plot_fixed_ifeval_lambda_heatmaps.R
   ```

2. Launch or resume the simulation:

   ```bash
   Rscript scripts/sample_size/run_fixed_ifeval_lambda_simulation.R
   ```

3. Plot interim or final results from the chunk checkpoints:

   ```bash
   Rscript scripts/sample_size/plot_fixed_ifeval_lambda_progress.R
   ```

## Script Roles

- `run_fixed_ifeval_lambda_simulation.R`: high-level launcher for the fixed
  IFEval-like design. It sets the simulation grid and worker allocation.
- `compare_original_simulation_joint_mfa_gibbs.R`: low-level fitting engine
  called by the launcher for each chunk.
- `plot_fixed_ifeval_lambda_heatmaps.R`: deterministic DGP visualization and
  Lambda table export.
- `plot_fixed_ifeval_lambda_progress.R`: reads checkpoint files from a running
  simulation and creates line/boxplot summaries.

