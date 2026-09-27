# Simulation and Rotation Scripts

## Recovery Study

- `run_three_arm_recovery_study.sh`: authoritative three-setting launcher.
- `run_fixed_ifeval_lambda_simulation.R`: resumable chunk scheduler.
- `compare_original_simulation_joint_mfa_gibbs.R`: simulation, fitting, alignment, and metric engine.
- `run_component_profile_recovery_addendum.sh`: matched rerun that records observation-level subtype recovery.
- `build_simulation_release_results.R`: combines the matched rerun with the large-`p` Product MAP extension.
- `validate_three_arm_recovery_study.R`: coverage, duplicate-key, and matched-seed validation.
- `plot_paper_dgp_figures.R`: deterministic DGP loading and mixture figures.
- `make_simulation_section_artifacts.R`: all paper and supplement recovery figures and LaTeX tables.

## Rotation Ablation

- `run_rotation_ablation_three_arms.sh`: authoritative three-setting ablation launcher.
- `compare_rotation_fastica_diagnostic.R`: estimated/oracle signal fitting engine.
- `summarize_rotation_ablation_three_arms.R`: validation, summaries, and paper figures.
- `plot_rotation_ablation_estimated_z_boxplots.R`: paneled estimated-signal error figure.

See the repository-level `REPRODUCE.md` for the fixed grid, worker settings, commands, and output locations.
