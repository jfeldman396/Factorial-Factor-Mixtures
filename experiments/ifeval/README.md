# IFEval Analysis

This folder documents the cleaned IFEval workflow.

The analysis starts from the binary item matrix and metadata in `data/ifeval`, then runs:

1. Singular-value shelf diagnostics.
2. Missing-aware held-out predictive likelihood tuning over rank, component-wise `G`, and sparse loading penalty.
3. Selected component-wise mixture refit with MAP refinement.
4. Loading heatmaps, cross-loading summaries, and exact item examples.
5. Factor-score and mixture-profile visualization.

The current component-wise interpretation uses `H = 4`,
`G = (3,3,1,3)`, and sparse-loading MAP refinement with
`lambda_l1_penalty = 4`. The rendered writeups are available at:

- `writeup/ifeval_componentwise_G3313.pdf`;
- `writeup/ifeval_analysis_writeup.pdf`.

See `REPLICATION.md` for the exact commands and `CODE_AUDIT.md` for the
current reproducibility audit notes.
