# IFEval Analysis

This folder documents the cleaned IFEval workflow.

The analysis starts from the binary item matrix and metadata in `data/ifeval`, then runs:

1. Missing-aware rank selection by held-out probit likelihood.
2. Predictive comparison between the mixture binary probit factor model and ordinary binary probit factor models.
3. Sparse loading penalty tuning for the selected mixture model.
4. Selected component-wise mixture interpretation.
5. Loading heatmaps, cross-loading summaries, and exact item examples.
6. 3D visualization of model factor scores and mixture-profile structure.

The current component-wise interpretation uses `H = 4`,
`G = (3,3,1,3)`, and sparse-loading MAP refinement with
`lambda_l1_penalty = 4`. The rendered writeups are available at:

- `writeup/ifeval_componentwise_G3313.pdf`;
- `writeup/ifeval_analysis_writeup.pdf`.

See `REPLICATION.md` for the exact commands and `CODE_AUDIT.md` for the
current reproducibility audit notes.
