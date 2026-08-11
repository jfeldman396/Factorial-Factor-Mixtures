# IFEval Analysis

This folder documents the cleaned IFEval workflow.

The analysis starts from the binary item matrix and metadata in `data/ifeval`, then runs:

1. Missing-aware rank selection by held-out probit likelihood.
2. Predictive comparison between the mixture binary probit factor model and ordinary binary probit factor models.
3. Sparse loading penalty tuning for the selected mixture model.
4. Selected `H = 3, G = 3` model interpretation.
5. Loading heatmaps, cross-loading summaries, and exact item examples.
6. 3D visualization of model factor scores and mixture-profile structure.

The rendered writeup is available at `writeup/ifeval_analysis_writeup.pdf`.
