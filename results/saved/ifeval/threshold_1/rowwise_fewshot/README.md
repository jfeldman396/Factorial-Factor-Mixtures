# Threshold-1 Row-Wise Few-Shot Results

These are the analysis-ready outputs for the paper's primary IFEval study.
The binary matrix contains 122 LLMs and 534 complete, non-degenerate IFEval
items. Five row folds hold out entire LLMs. Within each held-out row, 20% of
observed responses estimate the factor score and the remaining 80% evaluate
prediction, with three fixed reveal repetitions.

## Model Selection

`model_selection/` contains the complete 462-candidate P-IFA ranking, all
2,310 candidate-fold summaries, and the exact row folds and response indices.
The predictive winner is `H=4`, `G=(1,3,3,3)`, `lambda=4`. The near-tied
`G=(1,2,2,3)` component-count multiset is used for the paper's parsimonious
interpretation.

## Interpretation Model

`interpretation_model/` contains the full-data fixed-orientation
`H=4`, `G=(2,2,3,1)`, `lambda=4` summaries. The factor order is held fixed so
the component counts remain attached to the displayed factors. The item table
joins canonical loadings to prompts and instruction metadata.

## Gaussian Comparison

`gaussian_comparison/` contains the matched rank-4 Gaussian probit Gibbs
comparison using the exact P-IFA folds and reveal masks. The P-IFA winner has
mean log score `-0.30584`; the parsimonious interpretation configuration has
`-0.30617`; and the Gaussian model has `-0.30879`. Accuracy and Brier scores
are nearly identical, while P-IFA has the stronger predictive log score.

Large fitted objects and checkpoints are intentionally omitted. Use
`REPRODUCE.md` to regenerate them.
