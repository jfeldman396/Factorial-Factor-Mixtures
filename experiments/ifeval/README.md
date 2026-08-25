# IFEval Analysis

This folder documents the cleaned IFEval workflow.

The analysis starts from the binary item matrix and metadata in `data/ifeval`.
The binary entries are thresholded OpenEval `ifeval_strict_accuracy` scores:
the data formatter extracts numeric values from each nested OpenEval score
object, averages them within model-item pairs, and codes an entry as `1` when
that average is at least `0.5`.

The model-facing IFEval matrix is built by starting from the full OpenEval
model-by-item matrix, selecting columns with the `ifeval_` prefix, applying
coverage filters, requiring complete retained item columns, and retaining
nonconstant items.  The currently committed uncapped matrix starts from `541`
IFEval items and contains `122` models by `534` retained items.  See
`data/ifeval/README.md` and
`scripts/ifeval/create_ifeval_analysis_matrix.R` for the exact build step.
The retained item metadata includes raw prompts, evaluator instruction ids,
instruction families, instruction kwargs, and a long table with one row per
item-instruction pair.

The fitted workflow then runs:

1. Missing-aware held-out predictive likelihood tuning over rank, component-wise `G`, and sparse loading penalty.
2. Selected component-wise mixture refit with MAP refinement.
3. Loading heatmaps, cross-loading summaries, and exact item examples.
4. Factor-score and mixture-profile visualization.

The current component-wise interpretation uses `H = 4`,
`G = (3,3,1,3)`, and sparse-loading MAP refinement with
`lambda_l1_penalty = 4`. The rendered writeups are available at:

- `writeup/ifeval_componentwise_G3313.pdf`;
- `writeup/ifeval_analysis_writeup.pdf`.

See `REPLICATION.md` for the exact commands and `CODE_AUDIT.md` for the
current reproducibility audit notes.
