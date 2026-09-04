# Factorial Factor Mixtures

This repository collects the reproducible code, selected results, and writeups for factorial factor mixture analyses.

The current contents are organized around two analysis tracks:

1. `scripts/sample_size`: simulations comparing the proposed binary probit independent-mixture factor method against Viroli-style probit Gibbs baselines.
2. `scripts/ifeval`: IFEval rank, component-count, and sparsity selection by held-out predictive likelihood, followed by factor interpretation/visualization for the selected mixture model.

The main implementation lives in `R/`. Runnable scripts are kept under `scripts/` and source the shared implementation so the simulation and IFEval analyses use the same fitting code.

See `REPLICATION.md` for step-by-step instructions to rerun the IFEval
analysis and the sample-size simulation. See `CODE_AUDIT.md` for the latest
static audit notes and reproducibility caveats.

## Current Simulation

The active sample-size experiment uses:

- `n in {100, 200}`
- `p in {500, 1000, 1500, 2000}` for product MAP
- `p in {500, 1000}` for the Viroli Gibbs baselines
- `H in {5, 10, 15, 20}`
- `G in {2, 3}`, applied to every factor coordinate
- `25` Monte Carlo repetitions per setting
- Cross/IFEval-like loadings with randomly signed cross-loadings
- loading magnitudes `Uniform(1.25, 1.75)` or `Uniform(2.50, 3.00)`
- cross-loading probability `0.075` or `0.20`
- balanced and IFEval-like unbalanced item blocks
- item intercepts generated in an IFEval-like pattern
- product MAP with EM-SVD probit signal pretraining, sparse rotation, and MAP refinement
- Viroli-style probit Gibbs with either a Laplace loading prior or a diffuse Gaussian loading prior

The launcher is
`scripts/sample_size/run_final_product_viroli_simulation.R`.  The raw full
outputs are written to `results/full/signal_support_grid` and are ignored by
git.  Selected design heatmaps, progress plots, and summary CSVs are stored in
`results/selected_plots/sample_size` and `results/selected_tables/sample_size`.

## IFEval Analysis

The IFEval analysis includes:

- missing-aware rank/component/penalty selection by held-out predictive likelihood
- tuned sparse loading penalty for the mixture model
- selected-model loading interpretation, including column-specific mixture sizes when selected by CV
- cross-loading summaries
- 3D factor visualizations
- LaTeX writeup and rendered PDF in `writeup/`

### IFEval Data Files

The raw IFEval source is **not** a CSV in this repository.  It comes from the
OpenEval Hugging Face dataset, `human-centered-eval/OpenEval`, where the data
are stored as parquet shards:

```text
item/ifeval-00000-of-00001.parquet
response/ifeval-00000-of-00004.parquet
response/ifeval-00001-of-00004.parquet
response/ifeval-00002-of-00004.parquet
response/ifeval-00003-of-00004.parquet
```

The item parquet contains one row per prompt.  For IFEval, the raw item payload
contains the prompt text, the evaluator instruction ids, and instruction
kwargs.  The response parquet files contain one row per model response, with a
nested model descriptor and nested score object.  For this benchmark, the score
metric is `ifeval_strict_accuracy`, interpreted as the fraction of strict
instruction checks satisfied by that model on that prompt.

The formatter
[`scripts/data/format_openeval_binary_matrix.py`](scripts/data/format_openeval_binary_matrix.py)
converts the raw parquet response shards into CSV files by:

1. reading the requested OpenEval response shards;
2. extracting `model_name`, `item_id`, score metric values, and item metadata;
3. averaging score values within each model-item pair;
4. coding the binary response as `1(score >= threshold)`;
5. writing model-by-item matrices and prompt/instruction metadata.

The compact committed IFEval CSVs are here:

```text
data/ifeval/openeval_ifeval_only_binary_matrix.csv
data/ifeval/openeval_item_metadata.csv
data/ifeval/openeval_item_instruction_metadata_long.csv
data/ifeval/openeval_model_metadata.csv
data/ifeval/ifeval_analysis_matrix_build_summary.csv
```

`openeval_ifeval_only_binary_matrix.csv` is the model-by-item matrix used for
fitting.  Rows are LLMs, columns are retained IFEval prompts, and entries are
binary thresholded strict-accuracy scores.  The default threshold is
`score >= 0.5`.

`openeval_item_metadata.csv` is the wide item metadata table.  It includes the
raw OpenEval item payload plus parsed columns for `prompt`, `instruction_ids`,
`instruction_families`, `n_instructions`, `n_unique_instruction_ids`, and
`instruction_kwargs`.

`openeval_item_instruction_metadata_long.csv` is the long instruction table:
one row per retained prompt-instruction pair.  This is the table to use when
connecting the original IFEval instruction checks to the fitted low-rank
factor model.

The default committed matrix starts from `124` models and `541` IFEval prompts.
After removing two low-coverage models and seven constant prompts, it contains
`122` models and `534` retained prompts.  These `534` prompts contain `820`
prompt-instruction rows, `25` unique instruction ids, and `9` broader
instruction families.

Threshold sensitivity versions of the same CSVs are committed in:

```text
data/ifeval_threshold_0p5/
data/ifeval_threshold_0p67/
data/ifeval_threshold_1/
```

The `0p67` folder uses the exact two-thirds cutoff, $2/3$, but is labeled
`0p67` for readability.  This matters because the observed IFEval
strict-accuracy values are `0`, `1/3`, `1/2`, `2/3`, and `1`; using a literal
decimal threshold of `0.67` would exclude scores equal to `2/3` and would
therefore coincide with the threshold-1.0 matrix.

The retained matrix sizes differ slightly by threshold:

| Folder | Strict-accuracy rule | Models | Retained prompts | Mean binary score |
| --- | --- | ---: | ---: | ---: |
| `data/ifeval_threshold_0p5/` | `score >= 0.5` | 122 | 534 | 0.757 |
| `data/ifeval_threshold_0p67/` | `score >= 2/3` | 122 | 538 | 0.653 |
| `data/ifeval_threshold_1/` | `score >= 1` | 122 | 534 | 0.619 |

The current compact component-wise IFEval writeup is available as a rendered
PDF and as source Markdown:

- [Rendered IFEval component-wise PDF](writeup/ifeval_componentwise_G3313.pdf)
- [IFEval component-wise Markdown source](writeup/ifeval_componentwise_G3313.md)

It summarizes the `H = 4`, `G = (3,3,1,3)` fit, representative solo-loading
and cross-loading items, joint LLM ability profiles, the loading heatmap,
marginal mixture fits, and factor-score visualizations.

The cleaned IFEval data used by the scripts live in `data/ifeval`.
The exact OpenEval-to-IFEval matrix construction is documented in
`data/ifeval/README.md` and `REPRODUCE.md`; entries are thresholded OpenEval
numeric scores, then the IFEval analysis retains complete, nonconstant item
columns.  The same folder includes parsed prompt and instruction metadata for
each retained item, including a long item-instruction table for comparing the
25 retained IFEval instruction ids to the fitted lower-rank factors.

## Remote GitHub Repository

The configured remote is:

```sh
https://github.com/jfeldman396/Factorial-Factor-Mixtures.git
```
