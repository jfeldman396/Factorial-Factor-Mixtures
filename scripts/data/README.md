# Data Formatting Scripts

## OpenEval To Binary Matrix

`format_openeval_binary_matrix.py` downloads or reads a local snapshot of
`human-centered-eval/OpenEval` and writes model-by-item binary matrices plus
item metadata.

Install the Python data dependencies:

```sh
python3 -m pip install pandas pyarrow huggingface_hub
```

Build the uncapped IFEval source matrix from Hugging Face:

```sh
python3 scripts/data/format_openeval_binary_matrix.py \
  --benchmarks ifeval \
  --binary-threshold 0.5 \
  --out-dir data/openeval_ifeval_formatted_uncapped \
  --min-item-response-prop 0.25 \
  --min-model-response-prop 0.25
```

If the dataset snapshot is already cached locally, add:

```sh
--local-snapshot-dir "$HOME/.cache/huggingface/hub/datasets--human-centered-eval--OpenEval/snapshots/<snapshot-id>"
```

The key outputs are:

- `openeval_binary_matrix_raw.csv`: raw thresholded model-by-item matrix before the final IFEval analysis filter.
- `openeval_response_long_scores.csv`: one row per retained OpenEval response score before binary thresholding.
- `openeval_item_metadata.csv`: item metadata with raw payload, parsed prompt, instruction ids, instruction families, and instruction kwargs.
- `openeval_item_instruction_metadata_long.csv`: one row per item-instruction pair.

For IFEval, the formatter extracts OpenEval `ifeval_strict_accuracy` values,
which are fractions of strict instruction checks satisfied.  The current
snapshot has values `0`, `1/3`, `1/2`, `2/3`, and `1`.
