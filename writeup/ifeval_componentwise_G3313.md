# IFEval Component-Wise Mixture Fit: H=4, G=(3,3,1,3)

This note summarizes a focused IFEval fit of the binary probit independent-mixture factor model with

```text
H = 4
G = (3, 3, 1, 3)
lambda_l1_penalty = 4
pretraining iterations = 20 maximum
MAP refinement iterations = 20 maximum
```

The fitted full-data model stopped after 15 pretraining iterations and 12 MAP refinement iterations. The full-data binary probit log likelihood was `-0.21897` per response. In the completed cross-validation table, this specification has held-out log likelihood `-0.26978` per response, very close to the current best component-wise fit `G=(3,3,2,2)`, which has held-out log likelihood `-0.26944` per response.

## Component-Wise CV Context

The current component-wise CV screen strongly favors `H=4` and `lambda=4`. Among these fits, `G=(3,3,1,3)` is near the top while keeping the third coordinate Gaussian.

![Component-wise CV metrics for H=4, lambda=4](ifeval_componentwise_G3313_assets/componentwise_cv_H4_lambda4.png)

This is useful substantively: F3 is retained as a continuous axis, while F1, F2, and F4 are allowed to form discrete model subtypes.

## Marginal Factor Distributions

![Marginal mixture fits](ifeval_componentwise_G3313_assets/marginal_mixtures_G3313.png)

The fitted marginal mixture structure is:

| Factor | Components | Interpretation |
| --- | ---: | --- |
| F1 | 3 | Broad general instruction-following axis with high/mid/low model groups. |
| F2 | 3 | A discrete axis for language/case/script constraints, especially all-lowercase and single-language requirements. |
| F3 | 1 | A continuous axis rather than a discrete subtype. This is consistent with the earlier finding that extra F3 components mostly model central concentration rather than a meaningful left-tail group. |
| F4 | 3 | A discrete formatting/structured-output axis, including JSON wrapping, sentence-count limits, and local formatting constraints. |

The fitted component sizes were:

| Factor | Component counts |
| --- | --- |
| F1 | 27 / 16 / 79 |
| F2 | 21 / 61 / 40 |
| F3 | 122 |
| F4 | 20 / 50 / 52 |

## Loading Matrix

![Lambda ordered by strongest factor](ifeval_componentwise_G3313_assets/lambda_ordered_by_strongest_factor_G3313.png)

The heatmap is a row permutation of the fitted `500 x 4` loading matrix. For each item `j`, define

```text
strongest_factor_j = argmax_h |lambda_jh|
max_abs_loading_j = max_h |lambda_jh|.
```

Rows are sorted by `strongest_factor_j`, then by decreasing `max_abs_loading_j`. Thus, the apparent blocks are not imposed by the model; they are a diagnostic ordering of the learned loadings. The fitted structure is dominated by a broad F1 factor, with smaller private-factor blocks for F2, F3, and F4.

Among active items with `|loading| > 0.1`, the strongest-factor counts are:

| Strongest factor | Active items |
| --- | ---: |
| F1 | 410 |
| F2 | 37 |
| F3 | 26 |
| F4 | 23 |

The numeric loading summary also shows this structure:

| Factor | Mean abs loading | `|lambda| > 0.5` | `|lambda| > 1.0` |
| --- | ---: | ---: | ---: |
| F1 | 0.863 | 396 | 199 |
| F2 | 0.240 | 38 | 28 |
| F3 | 0.203 | 63 | 13 |
| F4 | 0.222 | 69 | 18 |

## Factor Interpretation

### F1: Broad Instruction Following

F1 is the dominant general axis. High-loading items often combine ordinary task completion with explicit constraints, such as repeating the request exactly, ending with a required phrase, using all caps, or wrapping text in prescribed delimiters.

Representative high-F1 items include:

| Item | Loading | Example constraint |
| --- | ---: | --- |
| `ifeval_20260421T021146Z_127` | 1.785 | End the response with a specified sentence and no words after it. |
| `ifeval_20260421T021146Z_242` | 1.758 | Repeat the request exactly, then answer with a title in double angle brackets. |
| `ifeval_20260421T021146Z_362` | 1.742 | Repeat the request word-for-word before answering. |
| `ifeval_20260421T021146Z_6` | 1.701 | Write two labeled all-caps paragraphs. |

### F2: Language, Case, and Script Control

F2 is more specialized. The strongest positive F2 items emphasize lowercasing, single-language output, and constraints on allowed alphabets or scripts.

Representative high-F2 items include:

| Item | Loading | Example constraint |
| --- | ---: | --- |
| `ifeval_20260421T021146Z_293` | 3.071 | Product description entirely in lowercase. |
| `ifeval_20260421T021146Z_409` | 2.967 | English-only response, all lowercase. |
| `ifeval_20260421T021146Z_280` | 2.963 | Hindi poem, only Hindi, no other language. |
| `ifeval_20260421T021146Z_292` | 2.822 | Haiku in lowercase with no capital letters. |

### F3: Creative Lexical Constraint Axis

F3 is best treated as continuous in this fit. Its high-loading items often require creative generation under severe lexical constraints, such as avoiding or limiting a common letter while still producing a poem, song, quiz, or long essay.

Representative high-F3 items include:

| Item | Loading | Example constraint |
| --- | ---: | --- |
| `ifeval_20260421T021146Z_277` | 1.637 | Write a rubric as a poem while using the letter `w` fewer than two times. |
| `ifeval_20260421T021146Z_477` | 1.565 | Write a song while using the letter `a` at most once. |
| `ifeval_20260421T021146Z_23` | 1.516 | Write a logic quiz while using the letter `t` at most once. |
| `ifeval_20260421T021146Z_147` | 1.324 | Write a long essay of at least 50 sentences. |

The Gaussian choice for F3 is important: it suggests this ability varies continuously across LLMs rather than forming a stable discrete subtype in the current data.

### F4: Structured Output and Local Formatting

F4 captures a smaller but interpretable formatting axis. High-loading items frequently require JSON wrapping, sentence-count limits, required endings, highlighted sections, or other local structural constraints.

Representative high-F4 items include:

| Item | Loading | Example constraint |
| --- | ---: | --- |
| `ifeval_20260421T021146Z_269` | 1.219 | Wrap the entire response in JSON format. |
| `ifeval_20260421T021146Z_58` | 1.213 | Answer in JSON format, with markdown ticks allowed. |
| `ifeval_20260421T021146Z_359` | 1.207 | Use fewer than 10 sentences and end with a required sentence. |
| `ifeval_20260421T021146Z_33` | 1.173 | Use a repeated keyword, fewer than six sentences, and highlighted text sections. |

## Working Interpretation

The `G=(3,3,1,3)` fit gives a clean compromise between predictive performance and interpretability. F1 is the broad general IFEval ability. F2 and F4 describe discrete model subtypes for specialized constraint-following modes. F3 remains a continuous creative-control axis, especially for generative tasks with severe lexical restrictions.

This structure is more interpretable than forcing F3 to have multiple mixture components: when F3 is given two or three components, the extra components tend to model central concentration rather than a meaningful left-tail or high-skill subtype.
