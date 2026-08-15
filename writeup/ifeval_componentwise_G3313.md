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

## Fitted Mixture Structure

![Marginal mixture fits](ifeval_componentwise_G3313_assets/marginal_mixtures_G3313.png)

The selected fit uses component counts `G=(3,3,1,3)`. Thus F1, F2, and F4 have discrete low/middle/high model groups, while F3 is deliberately modeled as one continuous Gaussian coordinate. The group labels below are interpreted using the reported factor scores and MAP group assignments in `openeval_model_factor_scores_profiles.csv`.

| Factor | Groups | Score interpretation | Substantive interpretation |
| --- | --- | --- | --- |
| F1 | 1/2/3: 27 / 16 / 79 models | low / middle / high | Broad instruction-following ability. |
| F2 | 1/2/3: 21 / 61 / 40 models | low / middle / high | Language, case, and script-control ability. |
| F3 | 1: 122 models | continuous | Creative lexical-control ability, not a discrete subtype. |
| F4 | 1/2/3: 20 / 50 / 52 models | low / middle / high | Structured-output and local-formatting ability. |

Empirically, the MAP groups have the following reported factor-score summaries. The weights in this table are empirical occupancies, so they are directly tied to the plotted model profiles.

| Factor | Group | Weight | Mean score | SD score |
| --- | --- | ---: | ---: | ---: |
| F1 | low | 0.221 | -1.085 | 0.446 |
| F1 | middle | 0.131 | 0.336 | 0.193 |
| F1 | high | 0.648 | 1.221 | 0.287 |
| F2 | low | 0.172 | -2.004 | 1.164 |
| F2 | middle | 0.500 | -0.545 | 0.346 |
| F2 | high | 0.328 | 0.593 | 0.237 |
| F3 | continuous | 1.000 | -0.138 | 0.960 |
| F4 | low | 0.164 | -2.169 | 0.770 |
| F4 | middle | 0.410 | -0.235 | 0.284 |
| F4 | high | 0.426 | 0.516 | 0.231 |

The fitted mixture is therefore not just a ranking model. It separates broad IFEval strength from more specific forms of constraint-following. For example, a model can be in the high-F1 group while only middle on F2 or F4, meaning it generally follows instructions well but is less reliably specialized for language/case/script constraints or structured-output constraints.

## LLM Scores and Mixture Profiles

The next two plots show the fitted LLM coordinates. Models are ordered from left to right by empirical IFEval accuracy.

![Refined factor scores by LLM](ifeval_componentwise_G3313_assets/llm_factor_scores_G3313_cropped.png)

The continuous score heatmap shows the broad role of F1: the highest-accuracy models tend to have high F1 scores, while weaker models tend to move downward on F1. F2, F3, and F4 are less monotone in overall accuracy. They capture different styles of constraint-following rather than just "better versus worse" performance.

The group heatmap translates the continuous scores into factor-wise MAP mixture assignments. F3 is intentionally a single Gaussian coordinate, so every model has group 1 on F3. F1, F2, and F4 show discrete heterogeneity: models can be strong on the broad instruction-following axis while belonging to different F2 or F4 subtypes.

![MAP mixture groups by LLM](ifeval_componentwise_G3313_assets/llm_mixture_groups_G3313_cropped.png)

### LLM Ability Profiles

The marginal factor coordinates give a more nuanced model profile than overall accuracy alone. The joint profiles are interpreted in detail after the loading matrix because the meaning of each coordinate is clearer once the item loadings are inspected.

- **F1: broad instruction-following ability.** High F1 models include `claude-sonnet-4-5-20250929`, `llama-4-scout-17b-16e-instruct`, `gpt-5.4-pro-2026-03-05`, and `gpt-5.2-pro-2025-12-11`. Low F1 models include small or older chat/instruct models such as `qwen-1.5-0.5b-chat`, `falcon-7b-instruct`, and `redpajama-incite` variants. This is the closest coordinate to a broad IFEval ability axis.

- **F2: language, case, and script-control ability.** High F2 models include `qwen-2.5-32b-instruct`, `granite-4.0-micro`, `gpt-4.1-nano-2025-04-14`, `gpt-5.1-2025-11-13`, and `gpt-4o-2024-08-06`. Low F2 models include `gemma-3-1b-it`, `gemma-3-4b-it`, and `gemma-2-2b-it`. This axis separates models on lowercasing, single-language output, and script-control constraints.

- **F3: continuous creative lexical-control ability.** High F3 models include `gpt-5-2025-08-07`, `gpt-5.4-pro-2026-03-05`, `gpt-5-nano-2025-08-07`, `gpt-5-mini-2025-08-07`, and `o3-2025-04-16`. Low F3 models include several Qwen and Mistral-family entries such as `qwen-1.5-32b-chat`, `qwen-2.5-0.5b-instruct`, and `mistral-7b-instruct-v0.3`. Because F3 is modeled with one Gaussian component, these are high/low positions on a continuous axis rather than discrete subtypes.

- **F4: structured-output and local-formatting ability.** High F4 models include `nova-pro-v1:0`, `claude-3-5-sonnet-20241022`, `gemini-1.5-flash-002`, `qwen-2.5-14b-instruct`, and `gpt-5-nano-2025-08-07`. Low F4 models include `llama-2` base models, `gemma-3-12b-it`, `redpajama-incite-7b-chat`, and `qwen-1.5-32b-chat`. This coordinate appears to capture JSON wrapping, exact endings, sentence-count limits, and similar formatting demands.

## Loading Matrix

![Lambda ordered by strongest factor](ifeval_componentwise_G3313_assets/lambda_ordered_by_strongest_factor_G3313.png)

The heatmap displays the fitted item loading matrix, `Lambda`, after permuting the 500 item rows. Each row is one IFEval item and each column is one factor. Red entries are positive loadings, blue entries are negative loadings, and white entries are near zero. Larger absolute values mean that the item's success probability changes more strongly with that factor score.

The row ordering is diagnostic, not part of the fitted model. For each item `j`, define

```text
strongest_factor_j = argmax_h |lambda_jh|
max_abs_loading_j = max_h |lambda_jh|.
```

Rows are sorted by `strongest_factor_j`, then by decreasing `max_abs_loading_j`. Thus, the apparent blocks are not imposed by the model; they are a visualization of the learned loading pattern. The large upper block consists of items whose strongest absolute loading is on F1. The smaller lower blocks consist of items whose strongest loading is on F2, F3, or F4. The fact that these blocks become visible after a genuine row permutation is evidence that the fitted loadings have a broad-plus-specific structure: one dominant general IFEval factor and several smaller private axes.

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

The loadings give the item-level interpretation of the factors. A strong solo-loading item has one dominant loading and small loadings on the other coordinates, so it is useful for naming a factor. A cross-loading item has large loadings on more than one coordinate, so it reveals tasks that require a combination of abilities. This section lists exact IFEval item ids, fitted loading values, and the item text or constraint being measured.

### Specific Strong Solo-Loading Items

These items have one clearly dominant loading. They are the cleanest examples for interpreting each coordinate.

- **F1 solo item:** `ifeval_20260421T021146Z_242`, with loadings `(1.758, 0.000, 0.000, 0.000)`. The item asks the model to first repeat the request word for word, then name a black dog, and include a title wrapped in double angular brackets. This is a clean broad instruction-following item because success requires remembering and executing several explicit instructions in sequence.

- **F2 solo item:** `ifeval_20260421T021146Z_293`, with loadings `(0.137, 3.071, 0.105, 0.000)`. The item asks: "Create a product description for a product that will help me to stop snoring. Use all lowercase letters." This is a clean language/case-control item.

- **F3 solo item:** `ifeval_20260421T021146Z_277`, with loadings `(0.167, -0.087, 1.637, 0.343)`. The item asks for a rubric written as a poem while using the letter `w` fewer than two times. This is a clean creative lexical-control item because the response must be generative and stylistic while obeying a severe letter constraint.

- **F4 solo item:** `ifeval_20260421T021146Z_359`, with loadings `(0.479, 0.248, 0.000, 1.207)`. The item asks how to repair a damaged ship hull, requires fewer than 10 sentences, and requires the response to end with the exact sentence "That is all you need!" This is a clean local-formatting item because the main difficulty is satisfying the sentence-count and required-ending constraints.

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

### Specific Cross-Loading Items

Cross-loading items reveal tasks that mix multiple skill demands. Using the rule that two factors both have `|loading| >= 0.5`, the observed pair types are F1/F2, F1/F3, F1/F4, and F3/F4. No F2/F3 or F2/F4 items met this threshold in this fit.

- **F1/F2 cross-loading:** `ifeval_20260421T021146Z_245`, with loadings `(1.001, 0.707, 0.000, 0.416)`. The item asks for exactly three names for a black and white dog using markdown bullet points. It combines broad instruction-following with enumeration and list-format control.

- **F1/F3 cross-loading:** `ifeval_20260421T021146Z_80`, with loadings `(1.300, 0.198, 1.301, -0.018)`. The item asks for a 30-line poem with short sentences, no commas, exactly one sentence per line, correct punctuation at the end of each line, and no extra content beyond the poem. It combines broad constraint-following with creative poetic generation under punctuation and line-structure restrictions.

- **F1/F4 cross-loading:** `ifeval_20260421T021146Z_269`, with loadings `(1.199, 0.000, -0.177, 1.219)`. The item asks what can be said about Radcliffe after he gets past guards, but requires the entire response to be wrapped in JSON format. It combines semantic interpretation with structured-output formatting.

- **F3/F4 cross-loading:** `ifeval_20260421T021146Z_477`, with loadings `(0.000, 0.000, 1.565, 0.816)`. The item asks for a song about Layton while using the letter `a` at most once. It combines creative generation with severe local lexical control.

The cross-loading pattern reinforces the heatmap interpretation. The learned matrix is not pure block diagonal. It is closer to a broad general factor with smaller specialized factors that partially overlap when an item combines content, language/case, lexical creativity, and output-format constraints.

## Joint LLM Evaluation by Mixture Profile

The joint MAP profile has the form `F1-F2-F3-F4`. Since F3 has one Gaussian component, the third entry is always `1`; the continuous F3 score still matters, but it does not create a discrete cluster. For F1, F2, and F4, group `3` is the high-score group, group `2` is the middle group, and group `1` is the low-score group. The profiles below jointly evaluate LLMs by their combination of abilities, not by any single coordinate.

- **Profile `3-3-1-3`: strong on all discrete axes.** This profile contains 14 models, with mean accuracy `0.948` and range `0.910-0.992`. Examples include `gpt-5.4-pro-2026-03-05`, with scores `(1.817, 0.828, 1.652, 0.413)` and accuracy `0.992`; `gpt-5.2-pro-2025-12-11`, with `(1.783, 0.610, 0.960, 0.595)` and accuracy `0.984`; and `grok-3-mini-beta`, with `(1.525, 0.603, 1.153, 0.277)` and accuracy `0.978`. These models excel broadly and also retain strength on language/case/script and structured-output constraints. Their remaining weaknesses are not captured by a low discrete cluster; they are more item-specific or tied to the continuous F3 creative lexical axis.

- **Profile `3-3-1-2`: strong broad and language/case/script ability, middle structured-output ability.** This profile contains 14 models, with mean accuracy `0.919` and range `0.794-0.974`. Examples include `grok-4-0709`, with scores `(1.465, 0.625, 0.876, -0.254)` and accuracy `0.974`; `gpt-5.1-2025-11-13`, with `(1.521, 0.856, 0.815, -0.064)` and accuracy `0.972`; and `qwen-3-80b-instruct`, with `(1.446, 0.711, 0.197, -0.443)` and accuracy `0.946`. These models excel on F1 and F2, but their relative weakness is F4: JSON wrapping, exact endings, sentence-count limits, and similar local formatting demands.

- **Profile `3-2-1-3`: strong broad and structured-output ability, middle language/case/script ability.** This profile contains 17 models, with mean accuracy `0.895` and range `0.842-0.930`. Examples include `grok-3-beta`, with scores `(1.586, -0.686, 1.033, 0.335)` and accuracy `0.930`; `gpt-5.4-2026-03-05`, with `(1.478, -0.479, 1.168, 0.340)` and accuracy `0.926`; and `gpt-5-2025-08-07`, with `(1.335, -0.731, 2.050, 0.484)` and accuracy `0.922`. These models excel on broad instruction-following and structured output, but they are less specialized for lowercase, single-language, and script-restricted tasks.

- **Profile `3-1-1-2`: high broad ability, weak language/case/script ability, middle structured-output ability.** This profile contains 10 models, with mean accuracy `0.865` and range `0.802-0.908`. Examples include `claude-sonnet-4-5-20250929`, with scores `(1.832, -1.871, 0.544, -0.212)` and accuracy `0.908`; `claude-sonnet-4-20250514`, with `(1.445, -1.596, 0.292, 0.055)` and accuracy `0.896`; and `palmyra-x5`, with `(1.600, -1.490, -0.196, -0.158)` and accuracy `0.886`. These models are genuinely strong on broad instruction-following, but the profile says where they suffer: language/case/script control is low and structured-output control is only middle.

- **Profile `1-2-1-1`: low broad and structured-output ability.** This profile contains 8 models, with mean accuracy `0.318` and range `0.244-0.392`. Examples include `llama-2-70b-hf`, with scores `(-0.746, -0.187, 0.110, -3.335)` and accuracy `0.392`; `llama-2-7b-hf`, with `(-0.942, -0.344, 0.691, -3.208)` and accuracy `0.376`; and `mistral-7b-v0.1`, with `(-1.190, -0.368, 0.961, -2.228)` and accuracy `0.334`. These models may have moderate F2 or F3 scores, but low F1 and low F4 produce broad instruction-following failures plus weak local formatting control.

- **Profile `1-1-1-1`: jointly weak on the discrete axes.** This profile contains 4 models, with mean accuracy `0.365` and range `0.208-0.452`. Examples include `gemma-2-9b-it`, with scores `(-0.425, -2.172, -0.570, -1.243)` and accuracy `0.452`; `gemma-2-2b-it`, with `(-0.607, -3.114, -0.181, -1.295)` and accuracy `0.412`; and `redpajama-incite-7b-chat`, with `(-1.671, -1.622, 0.275, -2.854)` and accuracy `0.208`. These models suffer jointly: low broad instruction-following, low language/case/script control, and low structured-output control.

This joint-profile view is the main practical advantage of the fitted mixture. Overall accuracy says which models are better on average, while the profile identifies how they are better or worse: broad instruction-following, language/case/script control, creative lexical control, or structured-output formatting.

## Where Weaker Profiles Fail

The fitted factors also identify concrete item types that expose weaker models. The examples below compare the weak or targeted-weakness profile to the strong `3-3-1-3` profile. The accuracies are empirical accuracies within that profile on the named item.

- **Low-F1 broad instruction-following failures:** The weak `1-2-1-1` profile includes older/smaller models such as `llama-2-70b-hf`, `llama-2-7b-hf`, `llama-2-13b-hf`, and `mistral-7b-v0.1`. On `ifeval_20260421T021146Z_242`, with loadings `(1.758, 0.000, 0.000, 0.000)`, the item asks the model to repeat the request word for word, then name a black dog, and include a title in double angular brackets. The `1-2-1-1` profile has accuracy `0.000`, while the strong `3-3-1-3` profile has accuracy `1.000`. The same profile also fails `ifeval_20260421T021146Z_127`, which requires ending with an exact sentence and no words after it; again the profile accuracy is `0.000` versus `1.000` for `3-3-1-3`. These are broad F1 failures: the model loses track of explicit multi-part instructions.

- **Low-F2 language/case/script failures:** The `3-1-1-2` profile contains models that are not weak overall, including `claude-sonnet-4-5-20250929`, `claude-sonnet-4-20250514`, and `palmyra-x5`, but their low F2 score predicts failures on language and case-control items. On `ifeval_20260421T021146Z_293`, with loadings `(0.137, 3.071, 0.105, 0.000)`, the prompt asks for a product description using all lowercase letters; this profile has accuracy `0.000`, compared with `1.000` for `3-3-1-3`. It also has accuracy `0.000` on `ifeval_20260421T021146Z_409`, which requires an English-only answer with only lowercase letters, and on `ifeval_20260421T021146Z_280`, which asks for a Hindi poem with no other language. This is exactly the F2 weakness: strong general ability does not guarantee reliable case, language, or script control.

- **Low-F4 structured-output failures:** The `1-2-1-1` profile is also weak on F4. On `ifeval_20260421T021146Z_269`, with loadings `(1.199, 0.000, -0.177, 1.219)`, the item asks a simple semantic question but requires the entire response to be wrapped in JSON; the weak profile has accuracy `0.000`, while `3-3-1-3` has accuracy `1.000`. The same profile has accuracy `0.000` on `ifeval_20260421T021146Z_359`, which requires fewer than 10 sentences and an exact final sentence. These failures are not about domain knowledge; they are local formatting and structured-output failures.

- **Creative lexical-control failures:** The jointly weak `1-1-1-1` profile includes models such as `gemma-2-9b-it`, `gemma-2-2b-it`, `gemma-3-4b-it`, and `redpajama-incite-7b-chat`. On `ifeval_20260421T021146Z_80`, with loadings `(1.300, 0.198, 1.301, -0.018)`, the prompt asks for a 30-line poem with no commas, exactly one sentence per line, correct punctuation, and no extra content. The `1-1-1-1` profile has accuracy `0.000`, compared with `1.000` for `3-3-1-3`. On `ifeval_20260421T021146Z_277`, which asks for a rubric as a poem while using the letter `w` fewer than two times, the `1-1-1-1` profile again has accuracy `0.000`. These items combine F1 and F3: the model must generate creative text while obeying exact lexical constraints.

These examples show why the factor mixture is useful for evaluation. A weak model is not merely "bad"; it fails in recognizable ways. Some profiles fail broad instruction-following, some fail language/case/script control despite strong general ability, some fail local formatting, and some fail creative lexical control.

## Working Interpretation

The `G=(3,3,1,3)` fit gives a clean compromise between predictive performance and interpretability. F1 is the broad general IFEval ability. F2 and F4 describe discrete model subtypes for specialized constraint-following modes. F3 remains a continuous creative-control axis, especially for generative tasks with severe lexical restrictions.

This structure is more interpretable than forcing F3 to have multiple mixture components: when F3 is given two or three components, the extra components tend to model central concentration rather than a meaningful left-tail or high-skill subtype.
