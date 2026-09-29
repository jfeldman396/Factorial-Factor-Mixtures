# IFEval Results

The primary paper analysis uses the strict threshold-1 binary matrix and
five-fold row-wise few-shot cross-validation. The release under
`threshold_1/rowwise_fewshot/` contains:

- the complete 462-candidate P-IFA ranking and fold scores;
- exact row-fold and reveal/prediction indices;
- the fixed-orientation `H=4`, `G=(2,2,3,1)`, `lambda=4` interpretation tables;
- the matched rank-4 Gaussian Gibbs comparison.

Large fitted R objects, checkpoints, and optimization diagnostics remain in
the ignored `results/full/` tree. Regenerate them with the commands in
`REPRODUCE.md`. `artifact_manifest.csv` records the released files, sizes, and
MD5 checksums.
