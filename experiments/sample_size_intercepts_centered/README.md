# Sample-Size Simulation: Probit IMFM vs Joint-Mixture Gibbs

This experiment compares:

- Proposed method: binary probit independent-mixture factor model, initialized by centered augmented-`Z` SVD and refined by MAP updates.
- Comparator: full Gibbs sampler for a joint-mixture factor model with `G^H` latent profiles and diagonal within-profile factor covariance.

Each repetition simulates a new dataset, including new latent factor draws, class labels, probit augmentation noise, binary responses, and item intercepts/loadings generated from the design.

Post-processing uses the same alignment idea for both methods: choose the best factor permutation and sign pattern to match the data-generating loadings, then compare flattened parameter vectors and parameter blocks.

The current long run is launched by:

```sh
Rscript scripts/sample_size/run_moderate_crossloading_sample_size_MAP_intercepts.R
```

Checkpoint plots can be regenerated with:

```sh
Rscript scripts/sample_size/plot_moderate_crossloading_checkpoint.R
```
