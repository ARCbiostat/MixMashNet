# Compute community scores from a fitted MixMashNet model

Computes subject-level community scores. Community scores are obtained
as weighted sums of the variables belonging to each detected community,
where weights correspond to the standardized community loadings
estimated via
[`EGAnet::net.loads`](https://r-ega.net/reference/net.loads.html) and
stored in the fitted `mixMN_fit` object. Scores are computed using the
dataset provided via the `data` argument. If `data = NULL`, the original
dataset used to fit the model (`fit$model$data`) is used by default.
Optionally, percentile bootstrap quantile regions for the community
scores can be computed if bootstrap community loadings are available in
`fit$community_loadings$boot`. Community scores are only available if
community loadings were computed in the fitted model. This requires that
all variables in the community subgraph are of MGM type Gaussian
(`"g"`), Poisson (`"p"`), or binary categorical (`"c"` with
`level == 2`).

## Usage

``` r
community_scores(
  fit,
  data = NULL,
  layer = NULL,
  scale = TRUE,
  quantile_level = NULL,
  return_quantile_region = FALSE,
  na_action = c("stop", "omit")
)
```

## Arguments

- fit:

  A fitted object of class
  `c("mixmashnet","mixMN_fit", "multimixMN_fit")` returned by
  [`mixMN()`](https://arcbiostat.github.io/MixMashNet/reference/mixMN.md)
  or
  [`multimixMN()`](https://arcbiostat.github.io/MixMashNet/reference/multimixMN.md).

- data:

  Optional data.frame with variables in columns. If `NULL`, uses
  `fit$model$data`. Errors if both are `NULL`.

- layer:

  Optional. If fit is a multimixMN_fit, specify which layer to score
  (name or index). If NULL, scores are computed for all layers and
  returned as a named list.

- scale:

  Logical; if `TRUE` (default), z-standardize variables used for
  scoring, using the mean/SD computed from the dataset used for scoring.

- quantile_level:

  Optional numeric from 0 to 1, e.g. 0.95 or 0.99. If provided,
  percentile bootstrap quantile regions are computed for community
  scores (requires `fit$community_loadings$boot`).

- return_quantile_region:

  Logical; if `TRUE`, return quantile regions.

- na_action:

  Character. How to handle missing values in the scoring data: `"stop"`
  (default) stops if any missing value is present in the required
  variables; `"omit"` computes scores using row-wise omission within
  each community (i.e., uses available variables only, re-normalizing
  weights within community for that row).

## Value

A list with class `c("mixmashnet","community_scores")` containing:

- `call`:

  The matched call.

- `settings`:

  List with `scale`, `quantile_level`, and `na_action`.

- `ids`:

  Character vector of subject IDs (rownames of `data`).

- `communities`:

  Character vector of community score names.

- `scores`:

  Numeric matrix of scores (n × K).

- `quantile_region`:

  If requested and available, a list with `lower` and `upper` matrices
  (n × K) for percentile bootstrap quantile regions; otherwise `NULL`.

- `details`:

  List containing `nodes_used`, `loadings_true`,
  `loadings_boot_available`, and scaling parameters (`center`, `scale`).

If `fit` is a `mixMN_fit` (or a `multimixMN_fit` with `layer`
specified), returns a `c("mixmashnet","community_scores")` object. If
`fit` is a `multimixMN_fit` and `layer = NULL`, returns a named list of
`community_scores` objects (one per layer).

## Details

The function requires that `fit$community_loadings$true` exists and that
the input `data` contains all required variables in
`fit$community_loadings$nodes`. It errors otherwise.

## References

Christensen, A. P., Golino, H., Abad, F. J., & Garrido, L. E. (2025).
Revised network loadings. *Behavior Research Methods*, 57(4), 114.
[doi:10.3758/s13428-025-02640-3](https://doi.org/10.3758/s13428-025-02640-3)
