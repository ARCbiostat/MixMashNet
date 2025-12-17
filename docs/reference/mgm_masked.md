# Masked MGM (modified from mgm)

A modified version of the model-fitting routine from mgm that adds
support for a per-node predictor mask via `mask_list`. When provided,
each node is estimated using only the allowed predictors specified for
that node. All other functionality mirrors the original mgm nodewise
estimation and graph assembly.

## Usage

``` r
mgm_masked(
  data,
  type,
  level,
  regularize,
  lambdaSeq,
  lambdaSel,
  lambdaFolds,
  lambdaGam,
  alphaSeq,
  alphaSel,
  alphaFolds,
  alphaGam,
  k,
  moderators,
  ruleReg,
  weights,
  threshold,
  method,
  binarySign,
  scale,
  verbatim,
  pbar,
  warnings,
  saveModels,
  saveData,
  overparameterize,
  thresholdCat,
  signInfo,
  mask_list = NULL,
  ...
)
```

## Arguments

- data:

  Numeric matrix (n × p). No missing values allowed.

- type:

  Character vector of length `p` with variable types as in mgm (e.g.,
  `"g"` for Gaussian, `"c"` for categorical, `"p"` for Poisson).

- level:

  Integer vector of length `p` with variable levels as in mgm.

- regularize:

  Logical; if `FALSE`, equivalent to no regularization
  (`lambdaSel = "EBIC"`, `lambdaSeq = 0`, `threshold = "none"`).

- lambdaSeq, lambdaSel, lambdaFolds, lambdaGam:

  Lambda grid/selection settings (see mgm).

- alphaSeq, alphaSel, alphaFolds, alphaGam:

  Mixing parameter grid/selection (see mgm).

- k:

  Interaction order (default `2`); same meaning as in mgm.

- moderators:

  Optional moderators specification (as in mgm).

- ruleReg:

  Regularization rule, default `"AND"` (as in mgm).

- weights:

  Observation weights (length `n`); internally normalized.

- threshold:

  Thresholding rule (default `"LW"`) as in mgm.

- method:

  Fitting backend (default `"glm"`) as in mgm.

- binarySign:

  Logical; if `TRUE`, store sign information.

- scale:

  Logical; if `TRUE`, standardize Gaussian variables.

- verbatim, pbar, warnings, saveModels, saveData:

  Overhead/UX flags as in mgm.

- overparameterize:

  Logical; use overparameterized design matrix (as in mgm).

- thresholdCat:

  Logical; categorical thresholding (as in mgm).

- signInfo:

  Logical; store sign information (as in mgm).

- mask_list:

  Optional list of length `p`; each element is an integer vector of
  column indices (1..p) indicating which predictors are allowed for the
  nodewise regression of the corresponding target node. If `NULL`, no
  masking is applied.

- ...:

  Further arguments forwarded consistently with mgm.

## Value

An object of class `c("mgm","core")` with fields analogous to mgm:

- `$pairwise$`: `wadj`, `signs`, and nodewise variants.

- `$interactions$`, `$intercepts$`, `$nodemodels$`: as in mgm.

- `$call$`: input metadata (types, levels, thresholds, etc.).

## Details

This function is adapted from the mgm workflow and preserves its
nodewise estimation, cross-validation / EBIC selection, thresholding,
and graph assembly. The only substantive extension is `mask_list`, which
filters each node's design matrix to the allowed predictors before
fitting. This enables constrained structures (e.g., multilayer masks)
while keeping mgm’s estimation logic intact.

## See also

[`mgm`](https://rdrr.io/pkg/mgm/man/mgm.html) for the original
implementation and argument semantics.
