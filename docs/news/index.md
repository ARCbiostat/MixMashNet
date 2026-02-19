# Changelog

## MixMashNet 0.4.0

### Major updates

- Added [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) for
  [`community_scores()`](https://arcbiostat.github.io/MixMashNet/reference/community_scores.md).

## MixMashNet 0.3.0

### Major updates

- Simplified the interface of
  [`mixMN()`](https://arcbiostat.github.io/MixMashNet/reference/mixMN.md)
  and
  [`multimixMN()`](https://arcbiostat.github.io/MixMashNet/reference/multimixMN.md)
  by removing the `type` and `level` arguments. Variable types are now
  inferred automatically.
- Community score computation is now restricted to gaussian, poisson and
  binary categorical variables.

## MixMashNet 0.1.0

- Development version.
