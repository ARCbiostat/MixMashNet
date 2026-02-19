# Print method for MixMashNet objects

Compact textual summary for objects returned by
[`mixMN()`](https://arcbiostat.github.io/MixMashNet/reference/mixMN.md)
and
[`multimixMN()`](https://arcbiostat.github.io/MixMashNet/reference/multimixMN.md).
The method reports:

- whether the fit is single layer (`mixMN`) or multilayer
  (`multimixMN`);

- number of subjects (if available) and variables;

- for multilayer fits, the number of nodes and non-zero edges per layer
  and, if present, per interlayer pair;

- size of the global graph (nodes and edges);

- number of communities (single layer) or communities per layer
  (multilayer);

- covariates used for adjustment and nodes excluded from the graph
  and/or clustering;

- main settings for community detection and bootstrap;

- data info.

## Usage

``` r
# S3 method for class 'mixmashnet'
print(x, ...)
```

## Arguments

- x:

  An object of class `mixmashnet`, returned by
  [`mixMN()`](https://arcbiostat.github.io/MixMashNet/reference/mixMN.md)
  or
  [`multimixMN()`](https://arcbiostat.github.io/MixMashNet/reference/multimixMN.md).

- ...:

  Additional arguments.

## Value

The input object `x`, returned invisibly.

## See also

[`mixMN`](https://arcbiostat.github.io/MixMashNet/reference/mixMN.md),
[`multimixMN`](https://arcbiostat.github.io/MixMashNet/reference/multimixMN.md)
