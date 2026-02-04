# Print method for MixMashNet objects

Compact textual summary for objects returned by [`mixMN()`](mixMN.md)
and [`multimixMN()`](multimixMN.md). The method reports:

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

- main settings for community detection and bootstrap.

## Usage

``` r
# S3 method for class 'mixmashnet'
print(x, ...)
```

## Arguments

- x:

  An object of class `mixmashnet`, typically returned by
  [`mixMN()`](mixMN.md) or [`multimixMN()`](multimixMN.md).

- ...:

  Additional arguments.

## Value

The input object `x`, returned invisibly.

## See also

[`mixMN`](mixMN.md), [`multimixMN`](multimixMN.md)
