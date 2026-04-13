# Extract a single layer from a multilayer MixMashNet object

Extracts one layer from a fitted multilayer `multimixMN_fit` object
returned by
[`multimixMN()`](https://arcbiostat.github.io/MixMashNet/reference/multimixMN.md).

The selected layer is returned as the corresponding single layer
`mixMN_fit` object stored in `layer_fits`.

## Usage

``` r
layer_slice(object, ...)

# S3 method for class 'multimixMN_fit'
layer_slice(object, layer, ...)
```

## Arguments

- object:

  An object of class `"multimixMN_fit"` returned by
  [`multimixMN()`](https://arcbiostat.github.io/MixMashNet/reference/multimixMN.md).

- ...:

  Further arguments passed to methods.

- layer:

  Character string giving the layer to extract.

## Value

An object of class `"mixMN_fit"` corresponding to the selected layer.
