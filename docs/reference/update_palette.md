# Update community and layer color palettes in MixMashNet objects

Updates the color palettes associated with communities and/or layers in
`mixMN_fit` and `multimixMN_fit` objects. The function replaces only the
colors corresponding to the provided names, leaving all other colors
unchanged.

If colors are provided for names that do not exist in the object (e.g.,
unknown community labels or layer names), a warning is issued and those
entries are ignored. If some communities or layers are not specified,
their original colors are preserved.

## Usage

``` r
update_palette(fit, community_colors = NULL, layer_colors = NULL)
```

## Arguments

- fit:

  An object of class `mixMN_fit` or `multimixMN_fit`.

- community_colors:

  Optional named character vector specifying new colors for communities.
  Names must correspond to existing community labels (as stored in
  `communities$palette`). Missing names are ignored.

- layer_colors:

  Optional named character vector specifying new colors for layers.
  Names must correspond to existing layer names (as stored in
  `layers$palette`). Only applicable to `multimixMN_fit` objects.

## Value

The input object `fit`, with updated community and/or layer palettes.

## Details

For `mixMN_fit` objects, community colors are updated in
`fit$communities$palette`.

For `multimixMN_fit` objects, community colors are updated separately
within each layer (i.e., in `fit$layer_fits[[L]]$communities$palette`),
while layer colors are updated in `fit$layers$palette`.

The function performs in-place modification of the palettes and returns
the updated object.
