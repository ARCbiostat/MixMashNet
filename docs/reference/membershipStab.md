# Node stability from bootstrap community assignments

Computes per-node stability given the empirical community structure and
the homogenized bootstrap memberships contained in a `mixMN_fit` object.
This function is used internally by [`mixMN()`](mixMN.md) and
[`multimixMN()`](multimixMN.md). Stability is expressed as the
proportion of bootstrap replications that assign each node to its
empirical (original) community.

## Usage

``` r
membershipStab(fit, IS.plot = FALSE)
```

## Arguments

- fit:

  An object returned by [`mixMN()`](mixMN.md) (class `mixMN_fit`),
  containing `$communities$original_membership` and
  `$communities$boot_memberships`. Bootstrap memberships must be
  available, i.e. `reps > 0` and `"community" %in% boot_what`.

- IS.plot:

  Logical; if `TRUE`, prints a stability plot via the internal helper
  `membershipStab_plot()`.

## Value

An object of class `c("membershipStab")`, with components:

- `membership`:

  List with:

  `empirical`

  :   Named integer vector of empirical community labels

  `bootstrap`

  :   Matrix of homogenized bootstrap labels (`reps × p`)

- `membership.stability`:

  List with:

  `empirical.dimensions`

  :   Named numeric vector of node-level stability (proportion assigned
      to empirical community)

  `all.dimensions`

  :   Matrix (`p × K`) with proportions of assignment to each community

- `community_palette`:

  Named vector of colors for communities, if available

## Details

Bootstrap community labels are first aligned to the empirical solution
using
[`EGAnet::community.homogenize()`](https://r-ega.net/reference/community.homogenize.html).
Stability is then computed node-wise as the proportion of bootstrap runs
in which the node's community matches its empirical assignment.

## References

Christensen, A. P., & Golino, H. (2021). Estimating the Stability of
Psychological Dimensions via Bootstrap Exploratory Graph Analysis: A
Monte Carlo Simulation and Tutorial. *Psych*, 3(3), 479–500.
[doi:10.3390/psych3030032](https://doi.org/10.3390/psych3030032)
