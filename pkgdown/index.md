# MixMashNet

**MixMashNet** is an R package that provides a unified framework for estimating
and analysing single-layer and multilayer networks using **Mixed Graphical
Models (MGMs)**.

The package is designed for the analysis of complex biomedical and
epidemiological data, allowing the joint modelling of conditional dependencies
between variables of mixed types (continuous, binary, and categorical), both
within and across multiple layers.

---

## Main features

MixMashNet provides tools for:

- estimation of homogeneous and heterogeneous MGM networks;
- specification and estimation of multilayer network structures;
- bootstrap-based confidence intervals for edges and node-level indices;
- computation of centrality and bridge metrics;
- assessment of membership stability and community scores;
- visualisation of networks and associated metrics, including interactive
  exploration via **Shiny** applications.

---

## Installation

The development version of **MixMashNet** can be installed from GitHub:

```r
# install.packages("devtools")
devtools::install_github("ARCbiostat/MixMashNet")
```
