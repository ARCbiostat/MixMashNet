#' MixMashNet: Multilayer and Single Layer Network Modeling
#'
#' @description
#' Tools for estimating and analyzing single layer and multilayer networks
#' using Mixed Graphical Models (MGMs), accommodating continuous, count, and
#' categorical variables.
#' In the multilayer setting, layers may comprise different types and numbers
#' of variables, and users can explicitly impose a predefined multilayer
#' topology to constrain the estimation of inter and intralayer connections.
#' The package implements bootstrap procedures to derive confidence intervals
#' for edge weights and node-level centrality and bridge metrics, and provides
#' tools to assess the stability of node community membership. In addition,
#' subject-level community scores can be computed to summarize the latent
#' dimensions identified through network clustering, facilitating downstream
#' statistical analyses.
#'
#' @references
#'
#' Christensen, A. P., & Golino, H. (2021).
#' Estimating the Stability of Psychological Dimensions via Bootstrap Exploratory Graph Analysis:
#' A Monte Carlo Simulation and Tutorial. \emph{Psych}, 3(3), 479–500.
#' \doi{10.3390/psych3030032}
#'
#' Christensen, A. P., Golino, H., Abad, F. J., & Garrido, L. E. (2025).
#' Revised network loadings. \emph{Behavior Research Methods}, 57(4), 114.
#' \doi{10.3758/s13428-025-02640-3}
#'
#' Epskamp, S., Borsboom, D., & Fried, E. I. (2018).
#' Estimating psychological networks and their accuracy: A tutorial paper.
#' \emph{Behavior Research Methods}, 50(1), 195–212.
#' \doi{10.3758/s13428-017-0862-1}
#'
#' Haslbeck, J. M. B., & Waldorp, L. J. (2020).
#' mgm: Estimating Time-Varying Mixed Graphical Models in High-Dimensional Data.
#' \emph{Journal of Statistical Software}, 93(8).
#' \doi{10.18637/jss.v093.i08}
#'
#' Jones, P. J., Ma, R., & McNally, R. J. (2021).
#' Bridge Centrality: A Network Approach to Understanding Comorbidity.
#' \emph{Multivariate Behavioral Research}, 56(2), 353–367.
#' \doi{10.1080/00273171.2019.1614898}
#'
#' @name MixMashNet-package
#' @docType package
#' @author Maria De Martino, Caterina Gregorio, Adrien Perigord
#'
#' \email{maria.demartino@uniud.it}
#' @keywords package
"_PACKAGE"
