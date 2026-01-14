#' Global variables used in non-standard evaluation
#'
#' These names are referenced in dplyr/ggplot2 code paths and are declared
#' to silence R CMD check notes about undefined globals.
#'
#' @name MixMashNet-globals
#' @keywords internal
#' @noRd
NULL

utils::globalVariables(c(
  "observed","abs_obs","lower","upper","node","community","community_factor",
  "node_order_value","node_order_alpha","node_order_comm",
  "order_reversed","label_colored","includes_zero","item",
  "sum_abs_w","sum_signed_w","sum_signed_w2",
  "inv_mean_dist","pair","hits"
))
