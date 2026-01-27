#' Global variables used in non-standard evaluation
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
  "inv_mean_dist","pair","hits",
  "n_g","m","s"
))
