#' Plot node stability per community (internal helper)
#'
#' @description
#' Internal helper used by \code{plot.mixmashnet()}
#' to visualize per-node stability by community using a horizontal barplot.
#'
#' @param stab_obj Stability object as returned by \code{membershipStab()},
#'   containing node-level stability and empirical community assignments.
#' @param title Plot title (character). Default:
#'   \code{"Node Stability by Community"}.
#' @param cutoff Numeric from 0 to 1. Optional vertical reference line for a
#'   stability threshold. Use \code{NULL} to hide the line. Default: 0.7.
#'
#' @return A \code{ggplot2} object showing node stability per community.
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom ggplot2 ggplot aes geom_col geom_text geom_vline scale_x_continuous
#' @importFrom ggplot2 scale_color_manual scale_fill_manual theme_minimal labs theme
#' @importFrom ggplot2 element_text element_line element_blank
membershipStab_plot <- function(stab_obj, title = "Node Stability by Community",
                                cutoff = 0.7) {
  stability  <- stab_obj$membership.stability$empirical.dimensions
  membership <- stab_obj$membership$empirical

  palette <- attr(stab_obj, "palette")
  if (is.null(palette)) palette <- stab_obj$community_palette
  if (is.null(palette)) stop("Missing palette: add 'community_palette' to the object.")

  df <- data.frame(
    node      = names(stability),
    stability = as.numeric(stability),
    community = as.factor(membership),
    stringsAsFactors = FALSE
  )

  # order by community, then decreasing stability
  df <- df[order(df$community, -df$stability), ]
  df$node <- factor(df$node, levels = df$node)

  # --- Palette mapping ---
  lv <- levels(df$community)
  pal_vals <- if (!is.null(names(palette))) unname(palette[as.character(lv)]) else palette[seq_along(lv)]
  pal_vals[is.na(pal_vals)] <- "grey80"
  pal_vals <- ifelse(is.na(pal_vals), "grey80", pal_vals)

  # Validate cutoff
  if (!is.null(cutoff)) {
    if (!is.numeric(cutoff) || length(cutoff) != 1L || is.na(cutoff) || cutoff < 0 || cutoff > 1) {
      stop("'cutoff' must be a single numeric value in [0, 1], or NULL.")
    }
  }

  # Auto x-max so labels are not clipped
  x_max <- max(1.1, max(df$stability, na.rm = TRUE) + 0.15)
  if (!is.null(cutoff)) x_max <- max(x_max, cutoff + 0.15)

  p <- ggplot2::ggplot(df, ggplot2::aes(y = node, x = stability, fill = community, color = community)) +
    ggplot2::geom_col(width = 0.7, fill = "white", linewidth = 1.2) +
    ggplot2::geom_text(
      ggplot2::aes(label = round(stability, 2)),
      hjust = -0.2, color = "black", size = 3.2, fontface = "bold"
    ) +
    ggplot2::scale_x_continuous(limits = c(0, x_max), expand = c(0, 0)) +
    ggplot2::scale_fill_manual(values = pal_vals, breaks = lv, labels = lv, name = "Community") +
    ggplot2::scale_color_manual(values = pal_vals, breaks = lv, labels = lv, guide = "none") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::labs(title = title, x = "Node Stability", y = NULL) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 10),
      axis.text.x = ggplot2::element_text(size = 10),
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(color = "grey90", linewidth = 0.4)
    )

  # Optional cutoff line
  if (!is.null(cutoff)) {
    p <- p + ggplot2::geom_vline(xintercept = cutoff, linetype = "dashed", linewidth = 0.6)
  }

  p
}
