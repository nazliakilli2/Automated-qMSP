#' Plot amplification curves
#'
#' @param x A `qmsp_run` or `qmsp_result`.
#' @param run Run name(s) to show. Defaults to the first run.
#' @param targets Genes to show (default: all).
#' @param samples Samples to show (default: all).
#' @param log Use a log10 delta Rn axis.
#' @param threshold Draw the instrument threshold for each gene.
#' @return A ggplot object, one panel per gene.
#' @export
plot_amplification <- function(x, run = NULL, targets = NULL, samples = NULL,
                               log = FALSE, threshold = TRUE) {
  run <- run %||% unique(x$wells$run)[1]
  w <- x$wells[x$wells$run %in% run, ]
  if (!is.null(targets)) w <- w[w$target %in% targets, ]
  if (!is.null(samples)) w <- w[w$sample %in% samples, ]
  if (!nrow(w)) stop("Nothing to plot for this selection.", call. = FALSE)

  cv <- x$curves[x$curves$run %in% run, ]
  cv <- merge(cv, w[, intersect(c("run", "well_index", "target", "well",
                                  "sample", "result"), names(w))],
              by = c("run", "well_index", "target"))
  cv$label <- paste0(cv$sample, " (", cv$well, ")")
  if (log) cv$delta_rn <- pmax(cv$delta_rn, 1)

  p <- ggplot2::ggplot(cv, ggplot2::aes(.data$cycle, .data$delta_rn,
                                        group = .data$label,
                                        colour = .data$sample)) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::facet_wrap(~target, scales = "free_y") +
    ggplot2::labs(x = "Cycle", y = expression(Delta * Rn), colour = "Sample",
                  title = paste(run, collapse = ", ")) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
  if ("result" %in% names(cv) && any(cv$result == "Review")) {
    p <- p + ggplot2::aes(linetype = .data$result == "Review") +
      ggplot2::scale_linetype_manual(values = c(`FALSE` = "solid",
                                                `TRUE` = "dashed"),
                                     name = "Needs review")
  }
  th <- x$thresholds
  if (threshold && !is.null(th) && nrow(th)) {
    th <- th[th$run %in% run & th$target %in% unique(cv$target), ]
    if (nrow(th)) {
      p <- p + ggplot2::geom_hline(data = th,
                                   ggplot2::aes(yintercept = .data$threshold),
                                   linetype = "dotted", colour = "grey30")
    }
  }
  if (log) p <- p + ggplot2::scale_y_log10()
  p
}

#' Plot the plate layout
#'
#' @param x A `qmsp_run` or `qmsp_result`.
#' @param run Run to show. Defaults to the first run.
#' @param fill What to colour wells by: `"ct"`, `"target"` or `"result"`
#'   (`"result"` needs a `qmsp_result`).
#' @return A ggplot object.
#' @export
plot_plate <- function(x, run = NULL, fill = c("ct", "target", "result")) {
  fill <- match.arg(fill)
  run <- run %||% unique(x$wells$run)[1]
  w <- x$wells[x$wells$run == run, ]
  if (fill == "result" && !"result" %in% names(w)) {
    stop("fill = \"result\" needs the output of analyze_qmsp().", call. = FALSE)
  }
  n_rows <- max(8L, match(max(w$row), LETTERS))
  n_cols <- max(12L, max(w$col))
  w$label <- paste0(substr(w$sample, 1, 10), "\n", substr(w$target, 1, 10),
                    "\n", ifelse(is.na(w$ct), "Undet.", sprintf("%.1f", w$ct)))
  w$row <- factor(w$row, levels = rev(LETTERS[seq_len(n_rows)]))
  w$col <- factor(w$col, levels = seq_len(n_cols))

  p <- ggplot2::ggplot(w, ggplot2::aes(.data$col, .data$row)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data[[fill]]), colour = "white",
                       linewidth = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = .data$label), size = 2.2,
                       lineheight = 0.9) +
    ggplot2::scale_x_discrete(position = "top", drop = FALSE) +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = NULL, y = NULL, title = run) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  if (fill == "ct") {
    p <- p + ggplot2::scale_fill_gradient(low = "#d94801", high = "#fee6ce",
                                          na.value = "grey92", name = "Ct")
  } else if (fill == "result") {
    p <- p + ggplot2::scale_fill_manual(
      values = c(Positive = "#fdae6b", Negative = "grey90", Review = "#fdd835"),
      name = "Result")
  }
  p
}

#' Methylation heatmap: samples x genes
#'
#' @param x A `qmsp_result` from [analyze_qmsp()].
#' @param value `"call"` (categorical) or `"beta"`, `"pmr"`, `"ratio"`, `"delta_ct"`,
#'   `"ct"` (numeric).
#' @param controls Include control samples?
#' @return A ggplot object.
#' @export
plot_methylation <- function(x, value = c("call", "beta", "pmr", "ratio",
                                          "delta_ct", "ct"),
                             controls = FALSE) {
  value <- match.arg(value)
  if (!inherits(x, "qmsp_result")) {
    stop("`x` must come from analyze_qmsp().", call. = FALSE)
  }
  r <- x$results
  if (!controls) r <- r[r$role == "sample", ]
  if (!nrow(r)) stop("No samples to plot.", call. = FALSE)
  multi_run <- length(unique(r$run)) > 1
  r$sample_label <- if (multi_run) paste0(r$sample, "  [", r$run, "]") else
    r$sample
  r$sample_label <- factor(r$sample_label,
                           levels = rev(sort(unique(r$sample_label))))
  r$fill <- if (value == "call") r$call else r[[value]]
  r$text <- if (value == "call") {
    ifelse(is.na(r$pmr) | r$call != "Methylated", "",
           sprintf("%.0f", r$pmr))
  } else {
    ifelse(is.na(r[[value]]), "", formatC(r[[value]], digits = 3,
                                          format = "g"))
  }

  p <- ggplot2::ggplot(r, ggplot2::aes(.data$target, .data$sample_label)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$fill), colour = "white",
                       linewidth = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = .data$text), size = 3) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 0))
  if (value == "call") {
    p + ggplot2::scale_fill_manual(
      values = c(Methylated = "#c0392b", Unmethylated = "#d6e4f0",
                 Review = "#f5c542", Invalid = "grey70"),
      drop = FALSE, name = "Call (PMR shown)")
  } else {
    p + ggplot2::scale_fill_gradient(low = "#fff5f0", high = "#a50f15",
                                     na.value = "grey92", name = value)
  }
}
