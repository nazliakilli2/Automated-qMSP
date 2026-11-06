#' Find a Ct cutoff for each gene from labelled samples
#'
#' Uses samples whose diagnosis is known to choose, for every gene, the Ct
#' cutoff that best separates cancer samples from healthy controls. A sample
#' counts as methylated for a gene when its Ct is at or below the cutoff, so
#' a good cutoff is below most cancer Cts and above most control Cts.
#'
#' The steps are:
#' \enumerate{
#'   \item the Ct of every well is read where its amplification curve crosses
#'     `threshold` (delta Rn 10,000 by default, see [analyze_qmsp()]);
#'   \item every sample gets one Ct per gene: the mean Ct of its wells when at
#'     least half of them have a Ct, otherwise no Ct. Wells from different
#'     runs with the same sample name are pooled;
#'   \item samples whose reference gene has no Ct at or below `ref_ct_max`
#'     (failed DNA) are left out, unless `exclude_failed_reference = FALSE`;
#'   \item every possible cutoff is tried. For each, sensitivity (cancer
#'     samples at or below it) and specificity (controls above it, or with no
#'     Ct) are computed. The chosen cutoff is halfway between the last Ct it
#'     includes and the next Ct up.
#' }
#'
#' Methods:
#' \describe{
#'   \item{`"youden"`}{best balance: the highest sensitivity + specificity.}
#'   \item{`"specificity"`}{the highest sensitivity among cutoffs with a
#'     specificity of at least `target` (few false "cancer" results).}
#'   \item{`"sensitivity"`}{the highest specificity among cutoffs with a
#'     sensitivity of at least `target` (few missed cancers).}
#' }
#' Ties go to the higher specificity, then to the lower cutoff.
#'
#' With few samples the cutoff can move a lot when one sample is added; the
#' `note` column warns when a group has fewer than 10 samples.
#'
#' @param x A `qmsp_run` from [read_eds()] / [read_eds_files()].
#' @param groups The group of every sample: a data frame with columns
#'   `sample` and `group`, or a named character vector (`c(KAS48 = "cancer")`).
#'   Groups are `cancer`, `control` (healthy, non-cancer), `positive`
#'   (positive control), `ntc` (no-template / NK water control) or `exclude`.
#'   Samples that are not listed are left out, except that names matching
#'   `ntc` or `positive` are taken as those controls.
#' @param method How to choose the cutoff: `"youden"`, `"specificity"` or
#'   `"sensitivity"`.
#' @param target The minimum specificity or sensitivity for those methods.
#' @param threshold Fluorescence (delta Rn) threshold at which the Ct is read.
#' @param reference,ntc,positive Regular expressions for the reference gene
#'   and the control sample names, as in [analyze_qmsp()].
#' @param ref_ct_max Reference gene Ct above this (or none) means failed DNA.
#' @param exclude_failed_reference Leave out samples with failed DNA?
#' @param cutoffs Optional named vector of cutoffs chosen by hand, e.g.
#'   `c(TAC1 = 36.5)`. Those genes keep that cutoff (the statistics are
#'   computed for it); the others get the suggested one.
#'
#' @return A `qmsp_cutoffs` list with
#'   \describe{
#'     \item{cutoffs}{one row per gene: `cutoff` (the one used),
#'       `suggested_cutoff`, `source` ("suggested" or "chosen"), `sensitivity`,
#'       `specificity`, `auc` (area under the ROC curve: 1 = perfect
#'       separation, 0.5 = no better than chance), the counts `tp`, `fn`,
#'       `fp`, `tn`, the positive and NK control Cts, and a `note`.}
#'     \item{samples}{one row per sample and gene with its group and Ct.}
#'     \item{roc}{every cutoff tried, with its sensitivity and specificity.}
#'     \item{settings}{the settings used.}
#'   }
#' @export
find_cutoffs <- function(x, groups,
                         method = c("youden", "specificity", "sensitivity"),
                         target = 0.9,
                         threshold = 10000,
                         reference = "ACTB|B.?ACTIN|BETA.?ACTIN",
                         ntc = "NTC|dH2O|dH20|water|blank|^NK",
                         positive = "H460|A549|HT29|positive|^PC\\b",
                         ref_ct_max = 40,
                         exclude_failed_reference = TRUE,
                         cutoffs = NULL) {
  if (!inherits(x, "qmsp_run")) {
    stop("`x` must come from read_eds() or read_eds_files().", call. = FALSE)
  }
  method <- match.arg(method)
  settings <- list(method = method, target = target, threshold = threshold,
                   reference = reference, ref_ct_max = ref_ct_max,
                   exclude_failed_reference = exclude_failed_reference)

  wells <- apply_thresholds(x$wells, x$curves, threshold)
  wells$group <- sample_groups(wells$sample, groups, ntc, positive)
  wells$is_reference <- if (is.null(reference)) rep(FALSE, nrow(wells)) else
    matches(wells$target, reference)

  samples <- per_sample_ct(wells, ref_ct_max)
  samples$used <- samples$group %in% c("cancer", "control") &
    !(exclude_failed_reference & samples$ref_failed %in% TRUE)

  genes <- unique(samples$target)
  n_cycles <- if (!is.null(x$curves) && nrow(x$curves))
    max(x$curves$cycle, na.rm = TRUE) else 50
  rows <- lapply(genes, function(g) {
    s <- samples[samples$target == g, ]
    fixed <- if (g %in% names(cutoffs)) unname(cutoffs[[g]]) else NULL
    gene_cutoff(g, s, method, target, n_cycles, fixed)
  })
  table <- do.call(rbind, lapply(rows, `[[`, "cutoff"))
  roc <- do.call(rbind, lapply(rows, `[[`, "roc"))
  rownames(table) <- NULL
  if (!is.null(roc)) rownames(roc) <- NULL

  structure(list(cutoffs = table, samples = samples, roc = roc,
                 settings = settings),
            class = "qmsp_cutoffs")
}

#' @export
print.qmsp_cutoffs <- function(x, ...) {
  s <- x$samples[x$samples$used, ]
  cat("<qmsp_cutoffs> ", length(unique(s$sample[s$group == "cancer"])),
      " cancer and ", length(unique(s$sample[s$group == "control"])),
      " control samples, method: ", x$settings$method, "\n", sep = "")
  c <- x$cutoffs
  out <- data.frame(gene = c$target, cutoff = c$cutoff,
                    sensitivity = round(c$sensitivity, 2),
                    specificity = round(c$specificity, 2),
                    auc = round(c$auc, 2), note = c$note)
  print(out, row.names = FALSE)
  invisible(x)
}

# Normalise the user's labels; unlabelled controls are recognised by name.
sample_groups <- function(sample, groups, ntc, positive) {
  if (is.data.frame(groups)) {
    if (!all(c("sample", "group") %in% names(groups))) {
      stop("`groups` needs the columns `sample` and `group`.", call. = FALSE)
    }
    groups <- stats::setNames(as.character(groups$group),
                              as.character(groups$sample))
  }
  if (is.null(names(groups))) {
    stop("`groups` must be named by sample.", call. = FALSE)
  }
  g <- tolower(trimws(groups[sample]))
  g <- ifelse(grepl("^(cancer|tumou?r|case|patient)", g), "cancer",
       ifelse(grepl("^(control|healthy|normal|non.?cancer|benign)", g), "control",
       ifelse(grepl("^(positive|pc|pos)", g), "positive",
       ifelse(grepl("^(ntc|nk|no.?template|water|blank)", g), "ntc",
       ifelse(is.na(g), NA_character_, "exclude")))))
  guess <- ifelse(matches(sample, ntc), "ntc",
                  ifelse(matches(sample, positive), "positive", "exclude"))
  ifelse(is.na(g), guess, g)
}

# One Ct per sample and gene (wells pooled over runs), plus the reference Ct.
per_sample_ct <- function(wells, ref_ct_max) {
  key <- paste(wells$sample, wells$target, sep = "\r")
  idx <- split(seq_len(nrow(wells)), key)
  rows <- lapply(idx, function(i) {
    w <- wells[i, ]
    n_ct <- sum(!is.na(w$ct))
    data.frame(sample = w$sample[1], group = w$group[1], target = w$target[1],
               is_reference = w$is_reference[1],
               ct = if (n_ct > 0 && n_ct >= nrow(w) / 2) mean(w$ct, na.rm = TRUE)
                    else NA_real_,
               n_wells = nrow(w), n_with_ct = n_ct,
               runs = paste(unique(w$run), collapse = " | "),
               stringsAsFactors = FALSE)
  })
  s <- do.call(rbind, rows)
  ref <- s[s$is_reference, ]
  s <- s[!s$is_reference, ]
  s$ref_ct <- ref$ct[match(s$sample, ref$sample)]
  has_ref <- s$sample %in% ref$sample
  s$ref_failed <- ifelse(has_ref, is.na(s$ref_ct) | s$ref_ct > ref_ct_max, NA)
  s$is_reference <- NULL
  s <- s[order(s$target, match(s$group, c("cancer", "control", "positive",
                                          "ntc", "exclude")), s$sample), ]
  rownames(s) <- NULL
  s
}

gene_cutoff <- function(gene, s, method, target, n_cycles, fixed = NULL) {
  used <- s[s$used, ]
  ca <- used$ct[used$group == "cancer"]
  co <- used$ct[used$group == "control"]
  pos_ct <- safe_mean(s$ct[s$group == "positive"])
  ntc_ct <- suppressWarnings(min(c(s$ct[s$group == "ntc"], Inf), na.rm = TRUE))
  if (is.infinite(ntc_ct)) ntc_ct <- NA_real_

  empty <- data.frame(target = gene, cutoff = NA_real_,
                      suggested_cutoff = NA_real_, source = "suggested",
                      sensitivity = NA_real_,
                      specificity = NA_real_, youden = NA_real_, auc = NA_real_,
                      n_cancer = length(ca), n_control = length(co),
                      tp = NA_integer_, fn = NA_integer_, fp = NA_integer_,
                      tn = NA_integer_, positive_ct = pos_ct, ntc_ct = ntc_ct,
                      note = "", stringsAsFactors = FALSE)
  if (!length(ca) || !length(co)) {
    empty$note <- "needs at least one cancer and one control sample"
    if (!is.null(fixed)) {
      empty$cutoff <- fixed
      empty$source <- "chosen"
    }
    return(list(cutoff = empty, roc = NULL))
  }

  roc <- roc_points(ca, co, n_cycles)
  best <- pick_cutoff(roc, method, target)
  r <- roc[best, ]
  out <- empty
  out$suggested_cutoff <- r$cutoff
  if (!is.null(fixed)) {
    r <- data.frame(cutoff = fixed,
                    sensitivity = mean(!is.na(ca) & ca <= fixed),
                    specificity = mean(is.na(co) | co > fixed))
    out$source <- "chosen"
  }
  out$cutoff <- r$cutoff
  out$sensitivity <- r$sensitivity
  out$specificity <- r$specificity
  out$youden <- r$sensitivity + r$specificity - 1
  out$auc <- auc(ca, co)
  out$tp <- sum(!is.na(ca) & ca <= r$cutoff)
  out$fn <- length(ca) - out$tp
  out$fp <- sum(!is.na(co) & co <= r$cutoff)
  out$tn <- length(co) - out$fp
  out$note <- join_flags(
    if (min(length(ca), length(co)) < 10)
      sprintf("few samples (%d cancer, %d control): cutoff is uncertain",
              length(ca), length(co)) else "",
    if (is.na(r$cutoff)) "no cutoff separates the groups: no Ct is called methylated" else "",
    if (!is.na(r$cutoff) && !is.na(pos_ct) && pos_ct > r$cutoff)
      sprintf("positive control (Ct %.1f) is above the cutoff", pos_ct) else "",
    if (!is.na(r$cutoff) && !is.na(ntc_ct) && ntc_ct <= r$cutoff)
      sprintf("NK control (Ct %.1f) is at or below the cutoff", ntc_ct) else "",
    if (!is.na(out$auc) && out$auc < 0.6) "the gene hardly separates cancer from control" else ""
  )
  roc$target <- gene
  list(cutoff = out, roc = roc[, c("target", "cutoff", "sensitivity",
                                   "specificity")])
}

# Every distinct cutoff: "methylated if Ct <= cutoff". The first row (NA)
# calls nothing methylated. Cutoffs sit halfway to the next observed Ct.
roc_points <- function(ca, co, n_cycles) {
  v <- sort(unique(c(ca, co)[!is.na(c(ca, co))]))
  nxt <- c(v[-1], NA)
  cut <- ifelse(is.na(nxt), pmin(v + 1, pmax(n_cycles, v)), (v + nxt) / 2)
  cut <- round(cut, 2)
  sens <- vapply(v, function(k) mean(!is.na(ca) & ca <= k), numeric(1))
  spec <- vapply(v, function(k) mean(is.na(co) | co > k), numeric(1))
  data.frame(cutoff = c(NA_real_, cut), sensitivity = c(0, sens),
             specificity = c(1, spec))
}

pick_cutoff <- function(roc, method, target) {
  ok <- switch(method,
               youden = rep(TRUE, nrow(roc)),
               specificity = roc$specificity >= target - 1e-9,
               sensitivity = roc$sensitivity >= target - 1e-9)
  if (!any(ok)) ok <- roc$sensitivity == max(roc$sensitivity)
  score <- switch(method,
                  youden = roc$sensitivity + roc$specificity,
                  specificity = roc$sensitivity,
                  sensitivity = roc$specificity)
  score[!ok] <- -Inf
  cand <- which(score >= max(score) - 1e-9)
  cand <- cand[roc$specificity[cand] == max(roc$specificity[cand])]
  cand[1]
}

# Probability that a cancer sample has a lower Ct than a control (no Ct
# counts as the highest Ct; ties count half).
auc <- function(ca, co) {
  a <- ifelse(is.na(ca), Inf, ca)
  b <- ifelse(is.na(co), Inf, co)
  mean(outer(a, b, function(p, q) (p < q) + 0.5 * (p == q)))
}

#' Sensitivity and specificity of a gene panel
#'
#' Calls a sample "cancer" when at least `min_methylated_genes` genes have a
#' Ct at or below their cutoff, and compares that with the known groups.
#'
#' @param x A `qmsp_cutoffs` from [find_cutoffs()].
#' @param cutoffs Named vector of cutoffs per gene (default: the ones found).
#' @param min_methylated_genes Methylated genes needed to call cancer.
#' @return A list with `summary` (sensitivity, specificity and counts) and
#'   `samples` (one row per sample with its methylated genes).
#' @export
panel_performance <- function(x, cutoffs = NULL, min_methylated_genes = 1) {
  if (is.null(cutoffs)) {
    cutoffs <- stats::setNames(x$cutoffs$cutoff, x$cutoffs$target)
  }
  s <- x$samples[x$samples$used & x$samples$target %in% names(cutoffs), ]
  s$methylated <- !is.na(s$ct) & !is.na(cutoffs[s$target]) &
    s$ct <= cutoffs[s$target]
  ids <- unique(s[, c("sample", "group")])
  ids$methylated_genes <- vapply(ids$sample, function(id) {
    paste(s$target[s$sample == id & s$methylated], collapse = ", ")
  }, character(1))
  ids$n_methylated <- vapply(ids$sample, function(id) {
    sum(s$methylated[s$sample == id])
  }, numeric(1))
  ids$called_cancer <- ids$n_methylated >= min_methylated_genes
  ca <- ids$called_cancer[ids$group == "cancer"]
  co <- ids$called_cancer[ids$group == "control"]
  list(summary = data.frame(
         min_methylated_genes = min_methylated_genes,
         sensitivity = if (length(ca)) mean(ca) else NA_real_,
         specificity = if (length(co)) mean(!co) else NA_real_,
         tp = sum(ca), fn = sum(!ca), fp = sum(co), tn = sum(!co)),
       samples = ids)
}

#' Save and load Ct cutoffs
#'
#' `write_cutoffs()` writes a CSV with one row per gene (`gene`, `ct_cutoff`
#' and the statistics). Genes without cancer or control samples are left out. `read_cutoffs()` reads such a file back as the
#' `ct_cutoff` argument of [analyze_qmsp()]: `c(40, TAC1 = 36.5, ...)`, where
#' genes not in the file use `default`. A gene listed with an empty cutoff
#' (no cutoff separated the groups) is never called methylated.
#'
#' @param x A `qmsp_cutoffs` from [find_cutoffs()], or a data frame with
#'   columns `gene` (or `target`) and `ct_cutoff` (or `cutoff`).
#' @param file Path of the CSV file.
#' @param default Cutoff for genes that are not in the file.
#' @return `write_cutoffs()`: `file`, invisibly. `read_cutoffs()`: a named
#'   numeric vector.
#' @export
write_cutoffs <- function(x, file = "qmsp_cutoffs.csv") {
  c <- if (inherits(x, "qmsp_cutoffs")) x$cutoffs else x
  # Genes without enough labelled samples are left out (they keep the default).
  if (all(c("n_cancer", "n_control") %in% names(c))) {
    c <- c[(c$n_cancer > 0 & c$n_control > 0) | !is.na(c$cutoff), , drop = FALSE]
  }
  names(c)[names(c) == "target"] <- "gene"
  names(c)[names(c) == "cutoff"] <- "ct_cutoff"
  utils::write.csv(c, file, row.names = FALSE, na = "")
  invisible(file)
}

#' @rdname write_cutoffs
#' @export
read_cutoffs <- function(file, default = 40) {
  d <- utils::read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  gene <- d[[intersect(c("gene", "target"), names(d))[1]]]
  cut <- suppressWarnings(as.numeric(
    d[[intersect(c("ct_cutoff", "cutoff"), names(d))[1]]]))
  if (is.null(gene) || is.null(cut)) {
    stop("The file needs the columns `gene` and `ct_cutoff`.", call. = FALSE)
  }
  # A gene without a cutoff never counts as methylated.
  cut[is.na(cut)] <- -Inf
  c(default, stats::setNames(cut, gene))
}

#' Export cutoffs to Excel
#'
#' Sheets: the cutoff per gene, every sample's Ct and group, every cutoff
#' tried (ROC), and the settings. CSV files are written when `writexl` is not
#' installed.
#'
#' @param x A `qmsp_cutoffs` from [find_cutoffs()].
#' @param path Output `.xlsx` path.
#' @param chosen Optional named vector of the cutoffs finally chosen per gene
#'   (e.g. after adjusting them by hand); added as column `chosen_cutoff`.
#' @return The path(s) written, invisibly.
#' @export
export_cutoffs <- function(x, path = "qmsp_cutoffs.xlsx", chosen = NULL) {
  c <- x$cutoffs
  if (!is.null(chosen)) c$chosen_cutoff <- unname(chosen[c$target])
  sheets <- list(
    Cutoffs = c,
    Samples = x$samples,
    ROC = x$roc %||% data.frame(),
    Settings = data.frame(setting = names(x$settings),
                          value = vapply(x$settings, function(v)
                            if (is.null(v)) "" else paste(v, collapse = ", "),
                            character(1)))
  )
  if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(sheets, path)
    return(invisible(path))
  }
  dir <- sub("\\.xlsx$", "", path, ignore.case = TRUE)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  files <- file.path(dir, paste0(names(sheets), ".csv"))
  Map(function(d, f) utils::write.csv(d, f, row.names = FALSE), sheets, files)
  invisible(files)
}

#' Plot sample Cts per gene with the cutoffs
#'
#' One panel per gene: every sample's Ct by group, samples without a Ct at the
#' top, and the cutoff as a dashed line.
#'
#' @param x A `qmsp_cutoffs` from [find_cutoffs()].
#' @param cutoffs Named vector of cutoffs to draw (default: the ones found).
#' @return A ggplot object.
#' @export
plot_cutoffs <- function(x, cutoffs = NULL) {
  if (is.null(cutoffs)) {
    cutoffs <- stats::setNames(x$cutoffs$cutoff, x$cutoffs$target)
  }
  s <- x$samples[x$samples$group != "exclude", ]
  top <- max(c(45, s$ct), na.rm = TRUE) + 2
  s$y <- ifelse(is.na(s$ct), top, s$ct)
  s$group <- factor(s$group, levels = c("cancer", "control", "positive", "ntc"),
                    labels = c("Cancer", "Control", "Positive control", "NK"))
  s$shape <- ifelse(s$used | !s$group %in% c("Cancer", "Control"), "used",
                    "left out")
  cut_df <- data.frame(target = names(cutoffs), cutoff = unname(cutoffs))
  cut_df <- cut_df[!is.na(cut_df$cutoff), ]
  ggplot2::ggplot(s, ggplot2::aes(.data$group, .data$y, colour = .data$group)) +
    ggplot2::geom_jitter(ggplot2::aes(shape = .data$shape), width = 0.15,
                         height = 0, size = 2.4) +
    ggplot2::geom_hline(data = cut_df, ggplot2::aes(yintercept = .data$cutoff),
                        linetype = "dashed", colour = "#c0392b") +
    ggplot2::facet_wrap(~target) +
    ggplot2::scale_y_reverse(breaks = c(seq(15, top - 2, 5), top),
                             labels = c(seq(15, top - 2, 5), "no Ct")) +
    ggplot2::scale_colour_manual(values = c(Cancer = "#c0392b",
                                            Control = "#1f6fb4",
                                            "Positive control" = "#7d3c98",
                                            NK = "#555555"),
                                 drop = FALSE) +
    ggplot2::scale_shape_manual(values = c(used = 16, "left out" = 1),
                                name = NULL) +
    ggplot2::labs(x = NULL, y = "Ct", colour = NULL) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "bottom",
                   axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}
