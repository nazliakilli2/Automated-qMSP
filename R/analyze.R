#' Analyse a qMSP run
#'
#' Turns raw well results into qMSP calls. The steps are:
#' \enumerate{
#'   \item every well is labelled as no-template control (NTC), positive
#'     control or sample, using the sample name (or Task = "NTC");
#'   \item every well gets a result -- `Positive`, `Negative` or `Review` --
#'     from its Ct, the instrument's amplification status and Cq confidence,
#'     and the height of its amplification curve compared with the positive
#'     controls of the same gene (low, linear "drift" curves are not trusted);
#'   \item NTC and positive controls are checked per run and gene;
#'   \item replicates are merged per run, sample and gene, the reference gene
#'     (e.g. beta-actin) is checked, and the methylation call, delta Ct,
#'     methylation ratio and PMR are computed.
#' }
#'
#' Methylation values:
#' \itemize{
#'   \item `delta_ct` = Ct(gene) - Ct(reference)
#'   \item `ratio` = 2^-delta_ct (0 for unmethylated samples with a valid
#'     reference)
#'   \item `pmr` = percentage of methylated reference: 100 x ratio(sample) /
#'     mean ratio(positive controls) for the same gene in the same run.
#' }
#'
#' @param x A `qmsp_run` from [read_eds()] / [read_eds_files()].
#' @param reference Regular expression matching the reference gene name
#'   (case-insensitive). Set to `NULL` when the run has no reference gene;
#'   calls are then made from the gene Ct alone.
#' @param ntc Regular expression matching no-template control sample names.
#' @param positive Regular expression matching positive control sample names
#'   (e.g. fully methylated / bisulfite-converted cell line DNA).
#' @param ct_cutoff Genes with Ct above this are called unmethylated.
#' @param ref_ct_max Samples whose reference Ct is above this (or undetermined)
#'   are invalid.
#' @param ref_ct_warn Reference Ct above this is flagged as low DNA input.
#' @param min_cq_conf Amplified wells with an instrument Cq confidence below this
#'   are sent to review.
#' @param min_plateau Amplified wells whose final delta Rn is below this
#'   fraction of the positive controls' (same gene and run) are sent to review.
#'
#' @return A `qmsp_result` list with `results` (one row per run, sample and
#'   gene), `controls` (NTC / positive control status per run and gene),
#'   `wells` (every well with its role, flags and result), `curves`,
#'   `thresholds` and `settings`.
#' @export
analyze_qmsp <- function(x,
                         reference = "ACTB|B.?ACTIN|BETA.?ACTIN",
                         ntc = "NTC|dH2O|dH20|water|blank|^NK",
                         positive = "H460|A549|HT29|positive|^PC\\b",
                         ct_cutoff = 40,
                         ref_ct_max = 40,
                         ref_ct_warn = 35,
                         min_cq_conf = 0.5,
                         min_plateau = 0.2) {
  if (!inherits(x, "qmsp_run")) {
    stop("`x` must come from read_eds() or read_eds_files().", call. = FALSE)
  }
  settings <- list(reference = reference, ntc = ntc, positive = positive,
                   ct_cutoff = ct_cutoff, ref_ct_max = ref_ct_max,
                   ref_ct_warn = ref_ct_warn, min_cq_conf = min_cq_conf,
                   min_plateau = min_plateau)

  wells <- classify_wells(x$wells, x$curves, settings)
  controls <- control_status(wells)
  results <- sample_results(wells, controls, settings)

  thresholds <- x$thresholds
  thresholds <- thresholds[!grepl("DEFAULT_SETTINGS", thresholds$target), ,
                           drop = FALSE]

  structure(
    list(results = results, controls = controls, wells = wells,
         curves = x$curves, thresholds = thresholds, meta = x$meta,
         settings = settings),
    class = "qmsp_result"
  )
}

#' @export
print.qmsp_result <- function(x, ...) {
  r <- x$results[x$results$role == "sample", ]
  cat("<qmsp_result> ", length(unique(x$results$run)), " run(s), ",
      length(unique(paste(r$run, r$sample))), " samples, ",
      length(unique(r$target)), " genes\n", sep = "")
  tab <- table(factor(r$call, levels = call_levels()))
  cat("  calls: ", paste(names(tab), tab, sep = " = ", collapse = ", "), "\n",
      sep = "")
  bad <- x$controls[x$controls$ntc_status %in% c("Fail", "Review") |
                      x$controls$positive_status %in% c("Fail"), ]
  if (nrow(bad)) {
    cat("  control problems:\n")
    for (i in seq_len(nrow(bad))) {
      cat("    ", bad$run[i], " / ", bad$target[i], ": NTC ",
          bad$ntc_status[i], ", positive control ", bad$positive_status[i],
          "\n", sep = "")
    }
  }
  n_review <- sum(x$wells$result == "Review")
  if (n_review) cat("  ", n_review, " well(s) need manual review (see $wells)\n",
                    sep = "")
  invisible(x)
}

call_levels <- function() {
  c("Methylated", "Unmethylated", "Review", "Invalid")
}

# ---- step 1 + 2: wells ------------------------------------------------------

classify_wells <- function(wells, curves, s) {
  is_ref <- if (is.null(s$reference)) rep(FALSE, nrow(wells)) else
    matches(wells$target, s$reference)
  role <- ifelse(tolower(wells$task) == "ntc" | matches(wells$sample, s$ntc),
                 "ntc",
                 ifelse(matches(wells$sample, s$positive), "positive", "sample"))
  wells$role <- role
  wells$is_reference <- is_ref

  wells$plateau <- final_delta_rn(wells, curves)
  wells$plateau_ratio <- plateau_ratio(wells, s$min_cq_conf)

  amplified <- !is.na(wells$ct)
  cutoff <- ifelse(is_ref, s$ref_ct_max, s$ct_cutoff)
  within <- amplified & wells$ct <= cutoff

  f_conf <- within & !is.na(wells$cq_conf) & wells$cq_conf < s$min_cq_conf
  f_status <- within & !is.na(wells$amp_status) & wells$amp_status != 1
  f_plateau <- within & !is.na(wells$plateau_ratio) &
    wells$plateau_ratio < s$min_plateau

  wells$flags <- join_flags(
    ifelse(f_conf, sprintf("low Cq confidence (%.2f)", wells$cq_conf), ""),
    ifelse(f_status, ifelse(wells$amp_status == 0,
                            "instrument: inconclusive amplification",
                            "instrument: no amplification"), ""),
    ifelse(f_plateau, sprintf("weak curve (%.0f%% of positive control)",
                              100 * wells$plateau_ratio), ""),
    ifelse(amplified & !within, sprintf("Ct above cutoff (%.1f)", wells$ct), "")
  )
  wells$result <- ifelse(!within, "Negative",
                         ifelse(f_conf | f_status | f_plateau, "Review",
                                "Positive"))
  wells
}

# Height of each curve: highest delta Rn over the last 5 cycles.
final_delta_rn <- function(wells, curves) {
  if (is.null(curves) || !nrow(curves)) return(rep(NA_real_, nrow(wells)))
  key <- paste(curves$run, curves$well_index, curves$target, sep = "\r")
  last <- stats::ave(curves$cycle, key, FUN = max)
  tail <- curves[curves$cycle > last - 5, ]
  tail_key <- paste(tail$run, tail$well_index, tail$target, sep = "\r")
  top <- tapply(tail$delta_rn, tail_key, max, na.rm = TRUE)
  unname(top[paste(wells$run, wells$well_index, wells$target, sep = "\r")])
}

# Curve height relative to the typical good curve of the same gene and run:
# positive controls if there are any, otherwise all confident wells.
plateau_ratio <- function(wells, min_cq_conf) {
  good <- !is.na(wells$ct) & !is.na(wells$cq_conf) &
    wells$cq_conf >= min_cq_conf & wells$amp_status %in% 1
  group <- paste(wells$run, wells$target, sep = "\r")
  ref <- vapply(seq_len(nrow(wells)), function(i) {
    same <- group == group[i] & good
    pc <- same & wells$role == "positive"
    use <- if (any(pc)) pc else same
    if (!any(use)) NA_real_ else stats::median(wells$plateau[use], na.rm = TRUE)
  }, numeric(1))
  ifelse(is.na(ref) | ref <= 0, NA_real_, wells$plateau / ref)
}

# ---- step 3: controls -------------------------------------------------------

control_status <- function(wells) {
  groups <- unique(wells[, c("run", "target")])
  rows <- lapply(seq_len(nrow(groups)), function(i) {
    w <- wells[wells$run == groups$run[i] & wells$target == groups$target[i], ]
    ntc <- w[w$role == "ntc", ]
    pos <- w[w$role == "positive", ]
    data.frame(
      run = groups$run[i],
      target = groups$target[i],
      is_reference = any(w$is_reference),
      ntc_status = if (!nrow(ntc)) "Not run" else
        if (any(ntc$result == "Positive")) "Fail" else
          if (any(ntc$result == "Review")) "Review" else "Pass",
      ntc_min_ct = suppressWarnings(min(c(ntc$ct, Inf), na.rm = TRUE)),
      positive_status = if (!nrow(pos)) "Not run" else
        if (any(pos$result == "Positive")) "Pass" else "Fail",
      positive_mean_ct = safe_mean(pos$ct[pos$result == "Positive"]),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$ntc_min_ct[is.infinite(out$ntc_min_ct)] <- NA_real_
  rownames(out) <- NULL
  out
}

# ---- step 4: samples --------------------------------------------------------

sample_results <- function(wells, controls, s) {
  # Reference gene per run and sample.
  ref_w <- wells[wells$is_reference, ]
  ref_key <- paste(ref_w$run, ref_w$sample, sep = "\r")
  ref_ct <- tapply(ifelse(ref_w$result == "Positive", ref_w$ct, NA_real_),
                   ref_key, safe_mean)

  gene_w <- wells[!wells$is_reference, ]
  if (!nrow(gene_w)) {
    stop("No gene wells found besides the reference gene.", call. = FALSE)
  }
  groups <- unique(gene_w[, c("run", "sample", "role", "target")])
  runs_with_ref <- unique(ref_w$run)

  rows <- lapply(seq_len(nrow(groups)), function(i) {
    g <- groups[i, ]
    w <- gene_w[gene_w$run == g$run & gene_w$sample == g$sample &
                  gene_w$target == g$target, ]
    n_pos <- sum(w$result == "Positive")
    n_rev <- sum(w$result == "Review")
    ct_used <- w$ct[w$result == "Positive"]
    gene_call <- if (n_pos > 0 && n_pos >= nrow(w) / 2 && n_rev == 0) "Positive"
      else if (n_pos + n_rev == 0) "Negative" else "Review"

    key <- paste(g$run, g$sample, sep = "\r")
    has_ref_well <- key %in% ref_key
    r_ct <- if (has_ref_well) unname(ref_ct[key]) else NA_real_
    ref_status <- if (g$role == "ntc") "n/a"
      else if (!g$run %in% runs_with_ref) "Not run"
      else if (!has_ref_well) "Missing"
      else if (is.na(r_ct)) "Failed"
      else if (r_ct > s$ref_ct_warn) "Low input" else "OK"

    data.frame(
      run = g$run, sample = g$sample, role = g$role, target = g$target,
      n_wells = nrow(w), n_positive = n_pos, n_review = n_rev,
      ct = safe_mean(ct_used), ct_sd = safe_sd(ct_used),
      ref_ct = r_ct, ref_status = ref_status, gene_call = gene_call,
      wells = paste(w$well, collapse = ","),
      well_flags = paste(unique(w$flags[nzchar(w$flags)]), collapse = " | "),
      stringsAsFactors = FALSE
    )
  })
  res <- do.call(rbind, rows)

  ctrl <- controls[match(paste(res$run, res$target),
                         paste(controls$run, controls$target)), ]
  res$ntc_status <- ctrl$ntc_status
  res$positive_status <- ctrl$positive_status

  ref_ok <- res$ref_status %in% c("OK", "Low input")
  # A contaminated NTC puts sample positives in doubt, not the controls.
  ntc_doubt <- res$role == "sample" & res$ntc_status %in% c("Fail", "Review")
  res$call <- ifelse(res$ref_status %in% c("Failed", "Missing"), "Invalid",
              ifelse(res$gene_call == "Negative", "Unmethylated",
              ifelse(res$gene_call == "Review" | ntc_doubt, "Review",
                     "Methylated")))

  res$delta_ct <- ifelse(ref_ok, res$ct - res$ref_ct, NA_real_)
  res$ratio <- ifelse(ref_ok & res$gene_call == "Negative", 0,
                      2^-res$delta_ct)

  # PMR relative to the positive controls of the same run and gene.
  pc <- res$role == "positive" & res$gene_call == "Positive" & ref_ok &
    !is.na(res$ratio)
  pc_ratio <- tapply(res$ratio[pc], paste(res$run, res$target)[pc], mean)
  res$pmr <- 100 * res$ratio / unname(pc_ratio[paste(res$run, res$target)])

  res$notes <- join_flags(
    ifelse(res$ref_status == "Failed",
           "reference gene failed - repeat sample", ""),
    ifelse(res$ref_status == "Missing", "no reference gene well for sample", ""),
    ifelse(res$ref_status == "Low input",
           sprintf("low DNA input (reference Ct %.1f)", res$ref_ct), ""),
    ifelse(res$ref_status == "Not run", "no reference gene in run", ""),
    ifelse(res$ntc_status == "Fail", "NTC amplified for this gene", ""),
    ifelse(res$ntc_status == "Review", "NTC signal needs review", ""),
    ifelse(res$positive_status == "Fail", "positive control failed", ""),
    ifelse(res$n_wells > 1 & res$n_positive > 0 & res$n_positive < res$n_wells,
           "replicates disagree", ""),
    res$well_flags
  )
  res$well_flags <- NULL

  res$call <- factor(res$call, levels = call_levels())
  res <- res[order(res$run, res$role != "sample", res$sample, res$target), ]
  rownames(res) <- NULL
  res
}

#' Wide summary table: one row per sample, one column per gene
#'
#' @param x A `qmsp_result` from [analyze_qmsp()].
#' @param value Which value to show: `"call"`, `"ct"`, `"delta_ct"`,
#'   `"ratio"` or `"pmr"`.
#' @param controls Include control samples?
#' @return A data frame.
#' @export
results_wide <- function(x, value = c("call", "ct", "delta_ct", "ratio", "pmr"),
                         controls = FALSE) {
  value <- match.arg(value)
  r <- x$results
  if (!controls) r <- r[r$role == "sample", ]
  r$value <- if (value == "call") as.character(r$call) else
    signif(r[[value]], 4)
  ids <- unique(r[, c("run", "sample")])
  genes <- unique(r$target)
  out <- ids
  for (g in genes) {
    m <- r[r$target == g, ]
    out[[g]] <- m$value[match(paste(ids$run, ids$sample),
                              paste(m$run, m$sample))]
  }
  ref <- unique(r[, c("run", "sample", "ref_ct")])
  out$reference_ct <- signif(ref$ref_ct[match(paste(ids$run, ids$sample),
                                              paste(ref$run, ref$sample))], 4)
  rownames(out) <- NULL
  out
}
