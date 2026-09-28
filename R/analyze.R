#' Analyse a qMSP run
#'
#' Turns raw well results into qMSP calls with one fixed rule:
#' \strong{a gene is methylated when its Ct is 40 or less}, where the Ct is
#' read where the amplification curve crosses a delta Rn of \strong{10,000}.
#'
#' The steps are:
#' \enumerate{
#'   \item every well is labelled as no-template control (NTC), positive
#'     control or sample, using the sample name (or Task = "NTC");
#'   \item the Ct of every well is read from its amplification curve at the
#'     threshold (delta Rn 10,000): the cycle where the curve crosses the
#'     threshold and stays above it (log-linear interpolation). Curves that
#'     never reach it have no Ct;
#'   \item a well is `Positive` when its Ct is at or below the Ct cutoff (40),
#'     otherwise `Negative`;
#'   \item NTC and positive controls are checked per run and gene;
#'   \item replicates are merged per run, sample and gene, and every gene gets
#'     a call: `Methylated`, `Unmethylated` or `Not determined`.
#' }
#'
#' A gene is `Not determined` when the call cannot be trusted: the sample's
#' reference gene (e.g. beta-actin) has no Ct at or below `ref_ct_max`, the
#' no-template control of that gene amplified (for a methylated call), the
#' positive control of that gene did not amplify (for an unmethylated call),
#' or fewer than half of the replicates agree.
#'
#' Every sample then gets a verdict (see [risk_report()]): `Potential cancer`,
#' `Not determined` or `Low risk`.
#'
#' Delta Ct, PMR and beta are also computed for information; they do not
#' change the call.
#'
#' @param x A `qmsp_run` from [read_eds()] / [read_eds_files()].
#' @param reference Regular expression matching the reference gene name
#'   (case-insensitive). Set to `NULL` when the run has no reference gene.
#' @param ntc Regular expression matching no-template control sample names.
#' @param positive Regular expression matching positive control sample names
#'   (e.g. fully methylated / bisulfite-converted cell line DNA).
#' @param ct_cutoff A gene is methylated when its Ct is at or below this.
#'   Lab standard: 40.
#' @param threshold Fluorescence (delta Rn) threshold at which the Ct is read.
#'   Lab standard: 10,000, for every gene including the reference gene.
#'   `NULL` uses the Ct values of the instrument software instead.
#' @param ref_ct_max A sample whose reference gene Ct is above this (or
#'   undetermined) is not determined.
#' @param ref_ct_warn Reference Ct above this is flagged as low DNA input.
#' @param min_methylated_genes Number of methylated panel genes needed for a
#'   `Potential cancer` verdict.
#' @param panel Genes that count for the verdict (default: all genes except the
#'   reference).
#'
#' @return A `qmsp_result` list with `report` (one verdict per run and
#'   sample), `results` (one row per run, sample and gene), `controls` (NTC /
#'   positive control status per run and gene), `wells` (every well with its
#'   role, Ct and result), `curves`, `thresholds` and `settings`.
#' @export
analyze_qmsp <- function(x,
                         reference = "ACTB|B.?ACTIN|BETA.?ACTIN",
                         ntc = "NTC|dH2O|dH20|water|blank|^NK",
                         positive = "H460|A549|HT29|positive|^PC\\b",
                         ct_cutoff = 40,
                         threshold = 10000,
                         ref_ct_max = 40,
                         ref_ct_warn = 35,
                         min_methylated_genes = 1,
                         panel = NULL) {
  if (!inherits(x, "qmsp_run")) {
    stop("`x` must come from read_eds() or read_eds_files().", call. = FALSE)
  }
  settings <- list(reference = reference, ntc = ntc, positive = positive,
                   ct_cutoff = ct_cutoff, threshold = threshold,
                   ref_ct_max = ref_ct_max, ref_ct_warn = ref_ct_warn,
                   min_methylated_genes = min_methylated_genes, panel = panel)

  wells <- apply_thresholds(x$wells, x$curves, threshold)
  wells <- classify_wells(wells, settings)
  controls <- control_status(wells)
  results <- sample_results(wells, controls, settings)

  thresholds <- x$thresholds
  thresholds <- thresholds[!grepl("DEFAULT_SETTINGS", thresholds$target), ,
                           drop = FALSE]
  thresholds <- used_thresholds(thresholds, wells, threshold)

  out <- structure(
    list(report = NULL, results = results, controls = controls, wells = wells,
         curves = x$curves, thresholds = thresholds, meta = x$meta,
         settings = settings),
    class = "qmsp_result"
  )
  out$report <- risk_report(out, min_methylated_genes = min_methylated_genes,
                            panel = panel)
  out
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
  v <- table(x$report$verdict)
  cat("  verdicts: ", paste(names(v), v, sep = " = ", collapse = ", "), "\n",
      sep = "")
  bad <- x$controls[x$controls$ntc_status == "Fail" |
                      x$controls$positive_status == "Fail", ]
  if (nrow(bad)) {
    cat("  control problems:\n")
    for (i in seq_len(nrow(bad))) {
      cat("    ", bad$run[i], " / ", bad$target[i], ": NTC ",
          bad$ntc_status[i], ", positive control ", bad$positive_status[i],
          "\n", sep = "")
    }
  }
  invisible(x)
}

call_levels <- function() c("Methylated", "Unmethylated", "Not determined")

# ---- step 0: own fluorescence thresholds ------------------------------------

# Per-gene threshold lookup; genes without a value get NA (= instrument Ct).
threshold_value <- function(threshold, genes) {
  if (is.null(threshold) || !length(threshold)) {
    return(rep(NA_real_, length(genes)))
  }
  nm <- names(threshold)
  if (!is.null(nm) && all(nzchar(nm))) threshold <- c(NA_real_, threshold)
  as.numeric(gene_value(threshold, genes))
}

# Cycle where a delta Rn curve crosses `thr` for the last time from below
# (so early noise spikes are ignored), interpolated on the log scale.
# NA when the curve ends below the threshold.
ct_at_threshold <- function(cycle, delta_rn, thr) {
  o <- order(cycle)
  cycle <- cycle[o]
  d <- delta_rn[o]
  n <- length(d)
  if (!n || is.na(thr) || is.na(d[n]) || d[n] < thr) return(NA_real_)
  below <- which(is.na(d) | d < thr)
  if (!length(below)) return(cycle[1])
  j <- max(below)
  if (is.na(d[j])) return(cycle[j + 1])
  step <- cycle[j + 1] - cycle[j]
  frac <- if (d[j] > 0) (log(thr) - log(d[j])) / (log(d[j + 1]) - log(d[j]))
    else (thr - d[j]) / (d[j + 1] - d[j])
  cycle[j] + step * frac
}

# Recompute Ct for genes that have their own threshold; keep the instrument's
# value in `ct_instrument`.
apply_thresholds <- function(wells, curves, threshold) {
  wells$ct_instrument <- wells$ct
  wells$threshold <- threshold_value(threshold, wells$target)
  todo <- which(!is.na(wells$threshold))
  if (!length(todo) || is.null(curves) || !nrow(curves)) return(wells)
  curve_key <- paste(curves$run, curves$well_index, curves$target, sep = "\r")
  by_well <- split(seq_len(nrow(curves)), curve_key)
  well_key <- paste(wells$run, wells$well_index, wells$target, sep = "\r")
  for (i in todo) {
    rows <- by_well[[well_key[i]]]
    if (is.null(rows)) next
    wells$ct[i] <- ct_at_threshold(curves$cycle[rows], curves$delta_rn[rows],
                                   wells$threshold[i])
  }
  wells
}

# Threshold table per run and gene: the instrument's and the one used.
used_thresholds <- function(thresholds, wells, threshold) {
  groups <- unique(wells[, c("run", "target")])
  inst <- thresholds$threshold[match(paste(groups$run, groups$target),
                                     paste(thresholds$run, thresholds$target))]
  auto <- thresholds$auto_threshold[match(paste(groups$run, groups$target),
                                          paste(thresholds$run,
                                                thresholds$target))]
  own <- threshold_value(threshold, groups$target)
  out <- data.frame(run = groups$run, target = groups$target,
                    threshold = ifelse(is.na(own), inst, own),
                    instrument_threshold = inst,
                    auto_threshold = auto,
                    source = ifelse(is.na(own), "instrument", "set"),
                    stringsAsFactors = FALSE)
  out[!is.na(out$threshold), , drop = FALSE]
}

# ---- step 1 + 2: wells ------------------------------------------------------

classify_wells <- function(wells, s) {
  is_ref <- if (is.null(s$reference)) rep(FALSE, nrow(wells)) else
    matches(wells$target, s$reference)
  role <- ifelse(tolower(wells$task) == "ntc" | matches(wells$sample, s$ntc),
                 "ntc",
                 ifelse(matches(wells$sample, s$positive), "positive", "sample"))
  wells$role <- role
  wells$is_reference <- is_ref

  amplified <- !is.na(wells$ct)
  cutoff <- ifelse(is_ref, s$ref_ct_max, gene_value(s$ct_cutoff, wells$target))
  within <- amplified & wells$ct <= cutoff
  never <- !amplified & !is.na(wells$ct_instrument) & !is.na(wells$threshold)
  wells$flags <- ifelse(amplified & !within,
                        sprintf("Ct above cutoff (%.1f)", wells$ct),
                 ifelse(never, sprintf(paste("curve does not reach the threshold",
                                             "(instrument Ct %.1f)"),
                                       wells$ct_instrument), ""))
  wells$result <- ifelse(within, "Positive", "Negative")
  wells
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
        if (any(ntc$result == "Positive")) "Fail" else "Pass",
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
    # Ct of the positive wells; if none, of the wells that crossed after the
    # cutoff (shown for information).
    ct_used <- if (n_pos) w$ct[w$result == "Positive"] else w$ct
    gene_call <- if (n_pos == 0) "Negative"
      else if (n_pos >= nrow(w) / 2) "Positive" else "Mixed"

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
      n_wells = nrow(w), n_positive = n_pos,
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
  ref_bad <- res$ref_status %in% c("Failed", "Missing")
  # A contaminated NTC puts methylated calls in doubt; a failed positive
  # control puts unmethylated calls in doubt. Controls are judged as they are.
  ntc_doubt <- res$role == "sample" & res$gene_call == "Positive" &
    res$ntc_status == "Fail"
  pc_doubt <- res$role == "sample" & res$gene_call == "Negative" &
    res$positive_status == "Fail"
  res$call <- ifelse(ref_bad | res$gene_call == "Mixed" | ntc_doubt | pc_doubt,
                     "Not determined",
                     ifelse(res$gene_call == "Positive", "Methylated",
                            "Unmethylated"))

  res$delta_ct <- ifelse(ref_ok, res$ct - res$ref_ct, NA_real_)
  res$ratio <- ifelse(ref_ok & res$gene_call == "Negative", 0,
                      2^-res$delta_ct)

  # PMR relative to the positive controls of the same run and gene.
  pc <- res$role == "positive" & res$gene_call == "Positive" & ref_ok &
    !is.na(res$ratio)
  pc_ratio <- tapply(res$ratio[pc], paste(res$run, res$target)[pc], mean)
  res$pmr <- 100 * res$ratio / unname(pc_ratio[paste(res$run, res$target)])

  # Beta: methylation level 0-1 relative to the fully methylated control.
  no_ref <- res$ref_status %in% c("Not run", "n/a")
  beta_pc <- ifelse(res$gene_call == "Negative", 0,
                    2^-(res$ct - ctrl$positive_mean_ct))
  res$beta <- pmin(1, ifelse(ref_ok, res$pmr / 100,
                             ifelse(no_ref, beta_pc, NA_real_)))
  res$beta_method <- ifelse(is.na(res$beta), "",
                            ifelse(ref_ok, "PMR/100 (reference-normalised)",
                                   "vs positive control Ct (no reference)"))
  res$notes <- join_flags(
    ifelse(res$ref_status == "Failed",
           sprintf("reference gene failed (no Ct at or below %g) - repeat sample",
                   s$ref_ct_max), ""),
    ifelse(res$ref_status == "Missing", "no reference gene well for sample", ""),
    ifelse(res$ref_status == "Low input",
           sprintf("low DNA input (reference Ct %.1f)", res$ref_ct), ""),
    ifelse(res$ref_status == "Not run", "no reference gene in run", ""),
    ifelse(res$ntc_status == "Fail", "NTC amplified for this gene", ""),
    ifelse(res$positive_status == "Fail", "positive control failed", ""),
    ifelse(res$n_wells > 1 & res$n_positive > 0 & res$n_positive < res$n_wells,
           sprintf("replicates disagree (%d of %d positive)", res$n_positive,
                   res$n_wells), ""),
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
#'   `"ratio"`, `"pmr"` or `"beta"`.
#' @param controls Include control samples?
#' @return A data frame.
#' @export
results_wide <- function(x, value = c("call", "ct", "delta_ct", "ratio",
                                      "pmr", "beta"),
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
