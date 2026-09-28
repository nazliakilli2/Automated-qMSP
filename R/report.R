#' Per-sample verdict: potential cancer, not determined or low risk
#'
#' Counts, for every sample, how many genes of the panel are methylated
#' (Ct at or below 40 at a delta Rn threshold of 10,000).
#' \itemize{
#'   \item `Potential cancer`: at least `min_methylated_genes` genes are
#'     methylated.
#'   \item `Not determined`: fewer are methylated, but genes that could not be
#'     determined (failed reference gene, failed controls, replicates that
#'     disagree) could change that. Repeat the sample.
#'   \item `Low risk`: fewer than `min_methylated_genes` genes are methylated,
#'     even counting the ones that could not be determined.
#' }
#' This is a research tool, not a diagnosis.
#'
#' @param x A `qmsp_result` from [analyze_qmsp()].
#' @param min_methylated_genes Number of methylated genes needed for
#'   `Potential cancer`.
#' @param panel Genes that count (default: all genes except the reference).
#' @return A data frame with one row per run and sample.
#' @export
risk_report <- function(x, min_methylated_genes = 1, panel = NULL) {
  r <- x$results[x$results$role == "sample", ]
  if (length(panel)) r <- r[r$target %in% panel, ]
  if (!nrow(r)) {
    return(data.frame(run = character(), sample = character(),
                      verdict = factor(character(), levels = verdict_levels()),
                      reason = character(), stringsAsFactors = FALSE))
  }
  ids <- unique(r[, c("run", "sample")])
  rows <- lapply(seq_len(nrow(ids)), function(i) {
    s <- r[r$run == ids$run[i] & r$sample == ids$sample[i], ]
    call <- as.character(s$call)
    meth <- s$target[call == "Methylated"]
    n_meth <- length(meth)
    nd <- s$target[call == "Not determined"]
    verdict <- if (n_meth >= min_methylated_genes) "Potential cancer"
      else if (n_meth + length(nd) >= min_methylated_genes) "Not determined"
      else "Low risk"
    meth_ct <- s$ct[call == "Methylated"]
    reason <- join_flags(
      sprintf("%d of %d genes methylated%s", n_meth, nrow(s),
              if (n_meth) paste0(" (", paste(sprintf("%s Ct %.1f", meth, meth_ct),
                                             collapse = ", "), ")") else ""),
      if (length(nd)) paste0("not determined: ", paste(nd, collapse = ", ")) else "",
      if (any(s$ref_status == "Failed")) "reference gene failed - repeat sample" else "",
      if (any(s$ref_status == "Missing")) "no reference gene well" else "",
      if (any(s$ref_status == "Low input")) "low DNA input" else ""
    )
    data.frame(
      run = ids$run[i], sample = ids$sample[i], verdict = verdict,
      reason = reason, n_genes = nrow(s), n_methylated = n_meth,
      methylated_genes = paste(meth, collapse = ", "),
      n_not_determined = length(nd),
      lowest_ct = if (all(is.na(s$ct))) NA_real_ else min(s$ct, na.rm = TRUE),
      ref_ct = s$ref_ct[1],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$verdict <- factor(out$verdict, levels = verdict_levels())
  out <- out[order(out$verdict, out$run, out$sample), ]
  rownames(out) <- NULL
  out
}

verdict_levels <- function() c("Potential cancer", "Not determined", "Low risk")

#' Write a printable HTML report
#'
#' A single self-contained HTML file with the verdict for every sample, the
#' Ct of every gene, the control checks and the settings used. Open it in a browser. To get a PDF, print it
#' and choose "Save as PDF".
#'
#' @param x A `qmsp_result` from [analyze_qmsp()].
#' @param file Output `.html` path.
#' @param title Report title.
#' @return `file`, invisibly.
#' @export
write_report <- function(x, file = "qmsp_report.html", title = "qMSP report") {
  writeLines(report_html(x, title), file, useBytes = TRUE)
  invisible(file)
}

report_html <- function(x, title = "qMSP report") {
  if (!inherits(x, "qmsp_result")) {
    stop("`x` must come from analyze_qmsp().", call. = FALSE)
  }
  rep <- x$report
  counts <- table(factor(rep$verdict, levels = verdict_levels()))
  multi_run <- length(unique(rep$run)) > 1

  verdict_rows <- vapply(seq_len(nrow(rep)), function(i) {
    v <- as.character(rep$verdict[i])
    tr(c(esc(rep$sample[i]), if (multi_run) esc(rep$run[i]),
         sprintf('<span class="badge %s">%s</span>', verdict_class(v), esc(v)),
         esc(rep$reason[i]), fmt_fixed(rep$ref_ct[i])))
  }, character(1))

  s <- x$settings
  settings <- data.frame(
    Setting = c("Rule", "Fluorescence threshold (delta Rn)",
                "Ct cutoff (methylated if at or below)", "Reference gene",
                "Reference gene Ct max / low-input warning",
                "NTC names", "Positive control names",
                "Methylated genes needed for 'Potential cancer'",
                "Genes in panel"),
    Value = c("A gene is methylated when its Ct is at or below the Ct cutoff",
              if (length(s$threshold)) fmt_setting(s$threshold)
              else "instrument software",
              fmt_setting(s$ct_cutoff), s$reference %||% "none",
              paste(s$ref_ct_max, "/", s$ref_ct_warn), s$ntc, s$positive,
              s$min_methylated_genes,
              if (length(s$panel)) paste(s$panel, collapse = ", ") else
                "all genes except the reference"),
    stringsAsFactors = FALSE
  )

  ctrl <- x$controls
  ctrl_rows <- vapply(seq_len(nrow(ctrl)), function(i) {
    tr(c(if (multi_run) esc(ctrl$run[i]), esc(ctrl$target[i]),
         status_cell(ctrl$ntc_status[i]), status_cell(ctrl$positive_status[i]),
         fmt_fixed(ctrl$positive_mean_ct[i])))
  }, character(1))

  runs <- if (is.data.frame(x$meta)) x$meta else as.data.frame(x$meta)
  run_list <- paste0("<li>", esc(runs$run),
                     ifelse(is.na(runs$run_start), "",
                            paste0(" &middot; ", format(runs$run_start,
                                                        "%d %b %Y"))),
                     "</li>", collapse = "")

  paste0(
    '<!doctype html><html><head><meta charset="utf-8">',
    '<meta name="viewport" content="width=device-width, initial-scale=1">',
    "<title>", esc(title), "</title><style>", report_css(), "</style></head>",
    "<body><main>",
    "<h1>", esc(title), "</h1>",
    '<p class="muted">Generated ', format(Sys.time(), "%d %b %Y %H:%M"),
    " with autoqmsp</p>",
    "<ul class=\"runs\">", run_list, "</ul>",
    '<p class="disclaimer">For research use only. These results are not a ',
    "diagnosis. A &ldquo;Potential cancer&rdquo; result means that ",
    "methylation of cancer-associated genes was detected. It must be confirmed ",
    "by a clinician with standard diagnostic tests.</p>",
    '<p class="rule">A gene is <b>methylated</b> when its Ct is <b>',
    esc(fmt_setting(s$ct_cutoff)), " or less</b>, with the Ct read where the ",
    "amplification curve crosses <b>&Delta;Rn ",
    if (length(s$threshold)) esc(fmt_thr(s$threshold)) else "(instrument threshold)",
    "</b>.</p>",
    '<div class="tiles">',
    tile(counts[["Potential cancer"]], "Potential cancer", "bad"),
    tile(counts[["Not determined"]], "Not determined (repeat)", "warn"),
    tile(counts[["Low risk"]], "Low risk", "good"),
    "</div>",
    "<h2>Result per sample</h2>",
    table_html(c("Sample", if (multi_run) "Run", "Result", "Why",
                 "Reference Ct"), verdict_rows),
    "<h2>Ct per gene</h2>",
    '<p class="muted">Red = methylated (Ct at or below the cutoff); ',
    "&ndash; = the curve did not reach the threshold; ",
    "yellow (n.d.) = not determined.</p>",
    ct_table(x, multi_run),
    "<h2>Control checks</h2>",
    table_html(c(if (multi_run) "Run", "Gene", "No-template control",
                 "Positive control", "Positive control Ct"), ctrl_rows),
    "<h2>Settings used</h2>",
    table_html(c("Setting", "Value"),
               vapply(seq_len(nrow(settings)), function(i) {
                 tr(c(esc(settings$Setting[i]), esc(settings$Value[i])))
               }, character(1))),
    "</main></body></html>"
  )
}

ct_table <- function(x, multi_run) {
  r <- x$results[x$results$role == "sample", ]
  if (length(x$settings$panel)) r <- r[r$target %in% x$settings$panel, ]
  ids <- unique(r[, c("run", "sample")])
  genes <- sort(unique(r$target))
  rows <- vapply(seq_len(nrow(ids)), function(i) {
    cells <- vapply(genes, function(g) {
      m <- r[r$run == ids$run[i] & r$sample == ids$sample[i] & r$target == g, ]
      if (!nrow(m)) return('<td class="na"></td>')
      call <- as.character(m$call[1])
      cls <- switch(call, Methylated = "meth", "Not determined" = "nd", "")
      label <- if (call == "Not determined") {
        paste0(if (is.na(m$ct[1])) "" else paste0(fmt_fixed(m$ct[1]), " "),
               "n.d.")
      } else fmt_fixed(m$ct[1])
      sprintf('<td class="%s" title="%s">%s</td>', cls,
              esc(paste0(call, if (nzchar(m$notes[1])) paste0(": ", m$notes[1]))),
              label)
    }, character(1))
    paste0("<tr><td>", esc(ids$sample[i]), "</td>",
           if (multi_run) paste0("<td>", esc(ids$run[i]), "</td>"),
           paste(cells, collapse = ""), "</tr>")
  }, character(1))
  table_html(c("Sample", if (multi_run) "Run", genes), rows)
}

report_css <- function() {
  paste(
    "body{margin:0;background:#f7f7f5;color:#1f2328;",
    "font:15px/1.5 -apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif}",
    "main{max-width:1100px;margin:0 auto;padding:24px 16px 48px}",
    "h1{font-size:26px;margin:0 0 4px}h2{font-size:18px;margin:32px 0 8px}",
    ".muted{color:#656d76;font-size:13px;margin:4px 0}",
    ".runs{font-size:13px;color:#444;padding-left:18px}",
    ".disclaimer{background:#fff8e1;border-left:4px solid #f5b400;",
    "padding:10px 12px;font-size:13px}",
    ".tiles{display:flex;gap:12px;flex-wrap:wrap;margin:16px 0}",
    ".tile{flex:1 1 180px;background:#fff;border-radius:8px;padding:14px;",
    "border-top:5px solid #ccc}.tile b{display:block;font-size:30px}",
    ".tile.bad{border-color:#c0392b}.tile.warn{border-color:#f5b400}",
    ".tile.good{border-color:#2e8b57}",
    ".scroll{overflow-x:auto}",
    "table{border-collapse:collapse;width:100%;background:#fff;font-size:13px}",
    "th,td{border:1px solid #e3e3e0;padding:6px 8px;text-align:left;",
    "vertical-align:top}th{background:#efefec}",
    ".rule{background:#fff;border-left:4px solid #1f6fb4;padding:10px 12px}",
    "td.meth{font-weight:700;color:#7b1d13;background:#fbe3e0}",
    "td.nd{background:#fdf1c4}td.na{background:#fafafa}",
    ".badge{display:inline-block;padding:2px 8px;border-radius:10px;",
    "font-weight:600;white-space:nowrap}",
    ".badge.bad{background:#fbe3e0;color:#9b2c1f}",
    ".badge.warn{background:#fdf1c4;color:#7a5b00}",
    ".badge.good{background:#dff3e7;color:#1e6b41}",
    "@media print{body{background:#fff}.tile{border:1px solid #ccc}}"
  )
}

verdict_class <- function(v) {
  switch(v, "Potential cancer" = "bad", "Not determined" = "warn",
         "Low risk" = "good", "")
}

status_cell <- function(status) {
  cls <- switch(status, Pass = "good", Fail = "bad", Review = "warn", "")
  sprintf('<td><span class="badge %s">%s</span></td>', cls, esc(status))
}

tile <- function(n, label, cls) {
  sprintf('<div class="tile %s"><b>%d</b>%s</div>', cls, as.integer(n),
          esc(label))
}

table_html <- function(header, rows) {
  if (!length(rows)) return('<p class="muted">None.</p>')
  paste0('<div class="scroll"><table><thead><tr>',
         paste0("<th>", esc(header), "</th>", collapse = ""),
         "</tr></thead><tbody>", paste(rows, collapse = ""),
         "</tbody></table></div>")
}

# Cells that are already <td> stay as they are; others are wrapped.
tr <- function(cells) {
  cells <- ifelse(startsWith(cells, "<td"), cells,
                  paste0("<td>", cells, "</td>"))
  paste0("<tr>", paste(cells, collapse = ""), "</tr>")
}

esc <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub('"', "&quot;", x, fixed = TRUE)
}

fmt_num <- function(x, digits) {
  ifelse(is.na(x), "&ndash;", formatC(x, format = "fg", digits = digits))
}

fmt_fixed <- function(x) ifelse(is.na(x), "&ndash;", sprintf("%.1f", x))

fmt_thr <- function(v) {
  if (!is.null(names(v))) return(fmt_setting(v))
  formatC(v, format = "fg", big.mark = ",", digits = 6)
}

fmt_setting <- function(v) {
  nm <- names(v)
  if (is.null(nm)) return(paste(format(v), collapse = ", "))
  paste(ifelse(nzchar(nm), paste0(nm, " = ", v), paste0("default ", v)),
        collapse = ", ")
}
