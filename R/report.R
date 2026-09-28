#' Per-sample verdict: potentially cancer, inconclusive or not risky
#'
#' Counts, for every sample, how many genes of the panel are methylated.
#' \itemize{
#'   \item `Potentially cancer`: at least `min_methylated_genes` genes are
#'     methylated.
#'   \item `Inconclusive`: fewer are methylated, but enough genes need review or
#'     are invalid that the verdict could change. Repeat the sample.
#'   \item `Not risky`: fewer than `min_methylated_genes` genes are methylated,
#'     even counting the ones that need review.
#' }
#' This is a research tool, not a diagnosis.
#'
#' @param x A `qmsp_result` from [analyze_qmsp()].
#' @param min_methylated_genes Number of methylated genes needed for
#'   `Potentially cancer`.
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
    n_unsure <- sum(call %in% c("Review", "Invalid"))
    verdict <- if (n_meth >= min_methylated_genes) "Potentially cancer"
      else if (n_meth + n_unsure >= min_methylated_genes) "Inconclusive"
      else "Not risky"
    reason <- join_flags(
      sprintf("%d of %d genes methylated%s", n_meth, nrow(s),
              if (n_meth) paste0(" (", paste(meth, collapse = ", "), ")") else ""),
      if (any(call == "Review"))
        paste0("needs review: ", paste(s$target[call == "Review"],
                                       collapse = ", ")) else "",
      if (any(call == "Invalid")) "reference gene failed - repeat sample" else "",
      if (any(s$ref_status == "Low input")) "low DNA input" else ""
    )
    data.frame(
      run = ids$run[i], sample = ids$sample[i], verdict = verdict,
      reason = reason, n_genes = nrow(s), n_methylated = n_meth,
      methylated_genes = paste(meth, collapse = ", "),
      n_review = sum(call == "Review"), n_invalid = sum(call == "Invalid"),
      max_beta = if (all(is.na(s$beta))) NA_real_ else max(s$beta, na.rm = TRUE),
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

verdict_levels <- function() c("Potentially cancer", "Inconclusive", "Not risky")

#' Write a printable HTML report
#'
#' A single self-contained HTML file with the verdict for every sample, the
#' beta value of every gene, the control checks, the wells that need a manual
#' look, and the settings used. Open it in a browser. To get a PDF, print it
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
         esc(rep$reason[i]), fmt_num(rep$max_beta[i], 3),
         fmt_fixed(rep$ref_ct[i])))
  }, character(1))

  s <- x$settings
  settings <- data.frame(
    Setting = c("Reference gene", "NTC names", "Positive control names",
                "Gene Ct cutoff (methylated if at or below)",
                "Beta cutoff (methylated if at or above)",
                "Fluorescence threshold (delta Rn)",
                "Reference Ct max / low-input warning",
                "Minimum Cq confidence", "Minimum curve height",
                "Methylated genes needed for 'Potentially cancer'",
                "Genes in panel"),
    Value = c(s$reference %||% "none", s$ntc, s$positive,
              fmt_setting(s$ct_cutoff), fmt_setting(s$beta_cutoff),
              if (length(s$threshold)) paste(fmt_setting(s$threshold),
                                             "(other genes: instrument)")
              else "instrument software",
              paste(s$ref_ct_max, "/", s$ref_ct_warn), s$min_cq_conf,
              paste0(100 * s$min_plateau, "% of positive control"),
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

  review <- x$wells[x$wells$result == "Review", ]
  review_rows <- vapply(seq_len(nrow(review)), function(i) {
    tr(c(if (multi_run) esc(review$run[i]), esc(review$well[i]),
         esc(review$sample[i]), esc(review$target[i]),
         fmt_fixed(review$ct[i]), esc(review$flags[i])))
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
    "diagnosis. A &ldquo;Potentially cancer&rdquo; result means that ",
    "methylation of cancer-associated genes was detected. It must be confirmed ",
    "by a clinician with standard diagnostic tests.</p>",
    '<div class="tiles">',
    tile(counts[["Potentially cancer"]], "Potentially cancer", "bad"),
    tile(counts[["Inconclusive"]], "Inconclusive (repeat)", "warn"),
    tile(counts[["Not risky"]], "Not risky", "good"),
    "</div>",
    "<h2>Result per sample</h2>",
    table_html(c("Sample", if (multi_run) "Run", "Result", "Why",
                 "Highest beta", "Reference Ct"), verdict_rows),
    "<h2>Methylation level (beta) per gene</h2>",
    '<p class="muted">0 = unmethylated, 1 = as methylated as the positive ',
    "control. Bold red = called methylated; yellow = needs review; ",
    "grey = invalid.</p>",
    beta_table(x, multi_run),
    "<h2>Control checks</h2>",
    table_html(c(if (multi_run) "Run", "Gene", "No-template control",
                 "Positive control", "Positive control Ct"), ctrl_rows),
    if (nrow(review)) paste0(
      "<h2>Wells to check by eye</h2>",
      table_html(c(if (multi_run) "Run", "Well", "Sample", "Gene", "Ct",
                   "Why"), review_rows)),
    "<h2>Settings used</h2>",
    table_html(c("Setting", "Value"),
               vapply(seq_len(nrow(settings)), function(i) {
                 tr(c(esc(settings$Setting[i]), esc(settings$Value[i])))
               }, character(1))),
    "</main></body></html>"
  )
}

beta_table <- function(x, multi_run) {
  r <- x$results[x$results$role == "sample", ]
  if (length(x$settings$panel)) r <- r[r$target %in% x$settings$panel, ]
  ids <- unique(r[, c("run", "sample")])
  genes <- sort(unique(r$target))
  rows <- vapply(seq_len(nrow(ids)), function(i) {
    cells <- vapply(genes, function(g) {
      m <- r[r$run == ids$run[i] & r$sample == ids$sample[i] & r$target == g, ]
      if (!nrow(m)) return('<td class="na"></td>')
      call <- as.character(m$call[1])
      b <- m$beta[1]
      cls <- switch(call, Methylated = "meth", Review = "review",
                    Invalid = "invalid", "")
      shade <- if (!is.na(b) && call == "Methylated")
        sprintf(' style="background:rgba(192,57,43,%.2f)"', 0.15 + 0.6 * b) else ""
      label <- if (call == "Invalid") "invalid" else
        if (call == "Review") paste0(fmt_num(b, 3), " ?") else fmt_num(b, 3)
      sprintf('<td class="%s"%s title="%s">%s</td>', cls, shade, esc(call),
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
    "td.meth{font-weight:700;color:#7b1d13}td.review{background:#fdf1c4}",
    "td.invalid{background:#e0e0e0;color:#555}td.na{background:#fafafa}",
    ".badge{display:inline-block;padding:2px 8px;border-radius:10px;",
    "font-weight:600;white-space:nowrap}",
    ".badge.bad{background:#fbe3e0;color:#9b2c1f}",
    ".badge.warn{background:#fdf1c4;color:#7a5b00}",
    ".badge.good{background:#dff3e7;color:#1e6b41}",
    "@media print{body{background:#fff}.tile{border:1px solid #ccc}}"
  )
}

verdict_class <- function(v) {
  switch(v, "Potentially cancer" = "bad", "Inconclusive" = "warn",
         "Not risky" = "good", "")
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

fmt_setting <- function(v) {
  nm <- names(v)
  if (is.null(nm)) return(paste(format(v), collapse = ", "))
  paste(ifelse(nzchar(nm), paste0(nm, " = ", v), paste0("default ", v)),
        collapse = ", ")
}
