#' Export results to Excel (or CSV files)
#'
#' Writes one sheet per table: the verdict per sample, wide tables of calls
#' and Ct values, the full per-gene results (with delta Ct, PMR and beta for
#' information), the control checks, the fluorescence thresholds, every well, and
#' the settings used. If the
#' `writexl` package is not installed, CSV files are written instead.
#'
#' @param x A `qmsp_result` from [analyze_qmsp()].
#' @param path Output `.xlsx` path. For CSV output, a folder.
#' @return The path(s) written, invisibly.
#' @export
export_results <- function(x, path = "qmsp_results.xlsx") {
  if (!inherits(x, "qmsp_result")) {
    stop("`x` must come from analyze_qmsp().", call. = FALSE)
  }
  wells <- x$wells
  report <- x$report
  report$verdict <- as.character(report$verdict)
  sheets <- list(
    Report = report,
    Calls = results_wide(x, "call"),
    Ct = results_wide(x, "ct"),
    Results = within_results(x$results),
    Controls = x$controls,
    Thresholds = x$thresholds,
    Wells = wells[, setdiff(names(wells), c("row", "col"))],
    Settings = data.frame(
      setting = names(x$settings),
      value = vapply(x$settings, function(v) {
        if (is.null(v)) "" else fmt_setting(v)
      }, character(1)),
      stringsAsFactors = FALSE
    )
  )
  if (requireNamespace("writexl", quietly = TRUE) &&
      grepl("\\.xlsx$", path, ignore.case = TRUE)) {
    writexl::write_xlsx(sheets, path)
    return(invisible(path))
  }
  dir <- if (grepl("\\.xlsx$", path, ignore.case = TRUE)) {
    sub("\\.xlsx$", "", path, ignore.case = TRUE)
  } else path
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  files <- file.path(dir, paste0(names(sheets), ".csv"))
  Map(function(d, f) utils::write.csv(d, f, row.names = FALSE), sheets, files)
  message("writexl not installed: wrote CSV files to ", dir)
  invisible(files)
}

within_results <- function(r) {
  r$call <- as.character(r$call)
  num <- vapply(r, is.numeric, logical(1))
  r[num] <- lapply(r[num], function(v) signif(v, 5))
  r
}
