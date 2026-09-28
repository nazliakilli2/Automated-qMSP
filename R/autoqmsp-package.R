#' autoqmsp: automated qMSP analysis
#'
#' Typical workflow:
#' ```
#' runs <- read_eds_files(c("run1.eds", "run2.eds"))
#' res  <- analyze_qmsp(runs, reference = "B ACTIN")
#' res$results
#' plot_methylation(res)
#' export_results(res, "results.xlsx")
#' ```
#' Or start the point-and-click app with [run_app()].
#'
#' @keywords internal
#' @importFrom ggplot2 .data
"_PACKAGE"

#' Start the autoqmsp app
#'
#' Opens a browser app where you upload `.eds` files, adjust settings and
#' download an Excel report. Needs the `shiny` package.
#'
#' @param ... Passed to [shiny::runApp()].
#' @export
run_app <- function(...) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Install the 'shiny' package first: install.packages(\"shiny\")",
         call. = FALSE)
  }
  shiny::runApp(system.file("shiny", package = "autoqmsp"), ...)
}
