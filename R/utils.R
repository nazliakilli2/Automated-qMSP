`%||%` <- function(a, b) if (is.null(a) || !length(a) || identical(a, "")) b else a

`%|na|%` <- function(a, b) if (is.na(a)) b else a

xml_text1 <- function(doc, xpath) {
  if (is.null(doc)) return(NULL)
  node <- xml2::xml_find_first(doc, xpath)
  if (inherits(node, "xml_missing")) return(NULL)
  xml2::xml_text(node)
}

ms_to_time <- function(ms) {
  ms <- suppressWarnings(as.numeric(ms))
  if (!length(ms) || is.na(ms) || ms <= 0) return(as.POSIXct(NA))
  as.POSIXct(ms / 1000, origin = "1970-01-01", tz = "UTC")
}

safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

safe_sd <- function(x) if (sum(!is.na(x)) < 2) NA_real_ else stats::sd(x, na.rm = TRUE)

# Case-insensitive match of names against a regular expression; NULL/"" = none.
matches <- function(x, pattern) {
  if (is.null(pattern) || !nzchar(pattern)) return(rep(FALSE, length(x)))
  grepl(pattern, x, ignore.case = TRUE, perl = TRUE)
}

join_flags <- function(...) {
  parts <- list(...)
  out <- do.call(paste, c(parts, sep = "; "))
  out <- gsub("(; )+", "; ", out)
  out <- gsub("^; |; $", "", out)
  out
}
