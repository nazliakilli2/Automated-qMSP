#' Read a QuantStudio / Applied Biosystems experiment file (.eds)
#'
#' An `.eds` file is a zip archive written by the QuantStudio Design &
#' Analysis software. It already contains the analysed results (Ct, amplification
#' status, Cq confidence) and the amplification curves for every well, so no
#' manual export is needed.
#'
#' @param path Path to an `.eds` file. A folder containing an unzipped `.eds`
#'   (with `apldbio/sds/analysis_result.txt` inside) or a bare
#'   `analysis_result.txt` also work.
#' @param name Optional run name. Defaults to the experiment name stored in the
#'   file, or the file name.
#'
#' @return A list of class `qmsp_run` with elements
#'   \describe{
#'     \item{wells}{one row per well and target: `run`, `well_index`, `well`
#'       (e.g. "C1"), `row`, `col`, `sample`, `target`, `task`, `ct` (`NA` when
#'       undetermined), `ct_raw`, `amp_status` (1 = amplified, 0 = inconclusive,
#'       -1 = no amplification) and `cq_conf`.}
#'     \item{curves}{long table of `run`, `well_index`, `target`, `cycle`, `rn`,
#'       `delta_rn`.}
#'     \item{thresholds}{threshold used by the instrument software per target.}
#'     \item{meta}{run name, file, start / end time, number of cycles.}
#'   }
#' @export
#' @examples
#' \dontrun{
#' run <- read_eds("my_run.eds")
#' head(run$wells)
#' }
read_eds <- function(path, name = NULL) {
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
  src <- eds_source(path)

  result_lines <- src$read("apldbio/sds/analysis_result.txt")
  if (is.null(result_lines)) {
    stop("No analysis results found in '", basename(path), "'. ",
         "Open and analyse the run in the QuantStudio software, save it, and ",
         "try again (template .edt files contain no results).", call. = FALSE)
  }

  experiment <- src$read_xml("apldbio/sds/experiment.xml")
  plate <- src$read_xml("apldbio/sds/plate_setup.xml")
  protocol <- src$read_xml("apldbio/sds/tcprotocol.xml")
  analysis <- src$read_xml("apldbio/sds/analysis_protocol.xml")

  run_name <- name %||% xml_text1(experiment, "/Experiment/Name") %||%
    tools::file_path_sans_ext(basename(path))
  n_cols <- as.integer(xml_text1(plate, "/Plate/Columns") %||% 12L)
  n_cycles <- max_cycles(protocol)

  parsed <- parse_analysis_result(result_lines)
  wells <- parsed$wells
  curves <- parsed$curves
  if (is.na(n_cycles) && nrow(curves)) n_cycles <- max(curves$cycle)

  wells$run <- rep(run_name, nrow(wells))
  wells$row <- LETTERS[wells$well_index %/% n_cols + 1L]
  wells$col <- wells$well_index %% n_cols + 1L
  wells$well <- paste0(wells$row, wells$col)
  # The software reports the last cycle (e.g. 50) as Ct for "Undetermined".
  undetermined <- is.na(wells$ct_raw) |
    (!is.na(n_cycles) & wells$ct_raw >= n_cycles)
  wells$ct <- ifelse(undetermined, NA_real_, wells$ct_raw)
  wells <- wells[, c("run", "well_index", "well", "row", "col", "sample",
                     "target", "task", "ct", "ct_raw", "amp_status", "cq_conf")]

  curves$run <- rep(run_name, nrow(curves))
  curves <- curves[, c("run", "well_index", "target", "cycle", "rn", "delta_rn")]

  thresholds <- detector_thresholds(analysis)
  if (nrow(thresholds)) thresholds$run <- run_name

  start <- xml_text1(experiment, "/Experiment/RunStartTime")
  end <- xml_text1(experiment, "/Experiment/RunEndTime")

  structure(
    list(
      wells = wells,
      curves = curves,
      thresholds = thresholds,
      meta = list(
        run = run_name,
        file = basename(path),
        n_cycles = n_cycles,
        run_start = ms_to_time(start),
        run_end = ms_to_time(end)
      )
    ),
    class = "qmsp_run"
  )
}

#' Read several .eds files into one dataset
#'
#' @param paths Character vector of `.eds` file paths (or a single folder, in
#'   which case every `.eds` file inside it is read).
#' @param names Optional run names, one per file.
#' @return A `qmsp_run` object with the wells, curves and thresholds of all
#'   runs stacked; `meta` becomes a data frame with one row per run.
#' @export
read_eds_files <- function(paths, names = NULL) {
  if (length(paths) == 1L && dir.exists(paths)) {
    paths <- list.files(paths, pattern = "\\.eds$", full.names = TRUE,
                        ignore.case = TRUE)
  }
  if (!length(paths)) stop("No .eds files given.", call. = FALSE)
  runs <- lapply(seq_along(paths), function(i) {
    read_eds(paths[[i]], name = if (!is.null(names)) names[[i]])
  })
  run_names <- vapply(runs, function(r) r$meta$run, character(1))
  if (anyDuplicated(run_names)) {
    dup <- duplicated(run_names) | duplicated(run_names, fromLast = TRUE)
    run_names[dup] <- make.unique(run_names[dup], sep = " #")
    runs <- Map(rename_run, runs, run_names)
  }
  combine_runs(runs)
}

#' @export
print.qmsp_run <- function(x, ...) {
  runs <- unique(x$wells$run)
  cat("<qmsp_run> ", length(runs), " run(s), ", nrow(x$wells), " wells\n",
      sep = "")
  for (r in runs) {
    w <- x$wells[x$wells$run == r, ]
    cat("  ", r, "\n    samples: ", paste(unique(w$sample), collapse = ", "),
        "\n    targets: ", paste(unique(w$target), collapse = ", "), "\n",
        sep = "")
  }
  invisible(x)
}

# ---- internals --------------------------------------------------------------

# Abstracts over "zip file", "unzipped folder" and "bare analysis_result.txt".
eds_source <- function(path) {
  if (dir.exists(path)) {
    get <- function(member) {
      f <- file.path(path, member)
      if (file.exists(f)) f else NULL
    }
    read <- function(member) {
      f <- get(member)
      if (is.null(f)) NULL else readLines(f, warn = FALSE, encoding = "UTF-8")
    }
  } else if (grepl("\\.txt$", path, ignore.case = TRUE)) {
    read <- function(member) {
      if (basename(member) == "analysis_result.txt") {
        readLines(path, warn = FALSE, encoding = "UTF-8")
      }
    }
  } else {
    members <- tryCatch(utils::unzip(path, list = TRUE)$Name,
                        error = function(e) {
      stop("'", basename(path), "' is not a valid .eds (zip) file.",
           call. = FALSE)
    })
    read <- function(member) {
      if (!member %in% members) return(NULL)
      con <- unz(path, member, encoding = "UTF-8")
      on.exit(close(con))
      readLines(con, warn = FALSE)
    }
  }
  read_xml <- function(member) {
    txt <- read(member)
    if (is.null(txt)) return(NULL)
    tryCatch(xml2::read_xml(paste(txt, collapse = "\n")),
             error = function(e) NULL)
  }
  list(read = read, read_xml = read_xml)
}

# analysis_result.txt layout (tab separated):
#   Session Name
#   Well  Sample Name  Detector  Task  Ct  ...  Amp Status  Cq Conf
#   24    NK-dH20      TAC1      Target 50.0 ...
#   Rn values        v1 v2 ...
#   Delta Rn values  v1 v2 ...
parse_analysis_result <- function(lines) {
  fields <- strsplit(lines, "\t", fixed = TRUE)
  first <- vapply(fields, function(f) if (length(f)) f[[1]] else "",
                  character(1))

  header_i <- which(first == "Well")[1]
  if (is.na(header_i)) stop("Unrecognised analysis_result.txt format.",
                            call. = FALSE)
  header <- trimws(fields[[header_i]])
  col_of <- function(nm) match(nm, header)
  i_sample <- col_of("Sample Name")
  i_target <- col_of("Detector") %|na|% col_of("Target Name")
  i_task <- col_of("Task")
  i_ct <- col_of("Ct") %|na|% col_of("CT") %|na|% col_of("Cq")
  i_amp <- col_of("Amp Status")
  i_conf <- col_of("Cq Conf")

  well_rows <- which(grepl("^[0-9]+$", first) & seq_along(first) > header_i)
  pick <- function(f, i) if (is.na(i) || i > length(f)) NA_character_ else f[[i]]
  wf <- fields[well_rows]
  wells <- data.frame(
    well_index = as.integer(first[well_rows]),
    sample = trimws(vapply(wf, pick, "", i = i_sample)),
    target = trimws(vapply(wf, pick, "", i = i_target)),
    task = trimws(vapply(wf, pick, "", i = i_task)),
    ct_raw = suppressWarnings(as.numeric(vapply(wf, pick, "", i = i_ct))),
    amp_status = suppressWarnings(as.integer(vapply(wf, pick, "", i = i_amp))),
    cq_conf = suppressWarnings(as.numeric(vapply(wf, pick, "", i = i_conf))),
    stringsAsFactors = FALSE
  )

  # Curve lines belong to the closest well row above them.
  owner <- findInterval(seq_along(first), well_rows)
  curve_block <- function(label) {
    idx <- which(first == label & owner > 0)
    out <- lapply(idx, function(i) {
      v <- suppressWarnings(as.numeric(fields[[i]][-1]))
      v <- v[!is.na(v)]
      data.frame(key = owner[i], cycle = seq_along(v), value = v)
    })
    do.call(rbind, c(out, list(data.frame(key = integer(), cycle = integer(),
                                          value = numeric()))))
  }
  rn <- curve_block("Rn values")
  drn <- curve_block("Delta Rn values")
  names(rn)[3] <- "rn"
  names(drn)[3] <- "delta_rn"
  curves <- merge(rn, drn, by = c("key", "cycle"), all = TRUE)
  curves$well_index <- wells$well_index[curves$key]
  curves$target <- wells$target[curves$key]
  curves <- curves[order(curves$key, curves$cycle),
                   c("well_index", "target", "cycle", "rn", "delta_rn")]
  rownames(curves) <- NULL

  list(wells = wells, curves = curves)
}

max_cycles <- function(protocol) {
  if (is.null(protocol)) return(NA_integer_)
  reps <- as.integer(xml2::xml_text(
    xml2::xml_find_all(protocol, ".//NumOfRepetitions")))
  if (!length(reps) || all(is.na(reps))) NA_integer_ else max(reps, na.rm = TRUE)
}

detector_thresholds <- function(analysis) {
  empty <- data.frame(target = character(), threshold = numeric(),
                      auto_threshold = logical(), stringsAsFactors = FALSE)
  if (is.null(analysis)) return(empty)
  blocks <- xml2::xml_find_all(
    analysis,
    "//JaxbAnalysisSettings[contains(Type, 'IDetectorSettings')]"
  )
  rows <- lapply(blocks, function(b) {
    settings <- xml2::xml_find_all(b, "./JaxbSettingValue")
    keys <- xml2::xml_text(xml2::xml_find_first(settings, "./Name"))
    vals <- xml2::xml_text(xml2::xml_find_first(settings, "./JaxbValueItem/*"))
    s <- stats::setNames(vals, keys)
    if (is.na(s["ObjectName"])) return(NULL)
    data.frame(target = unname(s["ObjectName"]),
               threshold = as.numeric(s["Threshold"]),
               auto_threshold = identical(unname(s["AutoThreshold"]), "true"),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, c(rows, list(empty)))
  rownames(out) <- NULL
  out
}

rename_run <- function(run, new) {
  run$wells$run <- rep(new, nrow(run$wells))
  run$curves$run <- rep(new, nrow(run$curves))
  if (nrow(run$thresholds)) run$thresholds$run <- new
  run$meta$run <- new
  run
}

combine_runs <- function(runs) {
  meta <- do.call(rbind, lapply(runs, function(r) {
    data.frame(run = r$meta$run, file = r$meta$file,
               n_cycles = r$meta$n_cycles, run_start = r$meta$run_start,
               run_end = r$meta$run_end, stringsAsFactors = FALSE)
  }))
  structure(
    list(
      wells = do.call(rbind, lapply(runs, `[[`, "wells")),
      curves = do.call(rbind, lapply(runs, `[[`, "curves")),
      thresholds = do.call(rbind, lapply(runs, function(r) {
        t <- r$thresholds
        if (!nrow(t)) t$run <- character()
        t
      })),
      meta = meta
    ),
    class = "qmsp_run"
  )
}
