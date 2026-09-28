# Builds a small synthetic QuantStudio analysis_result.txt.
sigmoid <- function(ct, height = 30000, n = 50) {
  if (is.na(ct)) return(stats::rnorm(n, 0, 50))
  height / (1 + exp(-(seq_len(n) - ct - 2) / 1.2))
}
drift <- function(slope = 60, n = 50) slope * pmax(seq_len(n) - 15, 0)

fixture_wells <- function() {
  data.frame(
    well = c(0, 1, 2, 3, 12, 13, 14, 15, 24, 25, 26, 27, 36, 37),
    sample = c("PC", "PC", "S1", "S1", "S2", "S2", "S3", "S3",
               "NTC", "NTC", "S4", "S4", "S5", "S5"),
    target = c("GENE1", "B ACTIN", "GENE1", "B ACTIN", "GENE1", "B ACTIN",
               "GENE1", "B ACTIN", "GENE1", "B ACTIN", "GENE1", "B ACTIN",
               "GENE1", "B ACTIN"),
    ct = c(28, 27, 30, 26, NA, 25, 32, NA, NA, NA, 29, 28, 30, 27),
    amp = c(1, 1, 1, 1, -1, 1, 1, -1, -1, -1, 1, 1, 1, 1),
    conf = c(.98, .97, .95, .96, 0, .97, .9, 0, 0, 0, .2, .95, .9, .96),
    weak = c(F, F, F, F, F, F, F, F, F, F, T, F, F, F),
    stringsAsFactors = FALSE
  )
}

write_fixture <- function(path = tempfile(fileext = ".txt"), wells = fixture_wells()) {
  set.seed(1)
  lines <- c("Session Name\t",
             paste("Well", "Sample Name", "Detector", "Task", "Ct", "Avg Ct",
                   "Ct SD", "Delta Ct", "Qty", "Avg Qty", "Qty SD",
                   "Amp Status", "Cq Conf", sep = "\t"))
  for (i in seq_len(nrow(wells))) {
    w <- wells[i, ]
    ct <- if (is.na(w$ct)) "50.0" else format(w$ct)
    lines <- c(lines, paste(w$well, w$sample, w$target, "Target", ct, ct, "",
                            "", "", "", "", w$amp, w$conf, sep = "\t"))
    d <- if (w$weak) drift() else sigmoid(w$ct)
    lines <- c(lines,
               paste(c("Rn values", format(d + 50000)), collapse = "\t"),
               paste(c("Delta Rn values", format(d)), collapse = "\t"))
  }
  writeLines(lines, path)
  path
}
