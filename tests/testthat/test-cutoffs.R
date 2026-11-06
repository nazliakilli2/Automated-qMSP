# Cancers cross the instrument threshold at 30, 32, 35; controls at 34, 38
# and never. At delta Rn 10000 the fixture curves cross `shift` later.
shift <- 2 - 1.2 * log(2)

cutoff_fixture <- function() {
  s <- c("C1", "C2", "C3", "H1", "H2", "H3", "BAD", "PC", "NK")
  gene <- c(30, 32, 35, 34, 38, NA, 25, 28, NA)
  ref <- c(25, 25, 25, 25, 25, 25, NA, 25, NA)
  n <- length(s)
  wells <- data.frame(
    well = seq_len(2 * n) - 1,
    sample = rep(s, 2),
    target = rep(c("GENE1", "B ACTIN"), each = n),
    ct = c(gene, ref),
    amp = ifelse(is.na(c(gene, ref)), -1, 1),
    conf = ifelse(is.na(c(gene, ref)), 0, 0.95),
    weak = FALSE, stringsAsFactors = FALSE)
  read_eds(write_fixture(wells = wells), name = "train")
}
labels <- c(C1 = "cancer", C2 = "Cancer", C3 = "cancer", BAD = "cancer",
            H1 = "control", H2 = "healthy", H3 = "control",
            PC = "positive control", NK = "NK")

test_that("the best cutoff separates cancer from control", {
  r <- find_cutoffs(cutoff_fixture(), labels)
  expect_s3_class(r, "qmsp_cutoffs")
  c <- r$cutoffs
  expect_equal(c$target, "GENE1")
  expect_equal(c$n_cancer, 3)               # BAD has failed DNA: left out
  expect_equal(c$n_control, 3)
  expect_equal(c$cutoff, round(33 + shift, 2), tolerance = 0.01)
  expect_equal(c$sensitivity, 2 / 3)
  expect_equal(c$specificity, 1)
  expect_equal(c$auc, 8 / 9)
  expect_equal(c(c$tp, c$fn, c$fp, c$tn), c(2, 1, 0, 3))
  expect_match(c$note, "few samples")
  expect_equal(c$positive_ct, 28 + shift, tolerance = 0.01)
  expect_true(is.na(c$ntc_ct))
  expect_false(r$samples$used[r$samples$sample == "BAD"])
  expect_output(print(r), "qmsp_cutoffs")
})

test_that("methods favour specificity or sensitivity", {
  x <- cutoff_fixture()
  spec <- find_cutoffs(x, labels, method = "specificity", target = 0.9)$cutoffs
  expect_equal(spec$specificity, 1)
  sens <- find_cutoffs(x, labels, method = "sensitivity", target = 0.9)$cutoffs
  expect_equal(sens$sensitivity, 1)
  expect_equal(sens$specificity, 2 / 3)
  expect_equal(sens$cutoff, round(36.5 + shift, 2), tolerance = 0.01)
  hand <- find_cutoffs(x, labels, cutoffs = c(GENE1 = 37))$cutoffs
  expect_equal(hand$cutoff, 37)
  expect_equal(hand$source, "chosen")
  expect_equal(hand$suggested_cutoff, round(33 + shift, 2), tolerance = 0.01)
  expect_equal(c(hand$sensitivity, hand$specificity), c(1, 2 / 3))
  keep <- find_cutoffs(x, labels, exclude_failed_reference = FALSE)$cutoffs
  expect_equal(keep$n_cancer, 4)
})

test_that("labels can be a data frame; missing groups give no cutoff", {
  x <- cutoff_fixture()
  df <- data.frame(sample = names(labels), group = unname(labels))
  expect_equal(find_cutoffs(x, df)$cutoffs$cutoff,
               find_cutoffs(x, labels)$cutoffs$cutoff)
  none <- find_cutoffs(x, labels[!grepl("control|healthy", labels)])$cutoffs
  expect_true(is.na(none$cutoff))
  expect_match(none$note, "at least one cancer and one control")
})

test_that("cutoffs are saved, read back and used by analyze_qmsp", {
  x <- cutoff_fixture()
  r <- find_cutoffs(x, labels)
  f <- write_cutoffs(r, tempfile(fileext = ".csv"))
  cut <- read_cutoffs(f)
  expect_equal(unname(cut[1]), 40)
  expect_equal(unname(cut["GENE1"]), r$cutoffs$cutoff)
  res <- analyze_qmsp(x, positive = "^PC$", ct_cutoff = cut)
  call <- as.character(res$results$call[res$results$target == "GENE1"])
  names(call) <- res$results$sample[res$results$target == "GENE1"]
  expect_equal(unname(call[c("C2", "C3", "H1")]),
               c("Methylated", "Unmethylated", "Unmethylated"))

  p <- panel_performance(r)
  expect_equal(p$summary$sensitivity, 2 / 3)
  expect_equal(p$summary$specificity, 1)
  expect_s3_class(plot_cutoffs(r), "ggplot")
  out <- tempfile(fileext = ".xlsx")
  export_cutoffs(r, out, chosen = c(GENE1 = 37))
  expect_true(file.exists(out) || dir.exists(sub("\\.xlsx$", "", out)))
  html <- report_html(res)
  expect_match(html, "GENE1 [0-9.]+; other genes 40")
})
