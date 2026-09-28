res <- function(...) {
  analyze_qmsp(read_eds(write_fixture(), name = "run1"), positive = "^PC$", ...)
}
call_of <- function(r, s) as.character(r$results$call[r$results$sample == s])

# Fixture curves (height 30000) cross delta Rn 10000 this many cycles after
# the instrument Ct.
shift <- 2 - 1.2 * log(2)

test_that("samples get the expected calls", {
  r <- res()
  expect_s3_class(r, "qmsp_result")
  expect_equal(call_of(r, "S1"), "Methylated")
  expect_equal(call_of(r, "S2"), "Unmethylated")
  expect_equal(call_of(r, "S3"), "Not determined")  # reference failed
  expect_equal(call_of(r, "S4"), "Unmethylated")    # drift never reaches 10000
  expect_equal(call_of(r, "PC"), "Methylated")
})

test_that("Ct is read at delta Rn 10000 and methylated means Ct <= 40", {
  r <- res()
  expect_equal(r$settings$threshold, 10000)
  expect_equal(r$settings$ct_cutoff, 40)
  w <- r$wells
  s1 <- w$sample == "S1" & w$target == "GENE1"
  expect_equal(w$ct[s1], 30 + shift, tolerance = 0.001)
  expect_equal(w$ct_instrument[s1], 30)
  expect_true(all(r$thresholds$threshold == 10000))
  expect_true(all(r$thresholds$source == "set"))
  s4 <- w$sample == "S4" & w$target == "GENE1"
  expect_true(is.na(w$ct[s4]))
  expect_match(w$flags[s4], "does not reach the threshold")

  # a sample crossing at Ct 40.5 is not methylated, at 39.5 it is
  wl <- fixture_wells()
  wl$ct[wl$sample == "S1" & wl$target == "GENE1"] <- 40.5 - shift
  wl$ct[wl$sample == "S5" & wl$target == "GENE1"] <- 39.5 - shift
  r2 <- analyze_qmsp(read_eds(write_fixture(wells = wl)), positive = "^PC$")
  expect_equal(call_of(r2, "S1"), "Unmethylated")
  expect_equal(call_of(r2, "S5"), "Methylated")

  inst <- res(threshold = NULL)
  expect_equal(inst$wells$ct[s1], 30)
})

test_that("delta Ct, ratio and PMR are computed against the positive control", {
  r <- res()$results
  s1 <- r[r$sample == "S1", ]
  pc <- r[r$sample == "PC", ]
  expect_equal(s1$delta_ct, 4)
  expect_equal(s1$ratio, 2^-4)
  expect_equal(pc$pmr, 100)
  expect_equal(s1$pmr, 100 * 2^-4 / 2^-1)
  expect_equal(r$ratio[r$sample == "S2"], 0)
})

test_that("Ct cutoff and reference cutoffs are applied", {
  expect_equal(call_of(res(ct_cutoff = 31), "S1"), "Unmethylated")
  # cutoff below the positive control: the negative call cannot be trusted
  expect_equal(call_of(res(ct_cutoff = 29), "S1"), "Not determined")
  r <- res(ref_ct_warn = 25.5)
  expect_equal(r$results$ref_status[r$results$sample == "S1"], "Low input")
  expect_equal(call_of(res(ref_ct_max = 25.5), "S1"), "Not determined")
})

test_that("a contaminated NTC makes methylated sample calls not determined", {
  w <- fixture_wells()
  w$ct[w$sample == "NTC" & w$target == "GENE1"] <- 33
  w$amp[w$sample == "NTC" & w$target == "GENE1"] <- 1
  w$conf[w$sample == "NTC" & w$target == "GENE1"] <- 0.9
  r <- analyze_qmsp(read_eds(write_fixture(wells = w)), positive = "^PC$")
  expect_equal(r$controls$ntc_status[r$controls$target == "GENE1"], "Fail")
  expect_equal(call_of(r, "S1"), "Not determined")
  expect_equal(call_of(r, "S2"), "Unmethylated")
  expect_equal(call_of(r, "PC"), "Methylated")
})

test_that("a failed positive control makes unmethylated calls not determined", {
  w <- fixture_wells()
  w$ct[w$sample == "PC" & w$target == "GENE1"] <- NA
  r <- analyze_qmsp(read_eds(write_fixture(wells = w)), positive = "^PC$")
  expect_equal(r$controls$positive_status[r$controls$target == "GENE1"], "Fail")
  expect_equal(call_of(r, "S2"), "Not determined")
  expect_equal(call_of(r, "S1"), "Methylated")
})

test_that("replicates that mostly disagree are not determined", {
  w <- fixture_wells()
  extra <- w[w$sample == "S5" & w$target == "GENE1", ][c(1, 1), ]
  extra$well <- c(38, 39)
  extra$ct <- NA
  r <- analyze_qmsp(read_eds(write_fixture(wells = rbind(w, extra))),
                    positive = "^PC$")
  expect_equal(call_of(r, "S5"), "Not determined")
  expect_match(r$results$notes[r$results$sample == "S5"], "1 of 3")
})

test_that("runs without a reference gene are called from gene Ct alone", {
  r <- res(reference = NULL)
  expect_true(all(r$results$ref_status[r$results$role == "sample"] == "Not run"))
  expect_true(all(is.na(r$results$delta_ct)))
  expect_equal(call_of(r, "S1")[r$results$target[r$results$sample == "S1"] == "GENE1"],
               "Methylated")
})

test_that("wide table, export and plots work", {
  r <- res()
  wide <- results_wide(r)
  expect_equal(wide$GENE1[wide$sample == "S1"], "Methylated")
  expect_false("PC" %in% wide$sample)

  out <- tempfile(fileext = ".xlsx")
  export_results(r, out)
  if (requireNamespace("writexl", quietly = TRUE)) {
    expect_true(file.exists(out))
  } else {
    expect_true(dir.exists(sub("\\.xlsx$", "", out)))
  }

  expect_s3_class(plot_methylation(r), "ggplot")
  expect_s3_class(plot_methylation(r, "pmr", controls = TRUE), "ggplot")
  expect_s3_class(plot_amplification(r), "ggplot")
  expect_s3_class(plot_plate(r, fill = "result"), "ggplot")
  expect_output(print(r), "qmsp_result")
})

test_that("Ct is recomputed from the curves at a threshold", {
  expect_equal(ct_at_threshold(1:5, c(0, 1, 10, 100, 1000), 10), 3)
  expect_equal(ct_at_threshold(1:5, c(0, 1, 10, 100, 1000), 31.6227766), 3.5,
               tolerance = 1e-6)
  expect_true(is.na(ct_at_threshold(1:5, c(0, 1, 10, 100, 50), 80)))
  # an early noise spike is ignored: the last crossing counts
  expect_equal(ct_at_threshold(1:6, c(0, 20, 1, 10, 100, 1000), 10), 4)

  thr_exact <- 30000 / (1 + exp(2 / 1.2))  # the fixture curves reach it at Ct
  r <- res(threshold = c(GENE1 = thr_exact))
  w <- r$wells
  s1 <- w$sample == "S1" & w$target == "GENE1"
  expect_equal(w$ct[s1], 30, tolerance = 0.2)
  expect_equal(w$ct_instrument[s1], 30)
  ref <- w$target == "B ACTIN" & w$sample == "S1"
  expect_equal(w$ct[ref], w$ct_instrument[ref])  # other genes untouched
  expect_equal(r$thresholds$source[r$thresholds$target == "GENE1"], "set")
  expect_false("B ACTIN" %in% r$thresholds$target)  # fixture has none
  expect_equal(r$thresholds$threshold[r$thresholds$target == "GENE1"], thr_exact)

  high <- res(threshold = c(GENE1 = 1e6))
  expect_true(all(is.na(high$wells$ct[high$wells$target == "GENE1"])))
  expect_equal(call_of(high, "S1")[1], "Not determined")  # PC fails too

  expect_s3_class(plot_amplification(r, ct_cutoff = c(38, GENE1 = 35)), "ggplot")
})
