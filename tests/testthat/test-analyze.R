res <- function(...) {
  analyze_qmsp(read_eds(write_fixture(), name = "run1"), positive = "^PC$", ...)
}
call_of <- function(r, s) as.character(r$results$call[r$results$sample == s])

test_that("samples get the expected calls", {
  r <- res()
  expect_s3_class(r, "qmsp_result")
  expect_equal(call_of(r, "S1"), "Methylated")
  expect_equal(call_of(r, "S2"), "Unmethylated")
  expect_equal(call_of(r, "S3"), "Invalid")     # reference failed
  expect_equal(call_of(r, "S4"), "Review")      # low confidence drift curve
  expect_equal(call_of(r, "PC"), "Methylated")
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
  expect_equal(call_of(res(ct_cutoff = 29), "S1"), "Unmethylated")
  r <- res(ref_ct_warn = 25.5)
  expect_equal(r$results$ref_status[r$results$sample == "S1"], "Low input")
  expect_equal(call_of(res(ref_ct_max = 25.5), "S1"), "Invalid")
})

test_that("weak curves are sent to review", {
  r <- res()
  w <- r$wells[r$wells$sample == "S4" & r$wells$target == "GENE1", ]
  expect_equal(w$result, "Review")
  expect_match(w$flags, "low Cq confidence")
  expect_match(w$flags, "weak curve")
})

test_that("a contaminated NTC is reported and puts sample positives in review", {
  w <- fixture_wells()
  w$ct[w$sample == "NTC" & w$target == "GENE1"] <- 33
  w$amp[w$sample == "NTC" & w$target == "GENE1"] <- 1
  w$conf[w$sample == "NTC" & w$target == "GENE1"] <- 0.9
  r <- analyze_qmsp(read_eds(write_fixture(wells = w)), positive = "^PC$")
  expect_equal(r$controls$ntc_status[r$controls$target == "GENE1"], "Fail")
  expect_equal(call_of(r, "S1"), "Review")
  expect_equal(call_of(r, "PC"), "Methylated")
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
