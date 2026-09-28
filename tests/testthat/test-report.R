res <- function(...) {
  analyze_qmsp(read_eds(write_fixture(), name = "run1"), positive = "^PC$", ...)
}
verdict_of <- function(r, s) as.character(r$report$verdict[r$report$sample == s])

test_that("beta is PMR / 100, capped at 1, and 0 when unmethylated", {
  r <- res()$results
  expect_equal(r$beta[r$sample == "S1"], 2^-4 / 2^-1)
  expect_equal(r$beta[r$sample == "PC"], 1)
  expect_equal(r$beta[r$sample == "S2"], 0)
  expect_true(is.na(r$beta[r$sample == "S3"]))
})

test_that("beta without a reference gene is relative to the positive control Ct", {
  r <- res(reference = NULL)$results
  s1 <- r[r$sample == "S1" & r$target == "GENE1", ]
  expect_equal(s1$beta, 2^-(30 - 28))
  expect_match(s1$beta_method, "no reference")
})

test_that("beta is information only: it does not change the call", {
  r <- res()$results
  s1 <- r[r$sample == "S1", ]
  expect_lt(s1$beta, 0.2)
  expect_equal(as.character(s1$call), "Methylated")
  expect_error(res(beta_cutoff = 0.2), "unused argument")
})

test_that("cutoffs can be set per gene", {
  r <- res(ct_cutoff = c(40, GENE1 = 31))
  expect_equal(as.character(r$results$call[r$results$sample == "S1"]),
               "Unmethylated")
  expect_error(res(ct_cutoff = c(GENE1 = 29)), "unnamed default")
})

test_that("every sample gets a verdict", {
  r <- res()
  expect_equal(verdict_of(r, "S1"), "Potential cancer")
  expect_equal(verdict_of(r, "S2"), "Low risk")
  expect_equal(verdict_of(r, "S3"), "Not determined")   # reference failed
  expect_equal(verdict_of(r, "S4"), "Low risk")         # drift curve
  expect_false(any(r$report$sample %in% c("PC", "NTC")))

  r2 <- res(min_methylated_genes = 2)
  expect_equal(verdict_of(r2, "S1"), "Low risk")
  expect_equal(verdict_of(res(panel = "OTHER"), "S1"), character())
})

test_that("the HTML report is written", {
  f <- write_report(res(), tempfile(fileext = ".html"))
  html <- paste(readLines(f), collapse = "\n")
  expect_match(html, "Potential cancer")
  expect_match(html, "10,000")
  expect_match(html, "Ct per gene")
  expect_match(html, "not a diagnosis")
  expect_match(html, "<table")
})
