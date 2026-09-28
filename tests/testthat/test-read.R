test_that("analysis_result.txt is parsed into wells and curves", {
  x <- read_eds(write_fixture(), name = "run1")
  expect_s3_class(x, "qmsp_run")
  expect_equal(nrow(x$wells), 14)
  expect_equal(x$wells$well[1:3], c("A1", "A2", "A3"))
  expect_equal(x$wells$well[x$wells$well_index == 12][1], "B1")
  # "50.0" means undetermined
  expect_true(is.na(x$wells$ct[x$wells$sample == "S2" & x$wells$target == "GENE1"]))
  expect_equal(x$wells$ct_raw[x$wells$sample == "S2" & x$wells$target == "GENE1"], 50)
  expect_equal(nrow(x$curves), 14 * 50)
  expect_equal(unique(x$wells$run), "run1")
})

test_that("unzipped folders and zipped .eds files are read", {
  dir <- tempfile()
  dir.create(file.path(dir, "apldbio", "sds"), recursive = TRUE)
  write_fixture(file.path(dir, "apldbio", "sds", "analysis_result.txt"))
  writeLines(c("<Experiment><Name>My run</Name>",
               "<RunStartTime>1790615284603</RunStartTime></Experiment>"),
             file.path(dir, "apldbio", "sds", "experiment.xml"))
  writeLines("<Protocol><Stage><NumOfRepetitions>50</NumOfRepetitions></Stage></Protocol>",
             file.path(dir, "apldbio", "sds", "tcprotocol.xml"))
  x <- read_eds(dir)
  expect_equal(x$meta$run, "My run")
  expect_equal(x$meta$n_cycles, 50)

  skip_if(Sys.which("zip") == "", "zip not available")
  eds <- tempfile(fileext = ".eds")
  withr_dir <- getwd()
  on.exit(setwd(withr_dir))
  setwd(dir)
  utils::zip(eds, "apldbio", flags = "-rq")
  setwd(withr_dir)
  y <- read_eds(eds)
  expect_equal(y$wells, x$wells)
})

test_that("files without results give a helpful error", {
  f <- tempfile(fileext = ".eds")
  writeLines("not a zip", f)
  expect_error(read_eds(f), "not a valid")
})

test_that("several files are combined and duplicate names made unique", {
  f <- write_fixture()
  x <- read_eds_files(c(f, f), names = c("a", "a"))
  expect_equal(sort(unique(x$wells$run)), c("a", "a #1"))
  expect_equal(nrow(x$wells), 28)
})
