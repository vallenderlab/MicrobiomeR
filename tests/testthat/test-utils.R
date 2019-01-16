library(MicrobiomeR)
library(testthat)

context("Testing of utils.R")

# Use existing data for test.
#data <- as_MicrobiomeR_format(raw_silva, format = "analyzed_format")

test_that("get_output_dir function works", {
  expect_true(dir.exists(get_output_dir(start_path="output", experiment="test", mkdir=TRUE)))
  expect_equal(print(get_output_dir(start_path="output", experiment="test", mkdir=FALSE)), paste0(getwd(), "/output/test/", format(Sys.time(), "%Y-%m-%d_%s")))
})


# Remove direcotry created by test
if (dir.exists("output/test")) {
  unlink("output", recursive = TRUE)
}
