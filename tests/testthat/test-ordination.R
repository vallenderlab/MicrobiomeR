library(MicrobiomeR)
library(testthat)

context("Testing of ordination plots and methods")

# Use existing data for test.
data <- as_MicrobiomeR_format(raw_silva, format = "analyzed_format")


test_that("basic ordination plot works", {
  expect_true(!is.null(ordination_plot(obj = data)))
})

# Remove file created by test
if (file.exists("Rplots.pdf")) {
  file.remove("Rplots.pdf")
}
