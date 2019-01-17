library(MicrobiomeR)
library(testthat)

context("Testing of ordination plots and methods")

# Use existing data for test.
data <- as_MicrobiomeR_format(raw_silva_2, format = "analyzed_format")


test_that("basic ordination plot works", {
  expect_true(!is.null(ordination_plot(obj = data)))
})

test_that("basic ordination data return works", {
  expect_true(!is.null(ordination_plot(obj = data, only_data = TRUE)))
})

# Remove file created by test
if (file.exists("Rplots.pdf")) {
  file.remove("Rplots.pdf")
}
