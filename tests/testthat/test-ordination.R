library(MicrobiomeR)
library(testthat)

context("Testing of ordination plots and methods")

# Use existing data for test.
data <- as_MicrobiomeR_format(raw_silva_2, format = "analyzed_format")


test_that("basic ordination plot works", {
  expect_true(!is.null(ordination_plot(obj = data)))
  expect_true(!is.null(ordination_plot(obj = data, title = "Test Plot")))
})

test_that("basic ordination data return works", {
  expect_true(!is.null(ordination_plot(obj = data, only_data = TRUE)))
  expect_true(!is.null(ordination_plot(obj = data, only_data = TRUE, title = "Test Plot")))
})

test_that("creating multiple ordination plots works", {
  expect_true(!is.null(ordination_plots(obj = data)))
  expect_true(!is.null(ordination_plots(obj = data, methods = NULL, distances = c("wunifrac", "unifrac", "bray"))))
  expect_true(!is.null(ordination_plots(obj = data, methods = c("PCoA", "NMDS"), distances = NULL)))
})

# Remove file created by test
if (file.exists("Rplots.pdf")) {
  file.remove("Rplots.pdf")
}
