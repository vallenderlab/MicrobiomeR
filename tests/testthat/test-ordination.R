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

plots <- ordination_plots(obj = data)
save_ordination_plots(ord_plots = plots, custom_path = "output/")

test_that("ordination plots exist", {
  expect_true(file.exists("output/pcoa_wunifrac_ordination.tiff"))
  expect_true(file.exists("output/pcoa_unifrac_ordination.tiff"))
  expect_true(file.exists("output/pcoa_bray_ordination.tiff"))
  expect_true(file.exists("output/nmds_wunifrac_ordination.tiff"))
  expect_true(file.exists("output/nmds_unifrac_ordination.tiff"))
  expect_true(file.exists("output/nmds_bray_ordination.tiff"))
})

# Remove file created by test
if (file.exists("Rplots.pdf")) {
  file.remove("Rplots.pdf")
}

# Remove directory created by test
if (dir.exists("output")) {
  unlink("output", recursive = TRUE)
}
