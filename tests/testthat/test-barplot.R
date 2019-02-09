library(MicrobiomeR)
library(testthat)

context("Testing of stacked barplot")

# Use existing data for test.
data <- as_MicrobiomeR_format(raw_silva_2, format = "analyzed_format")
test_palette <- get_color_palette(color_no = 12)


test_that("basic stacked barplot works", {
  expect_true(!is.null(stacked_barplot(obj = data)))
})

test_that("creating multiple stacked barplots works", {
  expect_true(!is.null(stacked_barplots(obj = data)))
})

plots <- stacked_barplots(obj = data)
save_stacked_barplots(plots = plots, custom_path = "output/")

test_that("stacked barplots exist", {
  expect_true(file.exists("output/phylum_stacked_barplot.tiff"))
  expect_true(file.exists("output/class_stacked_barplot.tiff"))
  expect_true(file.exists("output/order_stacked_barplot.tiff"))
})

# Remove file created by test
if (file.exists("Rplots.pdf")) {
  file.remove("Rplots.pdf")
}

# Remove directory created by test
if (dir.exists("output")) {
  unlink("output", recursive = TRUE)
}
