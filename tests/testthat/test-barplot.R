library(MicrobiomeR)
library(testthat)

context("Testing of stacked barplot")

# Use existing data for test.
data <- as_MicrobiomeR_format(raw_silva_2, format = "analyzed_format")
test_palette <- get_color_palette(color_no = 12)


test_that("basic stacked barplot works", {
  expect_true(!is.null(stacked_barplot(obj = data, palette_values = test_palette)))
})

test_that("creating multiple stacked barplots works", {
  expect_true(!is.null(stacked_barplots(obj = data)))
})

# Remove file created by test
if (file.exists("Rplots.pdf")) {
  file.remove("Rplots.pdf")
}
