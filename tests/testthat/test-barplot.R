library(MicrobiomeR)
library(testthat)

context("Testing of stacked barplot")

# Create data for the test
data <- analyzed_silva
test_palette <- get_color_palette(color_no = 12)


test_that("basic barplot works", {
  expect_true(!is.null(stacked_barplot(obj = data, palette_values = test_palette)))
})

# Remove file created by test
if (file.exists("Rplots.pdf")) {
  file.remove("Rplots.pdf")
}
