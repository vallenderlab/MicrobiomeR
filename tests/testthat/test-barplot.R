library(MicrobiomeR)
library(testthat)

context("Testing of stacked barplot")

# Use a deterministic subset so routine package checks stay fast.
data <- small_taxmap_fixture("analyzed_format")
test_palette <- get_color_palette(color_no = 12)


test_that("basic stacked barplot works", {
  expect_true(!is.null(stacked_barplot(obj = data)))
})

test_that("stacked barplot displays sample labels", {
  plot <- stacked_barplot(obj = data)

  expect_s3_class(plot$theme$axis.text.x, "element_text")
  expect_equal(plot$theme$axis.text.x$angle, 90)
})

test_that("stacked barplot removes y-axis baseline padding", {
  plot <- stacked_barplot(obj = data)
  y_scale <- plot$scales$get_scales("y")

  expect_equal(y_scale$expand, ggplot2::expansion(mult = c(0, 0)))
})

test_that("creating multiple stacked barplots works", {
  expect_true(!is.null(stacked_barplots(obj = data)))
})

plots <- stacked_barplots(obj = data)
save_stacked_barplots(sb_plots = plots, custom_path = "output/")

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
