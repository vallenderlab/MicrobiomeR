library(MicrobiomeR)
library(testthat)

context("Testing of alpha diversity measures and plot")

# Use existing data for test.
data <- as_MicrobiomeR_format(raw_silva_2, format = "analyzed_format")


test_that("alpha diversity measures all exist", {
  expect_equal(length(alpha_diversity_measures(data)), 10)
  expect_true(!is.null(alpha_diversity_measures(obj = data)$Shannon))
  expect_true(!is.null(alpha_diversity_measures(obj = data)$InverseSimpson))
  expect_true(!is.null(alpha_diversity_measures(obj = data)$GiniSimpson))
  expect_true(!is.null(alpha_diversity_measures(obj = data)$Fisher))
  expect_true(!is.null(alpha_diversity_measures(obj = data)$Coverage))
  expect_true(!is.null(alpha_diversity_measures(obj = data)$group.pairs))
})

test_that("base alpha diversity plot works", {
  expect_true(!is.null(alpha_diversity_plot(obj = data, measure = "Shannon")))
})

test_that("creating multiple alpha diversity plots works", {
  expect_true(!is.null(alpha_diversity_plots(obj = data)))
})

plots <- alpha_diversity_plots(obj = data)
save_alpha_diversity_plots(alpha_div_plots = plots, custom_path = "output/")

test_that("alpha diversity plots exist", {
  expect_true(file.exists("output/shannon_alpha_diversity.tiff"))
  expect_true(file.exists("output/ginisimpson_alpha_diversity.tiff"))
  expect_true(file.exists("output/inversesimpson_alpha_diversity.tiff"))
})

# Remove file created by test
if (file.exists("Rplots.pdf")) {
  file.remove("Rplots.pdf")
}

# Remove directory created by test
if (dir.exists("output")) {
  unlink("output", recursive = TRUE)
}
