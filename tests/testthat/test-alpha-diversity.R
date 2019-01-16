library(MicrobiomeR)
library(testthat)

context("Testing of alpha diversity measures and plot")

# Use existing data for test.
data <- as_MicrobiomeR_format(raw_silva, format = "analyzed_format")


test_that("alpha diversity measures all exist", {
  expect_equal(length(get_alpha_diversity_measures(data)), 10)
  expect_true(!is.null(get_alpha_diversity_measures(obj = data)$Shannon))
  expect_true(!is.null(get_alpha_diversity_measures(obj = data)$InverseSimpson))
  expect_true(!is.null(get_alpha_diversity_measures(obj = data)$GiniSimpson))
  expect_true(!is.null(get_alpha_diversity_measures(obj = data)$Fisher))
  expect_true(!is.null(get_alpha_diversity_measures(obj = data)$Coverage))
  expect_true(!is.null(get_alpha_diversity_measures(obj = data)$group.pairs))
})

test_that("base alpha diversity plot works", {
  expect_true(!is.null(alpha_diversity_plot(obj = data, measure = "shannon", select_otu_table = "otu_proportions")))
})

# Remove file created by test
if (file.exists("Rplots.pdf")) {
  file.remove("Rplots.pdf")
}
