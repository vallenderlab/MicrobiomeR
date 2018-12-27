library(MicrobiomeR)
library(testthat)

context("Testing of alpha diversity measures and plot")

# Use existing data for test.
data <- as_MicrobiomeR_format(raw_silva, format = "analyzed_format")


test_that("basic alpha diversity plot works", {
  expect_true(!is.null(alpha_diversity_plot(obj = data, measure = "shannon", select_otu_table = "otu_proportions", save = FALSE)))
})

# Remove file created by test
if (file.exists("Rplots.pdf")) {
  file.remove("Rplots.pdf")
}
