library(MicrobiomeR)
library(testthat)

context("Testing of alpha diversity measures and plot")

# Create data for the test
data <- analyzed_silva


test_that("basic alpha diversity plot works", {
  expect_true(!is.null(alpha_diversity_plot(obj = data, measure = "shannon", select_otu_table = "otu_proportions", save = FALSE)))
})

# Remove file created by test
if (file.exists("Rplots.pdf")) {
  file.remove("Rplots.pdf")
}
