context("test-heat-tree")

# Use existing data for test.
data <- as_MicrobiomeR_format(raw_silva_2, format = "analyzed_format")


test_that("basic heat-tree plot works", {
  expect_true(!is.null(get_heat_tree_plots(obj = data, rank_list = c("Phylum", "Class"))))
})
