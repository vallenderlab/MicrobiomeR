library(MicrobiomeR)
library(testthat)

context("test-heat-tree")

# Use existing data for test.
data <- as_MicrobiomeR_format(raw_silva_2, format = "analyzed_format")


test_that("basic heat-tree plot works", {
  expect_true(!is.null(heat_tree_plots(obj = data, rank_list = c("Phylum", "Class"))))
})

test_that("heat tree node counts are preserved for retained taxa", {
  filtered_data <- MicrobiomeR:::filter_taxa_compat(
    obj = data,
    subset = n_supertaxa < pkg.private$rank_index[["Phylum"]],
    supertaxa = TRUE,
    reassign_obs = FALSE
  )
  otu_counts <- MicrobiomeR:::get_heat_tree_otu_counts(
    filtered_obj = filtered_data,
    source_obj = data
  )

  expect_true(max(otu_counts) > 0)
})

test_that("heat_tree_plots is reproducible and preserves RNG state", {
  set.seed(2026)
  seed_before <- .Random.seed

  heat_tree_a <- heat_tree_plots(
    obj = data,
    rank_list = "Phylum",
    seed = 17
  )$heat_trees$Phylum
  seed_after_first <- .Random.seed

  heat_tree_b <- heat_tree_plots(
    obj = data,
    rank_list = "Phylum",
    seed = 17
  )$heat_trees$Phylum
  seed_after_second <- .Random.seed

  built_a <- ggplot2::ggplot_build(heat_tree_a)$data
  built_b <- ggplot2::ggplot_build(heat_tree_b)$data

  expect_equal(built_a, built_b)
  expect_equal(seed_after_first, seed_before)
  expect_equal(seed_after_second, seed_before)
})
