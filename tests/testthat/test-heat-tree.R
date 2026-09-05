library(MicrobiomeR)
library(testthat)

context("test-heat-tree")

# Build the reduced fixture once, then clone it within each test because Taxmap
# objects use reference semantics.
heat_tree_data <- small_taxmap_fixture("analyzed_format")


test_that("basic heat-tree plot works", {
  # Arrange
  input_data <- heat_tree_data$clone()

  # Act
  result <- heat_tree_plots(obj = input_data, rank_list = c("Phylum", "Class"))

  # Assert
  expect_true(!is.null(result))
})

test_that("heat tree node sizes count retained descendant OTUs", {
  # Arrange
  source_data <- heat_tree_data$clone()
  filtered_data <- MicrobiomeR:::filter_taxa_compat(
    obj = source_data$clone(),
    subset = n_supertaxa < pkg.private$rank_index[["Phylum"]],
    supertaxa = TRUE,
    reassign_obs = FALSE
  )

  # Act
  otu_counts <- MicrobiomeR:::get_heat_tree_otu_counts(
    filtered_obj = filtered_data,
    source_obj = source_data
  )

  # Assert
  retained_taxon_ids <- filtered_data$taxon_ids()
  source_taxon_ids <- source_data$taxon_ids()
  expected_counts <- source_data$n_obs("otu_abundance")[
    match(retained_taxon_ids, source_taxon_ids)
  ]

  expect_equal(unname(otu_counts), unname(expected_counts))
  expect_true(max(otu_counts) > 0)
})

test_that("heat tree parameters label descendant OTU counts", {
  # Arrange
  input_data <- heat_tree_data$clone()
  treatment_counts <- c(2L, 3L)

  # Act
  node_size_labels <- vapply(treatment_counts, function(treatment_no) {
    parameters <- suppressWarnings(
      heat_tree_parameters(
        obj = input_data,
        title = "Test heat tree",
        treatment_no = treatment_no
      )
    )
    parameters$node_size_axis_label
  }, character(1))

  # Assert
  expected_label <- "Number of retained descendant OTUs"
  expect_identical(unname(node_size_labels), rep(expected_label, 2L))
})

test_that("heat_tree_plots preserves its input Taxmap", {
  # Arrange
  input_data <- heat_tree_data$clone()
  taxon_ids_before <- input_data$taxon_ids()
  otu_abundance_before <- input_data$data$otu_abundance

  # Act
  heat_tree_plots(obj = input_data, rank_list = c("Phylum", "Class"))

  # Assert
  expect_identical(input_data$taxon_ids(), taxon_ids_before)
  expect_identical(input_data$data$otu_abundance, otu_abundance_before)
})

test_that("heat_tree_plots filters each rank independently", {
  # Arrange
  input_data <- heat_tree_data$clone()

  # Act
  result <- heat_tree_plots(
    obj = input_data,
    rank_list = c("Phylum", "Class", "Order")
  )

  # Assert
  rank_taxon_counts <- vapply(result$taxmaps, function(obj) {
    length(obj$taxon_ids())
  }, integer(1))

  expect_gt(rank_taxon_counts[["Class"]], rank_taxon_counts[["Phylum"]])
  expect_gt(rank_taxon_counts[["Order"]], rank_taxon_counts[["Class"]])
})

test_that("heat_tree_plots retains source OTU counts for every rank", {
  # Arrange
  input_data <- heat_tree_data$clone()

  # Act
  result <- heat_tree_plots(
    obj = input_data,
    rank_list = c("Phylum", "Class", "Order")
  )
  otu_counts_by_rank <- lapply(result$taxmaps, function(filtered_obj) {
    MicrobiomeR:::get_heat_tree_otu_counts(
      filtered_obj = filtered_obj,
      source_obj = result$metacoder_object
    )
  })

  # Assert
  all_nodes_have_counts <- vapply(otu_counts_by_rank, function(otu_counts) {
    all(otu_counts > 0)
  }, logical(1))

  expect_true(all(all_nodes_have_counts))
})

test_that("heat_tree_plots is reproducible and preserves RNG state", {
  # Arrange
  input_data <- heat_tree_data$clone()
  set.seed(2026)
  seed_before <- .Random.seed

  # Act
  heat_tree_a <- heat_tree_plots(
    obj = input_data,
    rank_list = "Phylum",
    seed = 17
  )$heat_trees$Phylum
  seed_after_first <- .Random.seed

  heat_tree_b <- heat_tree_plots(
    obj = input_data,
    rank_list = "Phylum",
    seed = 17
  )$heat_trees$Phylum
  seed_after_second <- .Random.seed

  built_a <- ggplot2::ggplot_build(heat_tree_a)$data
  built_b <- ggplot2::ggplot_build(heat_tree_b)$data

  # Assert
  expect_equal(built_a, built_b)
  expect_equal(seed_after_first, seed_before)
  expect_equal(seed_after_second, seed_before)
})
