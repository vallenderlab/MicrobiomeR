library(MicrobiomeR)
library(testthat)
library(phyloseq)

context("Testing the phyloseq functions.")

# Get internal data from package
biom_file <- system.file("extdata", "silva_OTU.biom", package = "MicrobiomeR")
tree_file <- system.file("extdata", "silva.tre", package = "MicrobiomeR")
md_file <- system.file("extdata", "nephele_metadata2.txt", package = "MicrobiomeR")

test_that("data imported as phyloseq", {
  expect_is(
    create_phyloseq(biom_file = biom_file, tree_file = tree_file, metadata_file = md_file, parse_func = parse_taxonomy_silva_128, save_rooted_tree = TRUE, recursive_save = TRUE),
    "phyloseq")
  expect_true(ape::is.rooted(phy_tree(phy_obj)))
  expect_true(file.exists(system.file("extdata", "rooted_silva.tre", package = "MicrobiomeR")))
})
