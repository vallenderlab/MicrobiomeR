library(MicrobiomeR)
library(testthat)
library(metacoder)
library(taxa)


context("test-metacoder")

test_that("basic sample_id_filter works", {
  result <- sample_id_filter(
    obj = small_taxmap_fixture("raw_format"),
    .f_filter = ~sum(.),
    .f_condition = ~. >= 20,
    validated = FALSE
  )

  expect_is(result, "Taxmap")
})

test_that("basic taxon_id_filter works", {
  result <- taxon_id_filter(
    obj = small_taxmap_fixture("raw_format"),
    .f_filter = ~sum(.),
    .f_condition = ~. >= 2000,
    validated = FALSE
  )

  expect_is(result, "Taxmap")
})

test_that("basic otu_id_filter works", {
  result <- otu_id_filter(
    obj = small_taxmap_fixture("raw_format"),
    .f_transform = ~. / sum(.),
    .f_filter = ~mean(.),
    .f_condition = ~. > 0.00005
  )

  expect_is(result, "Taxmap")
})

test_that("bimodality coefficient helper returns finite numeric output", {
  bimodality_value <- microbiomer_bimodality_coefficient(c(1, 2, 3, 4, 5, 6))
  expect_type(bimodality_value, "double")
  expect_true(is.finite(bimodality_value))
})
test_that("basic agglomerate_taxmap works", {
  agglomerated_data <- agglomerate_taxmap(
    small_taxmap_fixture("raw_format"),
    "Phylum"
  )
  lower_ranks <- c("Class", "Order", "Family", "Genus", "Species")

  expect_false(any(unique(taxon_ranks(agglomerated_data)) %in% lower_ranks))
})

test_that("basic otu_proportion_filter works", {
  result <- otu_proportion_filter(
    obj = small_taxmap_fixture("raw_format"),
    otu_percentage = 0.00001
  )

  expect_is(result, "Taxmap")
})

test_that("basic otu_prevalence_filter works", {
  result <- otu_prevalence_filter(
    obj = small_taxmap_fixture("raw_format"),
    validated = FALSE
  )

  expect_is(result, "Taxmap")
})

test_that("basic taxa_prevalence_filter works", {
  result <- taxa_prevalence_filter(
    obj = small_taxmap_fixture("raw_format"),
    rank = "Class",
    validated = FALSE
  )

  expect_is(result, "Taxmap")
})

test_that("basic cov_filter works", {
  result <- cov_filter(
    obj = small_taxmap_fixture("raw_format"),
    coefficient_of_variation = 3,
    validated = FALSE
  )

  expect_is(result, "Taxmap")
})
