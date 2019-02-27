library(MicrobiomeR)
library(testthat)
library(metacoder)
library(taxa)


context("test-metacoder")

test_that("basic sample_id_filter works", {
  expect_is(sample_id_filter(obj = raw_silva_2, .f_filter = ~sum(.), .f_condition = ~.>= 20, validated = FALSE),
            "Taxmap")
})

test_that("basic taxon_id_filter works", {
  expect_is(taxon_id_filter(obj = raw_silva_2, .f_filter    = ~sum(.), .f_condition = ~.>= 2000, validated = FALSE),
            "Taxmap")
})

test_that("basic otu_id_filter works", {
  expect_is(otu_id_filter(obj = raw_silva_2, .f_transform = ~./sum(.), .f_filter = ~mean(.), .f_condition = ~.> 0.00005),
            "Taxmap")
})
test_that("basic agglomerate_taxmap works", {
  expect_false(any(unique(taxon_ranks(agglomerate_taxmap(raw_silva_2, "Phylum"))) %in% c("Class", "Order", "Family", "Genus", "Species")))
})

test_that("basic otu_proportion_filter works", {
  expect_is(otu_proportion_filter(obj = raw_silva_2, otu_percentage = 0.00001),
            "Taxmap")
})

test_that("basic otu_prevalence_filter works", {
  expect_is(otu_prevalence_filter(obj = raw_silva_2, validated = FALSE),
            "Taxmap")
})

test_that("basic taxa_prevalence_filter works", {
  expect_is(taxa_prevalence_filter(obj = raw_silva_2, rank = "Class", validated = FALSE),
            "Taxmap")
})

test_that("basic cov_filter works", {
  expect_is(cov_filter(obj = raw_silva_2, coefficient_of_variation = 3, validated = FALSE),
            "Taxmap")
})

