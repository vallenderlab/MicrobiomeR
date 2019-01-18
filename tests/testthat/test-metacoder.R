library(MicrobiomeR)
library(testthat)
library(metacoder)
library(taxa)


context("test-metacoder")

test_that("basic sample_id_filter works", {
  expect_is(sample_id_filter(obj = raw_silva_2, .f_filter = ~sum(.), .f_condition = ~.>= 20, validated = TRUE),
            "Taxmap")
})

test_that("basic taxon_id_filter works", {
  expect_is(taxon_id_filter(obj = raw_silva_2, .f_filter    = ~sum(.), .f_condition = ~.>= 2000, validated = TRUE),
            "Taxmap")
})

test_that("basic otu_id_filter works", {
  expect_is(otu_id_filter(obj = raw_silva_2, .f_transform = ~./sum(.), .f_filter = ~mean(.), .f_condition = ~.> 0.00005),
            "Taxmap")
})
