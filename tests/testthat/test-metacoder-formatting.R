library(MicrobiomeR)
library(testthat)

context("Testing of metacoder formatting functions")

test_that("which_format works", {
  expect_equal(which_format(create_taxmap(phyloseq_silva_2)), "phyloseq_format")
  expect_equal(which_format(raw_silva_2), "raw_format")
  expect_equal(which_format(as_basic_format(raw_silva_2)), "basic_format")
  expect_equal(which_format(as_analyzed_format(raw_silva_2)), "analyzed_format")
})

test_that("as_*_format works", {
  expect_true(is_raw_format(as_raw_format(create_taxmap(phyloseq_silva_2))))
})

test_that("is_*_format works", {
  expect_true(is_phyloseq_format(create_taxmap(phyloseq_silva_2)))
  expect_false(is_phyloseq_format(raw_silva_2))
  expect_true(is_raw_format(raw_silva_2))
  expect_false(is_raw_format(create_taxmap(phyloseq_silva_2)))
  expect_true(is_basic_format(as_basic_format(raw_silva_2)))
  expect_false(is_basic_format(raw_silva_2))
  expect_true(is_analyzed_format(as_analyzed_format(raw_silva_2)))
  expect_false(is_analyzed_format(raw_silva_2))
})

test_that("as_custom_format rejects unknown rename sources", {
  expect_error(
    as_custom_format(raw_silva_2, format = "raw_format", change_name_list = list(not_a_table = "otu_abundance")),
    "None of the supplied table names were found"
  )
})
