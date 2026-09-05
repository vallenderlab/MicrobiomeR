library(MicrobiomeR)
library(testthat)

context("Testing of metacoder formatting functions")

test_that("which_format works", {
  expect_equal(which_format(small_taxmap_fixture("phyloseq_format")), "phyloseq_format")
  expect_equal(which_format(small_taxmap_fixture("raw_format")), "raw_format")
  expect_equal(which_format(small_taxmap_fixture("basic_format")), "basic_format")
  expect_equal(which_format(small_taxmap_fixture("analyzed_format")), "analyzed_format")
})

test_that("as_*_format works", {
  expect_true(is_raw_format(as_raw_format(small_taxmap_fixture("phyloseq_format"))))
})

test_that("is_*_format works", {
  expect_true(is_phyloseq_format(small_taxmap_fixture("phyloseq_format")))
  expect_false(is_phyloseq_format(small_taxmap_fixture("raw_format")))
  expect_true(is_raw_format(small_taxmap_fixture("raw_format")))
  expect_false(is_raw_format(small_taxmap_fixture("phyloseq_format")))
  expect_true(is_basic_format(small_taxmap_fixture("basic_format")))
  expect_false(is_basic_format(small_taxmap_fixture("raw_format")))
  expect_true(is_analyzed_format(small_taxmap_fixture("analyzed_format")))
  expect_false(is_analyzed_format(small_taxmap_fixture("raw_format")))
})

test_that("as_custom_format rejects unknown rename sources", {
  expect_error(
    as_custom_format(
      small_taxmap_fixture("raw_format"),
      format = "raw_format",
      change_name_list = list(not_a_table = "otu_abundance")
    ),
    "None of the supplied table names were found"
  )
})
