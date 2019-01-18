library(MicrobiomeR)
library(testthat)

test_that("which_format works", {
  expect_equal(which_format(object_handler(phyloseq_silva_2)), "phyloseq_format")
  expect_equal(which_format(raw_silva_2), "raw_format")
  expect_equal(which_format(as_basic_format(raw_silva_2)), "basic_format")
  expect_equal(which_format(as_analyzed_format(raw_silva_2)), "analyzed_format")
})

test_that("is_*_format works", {
  expect_true(is_phyloseq_format(object_handler(phyloseq_silva_2)))
  expect_false(is_phyloseq_format(raw_silva_2))
  expect_true(is_raw_format(raw_silva_2))
  expect_false(is_raw_format(object_handler(phyloseq_silva_2)))
  expect_true(is_basic_format(as_basic_format(raw_silva_2)))
  expect_false(is_basic_format(raw_silva_2))
  expect_true(is_analyzed_format(as_analyzed_format(raw_silva_2)))
  expect_false(is_analyzed_format(raw_silva_2))
})

