library(MicrobiomeR)
library(testthat)

context("Testing of palettes")

test_that("viridis palette works", {
  expect_true(!is.null(viridis::viridis(n=25, option="magma")))
  expect_equal(length(viridis::viridis(n=25, option="magma")), 25)
})

test_that("scico palette works", {
  expect_true(!is.null(scico_palette(scico_palette="hawaii")(25)))
  expect_equal(length(scico_palette(scico_palette="hawaii")(25)), 25)
})
