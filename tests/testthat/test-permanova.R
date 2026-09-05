library(MicrobiomeR)
library(testthat)

context("Testing of permanova")

# Use a deterministic subset so routine package checks stay fast.
data <- small_taxmap_fixture("analyzed_format")

test_that("default permanova works", {
  expect_equal(length(permanova(data)), 4)
  expect_true(!is.null(permanova(data)$coefficients))
  expect_true(!is.null(permanova(data)$top_coefficients))
})

test_that("permanova using unifrac methods works", {
  expect_warning(permanova(data, distance_method = "wunifrac"), "Coefficients were not able to be generated using this distance method.")
  expect_warning(permanova(data, distance_method = "unifrac"), "Coefficients were not able to be generated using this distance method.")
})

p <- permanova(data)

test_that("top coefficients barplot works", {
  p <- permanova(data)
  expect_true(!is.null(top_coefficients_barplot(top_coefficients = p$top_coefficients)))
})

test_that("permanova is reproducible and preserves RNG state", {
  set.seed(2026)
  seed_before <- .Random.seed

  permanova_a <- permanova(data, seed = 23)
  seed_after_first <- .Random.seed

  permanova_b <- permanova(data, seed = 23)
  seed_after_second <- .Random.seed

  expect_equal(permanova_a$permanova$aov.tab, permanova_b$permanova$aov.tab)
  expect_equal(permanova_a$anova, permanova_b$anova)
  expect_equal(permanova_a$coefficients, permanova_b$coefficients)
  expect_equal(permanova_a$top_coefficients, permanova_b$top_coefficients)
  expect_equal(seed_after_first, seed_before)
  expect_equal(seed_after_second, seed_before)
})


# Remove file created by test
if (file.exists("Rplots.pdf")) {
  file.remove("Rplots.pdf")
}
