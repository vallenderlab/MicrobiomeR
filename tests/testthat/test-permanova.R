library(MicrobiomeR)
library(testthat)

context("Testing of permanova")

# Use existing data for test.
data <- as_MicrobiomeR_format(raw_silva_2, format = "analyzed_format")

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


# Remove file created by test
if (file.exists("Rplots.pdf")) {
  file.remove("Rplots.pdf")
}
