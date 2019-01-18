context("test-correlation-plot")


# Use existing data for test.
data <- as_MicrobiomeR_format(raw_silva_2, format = "analyzed_format")


test_that("basic correlation plot works", {
  expect_true(!is.null(get_correlation_plots(obj = data, primary_ranks = c("Class", "Order"))))
})
