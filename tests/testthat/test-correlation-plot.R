context("test-correlation-plot")


# Use a deterministic subset so routine package checks stay fast.
data <- small_taxmap_fixture("analyzed_format")


test_that("basic correlation plot works", {
  expect_true(!is.null(correlation_plots(obj = data, primary_ranks = c("Class", "Order"))))
})
