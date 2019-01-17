library(MicrobiomeR)
library(testthat)

context("Testing of utils.R")

# Use existing data for test.
basic_silva <- as_MicrobiomeR_format(raw_silva, format = "basic_format")

test_that("get_output_dir function works", {
  expect_true(!dir.exists(get_output_dir(start_path="output", experiment="test", mkdir=FALSE)))
  expect_true(dir.exists(get_output_dir(start_path="output", experiment="test", mkdir=TRUE)))
  expect_error(get_output_dir(start_path="output", experiment="test", mkdir=TRUE, overwrite=FALSE), paste0("The directory ", getwd(), "/output/test/", format(Sys.time(), "%Y-%m-%d_%s"), " already exists. And you don't want to overwrite the directory."))
  expect_equal(get_output_dir(start_path="output", experiment="test", mkdir=FALSE), paste0(getwd(), "/output/test/", format(Sys.time(), "%Y-%m-%d_%s")))
  expect_equal(get_output_dir(start_path="output", experiment="test", mkdir=FALSE, end_path = "cool"), paste0(getwd(), "/output/test/cool"))
  expect_equal(get_output_dir(mkdir=FALSE, plot_type = "cool"), paste0(getwd(), "/cool"))
  expect_equal(get_output_dir(mkdir=FALSE), paste0(getwd(), "/output"))
  expect_equal(get_output_dir(start_path="output", experiment="test", mkdir=FALSE, plot_type = "plot"), paste0(getwd(), "/output/test/plot/", format(Sys.time(), "%Y-%m-%d_%s")))
})

test_that("object handler works", {
  expect_true(!is.null(object_handler(phyloseq_silva)))
  expect_error(object_handler(obj = "data/raw_silva.rda"), "object 'metacoder_object' not found")
  expect_error(object_handler(obj = NULL), "Please use a metacoder/phyloseq object or an rdata file.")

})

test_that("transposer works", {
  expect_true(!is.null(transposer(basic_silva$data$taxa_abundance, ids = "taxon_id", header_name = "samples")))
  expect_true(!is.null(transposer(basic_silva$data$taxa_abundance, ids = "taxon_id", header_name = "samples", preserved_categories = FALSE)))
  expect_true(!is.null(transposer(basic_silva$data$taxa_abundance, header_name = "samples")))
  expect_error(transposer(basic_silva$data, header_name = "samples"), "Data not transposable.")
})

test_that("transformer works", {
  expect_true(!is.null(transformer(basic_silva$data$taxa_abundance, function(x)x/sum(x))))
  expect_true(!is.null(transformer(basic_silva$data$taxa_abundance, function(x)x/sum(x), preserved_categories = FALSE)))
})

test_that("create_pub_table works", {
  expect_true(!is.null(create_pub_table(basic_silva$data$taxa_abundance[1:5, 1:5])))
})

test_that("vlookup works", {
  expect_true(!is.null(vlookup(lookup_vector = c("ab", "ac"), df = basic_silva$data$taxa_abundance, match_var = "taxon_id", return_var = "Sample_1")))
})

# Remove direcotry created by test
if (dir.exists("output/test")) {
  unlink("output", recursive = TRUE)
}

# Remove file created by test
if (file.exists("Rplots.pdf")) {
  file.remove("Rplots.pdf")
}
