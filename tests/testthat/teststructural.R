context("quad.rarefied")
library(testthat)
library(lotus)

# Test whether the output is a data frame
test_that("structural.specificity() returns a data frame", {
  output_table <- structural.specificity(quad.rarefied, abundance.weighted = FALSE, trim = TRUE)
  expect_is(output_table, "data.frame")
})
