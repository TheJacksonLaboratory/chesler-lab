
library(testthat)

# Sample data for testing
gs_valid <- data.frame(gene_id = c("ENSMUSG0000001", "ENSMUSG0000002"),
                       score = c(1.5, 2.3))

gs_invalid_colnames <- data.frame(gene_id = c("ENSMUSG0000001"),
                                  value = c(1.5))

gs_invalid_gene_ids <- data.frame(gene_id = c("INVALID_ID", "ENSMUSG0000002"),
                                  score = c(1.5, 2.3))

gs_with_duplicates <- data.frame(gene_id = c("ENSMUSG0000001", "ENSMUSG0000001", "ENSMUSG0000002"),
                                 score = c(1.5, -2.0, 2.3))

gs_non_numeric_scores <- data.frame(gene_id = c("ENSMUSG0000001", "ENSMUSG0000002"),
                                    score = c("high", "low"))

# Test cases
test_that("Function checks for required column names", {
  expect_error(preprocess_gs(gs_invalid_colnames),
               "ERROR: column names must contain 'gene_id' and 'score'")
})

test_that("Function removes duplicates and keeps highest absolute score", {
  result <- preprocess_gs(gs_with_duplicates, verbose = TRUE)
  expect_equal(nrow(result), 2)  # Should keep two unique gene_ids
  expect_equal(result$score[result$gene_id == "ENSMUSG0000001"], -2.0)  # Should keep the higher absolute score
})

test_that("Function checks for numeric scores", {
  expect_error(preprocess_gs(gs_non_numeric_scores),
               "ERROR: scores are not numeric. Scores are of class character")
})

test_that("Function returns a data frame with correct columns", {
  result <- preprocess_gs(gs_valid)
  expect_true(is.data.frame(result))
  expect_true(all(c("gene_id", "score") %in% colnames(result)))
})

test_that("Function prints the correct number of genes", {
  expect_output(preprocess_gs(gs_valid, verbose = TRUE),
                "number of genes in geneset: 2")
})

test_that("Function prints the correct range of scores", {
  expect_output(preprocess_gs(gs_valid, verbose = TRUE),
                "Range of scores: 1.5, 2.3")
})
