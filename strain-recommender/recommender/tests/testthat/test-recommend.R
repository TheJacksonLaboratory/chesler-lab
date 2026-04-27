
library(testthat)
library(dplyr)

# Sample data for testing
gs_test <- data.frame(gene_id = c("ENSMUSG0000001", "ENSMUSG0000002"),
                      score=c(1, -4))
meta_test <- data.frame(
  mouse_id = c("sample1", "sample2"),
  sex      = c("both", "both"),
  strain   = c("strain1", "strain2")
)
df_test <- data.frame(gene_id = c("ENSMUSG0000001", "ENSMUSG0000002"),
                      sample1 = c(1.2, 2.3),
                      sample2 = c(3.4, 4.5))

# Test cases
test_that("Function handles incorrect input types", {
  expect_error(recommend(gs = "not a dataframe", meta = meta_test, df = df_test))
})
test_that("Function checks for missing gene_id column", {
  df_missing_gene_id <- df_test[, -1]  # Remove gene_id column
  expect_error(recommend(gs = gs_test, meta = meta_test, df = df_missing_gene_id),
               "ERROR: 'gene_id' must be a column name")
})
test_that("Function checks for non-Ensembl gene IDs", {
  df_invalid_gene_id <- data.frame(gene_id = c("INVALID_ID"),
                                   sample1 = c(1.2),
                                   sample2 = c(3.4))
  expect_error(recommend(gs = gs_test, meta = meta_test, df = df_invalid_gene_id, verbose = TRUE))
})
test_that("Function checks for NA values in df", {
  df_with_na <- data.frame(gene_id = c("ENSMUSG0000001", "ENSMUSG0000002"),
                           sample1 = c(1.2, NA),
                           sample2 = c(3.4, 4.5))
  expect_error(recommend(gs = gs_test, meta = meta_test, df = df_with_na),
               "ERROR: no NA values for expression are allowed")
})
test_that("Function checks for numeric entries in df", {
  df_non_numeric <- data.frame(gene_id = c("ENSMUSG0000001", "ENSMUSG0000002"),
                               sample1 = c("a", 2.3),
                               sample2 = c(3.4, 4.5))
  expect_error(recommend(gs = gs_test, meta = meta_test, df = df_non_numeric),
               "ERROR: not all entries in df are numeric.")
})
test_that("Function checks for required columns in meta", {
  meta_missing_columns <- data.frame(mouse_id = c("sample1", "sample2"))
  expect_error(recommend(gs = gs_test, meta = meta_missing_columns, df = df_test),
               "ERROR: 'sex' and 'strain' must columns in meta")
})
test_that("Function returns a concordance data frame", {
  result <- recommend(gs = gs_test, meta = meta_test, df = df_test)
  expect_true(is.data.frame(result))
  expect_true("fpr" %in% colnames(result))
})
