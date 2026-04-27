
library(testthat)

# Sample data for testing
full_df <- data.frame(
  gene_id = c("ENSMUSG0000001", "ENSMUSG0000002", "ENSMUSG0000003"),
  sample1 = c(1.2, 2.3, 3.1),
  sample2 = c(1.5, 2.1, 3.3)
)
rownames(full_df) <- full_df$gene_id
full_df <- full_df[, -1]

gs_test <- data.frame(
  gene_id = c("ENSMUSG0000001", "ENSMUSG0000002"),
  score = c(0.5, 0.8)
)

# Test cases
test_that("Function returns a numeric vector", {
  result <- generate_concordance_null_distribution(full_df, gs = gs_test, B = 100)
  expect_true(is.numeric(result))
})

test_that("Function returns the correct length of output", {
  B <- 100
  result <- generate_concordance_null_distribution(full_df, gs_test, B = B)
  expect_equal(length(result), B)  # Output length should match B
})

test_that("Function produces consistent results with the same seed", {
  result1 <- generate_concordance_null_distribution(full_df, gs_test, B = 100, seed=421)
  result2 <- generate_concordance_null_distribution(full_df, gs_test, B = 100, seed=421)
  expect_equal(result1, result2)  # Results should be identical
})

test_that("Function handles large B values", {
  B <- 10000
  result <- generate_concordance_null_distribution(full_df, gs_test, B = B)
  expect_equal(length(result), B)  # Output length should match B
})
