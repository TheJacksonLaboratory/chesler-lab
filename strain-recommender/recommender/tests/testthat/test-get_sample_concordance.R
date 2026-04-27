
library(testthat)

# Sample data for testing
df_test <- data.frame(
  gene_id = c("ENSMUSG0000001", "ENSMUSG0000002", "ENSMUSG0000003"),
  sample1 = c(1.2, 2.3, 3.1),
  sample2 = c(1.5, 2.1, 3.3)
)
rownames(df_test) <- df_test$gene_id

gs_test <- data.frame(
  gene_id = c("ENSMUSG0000001", "ENSMUSG0000002"),
  score = c(0.5, 0.8)
)

meta_test <- data.frame(
  mouse_id = c("sample1", "sample2"),
  strain = c("strainA", "strainB"),
  sex = c("both", "both")
)

# Test cases
test_that("Function returns a data frame", {
  result <- get_sample_concordance(df_test, gs_test, meta_test, samples = c("sample1", "sample2"), sex = "both")
  expect_true(is.data.frame(result))
})

test_that("Function calculates correlation correctly", {
  result <- get_sample_concordance(df_test, gs_test, meta_test, samples = c("sample1", "sample2"), sex = "both")
  expect_equal(nrow(result), 2)  # Should return one row per strain
  expect_true(all(c("strain", "sex", "corr", "n") %in% colnames(result)))
})

test_that("Function handles missing samples gracefully", {
  result <- get_sample_concordance(df_test, gs_test, meta_test, samples = c("sample1"), sex = "both")
  expect_equal(nrow(result), 1)  # Should return one row for the existing strain
})

test_that("Function handles no matching gene_ids", {
  gs_no_match <- data.frame(gene_id = c("ENSMUSG0000004"), score = c(0.9))
  expect_error(get_sample_concordance(df_test, gs_no_match, meta_test, samples = c("sample1", "sample2"), sex = "both"),
               "ERROR: There are no genes in both gs and df...stopping")  # No matching gene_ids should return an empty data frame
})

test_that("Function handles NA values in df", {
  df_with_na <- data.frame(
    gene_id = c("ENSMUSG0000001", "ENSMUSG0000002", "ENSMUSG0000003"),
    sample1 = c(1.2, NA, 3.1),
    sample2 = c(1.5, 2.1, 3.3)
  )
  rownames(df_with_na) <- df_with_na$gene_id
  result <- get_sample_concordance(df_with_na, gs_test, meta_test, samples = c("sample1", "sample2"), sex = "both")
  expect_true(all(is.na(result$score)))
  # Correlation should be calculated even with NA
})
