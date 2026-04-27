#' Generate Concordance Null Distribution
#'
#' This function generates a null distribution of concordance scores by bootstrapping.
#'
#' @param full_df A data frame where rows are genes and columns are samples.
#' @param gs A data frame containing gene sets with `gene_id` and `score` columns.
#' @param B An integer specifying the number of bootstrap samples. Default is 4000.
#' @param seed An integer specifying the random seed for reproducibility.
#'
#' @return A numeric vector of concordance scores.
#'
#' @examples
#' \dontrun{
#' full_df <- data.frame(matrix(rnorm(1000), nrow=100, ncol=10))
#' rownames(full_df) <- paste0("gene", 1:100)
#' colnames(full_df) <- paste0("sample", 1:10)
#' gs <- data.frame(gene_id=paste0("gene", sample(1:100, 50)), score=rnorm(50))
#' generate_concordance_null_distribution(full_df, gs, B=1000)
#' }
#'
generate_concordance_null_distribution <- function(full_df, gs, B=4000, seed=303) {
  gs <- gs[which(gs$gene_id %in% rownames(full_df)), ]
  size <- nrow(gs)
  rho_boot <- rep(NA, B)
  set.seed(seed)
  for (b in 1:B) {
    genes <- sample(rownames(full_df), size=size, replace=FALSE)
    s <- sample(colnames(full_df), size=1)
    df_boot <- full_df[genes, s]
    rho_boot[b] <- cor(as.vector(df_boot, "numeric"), as.vector(gs$score, "numeric"),
                       method="spearman", use="pairwise.complete.obs")
  }
  rho_boot
}
