#' Combine Gene Sets
#'
#' Normalizes scores within each gene set using signed normalization.
#' This can be called on its own so a user can see what the gene set would
#' look like before using as input to `recommend()`. If using the combined
#' data frame to then run `recommend()`, set `combine = FALSE` in the
#' `recommend()` function. Set `combine = TRUE` in `recommend()` function
#' to combine gene sets and run the entire function at once.
#'
#' @param gs A data frame containing gene scores with columns 'gene_id' and
#'   'score', and optionally 'geneset'. If no 'geneset' column is present,
#'   the function passes directly to [preprocess_gs()].
#'
#' @return A cleaned data frame with 'gene_id' and 'score' columns, ready for
#'   viewing and/or use as input to [recommend()].
#'
#' @details
#' If a 'geneset' column is present, the function loops through each unique
#' geneset and normalizes positive and negative scores independently:
#' - Positive scores are divided by the maximum positive score, scaling them
#'   to (0, 1].
#' - Negative scores are divided by the absolute minimum negative score,
#'   scaling them to [-1, 0).
#'
#' This ensures each direction reaches its own maximum of +1 or -1, preserving
#' the sign while placing all scores on a common scale across gene sets.
#'
#' @examples
#' \dontrun{
#' gs <- data.frame(
#'   gene_id = c("ENSMUSG00000000001", "ENSMUSG00000000003",
#'               "ENSMUSG00000000010", "ENSMUSG00000000028"),
#'   score = c(0.5, 2.0, -1.0, -3.0),
#'   geneset = c("setA", "setA", "setB", "setB")
#' )
#' result <- combine_genesets(gs)
#' }
#'
#' @export
combine_genesets <- function(gs) {

  if ("geneset" %in% colnames(gs)) {
    genesets <- unique(gs$geneset)

    for (geneset in genesets) {
      is_geneset <- which(gs$geneset == geneset)

      pos <- which(gs$score[is_geneset] >= 0)
      neg <- which(gs$score[is_geneset] < 0)

      if (length(pos) > 0)
        gs$score[is_geneset[pos]] <- gs$score[is_geneset[pos]] / max(gs$score[is_geneset[pos]])
      if (length(neg) > 0)
        gs$score[is_geneset[neg]] <- gs$score[is_geneset[neg]] / abs(min(gs$score[is_geneset[neg]]))
    }
  }

  preprocess_gs(gs)
}
