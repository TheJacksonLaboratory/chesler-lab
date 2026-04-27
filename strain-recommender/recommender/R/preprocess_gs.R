#' Preprocess Gene Scores
#'
#' This function preprocesses a data frame of gene scores by performing several checks and cleaning steps.
#'
#' @param gs A data frame containing gene scores with at least two columns: 'gene_id' and 'score'.
#' @param verbose A logical value indicating whether to print detailed messages during processing. Default is TRUE.
#'
#' @return A cleaned data frame for the input geneset containing only 'gene_id' and 'score' columns of mouse Ensembl genes.
#'
#' @details
#' The function performs the following steps:
#' - Checks if the input data frame contains the required columns 'gene_id' and 'score'.
#' - Prints the number of unique genes in the gene set if `verbose` is TRUE.
#' - Removes rows where 'gene_id' does not match the mouse Ensembl gene ID pattern.
#' - Removes duplicate rows, keeping the score with the highest absolute value for duplicated 'gene_id's.
#' - Checks if the 'score' column is numeric.
#' - Prints a summary of the score range if `verbose` is TRUE.
#' - Returns the cleaned data frame with only 'gene_id' and 'score' columns.
#'
#' @examples
#' \dontrun{
#' gs <- data.frame(gene_id = c("ENSMUSG00000000001", "ENSMUSG00000000003", "INVALID_ID"),
#'                  score = c(1.5, 2.3, 3.1))
#' cleaned_gs <- preprocess_gs(gs, verbose=TRUE)
#' }
#'
#'@keywords internal
#'
preprocess_gs <- function(gs, verbose=TRUE) {
  # Check for gene_id errors
  if (!all(c('gene_id', 'score') %in% colnames(gs))) {
    stop("ERROR: column names must contain 'gene_id' and 'score'")
  }
  # Print the number of unique genes in the geneset if verbose is TRUE
  if (verbose) print(paste0("number of genes in geneset: ", length(unique(gs$gene_id))))
  if (!all(grepl("^ENSMUSG", gs$gene_id))) {
    print("WARNING: not all gene_ids are mouse Ensembl gene ids. Removing...")
    remove <- which(!grepl("^ENSMUSG", gs$gene_id))
    gs <- gs[-remove,]
    if (verbose) print(paste0("Removed ", length(remove), " that were not mouse Ensembl genes"))
    if (nrow(gs) == 0) {
      print("ERROR: no mouse Ensembl ids in the geneset")
      stop()
    }
  }
  # Remove duplicate rows
  gs <- unique(gs)
  dups <- unique(gs$gene_id[which(duplicated(gs$gene_id))])
  if (length(dups) > 0) {
    print(paste0("WARNING: There are ", length(dups), " genes with duplicated gene_ids. Keeping the score with the highest absolute score."))
    for (dup in dups) {
      check <- which(gs$gene_id==dup)
      keep <- check[which.max(abs(gs$score[check]))][1]
      remove <- setdiff(check, keep)
      gs <- gs[-remove, ]
    }
    if (verbose) print(paste0("Removed ", length(dups), " duplicates."))
  }
  if (verbose) print(paste0("number of genes in geneset: ", length(gs$gene_id)))
  # Check for score errors
  if (!is.numeric(gs$score)) {
    c <- class(gs$score)
    stop(paste0("ERROR: scores are not numeric. Scores are of class ", c, "
                 Must be numeric."))
  }
  # Print summary of scores if verbose is TRUE
  if (verbose) print(paste0("Range of scores: ", paste(range(gs$score), collapse=", ")))
  # Export the data.frame with only 'gene_id' and 'score' columns
  gs[, c('gene_id', 'score')]
}
