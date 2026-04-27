#' Get Sample Concordance
#'
#' Calculates the Spearman correlation between gene expression data and gene
#' scores for each strain.
#'
#' @param df A data frame containing gene expression data with 'gene_id' as
#'   row names.
#' @param gs A data frame containing gene scores with at least two columns:
#'   'gene_id' and 'score'.
#' @param meta A data frame containing metadata with at least three columns:
#'   'mouse_id', 'sex', and 'strain'.
#' @param samples A vector of sample IDs to be included in the analysis.
#' @param sex A character string indicating the sex of the samples.
#'   Default is 'both'.
#'
#' @return A data frame containing the strain, sex, Spearman correlation
#'   coefficient, and the number of genes used in the calculation.
#'
#' @details
#' The function performs the following steps:
#' - Filters the gene scores (\code{gs}) to include only those expressed in
#'   the data (\code{df}).
#' - Subsets the gene expression data (\code{df}) to include only the genes
#'   present in the filtered gene scores.
#' - Iterates over each unique strain in the metadata (\code{meta}).
#' - For each strain, identifies the samples that match the given strain.
#' - Calculates the mean expression for each gene across the identified samples
#'   if there are multiple samples; otherwise, uses the single sample's expression.
#' - Computes the Spearman correlation between the mean gene expression and the
#'   gene scores.
#'
#' @keywords internal
get_sample_concordance <- function(df, gs, meta, samples, sex) {

  out <- data.frame()

  # Filter genes to only include those expressed in df
  gs <- gs[which(gs$gene_id %in% rownames(df)), ]

  # Make sure there are genes after filtering
  condition <- nrow(gs) == 0
  if (condition) { stop("ERROR: There are no genes in both gs and df...stopping")}

  df <- df[gs$gene_id, ]

  # remove NAs if there are any
  na_rows <- which(apply(df, 1, function(x) any(is.na(x))))
  condition <- length(na_rows) > 0
  if (condition) {
    # TODO fix this so it will still run per sample.
    df <- df[-na_rows, ]; gs <- gs[-na_rows, ]
    print(paste0("WARNING: Removed ", length(na_rows), " rows due to NAs for some samples."))
    if (nrow(df) == 0) stop("ERROR: After removing NAs, there are no genes in df.")
  }

  # Get strains for the samples from meta
  strains <- unique(meta$strain[meta$mouse_id %in% samples])


  for (strain in strains) {
    # strain_samples are the mouse_ids where the strain is the strain in the loop
    strain_samples <- meta$mouse_id[which(meta$mouse_id %in% samples & meta$strain==strain)]

    # If there is more than one mouse_id annotated to the strain,
    #   let x be the mean of imputed expression for those mouse_ids.
    # Otherwise, use the expression for the mouse_id annotated to the strain
    if (length(strain_samples) > 1) {
      x <- rowMeans(df[, strain_samples], na.rm=TRUE)
    } else {
      x <- df[, strain_samples]
    }

    # Calculate Spearman's rho correlation between
    #   (x) the strain's imputed gene expression genes in gs, and
    #   (y) the gs scores for the same genes (disease state)
    rho <- cor(x, gs$score, method="spearman", use="pairwise.complete.obs")

    # add results to the output
    out <- rbind(out, data.frame(strain=strain, sex=sex, corr=rho, n=nrow(df)))
  }
  out
}
