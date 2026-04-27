#' Recommend Vulnerability Scores
#'
#' This function recommends vulnerability scores for each sample based on gene expression data and gene scores.
#'
#' @param gs A data frame containing gene scores with at least two columns: 'gene_id' and 'score'.
#'   If a 'geneset' column is present and `combine=TRUE`, scores will be normalized within each
#'   geneset before processing.
#' @param meta A data frame containing metadata with at least three columns: 'mouse_id', 'sex', and 'strain'.
#' @param df A data frame containing gene expression data with 'gene_id' as one of the columns.
#' @param verbose A logical value indicating whether to print detailed messages. Default is TRUE.
#' @param combine_gs A logical value indicating if there are multiple gene sets to combine before running the function.
#' If so, gs data frame must contain genesets column to indicate different gene sets. Default is FALSE.
#' @param seed An integer specifying the random seed for reproducibility.
#'   Default is 303.
#' @param B An integer specifying the number of bootstrap samples for the null
#'   distribution. Default is 4000.
#'
#' @return A data frame containing the vulnerability scores and false positive rates (FPR) for each sample/strain.
#'
#' @importFrom dplyr select mutate
#' @importFrom magrittr %>%
#'
#' @details
#' The function performs several checks and preprocessing steps:
#' - Loads necessary libraries and sources required functions.
#' - Preprocesses the gene set using the `preprocess_gs` function.
#' - Ensures that `gs`, `meta`, and `df` are data frames.
#' - Checks that `df` contains a 'gene_id' column with valid mouse Ensembl gene IDs.
#' - Sets row names of `df` to 'gene_id' and removes the 'gene_id' column.
#' - Ensures that all column names of `df` are present in `meta$mouse_id`.
#' - Checks for the presence of 'sex' and 'strain' columns in `meta`.
#' - Subsets `meta` to include only samples present in `df`.
#'
#' The function then calculates vulnerability scores for each sample by:
#' - Calculating concordance for each sample using the `get_sample_concordance` function.
#' - Generating a null distribution of Spearman correlation coefficients using the `generate_concordance_null_distribution` function.
#' - Calculating the false positive rate (FPR) for each vulnerability score.
#'
#' @examples
#' gs <- data.frame(gene_id = c("ENSMUSG00000000001", "ENSMUSG00000000003"),
#'                  score = c(1.5, 2.3))
#' meta <- data.frame(mouse_id = c("sample1", "sample2"),
#'                    sex = c("both", "both"),
#'                    strain = c("strain1", "strain2"))
#' df <- data.frame(gene_id = c("ENSMUSG00000000001", "ENSMUSG00000000003"),
#'                  sample1 = c(5.2, 3.8), sample2 = c(4.1, 2.9))
#' result <- recommend(gs, meta, df, verbose=TRUE)
#'
#' @export
recommend <- function(gs, meta, df, seed=303, B=4000, combine_gs = FALSE, verbose=TRUE) {

  # Preprocess the geneset
  if (combine_gs) {
    gs <- combine_genesets(gs)
  } else {
    gs <- preprocess_gs(gs, verbose=verbose)
  }

  colnames(meta)[1] <- "mouse_id"

  # Check for errors in data types
  if (!all(c(class(df), class(gs), class(meta)) == "data.frame")) {
    stop(print("ERROR: gs, df, and meta must be dataframes."))
  }

  # Check for errors in df
  if (!'gene_id' %in% colnames(df)) {
    stop(print("ERROR: 'gene_id' must be a column name"))
  }
  if (!all(grepl("^ENSMUSG", df$gene_id))) {
    print("WARNING: not all gene_ids are mouse Ensembl gene ids. Removing...")
    remove <- which(!grepl("^ENSMUSG", df$gene_id))
    df <- df[-remove,]
    if (verbose) print(paste0("Removed ", length(remove), " gene ids that were not mouse Ensembl genes"))
    if (nrow(df) == 0) {
      stop(print("ERROR: no mouse Ensembl ids in the geneset"))
    }
  }

  # Set row names and remove 'gene_id' column from df
  rownames(df) <- df$gene_id
  df <- df[, -which(colnames(df)=='gene_id')]

  # Check if all column names of df are in meta
  if (!all(colnames(df) %in% meta$mouse_id)) {
    stop(print("ERROR: all column names of df except 'gene_id' must be samples in meta"))
  } else if (any(is.na(df))) {
    stop(print("ERROR: no NA values for expression are allowed"))
  } else if (!all(apply(df, 2, is.numeric))) {
    stop(print("ERROR: not all entries in df are numeric."))
  }

  # Check for errors in meta
  if (!all(c('sex','mouse_id') %in% colnames(meta))) {
    stop(print("ERROR: 'sex' and 'strain' must columns in meta"))
  } else if (!(all(meta$sex %in% c('m','f','both')))) {
    stop(print("ERROR: entries in meta$sex must be 'm', 'f', and/or 'both'"))
  }

  # Subset meta to only those in df
  meta <- meta[which(meta$mouse_id %in% colnames(df)), ]
  if (verbose) print("Input data passed all tests. Proceeding to recommendations...")

  # Function to calculate FPR from generated null distribution
  get_fpr <- function(x) {
    if (x >= 0) y <- mean(rho_null >= x)
    if (x < 0) y <- mean(rho_null <= x)
    y
  }

  # Calculate concordance
  if (verbose) print("Calculating vulnerability scores for each sample...")
  #sexes <- unique(meta$sex)
  sexes <- 'both'
  concordance <- data.frame()
  for (sex in sexes) {
    samples <- meta$mouse_id[which(meta$sex==sex)]
    concordance <- rbind(concordance,
                         get_sample_concordance(df=df, gs=gs, meta=meta, samples=samples, sex=sex))
  }

  # Calculate FPR
  if (verbose) print("Calculating uncertainty for each vulnerability score...")
  rho_null <- generate_concordance_null_distribution(full_df=df, gs=gs, B=B, seed=seed)
  fpr <- sapply(concordance$corr, get_fpr)
  concordance$fpr <- fpr
  if (verbose) print("Completed with no errors.")

  concordance <- concordance %>%
    dplyr::select(-"sex") %>%
    dplyr::mutate(score = corr / max(abs(corr)))

  concordance <- concordance %>% select(strain, score, fpr, spearman_rho = corr, n_genes = n)
  concordance
}

