# ==============================================================================
# process_data
# ==============================================================================

#' Process Strain Recommender results across tissues
#'
#' For each tissue in the global \code{tissues} vector, this function applies
#' exclusion filters, filters by false positive rate (FPR) threshold, annotates
#' rare diseases, computes overall/rare/common binomial success rates, and
#' derives per-disease sensitivity statistics.
#'
#' @param tissue Character vector of tissue name.
#'   Defaults to \code{"mESC"}.
#' @param rare A data frame with a column \code{"Disease name"} listing
#'   diseases classified as rare.
#' @param thresh Numeric scalar. FPR threshold; rows with \code{fprA > thresh}
#'   are dropped. Defaults to \code{0.05}.
#' @param df_all A data frame containing genotype-disease comparison records
#'   with at minimum the following columns:
#'   \describe{
#'     \item{disease_name}{Character. Name of the associated disease.}
#'     \item{fprA}{Numeric. False positive rate for strain A; rows with
#'       \code{fprA > thresh} are excluded.}
#'     \item{MP_term_name}{Character. Mammalian Phenotype term name; rows with
#'       an empty string are excluded.}
#'     \item{correct}{Numeric (0/1). Whether the comparison was correctly
#'       classified.}
#'     \item{gs_id}{Character. Genotype set identifier.}
#'     \item{strainA}{Character. Name/identifier of the first strain.}
#'     \item{strainB}{Character. Name/identifier of the second strain.}
#'     \item{jnum_id}{Character. J-number reference identifier.}
#'     \item{match_status_strainA}{Character. Match quality for strain A.}
#'     \item{match_status_strainB}{Character. Match quality for strain B.}
#'     \item{distance_pred}{Numeric. Predicted distance score.}
#'     \item{genes_used}{Numeric. Number of genes used in the comparison.}
#'     \item{tier}{Numeric. Tier classification.}
#'   }
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{dd}{The filtered data frame from the \emph{last} tissue
#'       processed, with an added \code{rare} column (0/1) indicating whether
#'       each disease is in the global \code{rare} data frame.}
#'     \item{sensitivities}{A data frame of per-disease sensitivity statistics
#'       for the last tissue processed, as returned by
#'       \code{\link{get_sensitivity}}.}
#'   }
#'
#' @section Side effects:
#' \itemize{
#'   \item Prints progress messages to the console for each tissue, including
#'     row counts, disease counts, and sensitivity summary statistics.
#'   \item Accumulates a \code{success} data frame internally (overall, rare,
#'     and common binomial estimates per tissue), but this is \strong{not}
#'     currently included in the return value.
#' }
#'
#' @section Global dependencies:
#' None. All previously global dependencies (\code{tissues}, \code{thresh},
#' \code{rare}) are now explicit function parameters.
#'
#' @note
#' \itemize{
#'   \item The \code{success} data frame (binomial success rates per tissue)
#'     is computed but not returned. Consider adding it to the return list if
#'     needed downstream.
#'   \item Only the results for the \strong{last} tissue in \code{tissues} are
#'     returned. If per-tissue results are needed, the return value should be
#'     refactored to a list keyed by tissue.
#'   \item The \code{df_all <- df_all} assignment inside the loop is a no-op
#'     and can be safely removed.
#' }
#'
#' @seealso \code{\link{check_exclusions}}, \code{\link{get_sensitivity}}
#'
#' @examples
#' \dontrun{
#' # Assumes tissue, thresh, and rare exist in the environment
#' tissue <- "kidney"
#' thresh  <- 0.05
#' rare    <- data.frame(`Disease name` = c("Disease X"), check.names = FALSE)
#'
#' result <- process_data(df_all)
#' results_filtered <- result[[1]]
#' sensitivities <- result[[2]]
#' }
process_data <- function(df_all, tissue="mESC", rare, thresh=0.05) {
  success <- data.frame()
  
  print(paste0("Processing ", tissue, "..."))
  print(paste0("number of rows = ", nrow(df_all)))
  print(paste0("number of diseases = ", length(unique(df_all$disease_name))))
  
  df <- check_exclusions(df_all)
  print(paste0("after exclusions, number of rows = ", nrow(df)))
  print(paste0("number of diseases = ", length(unique(df$disease_name))))
  
  # Filter by FPR threshold and non-empty MP term
  df %>%
    filter(fprA <= thresh & MP_term_name != "") -> dd
  dd_complete <- dd
  print(paste0("after filtering on FPR, number of rows = ", nrow(dd)))
  print(paste0("number of diseases = ", length(unique(dd$disease_name))))
  
  # Add rare disease indicator
  dd$rare <- as.numeric(dd$disease_name %in% rare$`Disease name`)
  
  # Verify that there is at least one row in dd
  condition <- nrow(dd) > 0
  if (!condition) {
    stop("Error: after exclusions, there are no results to assess. Stopping...")
  }
  
  # Overall success rate
  x <- sum(dd$correct); n <- nrow(dd)
  b <- binom.test(x = x, n = n)
  nd <- length(unique(dd$disease_name))
  success <- rbind(success, data.frame(
    tissue = tissue, correct = x, n = n,
    rate = b$estimate, lower = b$conf.int[1], upper = b$conf.int[2],
    n_disease = nd, label = "overall"
  ))
  
  # Rare disease success rate
  x  <- sum(dd$correct[which(dd$rare == 1)])
  n  <- length(which(dd$rare == 1))
  b  <- binom.test(x = x, n = n)
  nd <- length(unique(dd$disease_name[which(dd$rare == 1)]))
  success <- rbind(success, data.frame(
    tissue = tissue, correct = x, n = n,
    rate = b$estimate, lower = b$conf.int[1], upper = b$conf.int[2],
    n_disease = nd, label = "rare"
  ))
  
  # Non-rare disease success rate
  x  <- sum(dd$correct[which(dd$rare == 0)])
  n  <- length(which(dd$rare == 0))
  b  <- binom.test(x = x, n = n)
  nd <- length(unique(dd$disease_name[which(dd$rare == 0)]))
  success <- rbind(success, data.frame(
    tissue = tissue, correct = x, n = n,
    rate = b$estimate, lower = b$conf.int[1], upper = b$conf.int[2],
    n_disease = nd, label = "common"
  ))
  
  # Per-disease sensitivity
  sensitivities <- get_sensitivity(x = dd)
  print(paste0(
    "median ", median(sensitivities$sensitivity),
    "; IQR ", paste(quantile(sensitivities$sensitivity, probs = c(0.25, 0.75)), collapse = ", "),
    "; mean ", mean(sensitivities$sensitivity)
  ))
  
  return(list(dd, sensitivities))
}

# Unit tests in unit_tests/