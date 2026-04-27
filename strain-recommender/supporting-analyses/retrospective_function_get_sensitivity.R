# ==============================================================================
# get_sensitivity
# ==============================================================================

#' Compute per-disease sensitivity and summary statistics
#'
#' Iterates over each unique disease in the input data frame and computes a
#' set of classification performance and descriptive statistics. Optionally
#' computes a binomial confidence interval for diseases with more than one
#' observation (currently commented out in the returned columns).
#'
#' @param x A data frame with one row per genotype-disease comparison,
#'   containing at minimum the following columns:
#'   \describe{
#'     \item{disease_name}{Character. Name of the disease; used to group rows.}
#'     \item{rare}{Logical or integer. Whether the disease is considered rare;
#'       taken from the first row of each disease group.}
#'     \item{correct}{Numeric (0/1). Whether the comparison was correctly
#'       classified.}
#'     \item{distance_pred}{Numeric. Predicted distance score for each
#'       comparison; used to summarise missed (incorrect) cases.}
#'     \item{gs_id}{Character. Genotype set identifier.}
#'     \item{jnum_id}{Character. J-number reference identifier.}
#'     \item{genes_used}{Numeric. Number of genes used in the comparison.}
#'     \item{tier}{Numeric. Tier classification (e.g. 3 or 4).}
#'   }
#'
#' @return A data frame with one row per unique disease and the following
#'   columns:
#'   \describe{
#'     \item{rare}{Whether the disease is rare (from first row of group).}
#'     \item{disease}{Disease name.}
#'     \item{sensitivity}{Proportion of comparisons correctly classified
#'       (\code{mean(correct)}).}
#'     \item{n_correct}{Total number of correct comparisons.}
#'     \item{n_missed_mean_distance}{Mean absolute predicted distance for
#'       missed (incorrect) comparisons. \code{NaN} if no missed cases.}
#'     \item{n_comparisons}{Total number of comparisons for the disease.}
#'     \item{n_gs}{Number of unique genotype sets.}
#'     \item{median_n_genes}{Median number of genes used across comparisons.}
#'     \item{n_studies}{Number of unique J-number references.}
#'     \item{n_tier3}{Number of comparisons at tier 3.}
#'     \item{n_tier4}{Number of comparisons at tier 4.}
#'   }
#'
#' @note The binomial confidence interval (\code{lower}, \code{upper}) computed
#'   via \code{\link[stats]{binom.test}} is calculated internally for diseases
#'   with more than one comparison but is currently excluded from the returned
#'   data frame (commented out). Only diseases with a single comparison receive
#'   \code{NA} bounds.
#'
#' @examples
#' df <- data.frame(
#'   disease_name  = c("Disease X", "Disease X", "Disease Y"),
#'   rare          = c(TRUE, TRUE, FALSE),
#'   correct       = c(1, 0, 1),
#'   distance_pred = c(0.1, 0.8, 0.2),
#'   gs_id         = c("gs1", "gs2", "gs3"),
#'   jnum_id       = c("J001", "J002", "J003"),
#'   genes_used    = c(10, 12, 8),
#'   tier          = c(3, 4, 3),
#'   stringsAsFactors = FALSE
#' )
#' get_sensitivity(df)
get_sensitivity <- function(x) {
  keep <- data.frame()
  diseases <- unique(x$disease_name)
  
  for (disease in diseases) {
    is_d <- which(x$disease_name == disease)
    test <- x[is_d, ]
    is_rare <- test$rare[1]
    lower <- upper <- NA
    
    if (nrow(test) > 1) {
      b <- binom.test(x = sum(test$correct), n = nrow(test))
      lower <- b$conf.int[1]
      upper <- b$conf.int[2]
    }
    
    keep <- rbind(keep, data.frame(
      rare                  = is_rare,
      disease               = disease,
      sensitivity           = mean(test$correct),
      n_correct             = sum(test$correct),
      n_comparisons         = nrow(test),
      n_gs                  = length(unique(test$gs_id)),
      median_n_genes        = median(test$genes_used),
      n_studies             = length(unique(test$jnum_id)),
      n_tier3               = length(which(test$tier == 3)),
      n_tier4               = length(which(test$tier == 4))
    ))
  }
  
  keep
}


# Unit tests in unit_tests/