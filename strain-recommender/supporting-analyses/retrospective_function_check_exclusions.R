
# ==============================================================================
# check_exclusions
# ==============================================================================

#' Filter and deduplicate strain recommender results
#'
#' Applies a sequence of exclusion rules to remove invalid, incomplete, or
#' duplicate rows from the input data frame:
#'
#' 1. Removes rows where \code{disease_name} or \code{jnum_id} is empty.
#' 2. Removes rows where \code{match_status_strainA} or
#'    \code{match_status_strainB} is not \code{"exact"}.
#' 3. Removes duplicate rows based on the combination of:
#'    \code{gs_id}, \code{strainA}, \code{strainB}, \code{MP_term_name},
#'    \code{jnum_id}, and \code{disease_name}.
#' 4. Removes rows where \code{strainA} and \code{strainB} are identical.
#'
#' @param x A data frame containing at minimum the following columns:
#'   \describe{
#'     \item{gs_id}{Genotype set identifier.}
#'     \item{strainA}{Name/identifier of the first strain.}
#'     \item{strainB}{Name/identifier of the second strain.}
#'     \item{MP_term_name}{Mammalian Phenotype term name.}
#'     \item{jnum_id}{J-number identifier for the reference.}
#'     \item{disease_name}{Name of the associated disease.}
#'     \item{match_status_strainA}{Match quality for strain A; must be
#'       \code{"exact"} to be retained.}
#'     \item{match_status_strainB}{Match quality for strain B; must be
#'       \code{"exact"} to be retained.}
#'   }
#'
#' @return A filtered data frame with invalid, incomplete, and duplicate
#'   rows removed. Row indices are not reset.
#'
#' @examples
#' df <- data.frame(
#'   gs_id               = c("gs1", "gs2"),
#'   strainA             = c("A", "B"),
#'   strainB             = c("C", "D"),
#'   MP_term_name        = c("term1", "term2"),
#'   jnum_id             = c("J001", "J002"),
#'   disease_name        = c("Disease X", "Disease Y"),
#'   stringsAsFactors    = FALSE
#' )
#' check_exclusions(df)
check_exclusions <- function(x) {
  single <- c("gs_id", "strainA", "strainB", "MP_term_name", "jnum_id", "disease_name")
  
  # Rule 1: drop rows with empty disease_name or jnum_id
  exclude <- which(x$disease_name == "" | x$jnum_id == "")
  
  # Rule 2: drop rows where either strain match is not exact
  exclude <- union(exclude, which(x$match_status_strainA != "exact" |
                                    x$match_status_strainB != "exact"))
  if (length(exclude) > 0) x <- x[-exclude, ]
  
  # Rule 3: drop duplicate rows across key columns
  dups <- which(duplicated(x[, single]))
  if (length(dups) > 0) x <- x[-dups, ]
  
  # Rule 4: drop rows where strainA and strainB are the same
  remove <- which(x$strainA == x$strainB)
  if (length(remove) > 0) x <- x[-remove, ]
  
  x
}


# Unit tests in unit_tests/