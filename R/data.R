#' mzrollDB
#'
#' mzrollDB file for Metabolomics Workbench ST003519 dataset
#'
#' @details
#' The experiment contains 3 experimental conditions
#' \itemize{
#'   \item{control, ischemia, insuline_ischemia}
#' }
#'
#' @return path to the .mzrollDB dataset
#'
#' @family ST003519
#'
#' @export
ST003519_mzroll <- function() {
  path <- system.file("extdata", "ST003519.mzrollDB", package = "metstats")
  if (path == "" || !file.exists(path)) {
    stop("ST003519 mzrollDB file was not found")
  }

  return(path)
}


#' metadata
#'
#' metadata example for Metabolomics Workbench ST003519 dataset
#' \describe{
#'   \item{SampleName}{sub string of sample filename, must be a unique match for one sample}
#'   \item{Condition}{condition that the samle belongs to}
#'   \item{RefCondition}{reference condition for the condition that the samle belongs to}
#'   \item{Replicate}{replicate number for the sample}
#'   \item{treatment}{example of a variable, other dataset variables must be included as additional columns}
#' }
#' @family ST003519
#'
#'
#' @export
"ST003519_metadata"

#' id mapping
#'
#' KEGG id mapping for the features of Metabolomics Workbench ST003519 dataset
#'
#' @format A data frame:
#' \describe{
#'   \item{feature_name}{feature names in the dataset}
#'   \item{id}{KEGG id for features}
#'   \item{short_name}{short name for fatures}
#' }
#' @family ST003519
#'
#' @export
"ST003519_id_mapping"
