#' @title Sample of Occurrence Data of Species Group Two for Analysis
#'
#' @description  A \code{\link{data.frame}} with spatial coordinates of Species Group Two. Occurrences have been sampled using \code{\link[virtualspecies]{sampleOccurrences}} from generated niches using \code{\link[virtualspecies]{generateRandomSp}} of which generated niche is nested within the niche of their interactions. See \code{\link{int_matrix}} to see species interaction matrix.
#'
#' @format A \code{\link{data.frame}} with three columns, which are:
#' \describe{
#' \item{x}{Spatial longitud coordinates}
#' \item{y}{Spatial latitude coordinates}
#' \item{species}{species names belonging to each long lat coordinate}
#' }
"occ_data2"
