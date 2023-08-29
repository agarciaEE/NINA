#' Get distribution suitability
#'
#' @param x Ecological Niche Model of class NINA.
#'
#' @return Data.frame containing the longitude and latitude coordinates (first two columns) and the suitabilities of presence of species (columns) per location (rows).
#' @export
#'
#' @examples
#' EN <- EN_model(env_data, occ_data1, cluster = "env", n.clus = 5)
#' predicted_suitability <- getPredictions(EN)
getPredictions <- function(x){

  if (sum(class(x) %in% c("NINA", "ENmodel", "BCmodel", "ECmodel")) != 2) {
    stop("Argument 'x' is not a NINA species model")
  }
  return(x$pred.dis)
}

