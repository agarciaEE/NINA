#' Get species names from NINA model
#'
#' @param x Ecological Niche Model of class NINA.
#'
#' @return Vector of species names contained in a NINA model
#' @export
#'
#' @examples
#' EN <- EN_model(env_data, occ_data1, cluster = "env", n.clus = 5)
#' speciesNames <- getSpecies(EN)
getSpecies <- function(x){

  if (sum(class(x) %in% c("NINA", "ENmodel", "BCmodel", "ECmodel")) != 2) {
    stop("Argument 'x' is not a NINA species model")
  }

  spsNames <- colnames(x$tab)
  return(spsNames)
}
