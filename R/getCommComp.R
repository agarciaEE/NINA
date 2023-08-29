#' Get community composition
#'
#' @param x Ecological Niche Model of class NINA.
#'
#' @return data frame conatining the presence or absence of species (columns) in each region (rows) from a NINA model
#' @export
#'
#' @examples
#' EN <- EN_model(env_data, occ_data1, cluster = "env", n.clus = 5)
#' CommComp <- getCommComp(EN)
getCommComp <- function(x){

  if (sum(class(x) %in% c("NINA", "ENmodel", "BCmodel", "ECmodel")) != 2) {
    stop("Argument 'x' is not a NINA species model")
  }
  return(x$tab)
}
