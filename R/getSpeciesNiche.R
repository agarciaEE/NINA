#' Get species niches
#'
#' @param x Ecological Niche Model of class NINA.
#'
#' @return List of species niches. If regional models, each species will contain a list of niches per region.
#' @export
#'
#' @examples
#' EN <- EN_model(env_data, occ_data1, cluster = "env", n.clus = 5)
#' SpNiche_List <- getSpeciesNiche(EN)
getSpeciesNiche <- function(x){

  if (sum(class(x) %in% c("NINA", "ENmodel", "BCmodel", "ECmodel")) != 2) {
    stop("Argument 'x' is not a NINA species model")
  }
  if (!is.null(x$clus)){
    out <- reverse_list(x$z.mod)
  } else {
    out <- x$z.mod
  }
  return(out)
}
