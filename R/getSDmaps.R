#' Get Species Distribution Models
#'
#' @param x Ecological Niche Model of class NINA.
#'
#' @return Raster Stack containing Species Distribution Models from NINA models.
#' @export
#'
#' @examples
#' EN <- EN_model(env_data, occ_data1, cluster = "env", n.clus = 5)
#' SDMs <- getSDmaps(EN)
getSDmaps <- function(x){

  if (sum(class(x) %in% c("NINA", "ENmodel", "BCmodel", "ECmodel")) != 2) {
    stop("Argument 'x' is not a NINA species model")
  }
  return(x$maps)
}
