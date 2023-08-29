#' Get input data
#'
#' @param x Ecological Niche Model of class NINA.
#'
#' @return list of input data provided.
#' @export
#'
#' @examples
#' EN <- EN_model(env_data, occ_data1, cluster = "env", n.clus = 5)
#' EN_failed <- getData(EN)
getData<- function(x){

  if (sum(class(x) %in% c("NINA", "ENmodel", "BCmodel", "ECmodel")) != 2) {
    stop("Argument 'x' is not a NINA species model")
  }
  out <- list(env = x$pca$tab,
              sp = x$obs,
              clus = x$clus,
              A = x$A,
              C = x$C)

  return(out)
}
