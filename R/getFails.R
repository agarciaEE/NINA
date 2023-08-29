#' Get failed models
#'
#' @param x Ecological Niche Model of class NINA.
#'
#' @return Data frame containing failed models
#' @export
#'
#' @examples
#' EN <- EN_model(env_data, occ_data1, cluster = "env", n.clus = 5)
#' EN_failed <- getFails(EN)
getFails<- function(x){

  if (sum(class(x) %in% c("NINA", "ENmodel", "BCmodel", "ECmodel")) != 2) {
    stop("Argument 'x' is not a NINA species model")
  }
  return(x$fail)
}
