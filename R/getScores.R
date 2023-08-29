#' Get PCA data
#'
#' @param x Ecological Niche Model of class NINA.
#'
#' @return List containing coordinates and loadings from environmental pca, species corresponding  coordinates and loadings, name of predictors and pca object.
#' @export
#'
#' @examples
#' EN <- EN_model(env_data, occ_data1, cluster = "env", n.clus = 5)
#' EN_scores <- getScores(EN)
getScores<- function(x){

  if (sum(class(x) %in% c("NINA", "ENmodel", "BCmodel", "ECmodel")) != 2) {
    stop("Argument 'x' is not a NINA species model")
  }
  out <- list(env = x$env.scores,
              sp = x$sp.scores,
              predictors = x$predictors,
              pca = x$pca)

  return(out)
}


