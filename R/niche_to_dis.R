#' Niche to distribution
#'
#' @param env.scores Envirronmental PCA scores
#' @param z niche model
#' @param cor Boolean
#' @param cluster data frame with clustered data
#' @param w Weighting coefficients
#'
#' @examples
#' \dontrun{
#' accident_2015 <- fars_read("Project/data/accident_2015.csv.bz2")
#' }
#'
#' @importFrom raster extract
#'
#' @return data frame
#'
#' @keywords internal
#' @noRd
#'
niche_to_dis <- function(env.scores, z, cor = FALSE, cluster = NULL){

  if (is.null(cluster)){
    if (cor){
      vals = raster::extract(z$z, z$glob) / raster::extract(z$Z, z$glob)
      df <- cbind(env.scores[,1:2], vals)
    }
    else{
      df <- cbind(env.scores[,1:2], vals = raster::extract(z$z, z$glob))
    }
  }
  if (is.data.frame(cluster)){
    if (any(names(z) %in% levels(cluster[,3]))){
      if (any(!names(z) %in% levels(cluster[,3]))){
        warning("Some regions of z are not included in the provided cluster.")
      }
      df = cbind(env.scores[,1:2], vals = 0)
      for (e in names(z)){
        reg <- which(cluster[,3] == e)
        if (cor){
          df[reg,3] <- raster::extract(z[[e]]$z, z[[e]]$glob) / raster::extract(z[[e]]$Z, z[[e]]$glob)
        }
        else{
          df[reg,3] <- raster::extract(z[[e]]$z, z[[e]]$glob)
        }
      }
    }
  }
  df[,3] = df[,3] / max(df[,3], na.rm = T)
  df[is.na(df)] = 0
  return(df)
}
