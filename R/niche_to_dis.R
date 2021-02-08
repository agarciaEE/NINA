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
      df <- cbind(env.scores[,1:2], vals = raster::extract(z$z.cor, z[[e]]$glob))
    }
    else{
      df <- cbind(env.scores[,1:2], vals = raster::extract(z$z.uncor, z[[e]]$glob))
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
          df[reg,3] <- raster::extract(z[[e]]$z.cor, z[[e]]$glob)
        }
        else{
          df[reg,3] <- raster::extract(z[[e]]$z.uncor, z[[e]]$glob)
        }
      }
    }
  }
  df[is.na(df)] = 0
  return(df)
}
