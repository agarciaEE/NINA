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
#' @export
niche_to_dis <- function(env.scores, z, cor = FALSE, cluster = NULL, w = NULL){

  if (!is.null(w)){
    env.scores[,3] <- env.scores[,3]*w
    env.scores[,4] <- env.scores[,4]*w
  }
  if (is.null(cluster)){
    if (cor){
      df <- cbind(env.scores[,1:2], vals = raster::extract(z$z.cor, env.scores[,3:4]))
    }
    else{
      df <- cbind(env.scores[,1:2], vals = raster::extract(z$z.uncor, env.scores[,3:4]))
    }
  }
  if (is.data.frame(cluster)){
    if (any(names(z) %in% unique(cluster[,3]))){
      if (any(!names(z) %in% unique(cluster[,3]))){
        warning("Some regions of z are not included in the provided cluster.")
      }
      df = cbind(env.scores[,1:2], vals = 0)
      for (e in names(z)){
        reg <- which(cluster[,3] == e)
        if (cor){
          df[reg,3] <- raster::extract(z[[e]]$z.cor, env.scores[reg,3:4])
        }
        else{
          df[reg,3] <- raster::extract(z[[e]]$z.uncor, env.scores[reg,3:4])
        }
      }
    }
  }
  df[is.na(df)] = 0
  return(df)
}
