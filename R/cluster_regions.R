#' @title CLUSTERING FUNCTION
#'
#' @description Cluster a data frame with two first columns indicating long and lat coordinates
#'
#' @param res Resolution of the greographical grid
#' @param plot Boolean whether to plot the cluster as map
#' @param x data frame to be clustered
#' @param n.clus Number of clusters to perform. If null it will be computed
#' @param nstart burn up start value
#' @param K.max maximum number of clusters to estimate
#' @param B number of runs
#'
#' @return Clustered data frame.
#'
#'@details Returns an error if \code{filename} does not exist.
#'
#' @examples
#' \dontrun{
#' accident_2015 <- fars_read("Project/data/accident_2015.csv.bz2")
#' }
#'
#' @importFrom stats na.exclude kmeans
#' @importFrom NbClust NbClust
#' @importFrom raster rasterize raster
#'
#' @export
cluster_regions <- function(x, plot = TRUE, n.clus = NULL, nstart = 25, K.max = NULL, B = 100, res = 0.083){

  if (class(x) %in% c("raster", "RasterBrick", "RasterStack")){
    df <- na.exclude(raster::as.data.frame(x, xy = T))
  }
  if (is.data.frame(x)){
    df = na.exclude(x)
  }
  if (is.null(n.clus)){
    if (is.null(K.max)){
      K.max = ncol(df[,-c(1:2)])
    }
    message("\n- Parameters:", "\n    nstart: ", nstart, "\n    K max: ", K.max, "\n    B: ", B)
    if ( B < 100){
      warning("Low number of bootstraps selected, not recommended for final analyses")
    }
    nbclus <- NbClust(data = df, distance = "euclidean",
                      min.nc = 2, max.nc = K.max, method = "kmeans")
    n.clus <- as.numeric(names(table(nbclus$Best.nc[1,]))[which(table(nbclus$Best.nc[1,])
                                                                == max(table(nbclus$Best.nc[1,])))])
    message("\t- Number of clusters selected =", n.clus, "\n")
  }
  clus.df <- as.data.frame(cbind(df[,1:2], cluster = kmeans(df,n.clus)$cluster))
  if (plot == TRUE){
    nc = round(diff(range(clus.df[,1])) / res,0)
    nr = round(diff(range(clus.df[,2])) / res,0)
    plot(rasterize(clus.df[,1:2], raster(ncol = nc, nrow = nr), field = clus.df[,3]), col = viridis::viridis(n.clus))
  }
  return(clus.df)
}
