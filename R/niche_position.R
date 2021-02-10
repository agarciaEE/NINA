#' @title NICHE POSITION
#'
#' @param z niche model
#'
#' @param type type of estimation. Default is "unimodal"
#' @param method Method to estimate the niche centroid. Default is "max"
#' @param quantile Numeric quantile to filter the niche
#' @param cor Boolean to whether use cor or uncor estimation
#'
#' @description Vector
#'
#' @return Data frame.
#'
#'@details Returns an error if \code{filename} does not exist.
#'
#' @examples
#' \dontrun{
#' accident_2015 <- fars_read("Project/data/accident_2015.csv.bz2")
#' }
#'
#' @importFrom stats cov.wt median
#' @importFrom fpc dbscan
#'
#' @export
niche_position <- function(z, type = c("unimodal", "multimodal"),  method = c("max", "median", "mean"), quantile = 0.5, cor = FALSE){

  type = type[1]
  method = method[1]
  R = length(z$x)
  if (cor == TRUE){
    v = raster::as.data.frame(z$z.cor, xy = T)
  }
  if (cor == FALSE){
    v = raster::as.data.frame(z$z.uncor, xy = T)
  }
  v[is.na(v)] = 0
  v = v[v[,3] != 0,]
  qt <- quantile(v[,3], quantile, na.rm = T)
  opt = v[v[,3] >= qt,]
  if (type == "multimodal"){
    opt$cluster <- fpc::dbscan(opt[,1:2], eps = 0.15, MinPts = 5)$cluster
    if (method == "max"){ctr <- t(sapply(unique(opt$cluster), function(x) opt[which.max(opt[opt$cluster == x,3]),1:2]))}
    if (method == "mean"){ctr <- t(sapply(unique(opt$cluster), function(x) cov.wt(opt[opt$cluster == x,1:2], wt = opt[opt$cluster == x,3])$center))  }
    if (method == "median"){ctr <- t(sapply(unique(opt$cluster), function(x)  c(median(opt[opt$cluster == x,1]), median(opt[opt$cluster == x,2]))))}
    ctr = as.data.frame(ctr)
    colnames(ctr) = c("Axis1", "Axis2")
  }
  if (type == "unimodal"){
    if (method == "max"){ ctr <- opt[which.max(opt[,3]),1:2]}
    if (method == "mean"){ ctr <- t(as.matrix(cov.wt(opt[,1:2], wt = opt[,3])$center))}
    if (method == "median"){ ctr <- t(as.matrix(c(median(opt[,1]), median(opt[,2]))))}
    ctr = as.numeric(ctr)
    names(ctr) = c("Axis1", "Axis2")
  }
  return(ctr)
}
