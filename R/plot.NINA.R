#' @title Plot
#'
#' @description  Plots the summary of the output models of an object class NINA
#'
#' @param x An object class NINA
#' @param ... Additional arguments for S3 methods
#'
#' @return A \code{NINA} class object description
#'
#' @method plot NINA
#'
#' @examples
#' \dontrun{
#' EN = EN_model(env_data, occ_data1, cluster = "env")
#' plot(EN)
#' }
#'
#' @importFrom ecospat ecospat.plot.contrib
#' @importFrom plotrix addtable2plot
#' @importFrom raster rasterize maxValue
#' @importFrom sp SpatialPolygons Polygons Polygon SpatialPoints
#' @importFrom stats cov.wt
#' @importFrom car dataEllipse
#'
#' @export
plot.NINA <- function(x, ...){

  df = merge(x$env.scores, x$obs, by = c(1,2), all = T)
  pca = x$pca
  mode.region = if(!is.null(x$clus)){TRUE} else {FALSE}

  ellipse.env <- car::dataEllipse(df[,3], df[,4], levels = 0.9, xlab = "PC1", ylab = "PC2", draw = F)
  center.env <- stats::cov.wt(df[,3:4])$center
  ellipse.env = sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(list(sp::SpatialPoints(ellipse.env)))), 1)))
  ellipse.occ <- car::dataEllipse(df[!is.na(df[,5]),3], df[!is.na(df[,5]),4], levels = 0.9, xlab = "PC1", ylab = "PC2", draw = F)
  center.occ <- stats::cov.wt(df[!is.na(df[,5]),3:4])$center
  ellipse.occ = sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(list(sp::SpatialPoints(ellipse.occ)))), 1)))

  if (mode.region){
    layout(matrix(c(1,1,2,2,
                    1,1,3,4,
                    5,6,7,8), 3, 4, byrow = T))
    clus = raster::rasterize(x$clus[,1:2], x$maps[[1]], field = as.numeric(x$clus[,3]), fun = "last", na.rm = T)
    clusNames =levels(x$clus[,3])
    n.clus = length(clusNames)
    tab = x$tab
  } else{
    layout(matrix(c(1,1,2,2,
                    1,1,2,2,
                    3,4,5,6), 3, 4, byrow = T))
  }
  plot(df[,3:4], col = "green2", pch = 19)
  plot(ellipse.env, add = T,  border = "green4", lty = 2, lwd = 2)
  legend(center.env[1], center.env[2], "environment",
         xjust = 0.5,      # 0.5 means center adjusted
         yjust = 0.5,      # 0.5 means center adjusted
         x.intersp = -0.5, # adjust character interspacing as you like to effect box width
         y.intersp = 0.1)
  points(df[,3:4], col = rep("blue")[df[,5]], pch = 4)
  plot(ellipse.occ, add = T,  border = "blue4", lty = 2, lwd = 2)
  legend(center.occ[1], center.occ[2], "ocurrences",
         xjust = 0.5,      # 0.5 means center adjusted
         yjust = 0.5,      # 0.5 means center adjusted
         x.intersp = -0.5, # adjust character interspacing as you like to effect box width
         y.intersp = 0.1)
  ecospat::ecospat.plot.contrib(pca$co, pca$eig)
  if (mode.region){
    plot(clus, col = viridis::viridis(n.clus), legend = F)
    plot(clus, legend.only = T, breaks = seq(0.5,raster::maxValue(clus)+0.5,1), col = viridis::viridis(n.clus),
         axis.args=list(at = 1:raster::maxValue(clus),labels=clusNames))
    plot.new()
    plotrix::addtable2plot(0,0,tab,bty="o",display.rownames=T,hlines=F, cex=1.5)
  }
  if (raster::nlayers(x$maps) <= 4) {
    for (i in 1:raster::nlayers(x$maps)) plot(x$maps[[i]], main  = names(x$maps[[i]]))
  }
  if (raster::nlayers(x$maps) > 4) {
    for (i in 1:4) plot(x$maps[[i]], main  = names(x$maps[[i]]))
    warning(paste("Ploting only the first four species maps of a total of", raster::nlayers(x$maps)), immediate. = T)
  }
  par(mfrow=c(1,1))
}
