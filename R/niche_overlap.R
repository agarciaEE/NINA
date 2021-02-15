#' @title Niche Overlap
#'
#'
#' @description Estimates the proportion of niche overlap between two niche objects of class NINA. The function uses the Schoener's D metric with an optional implementation that weights niche differences in base to the proximity to the niche centroid of the first argument.
#'
#' @param x NINA niche object
#' @param y NINA niche object
#' @param cor Logical whether to use environmentally corrected densites
#' @param centroid.w logical whether to weight niche overlap by distance to niche centroid
#' @param type If \code{centroid.w} is TRUE, type of centroid to estimate. Default is 'unimodal'. See \code{\link[NINA]{niche_position}}
#' @param method If \code{centroid.w} is TRUE, statistic method to estimate the niche centroid. Default is 'median'. See \code{\link[NINA]{niche_position}}
#' @param quantile If \code{centroid.w} is TRUE, quantle threshold to filter niche densities. Default is 0.5. See \code{\link[NINA]{niche_position}}
#'
#' @return Numeric value indicating the proportion of niche overlap between two niches
#'
#'
#' @importFrom raster extent resample xyFromCell as.matrix raster
#'
#' @export
niche_overlap <- function(x, y, cor = F, centroid.w = F, type = "unimodal", method = "median", quantile = 0.95) {

  if (all(class(x) != c("NINA" , "niche"))){stop("'x' is not of a niche object of class NINA")}
  if (all(class(y) != c("NINA" , "niche"))){stop("'y' is not of a niche object of class NINA")}
  if (length(x$x) != length(y$x)) {stop("niche objects have different spatial grid size")}
  R = length(x$x)

  if (raster::compareRaster(list(x$Z, y$Z), stopiffalse = F) == F){
    ras = raster::extent(x$Z)
    rasterEx <- raster::extent(ras)
    ras.template <- raster::raster(nrow=R,ncol=R)
    raster::extent(ras.template) <- rasterEx

    for (ii in c("Z", "z", "z.uncor", "z.cor", "w")) {
      y[[ii]] <- raster::resample(y[[ii]], ras.template, method = "ngb")
      y[[ii]][is.na(y[[ii]])] <- 0
    }
  }
  if (cor) {
    p1 <- t(raster::as.matrix(x$z.cor)/sum(raster::as.matrix(x$z.cor)))
    p2 <- t(raster::as.matrix(y$z.cor)/sum(raster::as.matrix(y$z.cor)))
  } else {
    p1 <- t(raster::as.matrix(x$z.uncor)/sum(raster::as.matrix(x$z.uncor)))
    p2 <- t(raster::as.matrix(y$z.uncor)/sum(raster::as.matrix(y$z.uncor)))
  }

  if (centroid.w) {
    Cp = niche_position(x,  type = type, method = method, quantile = quantile, cor = cor)
    ## Compute Euclidean distances to the niche Centroid
    pos = raster::xyFromCell(x$Z, 1:R^2)
    dist = apply(Cp[1,], 1, function(ii) unlist(sqrt((ii[1]-pos[,1])^2 + (ii[2]-pos[,2])^2)))
    if(ncol(dist)> 1) {dist = rowSums(dist)} else {dist = dist[,1]}
    dist = 1 - dist / max(dist)

    D <- 1 - (0.5 * (sum(abs(p1 - p2)*dist)))
  } else {
    D <- 1 - (0.5 * (sum(abs(p1 - p2))))
  }

  return(D)
}
