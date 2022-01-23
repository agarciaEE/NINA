#' @title Spatial Abundance Density Index
#'
#' @description Estimates the occurrence density index of a set of occurrences over a 2-D gridded space
#'
#' @param S 2-D coordinates matrix
#' @param Sc smoothing parameter for the kernel estimation. Default is 'href'. Altenrtanively can be set to 'LSCV' or any given numeric value
#' @param So  the size of the square raster map to be created i.e. number of rows/columns
#' @param R extent of the raster to be created
#' @param h smoothing parameter for the kernel estimation. Default is 'href'. Altenrtanively can be set to 'LSCV' or any given numeric value
#' @param mask raster mask to resample the created kernel densiity grid raster
#' @param th.o numeric threshold to filter density values of occurrences
#' @param th.s numeric threshold to filter density values of environment
#' @param method "epanechnikov" or "bivnorm"
#'
#' @return list object class niche
#'
#' @details
#'
#' @examples
#' \dontrun{
#' env <- cbind(runif(1000, -1, 5), runif(1000, -2, 10))
#' sp <- cbind(runif(10, -1, 5), runif(10, -2, 10))
#' niche<- estimate_niche(sp, ext = c(-1,5,-2,10))
#' ecospat::ecospat.plot.niiche(kdg)
#' }
#'
#' @importFrom adehabitatMA ascgen
#' @importFrom sp SpatialPoints
#' @importFrom adehabitatHR kernelUD
#' @importFrom raster raster extract resample compareRaster
#' @importFrom stats quantile
#'
#'@export
estimate_niche <- function(glob, glob1, sp, R,  h = "href", mask = NULL,
                           th.o = NULL, th.s = NULL, method = c("epa", "bivnorm")) {

  method <- method[1]
  glob <- as.matrix(glob)
  glob1 <- as.matrix(glob1)
  sp <- as.matrix(sp)

  l <- list()
  if (any(sapply(list(glob, glob1, sp), ncol) > 2)) {
    warning("coordinates matrix have more than 2 dimensions...only the first 2 will be used!")
  }
  else if (any(sapply(list(glob, glob1, sp), ncol) > 2)) {
    stop("2 dimensional matrix required")
  }
  else {
    ext <- c(range(glob[,1]), range(glob[,2]))
    glob1.dens <- NINA:::kernel_density_grid(glob1, R = R, h = h, method = method, th = th.s, env.mask = mask, ext = ext)
    sp.dens <-  NINA:::kernel_density_grid(sp, R = R, h = h, method = method, th = th.o, env.mask = mask, ext = ext)
    l$x <- seq(from = ext[1], to = ext[2], length.out = R)
    l$y <- seq(from = ext[3], to = ext[4], length.out = R)
    l$glob <- glob
    l$glob1 <- glob1
    l$sp <- sp
    l$z <- sp.dens * nrow(sp)/raster::cellStats(sp.dens, "sum")
    l$Z <- glob1.dens * nrow(glob1)/raster::cellStats(glob1.dens,  "sum")
    l$z[l$Z == 0] = 0

    l$z.uncor <- l$z/raster::cellStats(l$z, "max")
    l$z.cor <- l$z/l$Z
    l$z.cor[is.na(l$z.cor)] <- 0
    l$z.cor <- l$z.cor/raster::cellStats(l$z.cor, "max")
    l$w <- l$z.uncor
    l$w[l$w > 0] <- 1
  }
  return(l)
}
