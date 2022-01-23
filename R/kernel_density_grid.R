#' @title Kernel density grid internal FUNCTION
#'
#' @description Estimates a density grid using epanechnikov or bivnorm method
#'
#' @param x 2-D coordinates matrix
#' @param h smoothing parameter for the kernel estimation. Default is 'href'. Altenrtanively can be set to 'LSCV' or any given numeric value
#' @param R  the size of the square raster map to be created (number of rows/columns)
#' @param ext extent of the raster to be created
#' @param env.mask raster mask to resample the created kernel densiity grid raster
#' @param th numeric threshold to filter density values
#' @param method "epa" or "bivnorm"
#'
#' @return raster object
#'
#' @details
#'
#' @examples
#' \dontrun{
#' sp <- cbind(runif(10, -1, 5), runif(10, -2, 10))
#' kdg<-kernel_density_grid(sp, ext = c(-1,5,-2,10))
#' plot(kdg)
#' }
#'
#' @importFrom adehabitatMA ascgen
#' @importFrom sp SpatialPoints
#' @importFrom adehabitatHR kernelUD
#' @importFrom raster raster extract resample compareRaster
#' @importFrom stats quantile
#'
#' @keywords internal
#' @noRd
#'
kernel_density_grid <- function(x,  h = "href", R = 100, ext = NULL, env.mask = NULL,
                                th = NULL,  method = c("epa", "bivnorm")){

  method = method[1]
  if (is.null(ext)){
    ext <- c(range(x[,1]), range(x[,2]))
  }
  xr <- data.frame(cbind((x[, 1] - ext[1])/abs(ext[2] -
                                                 ext[1]), (x[, 2] - ext[3])/abs(ext[4] - ext[3])))
  mask <- adehabitatMA::ascgen(sp::SpatialPoints(cbind((0:(R))/R,
                                                       (0:(R)/R))), nrcol = R - 2, count = FALSE)
  x.dens <- adehabitatHR::kernelUD(sp::SpatialPoints(xr[,
                                                        1:2]), h = h, grid = mask, kern = method)
  x.dens <- raster::raster(xmn = ext[1], xmx = ext[2],
                           ymn = ext[3], ymx = ext[4], matrix(x.dens$ud,
                                                              nrow = R))
  if (!is.null(th)) {
    th.value <- stats::quantile(raster::extract(x.dens,
                                         x), th)
    x.dens[x.dens < th.value] <- 0
  }
  if (!is.null(env.mask)) {
    if(raster::compareRaster(c(x.dens, env.mask), stopiffalse = F) == F){
      x.dens <- raster::resample(x.dens, env.mask, method = "ngb")
    }
    x.dens <- x.dens * env.mask
  }
  return(x.dens)
}
