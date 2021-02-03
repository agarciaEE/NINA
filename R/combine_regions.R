#' @title COMBIINE REGIONAL MODELS FUNCTION
#'
#' @param z.list List of regional niche models
#'
#' @param env.scores Environmental PCA scores
#' @param w.global template for global niche space
#' @param R niche grid
#'
#' @description Combine regional niche models
#'
#' @return Niche model
#'
#'@details Returns an error if \code{filename} does not exist.
#'
#' @examples
#' \dontrun{
#' accident_2015 <- fars_read("Project/data/accident_2015.csv.bz2")
#' }
#'
#' @importFrom stats na.exclude
#' @importFrom raster compareRaster
#' @importFrom ecospat ecospat.grid.clim.dyn
#'
#' @export
combine_regions <- function(z.list,  env.scores, w.global,  R = 100){
  z = reverse_list(z.list)
  z.comb= list()
  for (s in names(z)){
    if (length(z[[s]]) > 1){
      glob <- na.exclude(do.call(rbind, lapply(z[[s]], function(X) X$glob)))
      sp <- na.exclude(do.call(rbind, lapply(z[[s]], function(X) X$sp)))
      z.comb[[s]] <- ecospat.grid.clim.dyn(glob, glob, sp, R)
    }
    else {
      glob <- na.exclude(as.data.frame(lapply(z[[s]], function(X) X$glob)))
      sp <- na.exclude(as.data.frame(lapply(z[[s]], function(X) X$sp)))
      z.comb[[s]] <- ecospat.grid.clim.dyn(glob, glob, sp, R)
    }
    if(!missing(env.scores)){
      ras = c(min(env.scores[,3]),
              max(env.scores[,3]),
              min(env.scores[,4]),
              max(env.scores[,4]))
      rasterEx <- raster::extent(ras)
      ras.template <- raster::raster(nrow=R,ncol=R)
      raster::extent(ras.template) <- rasterEx
      z.comb[[s]]$x <- seq(ras[1], ras[2], length.out = R)
      z.comb[[s]]$y <- seq(ras[3], ras[4], length.out = R)
      z.comb[[s]]$Z = raster::resample(z.comb[[s]]$Z, ras.template, method = "ngb")
      z.comb[[s]]$z = raster::resample(z.comb[[s]]$z, ras.template, method = "ngb")
      z.comb[[s]]$Z[is.na(z.comb[[s]]$Z)] = 0
      z.comb[[s]]$z[is.na(z.comb[[s]]$z)] = 0
      z.comb[[s]]$z.uncor = z.comb[[s]]$z/cellStats(z.comb[[s]]$z, "max")
      z.comb[[s]]$w <- z.comb[[s]]$z.uncor
      z.comb[[s]]$w[z.comb[[s]]$w > 0] <- 1
      z.comb[[s]]$z.cor <- z.comb[[s]]$z/z.comb[[s]]$Z
      z.comb[[s]]$z.cor[is.na(z.comb[[s]]$z.cor)] <- 0
      z.comb[[s]]$z.cor <- z.comb[[s]]$z.cor/cellStats(z.comb[[s]]$z.cor, "max")
    }
    if(!missing(w.global)){
      if (compareRaster(z.comb[[s]]$Z, w.global[[s]]$Z, stopiffalse = F) == F){
        stop("rasters have different extent")
      }
      z.comb[[s]]$z.uncor =  z.comb[[s]]$z.uncor * w.global[[s]]$z.uncor
      z.comb[[s]]$z =  z.comb[[s]]$z.uncor * cellStats( z.comb[[s]]$z, "max")
      z.comb[[s]]$w <- z.comb[[s]]$z.uncor
      z.comb[[s]]$w[z.comb[[s]]$w > 0] <- 1
      z.comb[[s]]$z.cor <- z.comb[[s]]$z/z.comb[[s]]$Z
      z.comb[[s]]$z.cor[is.na(z.comb[[s]]$z.cor)] <- 0
      z.comb[[s]]$z.cor <- z.comb[[s]]$z.cor/cellStats(z.comb[[s]]$z.cor, "max")
    }
  }
  return(z.comb)
}
