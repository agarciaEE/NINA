#' @title COMBIINE REGIONAL MODELS FUNCTION
#'
#' @param z.list List of regional niche models
#'
#' @param env.scores Environmental PCA scores
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
combine_regions <- function(z.list, env.scores,  R = 100){

  z.comb= list()
  if(!missing(env.scores)){
    ras = c(min(env.scores[,3]),
            max(env.scores[,3]),
            min(env.scores[,4]),
            max(env.scores[,4]))
    rasterEx <- raster::extent(ras)
    ras.template <- raster::raster(nrow=R,ncol=R)
    raster::extent(ras.template) <- rasterEx
  } else {
    ras = c(min(sapply(z.list, function(x) { raster::extent(x[[1]]$Z)[1] })),
            max(sapply(z.list, function(x) { raster::extent(x[[1]]$Z)[2] })),
            min(sapply(z.list, function(x) { raster::extent(x[[1]]$Z)[3] })),
            max(sapply(z.list, function(x) { raster::extent(x[[1]]$Z)[4] })))
    rasterEx <- raster::extent(ras)
    ras.template <- raster::raster(nrow=R,ncol=R)
    raster::extent(ras.template) <- rasterEx
  }

  z = reverse_list(z.list)
  for (s in names(z)){
    if (length(z[[s]]) > 1){
      z.comb[[s]]$glob <- na.exclude(do.call(rbind, lapply(z[[s]], function(X) X$glob)))
      z.comb[[s]]$glob1 = na.exclude(do.call(rbind, lapply(z[[s]], function(X) X$glob1)))
      z.comb[[s]]$sp <- na.exclude(do.call(rbind, lapply(z[[s]], function(X) X$sp)))
#      z.comb[[s]] <- ecospat.grid.clim.dyn(glob, glob, sp, R)
    }
    else {
      z.comb[[s]]$glob <- na.exclude(as.data.frame(lapply(z[[s]], function(X) X$glob)))
      z.comb[[s]]$glob1 = na.exclude(as.data.frame(lapply(z[[s]], function(X) X$glob1)))
      z.comb[[s]]$sp <- na.exclude(as.data.frame(lapply(z[[s]], function(X) X$sp)))
#      z.comb[[s]] <- ecospat.grid.clim.dyn(glob, glob, sp, R)
    }

    z.comb[[s]]$x <- seq(ras[1], ras[2], length.out = R)
    z.comb[[s]]$y <- seq(ras[3], ras[4], length.out = R)
    z.comb[[s]]$Z = sum(raster::stack(sapply(z[[s]], function(i) raster::resample(i$Z, ras.template, method = "ngb"))))
    z.comb[[s]]$z = sum(raster::stack(sapply(z[[s]], function(i) raster::resample(i$z, ras.template, method = "ngb"))))
    z.comb[[s]]$Z[is.na(z.comb[[s]]$Z)] = 0
    z.comb[[s]]$z[is.na(z.comb[[s]]$z)] = 0
    z.comb[[s]]$z.uncor = z.comb[[s]]$z/cellStats(z.comb[[s]]$z, "max")
    z.comb[[s]]$w <- z.comb[[s]]$z.uncor
    z.comb[[s]]$w[z.comb[[s]]$w > 0] <- 1
    z.comb[[s]]$z.cor <- z.comb[[s]]$z/z.comb[[s]]$Z
    z.comb[[s]]$z.cor[is.na(z.comb[[s]]$z.cor)] <- 0
    z.comb[[s]]$z.cor <- z.comb[[s]]$z.cor/cellStats(z.comb[[s]]$z.cor, "max")

    if(!is.null(z[[s]][[1]][["betas"]])){
      betas = sapply(z[[s]], function(i) raster::unstack(raster::resample(i$betas, ras.template, method = "ngb")))
      names.betas = sapply(z[[s]], function(i) names(i$betas))
      for (i in 1:length(betas)){
        names(betas[[i]]) = names.betas[[i]]
      }
      betas = reverse_list(betas)
      betas = sapply(betas, function(i) sum(raster::stack(i), na.rm = T))
      maxV = sapply(betas, function(i) raster::maxValue(i))
      betas = sapply(betas, function(i) spatialEco::raster.gaussian.smooth(i[[1]], n = 5, type = mean))
      betas = sapply(1:raster::nlayers(betas), function(i) betas[[i]] * maxV[i] / raster::maxValue(betas[[i]]))

      z.comb[[s]]$betas = raster::stack(betas)
      z.comb[[s]]$alpha = mean(sapply(z[[s]], function(i) i$alpha), na.rm = T)
    }
  }
  return(z.comb)
}
