#' @title ESTIMATE EC MODEL internal function
#'
#' @param en NINA Niche model to transform the niche space
#' @param W Weighting coefficients
#' @param R Niche space grid
#' @param D Numeric value for independence of interactions
#'
#' @description Transform environmental niche space into ecological niche space
#'
#' @return Data frame.
#'
#'@details Returns an error if \code{filename} does not exist.
#'
#'
#' @importFrom raster maxValue rasterize stack
#' @importFrom spatialEco raster.gaussian.smooth
#'
#' @keywords internal
#' @noRd
#'
EC_model_ <- function(en, W, R = 100, D = 0){

  if(is.na(raster::maxValue(W$z.uncor)) || raster::maxValue(W$z.uncor) ==  0
     || is.na(raster::maxValue(en$z.uncor)) || raster::maxValue(en$z.uncor) ==  0){

    ec = en
    g = W
    message("\t...Empty niche.")
  }
  else{

    b = raster::as.data.frame(en$z, xy = T)
    w = raster::as.data.frame(W$z.uncor, xy = T)
    betas = raster::as.data.frame(W$betas, xy = T)
    D = 0
    b[,1:2] =  b[,1:2] * D + b[,1:2] * w[,3]
    w[,1:2] = w[,1:2] * w[,3]
    if( D == 0){
      b = b[w[,3] != 0,]
    }
    ras = c(min(c(w[w[,3] != 0,1],b[b[,3] != 0, 1]), na.rm = T),
            max(c(w[w[,3] != 0,1],b[b[,3] != 0, 1]), na.rm = T),
            min(c(w[w[,3] != 0,2],b[b[,3] != 0, 2]), na.rm = T),
            max(c(w[w[,3] != 0,2],b[b[,3] != 0, 2]), na.rm = T))
    rasterEx <- raster::extent(ras)
    ras.template <- raster::raster(nrow=R,ncol=R)
    raster::extent(ras.template) <- rasterEx
    ec = en
    g = W
    g$x = seq(ras[1], ras[2], length.out = R)
    g$y = seq(ras[3], ras[4], length.out = R)
    g$betas = raster::resample(g$betas, ras.template, method = "ngb")

    ec$glob = en$glob * D + en$glob * raster::extract(W$z.uncor, en$glob)
    ec$glob1 = en$glob1 * D + en$glob1 * raster::extract(W$z.uncor, en$glob1)
    ec$sp = en$sp* D + en$sp * raster::extract(W$z.uncor, en$sp)
    ec$x = seq(ras[1], ras[2], length.out = R)
    ec$y = seq(ras[3], ras[4], length.out = R)

    if(W$alpha == 1){
      g$Z = raster::resample(g$Z, ras.template, method = "ngb")
      g$z = raster::resample(g$z, ras.template, method = "ngb")

      ec$Z = g$z
      ec$z = raster::resample(ec$z, ras.template, method = "ngb")
    }
    else{
      g$Z = raster::resample(g$Z, ras.template, method = "ngb")
      g$z = raster::rasterize(w[,1:2], ras.template, field = w[,3], fun = mean)
      maxV = raster::maxValue(g$z)
      g$z = spatialEco::raster.gaussian.smooth(g$z, n = 5,type = mean)
      g$z = g$z * maxV / raster::maxValue(g$z)

      ec$Z = g$z * raster::cellStats(W$z, "max")
      ec$z = raster::rasterize(b[,1:2], ras.template, field = b[,3], fun = mean)
      maxV = raster::maxValue(ec$z)
      ec$z = spatialEco::raster.gaussian.smooth(ec$z, n = 5,type = mean)
      ec$z = ec$z * maxV / raster::maxValue(ec$z)
    }
    g$z[is.na(g$z)] <- 0
    g$Z[is.na(g$Z)] <- 0
    g$z.uncor = g$z * raster::cellStats(en$z, "max")
    g$z.uncor[is.na(g$z.uncor)] <- 0
    g$w <- g$z.uncor
    g$w[g$w > 0] <- 1
    g$z.cor <- g$z/g$Z
    g$z.cor[is.na(g$z.cor)] <- 0
    g$z.cor <- g$z.cor/raster::cellStats(g$z.cor, "max")

    ec$z[is.na(ec$z)] <- 0
    ec$Z[is.na(ec$Z)] <- 0
    ec$z.uncor = ec$z * raster::cellStats(en$z, "max")
    ec$z.uncor[is.na(ec$z.uncor)] <- 0
    ec$w <- ec$z.uncor
    ec$w[ec$w > 0] <- 1
    ec$z.cor <- ec$z/ec$Z
    ec$z.cor[is.na(ec$z.cor)] <- 0
    ec$z.cor <- ec$z.cor/raster::cellStats(ec$z.cor, "max")

    message("\t...Success!")
  }
  return(list(ec = ec, g = g))
}


