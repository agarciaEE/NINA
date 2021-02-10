#' @title ESTIMATE EC MODEL internnal function
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
#' @examples
#' \dontrun{
#' accident_2015 <- fars_read("Project/data/accident_2015.csv.bz2")
#' }
#'
#' @importFrom raster maxValue rasterize stack
#' @importFrom spatialEco raster.gaussian.smooth
#'
#' @keywords internal
#' @noRd
#'
EC_model_ <- function(en, W, R = 100, D = 0){

  if( !is.na(maxValue(en$z.uncor)) >  0){
    if(maxValue(W$z.uncor) ==  0){
      ec = en
      ec$Z = W$Z
      ec$z = W$z
      ec$w = W$w
      ec$z.uncor = W$z.uncor
      ec$z.cor = W$z.uncor
    }
    else{
      b = raster::as.data.frame(en$z.uncor, xy = T)
      w = raster::as.data.frame(W$z.uncor, xy = T)
      betas = raster::as.data.frame(W$betas, xy = T)
      b[,1:2] =  b[,1:2] * D + b[,1:2] * w[,3]
      if( D == 0){
        b = b[w[,3] != 0,]
      }
      ras = c(min(c(w[w[,3] != 0,1],b[b[,3] != 0, 1])),
              max(c(w[w[,3] != 0,1],b[b[,3] != 0, 1])),
              min(c(w[w[,3] != 0,2],b[b[,3] != 0, 2])),
              max(c(w[w[,3] != 0,2],b[b[,3] != 0, 2])))
      rasterEx <- raster::extent(ras)
      ras.template <- raster::raster(nrow=R,ncol=R)
      raster::extent(ras.template) <- rasterEx
      ec = en
      ec$glob = en$glob * raster::extract(W$z.uncor, en$glob)
      ec$glob1 = ec$glob
      ec$sp = en$sp * raster::extract(W$z.uncor, en$sp)
      ec$x = seq(ras[1], ras[2], length.out = R)
      ec$y = seq(ras[3], ras[4], length.out = R)
      if(W$alpha == 1){
        ec$z <- ec$z * D + ec$z * W$z.uncor
        ec$Z = raster::resample(W$z, ras.template, method = "ngb")
        ec$z = raster::resample(ec$z, ras.template, method = "ngb")
        ec$z.uncor = ec$z / raster::cellStats(ec$z, "max")
        ec$betas = raster::resample(round(W$betas), ras.template, method = "ngb")
      }
      else{
        ec$Z = raster::rasterize(w[,1:2]*w[,3], ras.template, field = w[,3], fun = mean)
        ec$Z= spatialEco::raster.gaussian.smooth(ec$Z, n = 3,type = mean)
        ec$Z = ec$Z * raster::cellStats(W$z, "max")
        ec$z.uncor = raster::rasterize(b[,1:2], ras.template, field = b[,3], fun = mean)
        ec$z.uncor= spatialEco::raster.gaussian.smooth(ec$z.uncor, n = 3,type = mean)
        ec$z.uncor = ec$z.uncor / raster::cellStats(ec$z.uncor, "max")
        ec$betas = raster::stack(sapply(names(W$betas), function(i) raster::rasterize(betas[,1:2]*w[,3], ras.template, field = betas[,i], fun = mean)))
      }
      ec$Z[is.na(ec$Z)] <- 0
      ec$z.uncor[is.na(ec$z.uncor)] <- 0
      ec$z = ec$z.uncor * raster::cellStats(en$z, "max")
      ec$w <- ec$z.uncor
      ec$w[ec$w > 0] <- 1
      ec$z.cor <- ec$z/ec$Z
      ec$z.cor[is.na(ec$z.cor)] <- 0
      ec$z.cor <- ec$z.cor/raster::cellStats(ec$z.cor, "max")
    }
    message("\t...Success!")
  }
  else{
    stop("input model has not positive predictions")
  }
  return(ec)
}


