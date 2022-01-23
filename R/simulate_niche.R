#' @title Simulate niche
#'
#' @param z m by n matrix indicating the association coefficient (-1 to 1). m are species to be modeled as rows and n interactions as columns
#' @param glob n by n matrix indicating the competition coefficient between interactions (0 to 1).
#' @param glob1 NINA EN model object for species group one
#' @param method NINA EN model object for species group two
#' @param sample.size String indicating whether to perform at a region or a global level. Note that if models have not been estimated at a region level and it is selected it will produce an error
#' @param R Numeric value indicating independence from biotic associations. Value must be comprised between 0 and 1.
#' @param use.whole.env Logical. Sample from the entire environmental space
#'
#' @description Simulates a niche object of class NINA from a former niche or a matrix of environmental coordinates. The function uses ecospat.grid.clim.dyn  from ecospat R package to generate the niche object
#'
#' @return niche object class NINA
#'
#' @details It needs either z or glob arguments to run, otherwise will throw an error
#'
#' @examples
#' \dontrun{
#' z <- simulate_niche(glob = cbind(rnorm(100), rnorm(100)))
#' }
#'
#' @importFrom ecospat ecospat.grid.clim.dyn
#' @importFrom raster extract cellStats res xyFromCell Which values extent
#'
#' @export
simulate_niche <- function(z, use.whole.env = F, method = 1, glob = NULL, glob1 = glob, sample.size = nrow(glob)*0.1, R = 100){

  if (missing(z)){
    if (!is.null(glob)){
      sp <- glob[sample(1:nrow(glob), size = sample.size),]
      z <- ecospat::ecospat.grid.clim.dyn(glob, glob1, sp, R = R)
    }
  } else {
    R = length(z$x)
    x.inc <- raster::res(z$Z)[1] # distance between blocks in x axis
    y.inc <- raster::res(z$Z)[2] # distance between blocks in y axis

    centroid.z <-niche_position(z)
    if(nrow(centroid.z) > 1){ # if more than one, choose one randomly
      centroid.z <- centroid.z[sample(1:nrow(centroid.z),1),]
    }
    if (use.whole.env){
      z.cells <- raster::Which(z$Z > 0, cells = T) # get available cells indexes
      weights.z <- raster::values(z$Z) # environmental weights
      weights.z <- weights.z[weights.z>0]
    } else {
      z.cells <- raster::Which(z$z > 0, cells = T) # get available cells indexes
      weights.z <- raster::values(z$z) # environmental weights
      weights.z <- weights.z[weights.z>0]
    }
    z.uncor <- z$z.uncor
    rand.centroid.z <- raster::xyFromCell(z$Z, sample(z.cells, size = 1, replace = FALSE, prob = weights.z))
    if (length(raster::Which(z$z.uncor > 0, cells = T)) > 5){
        if (method == 1){
        ## shift niche coordinates
        x <- z$x
        y <- z$y
        xshift.z <- as.numeric(trunc((rand.centroid.z[1] - centroid.z[1]) / x.inc))
        if (xshift.z != 0){
          if (xshift.z < 0){
            for (i in 1:xshift.z){
              x <- c(x[1]-x.inc, x[-length(x)])
            }
          } else {
            for (i in 1:xshift.z){
              x <- c(x[-1], rev(x)[1]+x.inc)
            }
          }
        }
        yshift.z <- as.numeric(trunc((rand.centroid.z[2] - centroid.z[2]) / y.inc))
        if (yshift.z != 0){
          if (yshift.z < 0){
            for (i in 1:yshift.z){
              y <- c(y[1]-y.inc, y[-length(y)])
            }
          } else {
            for (i in 1:yshift.z){
              y <- c(y[-1], rev(y)[1]+y.inc)
            }
          }
        }
        if (!is.null(intersect(extent(z.uncor), extent(z$Z)))){
          # re-adjust niche
          z.uncor <- raster::resample(z.uncor, z$Z, method = "ngb")
          # remove responses outside the environmental space
          Z <- z$Z
          Z[Z>0] = 1
          z.uncor <- z.uncor  * Z
          # get a sp shifted scores
          s.cells <- raster::Which(z.uncor > 0, cells = T)
          vals <- raster::values(z.uncor)[s.cells]
          vals[is.na(vals)] = 0
          vals <- vals[vals>0]
        }
      }
      if (nrow(z$sp) <= length(s.cells)){
        s.cells <- sample(s.cells, size = nrow(z$sp), replace = FALSE, prob = vals)
      }
      #s.cells <- s.cells[!duplicated(s.cells)]
      sp <- raster::xyFromCell(z$Z, s.cells)
      if(method == 2){
        sp <- cbind(Axis1 = z$sp[,1]+(rand.centroid.z[1] - centroid.z[1]),
                    Axis2 = z$sp[,2]+(rand.centroid.z[2] - centroid.z[2]))
        v <- raster::extract(z$w, sp)
        sp <- sp[which(v == 1),]
      }
      # estimate shifted niche
      z <- ecospat::ecospat.grid.clim.dyn(z$glob, z$glob1, sp, R = R)
    }
    else{
      s.cells <- raster::Which(z$z.uncor > 0, cells = T)
      sp <- raster::xyFromCell(z$z, s.cells)
      v <- raster::extract(z$z.uncor, sp)
      sp <- cbind(Axis1 = sp[,1]+(rand.centroid.z[1] - centroid.z[1]),
                  Axis2 = sp[,2]+(rand.centroid.z[2] - centroid.z[2]))
      Z <- z$Z
      Z[Z>0] = 1
      z$z.uncor <- rasterize(sp, Z, field = v, FUN = mean)
      z$z <- z$z.uncor * cellStats(z$z, "max")
      z$w <- z$z.uncor
      z$w[z$w > 0] = 1
      z$z.cor = z$z/z$Z
      z$z.cor[is.na(z$z.cor)] = 0
    }
  }
  class(z) <- c("NINA", "niche")
  return(z)
}
