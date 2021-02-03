#' @title CREATE RASTER PROJECTIONS
#'
#' @description Create a stack raster map with species model predictions
#'
#' @param spsNames species names to build maps from
#' @param df predicted model data frame
#' @param ras raster template in case extent, resolution and CRS are to be set at once
#'
#' @return Raster stack
#'
#'@details Returns an error if \code{filename} does not exist.
#'
#' @examples
#' \dontrun{
#' accident_2015 <- fars_read("Project/data/accident_2015.csv.bz2")
#' }
#'
#' @importFrom raster rasterFromXYZ compareRaster extent stack
#'
#' @export
raster_projection <- function(df, spsNames = NULL, ras = NULL) {
  ### inner elements
  sp.rasdis <- list()
  if(is.null(spsNames)){spsNames = colnames(df)[-c(1:2)]}
  if(!all(spsNames %in% colnames(df))){ stop("Some species selected are not present in the data frame.")}
  ### error messages
  ## give rasters same extent and resolution
  sp.rasdis <- sapply(spsNames, function(i) rasterFromXYZ(df[, c("x", "y", i)]))
  if (is.null(ras)){
    if (compareRaster(sp.rasdis, stopiffalse = F) == F){
      NR = max(sapply(sp.rasdis, function(x) { nrow(x) }))
      NC = max(sapply(sp.rasdis, function(x) { ncol(x) }))
      ras = c(min(sapply(sp.rasdis, function(x) { extent(x)[1] })),
              max(sapply(sp.rasdis, function(x) { extent(x)[2] })),
              min(sapply(sp.rasdis, function(x) { extent(x)[3] })),
              max(sapply(sp.rasdis, function(x) { extent(x)[4] })))
      ras.template <- raster::raster(nrow=NR,ncol=NC)
      raster::extent(ras.template) <- raster::extent(ras)
      sp.rasdis <- sapply(spsNames, function(i) raster::resample(sp.rasdis[[i]], ras.template, method='ngb'))
    }
  }
  if (!is.null(ras)){
    ras.template <- raster::raster(nrow=nrow(ras),ncol=ncol(ras))
    raster::extent(ras.template) <- raster::extent(ras)
    sp.rasdis <- sapply(spsNames, function(i) raster::resample(sp.rasdis[[i]], ras.template, method='ngb'))
  }
  sp.rasdis <- stack(sp.rasdis)
  message("Stacked species:")
  print(names(sp.rasdis))
  return(sp.rasdis)
}
