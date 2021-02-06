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
#' @importFrom raster rasterFromXYZ compareRaster extent stack crs
#'
#' @export
raster_projection <- function(df, spsNames = NULL, ras = NULL, crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") {
  ### inner elements
  sp.rasdis <- list()
  if(is.null(spsNames)){spsNames = colnames(df)[-c(1:2)]}
  if(!all(spsNames %in% colnames(df))){ stop("Some species selected are not present in the data frame.")}
  ### error messages
  if (!is.null(ras)){
    sp.rasdis <- sapply(spsNames, function(i) raster::rasterize(df[,1:2], ras, field = df[,i], fun = mean, na.rm = T, update = T))
  } else{
    sp.rasdis <- raster::rasterFromXYZ(df, crs = crs)
  }
  sp.rasdis <- stack(sp.rasdis)

  return(sp.rasdis)
}
