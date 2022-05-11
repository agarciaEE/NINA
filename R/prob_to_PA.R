#' @title Probabilities to Presence/Absence
#'
#' @param x raster or data frame object
#' @param th selected threshold to convert to P/A
#' @param output type of output objct wanted. data frame or raster object
#' @param ras User can provide a raster object as template to match resolution, grid size and projection
#'
#' @description Converts probability distributions into Presence/Absence
#'
#' @return Data frame or raster
#'
#' @details Returns an error if \code{x} does not exist.
#'
#' @importFrom raster as.data.frame
#' @export
prob_to_PA <- function(x, th = NULL, output = c("data.frame", "raster"), ras = NULL){

  output = output[1]
  if (class(x) %in% c("raster", "RasterBrick", "RasterStack")){
    df <- raster::as.data.frame(x, xy =T)
  }
  if(is.data.frame(x)){
    df <- x
    spsNames <- colnames(df)[-c(1:2)]
    if (is.null(th)){
      if (length(spsNames) == 1){
        th <- mean(df[,spsNames], na.rm = T)
        names(th) = spsNames
      }
      else {
        th <- apply(df[,spsNames], 2, function(i) mean(i, na.rm = T))
      }
    }
  }
  df[is.na(df)] = 0
  for (i in 1:length(spsNames)) {
    df[df[,i+2] < th[i],i+2] = 0
    df[df[,i+2] >= th[i],i+2] = 1
  }
  if(output == "raster"){
    df <- raster_projection(df, ras = ras)
  }
  return(df)
}
