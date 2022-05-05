#' @title  SAVE SPECIES NICHE
#'
#' @param filepath path to the file
#'
#' @description reads a species NINA niche into from disk
#'
#'
#' @examples
#' \dontrun{
#' EN1 <- EN_model(env_data, occ_data1, cluster = "env", n.clus = 5)
#' write_niche(EN1$z.mod$A$spB, path = "./", file = "ENtest")
#' read_niche(path = "./", file = "ENtest.txt")
#' }
#'
#' @importFrom raster raster extent
#' @importFrom utils read.table
#'
#' @export
read_niche <- function(filepath){

  x <- utils::read.table(filepath, header = F)

  l <- list()
  elements_idx <- which(grepl("\\/[a-z]*", x$V1))
  elements_nms <- gsub("\\/", "", x$V1[grepl("\\/[a-z]*", x$V1)])
i = 5
  for (i in 1:length(elements_idx)){
    idx = elements_idx[i]
    start = as.integer(idx+1)
    end = as.integer(elements_idx[i+1]-1)
    if (is.na(end)){end = nrow(x)}
    l[[elements_nms[i]]] <- as.numeric(x[start:end,])
  }
  for (d in elements_nms[3:5]){
    l[[d]] <- cbind.data.frame(split(l[[d]], rep(c("Axis1", "Axis2"),
                                                 times = length(l[[d]])/2)),stringsAsFactors=F)
  }
  ext <- raster::extent(c(range(l$x), range(l$y)))
  nr = length(l$y)
  nc = length(l$x)
  for (m in elements_nms[6:10]){
    r <- raster::raster(t(matrix(l[[m]], nc, nr)))
    raster::extent(r)  <-  ext
    l[[m]] <- r
  }

  return(l)
}

