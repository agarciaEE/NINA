#' @title ESTIMATE ENVIRONMENTAL ACCESSIBILITY
#'
#' @param y.list list of niche models
#'
#' @param id species id
#' @param A.matrix Association matrix
#' @param C.matrix Competition matrix
#'
#' @description Estimates the omega
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
#' @importFrom raster stack
#'
#' @export
estimate_w<- function(y.list, id,  A.matrix = NULL, C.matrix = NULL){

  if (is.null(A.matrix)){
    A.matrix =  matrix(1, nrow = 1, ncol = length(y.list))
    rownames(A.matrix) = id ; colnames(A.matrix) = names(y.list)
  }
  betas <- estimate_betas(y.list, C.matrix)
  Xvar = colnames(A.matrix)[A.matrix[id,] != 0]
  Xvar = Xvar[Xvar %in% names(y.list)]
  w <- y.list[[1]]
  w$sp <- do.call(rbind, lapply(y.list, function(i) i$sp))
  if(length(y.list) > 1){
    w$z.uncor = sum(stack(lapply(names(y.list), function(i) betas[[i]]*as.numeric(A.matrix[id,i]))), na.rm = T)
    w$z = sum(stack(lapply(y.list, function(i) i$z)), na.rm = T)
    w$z = w$z.uncor * cellStats(w$z, "max")
  }
  else{
    w$z.uncor = betas[[1]]*A.matrix[id,names(y.list)]
    w$z = w$z.uncor * cellStats(y.list[[1]]$z, "max")
  }
  w$z.uncor[is.na(w$z.uncor)] <- 0
  w$w <- w$z.uncor
  w$w[w$w > 0] <- 1
  w$z.cor <- w$z/w$Z
  w$z.cor[is.na(w$z.cor)] <- 0
  w$z.cor <- w$z.cor/cellStats(w$z.cor, "max")
  w$betas <- stack(betas)
  w$alpha <- length(betas[Xvar])/length(y.list)
  return(w)
}
