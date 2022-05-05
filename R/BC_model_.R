#' @title  ESTIMATE BC MODEL internal function
#'
#' @param z environmental niche model
#'
#' @param y.list list of environmental niche models of species interactors
#' @param id species name to be estimated. id must match the rownames in A.matrix
#' @param D Numeric indicating species independence on interactions. 0 for fully independent, 1 for species that cannot be present without the peresence of their interactors
#' @param A.matrix m by n matrix indicating the association coefficient (-1 to 1). m are species to be modelled as rows and n interactors as columns
#' @param C.matrix n by n matrix indicating the competition coefficient between interactors (0 to 1).
#' @param method Method; abundances or composition
#' @param cor Logical
#' @param K = Carrying capacity of each environmental cell
#'
#' @description Transforms environmental niche in base to species interactions
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
#' @importFrom plyr ldply
#' @importFrom tidyr spread
#' @importFrom raster cellStats
#'
#' @keywords internal
#' @noRd
#'
BC_model_ <- function(z, y.list, id, D = 1, A.matrix = NULL,
                      method = c("composition", "densities"),
                      cor = F,  K  = NULL, C.matrix = NULL ){

  out = list()
  Xvar = colnames(A.matrix)[A.matrix[id, ]  != 0]
  Xvar <- Xvar[Xvar %in% names(y.list)]
  if(length(Xvar) > 0){
    w <- estimate_w(y.list, id = id, A.matrix = A.matrix, cor = cor,
                        K = K, method = method, C.matrix = C.matrix)
    class(w) <- c("NINA", "niche")
    if (!is.null(w)){
      if (cor) {
        wc <- w$z.cor
      } else {
        wc <- w$z.uncor
      }
    } else{ wc = 0}
  }
  else{
    w = NULL
    wc = 0
  }
  if (method == "composition"){
    z$z <- z$z - z$z * D * (1-wc)
  }
  if (method == "densities"){
    if (D > 0){
      z$z <- z$z * wc / D
    }
  }
  z$z.uncor <- z$z / raster::cellStats(z$z, "max")
  z$z.uncor[is.na(z$z.uncor)] <- 0
  z$w <- z$z.uncor
  z$w[z$w > 0] <- 1
  z$z.cor <- z$z/z$Z
  z$z.cor[is.na(z$z.cor)] <- 0
  z$z.cor <- z$z.cor/raster::cellStats(z$z.cor, "max")
  class(z) <- c("NINA", "niche")

  out$w = w
  out$z = z

  message("\t...Success!")
  return(out)
}



