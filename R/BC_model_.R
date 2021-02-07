#' @title  ESTIMATE BC MODEL internal function
#'
#' @param z environmental niche model
#'
#' @param y.list list of environmental niche models of species interactors
#' @param id species name to be estimated. id must match the rownames in A.matrix
#' @param D Binary integer indicating species independence on interactions. 0 for fully dependent, 1 for species that can be present without interactors
#' @param A.matrix m by n matrix indicating the association coefficient (-1 to 1). m are species to be modelled as rows and n interactors as columns
#' @param C.matrix n by n matrix indicating the competition coefficient between interactors (0 to 1).
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
BC_model_ <- function(z, y.list, id, D = 0, A.matrix = NULL, C.matrix = NULL ){

  out = list()
  if(length(y.list) > 0){
    out$w <- estimate_w(y.list, id = id, A.matrix = A.matrix, C.matrix = C.matrix)
    wc <- out$w$z.uncor
  }
  else{ wc = 0}
  z$z.uncor <- z$z.uncor * D + z$z.uncor * wc
  z$z.uncor[is.na(z$z.uncor)] <- 0
  z$z.uncor <- z$z.uncor / cellStats(z$z.uncor, "max")
  z$z <- z$z.uncor * cellStats(z$z, "max")
  z$w <- z$z.uncor
  z$w[z$w > 0] <- 1
  z$z.cor <- z$z/z$Z
  z$z.cor[is.na(z$z.cor)] <- 0
  z$z.cor <- z$z.cor/cellStats(z$z.cor, "max")
  out$z = z
  message("Success!")
  return(out)
}




