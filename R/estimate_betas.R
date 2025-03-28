#' @title  ESTIMATE BETAS
#' @param y.list List of niche models
#' @param C.matrix Competition matrix
#' @param cor Logical
#' @param K = Carrying capacity of each environmental cell
#' @description Estimate the beta values of interactors
#'
#' @return list of lists
#'
#'@details Returns an error if \code{filename} does not exist.
#'
#' @examples
#' \dontrun{
#' accident_2015 <- fars_read("Project/data/accident_2015.csv.bz2")
#' }
#'
#' @importFrom raster stack maxValue
#'
#' @export
estimate_betas <- function(y.list, C.matrix = NULL, cor = F, K = NULL){

  if (is.null(C.matrix)){
    C.matrix = 1 - diag(length(y.list))
    rownames(C.matrix) = names(y.list) ; colnames(C.matrix) = names(y.list)
  }
  betas <- list()
  for (n in names(y.list)){
    if (cor){
      if (is.null(K)){
        #K <- length(y.list) # K equals to the nnumber of species predicted to be present in each cell
        K <- sum(stack(sapply(names(y.list), function(i) y.list[[i]]$z.cor))) # environment is compartmentalize by the species and is considered fully occupied
      }
      z <- y.list[[n]]$z.cor
      sum.co <- sum(stack(sapply(names(y.list), function(i) C.matrix[n, i]*y.list[[i]]$z.cor)))
    } else {
      if (is.null(K)){
        #K <- length(y.list) # K equals to the nnumber of species predicted to be present in each cell
        K <- sum(stack(sapply(names(y.list), function(i) y.list[[i]]$z.uncor))) # environment is compartmentalize by the species and is considered fully occupied
      }
      z <- y.list[[n]]$z.uncor
      sum.co <- sum(stack(sapply(names(y.list), function(i) C.matrix[n, i]*y.list[[i]]$z.uncor)))
    }
    betas[[n]] <- z * (1 - sum.co/K)
  }
  betas = betas[sapply(betas, function(i) maxValue(i) > 0)]

  return(betas)
}
