#' @title  ESTIMATE BETAS
#' @param y.list List of niche models
#'
#' @param C.matrix Competition matrix
#'
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
estimate_betas <- function(y.list, C.matrix = NULL){

  if (is.null(C.matrix)){
    C.matrix = 1 - diag(length(y.list))
    rownames(C.matrix) = names(y.list) ; colnames(C.matrix) = names(y.list)
  }
  betas <- list()
  for (n in names(y.list)){
    sum.co <- sum(stack(sapply(names(y.list), function(i) C.matrix[n, i]*y.list[[i]]$z.uncor)))
    K <- sum(stack(sapply(names(y.list), function(i) y.list[[i]]$z.uncor)))
    betas[[n]] <- 1 - sum.co/K
  }
  betas = betas[sapply(betas, function(i) maxValue(i) > 0)]

  return(betas)
}
