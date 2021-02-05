#' @title Print
#'
#' @description  Prints characteristics of an object class NINA
#'
#' @param x An object class NINA
#' @param ... Additional arguments for S3 methods
#'
#' @return A description of a \code{NINA} class object
#'
#' @method print NINA
#'
#' @examples
#' \dontrun{
#' EN = EN_model(env_data, occ_data1, cluster = "env")
#' print(EN)
#' }
#' @export
print.NINA <- function(x, ...){
  type = x$type
  if (type == "EN"){ type = "Environmental-only"}
  if (type == "BC"){ type = "Environmental-constrained"}
  if (type == "EC"){ type = "Ecological"}

  env.var = x$predictors

  mode.region = if(!is.null(x$clus)){TRUE} else {FALSE}
  mode.global = if(!is.null(x$z.mod.global)){TRUE} else {FALSE}
  fail = if(!is.null(x$fail)){TRUE} else {FALSE}
  eval = if(!is.null(x$eval)){TRUE} else {FALSE}

  if (mode.region){
    n.reg = nrow(x$tab)
    n.sps = ncol(x$tab)
  } else {
    n.reg = 1
    n.sps = length(x$tab)
  }

  cat("Object class: NINA\n\n")
  cat("Niche model type: "); message(type)
  cat("Predictors: "); message(paste(env.var, collapse = " "))
  cat("Spatially constrained: "); message(mode.region)
  cat("Geographical extents: "); message(n.reg)
  cat("Ensemble of regional models: "); message(mode.global)
  cat("Number of species: "); message(n.sps)
  cat("Failures or warnings: "); message(fail)
  cat("Models evaluation: "); message(eval)

}
