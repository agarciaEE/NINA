#' @title Summary
#'
#' @description  Summarizes the output models of an object class NINA
#'
#' @param x An object class NINA
#' @param ... Additional arguments for S3 methods
#'
#' @return A \code{NINA} class object description
#'
#' @method summary NINA
#'
#' @examples
#' \dontrun{
#' EN = EN_model(env_data, occ_data1, cluster = "env")
#' summary(EN)
#' }
#' @export
summary.NINA <- function(x, ...){
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

  cat("Niche model type: "); message(type, "\n")

  cat("Predictors: "); message(paste(env.var, collapse = " "))
  cat(x$pca$nf, "axis-components used\n")
  summary(x$pca)

  cat("Spatially constrained: "); message(mode.region)
  cat("Geographical extents: "); message(n.reg)
  print(x$tab)

  cat("Ensemble of regional models: "); message(mode.global)

  cat("Failures or warnings: "); message(fail)
  if(fail){print(x$fail)}

  cat("Models evaluation: "); message(eval)
  if(eval){print(x$eval$tab)}
}
