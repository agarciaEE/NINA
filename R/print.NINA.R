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
#' EN = EN_model(env_data, occ_data1, cluster = "env", n.clus = 5)
#' print(EN)
#' }
#' \dontrun{
#' EN = EN_model(env_data, occ_data2)
#' print(EN)
#' }
#' @export
print.NINA <- function(x, ...){

  type = class(x)[2]

  if (type %in% c("ENmodel", "BCmodel", "ECmodel")){
    if (type == "ENmodel"){ type = "Environmental-only"}
    else if (type == "BCmodel"){ type = "Environmental-constrained"}
    else if (type == "ECmodel"){ type = "Ecological"}

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
  else if (type == "eval"){
    cat("Object class: NINA\n\n")
    cat("Models evaluation: \n=========================")
    cat("\nPresences/Absences:\n"); print(x$n)
    cat("\nEvaluation:\n"); print(x$tab)
    cat("\nThresholds:\n"); print(x$threshold)
    cat("\nCases:\n"); print(x$cases)
  }
  else if (type == "metrics"){
    cat("Object class: NINA\n\n")
    cat("\nNiche metrics: \n=========================")
    cat("\nNiche overlap:\n"); print(x$D$obs); cat("\np-value:\n"); print(x$D$pvalue)
    cat("\nNiche position:\n"); print(x$np$obs)
    cat("\nCentrid Shift:\n"); print(x$CS$obs); cat("\np-value:\n"); print(x$CS$pvalue)
    cat("\nAxes coefficiients:\n"); print(x$CS$weights)
    cat("\nEnvironmental Relocation:\n"); print(x$ER$obs); cat("\np-value:\n"); print(x$ER$pvalue)
    cat("\nUSE metrics:\n"); print(x$USE$obs); cat("\np-value:\n"); print(x$USE$pvalue)
  }
}
