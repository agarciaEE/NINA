#' @title Summary
#'
#' @description  Summarizes the features of an object class NINA and prints a description
#'
#' @param object An object class NINA
#' @param ... Additional arguments for S3 methods
#'
#' @return A list summarizing the features of a \code{NINA} class object
#'
#' @method summary NINA
#'
#' @examples
#' \dontrun{
#' EN = EN_model(env_data, occ_data1, cluster = "env", n.clus = 5)
#' EN.summary = summary(EN)
#' }
#' \dontrun{
#' EN = EN_model(env_data, occ_data2)
#' summary(EN)
#' }
#' @export
summary.NINA <- function(object, ...){

  type = class(object)[2]

  if (type %in% c("ENmodel", "BCmodel", "ECmodel")){
    if (type == "ENmodel"){ type = "Environmental-only"}
    else if (type == "BCmodel"){ type = "Environmental-constrained"}
    else if (type == "ECmodel"){ type = "Ecological"}

    env.var = object$predictors

    mode.region = if(!is.null(object$clus)){TRUE} else {FALSE}
    mode.global = if(!is.null(object$z.mod.global)){TRUE} else {FALSE}
    fail = if(!is.null(object$fail)){TRUE} else {FALSE}
    eval = if(!is.null(object$eval)){TRUE} else {FALSE}

    if (mode.region){
      n.reg = nrow(object$tab)
      n.sps = ncol(object$tab)
    } else {
      n.reg = 1
      n.sps = length(object$tab)
    }

    cat("Niche model type: "); message(type, "\n")

    cat("Predictors: "); message(paste(env.var, collapse = " "))
    cat(object$pca$nf, "axis-components used\n")
    summary(object$pca)

    cat("Spatially constrained: "); message(mode.region)
    cat("Geographical extents: "); message(n.reg)
    print(object$tab)

    cat("Ensemble of regional models: "); message(mode.global)

    cat("Failures or warnings: "); message(fail)
    if(fail){print(object$fail)}

    cat("Models evaluation: "); message(eval)
    if(eval){print(object$eval$tab)}

    return(invisible(list(type = type,
                      predictors = env.var,
                      nf = object$pca$nf,
                      pca = object$pca,
                      mode.region = c(regions = mode.region, num.regions = n.reg),
                      mode.global = mode.global,
                      tab = object$tab,
                      failures = list(failures =fail, tab = object$fail),
                      eval = list(evaluation = eval, tab = object$eval))))
  }
  else if (type == "eval"){
    cat("Object class: NINA\n\n")
    cat("Models evaluation output: \n=========================")
    cat("Presences/Absences: "); print(object$n)
    cat("Evaluation: "); print(object$tab)
    cat("Thresholds: "); print(object$threshold)
    cat("Cases: "); print(object$cases)
  }
  else if (type == "metrics"){
    cat("Object class: NINA\n\n")
    cat("\nNiche metrics: \n=========================")
    cat("\nNiche overlap:\n"); print(object$D$obs); cat("\np-value:\n"); print(object$D$pvalue)
    cat("\nNiche position:\n"); print(object$np$obs)
    cat("\nCentrid Shift:\n"); print(object$CS$obs); cat("\np-value:\n"); print(object$CS$pvalue)
    cat("\nAxes coefficiients:\n"); print(object$CS$weights)
    cat("\nEnvironmental Relocation:\n"); print(object$ER$obs); cat("\np-value:\n"); print(object$ER$pvalue)
    cat("\nUSE metrics:\n"); print(object$USE$obs); cat("\np-value:\n"); print(object$USE$pvalue)
  }
}
