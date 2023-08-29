#' Get models evaluation
#'
#' @param x Ecological Niche Model of class NINA.
#'
#' @return List object containing models evaluation results.
#' @export
#'
#' @examples
#' EN <- EN_model(env_data, occ_data1, cluster = "env", n.clus = 5)
#' eval <- getEvaluation(EN)
getEval <- function(x){

  if (sum(class(x) %in% c("NINA", "ENmodel", "BCmodel", "ECmodel")) != 2) {
    stop("Argument 'x' is not a NINA species model")
  }
  if (!is.null(x$eval)){
    return(x$eval)
  } else {
    warning("Models evluation has not been performed.", immediate. = T)
    input <- readline(prompt="Do you want to perform models evaluation with standard arguments? (y/n)")
    if (input %in% c("y", "Y", "yes", "Yes", "YES")){
      eval <- models_evaluation(x)
      return(eval)
    }
  }
}
