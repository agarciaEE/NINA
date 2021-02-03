#' @title ASSEMBLING MODELS FUNCTION
#'
#' @param modelsList List of bootstrap niche models
#'
#' @param type String indicating the type of niche model
#' @param method String indicating method to weight bootstrap performance. Default "ACC"
#' @param threshold Threshold to select models cut-off. Default is 0.5
#' @param modelspath File path to the folder where bootstrap models can be found
#'
#' @description Assemble bootstrap models
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
#'
#' @export
assemble_snm_models <- function(modelsList, type = c("EN", "BC", "EC"), method = c("ACC", "Jaccard Similarity", "TSS", "AUC", "kappa"), threshold = 0.5,  modelspath = "./"){

  method = method[1]
  type = type[1]
  if(is.list(modelsList)){
    if(!is.null(modelsList[[1]]$eval)){
      bootstrap.eval = lapply(modelsList, function(x) x$eval$tab)
    }
    else{
      warning("Models not evaluated previously. Considering equal preformance...")
      eval = F
      score = 1
    }
    env.scores = modelsList[[1]]$env.scores
    sp.scores <- ldply(lapply(modelsList, function(x) x$sp.scores), data.frame, .id = "bootstrap")
    if(type == "EC"){
      t <- lapply(modelsList, function(x) x$t.mod)
      w <- lapply(modelsList, function(x) x$w)
    }
    if(type == "BC"){
      w <- lapply(modelsList, function(x) x$w)
    }
    if (!is.null(modelsList[[1]]$clus)){
      clus.df = modelsList[[1]]$clus
      cluster = T
    }
    z <- lapply(modelsList, function(x) x$z.mod)
  }
  else {
    bootstrap.eval <- list()
    z <- list()
    t <- list()
    w <- list()
    sp.scores <- list()
    for (i in 1:length(modelsList)){
      model <- readRDS(paste0(modelspath, modelsList[i]))
      if(!is.null(model$eval)){
        bootstrap.eval[[modelsList[i]]] = model$eval$tab
        eval = T
      }
      else{
        warning("Models not evaluated previously. Considering equal preformance...")
        eval = F
        score = 1
      }
      env.scores = model$env.scores
      sp.scores[[modelsList[i]]] <- model$sp.scores
      if (!is.null(model$clus)){
        clus.df = model$clus
        cluster = T
      }
      if(type == "EC"){
        t[[modelsList[i]]] <- model$t.mod
        w[[modelsList[i]]] <- model$w
      }
      if(type == "BC"){
        w[[modelsList[i]]] <- model$w
      }
      z[[modelsList[i]]] <- model$z.mod
      rm(model)
    }
  }
  model <- list()
  sp.scores <- ldply(sp.scores, data.frame, .id = "bootstrap")
  if(cluster){
    if(type == "BC"){
      w <- assemble_snm_bootstraps(w, env.scores, sp.scores = sp.scores, w = T,
                           bootstrap.eval = bootstrap.eval, eval = eval, threshold = threshold,
                           cluster = cluster, method = method)
      model$w = w$z.mod
    }
    if(type == "EC"){
      t <- assemble_snm_bootstraps(t, env.scores, sp.scores = sp.scores,
                           bootstrap.eval = bootstrap.eval, eval = eval, threshold = threshold,
                           cluster = cluster, method = method)
      model$t.mod = t$z.mod
      w <- assemble_snm_bootstraps(w, env.scores, sp.scores = sp.scores, w = T,
                           bootstrap.eval = bootstrap.eval, eval = eval, threshold = threshold,
                           cluster = cluster, method = method)
      model$w = w$z.mod
    }
    z <- assemble_snm_bootstraps(z, env.scores,sp.scores = sp.scores,
                         bootstrap.eval = bootstrap.eval, eval = eval, threshold = threshold,
                         cluster = cluster, method = method)
    z.mod = z$z.mod
    model$z.mod = z.mod
    tab = cbind(ldply(sapply(z.mod, function(x) names(x)), data.frame, .id = "region"), P = 1)
    tab =  spread(tab, "region", "P")
    rownames(tab) <- tab[,1]
    tab <- tab[,-1]
    tab[is.na(tab)] = 0
    z.mod = reverse_list(z.mod)
    mod.Val = sapply(names(z.mod), function(i) niche_to_dis(env.scores, z.mod[[i]], cluster = clus.df, cor = FALSE)[,3])
    mod.Val[is.na(mod.Val)] = 0
    model$pred.dis = cbind(env.scores[,1:2], mod.Val)
    model$tab = tab
    model$clus = clus.df
    model$env.scores = env.scores
    model$sp.scores = sp.scores
    if(!is.null(bootstrap.eval)){
      model$bootstrap.eval = bootstrap.eval
    }
    message("Model assembling completed.")
  }
  else{
    if(type == "BC"){
      w <- assemble_snm_bootstraps(w, env.scores, sp.scores = sp.scores, w = T,
                           bootstrap.eval = bootstrap.eval, eval = eval, threshold = threshold,
                           cluster = cluster, method = method)
      model$w = w$z.mod
    }
    if(type == "EC"){
      t <- assemble_snm_bootstraps(t, env.scores, sp.scores = sp.scores,
                           bootstrap.eval = bootstrap.eval, eval = eval, threshold = threshold,
                           cluster = cluster, method = method)
      model$t.mod = t$z.mod
      w <- assemble_snm_bootstraps(w, env.scores, sp.scores = sp.scores, w = T,
                           bootstrap.eval = bootstrap.eval, eval = eval, threshold = threshold,
                           cluster = cluster, method = method)
      model$w = w$z.mod
    }
    z <- assemble_snm_bootstraps(z, env.scores,sp.scores = sp.scores, w = NULL,
                         bootstrap.eval = bootstrap.eval, eval = eval, threshold = threshold,
                         cluster = cluster, method = method)
    z.mod = z$z.mod
    mod.Val = sapply(names(z.mod), function(i) niche_to_dis(env.scores, z.mod[[i]], cor = FALSE)[,3])
    mod.Val[is.na(mod.Val)] = 0
    model$pred.dis = cbind(env.scores[,1:2], mod.Val)
    model$z.mod = z.mod
    model$env.scores = env.scores
    model$sp.scores = sp.scores
    if(!is.null(bootstrap.eval)){
      model$bootstrap.eval = bootstrap.eval
    }
    message("Model assembling completed.")
  }
  return(model)
}
