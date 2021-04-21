#' @title ASSEMBLING MODELS FUNCTION
#'
#' @param modelsList List of bootstrap niche models
#' @param type String indicating the type of niche model
#' @param method String indicating method to weight bootstrap performance. Default "ACC"
#' @param threshold Threshold to select models cut-off. Default is 0.5
#' @param modelspath File path to the folder where bootstrap models can be found
#'
#' @description Assemble bootstrap models
#'
#' @return Ensemble NINA model
#'
#' @details Returns an error if \code{filename} does not exist.
#'
#'
#' @importFrom tidyr spread
#' @importFrom plyr ldply
#'
#' @export
assemble_models <- function(modelsList,
                                method = c("ACC", "Jaccard Similarity", "TSS", "AUC", "kappa"),
                                threshold = 0.5,  modelspath = "./"){

  method = method[1]
  type = NULL
  bootstrap.eval <- list()
  eval = F
  cluster = F
  if(is.list(modelsList)){
    if (any(sapply(modelsList, function(x) class(x)[1]) != "NINA")){stop("Some or all elements of the list are not NINA models")}
    if (any(!sapply(modelsList, function(x) class(x)[2]) %in% c("ENmodel", "BCmodel", "ECmodel"))){stop("Some or all elements of the list are not NINA models")}
    type = sapply(modelsList, function(x) class(x)[2])
    if (length(unique(type)) != 1){
      stop("Models on the list are different of different types")
    } else { type = unique(type)}
    if(!is.null(modelsList[[1]]$eval)){
      bootstrap.eval = lapply(modelsList, function(x) x$eval$tab)
      eval = T
    }
    else{
      warning("Models not evaluated previously. Considering equal preformance...")
      score = 1
    }
    env.var = modelsList[[1]]$predictors
    env.scores = modelsList[[1]]$env.scores
    sp.scores = lapply(modelsList, function(x) x$sp.scores)
    ras = modelsList[[1]]$maps[[1]]
    crs = modelsList[[1]]$crs
    names(sp.scores) = as.character(1:length(modelsList))
    if(type == "ECmodel"){
      t <- lapply(modelsList, function(x) x$t.mod)
      g <- lapply(modelsList, function(x) x$g)
    }
    if(type == "BCmodel"){
      w <- lapply(modelsList, function(x) x$w)
    }
    if (!is.null(modelsList[[1]]$clus)){
      clus.df = modelsList[[1]]$clus
      cluster = T
    }
    z <- lapply(modelsList, function(x) x$z.mod)
  }
  else {
    z <- t <- w <- g <- list()
    sp.scores <- list()
    for (i in 1:length(modelsList)){
      model <- readRDS(paste0(modelspath, modelsList[i]))
      if (class(model)[1] != "NINA"){stop(paste("Element", i, "of the list is not a NINA model"))}
      if (!class(model)[2] %in% c("ENmodel", "BCmodel", "ECmodel")){stop(paste("Element", i, "of the list is not a NINA model"))}
      if(!is.null(type) && class(model)[2] != type){stop(paste("Element", i, "of the list is of different type compared to previous element"))}
      type = class(model)[2]
      if(!is.null(model$eval)){
        bootstrap.eval[[modelsList[i]]] = model$eval$tab
        eval = T
      }
      else{
        warning("Models not evaluated previously. Considering equal preformance...")
        score = 1
      }
      env.var = model$predictors
      env.scores = model$env.scores
      ras = model$maps[[1]]
      crs = model$crs
      sp.scores[[modelsList[i]]] <- model$sp.scores
      if (!is.null(model$clus)){
        clus.df = model$clus
        cluster = T
      }
      if(type == "ECmodel"){
        t[[modelsList[i]]] <- model$t.mod
        g[[modelsList[i]]] <- model$g
      }
      if(type == "BCmodel"){
        w[[modelsList[i]]] <- model$w
      }
      z[[modelsList[i]]] <- model$z.mod
      rm(model)
    }
  }
  model <- list()
  if(length(bootstrap.eval) == 0){ bootstrap.eval = NULL}
  sp.scores <- plyr::ldply(sp.scores, .id = "bootstrap")
  if(cluster){
    if(type == "BCmodel"){
      w <- assemble_snm_bootstraps(w, env.scores, sp.scores = sp.scores,
                           bootstrap.eval = bootstrap.eval, eval = eval, threshold = threshold,
                           cluster = cluster, method = method)
      model$w = w$z.mod
    }
    if(type == "ECmodel"){
      t <- assemble_snm_bootstraps(t, env.scores, sp.scores = sp.scores,
                           bootstrap.eval = bootstrap.eval, eval = eval, threshold = threshold,
                           cluster = cluster, method = method)
      model$t.mod = t$z.mod
      g <- assemble_snm_bootstraps(g, env.scores, sp.scores = sp.scores,
                           bootstrap.eval = bootstrap.eval, eval = eval, threshold = threshold,
                           cluster = cluster, method = method)
      model$g = g$z.mod
    }
    z <- assemble_snm_bootstraps(z, env.scores,sp.scores = sp.scores,
                         bootstrap.eval = bootstrap.eval, eval = eval, threshold = threshold,
                         cluster = cluster, method = method)
    tab = cbind(ldply(sapply(z, function(x) names(x)), data.frame, .id = "region"), P = 1)
    tab =  spread(tab, "region", "P")
    rownames(tab) <- tab[,1]
    tab <- tab[,-1]
    tab[is.na(tab)] = 0
    z = reverse_list(z)
    mod.Val = sapply(names(z), function(i) niche_to_dis(env.scores, z[[i]], cluster = clus.df, cor = FALSE)[,3])
    mod.Val[is.na(mod.Val)] = 0
    model$pred.dis = cbind(env.scores[,1:2], mod.Val)
    model$tab = tab
    model$clus = clus.df
    if(!is.null(bootstrap.eval)){
      model$bootstrap.eval = bootstrap.eval
    }
  }
  else{
    if(type == "BCmodel"){
      w <- assemble_snm_bootstraps(w, env.scores, sp.scores = sp.scores,
                           bootstrap.eval = bootstrap.eval, eval = eval, threshold = threshold,
                           cluster = cluster, method = method)
      model$w = w$z.mod
    }
    if(type == "ECmodel"){
      t <- assemble_snm_bootstraps(t, env.scores, sp.scores = sp.scores,
                           bootstrap.eval = bootstrap.eval, eval = eval, threshold = threshold,
                           cluster = cluster, method = method)
      model$t.mod = t$z.mod
      g <- assemble_snm_bootstraps(g, env.scores, sp.scores = sp.scores,
                           bootstrap.eval = bootstrap.eval, eval = eval, threshold = threshold,
                           cluster = cluster, method = method)
      model$g = g$z.mod
    }
    z <- assemble_snm_bootstraps(z, env.scores,sp.scores = sp.scores,
                         bootstrap.eval = bootstrap.eval, eval = eval, threshold = threshold,
                         cluster = cluster, method = method)
    mod.Val = sapply(names(z), function(i) niche_to_dis(env.scores, z[[i]], cor = FALSE)[,3])
    mod.Val[is.na(mod.Val)] = 0
    model$tab = table(names(z))
    if(!is.null(bootstrap.eval)){
      model$bootstrap.eval = bootstrap.eval
    }
  }
  model$pred.dis = cbind(env.scores[,1:2], mod.Val)
  model$z.mod = z
  model$env.scores = env.scores
  model$sp.scores = sp.scores
  model$maps  = raster_projection(model$pred.dis, ras = ras, crs = crs)
  model$predictors = env.var
  model$crs = crs
  #model$type = type
  attr(model, "class") <- c("NINA", type)
  message("All models have been assembled.")
  return(model)
}
