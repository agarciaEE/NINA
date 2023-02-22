#' @title NICHE PARAMETERS
#'
#' @param x NINA model
#' @param y Optional. NINA model
#' @param type String indicating the scale of the analysis. Default is "region"
#' @param centroid.w Logial to indicate if niche overlap estimate is to be weighted in relation to niche centroid
#' @param np.type If \code{centroid.w} is TRUE, type of centroid to estimate. Default is 'unimodal'. See \code{\link[NINA]{niche_position}}
#' @param np.metric If \code{centroid.w} is TRUE, statistic method to estimate the niche centroid. Default is 'median'. See \code{\link[NINA]{niche_position}}
#' @param quantile If \code{centroid.w} is TRUE, quantle threshold to filter niche densities. Default is 0.5. See \code{\link[NINA]{niche_position}}
#' @param centroid.w logical whether to weight niche overlap by distance to niche centroid
#' @param cor Logical whether to use environmentally corrected densites
#' @param rep Number of random samples to carry out the statistic tests
#' @param rand Directions of the comparisons. 1 performs all tests with argument 'x' as base of all comparisons. 2 uses argument 'y'. Default is 1.
#' @param alternative String indicating the alternative hypothesis. Default is "greater".
#' @param rnd.test Logical to indicate if randomization tests are to be estimated
#'
#' @description Estimates niche parameters
#'
#' @return Data frame.
#'
#' @details Returns an error if \code{filename} does not exist.
#'
#'
#' @importFrom stats quantile
#' @importFrom ecospat ecospat.niche.overlap
#' @importFrom raster extent compareRaster maxValue values
#' @importFrom plyr ldply
#'
#' @export
niche_parameters <- function(x, y, type = c("region", "global"), centroid.w = F, rnd.test = F,
                             quantile = 0.75, np.type = "unimodal", np.metric = c("median", "mean", "max"),
                             cor = F, rep = 100, rand = 1, alternative = c("greater", "lower")){

  type = type[1]
  np.metric = np.metric[1]

  xW <- yW <- xG <- yG <- clus <-  F

  xmod = NULL
  ymod = NULL

  if (all(class(x) == c("NINA", "BCmodel"))) {
    X = list(EnvSuit = x$z.mod,
             EnvAcc = x$w)
    xW = T
    xmod = "BC"
  }
  else if (all(class(x) == c("NINA", "ECmodel"))) {
    X = list(EnvSuit = x$z.mod,
             EcoSuit = x$t.mod,
             EcoAcc = x$g)
    xG = T
    xmod = "EC"
  }
  else if (all(class(x) == c("NINA", "ENmodel"))) {
    X = list(EnvSuit = x$z.mod)
    xmod = "EN"
    if (missing(y)){ stop("Supplied model has nothing to compare to...\n\t...Provide another NINA model as 'y' argument")}
  }
  else{ stop("Object 'x' is not a NINA model")}
  if(!is.null(x$clus)){clus = T}

  if(!missing(y)){
    if (all(class(y) == c("NINA", "BCmodel"))) {
      Y = list(EnvSuit = y$z.mod,
               EnvAcc = y$w)
      yW = T
      ymod = "BC"
    }
    else if (all(class(y) == c("NINA", "ECmodel"))) {
      Y = list(EnvSuit = y$z.mod,
               EcoSuit = y$t.mod,
               EcoAcc = y$g)
      yG = T
      ymod = "EC"
    }
    else if (all(class(y) == c("NINA", "ENmodel"))) {
      Y = list(EnvSuit = y$z.mod)
      ymod = "EN"
    }
    else{ stop("Object 'y' is not a NINA model")}
    if(!is.null(y$clus) != clus) {stop("Model 'x' and model 'y' have been estimated in different scales") }

    if(!all(names(x$maps) %in% names(y$maps))) {
      if(all(!names(x$maps) %in% names(y$maps))) {
        stop("Species in model 'x' not found in model 'y'")
      }
      else {
        sp.absent = names(x$maps)[!names(x$maps) %in% names(y$maps)]
        warning(paste(sp.absent, "in model 'x' are not present in model 'y'"), immediate. = T)
      }
    }
  }
  out = list()
  np = list()
  p.values <- list()
  simList <- list()
  if (type == "region"){
    if (clus){
      for (e in names(X[[1]])){
        message(paste("Estimating niche metrics of species in region", e, "..."))
        np[[e]] = list()
        p.values[[e]] <- list()
        simList[[e]] <- list()
        for (i in names(X[[1]][[e]])){
          message(paste("Computing niche metrics of", i, "..."), appendLF = F)
          if (xmod == "BC" && is.null(ymod)){
            if(!is.null(X$EnvSuit[[e]][[i]]) && !is.null(X$EnvAcc[[e]][[i]])){
              message("\n\t...between environmental niche and environmental accessibillity....", appendLF = F)
              ncomp <-  niche_comparison(X$EnvSuit[[e]][[i]], X$EnvAcc[[e]][[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                         quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                         cor = cor, rep = rep, rand = rand, alternative = alternative)
              np[[e]][[i]] <- c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                         ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$CS$weights, ncomp$USE$obs)),
                                comparison = "Env2Acc")
              if (rnd.test){
                simList[[e]][[i]] <- cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                           ER = ncomp$ER$sim, ncomp$USE$sim,
                                           comparison = "Env2Acc")

                p.values[[e]][[i]] <- cbind(alternative = alternative, D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue,
                                            CS = ncomp$CS$pvalue, ER = ncomp$ER$pvalue,
                                            ncomp$USE$pvalue, comparison = "Env2Acc")
                p.values[[e]][[i]] <- as.data.frame(p.values[[e]][[i]], row.names = rep("", nrow(p.values[[e]][[i]])))
              }
              message("...done.")
            }
          }
          else if (xmod == "EC" && is.null(ymod)){
            if(!is.null(X$EnvSuit[[e]][[i]]) &&  !is.null(X$EcoSuit[[e]][[i]])){
              message("\n\t...between environmental niche and ecological niche...", appendLF = F)
              ncomp <-  niche_comparison(X$EnvSuit[[e]][[i]], X$EcoSuit[[e]][[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                         quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                         cor = cor, rep = rep, rand = rand, alternative = alternative)
              np[[e]][[i]] <- c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                         ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                                comparison = "Env2Eco")

              if (rnd.test){
                simList[[e]][[i]] <- cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                         ER = ncomp$ER$sim, ncomp$USE$sim,
                                         comparison = "Env2Eco")

                p.values[[e]][[i]] <- cbind(alternative = alternative, D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue,
                                          CS = ncomp$CS$pvalue, ER = ncomp$ER$pvalue,
                                          ncomp$USE$pvalue, comparison = "Env2Eco")
              }
              message("...done.")
            }
            if(!is.null(X$EcoSuit[[e]][[i]]) &&  !is.null(X$EcoAcc[[e]][[i]])){
              message("\n\t...between ecological niche and ecological accessibillity....", appendLF = F)
              ncomp <-  niche_comparison(X$EcoSuit[[e]][[i]], X$EcoAcc[[e]][[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                         quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                         cor = cor, rep = rep, rand = rand, alternative = alternative)

              np[[e]][[i]] <- rbind(np[[e]][[i]], c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                                             ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                                                    comparison = "Eco2Acc"))
              if (rnd.test){
                simList[[e]][[i]] <- rbind(simList[[e]][[i]], cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                                                  ER = ncomp$ER$sim, ncomp$USE$sim,
                                                                  comparison = "Eco2Acc"))

                p.values[[e]][[i]] <- rbind(p.values[[e]][[i]], cbind(alternative = alternative,
                                                                    D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue, CS = ncomp$CS$pvalue,
                                                                    ER = ncomp$ER$pvalue, ncomp$USE$pvalue,
                                                                    comparison = "Eco2Acc"))
                p.values[[e]][[i]] <- as.data.frame(p.values[[e]][[i]], row.names = rep("", nrow(p.values[[e]][[i]])))
              }
              message("...done.")
            }
          }
          else if (!is.null(Y[[1]][[e]][[i]])){
            if (xmod == "EN" && ymod == "EN"){
              if(!is.null(X$EnvSuit[[e]][[i]]) &&  !is.null(Y$EnvSuit[[e]][[i]])){
                message("\n\t...between environmental niches....", appendLF = F)
                ncomp <-  niche_comparison(X$EnvSuit[[e]][[i]], Y$EnvSuit[[e]][[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                           quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                           cor = cor, rep = rep, rand = rand, alternative = alternative)

                np[[e]][[i]] <- c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                           ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                                  comparison = "Env2Env")

                if (rnd.test){
                  simList[[e]][[i]] <- cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                           ER = ncomp$ER$sim, ncomp$USE$sim,
                                           comparison = "Env2Env")

                  p.values[[e]][[i]] <- cbind(alternative = alternative, D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue,
                                            CS = ncomp$CS$pvalue, ER = ncomp$ER$pvalue,
                                            ncomp$USE$pvalue, comparison = "Env2Env")
                  p.values[[e]][[i]] <- as.data.frame(p.values[[e]][[i]], row.names = rep("", nrow(p.values[[e]][[i]])))
                }
                message("...done.")
              }

            }
            else if (xmod == "EN" && ymod == "BC"){
              if(!is.null(X$EnvSuit[[e]][[i]]) &&  !is.null(Y$EnvSuit[[e]][[i]])){
                message("\n\t...between environmental niches....", appendLF = F)
                ncomp <-  niche_comparison(X$EnvSuit[[e]][[i]], Y$EnvSuit[[e]][[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                           quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                           cor = cor, rep = rep, rand = rand, alternative = alternative)
                np[[e]][[i]] <- c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                           ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                                  comparison = "Env2Env")

                if (rnd.test){
                  simList[[e]][[i]] <- cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                           ER = ncomp$ER$sim, ncomp$USE$sim,
                                           comparison = "Env2Env")

                  p.values[[e]][[i]] <- cbind(alternative = alternative, D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue,
                                            CS = ncomp$CS$pvalue, ER = ncomp$ER$pvalue,
                                            ncomp$USE$pvalue, comparison = "Env2Env")
                }
                message("...done.")
              }
              if(!is.null(X$EnvSuit[[e]][[i]]) &&  !is.null(Y$EnvAcc[[e]][[i]])){
                message("\n\t...between environmental niche and environmental accessibillity....", appendLF = F)
                ncomp <-  niche_comparison(X$EnvSuit[[e]][[i]], Y$EnvAcc[[e]][[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                           quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                           cor = cor, rep = rep, rand = rand, alternative = alternative)

                np[[e]][[i]] <- rbind(np[[e]][[i]], c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                                               ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                                                      comparison = "Env2Acc"))
                if (rnd.test){
                  simList[[e]][[i]] <- rbind(simList[[e]][[i]], cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                                                    ER = ncomp$ER$sim, ncomp$USE$sim,
                                                                    comparison = "Env2Acc"))

                  p.values[[e]][[i]] <- rbind(p.values[[e]][[i]], cbind(alternative = alternative,
                                                                      D = ncomp$D$pvalue, CS = ncomp$CS$pvalue, breadth.diff = ncomp$breadth.diff$pvalue,
                                                                      ER = ncomp$ER$pvalue, ncomp$USE$pvalue,
                                                                      comparison = "Env2Acc"))
                  p.values[[e]][[i]] <- as.data.frame(p.values[[e]][[i]], row.names = rep("", nrow(p.values[[e]][[i]])))
                }
                message("...done.")
              }
            }
            else if (xmod == "EN" && ymod == "EC"){
              if(!is.null(X$EnvSuit[[e]][[i]]) &&  !is.null(Y$EnvSuit[[e]][[i]])){
                message("\n\t...between environmental niches....", appendLF = F)
                ncomp <-  niche_comparison(X$EnvSuit[[e]][[i]], Y$EnvSuit[[e]][[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                           quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                           cor = cor, rep = rep, rand = rand, alternative = alternative)
                np[[e]][[i]] <- c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                           ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                                  comparison = "Env2Env")

                if (rnd.test){
                  simList[[e]][[i]] <- cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                           ER = ncomp$ER$sim, ncomp$USE$sim,
                                           comparison = "Env2Env")

                  p.values[[e]][[i]] <- cbind(alternative = alternative, D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue,
                                            CS = ncomp$CS$pvalue, ER = ncomp$ER$pvalue,
                                            ncomp$USE$pvalue, comparison = "Env2Env")
                  p.values[[e]][[i]] <- as.data.frame(p.values[[e]][[i]], row.names = rep("", nrow(p.values[[e]][[i]])))
                }
                message("...done.")
              }
            }
            else if (xmod == "BC" && ymod == "BC"){
              if(!is.null(X$EnvSuit[[e]][[i]]) &&  !is.null(Y$EnvSuit[[e]][[i]])){
                message("\n\t...between environmental niches....", appendLF = F)
                ncomp <-  niche_comparison(X$EnvSuit[[e]][[i]], Y$EnvSuit[[e]][[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                           quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                           cor = cor, rep = rep, rand = rand, alternative = alternative)
                np[[e]][[i]] <- c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                           ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                                  comparison = "Env2Env")

                if (rnd.test){
                  simList[[e]][[i]] <- cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                           ER = ncomp$ER$sim, ncomp$USE$sim,
                                           comparison = "Env2Env")

                  p.values[[e]][[i]] <- cbind(alternative = alternative, D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue,
                                            CS = ncomp$CS$pvalue, ER = ncomp$ER$pvalue,
                                            ncomp$USE$pvalue, comparison = "Env2Env")
                }
                message("...done.")
              }
              if(!is.null(X$EnvSuit[[e]][[i]]) &&  !is.null(Y$EnvAcc[[e]][[i]])){
                message("\n\t...between environmental accessibilities...", appendLF = F)
                ncomp <-  niche_comparison(X$EnvAcc[[e]][[i]], Y$EnvAcc[[e]][[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                           quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                           cor = cor, rep = rep, rand = rand, alternative = alternative)

                np[[e]][[i]] <- rbind(np[[e]][[i]], c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                                               ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                                                      comparison = "Acc2Acc"))
                if (rnd.test){
                  simList[[e]][[i]] <- rbind(simList[[e]][[i]], cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                                                    ER = ncomp$ER$sim, ncomp$USE$sim,
                                                                    comparison = "Acc2Acc"))

                  p.values[[e]][[i]] <- rbind(p.values[[e]][[i]], cbind(alternative = alternative,
                                                                      D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue, CS = ncomp$CS$pvalue,
                                                                      ER = ncomp$ER$pvalue, ncomp$USE$pvalue,
                                                                      comparison = "Acc2Acc"))
                  p.values[[e]][[i]] <- as.data.frame(p.values[[e]][[i]], row.names = rep("", nrow(p.values[[e]][[i]])))
                }
                message("...done.")
              }
            }
            else if (xmod == "BC" && ymod == "EC"){
              if(!is.null(X$EnvSuit[[e]][[i]]) &&  !is.null(Y$EnvSuit[[e]][[i]])){
                message("\n\t...between environmental niches....", appendLF = F)
                ncomp <-  niche_comparison(X$EnvSuit[[e]][[i]], Y$EnvSuit[[e]][[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                           quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                           cor = cor, rep = rep, rand = rand, alternative = alternative)
                np[[e]][[i]] <- c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                           ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                                  comparison = "Env2Env")

                if (rnd.test){
                  simList[[e]][[i]] <- cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                           ER = ncomp$ER$sim, ncomp$USE$sim,
                                           comparison = "Env2Env")

                  p.values[[e]][[i]] <- cbind(alternative = alternative, D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue,
                                            CS = ncomp$CS$pvalue, ER = ncomp$ER$pvalue,
                                            ncomp$USE$pvalue, comparison = "Env2Env")
                }

                message("...done.")
              }
              if(!is.null(X$EnvSuit[[e]][[i]]) &&  !is.null(Y$EcoSuit[[e]][[i]])){
                message("\n\t...between environmental niche and ecological niche...", appendLF = F)
                ncomp <-  niche_comparison(X$EnvSuit[[e]][[i]], Y$EcoSuit[[e]][[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                           quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                           cor = cor, rep = rep, rand = rand, alternative = alternative)

                np[[e]][[i]] <- rbind(np[[e]][[i]], c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                                               ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                                                      comparison = "Env2Eco"))
                if (rnd.test){
                  simList[[e]][[i]] <- rbind(simList[[e]][[i]], cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                                                    ER = ncomp$ER$sim, ncomp$USE$sim,
                                                                    comparison = "Env2Eco"))

                  p.values[[e]][[i]] <- rbind(p.values[[e]][[i]], cbind(alternative = alternative,
                                                                      D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue, CS = ncomp$CS$pvalue,
                                                                      ER = ncomp$ER$pvalue, ncomp$USE$pvalue,
                                                                      comparison = "Env2Eco"))
                  p.values[[e]][[i]] <- as.data.frame(p.values[[e]][[i]], row.names = rep("", nrow(p.values[[e]][[i]])))
                }
                message("...done.")
              }
              if(!is.null(X$EnvAcc[[e]][[i]]) &&  !is.null(Y$EcoAcc[[e]][[i]])){
                message("\n\t...between environmental accessibility and ecological accessibility", appendLF = F)
                ncomp <-  niche_comparison(X$EnvAcc[[e]][[i]], Y$EcoAcc[[e]][[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                           quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                           cor = cor, rep = rep, rand = rand, alternative = alternative)

                np[[e]][[i]] <- rbind(np[[e]][[i]], c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                                               ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                                                      comparison = "Acc2Acc"))
                if (rnd.test){
                  simList[[e]][[i]] <- rbind(simList[[e]][[i]], cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                                                    ER = ncomp$ER$sim, ncomp$USE$sim,
                                                                    comparison = "Acc2Acc"))

                  p.values[[e]][[i]] <- rbind(p.values[[e]][[i]], cbind(alternative = alternative,
                                                                      D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue, CS = ncomp$CS$pvalue,
                                                                      ER = ncomp$ER$pvalue, ncomp$USE$pvalue,
                                                                      comparison = "Acc2Acc"))
                  p.values[[e]][[i]] <- as.data.frame(p.values[[e]][[i]], row.names = rep("", nrow(p.values[[e]][[i]])))
                }
                message("...done.")
              }
            }
            else if (xmod == "EC" && ymod == "EC"){
              if(!is.null(X$EnvSuit[[e]][[i]]) &&  !is.null(Y$EnvSuit[[e]][[i]])){
                message("\n\t...between environmental niches....", appendLF = F)
                ncomp <-  niche_comparison(X$EnvSuit[[e]][[i]], Y$EnvSuit[[e]][[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                           quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                           cor = cor, rep = rep, rand = rand, alternative = alternative)
                np[[e]][[i]] <- c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                           ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                                  comparison = "Env2Env")

                if (rnd.test){
                  simList[[e]][[i]] <- cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                           ER = ncomp$ER$sim, ncomp$USE$sim,
                                           comparison = "Env2Env")

                  p.values[[e]][[i]] <- cbind(alternative = alternative, D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue,
                                            CS = ncomp$CS$pvalue, ER = ncomp$ER$pvalue,
                                            ncomp$USE$pvalue, comparison = "Env2Env")
                }
                message("...done.")
              }
              if(!is.null(X$EcoSuit[[e]][[i]]) &&  !is.null(Y$EcoSuit[[e]][[i]])){
                message("\n\t...between ecological niches....", appendLF = F)
                ncomp <-  niche_comparison(X$EcoSuit[[e]][[i]], Y$EcoSuit[[e]][[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                           quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                           cor = cor, rep = rep, rand = rand, alternative = alternative)

                np[[e]][[i]] <- rbind(np[[e]][[i]], c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                                               ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                                                      comparison = "Eco2Eco"))
                if (rnd.test){
                  simList[[e]][[i]] <- rbind(simList[[e]][[i]], cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                                                    ER = ncomp$ER$sim, ncomp$USE$sim,
                                                                    comparison = "Eco2Eco"))

                  p.values[[e]][[i]] <- rbind(p.values[[e]][[i]], cbind(alternative = alternative,
                                                                      D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue, CS = ncomp$CS$pvalue,
                                                                      ER = ncomp$ER$pvalue, ncomp$USE$pvalue,
                                                                      comparison = "Eco2Eco"))
                  p.values[[e]][[i]] <- as.data.frame(p.values[[e]][[i]], row.names = rep("", nrow(p.values[[e]][[i]])))
                }
                message("...done.")
              }
              if(!is.null(X$EcoAcc[[e]][[i]]) &&  !is.null(Y$EnvAcc[[e]][[i]])){
                message("\n\t...between environmental acessibilities...", appendLF = F)
                ncomp <-  niche_comparison(X$EcoAcc[[e]][[i]], Y$EnvAcc[[e]][[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                           quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                           cor = cor, rep = rep, rand = rand, alternative = alternative)

                np[[e]][[i]] <- rbind(np[[e]][[i]], c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                                               ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                                                      comparison = "Acc2Acc"))
                if (rnd.test){
                  simList[[e]][[i]] <- rbind(simList[[e]][[i]], cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                                                    ER = ncomp$ER$sim, ncomp$USE$sim,
                                                                    comparison = "Acc2Acc"))

                  p.values[[e]][[i]] <- rbind(p.values[[e]][[i]], cbind(alternative = alternative,
                                                                      D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue, CS = ncomp$CS$pvalue,
                                                                      ER = ncomp$ER$pvalue, ncomp$USE$pvalue,
                                                                      comparison = "Acc2Acc"))
                  p.values[[e]][[i]] <- as.data.frame(p.values[[e]][[i]], row.names = rep("", nrow(p.values[[e]][[i]])))
                }
                message("...done.")
              }
            }
            else if (xmod == "BC" && ymod == "EN"){
              if(!is.null(X$EnvSuit[[e]][[i]]) &&  !is.null(Y$EnvSuit[[e]][[i]])){
                message("\n\t...between environmental niches....", appendLF = F)
                ncomp <-  niche_comparison(X$EnvSuit[[e]][[i]], Y$EnvSuit[[e]][[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                           quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                           cor = cor, rep = rep, rand = rand, alternative = alternative)
                np[[e]][[i]] <- c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                           ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                                  comparison = "Env2Env")

                if (rnd.test){
                  simList[[e]][[i]] <- cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                           ER = ncomp$ER$sim, ncomp$USE$sim,
                                           comparison = "Env2Env")

                  p.values[[e]][[i]] <- cbind(alternative = alternative, D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue,
                                            CS = ncomp$CS$pvalue, ER = ncomp$ER$pvalue,
                                            ncomp$USE$pvalue, comparison = "Env2Env")
                }
                message("...done.")
              }
              if(!is.null(X$EnvAcc[[e]][[i]]) &&  !is.null(Y$EnvSuit[[e]][[i]])){
                message("\n\t...between environmental acccessiibiliity and environmental niche....", appendLF = F)
                ncomp <-  niche_comparison(X$EnvAcc[[e]][[i]], Y$EnvSuit[[e]][[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                           quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                           cor = cor, rep = rep, rand = rand, alternative = alternative)

                np[[e]][[i]] <- rbind(np[[e]][[i]], c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                                               ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                                                      comparison = "Acc2Env"))
                if (rnd.test){
                  simList[[e]][[i]] <- rbind(simList[[e]][[i]], cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                                                    ER = ncomp$ER$sim, ncomp$USE$sim,
                                                                    comparison = "Acc2Env"))

                  p.values[[e]][[i]] <- rbind(p.values[[e]][[i]], cbind(alternative = alternative,
                                                                      D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue, CS = ncomp$CS$pvalue,
                                                                      ER = ncomp$ER$pvalue, ncomp$USE$pvalue,
                                                                      comparison = "Acc2Env"))
                  p.values[[e]][[i]] <- as.data.frame(p.values[[e]][[i]], row.names = rep("", nrow(p.values[[e]][[i]])))
                }
                message("...done.")
              }
            }
            else if (xmod == "EC" && ymod == "EN"){
              if(!is.null(X$EnvSuit[[e]][[i]]) &&  !is.null(Y$EnvSuit[[e]][[i]])){
                message("\n\t...between environmental niches....", appendLF = F)
                ncomp <-  niche_comparison(X$EnvSuit[[e]][[i]], Y$EnvSuit[[e]][[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                           quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                           cor = cor, rep = rep, rand = rand, alternative = alternative)
                np[[e]][[i]] <- c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                           ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                                  comparison = "Env2Env")

                if (rnd.test){
                  simList[[e]][[i]] <- cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                           ER = ncomp$ER$sim, ncomp$USE$sim,
                                           comparison = "Env2Env")

                  p.values[[e]][[i]] <- cbind(alternative = alternative, D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue,
                                            CS = ncomp$CS$pvalue, ER = ncomp$ER$pvalue,
                                            ncomp$USE$pvalue, comparison = "Env2Env")
                  p.values[[e]][[i]] <- as.data.frame(p.values[[e]][[i]], row.names = rep("", nrow(p.values[[e]][[i]])))
                }
                message("...done.")
              }
            }
            else if (xmod == "EC" && ymod == "BC"){
              if(!is.null(X$EnvSuit[[e]][[i]]) &&  !is.null(Y$EnvSuit[[e]][[i]])){
                message("\n\t...between environmental niches....", appendLF = F)
                ncomp <-  niche_comparison(X$EnvSuit[[e]][[i]], Y$EnvSuit[[e]][[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                           quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                           cor = cor, rep = rep, rand = rand, alternative = alternative)
                np[[e]][[i]] <- c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                           ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                                  comparison = "Env2Env")

                if (rnd.test){
                  simList[[e]][[i]] <- cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                           ER = ncomp$ER$sim, ncomp$USE$sim,
                                           comparison = "Env2Env")

                  p.values[[e]][[i]] <- cbind(alternative = alternative, D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue,
                                            CS = ncomp$CS$pvalue, ER = ncomp$ER$pvalue,
                                            ncomp$USE$pvalue, comparison = "Env2Env")
                }
                message("...done.")
              }
              if(!is.null(X$EcoSuit[[e]][[i]]) &&  !is.null(Y$EnvSuit[[e]][[i]])){
                message("\n\t...between ecological niche and environmental niche....", appendLF = F)
                ncomp <-  niche_comparison(X$EcoSuit[[e]][[i]], Y$EnvSuit[[e]][[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                           quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                           cor = cor, rep = rep, rand = rand, alternative = alternative)

                np[[e]][[i]] <- rbind(np[[e]][[i]], c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                                               ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                                                      comparison = "Eco2Env"))
                if (rnd.test){
                  simList[[e]][[i]] <- rbind(simList[[e]][[i]], cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                                                    ER = ncomp$ER$sim, ncomp$USE$sim,
                                                                    comparison = "Eco2Env"))

                  p.values[[e]][[i]] <- rbind(p.values[[e]][[i]], cbind(alternative = alternative,
                                                                      D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue, CS = ncomp$CS$pvalue,
                                                                      ER = ncomp$ER$pvalue, ncomp$USE$pvalue,
                                                                      comparison = "Eco2Env"))
                  p.values[[e]][[i]] <- as.data.frame(p.values[[e]][[i]], row.names = rep("", nrow(p.values[[e]][[i]])))
                }
                message("...done.")
              }
              if(!is.null(X$EcoAcc[[e]][[i]]) &&  !is.null(Y$EnvAcc[[e]][[i]])){
                message("\n\t...between ecological accessibillity and environmental accessibillity....", appendLF = F)
                ncomp <-  niche_comparison(X$EcoAcc[[e]][[i]], Y$EnvAcc[[e]][[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                           quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                           cor = cor, rep = rep, rand = rand, alternative = alternative)

                np[[e]][[i]] <- rbind(np[[e]][[i]], c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                                               ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                                                      comparison = "Acc2Acc"))
                if (rnd.test){
                  simList[[e]][[i]] <- rbind(simList[[e]][[i]], cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                                                    ER = ncomp$ER$sim, ncomp$USE$sim,
                                                                    comparison = "Acc2Acc"))

                  p.values[[e]][[i]] <- rbind(p.values[[e]][[i]], cbind(alternative = alternative,
                                                                      D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue, CS = ncomp$CS$pvalue,
                                                                      ER = ncomp$ER$pvalue, ncomp$USE$pvalue,
                                                                      comparison = "Acc2Acc"))
                  p.values[[e]][[i]] <- as.data.frame(p.values[[e]][[i]], row.names = rep("", nrow(p.values[[e]][[i]])))
                }
                message("...done.")
              }
            }
          }
          else{
            message(paste("...nothing to analyze."))
          }
        }
        np[[e]] <- plyr::ldply(np[[e]], .id = "species")
      }
      if (rnd.test){
        simList <- reverse_list(simList)
        p.values <- reverse_list(p.values)
        out$np <- plyr::ldply(np, .id = "region")
        out$sim <- lapply(simList, function(i) ldply(i, .id = "region"))
        out$pvalue <- lapply(p.values, function(i) ldply(i, .id = "region"))
      } else {
        out <- plyr::ldply(np, .id = "region")
      }
    }
    else {
      stop("Regional models not found")
    }
  }
  if (type == "global"){
    if(clus){
      if(!is.null(x$z.mod.global)){
        if(xmod == "EN"){
          X = list(EnvSuit = x$z.mod.global)
        } else if (xmod == "BC") {
          X = list(EnvSuit = x$z.mod.global,
                   EnvAcc = x$w.global)
        } else if (xmod == "EC") {
          X = list(EnvSuit = x$z.mod.global,
                   EcoSuit = x$t.mod.global,
                   EcoAcc = x$g.global)
        }
      } else {
        stop("Global models of argument 'x' not found")
      }
      if (!missing(y)){
        if(!is.null(y$z.mod.global)){
          if(ymod == "EN"){
            Y = list(EnvSuit = y$z.mod.global)
          } else if (ymod == "BC") {
            Y = list(EnvSuit = y$z.mod.global,
                     EnvAcc = y$w.global)
          } else if (ymod == "EC") {
            Y = list(EnvSuit = y$z.mod.global,
                     EcoSuit = y$t.mod.global,
                     EcoAcc = y$g.global)
          }
        } else {
          stop("Global models of argument 'y' not found")
        }
      }
    }
    np = list()
    p.values <- list()
    simList <- list()
    for (i in names(X[[1]])){
      message(paste("Computing niche metrics of", i, "..."), appendLF = F)
      if (xmod == "BC" && is.null(ymod)){
        if(!is.null(X$EnvSuit[[i]]) &&  !is.null(X$EnvAcc[[i]])){
          message("\n\t...between environmental niche and environmental accessibillity....", appendLF = F)
          ncomp <-  niche_comparison(X$EnvSuit[[i]], X$EnvAcc[[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                     quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                     cor = cor, rep = rep, rand = rand, alternative = alternative)
          np[[i]] <- c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                     ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$CS$weights, ncomp$USE$obs)),
                            comparison = "Env2Acc")
          if (rnd.test){
            simList[[i]] <- cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                       ER = ncomp$ER$sim, ncomp$USE$sim,
                                       comparison = "Env2Acc")

            p.values[[i]] <- cbind(alternative = alternative, D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue,
                                        CS = ncomp$CS$pvalue, ER = ncomp$ER$pvalue,
                                        ncomp$USE$pvalue, comparison = "Env2Acc")
            p.values[[i]] <- as.data.frame(p.values[[i]], row.names = rep("", nrow(p.values[[i]])))
          }
          message("...done.")
        }
      }
      else if (xmod == "EC" && is.null(ymod)){
        if(!is.null(X$EnvSuit[[i]]) &&  !is.null(X$EcoSuit[[i]])){
          message("\n\t...between environmental niche and ecological niche...", appendLF = F)
          ncomp <-  niche_comparison(X$EnvSuit[[i]], X$EcoSuit[[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                     quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                     cor = cor, rep = rep, rand = rand, alternative = alternative)
          np[[i]] <- c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                     ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                            comparison = "Env2Eco")

          if (rnd.test){
            simList[[i]] <- cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                       ER = ncomp$ER$sim, ncomp$USE$sim,
                                       comparison = "Env2Eco")

            p.values[[i]] <- cbind(alternative = alternative, D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue,
                                        CS = ncomp$CS$pvalue, ER = ncomp$ER$pvalue,
                                        ncomp$USE$pvalue, comparison = "Env2Eco")
          }
          message("...done.")
        }
        if(!is.null(X$EcoSuit[[i]]) &&  !is.null(X$EcoAcc[[i]])){
          message("\n\t...between ecological niche and ecological accessibillity....", appendLF = F)
          ncomp <-  niche_comparison(X$EcoSuit[[i]], X$EcoAcc[[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                     quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                     cor = cor, rep = rep, rand = rand, alternative = alternative)

          np[[i]] <- rbind(np[[i]], c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                                         ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                                                comparison = "Eco2Acc"))
          if (rnd.test){
            simList[[i]] <- rbind(simList[[i]], cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                                                ER = ncomp$ER$sim, ncomp$USE$sim,
                                                                comparison = "Eco2Acc"))

            p.values[[i]] <- rbind(p.values[[i]], cbind(alternative = alternative,
                                                                  D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue, CS = ncomp$CS$pvalue,
                                                                  ER = ncomp$ER$pvalue, ncomp$USE$pvalue,
                                                                  comparison = "Eco2Acc"))
            p.values[[i]] <- as.data.frame(p.values[[i]], row.names = rep("", nrow(p.values[[i]])))
          }
          message("...done.")
        }
      }
      else if (!is.null(Y[[1]][[i]])){
        if (xmod == "EN" && ymod == "EN"){
          if(!is.null(X$EnvSuit[[i]]) &&  !is.null(Y$EnvSuit[[i]])){
            message("\n\t...between environmental niches....", appendLF = F)
            ncomp <-  niche_comparison(X$EnvSuit[[i]], Y$EnvSuit[[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                       quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                       cor = cor, rep = rep, rand = rand, alternative = alternative)

            np[[i]] <- c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                       ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                              comparison = "Env2Env")

            if (rnd.test){
              simList[[i]] <- cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                         ER = ncomp$ER$sim, ncomp$USE$sim,
                                         comparison = "Env2Env")

              p.values[[i]] <- cbind(alternative = alternative, D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue,
                                          CS = ncomp$CS$pvalue, ER = ncomp$ER$pvalue,
                                          ncomp$USE$pvalue, comparison = "Env2Env")
              p.values[[i]] <- as.data.frame(p.values[[i]], row.names = rep("", nrow(p.values[[i]])))
            }
            message("...done.")
          }

        }
        else if (xmod == "EN" && ymod == "BC"){
          if(!is.null(X$EnvSuit[[i]]) &&  !is.null(Y$EnvSuit[[i]])){
            message("\n\t...between environmental niches....", appendLF = F)
            ncomp <-  niche_comparison(X$EnvSuit[[i]], Y$EnvSuit[[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                       quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                       cor = cor, rep = rep, rand = rand, alternative = alternative)
            np[[i]] <- c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                       ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                              comparison = "Env2Env")

            if (rnd.test){
              simList[[i]] <- cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                         ER = ncomp$ER$sim, ncomp$USE$sim,
                                         comparison = "Env2Env")

              p.values[[i]] <- cbind(alternative = alternative, D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue,
                                          CS = ncomp$CS$pvalue, ER = ncomp$ER$pvalue,
                                          ncomp$USE$pvalue, comparison = "Env2Env")
            }
            message("...done.")
          }
          if(!is.null(X$EnvSuit[[i]]) &&  !is.null(Y$EnvAcc[[i]])){
            message("\n\t...between environmental niche and environmental accessibillity....", appendLF = F)
            ncomp <-  niche_comparison(X$EnvSuit[[i]], Y$EnvAcc[[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                       quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                       cor = cor, rep = rep, rand = rand, alternative = alternative)

            np[[i]] <- rbind(np[[i]], c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                                           ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                                                  comparison = "Env2Acc"))
            if (rnd.test){
              simList[[i]] <- rbind(simList[[i]], cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                                                  ER = ncomp$ER$sim, ncomp$USE$sim,
                                                                  comparison = "Env2Acc"))

              p.values[[i]] <- rbind(p.values[[i]], cbind(alternative = alternative,
                                                                    D = ncomp$D$pvalue, CS = ncomp$CS$pvalue, breadth.diff = ncomp$breadth.diff$pvalue,
                                                                    ER = ncomp$ER$pvalue, ncomp$USE$pvalue,
                                                                    comparison = "Env2Acc"))
              p.values[[i]] <- as.data.frame(p.values[[i]], row.names = rep("", nrow(p.values[[i]])))
            }
            message("...done.")
          }
        }
        else if (xmod == "EN" && ymod == "EC"){
          if(!is.null(X$EnvSuit[[i]]) &&  !is.null(Y$EnvSuit[[i]])){
            message("\n\t...between environmental niches....", appendLF = F)
            ncomp <-  niche_comparison(X$EnvSuit[[i]], Y$EnvSuit[[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                       quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                       cor = cor, rep = rep, rand = rand, alternative = alternative)
            np[[i]] <- c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                       ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                              comparison = "Env2Env")

            if (rnd.test){
              simList[[i]] <- cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                         ER = ncomp$ER$sim, ncomp$USE$sim,
                                         comparison = "Env2Env")

              p.values[[i]] <- cbind(alternative = alternative, D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue,
                                          CS = ncomp$CS$pvalue, ER = ncomp$ER$pvalue,
                                          ncomp$USE$pvalue, comparison = "Env2Env")
              p.values[[i]] <- as.data.frame(p.values[[i]], row.names = rep("", nrow(p.values[[i]])))
            }
            message("...done.")
          }
        }
        else if (xmod == "BC" && ymod == "BC"){
          if(!is.null(X$EnvSuit[[i]]) &&  !is.null(Y$EnvSuit[[i]])){
            message("\n\t...between environmental niches....", appendLF = F)
            ncomp <-  niche_comparison(X$EnvSuit[[i]], Y$EnvSuit[[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                       quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                       cor = cor, rep = rep, rand = rand, alternative = alternative)
            np[[i]] <- c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                       ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                              comparison = "Env2Env")

            if (rnd.test){
              simList[[i]] <- cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                         ER = ncomp$ER$sim, ncomp$USE$sim,
                                         comparison = "Env2Env")

              p.values[[i]] <- cbind(alternative = alternative, D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue,
                                          CS = ncomp$CS$pvalue, ER = ncomp$ER$pvalue,
                                          ncomp$USE$pvalue, comparison = "Env2Env")
            }
            message("...done.")
          }
          if(!is.null(X$EnvSuit[[i]]) &&  !is.null(Y$EnvAcc[[i]])){
            message("\n\t...between environmental accessibilities...", appendLF = F)
            ncomp <-  niche_comparison(X$EnvAcc[[i]], Y$EnvAcc[[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                       quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                       cor = cor, rep = rep, rand = rand, alternative = alternative)

            np[[i]] <- rbind(np[[i]], c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                                           ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                                                  comparison = "Acc2Acc"))
            if (rnd.test){
              simList[[i]] <- rbind(simList[[i]], cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                                                  ER = ncomp$ER$sim, ncomp$USE$sim,
                                                                  comparison = "Acc2Acc"))

              p.values[[i]] <- rbind(p.values[[i]], cbind(alternative = alternative,
                                                                    D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue, CS = ncomp$CS$pvalue,
                                                                    ER = ncomp$ER$pvalue, ncomp$USE$pvalue,
                                                                    comparison = "Acc2Acc"))
              p.values[[i]] <- as.data.frame(p.values[[i]], row.names = rep("", nrow(p.values[[i]])))
            }
            message("...done.")
          }
        }
        else if (xmod == "BC" && ymod == "EC"){
          if(!is.null(X$EnvSuit[[i]]) &&  !is.null(Y$EnvSuit[[i]])){
            message("\n\t...between environmental niches....", appendLF = F)
            ncomp <-  niche_comparison(X$EnvSuit[[i]], Y$EnvSuit[[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                       quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                       cor = cor, rep = rep, rand = rand, alternative = alternative)
            np[[i]] <- c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                       ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                              comparison = "Env2Env")

            if (rnd.test){
              simList[[i]] <- cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                         ER = ncomp$ER$sim, ncomp$USE$sim,
                                         comparison = "Env2Env")

              p.values[[i]] <- cbind(alternative = alternative, D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue,
                                          CS = ncomp$CS$pvalue, ER = ncomp$ER$pvalue,
                                          ncomp$USE$pvalue, comparison = "Env2Env")
            }

            message("...done.")
          }
          if(!is.null(X$EnvSuit[[i]]) &&  !is.null(Y$EcoSuit[[i]])){
            message("\n\t...between environmental niche and ecological niche...", appendLF = F)
            ncomp <-  niche_comparison(X$EnvSuit[[i]], Y$EcoSuit[[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                       quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                       cor = cor, rep = rep, rand = rand, alternative = alternative)

            np[[i]] <- rbind(np[[i]], c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                                           ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                                                  comparison = "Env2Eco"))
            if (rnd.test){
              simList[[i]] <- rbind(simList[[i]], cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                                                  ER = ncomp$ER$sim, ncomp$USE$sim,
                                                                  comparison = "Env2Eco"))

              p.values[[i]] <- rbind(p.values[[i]], cbind(alternative = alternative,
                                                                    D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue, CS = ncomp$CS$pvalue,
                                                                    ER = ncomp$ER$pvalue, ncomp$USE$pvalue,
                                                                    comparison = "Env2Eco"))
              p.values[[i]] <- as.data.frame(p.values[[i]], row.names = rep("", nrow(p.values[[i]])))
            }
            message("...done.")
          }
          if(!is.null(X$EnvAcc[[i]]) &&  !is.null(Y$EcoAcc[[i]])){
            message("\n\t...between environmental accessibility and ecological accessibility", appendLF = F)
            ncomp <-  niche_comparison(X$EnvAcc[[i]], Y$EcoAcc[[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                       quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                       cor = cor, rep = rep, rand = rand, alternative = alternative)

            np[[i]] <- rbind(np[[i]], c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                                           ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                                                  comparison = "Acc2Acc"))
            if (rnd.test){
              simList[[i]] <- rbind(simList[[i]], cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                                                  ER = ncomp$ER$sim, ncomp$USE$sim,
                                                                  comparison = "Acc2Acc"))

              p.values[[i]] <- rbind(p.values[[i]], cbind(alternative = alternative,
                                                                    D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue, CS = ncomp$CS$pvalue,
                                                                    ER = ncomp$ER$pvalue, ncomp$USE$pvalue,
                                                                    comparison = "Acc2Acc"))
              p.values[[i]] <- as.data.frame(p.values[[i]], row.names = rep("", nrow(p.values[[i]])))
            }
            message("...done.")
          }
        }
        else if (xmod == "EC" && ymod == "EC"){
          if(!is.null(X$EnvSuit[[i]]) &&  !is.null(Y$EnvSuit[[i]])){
            message("\n\t...between environmental niches....", appendLF = F)
            ncomp <-  niche_comparison(X$EnvSuit[[i]], Y$EnvSuit[[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                       quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                       cor = cor, rep = rep, rand = rand, alternative = alternative)
            np[[i]] <- c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                       ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                              comparison = "Env2Env")

            if (rnd.test){
              simList[[i]] <- cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                         ER = ncomp$ER$sim, ncomp$USE$sim,
                                         comparison = "Env2Env")

              p.values[[i]] <- cbind(alternative = alternative, D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue,
                                          CS = ncomp$CS$pvalue, ER = ncomp$ER$pvalue,
                                          ncomp$USE$pvalue, comparison = "Env2Env")
            }
            message("...done.")
          }
          if(!is.null(X$EcoSuit[[i]]) &&  !is.null(Y$EcoSuit[[i]])){
            message("\n\t...between ecological niches....", appendLF = F)
            ncomp <-  niche_comparison(X$EcoSuit[[i]], Y$EcoSuit[[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                       quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                       cor = cor, rep = rep, rand = rand, alternative = alternative)

            np[[i]] <- rbind(np[[i]], c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                                           ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                                                  comparison = "Eco2Eco"))
            if (rnd.test){
              simList[[i]] <- rbind(simList[[i]], cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                                                  ER = ncomp$ER$sim, ncomp$USE$sim,
                                                                  comparison = "Eco2Eco"))

              p.values[[i]] <- rbind(p.values[[i]], cbind(alternative = alternative,
                                                                    D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue, CS = ncomp$CS$pvalue,
                                                                    ER = ncomp$ER$pvalue, ncomp$USE$pvalue,
                                                                    comparison = "Eco2Eco"))
              p.values[[i]] <- as.data.frame(p.values[[i]], row.names = rep("", nrow(p.values[[i]])))
            }
            message("...done.")
          }
          if(!is.null(X$EcoAcc[[i]]) &&  !is.null(Y$EnvAcc[[i]])){
            message("\n\t...between environmental acessibilities...", appendLF = F)
            ncomp <-  niche_comparison(X$EcoAcc[[i]], Y$EnvAcc[[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                       quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                       cor = cor, rep = rep, rand = rand, alternative = alternative)

            np[[i]] <- rbind(np[[i]], c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                                           ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                                                  comparison = "Acc2Acc"))
            if (rnd.test){
              simList[[i]] <- rbind(simList[[i]], cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                                                  ER = ncomp$ER$sim, ncomp$USE$sim,
                                                                  comparison = "Acc2Acc"))

              p.values[[i]] <- rbind(p.values[[i]], cbind(alternative = alternative,
                                                                    D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue, CS = ncomp$CS$pvalue,
                                                                    ER = ncomp$ER$pvalue, ncomp$USE$pvalue,
                                                                    comparison = "Acc2Acc"))
              p.values[[i]] <- as.data.frame(p.values[[i]], row.names = rep("", nrow(p.values[[i]])))
            }
            message("...done.")
          }
        }
        else if (xmod == "BC" && ymod == "EN"){
          if(!is.null(X$EnvSuit[[i]]) &&  !is.null(Y$EnvSuit[[i]])){
            message("\n\t...between environmental niches....", appendLF = F)
            ncomp <-  niche_comparison(X$EnvSuit[[i]], Y$EnvSuit[[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                       quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                       cor = cor, rep = rep, rand = rand, alternative = alternative)
            np[[i]] <- c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                       ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                              comparison = "Env2Env")

            if (rnd.test){
              simList[[i]] <- cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                         ER = ncomp$ER$sim, ncomp$USE$sim,
                                         comparison = "Env2Env")

              p.values[[i]] <- cbind(alternative = alternative, D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue,
                                          CS = ncomp$CS$pvalue, ER = ncomp$ER$pvalue,
                                          ncomp$USE$pvalue, comparison = "Env2Env")
            }
            message("...done.")
          }
          if(!is.null(X$EnvAcc[[i]]) &&  !is.null(Y$EnvSuit[[i]])){
            message("\n\t...between environmental acccessiibiliity and environmental niche....", appendLF = F)
            ncomp <-  niche_comparison(X$EnvAcc[[i]], Y$EnvSuit[[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                       quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                       cor = cor, rep = rep, rand = rand, alternative = alternative)

            np[[i]] <- rbind(np[[i]], c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                                           ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                                                  comparison = "Acc2Env"))
            if (rnd.test){
              simList[[i]] <- rbind(simList[[i]], cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                                                  ER = ncomp$ER$sim, ncomp$USE$sim,
                                                                  comparison = "Acc2Env"))

              p.values[[i]] <- rbind(p.values[[i]], cbind(alternative = alternative,
                                                                    D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue, CS = ncomp$CS$pvalue,
                                                                    ER = ncomp$ER$pvalue, ncomp$USE$pvalue,
                                                                    comparison = "Acc2Env"))
              p.values[[i]] <- as.data.frame(p.values[[i]], row.names = rep("", nrow(p.values[[i]])))
            }
            message("...done.")
          }
        }
        else if (xmod == "EC" && ymod == "EN"){
          if(!is.null(X$EnvSuit[[i]]) &&  !is.null(Y$EnvSuit[[i]])){
            message("\n\t...between environmental niches....", appendLF = F)
            ncomp <-  niche_comparison(X$EnvSuit[[i]], Y$EnvSuit[[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                       quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                       cor = cor, rep = rep, rand = rand, alternative = alternative)
            np[[i]] <- c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                       ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                              comparison = "Env2Env")

            if (rnd.test){
              simList[[i]] <- cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                         ER = ncomp$ER$sim, ncomp$USE$sim,
                                         comparison = "Env2Env")

              p.values[[i]] <- cbind(alternative = alternative, D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue,
                                          CS = ncomp$CS$pvalue, ER = ncomp$ER$pvalue,
                                          ncomp$USE$pvalue, comparison = "Env2Env")
              p.values[[i]] <- as.data.frame(p.values[[i]], row.names = rep("", nrow(p.values[[i]])))
            }
            message("...done.")
          }
        }
        else if (xmod == "EC" && ymod == "BC"){
          if(!is.null(X$EnvSuit[[i]]) &&  !is.null(Y$EnvSuit[[i]])){
            message("\n\t...between environmental niches....", appendLF = F)
            ncomp <-  niche_comparison(X$EnvSuit[[i]], Y$EnvSuit[[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                       quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                       cor = cor, rep = rep, rand = rand, alternative = alternative)
            np[[i]] <- c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                       ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                              comparison = "Env2Env")

            if (rnd.test){
              simList[[i]] <- cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                         ER = ncomp$ER$sim, ncomp$USE$sim,
                                         comparison = "Env2Env")

              p.values[[i]] <- cbind(alternative = alternative, D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue,
                                          CS = ncomp$CS$pvalue, ER = ncomp$ER$pvalue,
                                          ncomp$USE$pvalue, comparison = "Env2Env")
            }
            message("...done.")
          }
          if(!is.null(X$EcoSuit[[i]]) &&  !is.null(Y$EnvSuit[[i]])){
            message("\n\t...between ecological niche and environmental niche....", appendLF = F)
            ncomp <-  niche_comparison(X$EcoSuit[[i]], Y$EnvSuit[[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                       quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                       cor = cor, rep = rep, rand = rand, alternative = alternative)

            np[[i]] <- rbind(np[[i]], c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                                           ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                                                  comparison = "Eco2Env"))
            if (rnd.test){
              simList[[i]] <- rbind(simList[[i]], cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                                                  ER = ncomp$ER$sim, ncomp$USE$sim,
                                                                  comparison = "Eco2Env"))

              p.values[[i]] <- rbind(p.values[[i]], cbind(alternative = alternative,
                                                                    D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue, CS = ncomp$CS$pvalue,
                                                                    ER = ncomp$ER$pvalue, ncomp$USE$pvalue,
                                                                    comparison = "Eco2Env"))
              p.values[[i]] <- as.data.frame(p.values[[i]], row.names = rep("", nrow(p.values[[i]])))
            }
            message("...done.")
          }
          if(!is.null(X$EcoAcc[[i]]) &&  !is.null(Y$EnvAcc[[i]])){
            message("\n\t...between ecological accessibillity and environmental accessibillity....", appendLF = F)
            ncomp <-  niche_comparison(X$EcoAcc[[i]], Y$EnvAcc[[i]], centroid.w = centroid.w, rnd.test = rnd.test,
                                       quantile = quantile, np.type = np.type,  np.metric = np.metric,
                                       cor = cor, rep = rep, rand = rand, alternative = alternative)

            np[[i]] <- rbind(np[[i]], c(unlist(c(D = ncomp$D$obs, breadth.diff = ncomp$breadth.diff$obs, CS = ncomp$CS$obs,
                                                           ER = ncomp$ER$obs, ncomp$CS$weights, ncomp$USE$obs)),
                                                  comparison = "Acc2Acc"))
            if (rnd.test){
              simList[[i]] <- rbind(simList[[i]], cbind(D  = ncomp$D$sim, breadth.diff = ncomp$breadth.diff$sim, CS = ncomp$CS$sim,
                                                                  ER = ncomp$ER$sim, ncomp$USE$sim,
                                                                  comparison = "Acc2Acc"))

              p.values[[i]] <- rbind(p.values[[i]], cbind(alternative = alternative,
                                                                    D = ncomp$D$pvalue, breadth.diff = ncomp$breadth.diff$pvalue, CS = ncomp$CS$pvalue,
                                                                    ER = ncomp$ER$pvalue, ncomp$USE$pvalue,
                                                                    comparison = "Acc2Acc"))
              p.values[[i]] <- as.data.frame(p.values[[i]], row.names = rep("", nrow(p.values[[i]])))
            }
            message("...done.")
          }
        }
      }
      else{
        message(paste("...nothing to analyze."))
      }
    }
    if (rnd.test){
      out$np <- plyr::ldply(np, .id = "species")
      out$sim <- simList
      out$pvalues <- p.values
    } else {
      out <- plyr::ldply(np, .id = "species")
    }
  }
  message("DONE")
  return(out)
}
