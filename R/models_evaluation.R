#' @title MODEL EVALUATION FUNCTION
#'
#' @description full models evaluation
#'
#' @param Pred Data frame with predicted niche models. Alternatively it can be an object class NINA, in which case \code{Obs} and  \code{predictors} is omitted.
#' @param Obs Occurrence dataset
#' @param predictors Environmental predictors
#' @param spsNames species names to evaluate
#' @param th threshol to perform cut off
#' @param ras raster to constrain pseudoabsences sampling
#' @param int.matrix interaction matrix between species and interactors
#' @param res spatial resolution
#' @param plot Boolean to whether plot the evaluation
#' @param sample.pseudoabsences Boolean to whether sample pseudo-absences
#' @param rep number of randomzation tests
#' @param best.th method to select the best thresholt. Default is "similarity"
#'
#' @return List
#'
#'@details Returns an error if \code{filename} does not exist.
#'
#' @examples
#' \dontrun{
#' EN = EN_model(env_data, occ_data1, cluster = "env", n.clus = 5)
#' eval = models_evaluation(EN$pred.dis, EN$obs, env_data)
#' print(eval)
#' }
#'
#' @importFrom raster stack rasterFromXYZ maxValue
#' @importFrom stats na.exclude
#' @importFrom ecospat ecospat.sample.envar
#' @importFrom tidyr gather
#' @import ggplot2
#'
#' @export
models_evaluation <- function(Pred, Obs, predictors, spsNames = NULL, th = NULL, ras = NULL, int.matrix = NULL, sample.pseudoabsences = TRUE,
                              res = NULL, plot = TRUE, rep = 100,
                              best.th = c("accuracy", "similarity") ){

  if (class(Pred) == "NINA"){
    Obs = Pred$obs
    predictors = Pred$env.scores
    pred.stack = Pred$maps
    if (Pred$type %in% c("BC", "EC")){
      if (Pred$type == "BC"){
        BioCons = sapply(reverse_list(Pred$z.mod), function(x) niche_to_dis(predictors, x, cluster = Pred$clus, cor = FALSE)[,3])
      }
      if (Pred$type == "EC"){
        BioCons = sapply(reverse_list(Pred$t.mod), function(x) niche_to_dis(predictors, x, cluster = Pred$clus, cor = FALSE)[,3])
      }
      BioCons[is.na(BioCons)] = 0
      BioCons = cbind(predictors[,1:2], BioCons)
      ras = raster_projection(BioCons, ras = Pred$maps[[1]])
      int.matrix = diag(nrow = raster::nlayers(Pred$maps))
      rownames(int.matrix) <- colnames(int.matrix) <- names(Pred$maps)
    }
    Pred = Pred$pred.dis
  }
  best.th = best.th[1]
  # check Obs argument
  if (ncol(Obs) == 3) { Obs$PA = 1  }
  message("Observations ... OK")
  # check predictors argument
  if (is.data.frame(predictors)){
    pred.var = colnames(predictors[,-c(1:2)]) # environmental variables
    pred.stack = raster_projection(predictors)
  }
  if (any(class(predictors) %in% c("raster", "RasterBrick", "RasterStack"))){
    pred.stack = predictors
    predictors <- na.exclude(raster::as.data.frame(pred.stack, xy = T))
    pred.var = colnames(predictors[,-c(1:2)]) # environmental variables
  }
  if (is.null(res)){
    res = max(raster::res(pred.stack))
  }
  message("Environmental predictors ... OK")
  # check Pred argument
  if (any(class(Pred) %in% c("raster", "RasterBrick", "RasterStack"))){
    Pred <- cbind(predictors[,1:2], raster::extract(Pred, predictors[,1:2]))
    #Pred <- ecospat.sample.envar(dfsp=predictors,colspxy=1:2,colspkept=1:2,dfvar=Pred,colvarxy=1:2,colvar= 3:ncol(Pred) ,resolution= res)
  }
  # check spsNames argument
  if (is.null(spsNames)){ spsNames = colnames(Pred)[-c(1:2)]}
  spsObs = levels(Obs[,3])
  if(any(!spsNames %in% spsObs)) {
    missing.sps <- spsNames[which(!spsNames %in% spsObs)]
    message("Following predicted species have no observations provided: ")
    print(missing.sps)
  }
  message("Model predictions ... OK")
  Obs <- Obs[Obs[,3] %in% spsNames,]
  if (sample.pseudoabsences == TRUE){
    message("Sampling pseudo_absences... ")
    Abs.samp <- sample_pseudoabsences(Obs,  predictors , spsNames = spsNames, ras = ras, int.matrix = int.matrix,  res = res)
    Obs <- rbind(Abs.samp$Presences, Abs.samp$Absences)
    occ.tab <- Abs.samp$tab
  }
  else{
    Obs <- na.exclude(ecospat.sample.envar(dfsp=predictors,colspxy=1:2,colspkept=1:2,dfvar=Obs,colvarxy=1:2,colvar= 3:4 ,resolution= res))
    Obs[,3] = as.factor(Obs[,3])
    levels(Obs[,3]) <- spsObs
    occ.tab <- t(sapply(spsObs, function(x) cbind(sum(Obs[,3] == x & Obs[,4] == 1),
                                                  sum(Obs$species == x & Obs[,4] == 0))))
    colnames(occ.tab) = c("n.occurrences", "n.absences")
  }
  tab <- NULL
  confusion <- NULL
  n <- data.frame(matrix(NA, length(spsNames), 2), row.names = spsNames)
  colnames(n) = c("np", "na")
  threshold <- NULL
  message("Performing models evaluation...")
  for (i in spsNames){
    message(paste("\t...Evaluating", i, "niche model..."))
    Obs.sp <- Obs[Obs$species == i,]
    pred.i <- which(colnames(Pred) %in% i)
    Pred.sp <- ecospat.sample.envar(dfsp=Obs.sp,colspxy=1:2,colspkept=1:2,dfvar=Pred,colvarxy=1:2,colvar= pred.i ,resolution= res)
    Obs. <- Obs.sp[,4]
    Fit. <- Pred.sp[,3]
    Fit.[is.na(Fit.)] = 0
    eval = evaluate_model(Fit., Obs., best.th = best.th, rep = rep, th = th, main = i, plot = plot)
    tab = rbind(tab, eval$tab)
    confusion <- rbind(confusion, cbind(eval$confusion, species = i))
    n[i,] <- eval$n
    threshold <-  rbind(threshold, eval$threshold)
  }
  tab <- as.data.frame(tab, row.names = spsNames)
  threshold <- as.data.frame(threshold, row.names = spsNames)
  if (plot == TRUE){
    z <- tidyr::gather(tab, "test", "value", -c(2,12) )
    z$test <- factor(z$test,levels = c("Pearson's correlation",  "Jaccard Similarity",
                                       "TPR" ,  "TNR", "TSS","ACC", "AUC", "kappa", "PPV", "NPV"))
    p <- ggplot(z, aes_string(x = "test", y = "value", group = "test")) +
      geom_boxplot() +
      ylim(0,1) +
      scale_x_discrete(labels = gsub('\\s','\n',levels(z$test))) +
      labs(title= "All models" ,x="", y = "Score") + theme_classic() +
      theme(axis.title.x = element_text( size=14, ),
            axis.title.y = element_text( size=14, ),
            axis.text.x =element_text( size=12, ),
            axis.text.y = element_text(size=12,))
    print(p)
  }
  out <- list(n = n, tab = tab, threshold = threshold, confusion = confusion, cases = occ.tab, type = "eval")
  attr(out, "class") <- "NINA"
  message("Models evaluation performed.")
  return(out)
}
