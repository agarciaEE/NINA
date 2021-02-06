#' @title MODEL EVALUATION FUNCTION
#'
#' @description full models evaluation
#'
#' @param Pred Predicted niche models
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
#' accident_2015 <- fars_read("Project/data/accident_2015.csv.bz2")
#' }
#'
#' @importFrom raster stack rasterFromXYZ maxValue
#' @importFrom stats na.exclude
#' @importFrom ecospat ecospat.sample.envar
#' @import ggplot2
#'
#' @export
models_evaluation <- function(Pred, Obs, predictors, spsNames = NULL, th = NULL, ras = NULL, int.matrix = NULL, sample.pseudoabsences = TRUE,
                              res = NULL, plot = TRUE, rep = 100,
                              best.th = c("accuracy", "similarity") ){

  best.th = best.th[1]
  # check Obs argument
  if (ncol(Obs) == 3) { Obs$PA = 1  }
  # check predictors argument
  if (is.data.frame(predictors)){
    pred.var = colnames(predictors[,-c(1:2)]) # environmental variables
    pred.stack = stack(sapply(pred.var, function(x) rasterFromXYZ(cbind(predictors[,1:2], predictors[,x]))))
  }
  if (class(predictors) %in% c("raster", "RasterBrick", "RasterStack")){
    pred.stack = predictors
    predictors <- na.exclude(raster::as.data.frame(pred.stack, xy = T))
    pred.var = colnames(predictors[,-c(1:2)]) # environmental variables
  }
  if (is.null(res)){
    res = max(res(pred.stack))
  }
  # check Pred argument
  if (class(Pred) %in% c("raster", "RasterBrick", "RasterStack")){
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
  Obs <- Obs[Obs[,3] %in% spsNames,]
  if (sample.pseudoabsences == TRUE){
    Abs.samp <- sample_pseudoabsences(Obs,  predictors , ras = ras, int.matrix = int.matrix,  res = res)
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
  n <- NULL
  threshold <- NULL
  for (i in spsNames){
    Obs.sp <- Obs[Obs$species == i,]
    pred.i <- which(colnames(Pred) %in% i)
    Pred.sp <- ecospat.sample.envar(dfsp=Obs.sp,colspxy=1:2,colspkept=1:2,dfvar=Pred,colvarxy=1:2,colvar= pred.i ,resolution= res)
    Obs. <- Obs.sp[,4]
    Fit. <- Pred.sp[,3]
    Fit.[is.na(Fit.)] = 0
    eval = evaluate_model(Fit., Obs., best.th = best.th, rep = rep, th = th, main = i, plot = plot)
    tab = rbind(tab, eval$tab)
    confusion <- rbind(confusion, cbind(eval$confusion, species = i))
    n <- rbind(n, eval$n)
    threshold <-  rbind(threshold, eval$threshold)
  }
  tab <- as.data.frame(tab, row.names = spsNames)
  threshold <- as.data.frame(threshold, row.names = spsNames)
  if (plot == TRUE){
    z <- tidyr::gather(tab, "test", "value", -c(2,12) )
    z$test <- factor(z$test,levels = c("Pearson's correlation",  "Jaccard Similarity",
                                       "TPR" ,  "TNR", "TSS","ACC", "AUC", "kappa", "PPV", "NPV"))
    p <- ggplot(z, aes(x = "test", y = "value", group = "test")) +
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
  out <- list(n = n, tab = tab, threshold = threshold, confusion = confusion, cases = occ.tab)
  return(out)
}
