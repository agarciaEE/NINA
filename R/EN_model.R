#' @title  EN MODELLING FUNCTION
#'
#' @description Samples pseudo-absences from a environmental extent based on euclidean distances to the observations in a environmental PCA
#'
#' @param env Environmental data frame
#' @param occ Occurrence dataset
#' @param sample.pseudoabsences  Boolean to indicate if pseudoabsences are to be sampled
#' @param extrapolate.niche Boolean wheter to allow niche extrapolation if needed
#' @param nstart Burn up start for clustering estimation
#' @param k.max Maximum number of cluseters
#' @param B Number of runs
#' @param combine.clusters Boolean whether to combine regions into global models
#' @param cluster data frame with clustered regions as levels
#' @param n.clus number of clusters to perform
#' @param R niche space grid
#' @param eval Boolean whether to evaluate the model
#' @param split.data Boolean whether to split data in case of model bootstrapping
#' @param split.percentage Split percentage, Default is 0.25
#' @param split.method String indicating the method to split occurrence data. Default is kmeans
#' @param res grid resolution of the spatial extent
#' @param path Directory path of partition models
#' @param project.name String indicating the project name. Default is "NINA_EN"
#' @param save.bootstraps Boolean whether to save partitions
#' @param bootstraps Integer indicating number of partitions to perform. Default is 1
#' @param assemble.models Boolean whether to assemble all partitions into ensemble model
#' @param assemble.method String indicating the evaluation parameter to weight and compute the assembling. Default is "ACC"
#' @param crs CRS object or a character string describing a projection and datum in PROJ.4 format
#'
#' @return List of elements
#'
#'@details Returns an error if \code{filename} does not exist.
#'
#' @examples
#' \dontrun{
#' EN<- EN_model(env_data, occ_data1, boot)
#' }
#' \dontrun{
#' EN<- EN_model(env_data, occ_data2, cluster = "env", n.clus = 5)
#' }
#'
#' @importFrom utils write.table capture.output
#' @importFrom stats na.exclude
#'
#' @export
EN_model <- function(env, occ, res = NULL, path = "./", project.name	= "NINA_EN",
                     nstart = 25, k.max = NULL, B = 100, crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0",
                     extrapolate.niche = FALSE, save.bootstraps = F, save.model = F,
                     combine.clusters = FALSE, cluster = NULL, n.clus = NULL, R = 100, sample.pseudoabsences = TRUE,
                     eval = FALSE, split.data = FALSE, split.percentage = 0.25, split.method  = c("kmeans", "Euclidean"),
                     bootstraps = 1, assemble.models = TRUE, assemble.method = c("ACC", "Jaccard Similarity", "TSS", "AUC", "kappa")){

  split.method = split.method[1]
  #write.table(capture.output(as.list(match.call(expand.dots=FALSE))), file = paste0(path,project.name, ".args.txt"))
  ### if bootstrap...
  ##################
  if (bootstraps > 1){
    if (save.bootstraps || save.model){
      pn = project.name
      i = 1
      while (dir.exists(file.path(path, pn))){
        pn = paste0(project.name,i)
        i = i+1
      }
      dir.create(paste0(file.path(path, pn),"/"), showWarnings = FALSE)
      path <- paste0(file.path(path, pn),"/")
    }
    modelsList <- list()
    for (n in 1:bootstraps){
      message("Carrying out EN model bootstrap n", n, "...")
      EN <- EN_model_(env, occ, res = res, nstart = nstart, k.max = k.max, B = B,
                      extrapolate.niche = extrapolate.niche, sample.pseudoabsences = sample.pseudoabsences,
                      combine.clusters = combine.clusters, cluster = cluster, n.clus = n.clus, R = R,
                      eval = eval, split.data = split.data, split.percentage = split.percentage, split.method = split.method)
      pca = EN$pca
      crs = EN$crs
      if (!is.null(cluster)){
        cluster = EN$clus
      }
      if (save.bootstraps){
        saveRDS(EN, file = paste0(path, project.name , "_bootstrap", n, ".RDS"))
      }
      else {
        modelsList[[n]]<- EN
      }
      rm(EN)
    }
    if (assemble.models == TRUE){
      message("\t Assembling bootstrap models...")
      if (save.bootstraps){
        modelsList = list.files(path, pattern = "\\.RDS$")
      }
      EN = assemble_snm_models(modelsList, method = assemble.method, modelspath = path)
      if (combine.clusters){
        EN$z.mod.global <- combine_regions(EN$z.mod, R = R)
      }
      if (eval){
        message("\t Carrying out assembled model evaluation...")
        EN$eval <- models_evaluation(EN$pred.dis, occ,
                                     predictors = env, sample.pseudoabsences = sample.pseudoabsences,  res = res)
      }
      EN$pca = pca
      EN$obs = occ
      if (save.model){
        saveRDS(EN, file = paste0(path, project.name ,"_ensemble.RDS"))
      }
    }
  }
  ### Single model
  ##################
  if (bootstraps == 1){
    message("Carrying out unique EN model...")
    EN <- EN_model_(env, occ, res = res, nstart = nstart, k.max = k.max, B = B,
                    extrapolate.niche = extrapolate.niche, sample.pseudoabsences = sample.pseudoabsences,
                    combine.clusters = combine.clusters, cluster = cluster, n.clus = n.clus, R = R,
                    eval = eval, split.data = split.data , split.percentage = split.percentage, split.method = split.method)
    EN$obs = occ
    if (save.model){
      saveRDS(EN, file = paste0(path, project.name ,".RDS"))
    }
  }
  return(EN)
}
