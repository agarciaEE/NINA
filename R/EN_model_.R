#' @title EN MODELLING internal FUNCTION
#'
#' @description Samples pseudo-absences from a environmental extent based on euclidean distances to the observations in a environmental PCA
#'
#' @param env Environmental data frame
#' @param occ Occurrence dataset
#' @param sample.pseudoabsences  Boolean to indicate if pseudoabsences are to be sampled
#' @param extrapolate.niche Boolean whether to allow niche extrapolation if needed
#' @param nstart Burn up start for clustering estimation
#' @param k.max Maximum number of clusters
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
#' @param crs CRS object or a character string describing a projection and datum in PROJ.4 format
#'
#' @return List of elements
#'
#' @details Returns an error if \code{filename} does not exist.
#'
#' @examples
#' \dontrun{
#' EN<- EN_model_(env_data, occ_data2, cluster = "env", n.clus = 5)
#' }
#'
#' @importFrom raster stack rasterFromXYZ maxValue res crs
#' @importFrom stats na.exclude
#' @importFrom ade4 dudi.pca
#' @importFrom ecospat ecospat.sample.envar ecospat.grid.clim.dyn
#' @importFrom tidyr spread
#' @importFrom plyr ldply
#'
#' @keywords internal
#' @noRd
#'
EN_model_ <- function(env, occ, res = NULL, sample.pseudoabsences = TRUE, crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0",
                                extrapolate.niche = FALSE, nstart = 25, k.max = NULL, B = 100,
                                combine.clusters = FALSE, cluster = NULL, n.clus = NULL, R = 100,
                                eval = FALSE, split.data = FALSE, split.percentage = 0.25, split.method  = c("kmeans", "Euclidean")){

  split.method = split.method[1]
  # check if there are enough occurrences per species
  ###################
  rm.sp = names(table(occ[,3]))[table(occ[,3]) < 4]
  spsNames = names(table(occ[,3]))[!names(table(occ[,3])) %in% rm.sp]
  occ = occ[occ[,3] %in% spsNames,]
  if (ncol(occ) == 3){
    occ$PA = 1
  }
  if (is.character(cluster) && cluster == "obs") { obs = occ }
  if (split.data == TRUE){
    split.df <- split_data(occ, method = split.method,  split.percentage = split.percentage)
    occ <- split.df$occ.train
    occ.test <- split.df$occ.test
  }
  occ = occ[!duplicated(occ),]
  occ.df = tidyr::spread(occ, 3, 4)
  occ.df[is.na(occ.df)] = 0
  if (length(rm.sp) > 0){
    warning(length(rm.sp), " species removed from the analysis... Not enough observations.")
    message("Following species were dropped out of the analyses")
    print(rm.sp)
  }
  message("Occurrence dataset ... OK")
  ###################
  # check env argument
  ###################
  if (is.data.frame(env)){
    env = stats::na.exclude(env)
    env.var = colnames(env[,-c(1:2)]) # environmental variables
    env.stack = raster_projection(env, crs = crs)
  }
  if (any(class(env) %in% c("raster", "RasterBrick", "RasterStack"))){
    env.stack = env
    crs = raster::crs(env.stack)
    env <- stats::na.exclude(raster::as.data.frame(env.stack, xy = T))
    env.var = colnames(env[,-c(1:2)]) # environmental variables
  }
  message("Environmental dataset ... OK")
  ###################
  # check res argument
  ###################
  if (is.null(res)){
    res = max(raster::res(env.stack))
  }
  ###################
  ### inner elements
  output = list()
  sps.scores = NULL
  z.mod = list()
  mod.Val = list()
  fail.m <- NULL
  ### error messages
  ### function
  message("\t- Conforming environmental space...")
  pca.cal <- dudi.pca(stats::na.exclude(env[,env.var]), center = T, scale = T, scannf = F, nf = 2) # the pca is calibrated on all the sites of the study area
  env.scores = cbind(env[,1:2], pca.cal$li) # environmental PCA scores
  occ.df <- ecospat::ecospat.sample.envar(dfsp=env,colspxy=1:2,colspkept=1:2,dfvar=occ.df,colvarxy=1:2,colvar= 3:ncol(occ.df) ,resolution=res) # match occ.df coordinates to env coordinates
  occ.df[is.na(occ.df)] = 0
  ## ENN method
  ### Whole niche
  ###################
  if (is.null(cluster)){
    message("\t- Carrying out species EN models...")
    sps.scores = lapply(spsNames, function(i) env.scores[occ.df[,i] != 0,])
    names(sps.scores) <-  spsNames
    for(i in spsNames) {
      message("\t\t- Modelling ", i, " Environmental Niche...")
      sp.scores <- sps.scores[[i]][,3:4]
      if (nrow(sp.scores) > 4) {
        z.mod[[i]] <- ecospat.grid.clim.dyn(env.scores[,3:4],env.scores[,3:4],sp.scores,R)
        mod.Val[[i]] <- cbind(env.scores[,1:2], response = raster::extract(z.mod[[i]]$z, z.mod[[i]]$glob))
      }
      else {
        fail.m <- c(fail.m, i)
      }
    }
    if (length(mod.Val) != 0) {
      tab  = table(names(mod.Val))
      mod.Val <- ldply(mod.Val, data.frame, .id = "species")
      sps.scores <- ldply(sps.scores, data.frame, .id = "species")
      output$tab <-tab
      if (length(fail.m) > 0){
        output$fail <- fail.m
        warning("Model failed to predict the following species due to lack of observations", immediate. = T)
        print(fail.m)
      }
    }
    else {
      stop("Modelling has failed to predict any species selected. \nProbable cause not enough scores inside the niche space.")
    }
  }
  ###################
  ### Cluster niche
  ###################
  if (!is.null(cluster)) {
    if (is.character(cluster) && cluster == "obs") {
      message("\t- Clustering species observations...")
      obs = spread(obs, 3, 4)
      obs <- ecospat.sample.envar(dfsp=env,colspxy=1:2,colspkept=1:2,dfvar=obs,colvarxy=1:2,colvar= 3:ncol(obs) ,resolution=res) # match occ.df coordinates to env coordinates
      obs[is.na(obs)] = 0
      clus.df = cluster_regions(obs, plot = TRUE, n.clus = n.clus, nstart = nstart, K.max = k.max, B = B)
      clus.df[,3] = as.factor(clus.df[,3])
      regions <- LETTERS[1:length(levels(clus.df[,3]))]
      levels(clus.df[,3]) = regions
    }
    if (is.character(cluster) && cluster == "env") {
      message("\t- Clustering environmental variables...")
      clus.df = cluster_regions(env, plot = TRUE, n.clus = n.clus, nstart = nstart, K.max = k.max, B = B)
      clus.df[,3] = as.factor(clus.df[,3])
      regions <- LETTERS[1:length(levels(clus.df[,3]))]
      levels(clus.df[,3]) = regions
    }
    if (is.data.frame(cluster)) {
      message("\t- Clustering based on provided regions...")
      cluster[,3] <- as.factor(cluster[,3])
      regions <- levels(cluster[,3])
      clus.df = na.exclude(ecospat.sample.envar(dfsp=env,colspxy=1:2,colspkept=1:2,dfvar=cluster,colvarxy=1:2,colvar= 3 ,resolution=res))
      clus.df[,3] = as.factor(clus.df[,3])
      levels(clus.df[,3]) <- regions
    }
    if (is.vector(cluster) && length(cluster) == nrow(env)) {
      message("\t- Clustering based on provided regions...")
      clus.df <- cbind(env[,1:2], cluster)
      clus.df[,3] = as.factor(clus.df[,3])
      regions <-levels(clus.df[,3])
    }
    message("\t- Regions:")
    print(regions)
    fail.m <- NULL
    for (e in regions) {
      message("\t- Carrying out species EN models in region ", e, "...")
      reg = rownames(clus.df)[which(clus.df[,3] == e)]
      occ.subset <- occ.df[reg,]
      # subseting environmental grid allows to have better resolution on each region
      env.scores.subset <- env.scores[reg,]
      if (nrow(occ.subset) > 0) {
        sps.subset <- names(which(sapply(spsNames, function(x) sum(occ.subset[,x]) > 0)))
        z.mod[[e]] = list()
        mod.Val[[e]] <- list()
        #fail.m[[e]] <- list()
        for(i in sps.subset) {
          message("\t\t- Modelling ", i, " Environmental Niche...")
          sp.scores <- merge(occ.df[occ.df[reg,i] != 0,1:2], env.scores.subset, by= 1:2)
          if (nrow(sp.scores) <= 4){
            if (extrapolate.niche == TRUE){
              sp.scores <- merge(occ.df[occ.df[,i] != 0,1:2], env.scores, by= 1:2)
              #sp.scores <- merge(sp.scores, env.scores.subset, by= 3:4)
              #sp.scores <- sp.scores[,c(5:6,1:2)]
              #colnames(sp.scores) <- colnames(env.scores.subset)
            }
          }
          if (nrow(sp.scores) > 4) {
            message("\t\t- Estimating ", i, " niche response in region ", e, "...")
            z.mod[[e]][[i]] <- ecospat.grid.clim.dyn(env.scores.subset[,3:4],env.scores.subset[,3:4],sp.scores[,3:4],R )
            mod.Val[[e]][[i]] <- cbind(env.scores.subset[,1:2], response = raster::extract(z.mod[[e]][[i]]$z, z.mod[[e]][[i]]$glob))
            sps.scores <- rbind(sps.scores, cbind(region = e, species = i, sp.scores))
          }
          else {
            warning("Not enough observations of ", i , " in region ", e, ".", immediate. = T)
            fail.m <- rbind(fail.m, cbind( region = e ,  species = i))
          }
        }
        mod.Val[[e]] <- ldply(mod.Val[[e]], data.frame, .id = "species")
      }
      else {
        warning("There are no observations in region ", e, ".", immediate. = T)
      }
    }
    z.mod <-  z.mod[unlist(lapply(z.mod, length) != 0)]
    z.mod <- lapply(z.mod, function(i) i[unlist(lapply(i, function(x) !is.null(x)))])
    z.mod <- lapply(z.mod, function(i) i[unlist(lapply(i, function(x) !is.na(maxValue(x$z))))])
    z.mod <- lapply(z.mod, function(i) i[unlist(lapply(i, function(x) length(x) > 0))])
    if (length(mod.Val) != 0) {
      if (combine.clusters == TRUE){
        message("\t- Assembling regions into global model...")
        output$z.mod.global = combine_regions(z.mod, R = R)
      }
      tab = cbind(ldply(sapply(z.mod, function(x) names(x)), data.frame, .id = "region"), P = 1)
      tab =  spread(tab, "region", "P")
      rownames(tab) <- tab[,1]
      tab <- tab[,-1]
      tab[is.na(tab)] = 0
      output$tab <-t(tab)
      mod.Val <- ldply(mod.Val, data.frame, .id = "region")
      mod.Val = mod.Val[,-1]
      output$clus <- clus.df
      if (nrow(fail.m) > 0){
        #fail.m <- plyr::ldply(fail.m, data.frame, .id = "region")
        #colnames(fail.m)[2] = "species"
        output$fail <- fail.m
        warning("Model failed to predict the following species due to lack of observations", immediate. = T)
        print(fail.m)
      }
      else {
        message("All species have been modelled correctly.")
      }
    }
    else {
      stop("Modelling has failed to predict any species selected. \nProbable cause not enough scores inside the niche space.")
    }
  }
  ###############
  mod.Val <- spread(mod.Val, "species", "response")
  mod.Val[is.na(mod.Val)] = 0
  mod.Val[,-c(1:2)] <- apply(mod.Val[,-c(1:2)], 2, function(i) i/max(i, na.rm = T))
  output$pca <- pca.cal
  output$sp.scores <- sps.scores
  output$env.scores <- env.scores
  output$z.mod <- z.mod
  output$pred.dis <- mod.Val
  if (split.data == TRUE){output$obs <- list(train = occ, test = occ.test)} else { output$obs <- occ}
  ## Model evaluation
  ###############
  if (eval == TRUE){
    message("\t- Carrying out models evaluations...")
    output$eval <- models_evaluation(mod.Val, if (split.data == TRUE){occ.test} else {occ}, predictors = env, sample.pseudoabsences = sample.pseudoabsences, res = res)
  }
  ###############
  output$maps  = raster_projection(mod.Val, ras = env.stack[[1]], crs = crs)
  output$predictors = env.var
  output$crs = crs
  output$type = "EN"
  attr(output, "class") <- "NINA"
  message("Species EN models succesfully completed!")
  return(output)
}


