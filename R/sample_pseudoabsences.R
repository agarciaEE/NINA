#' @title Sampling pseudo-absences
#'
#' @description Samples pseudo-absences from a environmental extent based on euclidean distances to the observations in a environmental PCA
#'
#' @param Obs a n x 3 data.frame class object indicating longitud, latitud and species in the three first columns
#' @param predictors a n x m data.frame class object indicating longitud, latitud and species in the two first columns. The rest of the columns must be the environmental variables from which the PCA will be computed. Alternatively a raster stack with environmental variables or a dudi.pca object. If dudi.pca, res argument must be provided, otherwhise it will be set to 1.
#' @param spsNames a string vector with selected species names to subset Obs. Default NULL will select all species present in the third column of paramteter Obs
#' @param th a numeric threshold to estimate pseudo-absence distances from observations environmental centroid. Default is 0.95
#' @param ras raster or stacked raster with the ecological contrains that will weight pseudo-absence probabilities. Default NULL
#' @param int.matrix m x o matrix defining the effect of the interactors as columns, provided in ras parameter, on each species to sample pseudo-absences, as rows.
#' @param res resolution of the environmental grid. Default NULL means that will be computed automatically.
#' @param plot Boolean to indicate if computed distances is to be plotted
#'
#' @return List
#'
#' @details Returns an error if \code{filename} does not exist.
#'
#' @examples
#' \dontrun{
#' EN = EN_model(env_data, occ_data1, cluster = "env", n.clus = 5)
#' PseudoA <- sample_pseudoabsences(EN)
#' }
#'
#' @importFrom ade4 dudi.pca
#' @importFrom stats na.exclude quantile
#' @importFrom raster stack rasterFromXYZ
#' @importFrom ecospat ecospat.sample.envar
#' @importFrom graphics points
#'
#' @export
sample_pseudoabsences <- function(Obs, predictors, spsNames = NULL, th = 0.95,
                                       ras = NULL, int.matrix, res = NULL, plot = F){

  # check spsNames argument
  if (is.null(spsNames)){spsNames = unique(Obs[,3])}
  # check predictors argument
  if (class(predictors) == "NINA"){
    pca_scores = predictors[,3:4]
    predictors = predictors[,1:2]
    if (is.null(res)){
      warning("res argument is not given and cannot be computed from argument predictors. It will be set to value 1")
      res = 1
    }
  }
  else if (is.data.frame(predictors)){
    pred.var = colnames(predictors[,-c(1:2)]) # environmental variables
    pred.stack = raster_projection(predictors)
    pca_scores<- ade4::dudi.pca(predictors[,pred.var], center = T, scale = T, scannf = F, nf = 2)$li # the pca is calibrated on all the sites of the study area
    if (is.null(res)){
      res = max(res(pred.stack))
    }
  }
  else if (any(class(predictors) %in% c("raster", "RasterBrick", "RasterStack"))){
    pred.stack = predictors
    predictors <- na.exclude(raster::as.data.frame(pred.stack, xy = T))
    pred.var = colnames(predictors[,-c(1:2)]) # environmental variables
    pca_scores <- ade4::dudi.pca(predictors[,pred.var], center = T, scale = T, scannf = F, nf = 2)$li # the pca is calibrated on all the sites of the study area
    if (is.null(res)){
      res = max(res(pred.stack))
    }
  }
  if (!is.null(ras)){
    if(is.null(int.matrix)){
      stop("Argument 'int.matrix' not provided. ")
    }
    else{
      if (is.data.frame(ras)){
        ras.var = colnames(ras[,-c(1:2)]) # environmental variables
        ras = raster_projection(ras)
      }
      if (any(class(ras) %in% c("raster", "RasterBrick", "RasterStack"))){
        #  ras.stack = ras
        #  ras <- na.exclude(raster::as.data.frame(ras.stack, xy = T))
        ras.var = names(ras) # environmental variables
      }
      if(all(spsNames %in% colnames(int.matrix)) && all(rownames(int.matrix) %in% ras.var)){
        int.matrix = t(int.matrix)
      }
      if(!all(spsNames %in% rownames(int.matrix))){
        stop("Argument 'int.matrix' is missing some of the species of study.")
      }
      if(!all(colnames(int.matrix) %in% ras.var)){
        stop("Argument 'ras' is missing some of the species interactions.")
      }
    }
  }
  out = list()
  occ.P <- NULL
  occ.A <- NULL

  def <- data.frame(matrix(NA, nrow = length(spsNames), ncol = 3))
  colnames(def) = c("occ/pred", "n.occurrences", "n.absences")
  rownames(def) = spsNames
  for (i in 1:length(spsNames)) {
    sp = spsNames[i]
    occ <- Obs[Obs$species == sp, 1:3]
    occ.inn <-na.exclude(ecospat::ecospat.sample.envar(dfsp=predictors,colspxy=1:2,colspkept=1:2,dfvar=occ,colvarxy=1:2,colvar= 3 ,resolution=res))
    occ.inn <- rownames(occ.inn)
    occ.P = rbind(occ.P, cbind(predictors[occ.inn,1:2], species = sp, PA = 1))
    if (!is.null(ras)){
      Xvar <- colnames(int.matrix)[which(int.matrix[sp,] != 0)]
      nvar<-length(Xvar)
      if (nvar > 1){
        m.ras <- sum(ras[[Xvar]], na.rm = T)
      } else {m.ras <- ras[[Xvar]]}
      m.ras <- cbind(predictors[,1:2], P = raster::extract(m.ras, predictors[,1:2]))
      m.ras = m.ras[m.ras[,3] != 0, 1:3]
      ras.inn <- rownames(m.ras)
    }
    dist.mat <- sapply(1:length(occ.inn), function(n) sqrt((pca_scores[,1] - pca_scores[occ.inn[n],1])^2 + (pca_scores[,2] - pca_scores[occ.inn[n],2])^2))
    dist.p = rowSums(dist.mat)
    Dmax = max(dist.p)
    dist.p = dist.p/Dmax
    rownames(dist.mat) = rownames(predictors)
    if(plot){ plot(rasterFromXYZ(cbind(predictors[,1:2], dist.p)))}
    dist.mat = as.data.frame(dist.mat, row.names = rownames(predictors))
    p.th = quantile(rowSums(dist.mat[occ.inn,])/Dmax, th, na.rm = T)
    inn_coord = rownames(dist.mat)[which(dist.p <= p.th)]
    if (!is.null(ras)){
      inn_coord = inn_coord[which(inn_coord %in% ras.inn)]
    }
    inn_coord = c(inn_coord, occ.inn)[!duplicated(c(inn_coord, occ.inn))]
    pred.inn <- predictors[inn_coord,1:2]
    def[i,1] = nrow(pred.inn)/nrow(predictors)
    def[i,2] = length(occ.inn)
    dist.p[which(rownames(predictors) %in% inn_coord)] <- 0
    pred.out <- predictors[!rownames(predictors) %in% occ.inn,1:2]
    dist.p <- dist.p[which(!rownames(predictors) %in% occ.inn)]
    if (length(occ.inn)/nrow(pred.inn) == 1) {
      def[i,3] = round(length(occ.inn)/nrow(pred.inn) * nrow(pred.out),0)
      abs =  as.data.frame(pred.out[sample(nrow(pred.out), def[i,3], replace = T, prob = dist.p),1:2])
      abs = abs[!duplicated(abs),]
    }
    else {
      def[i,3] = round(length(occ.inn)/nrow(pred.inn) * nrow(predictors),0)
      abs =  as.data.frame(pred.out[sample(nrow(pred.out), def[i,3], replace = F, prob = dist.p),1:2])
    }
    occ.A <- rbind(occ.A, cbind(abs, species = sp, PA = 0))
    if (plot){
      points(abs, pch = 19, col = "#FF000080")
      points(occ[,1:2], pch = 18, col = "#0000FF80")
    }
  }
  out$Presences = occ.P
  out$Absences = occ.A
  out$tab = def
  return(out)
}

