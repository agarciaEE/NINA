#' @title  ESTIMATE BC MODEL
#'
#' @param A.matrix m by n matrix indicating the association coefficient (-1 to 1). m are species to be modeled as rows and n interactions as columns
#' @param C.matrix n by n matrix indicating the competition coefficient between interactions (0 to 1).
#' @param x NINA EN model object for species group one
#' @param y NINA EN model object for species group two
#' @param type String indicating whether to perform at a region or a global level. Note that if models have not been estimated at a region level and it is selected it will produce an error
#' @param D Numeric value indicating independence from biotic associations. Value must be comprised between 0 and 1.
#' @param method Method; abundances or composition
#' @param cor Logical
#' @param R Integer. size of the grid for the niche space estimate. R represents number of columns and rows
#' @param combine.regions Logical. Whether to combine regional niche models into a global one
#' @param relative.niche Logical
#' @param eval Boolean whether to evaluate the model
#' @param K = Carrying capacity of each environmental cell
#' @param th threshold to perform cut off for model evaluation
#' @param ras raster to constrain pseudoabsences sampling in model evalluation
#' @param res spatial resolution
#' @param plot.eval Logical to whether plot the evaluation
#' @param sample.pseudoabsences Boolean to whether sample pseudo-absences
#' @param rep number of randomzation tests
#' @param best.th method to select the best thresholt. Default is "similarity"
#'
#' @description Transforms environmental niche in base to species interactions
#'
#' @return NINA model
#'
#'
#' @examples
#' \dontrun{
#' EN1 <- EN_model(env_data, occ_data1, cluster = "env", n.clus = 5)
#' EN2 <- EN_model(env_data, occ_data2, cluster = "env", n.clus = 5)
#' BC <- BC_model(EN1, EN2, A.matrix = int_matrix)
#' }
#'
#' @importFrom plyr ldply
#' @importFrom tidyr spread
#' @importFrom raster extract cellStats
#'
#' @export
BC_model <- function(x, y, A.matrix = NULL, C.matrix = NULL,
                     D = 1,  method = c("composition", "densities"), eval = TRUE,
                     relative.niche = T, K = NULL, sample.pseudoabsences = TRUE, R = 100,
                     res = NULL, plot.eval = FALSE, rep = 100, th = NULL, ras = NULL,
                     best.th = c("accuracy", "similarity"), combine.regions = F,
                     cor = F, type = c("region", "global")){

  type = type[1]
  method = method[1]
  BC = x
  w = F
  clus = F
  type = type[1]
  if(!is.null(x$clus)){clus = T}
  if (sum(class(x) %in% c("NINA", "ENmodel", "BCmodel")) != 2) {
    stop("Argument 'x' is not a NINA Environmental model")
  }
  if (sum(class(y) %in% c("NINA", "ENmodel", "BCmodel")) != 2) {
    stop("Argument 'y' is not a NINA Environmental model")
  }
  if (is.null(A.matrix)) { stop("Argument 'A.matrix' needed")}
  else {
    if(any(!names(x$maps) %in% rownames(A.matrix))){ stop("Some species in models 'x' are not present in argument 'A.matrix")}
    if(any(!colnames(A.matrix) %in% names(y$maps))){ stop("Some species in models 'y' are not present in argument 'A.matrix")}
  }
  if(!is.null(y$clus) != clus) {stop("Model x and model y have been estimated in different scales") }

  x.mod = x$z.mod
  y.mod = y$z.mod
  env.scores = x$env.scores
  if (type == "region"){
    if(!is.null(x$clus) && !is.null(y$clus)){
      z.mod = list()
      w.list = list()
      mod.Val = list()
      for (e in names(x.mod)){
        message(paste0("Estimating biotic constrains of species in region ", e, "..."))
        z.mod[[e]] = list()
        w.list[[e]] = list()
        mod.Val[[e]] = list()
        for (i in names(x.mod[[e]])){
          message(paste0("\tAdding biotic constrains to ", i, "..."), appendLF = F)
          z = x.mod[[e]][[i]]
          bc <- BC_model_(z, y.mod[[e]], id = i, D = D, K = K, A.matrix = A.matrix, method = method, cor = cor,  C.matrix = C.matrix)
          if (cor){
            mod.Val[[e]][[i]] <- cbind(env.scores[rownames(bc$z$glob),1:2], vals = raster::extract(bc$z$z/bc$z$Z, bc$z$glob))
          } else {
            mod.Val[[e]][[i]] <- cbind(env.scores[rownames(bc$z$glob),1:2], vals = raster::extract(bc$z$z, bc$z$glob))
          }
          z.mod[[e]][[i]] = bc$z
          w.list[[e]][[i]] = bc$w
        }
        z.mod[[e]] <-  z.mod[[e]][unlist(lapply(z.mod[[e]], length) != 0)]
        z.mod[[e]] <- z.mod[[e]][unlist(lapply(z.mod[[e]], function(x) !is.na(maxValue(x$z))))]
        z.mod[[e]] <- z.mod[[e]][unlist(lapply(z.mod[[e]], function(x) maxValue(x$z) != 0))]
        mod.Val[[e]] <- ldply(mod.Val[[e]], data.frame, .id = "species")
      }
      if (relative.niche){
        z.rev <- reverse_list(z.mod)
        for (ii in names(z.rev)){
          max.zdens = max(sapply(z.rev[[ii]], function(i)  raster::cellStats(i$z, "max")))
          max.Zdens = max(sapply(z.rev[[ii]], function(i)  raster::cellStats(i$Z, "max")))
          for (jj in names(z.rev[[ii]])){
            z = z.rev[[ii]][[jj]]
            z$z.uncor = z$z / max.zdens
            z$z.uncor[is.na(z$z.uncor)] <- 0
            z$w <- z$z.uncor
            z$w[z$w > 0] <- 1
            z$z.cor <- z$z/z$Z
            z$z.cor[is.na(z$z.cor)] <- 0
            z$z.cor <- z$z.cor/max.Zdens
            z.rev[[ii]][[jj]] = z
          }
        }
        z.mod <- reverse_list(z.rev)
      }
      z.mod <-  z.mod[unlist(lapply(z.mod, length) != 0)]
      tab = cbind(ldply(sapply(z.mod, function(x) names(x)), data.frame, .id = "region"), P = 1)
      tab =  spread(tab, "region", "P")
      rownames(tab) <- tab[,1]
      tab <- tab[,-1]
      tab[is.na(tab)] = 0
      mod.Val <- ldply(mod.Val, data.frame, .id = "region")
      mod.Val <- spread(mod.Val, "species", "vals")[,-1]
      mod.Val[is.na(mod.Val)] = 0
      mod.Val[,-c(1:2)] = apply(mod.Val[,-c(1:2)], 2, function(i) i/max(i, na.rm = T))
      BC$tab = t(tab)
      BC$pred.dis = mod.Val
      if (!is.null(BC$z.mod.global) && combine.regions){
        message("\t- Assembling regions into global model...")
        BC$z.mod.global = combine_regions(z.mod, env.scores = env.scores[,3:4], R = R)
      }
    }
    else {
      stop("Regional models not found")
    }
  }
  if (type == "global"){
    if(!is.null(x$clus)){x.mod = x$z.mod.global}
    if(!is.null(y$clus)){y.mod = y$z.mod.global}
    z.mod = list()
    w.list = list()
    mod.Val = list()
    for (i in names(x.mod)){
      message(paste0("\tAdding biotic constrains to ", i, "..."), appendLF = F)
      z = x.mod[[i]]
      bc <- BC_model_(z, y.mod, id = i, D = D, K = K, A.matrix = A.matrix, method = method, cor = cor, C.matrix = C.matrix)
      if (cor) {
        mod.Val[[i]] <- cbind(env.scores[,1:2], vals = raster::extract(bc$z$z, bc$z$glob))
      } else {
        mod.Val[[i]] <- cbind(env.scores[,1:2], vals = raster::extract(bc$z$z/bc$z$Z, bc$z$glob))
      }
      z.mod[[i]] = bc$z
      w.list[[i]] = bc$w
    }
    z.mod <-  z.mod[unlist(lapply(z.mod, length) != 0)]
    tab = names(z.mod)
    mod.Val <- ldply(mod.Val, data.frame, .id = "species")
    mod.Val <- spread(mod.Val, "species", "vals")
    mod.Val[is.na(mod.Val)] = 0
    mod.Val = cbind(mod.Val[,1:2], apply(mod.Val[,-c(1:2)], 2, function(i) i/max(i, na.rm = T)))
    BC$tab = tab
    BC$pred.dis = mod.Val
  }
  BC$z.mod = z.mod
  BC$w = w.list
  BC$maps  = raster_projection(mod.Val, ras = x$maps[[1]])
  if (eval == TRUE){
    message("\t- Carrying out models evaluations...")
    BC$eval <- models_evaluation(BC, sample.pseudoabsences = sample.pseudoabsences, res = res, plot = plot.eval,
                                int.matrix = A.matrix , rep = rep, th = th, best.th = best.th)
  }
  #BC$type = "BC"
  message("Models successfully corrected!")
  attr(BC, "class") <- c("NINA", "BCmodel")

  return(BC)
}

