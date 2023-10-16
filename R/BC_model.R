#' @title  ESTIMATE BC MODEL
#'
#' @param A.matrix m by n matrix indicating the association coefficient (-1 to 1). m are species to be modeled as rows and n interactions as columns
#' @param C.matrix n by n matrix indicating the competition coefficient between interactions (0 to 1).
#' @param x NINA EN model object for species group one
#' @param y NINA EN model object for species group two
#' @param type String indicating whether to perform at a "regional" or a "global" level. Note that if models have not been estimated at a region level and it is selected it will produce an error
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
#' @importFrom plyr plyr::ldply
#' @importFrom tidyr tidyr::spread
#' @importFrom raster extract cellStats
#'
#' @export
BC_model <- function(x, y, A.matrix = NULL, C.matrix = NULL,
                     D = 1,  method = c("composition", "densities"), eval = TRUE,
                     relative.niche = T, K = NULL, sample.pseudoabsences = TRUE, R = 100,
                     res = NULL, plot.eval = FALSE, rep = 100, th = NULL, ras = NULL,
                     best.th = c("accuracy", "similarity"), combine.regions = F,
                     cor = F, type = c("regional", "global")){

  type = type[1]
  method = method[1]

  if (!inherits(x, c("NINA", "ENmodel"))) {
    stop("Argument 'x' is not a NINA Environmental model")
  }
  if (!inherits(y, c("NINA", "ENmodel"))) {
    stop("Argument 'y' is not a NINA Environmental model")
  }

  BC = x
  w = F
  clus = F

  if(!identical(x$clus, y$clus)) { stop("Model x and model y have been estimated in different scales") }
  if(inherits(x$clus, "data.frame")){ clus = T }

  if (type == "regional" & !clus){ stop("No regions found to carry out regional models.") }

  xnames <- names(x$maps)
  ynames <- names(y$maps)
  env.scores = x$env.scores
  x.mod = x$z.mod
  y.mod = y$z.mod

  if (type == "regional"){

    if (inherits(A.matrix, "list")) {
      A.matrixList <- A.matrix
      if (any(!names(x.mod) %in% names(A.matrixList))) {
        x.mod <- x.mod[names(A.matrixList)]
        if (length(x.mod) == 0){
          stop("Names on A.matrix list does not correspond to regions names.")
        } else {
          warning("Some regions from models are missing on provided A.matrix list. Only provided regions will be estimated.")
        }
      }
      reg_x_int_match <- sapply(names(x.mod),
                               function(r) any(!sapply(x.mod, names, simplify = F)[[r]] %in%
                                                 sapply(A.matrixList, rownames, simplify = F)[[r]]))
      reg_y_int_match <- sapply(names(y.mod),
                                function(r) any(!sapply(y.mod, names, simplify = F)[[r]] %in%
                                                  sapply(A.matrixList, colnames, simplify = F)[[r]]))
      if(any(reg_x_int_match)){
        x.mod <- sapply(names(x.mod), function(r) x.mod[[r]][names(x.mod[[r]]) %in% rownames(A.matrixList[[r]])])
        if (all(sapply(x.mod, length) == 0)){
            stop("Species in model 'x' are not present in argument list 'A.matrix.")
        }
      }
      if(any(reg_y_int_match)){
        y.mod <- sapply(names(y.mod), function(r) y.mod[[r]][names(y.mod[[r]]) %in% colnames(A.matrixList[[r]])])
        if (all(sapply(y.mod, length) == 0)){
          stop("Species in model 'y' are not present in argument list 'A.matrix.")
        }
      }
    }
    else if (inherits(A.matrix, "data.frame")) {
      if(any(!xnames %in% rownames(A.matrix))){ stop("Some species in model 'x' are not present in argument 'A.matrix")}
      if(any(!colnames(A.matrix) %in% ynames)){ stop("Some species in model 'y' are not present in argument 'A.matrix")}
      if (clus){
        A.matrixList <- rep(list(A.matrix), length(x$z.mod))
        names(A.matrixList) <- names(x$z.mod)
      }
      warning("Only one interaction matrix provided as 'data.frame'. Assuming same interaction matrix in every region.")
    }
    else { stop("Argument 'A.matrix' of class 'data.frame' or 'list' has to be provided.") }

    z.mod = list()
    w.list = list()
    mod.Val = list()
    for (e in names(x.mod)){
      message(paste0("Estimating biotic constrains of species in region ", e, "..."))
      z.mod[[e]] = list()
      w.list[[e]] = list()
      mod.Val[[e]] = list()
      reg.env.scores <- env.scores[x$clus[,3] == e,]
      for (i in names(x.mod[[e]])){
        message(paste0("\tAdding biotic constrains to ", i, "..."), appendLF = F)
        z = x.mod[[e]][[i]]
        bc <- NINA:::BC_model_(z, y.mod[[e]], id = i, D = D, K = K, A.matrix = A.matrixList[[e]], method = method, cor = cor,  C.matrix = C.matrix)
        if (cor){
          mod.Val[[e]][[i]] <- cbind(reg.env.scores[,1:2], vals = raster::extract(bc$z$z/bc$z$Z, reg.env.scores[,3:4]))
        } else {
          mod.Val[[e]][[i]] <- cbind(reg.env.scores[,1:2], vals = raster::extract(bc$z$z, reg.env.scores[,3:4]))
        }
        z.mod[[e]][[i]] = bc$z
        w.list[[e]][[i]] = bc$w
      }
      z.mod[[e]] <-  z.mod[[e]][unlist(lapply(z.mod[[e]], length) != 0)]
      z.mod[[e]] <- z.mod[[e]][unlist(lapply(z.mod[[e]], function(x) !is.na(maxValue(x$z))))]
      z.mod[[e]] <- z.mod[[e]][unlist(lapply(z.mod[[e]], function(x) maxValue(x$z) != 0))]
      mod.Val[[e]] <- plyr::ldply(mod.Val[[e]], data.frame, .id = "species")
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
    #z.mod <-  z.mod[unlist(lapply(z.mod, length) != 0)]
    tab = cbind(plyr::ldply(sapply(z.mod, function(x) names(x)), data.frame, .id = "region"), P = 1)
    tab =  tidyr::spread(tab, "region", "P")
    rownames(tab) <- tab[,1]
    tab <- tab[,-1]
    tab[is.na(tab)] = 0
    mod.Val <- plyr::ldply(mod.Val, data.frame, .id = "region")
    mod.Val <- tidyr::spread(mod.Val, "species", "vals")[,-1]
    mod.Val[is.na(mod.Val)] = 0
    mod.Val[,-c(1:2)] = apply(mod.Val[,-c(1:2)], 2, function(i) i/max(i, na.rm = T))
    BC$tab = t(tab)
    BC$pred.dis = mod.Val
    if (!is.null(BC$z.mod.global) && combine.regions){
      message("\t- Assembling regions into global model...")
      BC$z.mod.global = combine_regions(z.mod, env.scores = env.scores[,3:4], R = R)
    }
  }
  if (type == "global"){

    if (inherits(A.matrix, "list")) {
      if (clus & type == "global"){
      stop("A global interaction matrix of class 'data.frame' has to be provided for global models.")
      }
    }
    else if (inherits(A.matrix, "data.frame")) {
      if(any(!xnames %in% rownames(A.matrix))){ stop("Some species in model 'x' are not present in argument 'A.matrix")}
      if(any(!colnames(A.matrix) %in% ynames)){ stop("Some species in model 'y' are not present in argument 'A.matrix")}
    }
    if(clus){
      x.mod = x$z.mod.global
      y.mod = y$z.mod.global
    }
    z.mod = list()
    w.list = list()
    mod.Val = list()
    for (i in names(x.mod)){
      message(paste0("\tAdding biotic constrains to ", i, "..."), appendLF = F)
      z = x.mod[[i]]
      bc <- NINA:::BC_model_(z, y.mod, id = i, D = D, K = K, A.matrix = A.matrix, method = method, cor = cor, C.matrix = C.matrix)
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
    mod.Val <- plyr::ldply(mod.Val, data.frame, .id = "species")
    mod.Val <- tidyr::spread(mod.Val, "species", "vals")
    mod.Val[is.na(mod.Val)] = 0
    mod.Val = cbind(mod.Val[,1:2], apply(mod.Val[,-c(1:2)], 2, function(i) i/max(i, na.rm = T)))
    BC$tab = tab
    BC$pred.dis = mod.Val
  }
  BC$z.mod = z.mod
  BC$w = w.list
  BC$A = A.matrix
  BC$C = C.matrix
  BC$maps  = raster_projection(mod.Val, ras = x$maps[[1]])
  if (eval == TRUE){
    message("\t- Carrying out models evaluations...")
    BC$eval <- models_evaluation(BC, sample.pseudoabsences = sample.pseudoabsences, res = res, plot = plot.eval,
                                int.matrix = A.matrix , rep = rep, th = th, best.th = best.th)
  }
  message("Models successfully corrected!")
  attr(BC, "class") <- c("NINA", "BCmodel")

  return(BC)
}

