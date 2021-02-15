#' @title ESTIMATE EC MODEL
#'
#' @param x NINA EN or BC model
#' @param y Optional. NINA EN Model
#' @param type String indicating whether to perform at a region or a global level. Note that if models have not been estimated at a region level and it is selected it will produce an error
#' @param D Numeric value for independence of interactions
#' @param A.matrix m by n matrix indicating the association coefficient (-1 to 1). m are species to be modeled as rows and n interactions as columns
#' @param C.matrix n by n matrix indicating the competition coefficient between interactions (0 to 1).
#'
#' @description Transform environmental niche space into ecological niche space
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
#' @importFrom raster maxValue rasterize stack
#' @importFrom spatialEco raster.gaussian.smooth
#'
#' @export
EC_model <- function(x, y,  D = 0,  A.matrix = NULL, C.matrix = NULL, type = c("region", "global")){

  w = F
  clus = F
  type = type[1]
  if(!is.null(x$clus)){clus = T}
  if (all(class(x) == c("NINA", "BCmodel"))) {
    EC = x[names(x) != "w"]
    z.mod = x$z.mod
    w.mod = x$w
    w = T
  }
  else if (all(class(x) == c("NINA", "ENmodel"))) {
    EC = x
    z.mod = x$z.mod
    if(missing(y)){ stop("Environmental-only model requires another Environmental-only or Environmental-constrained model as argument 'y'") }
    if (class(y) != "NINA") { stop("Argument 'y' is not a NINA model") }
    if (is.null(A.matrix)) { stop("Argument 'A.matrix' needed")}
    else {
      if(any(!names(x$maps) %in% rownames(A.matrix))){ stop("Some species in models 'x' are not present in argument 'A.matrix")}
      if(any(!colnames(A.matrix) %in% names(y$maps))){ stop("Some species in models 'y' are not present in argument 'A.matrix")}
    }
    if(!is.null(y$clus) != clus) {stop("Model x and model y have been estimated in different scales") }
    y.mod = y$z.mod
  }
  else if (all(class(x) == c("NINA", "ECmodel"))) {
    EC = x
    z.mod = x$t.mod
    if(missing(y)){ stop("Ecological model requires another Ecological model as argument 'y' to be reassessed") }
    if (class(y) != c("NINA", "ECmodel")) { stop("Argument 'y' is not a Ecological model of class NINA") }
    if (is.null(A.matrix)) { stop("Argument 'A.matrix' needed")} else {
      if(any(!names(x$maps) %in% rownames(A.matrix))){ stop("Some species are not present in argument 'A.matrix")}}
    if(sum(duplicated(rbind(x$env.scores, y$env.scores))) == nrow(x$env.scores)) {stop("Model x and model y have been estimated with different environmental spaces") }
    if(!is.null(y$clus) != clus) {stop("Model x and model y have been estimated in different scales") }
    y.mod = y$t.mod
  }

  env.scores = x$env.scores
  if (type == "region"){
    if(clus){
      clus.df = x$clus
      mod.Val = list()
      t.mod = list()
      g.mod = list()
      for (e in names(z.mod)){
        message(paste0("Transforming niche space of species in region ", e, "..."))
        t.mod[[e]] = list()
        g.mod[[e]] = list()
        mod.Val[[e]] = list()
        for (i in names(z.mod[[e]])){
          message(paste0("\tEstimating ecological niche of ", i, "..."), appendLF = if (w){F}else{T})
          en = z.mod[[e]][[i]]
          R = length(en$x)
          if (w) { W = w.mod[[e]][[i]] } else{
            message(paste0("\tComputing biotic constrains of ", i, "..."), appendLF = F)
            bc <- BC_model_(en, y.mod[[e]], id = i, D = D, A.matrix = A.matrix, C.matrix = C.matrix)
            W = bc$w
          }
          ec.mod = EC_model_(en, W, R = R, D = D)
          t.mod[[e]][[i]] = ec.mod$ec
          g.mod[[e]][[i]] = ec.mod$g

          mod.Val[[e]][[i]] <- cbind(env.scores[rownames(t.mod[[e]][[i]]$glob),1:2], vals = raster::extract(t.mod[[e]][[i]]$z, t.mod[[e]][[i]]$glob))
        }
        t.mod[[e]] <-  t.mod[[e]][unlist(lapply(t.mod[[e]], length) != 0)]
        t.mod[[e]] <- t.mod[[e]][unlist(lapply(t.mod[[e]], function(x) !is.na(maxValue(x$z))))]
        t.mod[[e]] <- t.mod[[e]][unlist(lapply(t.mod[[e]], function(x) maxValue(x$z) != 0))]
        mod.Val[[e]] <- ldply(mod.Val[[e]], data.frame, .id = "species")
      }
      t.mod <-  t.mod[unlist(lapply(t.mod, length) != 0)]
      mod.Val <- ldply(mod.Val, data.frame, .id = "region")
      mod.Val <- spread(mod.Val, "species", "vals")[,-1]
      mod.Val[is.na(mod.Val)] = 0
      mod.Val[,-c(1:2)] = apply(mod.Val[,-c(1:2)], 2, function(i) i/max(i, na.rm = T))
    }
    else {
      stop("Regional models not found")
    }
  }
  if (type == "global"){
    if(clus){
      if(!is.null(x$z.mod.global)){
        z.mod = x$z.mod.global
        if (x$type == "EN"){
          if(!is.null(y$z.mod.global)){
            y.mod = y$z.mod.global
          } else {stop("Global models of argument 'y' not found")}
        } else { w.mod = x$w.global}
      } else {stop("Global models of argument 'x' not found")}
    }
    t.mod = list()
    g.mod = list()
    mod.Val = list()
    for (i in names(z.mod)){
      message(paste0("Estimating ecological niche of ", i, "..."), appendLF = if (w){F}else{T})
      en = z.mod[[i]]
      R = length(en$x)
      if (w) { W = w.mod[[i]] } else{
        message(paste0("\tComputing biotic constrains of ", i, "..."), appendLF = F)
        bc <- BC_model_(en, y.mod, id = i, D = D, A.matrix = A.matrix, C.matrix = C.matrix)
        W = bc$w
      }
      if( !is.na(raster::maxValue(en$z.uncor)) && raster::maxValue(en$z.uncor) >  0){
        ec.mod = EC_model_(en, W, R = R, D = D)
        t.mod[[i]] = ec.mod$ec
        g.mod[[i]] = ec.mod$g

        mod.Val[[i]] <- cbind(env.scores[,1:2], vals = raster::extract(t.mod[[i]]$z, t.mod[[i]]$glob))
      }
    }
    t.mod <-  t.mod[unlist(lapply(t.mod, length) != 0)]
    mod.Val <- ldply(mod.Val, data.frame, .id = "species")
    mod.Val <- spread(mod.Val, "species", "vals")
    mod.Val[is.na(mod.Val)] = 0
    mod.Val = cbind(mod.Val[,1:2], apply(mod.Val[,-c(1:2)], 2, function(i) i/max(i, na.rm = T)))
  }
  EC$t.mod <- t.mod
  EC$g <- g.mod
  EC$pred.dis = mod.Val
  EC$maps  = raster_projection(mod.Val, ras = EC$maps[[1]])
  #EC$type = "EC"
  message("Models successfully transformed!")
  attr(EC, "class") <- c("NINA", "ECmodel")
  return(EC)
}

