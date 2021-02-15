#' @title  ESTIMATE BC MODEL
#'
#' @param A.matrix m by n matrix indicating the association coefficient (-1 to 1). m are species to be modeled as rows and n interactions as columns
#' @param C.matrix n by n matrix indicating the competition coefficient between interactions (0 to 1).
#' @param x NINA EN model object for species group one
#' @param y NINA EN model object for species group two
#' @param type String indicating whether to perform at a region or a global level. Note that if models have not been estimated at a region level and it is selected it will produce an error
#' @param D Numeric value indicating independence from biotic associations. Value must be comprised between 0 and 1.
#'
#' @description Transforms environmental niche in base to species interactions
#'
#' @return Data frame.
#'
#' @details Returns an error if \code{filename} does not exist.
#'
#' @examples
#' \dontrun{
#' accident_2015 <- fars_read("Project/data/accident_2015.csv.bz2")
#' }
#'
#' @importFrom plyr ldply
#' @importFrom tidyr spread
#' @importFrom raster extract cellStats
#'
#' @export
BC_model <- function(x, y, A.matrix = NULL, C.matrix = NULL, D = 0, type = c("region", "global")){

  type = type[1]
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
          bc <- BC_model_(z, y.mod[[e]], id = i, D = D, A.matrix = A.matrix, C.matrix = C.matrix)
          mod.Val[[e]][[i]] <- cbind(env.scores[rownames(bc$z$glob),1:2], vals = raster::extract(bc$z$z, bc$z$glob))
          z.mod[[e]][[i]] = bc$z
          w.list[[e]][[i]] = bc$w
        }
        z.mod[[e]] <-  z.mod[[e]][unlist(lapply(z.mod[[e]], length) != 0)]
        z.mod[[e]] <- z.mod[[e]][unlist(lapply(z.mod[[e]], function(x) !is.na(maxValue(x$z))))]
        z.mod[[e]] <- z.mod[[e]][unlist(lapply(z.mod[[e]], function(x) maxValue(x$z) != 0))]
        mod.Val[[e]] <- ldply(mod.Val[[e]], data.frame, .id = "species")
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
      bc <- BC_model_(z, y.mod, id = i, D = D, A.matrix = A.matrix, C.matrix = C.matrix)
      mod.Val[[i]] <- cbind(env.scores[,1:2], vals = raster::extract(bc$z$z, bc$z$glob))
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
  #BC$type = "BC"
  message("Models successfully corrected!")
  attr(BC, "class") <- c("NINA", "BCmodel")

  return(BC)
}

