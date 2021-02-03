#' @title  ESTIMATE BC MODEL
#'
#' @param z environmental niche model
#'
#' @param y.list list of environmental niche models of species interactors
#' @param id species name to be estimated. id must match the rownames in A.matrix
#' @param D Binary integer indicating species independence on interactions. 0 for fully dependent, 1 for species that can be present without interactors
#' @param A.matrix m by n matrix indicating the association coefficient (-1 to 1). m are species to be modelled as rows and n interactors as columns
#' @param C.matrix n by n matrix indicating the competition coefficient between interactors (0 to 1).
#'
#' @description Transforms environmental niche in base to species interactions
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
#' @importFrom plyr ldply
#' @importFrom tidyr spread
#' @importFrom raster cellStats
#'
#' @export
BC_model_ <- function(z, y.list, id, D = 0, A.matrix = NULL, C.matrix = NULL ){

  out = list()
  if(length(y.list) > 0){
    out$w <- estimate_w(y.list, id = id, A.matrix = A.matrix, C.matrix = C.matrix)
    wc <- out$w$z.uncor
  }
  else{ wc = 0}
  z$z.uncor <- z$z.uncor * D + z$z.uncor * wc
  z$z.uncor[is.na(z$z.uncor)] <- 0
  z$z.uncor <- z$z.uncor / cellStats(z$z.uncor, "max")
  z$z <- z$z.uncor * cellStats(z$z, "max")
  z$w <- z$z.uncor
  z$w[z$w > 0] <- 1
  z$z.cor <- z$z/z$Z
  z$z.cor[is.na(z$z.cor)] <- 0
  z$z.cor <- z$z.cor/cellStats(z$z.cor, "max")
  out$z = z
  return(out)
}



BC_model <- function(x, y, A.matrix = NULL, C.matrix = NULL, type = c("region", "global")){

  type = type[1]
  BC = x
  env.scores = x$env.scores
  if (type == "region"){
    if(!is.null(x$clus) && !is.null(y$clus)){
      x.mod = x$z.mod
      y.mod = y$z.mod
      z.mod = list()
      w.list = list()
      mod.Val = list()
      for (e in names(x.mod)){
        z.mod[[e]] = list()
        w.list[[e]] = list()
        mod.Val[[e]] = list()
        for (i in names(x.mod[[e]])){
          z = x.mod[[e]][[i]]
          bc <- BC_model_(z, y.mod[[e]], id = i, D = 0, A.matrix = A.matrix, C.matrix = C.matrix)
          mod.Val[[e]][[i]] <- cbind(env.scores[rownames(bc$z$glob),1:2], vals = raster::extract(bc$z$z.uncor, bc$z$glob))
          z.mod[[e]][[i]] = bc$z
          w.list[[e]][[i]] = bc$w
        }
        mod.Val[[e]] <- ldply(mod.Val[[e]], data.frame, .id = "species")
      }
      tab = cbind(ldply(sapply(z.mod, function(x) names(x)), data.frame, .id = "region"), P = 1)
      tab =  spread(tab, "region", "P")
      rownames(tab) <- tab[,1]
      tab <- tab[,-1]
      tab[is.na(tab)] = 0
      mod.Val <- ldply(mod.Val, data.frame, .id = "region")
      mod.Val <- spread(mod.Val, "species", "vals")[,-1]
      mod.Val[is.na(mod.Val)] = 0
      mod.Val[,-c(1:2)] = apply(mod.Val[,-c(1:2)], 2, function(i) i/max(i, na.rm = T))
      BC$tab = tab
      BC$pred.dis = mod.Val
    }
    else {
      stop("Regional models not found")
    }
  }
  if (type == "global"){
    if(!is.null(x$clus)){x.mod = x$z.mod.global}else{x.mod = x$z.mod}
    if(!is.null(y$clus)){y.mod = y$z.mod.global}else{y.mod = y$z.mod}
    z.mod = list()
    w.list = list()
    mod.Val = list()
    for (i in names(x.mod)){
      z = x.mod[[i]]
      bc <- BC_model_(z, y.mod, id = i, D = 0, A.matrix = A.matrix, C.matrix = C.matrix)
      mod.Val[[i]] <- cbind(env.scores[,1:2], vals = raster::extract(bc$z$z.uncor, bc$z$glob))
      z.mod[[i]] = bc$z
      w.list[[i]] = bc$w
    }
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
  BC$type = "BC"
  return(BC)
}
