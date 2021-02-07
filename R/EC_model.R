#' @title ESTIMATE EC MODEL
#'
#' @param BC NINA EN or BC model
#' @param EN Optional. NINA EN Model
#' @param type String indicating whether to perform at a region or a global level. Note that if models have not been estimated at a region level and it is selected it will produce an error
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
EC_model <- function(BC, EN, type = c("region", "global")){

  if (!missing(EN)){ EC = EN }
  else{ EC = BC }
  type = type[1]
  env.scores = EC$env.scores
  if (type == "region"){
    if(!is.null(EC$clus)){
      clus.df = EC$clus
      z.mod = EC$z.mod
      w.mod = BC$w
      mod.Val = list()
      t.mod = list()
      for (e in names(z.mod)){
        message(paste0("Transforming niche space of species in region ", e, "..."))
        t.mod[[e]] = list()
        mod.Val[[e]] = list()
        for (i in names(z.mod[[e]])){
          message(paste0("\tEstimating ecological niche of ", i, "..."))
          en = z.mod[[e]][[i]]
          W = w.mod[[e]][[i]]
          R = length(en$x)
          if( !is.na(maxValue(en$z.uncor)) >  0){
            t.mod[[e]][[i]] = EC_model_(en, W, R)
            mod.Val[[e]][[i]] <- cbind(env.scores[rownames(t.mod[[e]][[i]]$glob),1:2], vals = raster::extract(t.mod[[e]][[i]]$z.uncor, t.mod[[e]][[i]]$glob))
          }
        }
        mod.Val[[e]] <- ldply(mod.Val[[e]], data.frame, .id = "species")
      }
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
    if(!is.null(BC$clus)){
      z.mod = EC$z.mod.global
      w.mod = BC$w.global
    }
    else{
      z.mod = EC$z.mod
      w.mod  = BC$w
    }
    t.mod = list()
    mod.Val = list()
    for (i in names(z.mod)){
      message(paste0("\tEstimating ecological niche of ", i, "..."))
      en = z.mod[[i]]
      W = w.mod[[i]]
      R = length(en$x)
      if(!is.na(maxValue(en$z.uncor)) >  0){
        t.mod[[i]] = EC_model_(en, W, R)
        mod.Val[[i]] <- cbind(env.scores[,1:2], vals = raster::extract(t.mod[[i]]$z.uncor, t.mod[[i]]$glob))
      }
    }
    mod.Val <- ldply(mod.Val, data.frame, .id = "species")
    mod.Val <- spread(mod.Val, "species", "vals")
    mod.Val[is.na(mod.Val)] = 0
    mod.Val = cbind(mod.Val[,1:2], apply(mod.Val[,-c(1:2)], 2, function(i) i/max(i, na.rm = T)))
  }
  EC$pred.dis = mod.Val
  EC$w <- w.mod
  EC$t.mod <- t.mod
  BC$maps  = raster_projection(mod.Val, ras = BC$maps[[1]])
  EC$type = "EC"
  message("Models successfully transformed!")
  attr(EC, "class") <- "NINA"

  return(EC)
}

