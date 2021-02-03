#' @title ESTIMATE EC MODEL
#'
#' @param en Niche model ito transform niche space
#'
#' @param W Weighting coefficients
#' @param R Niche space grid
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
EC_model_ <- function(en, W, R){

  if( !is.na(maxValue(en$z.uncor)) >  0){
    if(maxValue(W$z.uncor) ==  0){
      ec = en
      ec$Z = W$Z
      ec$z = W$z
      ec$w = W$w
      ec$z.uncor = W$z.uncor
      ec$z.cor = W$z.uncor
    }
    else{
      b = raster::as.data.frame(en$z.uncor, xy = T)
      w = raster::as.data.frame(W$z.uncor, xy = T)
      betas = raster::as.data.frame(W$betas, xy = T)
      b[,1:2] = b[,1:2]*w[,3]
      b = b[w[,3] != 0,]
      ras = c(min(c(w[w[,3] != 0,1],b[b[,3] != 0, 1])),
              max(c(w[w[,3] != 0,1],b[b[,3] != 0, 1])),
              min(c(w[w[,3] != 0,2],b[b[,3] != 0, 2])),
              max(c(w[w[,3] != 0,2],b[b[,3] != 0, 2])))
      rasterEx <- raster::extent(ras)
      ras.template <- raster::raster(nrow=R,ncol=R)
      raster::extent(ras.template) <- rasterEx
      ec = en
      ec$glob = en$glob * raster::extract(W$z.uncor, en$glob)
      ec$glob1 = ec$glob
      ec$sp = en$sp * raster::extract(W$z.uncor, en$sp)
      ec$x = seq(ras[1], ras[2], length.out = R)
      ec$y = seq(ras[3], ras[4], length.out = R)
      if(W$alpha == 1){
        ec$Z = raster::resample(round(W$z.uncor), ras.template, method = "ngb")
        ec$z.uncor = raster::resample(en$z.uncor*W$z.uncor, ras.template, method = "ngb")
        ec$z.uncor = ec$z.uncor / cellStats(ec$z.uncor, "max")
        ec$betas = raster::resample(round(W$betas), ras.template, method = "ngb")
      }
      else{
        ec$Z = rasterize(w[,1:2]*w[,3], ras.template, field = w[,3], fun = mean)
        ec$Z= raster.gaussian.smooth(ec$Z, n = 3,type = mean)
        ec$z.uncor = rasterize(b[,1:2], ras.template, field = b[,3], fun = mean)
        ec$z.uncor= raster.gaussian.smooth(ec$z.uncor, n = 3,type = mean)
        ec$z.uncor[is.na(ec$z.uncor)] <- 0
        ec$z.uncor = ec$z.uncor / cellStats(ec$z.uncor, "max")
        ec$betas = stack(sapply(names(W$betas), function(i) rasterize(betas[,1:2]*w[,3], ras.template, field = betas[,i], fun = mean)))
      }
      ec$Z[is.na(ec$Z)] <- 0
      ec$z.uncor[is.na(ec$z.uncor)] <- 0
      ec$z = ec$z.uncor * cellStats(en$z, "max")
      ec$w <- ec$z.uncor
      ec$w[ec$w > 0] <- 1
      ec$z.cor <- ec$z/ec$Z
      ec$z.cor[is.na(ec$z.cor)] <- 0
      ec$z.cor <- ec$z.cor/cellStats(ec$z.cor, "max")
    }
  }
  else{
    stop("input model has not positive predictions")
  }
  return(ec)
}

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
        t.mod[[e]] = list()
        mod.Val[[e]] = list()
        for (i in names(z.mod[[e]])){
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
  EC$type = "EC"
  return(EC)
}
