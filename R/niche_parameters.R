#' @title NICHE PARAMETERS
#'
#' @param EN NINA object. Environmental niche model
#'
#' @param EC NINA object. Ecological niche model
#' @param type String indicating the scale of the analysis. Default is "region"
#'
#' @description Estimates niche parameters
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
#' @importFrom stats quantile
#' @importFrom ecospat ecospat.niche.overlap
#' @importFrom raster extent compareRaster maxValue getValues
#'
#' @export
niche_parameters <- function(EN, EC, type = c("region", "global")){

  type = type[1]
  niche.p = list()
  if (type == "region"){
    if(!is.null(EN$clus) && !is.null(EC$clus)){
      EN.mod = EN$z.mod
      BC.mod = EC$z.mod
      EC.mod = EC$t.mod
      W.mod = EC$w
      for (e in names(EN.mod)){
        niche.p[[e]] = list()
        for (i in names(EN.mod[[e]])){
          en = EN.mod[[e]][[i]]
          R = length(en$x)
          if(maxValue(BC.mod[[e]][[i]]$z.uncor) == 0){
            niche.p[[e]][[i]] <- c(CS = NA,
                                   CS.pval.lower = NA,
                                   CS.pval.greater = NA,
                                   alpha = NA,
                                   Denbc = 0,
                                   Denec = 0,
                                   Dbcec = 0,
                                   Unfilling = 1,
                                   Up = 1,
                                   Expansion = 0,
                                   Ep = 0,
                                   Stability = 0,
                                   Sp = 0)
          }
          else{
            bc = BC.mod[[e]][[i]]
            ec = EC.mod[[e]][[i]]
            w = W.mod[[e]][[i]]
            Denbc = ecospat.niche.overlap(en, bc, cor = F)$D
            if (compareRaster(list(en$z.uncor, ec$z.uncor), stopiffalse = F) == F){
              ras = extent(en$z.uncor)
              rasterEx <- raster::extent(ras)
              ras.template <- raster::raster(nrow=R,ncol=R)
              raster::extent(ras.template) <- rasterEx
              ec$z.uncor <- raster::resample(ec$z.uncor, ras.template, method = "ngb")
              ec$z.uncor[is.na(ec$z.uncor)] <- 0
            }
            Cen = niche_position(en,  type = "unimodal", method = "median", quantile = 0.95, cor = F)
            Cbc = niche_position(bc,  type = "unimodal", method = "median", quantile = 0.95, cor = F)
            Cw = niche_position(w,  type = "unimodal", method = "median", quantile = 0.95, cor = F)
            endf = raster::as.data.frame(en$z.uncor, xy = T)
            endf = endf[endf[,3] != 0,]
            Crnd = endf[endf[,3] >= quantile(endf[,3], 0.75, na.rm = T),]
            #Crnd = bindf[sample(nrow(bindf), 100, prob = bindf[,3], replace = T),1:2]
            drnd = sqrt((Cen[1]-Crnd[1])^2 + (Cen[2]-Crnd[2])^2)
            deb = sqrt((Cen[1]-Cbc[1])^2 + (Cen[2]-Cbc[2])^2)
            dew = sqrt((Cen[1]-Cw[1])^2 + (Cen[2]-Cw[2])^2)
            CS = as.numeric(deb/dew)
            CS.pval.l = (sum(drnd/dew <= deb) + 1)/(nrow(drnd) + 1)
            CS.pval.g = (sum(drnd/dew >= deb) + 1)/(nrow(drnd) + 1)
            alpha = Morpho::angle.calc(Cen, Cbc) * 180 / pi
            Denec = ecospat.niche.overlap(en, ec, cor = F)$D
            Dbcec = ecospat.niche.overlap(bc, ec, cor = F)$D
            U = 1 - ecospat.niche.overlap(en, bc, cor = F)$D
            E = 1 - ecospat.niche.overlap(w, bc, cor = F)$D
            S = ecospat.niche.overlap(en, w, cor = F)$D
            Up = 1 - sum(getValues(bc$z.uncor) > 0)/sum(getValues(en$z.uncor) > 0)
            Ep = 1 - sum(getValues(bc$z.uncor) > 0)/sum(getValues(w$z.uncor) > 0)
            Sp = sum(getValues(bc$z.uncor) > 0)/(sum(getValues(en$z.uncor) > 0) +  sum(getValues(w$z.uncor) > 0) - sum(getValues(bc$z.uncor) > 0))
            niche.p[[e]][[i]] <- c(CS = CS,
                                   CS.pval.lower = CS.pval.l,
                                   CS.pval.greater = CS.pval.g,
                                   alpha = alpha/180,
                                   Denbc = Denbc,
                                   Denec = Denec,
                                   Dbcec = Dbcec,
                                   Unfilling = U,
                                   Up = Up,
                                   Expansion = E,
                                   Ep = Ep,
                                   Stability = S,
                                   Sp = Sp)
          }
        }
        niche.p[[e]] <- ldply(niche.p[[e]], .id = "species")
      }
      niche.p <- ldply(niche.p, .id = "region")
      niche.p[is.na(niche.p$Dbcec),"Dbcec"] = 0
      niche.p[is.na(niche.p$Denec),"Denec"] = 0
    }
    else {
      stop("Regional models not found")
    }
  }
  if (type == "global"){
    if(!is.null(EN$clus)){
      EN.mod = EN$z.mod.global
    }
    else{
      EN.mod = EN$z.mod
    }
    if(!is.null(EC$clus)){
      BC.mod = EC$z.mod.global
      EC.mod = EC$t.mod.global
      W.mod = EC$w.global
    }
    else{
      BC.mod = EC$z.mod
      EC.mod = EC$t.mod
      W.mod = EC$w
    }
    for (i in names(EN.mod)){
      en = EN.mod[[i]]
      R = length(en$x)
      if(maxValue(BC.mod[[i]]$z.uncor) == 0){
        niche.p[[i]] <- c(CS = NA,
                          CS.pval.lower = NA,
                          CS.pval.greater = NA,
                          alpha = NA,
                          Denbc = 0,
                          Denec = 0,
                          Dbcec = 0,
                          Unfilling = 1,
                          Up = 1,
                          Expansion = 0,
                          Ep = 0,
                          Stability = 0,
                          Sp = 0)
      }
      else{
        bc = BC.mod[[i]]
        ec = EC.mod[[i]]
        w = W.mod[[i]]
        Denbc = ecospat.niche.overlap(en, bc, cor = F)$D
        if (compareRaster(list(en$z.uncor, ec$z.uncor), stopiffalse = F) == F){
          ras = extent(en$z.uncor)
          rasterEx <- raster::extent(ras)
          ras.template <- raster::raster(nrow=R,ncol=R)
          raster::extent(ras.template) <- rasterEx
          ec$z.uncor <- raster::resample(ec$z.uncor, ras.template, method = "ngb")
          ec$z.uncor[is.na(ec$z.uncor)] <- 0
        }
        Cen = niche_position(en,  type = "unimodal", method = "median", quantile = 0.95, cor = F)
        Cbc = niche_position(bc,  type = "unimodal", method = "median", quantile = 0.95, cor = F)
        Cw = niche_position(w,  type = "unimodal", method = "median", quantile = 0.95, cor = F)
        endf = raster::as.data.frame(en$z.uncor, xy = T)
        endf = endf[endf[,3] != 0,]
        Crnd = endf[endf[,3] >= quantile(endf[,3], 0.75, na.rm = T),]
        #Crnd = bindf[sample(nrow(bindf), 100, prob = bindf[,3], replace = T),1:2]
        drnd = sqrt((Cen[1]-Crnd[1])^2 + (Cen[2]-Crnd[2])^2)
        deb = sqrt((Cen[1]-Cbc[1])^2 + (Cen[2]-Cbc[2])^2)
        dew = sqrt((Cen[1]-Cw[1])^2 + (Cen[2]-Cw[2])^2)
        CS = as.numeric(deb/dew)
        CS.pval.l = (sum(drnd/dew <= deb) + 1)/(nrow(drnd) + 1)
        CS.pval.g = (sum(drnd/dew >= deb) + 1)/(nrow(drnd) + 1)
        alpha = Morpho::angle.calc(Cen, Cbc) * 180 / pi
        Denec = ecospat.niche.overlap(en, ec, cor = F)$D
        Dbcec = ecospat.niche.overlap(bc, ec, cor = F)$D
        U = 1 - ecospat.niche.overlap(en, bc, cor = F)$D
        E = 1 - ecospat.niche.overlap(w, bc, cor = F)$D
        S = ecospat.niche.overlap(en, w, cor = F)$D
        Up = 1 - sum(getValues(bc$z.uncor) > 0)/sum(getValues(en$z.uncor) > 0)
        Ep = 1 - sum(getValues(bc$z.uncor) > 0)/sum(getValues(w$z.uncor) > 0)
        Sp = sum(getValues(bc$z.uncor) > 0)/(sum(getValues(en$z.uncor) > 0) +  sum(getValues(w$z.uncor) > 0) - sum(getValues(bc$z.uncor) > 0))
        niche.p[[i]] <- c(CS = CS,
                          CS.pval.lower = CS.pval.l,
                          CS.pval.greater = CS.pval.g,
                          alpha = alpha/180,
                          Denbc = Denbc,
                          Denec = Denec,
                          Dbcec = Dbcec,
                          Unfilling = U,
                          Up = Up,
                          Expansion = E,
                          Ep = Ep,
                          Stability = S,
                          Sp = Sp)
      }
    }
    niche.p <- ldply(niche.p, .id = "species")
    niche.p[is.na(niche.p$Dbcec),"Dbcec"] = 0
    niche.p[is.na(niche.p$Denec),"Denec"] = 0
  }
  return(niche.p)
}
