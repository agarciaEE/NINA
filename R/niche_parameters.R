#' @title NICHE PARAMETERS
#'
#' @param EN NINA object. Environmental niche model
#'
#' @param EC NINA object. Ecological niche model
#' @param type String indicating the scale of the analysis. Default is "region"
#' @param quantile Numeric value defining quantile threshold to filter niche suitability
#' @param np.metric Statistic metric to compute the niche position. Defaullt is median
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
#' @importFrom raster extent compareRaster maxValue values
#'
#' @export
niche_parameters <- function(EN, EC, type = c("region", "global"), quantile = 0.75, np.metric = c("median", "mean", "max")){

  type = type[1]
  np.metric = np.metric[1]
  niche.p = list()
  if (type == "region"){
    if(!is.null(EN$clus) && !is.null(EC$clus)){
      EN.mod = EN$z.mod
      BC.mod = EC$z.mod
      EC.mod = EC$t.mod
      W.mod = EC$w
      for (e in names(EN.mod)){
        message(paste("Estimating niche parameters of species in region", e, "..."))
        niche.p[[e]] = list()
        for (i in names(EN.mod[[e]])){
          message(paste("Computing niche parameters of", i, "..."))
          en = EN.mod[[e]][[i]]
          R = length(en$x)
          if(!is.na(maxValue(BC.mod[[e]][[i]]$z)) && maxValue(BC.mod[[e]][[i]]$z) > 0){
            bc = BC.mod[[e]][[i]]
            ec = EC.mod[[e]][[i]]
            w = W.mod[[e]][[i]]
            Denbc = ecospat::ecospat.niche.overlap(en, bc, cor = F)$D
            if (compareRaster(list(en$z.uncor, ec$z.uncor), stopiffalse = F) == F){
              message("Resampling ecological niche not match the extent of the environmental niche..." )
              ras = extent(en$z.uncor)
              rasterEx <- raster::extent(ras)
              ras.template <- raster::raster(nrow=R,ncol=R)
              raster::extent(ras.template) <- rasterEx
              ec$z.uncor <- raster::resample(ec$z.uncor, ras.template, method = "ngb")
              ec$z.uncor[is.na(ec$z.uncor)] <- 0
              ec$Z <- raster::resample(ec$Z, ras.template, method = "ngb")
              ec$Z[is.na(ec$Z)] <- 0
              ec = as.list(ec)
            }
            g = ec
            g$z.uncor = g$Z / raster::cellStats(g$Z, "max")
            Cen = niche_position(en,  type = "unimodal", method = np.metric, quantile = quantile, cor = F)
            Cbc = niche_position(bc,  type = "unimodal", method = np.metric, quantile = quantile, cor = F)
            Cw = niche_position(w,  type = "unimodal", method = np.metric, quantile = quantile, cor = F)
            endf = raster::as.data.frame(en$z.uncor, xy = T)
            endf = endf[endf[,3] != 0,]
            endf = endf[endf[,3] >= quantile(endf[,3], quantile, na.rm = T),]
            Crnd = endf[sample(1:nrow(endf), 1000, prob = endf[,3], replace = T),1:2]
            Ernd = apply(Crnd[,1:2], 1, function(i) Morpho::angle.calc(Cen, i) / pi)
            drnd = unlist(sqrt((Cen[1]-Crnd[1])^2 + (Cen[2]-Crnd[2])^2))
            denbc = sqrt((Cen[1]-Cbc[1])^2 + (Cen[2]-Cbc[2])^2)
            denw = sqrt((Cen[1]-Cw[1])^2 + (Cen[2]-Cw[2])^2)
            CSrnd = drnd/denw
            CS = as.numeric(denbc/denw)
            CS.pval.l = (sum(CSrnd <= CS) + 1)/(length(CSrnd) + 1)
            CS.pval.g = (sum(CSrnd >= CS) + 1)/(length(CSrnd) + 1)
            ER = Morpho::angle.calc(Cen, Cbc) / pi
            Denec = ecospat::ecospat.niche.overlap(en, ec, cor = F)$D
            Dbcec = ecospat::ecospat.niche.overlap(bc, ec, cor = F)$D
            Denw = ecospat::ecospat.niche.overlap(en, w, cor = F)$D
            Dbcw = ecospat::ecospat.niche.overlap(bc, w, cor = F)$D
            Decg = ecospat::ecospat.niche.overlap(ec, g, cor = F)$D

            U = 1 - ecospat::ecospat.niche.overlap(en, bc, cor = F)$D
            E = 1 - ecospat::ecospat.niche.overlap(w, bc, cor = F)$D
            S = ecospat::ecospat.niche.overlap(en, w, cor = F)$D
            Up = 1 - sum(raster::values(bc$z.uncor) > 0)/sum(raster::values(en$z.uncor) > 0)
            Ep = 1 - sum(raster::values(bc$z.uncor) > 0)/sum(raster::values(w$z.uncor) > 0)
            Sp = sum(raster::values(bc$z.uncor) > 0)/(sum(raster::values(en$z.uncor) > 0) +  sum(raster::values(w$z.uncor) > 0) - sum(raster::values(bc$z.uncor) > 0))

            niche.p[[e]][[i]] <- c(CS = CS,
                               CS.pval.lower = CS.pval.l,
                               CS.pval.greater = CS.pval.g,
                               ER = ER,
                               ER.pval.l = (sum(Ernd <= ER) + 1)/(length(Ernd) + 1),
                               ER.pval.g = (sum(Ernd >= ER) + 1)/(length(Ernd) + 1),
                               EN.breadth = sum(raster::values(en$z.uncor) > 0)/sum(raster::values(en$Z) > 0),
                               BC.breadth = sum(raster::values(bc$z.uncor) > 0)/sum(raster::values(bc$Z) > 0),
                               EC.breadth = sum(raster::values(ec$z.uncor) > 0)/sum(raster::values(ec$Z) > 0),
                               W.breadth = sum(raster::values(w$z.uncor) > 0)/sum(raster::values(w$Z) > 0),
                               EP.breadth = sum(raster::values(ec$Z) > 0)/sum(raster::values(w$z.uncor) > 0),
                               Denbc = Denbc,
                               Denec = Denec,
                               BE = 1-Dbcec,
                               Denw = Denw,
                               Dbcw = Dbcw,
                               Decg = Decg,
                               Unfilling = U,
                               Up = Up,
                               Expansion = E,
                               Ep = Ep,
                               Stability = S,
                               Sp = Sp)
            message("Done")
          }
          else{
            message(paste(i, "has nothing to compare"))
            niche.p[[e]][[i]] <- c(CS = NA,
                               CS.pval.lower = NA,
                               CS.pval.greater = NA,
                               ER = NA,
                               ER.pval.l = NA,
                               ER.pval.g = NA,
                               EN.breadth = sum(raster::values(en$z.uncor) > 0)/sum(raster::values(en$Z) > 0),
                               BC.breadth = 0,
                               EC.breadth = 0,
                               W.breadth = 0,
                               EP.breadth = 0,
                               Denbc = 0,
                               Denec = 0,
                               BE = 1,
                               Denw = 0,
                               Dbcw = 0,
                               Decg = 0,
                               Unfilling = 1,
                               Up = 1,
                               Expansion = 0,
                               Ep = 0,
                               Stability = 0,
                               Sp = 0)
          }
        }
        niche.p[[e]] <- ldply(niche.p[[e]], .id = "species")
        print(niche.p[[e]])
      }
      niche.p <- ldply(niche.p, .id = "region")
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
      message(paste("Computing niche parameters of", i, "..."))
      en = EN.mod[[i]]
      R = length(en$x)
      if(!is.na(maxValue(BC.mod[[i]]$z)) && maxValue(BC.mod[[i]]$z) > 0){
        bc = BC.mod[[i]]
        ec = EC.mod[[i]]
        w = W.mod[[i]]
        Denbc = ecospat::ecospat.niche.overlap(en, bc, cor = F)$D
        if (compareRaster(list(en$z.uncor, ec$z.uncor), stopiffalse = F) == F){
          message("Resampling ecological niche not match the extent of the environmental niche..." )
          ras = raster::extent(en$z.uncor)
          rasterEx <- raster::extent(ras)
          ras.template <- raster::raster(nrow=R,ncol=R)
          raster::extent(ras.template) <- rasterEx
          ec$z.uncor <- raster::resample(ec$z.uncor, ras.template, method = "ngb")
          ec$z.uncor[is.na(ec$z.uncor)] <- 0
          ec$Z <- raster::resample(ec$Z, ras.template, method = "ngb")
          ec$Z[is.na(ec$Z)] <- 0

        }
        g = ec
        g$z.uncor = g$Z
        Cen = niche_position(en,  type = "unimodal", method = np.metric, quantile = quantile, cor = F)
        Cbc = niche_position(bc,  type = "unimodal", method = np.metric, quantile = quantile, cor = F)
        Cw = niche_position(w,  type = "unimodal", method = np.metric, quantile = quantile, cor = F)
        endf = raster::as.data.frame(en$z.uncor, xy = T)
        endf = endf[endf[,3] != 0,]
        endf = endf[endf[,3] >= quantile(endf[,3], quantile, na.rm = T),]
        Crnd = endf[sample(1:nrow(endf), 1000, prob = endf[,3], replace = T),1:2]
        Ernd = apply(Crnd[,1:2], 1, function(i) Morpho::angle.calc(Cen, i) / pi)
        drnd = unlist(sqrt((Cen[1]-Crnd[1])^2 + (Cen[2]-Crnd[2])^2))
        denbc = sqrt((Cen[1]-Cbc[1])^2 + (Cen[2]-Cbc[2])^2)
        denw = sqrt((Cen[1]-Cw[1])^2 + (Cen[2]-Cw[2])^2)
        CSrnd = drnd/denw
        CS = as.numeric(denbc/denw)
        CS.pval.l = (sum(CSrnd <= CS) + 1)/(length(CSrnd) + 1)
        CS.pval.g = (sum(CSrnd >= CS) + 1)/(length(CSrnd) + 1)
        ER = Morpho::angle.calc(Cen, Cbc) / pi
        Denec = ecospat::ecospat.niche.overlap(en, ec, cor = F)$D
        Dbcec = ecospat::ecospat.niche.overlap(bc, ec, cor = F)$D
        Denw = ecospat::ecospat.niche.overlap(en, w, cor = F)$D
        Dbcw = ecospat::ecospat.niche.overlap(bc, w, cor = F)$D
        Decg = ecospat::ecospat.niche.overlap(ec, g, cor = F)$D
        U = 1 - ecospat::ecospat.niche.overlap(en, bc, cor = F)$D
        E = 1 - ecospat::ecospat.niche.overlap(w, bc, cor = F)$D
        S = ecospat::ecospat.niche.overlap(en, w, cor = F)$D
        Up = 1 - sum(raster::values(bc$z.uncor) > 0)/sum(raster::values(en$z.uncor) > 0)
        Ep = 1 - sum(raster::values(bc$z.uncor) > 0)/sum(raster::values(w$z.uncor) > 0)
        Sp = sum(raster::values(bc$z.uncor) > 0)/(sum(raster::values(en$z.uncor) > 0) +  sum(raster::values(w$z.uncor) > 0) - sum(raster::values(bc$z.uncor) > 0))

        niche.p[[i]] <- c(CS = CS,
                               CS.pval.lower = CS.pval.l,
                               CS.pval.greater = CS.pval.g,
                               ER = ER,
                               ER.pval.l = (sum(Ernd <= ER) + 1)/(length(Ernd) + 1),
                               ER.pval.g = (sum(Ernd >= ER) + 1)/(length(Ernd) + 1),
                               EN.breadth = sum(raster::values(en$z.uncor) > 0)/sum(raster::values(en$Z) > 0),
                               BC.breadth = sum(raster::values(bc$z.uncor) > 0)/sum(raster::values(bc$Z) > 0),
                               EC.breadth = sum(raster::values(ec$z.uncor) > 0)/sum(raster::values(ec$Z) > 0),
                               W.breadth = sum(raster::values(w$z.uncor) > 0)/sum(raster::values(w$Z) > 0),
                               EP.breadth = sum(raster::values(ec$Z) > 0)/sum(raster::values(w$z.uncor) > 0),
                               Denbc = Denbc,
                               Denec = Denec,
                               BE = 1-Dbcec,
                               Denw = Denw,
                               Dbcw = Dbcw,
                               Decg = Decg,
                               Unfilling = U,
                               Up = Up,
                               Expansion = E,
                               Ep = Ep,
                               Stability = S,
                               Sp = Sp)
      }
      else{
        niche.p[[i]] <- c(CS = NA,
                                      CS.pval.lower = NA,
                                      CS.pval.greater = NA,
                                      ER = NA,
                                      ER.pval.l = NA,
                                      ER.pval.g = NA,
                                      EN.breadth = sum(values(en$z.uncor) > 0)/sum(values(en$Z) > 0),
                                      BC.breadth = 0,
                                      EC.breadth = 0,
                                      W.breadth = 0,
                                      EP.breadth = 0,
                                      Denbc = 0,
                                      Denec = 0,
                                      BE = 1,
                                      Denw = 0,
                                      Dbcw = 0,
                                      Decg = 0,
                                      Unfilling = 1,
                                      Up = 1,
                                      Expansion = 0,
                                      Ep = 0,
                                      Stability = 0,
                                      Sp = 0)
      }
    }
    niche.p <- ldply(niche.p, .id = "species")
    niche.p[is.na(niche.p$Denec),"Denec"] = 0
  }
  return(niche.p)
}
