#' @title NICHE COMPARISON
#'
#' @param x NINA niche object
#' @param y NINA niche object
#' @param np.type If \code{centroid.w} is TRUE, type of centroid to estimate. Default is 'unimodal'. See \code{\link[NINA]{niche_position}}
#' @param np.metric If \code{centroid.w} is TRUE, statistic method to estimate the niche centroid. Default is 'median'. See \code{\link[NINA]{niche_position}}
#' @param quantile If \code{centroid.w} is TRUE, quantle threshold to filter niche densities. Default is 0.5. See \code{\link[NINA]{niche_position}}
#' @param centroid.w logical whether to weight niche overlap by distance to niche centroid
#' @param cor Logical whether to use environmentally corrected densites
#' @param rep Number of random samples to carry out the statistic tests
#' @param rand Directions of the comparisons. 1 performs all tests with argument 'x' as base of all comparisons. 2 uses argument 'y'. Default is 1.
#' @param alternative String indicating the alternative hypothesis. Default is "greater".
#' @param rnd.test Logical to indicate if randomization tests are to be estimated
#'
#' @description Computes different metrics to analyze differences between two  niches
#'
#' @return An S3 object class NINA containing different niche anallyses.
#'
#' @importFrom raster extent resample xyFromCell as.matrix raster
#'
#' @export
niche_comparison <- function(x, y, centroid.w = F, rnd.test = F,
                             quantile = 0.75, np.type = c("unimodal", "multimodal"), np.metric = c("median", "mean", "max"),
                             cor = F, rep = 100, rand = 1, alternative = c("greater", "lower")){


  if (all(class(x) != c("NINA" , "niche"))){stop("'x' is not of a niche object of class NINA")}
  if (all(class(y) != c("NINA" , "niche"))){stop("'y' is not of a niche object of class NINA")}
  if (length(x$x) != length(y$x)) {stop("niche objects have different spatial grid size")}
  R = length(x$x)

  out <- list()
  np.metric = np.metric[1]
  np.type = np.type[1]
  if (raster::compareRaster(list(x$Z, y$Z), stopiffalse = F) == F){
    message("\nResampling niche 'y' to match the extent of the niche 'x'..." )
    ras = raster::extent(x$Z)
    rasterEx <- raster::extent(ras)
    ras.template <- raster::raster(nrow=R,ncol=R)
    raster::extent(ras.template) <- rasterEx

    for (ii in c("Z", "z", "z.uncor", "z.cor", "w")) {
      y[[ii]] <- raster::resample(y[[ii]], ras.template, method = "ngb")
      y[[ii]][is.na(y[[ii]])] <- 0
    }
  }
  ####           ####
  ## Niche overlap ##
  ####           ####
  message("\nAnalyzing niche overlap...")
  out$D <- list()
  ## Estimate niche overlap
  message("\t-Computing observed niche overlap...",  appendLF = F)
  D = niche_overlap(x, y, cor = cor,
                    centroid.w = centroid.w, type = np.type,
                    method = np.metric, quantile = quantile)
  ## Niche overlap output
  out$D$obs = as.numeric(D)
  message("...OK")
  ## Compute random niche overlaps
  if (rnd.test){
    message(paste("\t-Computing", rep, "random niche overlaps..."), appendLF = F)
    USE.sim = data.frame()
    breadth.sim <- NULL
    Drnd = NULL
    for (i in 1:rep){
      if (rand == 1){
        zrnd = ecospat::ecospat.grid.clim.dyn(x$glob, x$glob1, x$glob1[sample(1:nrow(x$glob1), nrow(x$sp), replace = F),], R = R)
        if (raster::compareRaster(list(x$Z, zrnd$Z), stopiffalse = F) == F){
          for (ii in c("Z", "z", "z.uncor", "z.cor", "w")) {
          zrnd[[ii]] <- raster::resample(zrnd[[ii]], x$Z, method = "ngb")
          zrnd[[ii]][is.na(zrnd[[ii]])] <- 0
          }
        }
        class(zrnd) <- c("NINA", "niche")
        Drnd = c(Drnd, niche_overlap(x, zrnd, cor = cor,
                                     centroid.w = centroid.w, type = np.type,
                                     method = np.metric, quantile = quantile))
        z = x$w
        zab = z*zrnd$w

        A = sum(raster::values(z) > 0)
        B = sum(raster::values(zrnd$w) > 0)
        AB = sum(raster::values(zab) > 0)

      }
      if (rand == 2){
        zrnd = ecospat::ecospat.grid.clim.dyn(y$glob, y$glob1, y$glob1[sample(1:nrow(y$glob1), nrow(y$sp), replace = F),], R = R)
        if (raster::compareRaster(list(y$Z, zrnd$Z), stopiffalse = F) == F){
          for (ii in c("Z", "z", "z.uncor", "z.cor", "w")) {
            zrnd[[ii]] <- raster::resample(zrnd[[ii]], y$Z, method = "ngb")
            zrnd[[ii]][is.na(zrnd[[ii]])] <- 0
          }
        }
        class(zrnd) <- c("NINA", "niche")
        Drnd = c(Drnd, niche_overlap(y, zrnd, cor = cor,
                                     centroid.w = centroid.w, type = np.type,
                                     method = np.metric, quantile = quantile))
        z = y$w
        zab = z*zrnd$w

        A = sum(raster::values(zrnd$w) > 0)
        B = sum(raster::values(z) > 0)
        AB = sum(raster::values(zab) > 0)
      }

      U = 1 - AB/A
      E = 1 - AB/B
      S = AB / (A+B-AB)

      USE.sim <- rbind(USE.sim, cbind(U, S, E))

      breadth.sim = c(breadth.sim, sum(raster::values(zrnd$z.uncor) > 0)/sum(raster::values(zrnd$Z) > 0))
    }
    message("...OK")

    ## Estimate D p-value
    message("\t-Estimating niche overlap p-values for alternative hypotheses...", appendLF = F)

    D.pval = NULL
    if ("greater" %in% alternative){
      D.pval = c(D.pval, greater = (sum(Drnd >= D) + 1)/(length(Drnd) + 1))
    }
    if ("lower" %in% alternative){
      D.pval = c(D.pval, lower = (sum(Drnd <= D) + 1)/(length(Drnd) + 1))
    }
    message("...OK")
    out$D$sim = Drnd
    out$D$pvalue = D.pval
  }
  ####           ####
  ## Niche breadth ##
  ####           ####
  message("Analyzing niche breadths...")
  out$breadth.diff <- list()
  message("\t-Computing observed niche breadth differences...", appendLF = F)
  ## Compute niche position
  x.breadth = sum(raster::values(x$z.uncor) > 0)/sum(raster::values(x$Z) > 0)
  y.breadth = sum(raster::values(y$z.uncor) > 0)/sum(raster::values(y$Z) > 0)
  if (rand == 1){ obs = x.breadth - y.breadth}
  if (rand == 2){ obs = y.breadth - x.breadth}
  out$breadth.diff$obs <- obs
  message("...OK")
  if (rnd.test){
    message("\t-Computing simulated niche breadths...", appendLF = F)
    ## Niche breqdth simulations
    if (rand == 1){ breadth.sim = x.breadth - breadth.sim}
    if (rand == 2){ breadth.sim = y.breadth - breadth.sim}
    out$breadth.diff$sim = breadth.sim
    message("...OK")
    ## Estimate CS p-value
    message("\t-Estimating niche breadth p-values for alternative hypotheses...", appendLF = F)

    breadth.pval = NULL
    if ("greater" %in% alternative){
      breadth.pval = c(breadth.pval, greater = (sum(breadth.sim >= obs) + 1)/(length(breadth.sim) + 1))
    }
    if ("lower" %in% alternative){
      breadth.pval = c(breadth.pval, lower = (sum(breadth.sim <= obs) + 1)/(length(breadth.sim) + 1))
    }
    ## niche breadth output test
    out$breadth.diff$pvalue <-  breadth.pval
    message("...OK")
  }
  ####            ####
  ## Niche position ##
  ####            ####
  message("Analyzing niche centroids...")
  out$np <- list()
  ## Compute niche position
  message("\t-Computing observed niche position...", appendLF = F)
  Cx = as.numeric(niche_position(x,  type = "unimodal", method = np.metric, quantile = quantile, cor = cor))
  Cy = as.numeric(niche_position(y,  type = "unimodal", method = np.metric, quantile = quantile, cor = cor))
  names(Cx) <- names(Cy) <- c("Axis1", "Axis2")
  out$np$obs = cbind(rbind(Cx, Cy), breadth = c(x.breadth, y.breadth))
  message("...OK")
  if (rnd.test){

    ## Compute random niche centroids
    rep.warning = F
    message(paste("\t-Computing", rep, "random niche centroids..."), appendLF = F)
    if (rand == 1){
      if (cor) {rndf = raster::as.data.frame(x$z.cor, xy = T)}
      else{rndf = raster::as.data.frame(x$z.uncor, xy = T)}
      rndf = rndf[rndf[,3] != 0,]
      rndf = rndf[rndf[,3] >= quantile(rndf[,3], quantile, na.rm = T),]
      if (rep > nrow(rndf)){
        rep = nrow(rndf)
        message(paste0("\nWARNING: Number of randomzations exceed the number of environments available.\n\t ...All environments (", rep, ") will be considered instead.\n\t ...To avoid this, either decrease the number of randomizations or lower the quantile threshold"))
        rep.warning = T
      }
      Crnd = rndf[sample(1:nrow(rndf), rep, prob = rndf[,3], replace = T),1:2]
      drnd = unlist(sqrt((Cx[1]-Crnd[1])^2 + (Cx[2]-Crnd[2])^2))
      ## Compute random alpha values within the selected niche
      Arnd = apply(Crnd[,1:2], 1, function(i) Morpho::angle.calc(Cx, i) / pi)
    }
    if (rand == 2){
      if (cor) {rndf = raster::as.data.frame(y$z.cor, xy = T)}
      else{rndf = raster::as.data.frame(y$z.uncor, xy = T)}
      rndf = rndf[rndf[,3] != 0,]
      rndf = rndf[rndf[,3] >= quantile(rndf[,3], quantile, na.rm = T),]
      if (rep > nrow(rndf)){
        rep = nrow(rndf)
        message(paste0("\nWARNING: Number of randomzations exceed the number of environments available.\n\t ...All environments (", rep, ") will be considered instead.\n\t ...To avoid this, either decrease the number of randomizations or lower the quantile threshold"))
        rep.warning = T
      }
      Crnd = rndf[sample(1:nrow(rndf), rep, prob = rndf[,3], replace = T),1:2]
      drnd = unlist(sqrt((Cy[1]-Crnd[1])^2 + (Cy[2]-Crnd[2])^2))
      ## Compute random alpha values within the selected niche
      Arnd = apply(Crnd[,1:2], 1, function(i) Morpho::angle.calc(Cy, i) / pi)
    }
    if(rep.warning){message("\n\t...OK")}
    else {message("...OK")}

    ## Niche position output
    out$np$sim = cbind(Crnd, rand)
    out$np$rep = data.frame(rep = rep, warning = rep.warning, row.names = "")
  }
  ####            ####
  ## Centroid Shift ##
  ####            ####
  message(paste("\t-Computing Centroid Shift..."), appendLF = F)
  out$CS <- list()
  ## Compute observed Centroid Shift
  CS = as.numeric(sqrt((Cx[1]-Cy[1])^2 + (Cx[2]-Cy[2])^2))
  out$CS$obs  <- CS
  message("...OK")
  if (rnd.test){
    ## Estimate CS p-value
    message("\t-Estimating Centroid Shift p-values for alternative hypotheses...", appendLF = F)
    CS.pval = NULL
    if ("greater" %in% alternative){
      CS.pval = c(CS.pval, greater = (sum(drnd >= CS) + 1)/(length(drnd) + 1))
    }
    if ("lower" %in% alternative){
      CS.pval = c(CS.pval, lower = (sum(drnd <= CS) + 1)/(length(drnd) + 1))
    }
    ## Centroid Shift output test
    out$CS$sim  <- drnd
    out$CS$pvalue <-  CS.pval
    message("...OK")
  }
  ## Estimate proportion of change on each Axis
  message(paste("\t-Computing Centroid Shift Axes weights..."), appendLF = F)
  dA1 = abs(Cx[1]-Cy[1])
  dA2 = abs(Cx[2]-Cy[2])

  ##
  Sp = (dA1 + dA2 + CS) / 2
  Area = sqrt(Sp*(Sp-dA1)*(Sp-dA2)*(Sp-CS))
  h = Area*2/CS

  Area1 = dA1*h*sin(45*pi/180) / 2
  Area2 = dA2*h*sin(45*pi/180) / 2
  pAx1 <- Area1 / Area
  pAx2 <- Area2 / Area

  out$CS$weights <- c(w = pAx1, w = pAx2)
  message("...OK")

  ####                      ####
  ## Environmental Relocation ##
  ####                      ####
  message(paste("\t-Estimating observed Environmental Relocation..."), appendLF = F)
  out$ER <- list()

  ## Compute observed angle alpha
  alpha = Morpho::angle.calc(Cx, Cy) / pi
  out$ER$obs  <- alpha
  message("...OK")

  if (rnd.test){
    ## Estimate alpha p-value
    message("\t-Estimating Environmental Relocation p-values for alternative hypotheses...", appendLF = F)
    alpha.pval = NULL
    if ("greater" %in% alternative){
      alpha.pval = c(alpha.pval, greater =(sum(Arnd >= alpha) + 1)/(length(Arnd) + 1))
    }
    if ("lower" %in% alternative){
      alpha.pval = c(alpha.pval, lower =(sum(Arnd <= alpha) + 1)/(length(Arnd) + 1))
    }
    message("...OK")

    ## Environmental relocation  output
    out$ER$sim  <- Arnd
    out$ER$pvalue <-  alpha.pval
  }

  ####         ####
  ## USE metrics ##
  ####         ####
  message(paste("Analyzing USE metrics..."))
  out$USE <- list()
  ## Compute USE values
  message("\t-Computing observed USE metrics...", appendLF = F)

  xb = x$w
  yb = y$w
  xyb = xb*yb

  A = sum(raster::values(xb) > 0)
  B = sum(raster::values(yb) > 0)
  AB = sum(raster::values(xyb) > 0)

  U = 1 - AB/A
  E = 1 - AB/B
  S = AB / (A+B-AB)
  out$USE$obs <- c(Unfilling = U, Stability = S, Expansion = E)
  message("...OK")

  if (rnd.test){
    ## Estimate alpha p-value
    message("\t-Estimating USE p-values for alternative hypotheses...", appendLF = F)
    USE.pval = data.frame(matrix(NA, length(alternative), 3), row.names = alternative)
    colnames(USE.pval) = c("Unfilling", "Stability", "Expansion")
    if ("greater" %in% alternative){
      USE.pval["greater", ] =  sapply(1:3, function(i) (sum(USE.sim[,i] >= c(U,S,E)[i]) + 1)/(nrow(USE.sim) + 1))
    }
    if ("lower" %in% alternative){
      USE.pval["lower", ] =  sapply(1:3, function(i) (sum(USE.sim[,i] <= c(U,S,E)[i]) + 1)/(nrow(USE.sim) + 1))
    }
    message("...OK")

    ## USE output
    out$USE$sim <- USE.sim
    out$USE$pvalue <- USE.pval
  }

  attr(out, "class") <- c("NINA", "metrics")

  return(out)
}


