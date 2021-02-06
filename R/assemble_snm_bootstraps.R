#' @title ASSEMBLING BOOTSTRAP MODELS
#' @param z List of niche models
#'
#' @param env.scores Envirronmental PCA scores
#' @param sp.scores Occurences environmental scores
#' @param w Boolean indicating the use of weigthing coefficients
#' @param bootstrap.eval List of model evaluations
#' @param eval Boolean indicating if use the model evaluations
#' @param threshold Numeric thesthold to select models cut-off. Default 0.5
#' @param cluster Boolean if models have been clustered into regions
#' @param method String indicating which parameter to use to weight models performance. Default "ACC"
#'
#' @description Assemble bootstrap niche models
#'
#' @return Niche model
#'
#'@details Returns an error if \code{filename} does not exist.
#'
#' @examples
#' \dontrun{
#' accident_2015 <- fars_read("Project/data/accident_2015.csv.bz2")
#' }
#'
#' @importFrom raster extent maxValue stack compareRaster
#'
#' @export
assemble_snm_bootstraps <- function(z, env.scores, sp.scores, w = FALSE,
                            bootstrap.eval = NULL, eval = F, threshold = 0.5,
                            cluster = F, method = "ACC"){

  z = reverse_list(z)
  z.mod <- list()
  #mod.Val <- list()
  score = 1
  if(cluster){
    for (e in names(z)){
      z[[e]] = reverse_list(z[[e]])
      z.mod[[e]] <- list()
      #mod.Val[[e]] <- list()
      for (sp in names(z[[e]])){
        zz.l <- list()
        Z.l <- list()
        z.mod[[e]][[sp]] <- z[[e]][[sp]][[1]]
        R = length(z[[e]][[sp]][[1]]$x)
        #mod.Val[[e]][[sp]] <- list()
        sc.vec <- NULL
        for (i in 1:length(z[[e]][[sp]])){
          if(eval){score = bootstrap.eval[[i]][sp,method]}
          if(score > threshold){
            sc.vec <- c(sc.vec, score)
            zz.l[[i]] <- z[[e]][[sp]][[i]]$z
            Z.l[[i]] <- z[[e]][[sp]][[i]]$Z
          }
        }
        zz.l = zz.l[which(sapply(zz.l, "maxValue") > 0)]
        Z.l = Z.l[which(sapply(Z.l, "maxValue") > 0)]
        if (length(zz.l) == 0){
          warning(paste("No partition of", sp, "in", e, "fits the assembling threshold standards."), immediate. = T)
          next()
        }
        if (length(zz.l) > 1){
          if (compareRaster(Z.l, stopiffalse = F) == F){
            ras = c(min(sapply(Z.l, function(x) { extent(x)[1] })),
                    max(sapply(Z.l, function(x) { extent(x)[2] })),
                    min(sapply(Z.l, function(x) { extent(x)[3] })),
                    max(sapply(Z.l, function(x) { extent(x)[4] })))
            rasterEx <- raster::extent(ras)
            ras.template <- raster::raster(nrow=R,ncol=R)
            raster::extent(ras.template) <- rasterEx
            zz.l <- lapply(zz.l, function(x) raster::resample(x, ras.template, method='ngb'))
            Z.l <- lapply(Z.l, function(x) raster::resample(x, ras.template, method='ngb'))
          }
          sc.vec <- sapply(sc.vec, function(i) i/sum(sc.vec))
          #sc.vec <- sc.vec/min(sc.vec)
          rasterEx <- raster::extent(zz.l[[1]])
          z.max <- max(maxValue(stack(zz.l)), na.rm = T)
          zz <- stack(sapply(1:length(zz.l), function(i) zz.l[[i]]*sc.vec[i]))
          Z <- stack(Z.l)
          z.mod[[e]][[sp]]$x = seq(rasterEx[1], rasterEx[2], length.out = R)
          z.mod[[e]][[sp]]$y = seq(rasterEx[3], rasterEx[4], length.out = R)
          sc = sp.scores[sp.scores[,"region"] == e & sp.scores[,"species"] == sp,]
          z.mod[[e]][[sp]]$sp =  sc[!duplicated(sc), c("Axis1", "Axis2")]
          zz <- sum(zz, na.rm=TRUE)
          Z <- mean(Z, na.rm=TRUE)
          z.uncor <- zz/cellStats(zz, "max")
          zz <- zz * z.max
          ww <- z.uncor
          ww[ww > 0] <- 1
          z.cor <- zz/Z
          z.cor[is.na(z.cor)] <- 0
          z.cor <- z.cor/cellStats(z.cor, "max")
          z.mod[[e]][[sp]]$z.uncor  <- z.uncor
          z.mod[[e]][[sp]]$z.cor <- z.cor
          z.mod[[e]][[sp]]$w <- ww
          z.mod[[e]][[sp]]$z <- zz
          z.mod[[e]][[sp]]$Z <- Z
          #sp.coords <- cbind(env.scores[,1:2], env.scores[,3:4])[rownames(z.mod[[e]][[sp]]$glob),]
          #mod.Val[[e]][[sp]] <- cbind(sp.coords[,1:2], vals = raster::extract(z.mod[[e]][[sp]]$z.uncor, z.mod[[e]][[sp]]$glob))
        }
        if (length(zz.l) == 1){
          zz <- zz.l[[1]]
          Z <- Z.l[[1]]
          if (compareRaster(zz,Z, stopiffalse = F) == F){
            rasterEx <- raster::extent(Z)
            ras.template <- raster::raster(nrow=R,ncol=R)
            raster::extent(ras.template) <- rasterEx
            zz <- raster::resample(zz, ras.template, method='ngb')
          }
          z.uncor <- zz/cellStats(zz, "max")
          ww <- z.uncor
          ww[ww > 0] <- 1
          z.cor <- zz/Z
          z.cor[is.na(z.cor)] <- 0
          z.cor <- z.cor/cellStats(z.cor, "max")
          z.mod[[e]][[sp]]$z.uncor  <- z.uncor
          z.mod[[e]][[sp]]$z.cor <- z.cor
          z.mod[[e]][[sp]]$w <- ww
          z.mod[[e]][[sp]]$z <- zz
          z.mod[[e]][[sp]]$Z <- Z
          #sp.coords <- cbind(env.scores[,1:2], env.scores[,3:4])[rownames(z.mod[[e]][[sp]]$glob),]
          #mod.Val[[e]][[sp]] <-  cbind(sp.coords[,1:2], vals = raster::extract(z.mod[[e]][[sp]]$z.uncor, z.mod[[e]][[sp]]$glob))
        }
        if(w) {
          z.mod[[e]][[sp]]$betas <- z[[e]][[sp]][[1]]$betas
          z.mod[[e]][[sp]]$alpha <- z[[e]][[sp]][[1]]$alpha
        }
        message("Assembling of ", sp, " in ", e, " completed.")
      }
      #mod.Val[[e]] <- ldply(mod.Val[[e]], data.frame, .id = "species")
      message(paste("Assembling models of region", e, "succesfully completed."))
    }
  }
  else {
    for (sp in names(z)){
      zz.l <- list()
      Z.l <- list()
      z.mod[[sp]] <- z[[sp]][[1]]
      #mod.Val[[sp]] <- list()
      for (i in names(z[[sp]])){
        if(eval){score = bootstrap.eval[[i]][sp,method]}
        zz.l[[i]] <- z[[sp]][[i]]$z * score
        Z.l[[i]] <- z[[sp]][[i]]$Z
      }
      zz.l = zz.l[which(sapply(zz.l, "maxValue") > 0)]
      Z.l = Z.l[which(sapply(Z.l, "maxValue") > 0)]
      if (length(zz.l) > 1){
        if (compareRaster(zz.l, stopiffalse = F) == F){
          ras = c(min(sapply(zz.l, function(x) { extent(x)[1] })),
                  max(sapply(zz.l, function(x) { extent(x)[2] })),
                  min(sapply(zz.l, function(x) { extent(x)[3] })),
                  max(sapply(zz.l, function(x) { extent(x)[4] })))
          rasterEx <- raster::extent(ras)
          ras.template <- raster::raster(nrow=R,ncol=R)
          raster::extent(ras.template) <- rasterEx
          zz.l <- lapply(zz.l, function(x) raster::resample(x, ras.template, method='ngb'))
          Z.l <- lapply(Z.l, function(x) raster::resample(x, ras.template, method='ngb'))
        }
        rasterEx <- raster::extent(zz.l[[1]])
        zz <- stack(zz.l)
        Z <- stack(Z.l)
        z.mod[[sp]]$x = seq(rasterEx[1], rasterEx[2], (rasterEx[2]-rasterEx[1])/R)
        z.mod[[sp]]$y = seq(rasterEx[3], rasterEx[4], (rasterEx[4]-rasterEx[3])/R)
        sc = sp.scores[sp.scores[,"species"] == sp,]
        z.mod[[sp]]$sp =  sc[!duplicated(sc), c("Axis1", "Axis2")]
        zz <- sum(zz, na.rm=TRUE)/sum(stack(lapply(z[[sp]], function(x) raster::resample(x, ras.template, method='ngb'))))
        Z <- mean(Z, na.rm=TRUE)
        z.uncor <- zz/cellStats(zz, "max")
        ww <- z.uncor
        ww[ww > 0] <- 1
        z.cor <- zz/Z
        z.cor[is.na(z.cor)] <- 0
        z.cor <- z.cor/cellStats(z.cor, "max")
        z.mod[[sp]]$z.uncor  <- z.uncor
        z.mod[[sp]]$z.cor <- z.cor
        z.mod[[sp]]$w <- ww
        z.mod[[sp]]$z <- zz
        z.mod[[sp]]$Z <- Z
        #mod.Val[[sp]] <- cbind(env.scores[,1:2], vals = raster::extract(z.mod[[sp]]$z.uncor, z.mod[[sp]]$glob))
      }
      if (length(zz.l) == 1){
        zz <- zz.l[[1]]
        Z <- Z.l[[1]]
        z.uncor <- zz/cellStats(zz, "max")
        ww <- z.uncor
        ww[ww > 0] <- 1
        z.cor <- zz/Z
        z.cor[is.na(z.cor)] <- 0
        z.cor <- z.cor/cellStats(z.cor, "max")
        z.mod[[sp]]$z.uncor  <- z.uncor
        z.mod[[sp]]$z.cor <- z.cor
        z.mod[[sp]]$w <- ww
        z.mod[[sp]]$z <- zz
        z.mod[[sp]]$Z <- Z
        #mod.Val[[sp]] <-  cbind(env.scores[,1:2], vals = raster::extract(z.mod[[sp]]$z.uncor, z.mod[[sp]]$glob))
      }
      if(w) {
        z.mod[[e]][[sp]]$betas <- z[[e]][[sp]][[1]]$betas
        z.mod[[e]][[sp]]$alpha <- z[[e]][[sp]][[1]]$alpha
      }
    }
  }
  return(list(z.mod = z.mod))
}
