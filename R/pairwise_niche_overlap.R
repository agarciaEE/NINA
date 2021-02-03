#' @title SPECIES MODEL PAIRWISE NICHE OVERLAP FUNCTION
#'
#' @param x Niche model. NINA object
#'
#' @param transformed Boolean whetehr use the environmental or ecoligcal niche
#' @param mode Scale of the estimation. Default is region
#' @param write.csv Boolean whether to write a csv file witht the results
#' @param filename Filename
#' @param alternative Type of statistic test. Default is "lower"
#' @param rand.type Type of randomization methos. Default is 1
#' @param cor Boolean whether to use cor or uncor estimation
#' @param rep number of randomization tests
#'
#' @description Estimate all pairwise niche overlap of modelled species
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
#' @importFrom utils combn write.table
#' @importFrom raster compareRaster extent
#'
#' @export
pairwise_niche_overlap <- function(x, transformed = F, mode = c("region", "global"),
                                        write.csv = T, filename = "pairwise.niche.overlap.csv",
                                        alternative = c("lower", "greater", "both"), rand.type = 1,
                                        cor = F, rep = 10){

  mode= mode[1]
  alternative = alternative[1]
  if (mode == "region"){
    if (is.null(x$clus)){
      mode = "global"
    }
    else {
      if (transformed){
        if(!is.null(x$t.mod)){ z.list <- x$t.mod}
        else {
          stop("There is not transformed niche")
        }
      }
      else {
        z.list <- x$z.mod
      }
      tab <- list()
      for (e in names(z.list)){
        z.list[[e]] <-  z.list[[e]][unlist(lapply(z.list[[e]], length) != 0)]
        if (length(z.list[[e]]) >= 2){
          R = length(z.list[[e]][[1]]$x)
          comb <- combn(names(z.list[[e]]),2)
          tab[[e]] <- data.frame(spa = comb[1,], spb = comb[2,])
          l <- sapply(z.list[[e]], function(x) x$Z)
          if (compareRaster(l, stopiffalse = F) == F){
            ras = c(min(sapply(l, function(x) { extent(x)[1] })),
                    max(sapply(l, function(x) { extent(x)[2] })),
                    min(sapply(l, function(x) { extent(x)[3] })),
                    max(sapply(l, function(x) { extent(x)[4] })))
            rasterEx <- raster::extent(ras)
            ras.template <- raster::raster(nrow=R,ncol=R)
            raster::extent(ras.template) <- rasterEx
            for (n in 1:length(z.list[[e]])){
              z.list[[e]][[n]]$x <- seq(ras[1], ras[2], length.out = R)
              z.list[[e]][[n]]$y <- seq(ras[3], ras[4], length.out = R)
              z.list[[e]][[n]]$Z <- raster::resample(z.list[[e]][[n]]$Z, ras.template, method='ngb')
              z.list[[e]][[n]]$Z[is.na(z.list[[e]][[n]]$Z)] = 0
              z.list[[e]][[n]]$z <- raster::resample(z.list[[e]][[n]]$z, ras.template, method='ngb')
              z.list[[e]][[n]]$z[is.na(z.list[[e]][[n]]$z)] = 0
              z.list[[e]][[n]]$z.uncor <- raster::resample(z.list[[e]][[n]]$z.uncor, ras.template, method='ngb')
              z.list[[e]][[n]]$z.uncor[is.na(z.list[[e]][[n]]$z.uncor)] = 0
              z.list[[e]][[n]]$z.cor <- raster::resample(z.list[[e]][[n]]$z.cor, ras.template, method='ngb')
              z.list[[e]][[n]]$z.cor[is.na(z.list[[e]][[n]]$z.cor)] = 0
              z.list[[e]][[n]]$w <- raster::resample(z.list[[e]][[n]]$w, ras.template, method='ngb')
              z.list[[e]][[n]]$w[is.na(z.list[[e]][[n]]$w)] = 0
            }
          }
          for (i in 1:ncol(comb)){
            name.spa = comb[1,i]
            name.spb = comb[2,i]
            tab[[e]][i,"D"] <- round(as.numeric(ecospat.niche.overlap(z.list[[e]][[name.spa]],z.list[[e]][[name.spb]],cor=cor)[1]),3)
            #equ <- ecospat.niche.equivalency.test2(z.list[[e]][[name.spa]],z.list[[e]][[name.spb]],rep=rep)
            #sim <- ecospat.niche.similarity.test2(z.list[[e]][[name.spa]],z.list[[e]][[name.spb]],rep=rep,rand.type = rand.type, alternative = "lower")
            #if (alternative %in% c("lower", "both")){
            #  tab[[e]][i,"equivalence.lower"] <- equ$p.D
            #  tab[[e]][i,"similarity.lower"] <- sim$p.D
            #}
            #if (alternative %in% c("greater", "both")){
            #  tab[[e]][i,"equivalence.greater"] <-(sum(equ$sim$D >= equ$obs$D) + 1)/(length(equ$sim$D) + 1)
            #  tab[[e]][i,"similarity.greater"] <- (sum(sim$sim$D >= sim$obs$D) + 1)/(length(sim$sim$D) + 1)
            #}
          }
          message("Species pairwise niche overlaps in ", e,":\n")
          print(tab[[e]])
          if (write.csv){
            write.table(as.data.frame(cbind(tab[[e]], region = e)),
                        col.names=F, row.names = F, sep = ",",
                        file = filename, append = T)
          }
        }
        else {
          next()
        }
      }
      tab = ldply(tab, data.frame, .id = "region")
    }
  }
  if (mode == "global"){
    if (transformed){
      if(!is.null(x$t.mod.global)){ z.list <- x$t.mod.global}
      else {
        stop("There is not transformed niche")
      }
    }
    else {
      z.list <- x$z.mod.global
    }
    R = length(z.list[[1]]$x)
    comb <- combn(names(z.list),2)
    tab <- data.frame(spa = comb[1,], spb = comb[2,])
    l <- sapply(z.list, function(x) x$Z)
    if (compareRaster(l, stopiffalse = F) == F){
      ras = c(min(sapply(l, function(x) { extent(x)[1] })),
              max(sapply(l, function(x) { extent(x)[2] })),
              min(sapply(l, function(x) { extent(x)[3] })),
              max(sapply(l, function(x) { extent(x)[4] })))
      rasterEx <- raster::extent(ras)
      ras.template <- raster::raster(nrow=R,ncol=R)
      raster::extent(ras.template) <- rasterEx
      for (n in 1:length(z.list)){
        z.list[[n]]$x <- seq(ras[1], ras[2], length.out = R)
        z.list[[n]]$y <- seq(ras[3], ras[4], length.out = R)
        z.list[[n]]$Z <- raster::resample(z.list[[n]]$Z, ras.template, method='ngb')
        z.list[[n]]$Z[is.na(z.list[[n]]$Z)] = 0
        z.list[[n]]$z <- raster::resample(z.list[[n]]$z, ras.template, method='ngb')
        z.list[[n]]$z[is.na(z.list[[n]]$z)] = 0
        z.list[[n]]$z.uncor <- raster::resample(z.list[[n]]$z.uncor, ras.template, method='ngb')
        z.list[[n]]$z.uncor[is.na(z.list[[n]]$z.uncor)] = 0
        z.list[[n]]$z.cor <- raster::resample(z.list[[n]]$z.cor, ras.template, method='ngb')
        z.list[[n]]$z.cor[is.na(z.list[[n]]$z.cor)] = 0
        z.list[[n]]$w <- raster::resample(z.list[[n]]$w, ras.template, method='ngb')
        z.list[[n]]$w[is.na(z.list[[n]]$w)] = 0
      }
    }
    for (i in 1:ncol(comb)){
      name.spa = comb[1,i]
      name.spb = comb[2,i]
      tab[i,"D"] <- round(as.numeric(ecospat.niche.overlap(z.list[[name.spa]],z.list[[name.spb]],cor=cor)[1]),3)
      #equ <- ecospat.niche.equivalency.test2(z.list[[name.spa]],z.list[[name.spb]],rep=rep, alternative = "lower")
      #sim <- ecospat.niche.similarity.test2(z.list[[name.spa]],z.list[[name.spb]],rep=rep,rand.type = rand.type, alternative = "lower")
      #if (alternative %in% c("lower", "both")){
      #  tab[i,"equivalence.lower"] <- equ$p.D
      #  tab[i,"similarity.lower"] <- sim$p.D
      #}
      #if (alternative %in% c("greater", "both")){
      #  tab[i,"equivalence.greater"] <-(sum(equ$sim$D >= equ$obs$D) + 1)/(length(equ$sim$D) + 1)
      #  tab[i,"similarity.greater"] <- (sum(sim$sim$D >= sim$obs$D) + 1)/(length(sim$sim$D) + 1)
      #}
    }
    cat(paste0("Species pairwise niche overlaps:\n"))
    if (write.csv){
      write.csv(tab, file = filename, append = TRUE)
    }
    print(tab)
  }
  return(tab)
}
