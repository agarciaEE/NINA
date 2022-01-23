#' @title ESTIMATE ENVIRONMENTAL ACCESSIBILITY
#'
#' @param y.list list of niche models
#' @param method Method; abundances or composition
#' @param id species id
#' @param A.matrix Association matrix
#' @param C.matrix Competition matrix
#' @param cor Logical
#' @param K = Carrying capacity of each environmental cell
#'
#' @description Estimates the omega
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
#' @importFrom raster stack
#'
#' @export
estimate_w <- function(y.list, id,  A.matrix = NULL, cor  = F, K = NULL,  method = c("composition", "densities"), C.matrix = NULL){

  method = method[1]
  if (is.null(A.matrix)){
    A.matrix =  matrix(1, nrow = 1, ncol = length(y.list))
    rownames(A.matrix) = id ; colnames(A.matrix) = names(y.list)
  }
  w <- NULL
  if (method == "composition"){
    betas <- estimate_betas(y.list, C.matrix = C.matrix, cor = cor, K = K)
    Xvar = colnames(A.matrix)[A.matrix[id,]  != 0]
    Xvar <- Xvar[Xvar %in% names(y.list)]
    if(length(Xvar) > 0){
      w <- y.list[[Xvar[1]]]
      if (cor){
        if(length(Xvar) == 1){
          w$sp <- as.data.frame(cbind(y.list[[Xvar]]$sp, species = Xvar))
          w$z.cor = betas[[Xvar]]*A.matrix[id,Xvar]
          w$z.cor[is.na(w$z.cor)] <- 0
        } else {
          w$sp <- do.call(rbind, lapply(y.list[Xvar], function(i) i$sp))
          w$z = sum(stack(lapply(y.list[Xvar], function(i) i$z)), na.rm = T)
          w$z.cor = sum(stack(lapply(names(y.list), function(i) betas[[i]]*as.numeric(A.matrix[id,i]))), na.rm = T)
          w$z.cor[is.na(w$z.cor)] <- 0
        }
        w$z <- w$z.cor * w$Z
        w$z.uncor <- w$z/raster::cellStats(w$z, "max")
        w$z.uncor[is.na(w$z.uncor)] <- 0
        w$w <- w$z.uncor
        w$w[w$w > 0] <- 1
        w$z.cor <- w$z.cor/raster::cellStats(w$z.cor, "max")
      } else {
        if(length(Xvar) == 1){
          w$sp <- as.data.frame(cbind(y.list[[Xvar]]$sp, species = Xvar))
          w$z.uncor = betas[[Xvar]]*A.matrix[id,Xvar]
        } else {
          w$sp <- do.call(rbind, lapply(y.list[Xvar], function(i) i$sp))
          w$z = sum(stack(lapply(y.list[Xvar], function(i) i$z)), na.rm = T)
          w$z.uncor = sum(stack(lapply(names(y.list), function(i) betas[[i]]*as.numeric(A.matrix[id,i]))), na.rm = T)
        }
        w$z = w$z.uncor * cellStats(w$z, "max")
        # w$Z <- sum(stack(sapply(y.list[Xvar], function(i) i$Z))) / length(Xvar) # REDUNDAT. Z env stays the same given by y.list[[1]]
        w$z.uncor[is.na(w$z.uncor)] <- 0
        w$w <- w$z.uncor
        w$w[w$w > 0] <- 1
        w$z.cor <- w$z/w$Z
        w$z.cor[is.na(w$z.cor)] <- 0
        w$z.cor <- w$z.cor/raster::cellStats(w$z.cor, "max")
      }
      w$betas <- stack(betas)
      w$alpha <- length(betas[Xvar])/length(y.list)
    }
  }
  if (method == "densities"){
    if(cor){
      betas <- sapply(y.list, function(i) i$z.cor) # probabilities of occurring in the environment based on occ density over environment density proportion
    } else {
      betas <- sapply(y.list, function(i) i$z.uncor) # Probability of occurring in environment based on occ density over the maximum estimated density.
    }
    #### POSIIVE INTERACTIONS
    ## Product of (1-P(positive species presence)) equals to Probability of Positive Absences
    Pvar = colnames(A.matrix)[A.matrix[id,]  > 0]
    Pvar <- Pvar[Pvar %in% names(y.list)]
    PPA <- 1
    for (i in Pvar)  {
      betas[[i]][is.na(betas[[i]])] = 0
      PPA <- (1-betas[[i]] * A.matrix[id, i]) * PPA #Pvar represent positive interactions, A.matrix weights the effect of the interaction (TO CHECK)
    }
    # PPA equals to Probability of all Positive species Absence
    PPP <- 1-PPA # PPP equals to Probability of at least one Positive species Presence
    ## Make positive effects output
    if(length(Pvar)>0) {
      ## Compute species independent relative effects
      Pw <- sapply(Pvar, function(i) betas[[i]]* A.matrix[id, i])
      r <- raster::cellStats(PPP, "max") / sum(stack(Pw))
      r[!is.finite(r)] = 0
      Pw <-  stack(Pw) * r
      w <- y.list[[Pvar[1]]]
      if(length(Pvar) == 1){
        w$sp <- as.data.frame(cbind(y.list[[Pvar]]$sp, species = Pvar))
      } else {
        w$sp <- do.call(rbind, lapply(y.list[Pvar], function(i) i$sp))
      }
      w$Z <- sum(stack(sapply(y.list[Pvar], function(i) i$Z)))[[1]]
      if (cor){
        w$z <- PPP * w$Z
        w$z.uncor <- w$z/raster::cellStats(w$z, "max")
        w$z.uncor[is.na(w$z.uncor)] <- 0
        w$w <- w$z.uncor
        w$w[w$w > 0] <- 1
        w$z.cor <- PPP
        w$z.cor[is.na(w$z.cor)] <- 0
        w$z.cor <- w$z.cor/raster::cellStats(w$z.cor, "max")
      } else{
        w$z <- PPP * raster::cellStats(sum(stack(sapply(y.list[Pvar], function(i) i$z))), "max")
        w$z.uncor <- PPP
        w$z.uncor[is.na(w$z.uncor)] <- 0
        w$w <- w$z.uncor
        w$w[w$w > 0] <- 1
        w$z.cor <- w$z/w$Z
        w$z.cor[is.na(w$z.cor)] <- 0
        w$z.cor <- w$z.cor/raster::cellStats(w$z.cor, "max")
      }
      w$Z <- w$Z / length(Pvar)
      w$betas <- Pw # species specific weights on omega to keep track
      w$alpha <- length(Pvar)/length(y.list) # proportion of positive interactions over total possible interactions
    }
    if(FALSE) { # Not yet implemented
      #### NEGATIVE INTERACTIONS
      ## Product of (1 - P(negative species presence)) equals to Probability of negative absence
      Nvar = colnames(A.matrix)[A.matrix[id,]  < 0]
      Nvar <- Nvar[Nvar %in% names(y.list)]
      ## Make negative effects output
      if(length(Nvar)>0) {
        PNA <- 1
        for (i in Nvar)  PNA <- (1-betas[[i]]* A.matrix[id, i]) * PNA #Nvar represent negative interactions. A.matrix weights the effect of the interaction (TO CHECK)
        # PNA equals to Probability of all Negative species Absence
        PNP <- 1-PNA  # PNP equals to Probability of at least one Negative species Presence
        ## Compute species independent relative effects
        Nw <- sapply(Nvar, function(i) betas[[i]]* A.matrix[id, i])
        r <- raster::cellStats(PNP, "max") / sum(stack(Nw))
        Nw <-  stack(Nw) * r
        Nz <- y.list[[Nvar[1]]]
        Nz$sp <- plyr::ldply(sapply(y.list[Nvar], function(i) i$sp), .id = "species")[,c(2:3,1)]
        Nz$Z <- sum(stack(sapply(y.list[Nvar], function(i) i$Z)))
        Nz$z <- PPP * Nz$Z
        Nz$z.uncor <- Nz$z/raster::cellStats(Nz$z, "max")
        Nz$z.uncor[is.na(Nz$z.uncor)] <- 0
        Nz$w <- Nz$z.uncor
        Nz$w[Nz$w > 0] <- 1
        Nz$z.cor <- Nz$z/Nz$Z
        Nz$z.cor[is.na(Nz$z.cor)] <- 0
        Nz$z.cor <- Nz$z.cor/cellStats(Nz$z.cor, "max")
        Nz$PNP <- PNP # Actual omega coefficients to weight species probabilities
        Nz$betas <- Nw # species independent relative effect to keep track
        Nz$alpha <- length(Pvar)/length(y.list) # proportion of negative interactions over total possible interactions
      }
    }
  }
  return(w)
}
