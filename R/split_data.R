#' @title SPLIT OCCURRENCES FOR CROSS-VALIDATION FUNCTION
#'
#' @description Samples pseudo-absences from a environmental extent based on euclidean distances to the observations in a environmental PCA
#'
#' @param Obs a nx3 data.frame class object indicating longitud, latitud and species in the three first columns
#' @param method Method used to split occurrences. Default is "Euclidean"
#' @param split.percentage split percentage between train adn test data. Default is 0.25
#' @param spsNames Species names to split occurrences from
#'
#' @return List
#'
#'@details Returns an error if \code{filename} does not exist.
#'
#' @examples
#' \dontrun{
#' accident_2015 <- fars_read("Project/data/accident_2015.csv.bz2")
#' }
#'
#'
#' @export
split_data <- function(Obs, spsNames = NULL, method = c("Euclidean", "kmeans"), split.percentage = 0.25){

  method= method[1]
  if (ncol(Obs) == 3) { Obs$PA = 1  }
  if (is.null(spsNames)){spsNames = unique(Obs[,3])}
  pres_train <- NULL
  pres_test <- NULL
  for (i in spsNames){
    df <- Obs[Obs[,3] == i,]
    ddf <- df[,c(1:2,4)]
    nr <- nrow(ddf)
    tr.p <- 1/split.percentage - 1
    rr <- round(nr*split.percentage,0)
    if (rr >= 1){
      if(method == "kmeans"){
        ddf.c <- cluster_regions(cbind(ddf[,1:2],ddf), plot = F, n.clus = rr)
        s <- unlist(lapply(split(ddf.c, ddf.c[,3]), function(r)
          if(nrow(r) >= tr.p) {rownames(r)[sample(nrow(r), round(nrow(r)*split.percentage,0))]}))
        ss <- unlist(sapply(split(ddf.c, ddf.c[,3]), function(r) if(nrow(r) < tr.p) {rownames(r)}))
        ss <- sample(ss, round(length(ss)*split.percentage,0))
        s <- c(s, ss)
        pres_train <- rbind(pres_train, df[!rownames(df) %in% s,])
        pres_test <- rbind(pres_test, df[s,])
        nrow(pres_test)/(nrow(pres_train)+nrow(pres_test))
      }
      if(method == "Euclidean"){
        t <- NULL
        s <- NULL
        ddf.d <- sapply(1:nrow(ddf), function(i) sapply(1:nrow(ddf), function(j) sqrt((ddf[i,1] - ddf[j,1])^2 + (ddf[i,2] - ddf[j,2])^2)))
        while(nr > tr.p){
          tt <- rownames(ddf[sample(nr, tr.p),])
          t <- c(t, tt)
          bck <- ddf[tt,]
          ddf <- ddf[!rownames(ddf) %in% tt,]
          bck.d <- sapply(1:nrow(ddf), function(i) sapply(1:nrow(bck), function(j) sqrt((ddf[i,1] - bck[j,1])^2 + (ddf[i,2] - bck[j,2])^2)))
          bck.d = 1-bck.d/max(ddf.d)
          if (tr.p > 1){ bck.d <- apply(bck.d, 2, function(x) max(x))}
          ss <- rownames(ddf[sample(nrow(ddf),1, prob = bck.d),])
          s <- c(s, ss)
          ddf <- ddf[!rownames(ddf) %in% ss,]
          nr = nrow(ddf)
        }
        if (nr <= tr.p) {
          t <- c(t, rownames(ddf))
        }
        pres_train <- rbind(pres_train, df[t,])
        pres_test <- rbind(pres_test, df[s,])
      }
    }
    else {
      warning("There is not enough occurrences of ", i, " to split. All observations will be used instead.")
      pres_train <- rbind(pres_train, df)
      pres_test <- rbind(pres_test, df)
    }
  }
  return(list(occ.train = pres_train,
              occ.test = pres_test))
}
