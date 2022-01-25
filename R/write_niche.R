#' @title  SAVE SPECIES NICHE
#'
#' @param x Niche model of class NINA
#' @param path directory path
#' @param file filename
#'
#' @description Writes a species NINA niche into the disk
#'
#'
#' @examples
#' \dontrun{
#' EN1 <- EN_model(env_data, occ_data1, cluster = "env", n.clus = 5)
#' write_niche(EN1$z.mod$A$spB, path = "./", file = "ENtest")
#' }
#'
#'
#' @export
write_niche <- function(x, path = "~", file = "niche"){

  extension = ".snm"
  file = paste0(file, extension)
  write.table("NINA niche", file.path(path, file), col.names = F, row.names = F)

  if (all(class(x) == c("NINA", "niche"))){
    xm <- lapply(x, as.matrix)
    for (l in names(xm)){
      write.table(paste0("/", l), file = file.path(path, file), quote = FALSE, row.names=FALSE, col.names=FALSE, append = T)
      write.table(xm[[l]], file = file.path(path, file), row.names=FALSE, col.names=FALSE, append = T)
    }
  }
  else {
    stop("x is not a niche object class NINA")
  }
}
