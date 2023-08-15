#' @title  LOAD NICHE MODEL
#'
#' @param path directory path
#' @param project.name name of the folder to be created
#'
#' @description Saves NINA niche model into the disk
#'
#'
#' @examples
#' \dontrun{
#' EN1 <- EN_model(env_data, occ_data1, cluster = "env", n.clus = 5)
#' save_model(EN1, path = "~/Desktop", project.name = "NINA1")
#' load_model("NINA1", path = "~/Desktop")
#' }
#'
#' @importFrom raster raster
#' @importFrom stringr str_detect
#' @importFrom utils read.table
#'
#' @export
load_model <- function(project.name, path = "./"){

  mpath = file.path(path, project.name)
  m <- list()
  f <- list.files(mpath)
  class(m) <- utils::read.table(file.path(mpath, "class.txt"))[,1]
  if (class(m)[1] == "NINA" && class(m)[2] == "modelsList") {
    for (i in f[f != "class.txt"]) {
      m[[i]] <- load_model(i, path = mpath)
    }
  }
  else if (class(m)[1] == "NINA" && class(m)[2] %in% c("ENmodel", "BCmodel", "ECmodel")){

    #### zmod
    if ("z" %in% f){
      m$z.mod = list()
      npath = file.path(mpath, "z")
      fn <- list.files(npath)
      if (all(stringr::str_detect(fn, ".snm"))){
        for (r in fn){
          rn <- gsub("\\.snm", "", r)
          m$z.mod[[rn]] = read_niche(filepath = file.path(npath,r))
          class(m$z.mod[[rn]]) <- c("NINA", "niche")
        }
      }
      else {
        for (r in fn){
          rpath = file.path(npath, r)
          rn <- list.files(rpath)
          if (all(stringr::str_detect(rn, ".snm"))){
            m$z.mod[[r]] = list()
            for (s in rn){
              sn <- gsub("\\.snm", "", s)
              m$z.mod[[r]][[sn]] = read_niche(filepath = file.path(rpath,s))
              class(m$z.mod[[r]][[sn]]) <- c("NINA", "niche")
            }
          }
          else {
            stop("files in folder ", r, " in a format not recognised")
          }
        }
      }
    }
    ####
    #### zmod global
    if ("zglobal" %in% f){
      m$z.mod.global = list()
      npath = file.path(mpath, "zglobal")
      fn <- list.files(npath)
      if (all(stringr::str_detect(fn, ".snm"))){
        for (r in fn){
          rn <- gsub("\\.snm", "", r)
          m$z.mod.global[[rn]] = read_niche(filepath = file.path(npath,r))
          class(m$z.mod.global[[rn]]) <- c("NINA", "niche")
        }
      }
      else {
        stop("files in folder 'zglobal' in a format not recognised")
      }
    }
    #### w
    if ("w" %in% f){
      m$w = list()
      npath = file.path(mpath, "w")
      fn <- list.files(npath)
      if (all(stringr::str_detect(fn, ".snm"))){
        for (r in fn){
          rn <- gsub("\\.snm", "", r)
          m$w[[rn]] = read_niche(filepath = file.path(npath,r))
          class(m$w[[rn]]) <- c("NINA", "niche")
        }
      }
      else {
        for (r in fn){
          rpath = file.path(npath, r)
          rn <- list.files(rpath)
          if (all(stringr::str_detect(rn, ".snm"))){
            m$w[[r]] = list()
            for (s in rn){
              sn <- gsub("\\.snm", "", s)
              m$w[[r]][[sn]] = read_niche(filepath = file.path(rpath,s))
              class(m$w[[r]][[sn]]) <- c("NINA", "niche")
            }
          }
          else {
            stop("files in folder ", r, " in a format not recognised")
          }
        }
      }
    }
    ####
    #### species distribution
    npath = file.path(mpath, "sd")
    m$pred.dis <- utils::read.table(file.path(npath, "species_distributions.txt"))
    map.files <- list.files(npath, pattern = "\\.tif$")
    m$maps <- stack(sapply(map.files , function(s) raster::raster(file.path(npath, s))))
    names(m$maps) <- gsub("\\.tif", "", names(m$maps))
    ####
    #### info
    npath = file.path(mpath, "info")
    m$tab <- utils::read.table(file.path(npath, "tab.txt"))
    if(file.exists(file.path(npath, "fail.txt"))){ m$fail <- utils::read.table(file.path(npath, "fail.txt")) }
    m$predictors <- as.vector(t(utils::read.table(file.path(npath, "predictors.txt"))))
    m$crs <- as.vector(t(utils::read.table(file.path(npath, "crs.txt"))))
    ####
    #### data
    npath = file.path(mpath, "data")
    if("regions.txt" %in% list.files(npath)){
      m$clus <- utils::read.table(file.path(npath, "regions.txt"))
    }
    m$obs <- utils::read.table(file.path(npath, "occurrences.txt"))
    m$sp.scores <- utils::read.table(file.path(npath, "sp_scores.txt"))
    m$env.scores <- utils::read.table(file.path(npath, "env_scores.txt"))
    m$pca <- readRDS(file.path(npath, "pca.RDS"))

    ####
    #### eval
    if("eval" %in% f){
      npath = file.path(mpath, "eval")
      ef <- list.files(npath)
      m$eval <- list()
      for (n in ef) {
        m$eval[[gsub("\\.txt", "", n)]] <- utils::read.table(file.path(npath, n))
      }
      attr(m$eval, "class") <- c("NINA", "eval")
    }
    if ("w" %in% f){
      attr(m, "class") <- c("NINA", "BCmodel")
    } else {
      attr(m, "class") <- c("NINA", "ENmodel")
    }
  }
  else {
    stop("Model not found.")
  }
  return(m)
}
