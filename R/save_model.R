#' @title  SAVE NICHE MODEL
#'
#' @param m Niche model of class NINA
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
#' }
#'
#' @importFrom raster writeRaster
#'
#' @export
save_model <- function(m, path = "~", project.name = "NINA1"){

  if (class(m)[1] == "NINA" && class(m)[2] == "modelsList") {
    path = file.path(path, project.name)
    dir.create(path)
    for (l in names(m)) {
      save_model(m[[l]], path = path, project.name = l)
    }
    write.table(class(m), file = paste0(path, "/class.txt"), row.names = F, col.names = F)
  }
  else if (class(m)[1] == "NINA" && class(m)[2] %in% c("ENmodel", "BCmodel", "ECmodel")){
    paths = list()
    ### initiate
    path = file.path(path, project.name)
    dir.create(path)
    write.table(class(m), file = paste0(path, "/class.txt"), row.names = F, col.names = F)
    #### zmod
    zmod = m$z.mod
    npath = file.path(path, "z")
    dir.create(npath)
    if (is.list(zmod)){
      for (r in names(zmod)){
        if (all(class(zmod[[r]]) == c("NINA", "niche"))){
          write_niche(zmod[[r]], path = npath, file = gsub("\\/", " ", r))
        }
        else if (is.list(zmod[[r]])){
          dir.create(path = file.path(npath, sub("\\/", " ", r)))
          for (s in names(zmod[[r]])){
            if (all(class(zmod[[r]][[s]]) == c("NINA", "niche"))){
              write_niche(zmod[[r]][[s]], path = file.path(npath,sub("\\/", " ", r)), file = s)
            }
          }
        }
      }
    }
    ####
    #### zmod global
    if (!is.null(m$z.mod.global)){
      zglob = m$z.mod.global
      npath = file.path(path, "zglobal")
      dir.create(npath)
      if (is.list(zglob)){
        for (r in names(zglob)){
          if (all(class(zglob[[r]]) == c("NINA", "niche"))){
            write_niche(zglob[[r]], path = npath, file = gsub("\\/", " ", r))
          }
          else if (is.list(zglob[[r]])){
            dir.create(path = file.path(npath, sub("\\/", " ", r)))
            for (s in names(zglob[[r]])){
              if (all(class(zglob[[r]][[s]]) == c("NINA", "niche"))){
                write_niche(zglob[[r]][[s]], path = file.path(npath,sub("\\/", " ", r)), file = s)
              }
            }
          }
        }
      }
    }
    ####
    #### w
    if(!is.null(m$w)){
      wmod = m$w
      npath = file.path(path, "w")
      dir.create(npath)
      if (is.list(m$w)){
        for (r in names(m$w)){
          if (all(class(m$w[[r]]) == c("NINA", "niche"))){
            write_niche(m$w[[r]], path = npath, file = gsub("\\/", " ", r))
          }
          else if (is.list(m$w[[r]])){
            dir.create(path = file.path(npath, sub("\\/", " ", r)))
            for (s in names(m$w[[r]])){
              if (all(class(m$w[[r]][[s]]) == c("NINA", "niche"))){
                write_niche(m$w[[r]][[s]], path = file.path(npath,sub("\\/", " ", r)), file = s)
              }
            }
          }
        }
      }
    }
    ####
    #### species distribution
    sd <- m$pred.dis
    npath = file.path(path, "sd")
    dir.create(npath)
    write.table(sd, file = paste0(npath, "/species_distributions.txt"))
    sapply(names(m$maps) , function(s) writeRaster(m$maps[[s]], filename = paste0(npath, "/", gsub("\\.", "_", s)), format = "GTiff", overwrite = T))
    ####
    #### info

    npath = file.path(path, "info")
    dir.create(npath)
    write.table(m$tab, file = paste0(npath, "/tab.txt"))
    if(!is.null(m$fail)){ write.table(m$fail, file = paste0(npath, "/fail.txt")) }
    write.table(m$predictors, file = paste0(npath, "/predictors.txt"), row.names = F, col.names = F)
    write.table(as.character(m$crs), file = paste0(npath, "/crs.txt"))
    ####
    #### data

    npath = file.path(path, "data")
    dir.create(npath)
    if (!is.null(m$clus)){
      write.table(m$clus, file = paste0(npath, "/regions.txt"))
    }
    if(is.data.frame(m$obs)) {
      write.table(m$obs, file = paste0(npath, "/occurrences.txt"))
    } else {
      write.table(do.call(rbind, m$obs), file = paste0(npath, "/occurrences.txt"))
    }
    write.table(m$sp.scores, file = paste0(npath, "/sp_scores.txt"))
    write.table(m$env.scores, file = paste0(npath, "/env_scores.txt"))
    saveRDS(m$pca, file = paste0(npath, "/pca.RDS"))

    ####
    #### eval
    if (!is.null(m$eval)){
      npath = file.path(path, "eval")
      dir.create(npath)
      for (n in names(m$eval)) {write.table(m$eval[[n]],  file = paste0(npath, "/", n, ".txt"))}
    }
    ####
  }
  else {
    stop("Object is not of class NINA")
  }
}
