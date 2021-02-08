#' @title Plot
#'
#' @description  Plots the summary of the output models of an object class NINA
#'
#' @param x An object class NINA
#' @param ... Additional arguments for S3 methods
#'
#' @return A \code{NINA} class object description
#'
#' @method plot NINA
#'
#' @examples
#' \dontrun{
#' EN = EN_model(env_data, occ_data1, cluster = "env", n.clus = 5)
#' plot(EN)
#' }
#' \dontrun{
#' EN = EN_model(env_data, occ_data2)
#' plot(EN)
#' }
#'
#' @importFrom ecospat ecospat.plot.contrib
#' @importFrom plotrix addtable2plot
#' @importFrom raster rasterize maxValue
#' @importFrom sp SpatialPolygons Polygons Polygon SpatialPoints
#' @importFrom stats cov.wt
#' @importFrom car dataEllipse
#' @import gridExtra
#' @importFrom graphics layout legend par plot.new
#'
#' @export
plot.NINA <- function(x, ...){

  type = x$type

  if (type %in% c("EN", "BC", "EC")){
    df = merge(x$env.scores, x$obs, by = c(1,2), all = T)
    pca = x$pca
    mode.region = if(!is.null(x$clus)){TRUE} else {FALSE}
    n.maps = raster::nlayers(x$maps)
    pos.maps = rep(1:n.maps, ceiling(4/n.maps))[1:4]
    ellipse.env <- car::dataEllipse(df[,3], df[,4], levels = 0.9, xlab = "PC1", ylab = "PC2", draw = F)
    center.env <- stats::cov.wt(df[,3:4])$center
    ellipse.env = sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(list(sp::SpatialPoints(ellipse.env)))), 1)))
    ellipse.occ <- car::dataEllipse(df[!is.na(df[,5]),3], df[!is.na(df[,5]),4], levels = 0.9, xlab = "PC1", ylab = "PC2", draw = F)
    center.occ <- stats::cov.wt(df[!is.na(df[,5]),3:4])$center
    ellipse.occ = sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(list(sp::SpatialPoints(ellipse.occ)))), 1)))

    if (mode.region){
      layout(matrix(c(1,1,2,2,
                      1,1,3,4,
                      pos.maps+4), 3, 4, byrow = T))
      clus = raster::rasterize(x$clus[,1:2], x$maps[[1]], field = as.numeric(x$clus[,3]), fun = "last", na.rm = T)
      clusNames =levels(x$clus[,3])
      n.clus = length(clusNames)
      tab = x$tab
    } else{
      layout(matrix(c(1,1,2,2,
                      1,1,2,2,
                      pos.maps+2), 3, 4, byrow = T))
    }
    par(mar=c(6,4,3,4))
    plot(df[,3:4], col = "green2", pch = 19)
    plot(ellipse.env, add = T,  border = "green4", lty = 2, lwd = 2)
    legend(center.env[1], center.env[2], "environment",
           xjust = 0.5,      # 0.5 means center adjusted
           yjust = 0.5,      # 0.5 means center adjusted
           x.intersp = -0.5, # adjust character interspacing as you like to effect box width
           y.intersp = 0.1)
    points(df[,3:4], col = rep("blue")[df[,5]], pch = 4)
    plot(ellipse.occ, add = T,  border = "blue4", lty = 2, lwd = 2)
    legend(center.occ[1], center.occ[2], "ocurrences",
           xjust = 0.5,      # 0.5 means center adjusted
           yjust = 0.5,      # 0.5 means center adjusted
           x.intersp = -0.5, # adjust character interspacing as you like to effect box width
           y.intersp = 0.1)
    par(mar=c(6,4,3,4))
    ecospat::ecospat.plot.contrib(pca$co, pca$eig)
    if (mode.region){
      par(mar=c(1,1,1,1))
      plot(clus, col = viridis::viridis(n.clus), legend = F)
      plot(clus, legend.only = T, breaks = seq(0.5,raster::maxValue(clus)+0.5,1), col = viridis::viridis(n.clus),
           axis.args=list(at = 1:raster::maxValue(clus),labels=clusNames))
      par(mar=c(4,4,4,4))
      plot.new()
      plotrix::addtable2plot(0,0,tab,bty="o",display.rownames=T,hlines=F, cex=1.5)
    }
    if(n.maps <= 4) {
      for (i in 1:n.maps) {
        par(mar=c(4,4,6,4))
        plot(x$maps[[i]], main  = names(x$maps[[i]]))
      }
    }
    if(n.maps > 4) {
      for (i in 1:4) {
        if (i == 1){
          par(mar=c(10,6,4,2))
        } else{ par(mar=c(5,3,2,2)) }
        plot(x$maps[[i]], sub  = names(x$maps[[i]]) )
      }
      message(paste("Ploting only the first four species maps of a total of", n.maps))
    }
    par(mfrow=c(1,1))
  }

  if (type == "eval"){
      grobList = list()
      for (i in unique(x$confusion$species)){
        res.n = x$confusion[x$confusion$test == "random" & x$confusion$species == i,]
        res. = x$confusion[x$confusion$test == "predicted" & x$confusion$species == i,]
        TNR = x$tab[i, "TNR"]
        TPR = x$tab[i, "TPR"]
        AUC = x$tab[i, "AUC"]
        main = i
        p.value = x$threshold[i,"p.value"]
        grobList[[i]] <- ggplot(res.n, aes(TNR, TPR)) +
                scale_y_continuous("sensitivity", limits = c(0,1), breaks = seq(0,1,0.25), expand = c(0.01,0.01)) +
                scale_x_reverse("specificity", limits = c(1,0), breaks = seq(0,1,0.25), expand = c(0.01,0.01)) +
                geom_tile(fill = "#132B42") +
                stat_density_2d(aes_string(fill = "..level..", col = "..level.."), geom = "polygon",
                                alpha = 0.1, bins = 10) +
                #geom_density_2d(col = "#E69F00" ) +
                #stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE, n = 100) +
                scale_color_gradient(low = "#132B42", high =  "#96B3C9") +
                geom_path(data = res., aes(TNR, TPR), col = "#E69F00", size = 1) +
                geom_point(data = data.frame(x = TNR, y = TPR), aes_string("x", "y"), col = "#E69F00", shape = 18, size = 4) +
                annotate("text", x = TNR, y = TPR, label =  paste("italic(p)==", format.pval(p.value)), parse = T, size = 2, hjust = -0.15, vjust = 1)+
                theme_classic() +
                coord_fixed() +
                labs(title= gsub("\\s\\s", "\n", paste0(gsub('\\.', ' ', main), "  AUC=", round(AUC,2))), parse = T) +
                theme(legend.position='none',
                      plot.margin = unit(c(1,1,1,1), "lines"),
                      panel.border=element_blank(),
                      panel.grid.major=element_blank(),
                      panel.grid.minor=element_blank(),
                      plot.title = element_text(color = "black", size = 10, vjust = 0.5, hjust = 0.5),
                      axis.title.x = element_text(color = "black", size = 10, vjust = 0.5, hjust = 0.5),
                      axis.title.y = element_text(color = "black", size = 10, vjust = 1, hjust = 0.5),
                      axis.text.x = element_text(color = "black", size = 8),
                      axis.text.y = element_text(color = "black", size = 8 ),
                      axis.line = element_blank())
      }
      tab <- tidyr::gather(x$tab, "test", "value", -c(2,12) )
      tab$test <- factor(tab$test,levels = c("Pearson's correlation",  "Jaccard Similarity",
                                         "TPR" ,  "TNR", "TSS","ACC", "AUC", "kappa", "PPV", "NPV"))
      bplot <- ggplot(x$tab, aes_string(x = 1, y = "AUC")) +
        geom_violin(fill = "#E69F00") +
        geom_boxplot(fill = "#E69F00", width = 0.05) +
        scale_y_continuous("Score", breaks = seq(0,1,0.25), expand = c(0,0), limits = c(-0.05,1.05)) +
        scale_x_discrete("AUC", breaks = c(0.5,1.5), expand = c(0,0)) +
        labs(title= "All models") + theme_classic() +
        coord_cartesian(xlim=c(0,2), ylim=c(-0.05,1.05)) +
        annotate(x=0, xend=0, y=0, yend=1, lwd = 0.75, colour="black", geom="segment") +
        annotate(x=0.5, xend=1.5, y=-0.05, yend=-0.05,  colour="black", lwd=1, geom="segment") +
        theme(panel.border=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              axis.line = element_blank(),
              axis.ticks = element_line(size = 0.75),
              axis.ticks.length=unit(.2, "cm"),
              plot.margin = unit(c(2, 0, 2, 0), "lines"),
              axis.title.x = element_text( size=14, ),
              axis.title.y = element_text( size=14, ),
              axis.text.x =element_text( size=12, ),
              axis.text.y = element_text(size=12,))
      grobList <- arrangeGrob(grobs = grobList, ncol=ceiling(sqrt(length(grobList))))
      grid.arrange(grobList, bplot, ncol=2, widths=c(3,1))
  }
}

