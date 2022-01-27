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
#' @importFrom raster rasterize maxValue image
#' @importFrom sp SpatialPolygons Polygons Polygon SpatialPoints plot
#' @importFrom stats cov.wt na.exclude
#' @importFrom car dataEllipse
#' @import gridExtra
#' @import ggplot2
#' @importFrom graphics layout legend par plot.new
#'
#' @export
plot.NINA <- function(x, ...){

  type = class(x)[2]

  if (type %in% c("ENmodel", "BCmodel", "ECmodel")){
    df.env = x$env.scores
    df.sp = x$sp.scores
    pca = x$pca
    mode.region = if(!is.null(x$clus)){TRUE} else {FALSE}
    n.maps = raster::nlayers(x$maps)
    pos.maps = rep(1:n.maps, ceiling(4/n.maps))[1:4]
    ellipse.env <- car::dataEllipse(df.env[,3], df.env[,4], levels = 0.9, xlab = "PC1", ylab = "PC2", draw = F)
    center.env <- stats::cov.wt(na.exclude(df.env[,3:4]))$center
    ellipse.env = sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(list(sp::SpatialPoints(ellipse.env)))), 1)))
    ellipse.occ <- car::dataEllipse(df.sp[,"Axis1"], df.sp[,"Axis2"], levels = 0.9, xlab = "PC1", ylab = "PC2", draw = F)
    center.occ <- stats::cov.wt(na.exclude(df.sp[,c("Axis1", "Axis2")]))$center
    ellipse.occ = sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(list(sp::SpatialPoints(ellipse.occ)))), 1)))
    if (mode.region){
      layout(matrix(c(5,5,6,6,
                      5,5,6,6,
                      1,2,3,4), 3, 4, byrow = T))
      clus = raster::rasterize(x$clus[,1:2], x$maps[[1]], field = as.numeric(as.factor(x$clus[,3])), fun = "last", na.rm = T)
      clusNames = levels(as.factor(x$clus[,3]))
      n.clus = length(clusNames)
      tab = x$tab
      df.env <- merge(df.env, x$clus, by = c("x", "y"))
      reg.cols <- colorspace::lighten(viridis::viridis(n.clus), 0.4)
      pca.cols <- rep(reg.cols)[as.factor(df.env[,5])]
    } else{
      layout(matrix(c(5,5,5,5,
                      5,5,5,5,
                      1,2,3,4), 3, 4, byrow = T))
      pca.cols = "#00EE0080"
    }
    if(n.maps <= 4) {
      for (i in 1:n.maps) {
        raster::image(x$maps[[i]], xlab  = names(x$maps[[i]]), ylab = "")
        par(mar=c(5,4,4,4))
      }
    }
    else {
      for (i in 1:4) {
        raster::image(x$maps[[i]], xlab  = names(x$maps[[i]]), ylab = "" )
        par(mar=c(5,4,4,4))
      }
      message(paste("Ploting only the first four species maps of a total of", n.maps))
    }
    par(mar=c(6,4,3,4))
    ## env
    plot(df.env[,3:4], col = pca.cols, pch = 20,
         xlab = paste0("Axis1 (", round(pca$eig[1]/sum(pca$eig) * 100, 2), "%)"),
         ylab = paste0("Axis2 (", round(pca$eig[2]/sum(pca$eig) * 100, 2), "%)"))
    sp::plot(ellipse.env, add = T,  border = "#089908", lty = 2, lwd = 2)
    ## occs
    points(df.sp[,c("Axis1", "Axis2")], col = "grey20", pch = 4)
    sp::plot(ellipse.occ, add = T,  border = "grey20", lty = 2, lwd = 2)
    ## env contributions
    co.mp <- floor(min(apply(apply(pca$li, 2, range), 2, diff))/max(apply(apply(pca$co, 2, range), 2, diff)))
    ade4::s.corcircle(pca$co*co.mp, clabel = 1.5,
                      box = F,  grid = F, fullcircle = F, add.plot = T)
    legend(center.env[1], center.env[2], "environment", box.lty = 0 , bg = "#089908", text.col = "white",
           xjust = 0.5, yjust = 0.5, x.intersp = -0.5, y.intersp = 0.1)
    legend(center.occ[1], center.occ[2], "occurrences", box.lty = 0 , bg = "grey20", text.col = "white",
           xjust = 0.5, yjust = 0.5, x.intersp = -0.5, y.intersp = 0.1)
    if (mode.region){
      #par(mar=c(1,1,1,1))
      raster::image(clus, col = reg.cols, ylab = "lat", xlab = "lon")
      points(df.sp[,c("x", "y")], col = "grey20", pch = 4)
      raster::plot(clus, legend.only = T, breaks = seq(0.5,raster::maxValue(clus)+0.5,1), col = viridis::viridis(n.clus),
           axis.args=list(at = 1:raster::maxValue(clus),labels=clusNames))
      #par(mar=c(4,4,4,4))
      #plot.new()
      #plotrix::addtable2plot(0,0,tab,bty="o",display.rownames=T,hlines=F, cex=1.5)
    }
  }
  else if (type == "eval"){
      grobList = list()
      for (i in unique(x$confusion$species)){
        res.n = x$confusion[x$confusion$test == "random" & x$confusion$species == i,]
        res. = x$confusion[x$confusion$test == "predicted" & x$confusion$species == i,]
        TNR = x$tab[i, "TNR"]
        TPR = x$tab[i, "TPR"]
        AUC = x$tab[i, "AUC"]
        main = i
        p.value = x$threshold[i,"p.value"]
        grobList[[i]] <- ggplot2::ggplot(res.n, ggplot2::aes(TNR, TPR)) +
          ggplot2::scale_y_continuous("sensitivity", limits = c(0,1), breaks = seq(0,1,0.25), expand = c(0.01,0.01)) +
          ggplot2::scale_x_reverse("specificity", limits = c(1,0), breaks = seq(0,1,0.25), expand = c(0.01,0.01)) +
          ggplot2::geom_tile(fill = "#132B42") +
          ggplot2::stat_density_2d(aes_string(fill = "..level..", col = "..level.."), geom = "polygon",
                                alpha = 0.1, bins = 10) +
                #geom_density_2d(col = "#E69F00" ) +
                #stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE, n = 100) +
          ggplot2::scale_color_gradient(low = "#132B42", high =  "#96B3C9") +
          ggplot2::geom_path(data = res., ggplot2::aes(TNR, TPR), col = "#E69F00", size = 1) +
          ggplot2::geom_point(data = data.frame(x = TNR, y = TPR), aes_string("x", "y"), col = "#E69F00", shape = 18, size = 4) +
          ggplot2::annotate("text", x = TNR, y = TPR, label =  paste("italic(p)==", format.pval(p.value)), parse = T, size = 2, hjust = -0.15, vjust = 1)+
          ggplot2::theme_classic() +
          ggplot2::coord_fixed() +
          ggplot2::labs(title= gsub("\\s\\s", "\n", paste0(gsub('\\.', ' ', main), "  AUC=", round(AUC,2))), parse = T) +
          ggplot2::theme(legend.position='none',
                      plot.margin = unit(c(1,1,1,1), "lines"),
                      panel.border=ggplot2::element_blank(),
                      panel.grid.major=ggplot2::element_blank(),
                      panel.grid.minor=ggplot2::element_blank(),
                      plot.title = ggplot2::element_text(color = "black", size = 10, vjust = 0.5, hjust = 0.5),
                      axis.title.x = ggplot2::element_text(color = "black", size = 10, vjust = 0.5, hjust = 0.5),
                      axis.title.y = ggplot2::element_text(color = "black", size = 10, vjust = 1, hjust = 0.5),
                      axis.text.x = ggplot2::element_text(color = "black", size = 8),
                      axis.text.y = ggplot2::element_text(color = "black", size = 8 ),
                      axis.line = ggplot2::element_blank())
      }
      tab <- tidyr::gather(x$tab, "test", "value", -c(2,12) )
      tab$test <- factor(tab$test,levels = c("Pearson's correlation",  "Jaccard Similarity",
                                         "TPR" ,  "TNR", "TSS","ACC", "AUC", "kappa", "PPV", "NPV"))
      bplot <- ggplot2::ggplot(x$tab, aes_string(x = 1, y = "AUC")) +
        ggplot2::geom_violin(fill = "#E69F00") +
        ggplot2::geom_boxplot(fill = "#E69F00", width = 0.05) +
        ggplot2::scale_y_continuous("Score", breaks = seq(0,1,0.25), expand = c(0,0), limits = c(-0.05,1.05)) +
        ggplot2::scale_x_discrete("AUC", breaks = c(0.5,1.5), expand = c(0,0)) +
        ggplot2::labs(title= "All models") + theme_classic() +
        ggplot2::coord_cartesian(xlim=c(0,2), ylim=c(-0.05,1.05)) +
        ggplot2::annotate(x=0, xend=0, y=0, yend=1, lwd = 0.75, colour="black", geom="segment") +
        ggplot2::annotate(x=0.5, xend=1.5, y=-0.05, yend=-0.05,  colour="black", lwd=1, geom="segment") +
        ggplot2::theme(panel.border=ggplot2::element_blank(),
              panel.grid.major=ggplot2::element_blank(),
              panel.grid.minor=ggplot2::element_blank(),
              axis.line = ggplot2::element_blank(),
              axis.ticks = ggplot2::element_line(size = 0.75),
              axis.ticks.length=unit(.2, "cm"),
              plot.margin = unit(c(2, 0, 2, 0), "lines"),
              axis.title.x = ggplot2::element_text( size=14, ),
              axis.title.y = ggplot2::element_text( size=14, ),
              axis.text.x =ggplot2::element_text( size=12, ),
              axis.text.y = ggplot2::element_text(size=12,))
      grobList <- gridExtra::arrangeGrob(grobs = grobList, ncol=ceiling(sqrt(length(grobList))))
          gridExtra::grid.arrange(grobList, bplot, ncol=2, widths=c(3,1))
  }
  par(mfrow=c(1,1))

}

