#' @title EVALUATE FUNCTION
#'
#' @description Evaluate independent model
#'
#' @param th threshold cut off
#' @param plot Boolean whetehr to plot the result
#' @param Fit. Vector of predicted values
#' @param Obs. Vector of observed values
#' @param rep number of randomization tests
#' @param best.th method to select the best threshold. Default is "similarity"
#' @param main If plot = TRUE. Indicate plot title
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
#' @importFrom stats cor.test rbinom
#' @importFrom ggplot2 ggplot stat_density_2d aes scale_y_continuous scale_x_reverse geom_tile scale_color_gradient geom_path geom_point annotate theme_classic element_blank element_text labs theme
#'
#' @keywords internal
#' @noRd
#'
evaluate_model <- function(Fit., Obs., th = NULL, rep = 1000, best.th = c("similarity", "accuracy"),
                          main = "", plot = T){

  best.th = best.th[1]

  p = Fit.[Fit. != 0]
  a = Fit.[Fit. == 0]
  np <- length(p)
  na <- length(a)
  N <- na + np

  eval <- list()
  cor <- try( cor.test(Fit., Obs.), silent=TRUE )
  if (class(cor) != 'try-error') {
    p.cor <- cor$p.value
    cor <- cor$estimate
  }
  if (is.null(th)) {
    th <- sort(unique(round(c(a,p), 5)))
    th <- c(th, th[length(th)] + 0.0001)
  }
  else {
    th <- sort(as.vector(th))
  }
  res. <- matrix(ncol=6, nrow=length(th))
  colnames(res.) <- c('tp', 'fp', 'fn', 'tn', 'jaccard.sim', 'threshold')
  jacc <- NULL
  res.n = list()
  for (n in 1:length(th)) {
    Fit.th = Fit.
    Fit.th[Fit.th >= th[n]] = 1
    Fit.th[Fit.th < th[n]] = 0
    res.[n,1] <- length(which(Fit.th[which(Obs. == 1)] == 1))  # a  true positives
    res.[n,2] <- length(which(Fit.th[which(Obs. == 0)] == 1))  # b  false positives
    res.[n,3] <- length(which(Fit.th[which(Obs. == 1)] == 0))     # c  false negatives
    res.[n,4] <- length(which(Fit.th[which(Obs. == 0)] == 0))     # d  true negatives
    res.[n,5] <- 1 - try( prabclus::jaccard(as.matrix(cbind(Fit.th, Obs.))) , silent=TRUE )[1,2]
    res.[n,6] <- th[n]

    res.n[[n]] <- matrix(ncol=6, nrow=rep)
    colnames(res.n[[n]]) <- c('tp', 'fp', 'fn', 'tn', 'jaccard.sim', 'threshold')
    for (r in 1:rep) {
      Fit.n <- rbinom(length(Fit.th), 1, length(Fit.th[Fit.th == 1])/length(Fit.th))
      res.n[[n]][r,1] <- length(which(Fit.n[which(Obs. == 1)] == 1))  # a  true positives
      res.n[[n]][r,2] <- length(which(Fit.n[which(Obs. == 0)] == 1))  # b  false positives
      res.n[[n]][r,3] <- length(which(Fit.n[which(Obs. == 1)] == 0))     # c  false negatives
      res.n[[n]][r,4] <- length(which(Fit.n[which(Obs. == 0)] == 0))     # d  true negatives
      res.n[[n]][r,5] <- 1 - try(prabclus::jaccard(as.matrix(cbind(Fit.n, Obs.))) , silent=TRUE )[1,2]
      res.n[[n]][r,6] <- th[n]
    }
  }
  res.n = ldply(res.n, data.frame)
  res. <- as.data.frame(res.)
  res.n[,"TPR"] = res.n[,1] / (res.n[,1] + res.n[,3])
  res.n[!is.finite(res.n[,"TPR"]),"TPR"] = 0
  res.n[,"TNR"] = res.n[,4] / (res.n[,4] + res.n[,2])
  res.n[!is.finite(res.n[,"TNR"]),"TNR"] = 0
  res.n[,"ACC"] = (res.n[,1] + res.n[,4]) / N
  res.n[!is.finite(res.n[,"ACC"]),"ACC"] = 0
  res.n$test = "random"

  res.[,"TPR"] = res.[,1] / (res.[,1] + res.[,3])
  res.[!is.finite(res.[,"TPR"]),"TPR"] = 0
  res.[,"TNR"] = res.[,4] / (res.[,4] + res.[,2])
  res.[!is.finite(res.[,"TNR"]),"TNR"] = 0
  res.[,"ACC"] = (res.[,1] + res.[,4]) / N
  res.[!is.finite(res.[,"ACC"]),"ACC"] = 0
  res.$test = "predicted"

  if( best.th == "similarity") {
    dist <- sapply(1:length(th), function(i) res.[res.$threshold == th[i], "jaccard.sim"] - mean(res.n[res.n$threshold == th[i], "jaccard.sim"], na.rm = T))
    dist <- (dist+1)/(max(dist)+1)
    p.value <- sapply(1:length(th), function(i) (sum(res.n[res.n$threshold == th[i], "jaccard.sim"] >=
                                                       res.[res.$threshold == th[i], "jaccard.sim"]) + 1)/(length(res.n[res.n$threshold == th[i]
                                                                                                                        , "jaccard.sim"]) + 1))
    score = rowSums(cbind(res.[,"jaccard.sim"], dist), na.rm = T)
    score[is.na(score)] = 0
    THR = th[which.max(score)]
    p.value <- p.value[which.max(score)]
  }
  if( best.th == "accuracy") {
    dist <- sapply(1:length(th), function(i) res.[res.$threshold == th[i], "ACC"] - mean(res.n[res.n$threshold == th[i], "ACC"], na.rm = T))
    dist <- (dist+1)/(max(dist)+1)
    p.value <- sapply(1:length(th), function(i) (sum(res.n[res.n$threshold == th[i], "ACC"] >=
                                                       res.[res.$threshold == th[i], "ACC"]) + 1)/(length(res.n[res.n$threshold == th[i]
                                                                                                                , "ACC"]) + 1))
    score = rowSums(cbind(res.[,"ACC"], dist), na.rm = T)
    score[is.na(score)] = 0
    THR = th[which.max(score)]
    p.value <- p.value[which.max(score)]
  }
  th.index = which.max(score)
  eval$confusion = rbind(res., res.n)
  A = res.[th.index,1]
  B = res.[th.index,2]
  C = res.[th.index,3]
  D = res.[th.index,4]
  eval$n <- c(np = A+B, na = C+D)
  # Accuracy
  ACC = (A + D) / N
  # True Positive Rate
  TPR = A / (A + C)
  if(is.na(TPR)){TPR = 0}
  # True Negative Rate
  TNR = D / (B + D)
  if(is.na(TNR)){TNR = 0}
  # Positive Predictive Value
  PPV = A/(A + B)
  if(is.na(PPV)){PPV = 0}
  # Negative Predictive Value
  NPV = D/(C + D)
  if(is.na(NPV)){NPV = 0}
  # OR
  OR = (A*D)/(C*B)
  if(is.na(OR)){OR = 0}
  # Jaccard similarity
  JACC = res.[th.index,5]
  # AUC
  tpr <- rev(res.[,"TPR"])
  tnr <- rev(1 - res.[,"TNR"])
  dFPR <- c(diff(tnr), 0)
  dTPR <- c(diff(tpr), 0)
  AUC =  sum(tpr * dFPR) + sum(dTPR * dFPR)/2
  # kappa
  prA = (A+D)/N
  prY = (A+B)/N * (A+C)/N
  prN = (C+D)/N * (B+D)/N
  prE = prY + prN
  kappa = (prA - prE) / (1-prE)
  # TSS
  TSS = TPR + TNR - 1

  if(TSS<0){TSS = 0}
  tab <- c( cor, p.cor, JACC, TPR, TNR, TSS, ACC, AUC, kappa, PPV, NPV, OR)
  names(tab) <- c("Pearson's correlation", "p.value", "Jaccard Similarity",
                  "TPR" ,  "TNR", "TSS","ACC", "AUC", "kappa", "PPV", "NPV", "OR")
  eval$tab <- tab
  eval$threshold <- c(threshold = THR, p.value = p.value)
  if (plot == TRUE){

    plot.eval <- function() {

      ggplot(res.n, aes(TNR, TPR)) +
      scale_y_continuous("sensitivity", limits = c(0,1), breaks = seq(0,1,0.25), expand = c(0.01,0.01)) +
      scale_x_reverse("specificity", limits = c(1,0), breaks = seq(0,1,0.25), expand = c(0.01,0.01)) +
      geom_tile(fill = "#132B42") +
      stat_density_2d(aes_string(fill = "..level..", col = "..level..", with = FALSE), geom = "polygon",
                      alpha = 0.1, bins = 10) +
      #geom_density_2d(col = "#E69F00" ) +
      #stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE, n = 100) +
      scale_color_gradient(low = "#132B42", high =  "#96B3C9") +
      geom_path(data = res., aes(TNR, TPR), col = "#E69F00", size = 1) +
      geom_point(data = data.frame(x = TNR, y = TPR), aes_string(x = "x", y = "y"), col = "#E69F00", shape = 18, size = 4) +
      annotate("text", x = TNR, y = TPR, label =  paste("italic(p)==", format.pval(p.value)), parse = T, size = 4, hjust = -0.15, vjust = 1)+
      theme_classic() +
      labs(title= gsub("\\s\\s", "\n", paste0(gsub('\\.', ' ', main), "  AUC=", round(AUC,2))), parse = T) +
      theme(legend.position='none',
            #panel.background = element_rect(fill = "#132B42",
            #                                colour = "#132B42",
            #                                size = 0.5, linetype = "solid"),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            title = element_text(color = "black", size = 10, vjust = 0.5, hjust = 0.5),
            axis.title.x = element_text(color = "black", size = 15, vjust = 0.5, hjust = 0.5),
            axis.title.y = element_text(color = "black", size = 15, vjust = 1, hjust = 0.5),
            axis.text.x = element_text(color = "black", size = 12),
            axis.text.y = element_text(color = "black", size = 12 ),
            axis.line = element_blank())
    }

    print(plot.eval())
    #plot(res.$TNR, res.$TPR, main = gsub("\\s", "\n", paste0(main, " AUC=", round(AUC,2))), asp = 0, type = "l",
    #     xlab = "specificity", ylab = "sensitivity", xlim = c(1,0))
    #points(TNR, TPR, pch = 18, bg = "red", col = "red", cex = 2, lwd = 2)
    #text(TNR + 0.05, TPR - 0.05, bquote(italic(p-value) == .(format.pval(p.value))) , cex = 0.8)
    message("\t...done.")
  }
  return(eval)
}


