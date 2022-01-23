context("Testing")
library(plyr)
test_that("Succes", {

  library(NINA)

  g1_EN = EN_model(env_data, occ_data1,  cluster = "env", n.clus = 5, combine.clusters = T)
  g2_EN = EN_model(env_data, occ_data2,  cluster = g1_EN$clus, combine.clusters = T)

  g2_BC <- BC_model(g2_EN, g1_EN, D = 0, A.matrix = int_matrix, C.matrix = NULL, type = "region")
  g2_EC <- EC_model(g2_BC, D = 0, type = "region")

  g2_NP <- niche_parameters(g2_BC, type = "region", rnd.test = T, rep = 10)
  test <- function(){
    print(g2_NP)
    expect_equal(class(g2_NP), c("NINA", "comp"))

    expect_equal(all(names(g2_BC$maps) %in% levels(g2_NP$np$species)), TRUE)

    expect_equal(class(g2_NP$np), "data.frame")
    expect_equal(class(g2_NP$sim), "list")
    expect_equal(class(g2_NP$pvalue), "list")

    expect_equal(names(g2_NP$sim), levels(g2_NP$np$species))
    expect_equal(names(g2_NP$pvalue), levels(g2_NP$np$species))
  }
  test()
  #g2_NP<- niche_parameters(g2_EC, type = "region", rnd.test = T, rep = 10)
  #test()
  #g2_NP<- niche_parameters(g2_EN, g2_EC, type = "region", rnd.test = T, rep = 10)
  #test()
  #g2_NP <- niche_parameters(g2_EN, g2_BC, type = "region", rnd.test = T, rep = 10)
  #test()
  #g2_NP <- niche_parameters(g2_BC, g2_EN, type = "region", rnd.test = T, rep = 10)
  #test()
  #g2_NP <- niche_parameters(g2_BC, g2_EC, type = "region", rnd.test = T, rep = 10)
  #test()
  #g2_NP <- niche_parameters(g2_EC, g2_BC, type = "region", rnd.test = T, rep = 10)
  #test()
  #g2_NP <- niche_parameters(g2_EC, g2_EN, type = "region", rnd.test = T, rep = 10)
  #test()

  #ncomp = niche_comparison(g2_BC$z.mod$A$sp1, g2_EC$t.mod$A$sp1, centroid.w = T, rnd.test = T, rep = 10)

  #expect_equal(class(ncomp), c("NINA", "metrics"))

  #expect_equal(unique(unlist(lapply(ncomp, class))), "list")

  #expect_equal(unique(unlist(lapply(ncomp, names))), c("obs", "sim", "pvalue", "rep", "weights"))

  #np = niche_position(g2_BC$z.mod[[1]][[1]], type = "multimodal", method = "mean", quantile = 0.95)

  #expect_equal(class(np), "data.frame")

  #expect_equal(ncol(np), 2)

  ##########
  #g1_EN = EN_model(env_data, occ_data1)
  #g2_EN = EN_model(env_data, occ_data2)

  #g2_BC <- BC_model(g2_EN, g1_EN, D = 0, A.matrix = int_matrix, C.matrix = NULL, type = "global")
  #g2_EC <- EC_model(g2_BC, D = 0, type = "global", method = "densisties")

  #g2_NP <- niche_parameters(g2_BC,  type = "global", rnd.test = T, rep = 10)
  #test()
  #g2_NP <- niche_parameters(g2_EC,  type = "global", rnd.test = T, rep = 10)
  #test()
  #g2_NP <- niche_parameters(g2_EN, g2_EC, type = "global", rnd.test = T, rep = 10)
  #test()
  #g2_NP <- niche_parameters(g2_EN, g2_BC, type = "global", rnd.test = T, rep = 10)
  #test()
  #g2_NP <- niche_parameters(g2_BC, g2_EN, type = "global", rnd.test = T, rep = 10)
  #test()
  #g2_NP <- niche_parameters(g2_BC, g2_EC, type = "global", rnd.test = T, rep = 10)
  #test()
  #g2_NP <- niche_parameters(g2_EC, g2_BC, type = "global", rnd.test = T, rep = 10)
  #test()
  #g2_NP <- niche_parameters(g2_EC, g2_EN, type = "global", rnd.test = T, rep = 10)
  #test()

})
