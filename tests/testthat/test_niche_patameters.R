context("Testing")

test_that("Succes", {

  library(NINA)

  g1_EN = EN_model(env_data, occ_data1,  cluster = "env", n.clus = 5)
  g2_EN = EN_model(env_data, occ_data2,  cluster = g1_EN$clus)

  g2_BC <- BC_model(g2_EN, g1_EN, D = 0, A.matrix = int_matrix, C.matrix = NULL, type = "region")
  g2_EC <- EC_model(g2_BC, D = 0, type = "region")

  g2_NP <- niche_parameters(g2_BC, g2_EC, type = "region")

  expect_equal(class(g2_NP), "data.frame")

  expect_equal(all(names(g2_BC$maps) %in% g2_NP$species), TRUE)

  np = niche_position(g2_BC$z.mod[[1]][[1]], type = "multimodal", method = "mean", quantile = 0.95)

  expect_equal(class(np), "data.frame")

  expect_equal(ncol(np), 2)

})
