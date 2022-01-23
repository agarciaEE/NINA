context("Testing")

test_that("Succes", {

  library(NINA)

  g1_EN = EN_model(env_data, occ_data1,  cluster = "env", n.clus = 5)
  g2_EN = EN_model(env_data, occ_data2,  cluster = g1_EN$clus)

  g2_BC <- BC_model(g2_EN, g1_EN, A.matrix = int_matrix, C.matrix = NULL, method = "densities", type = "region")

  eval = models_evaluation(g1_EN$pred.dis, g1_EN$obs, env_data, plot = F)
  print(eval)
  plot(eval)
  expect_equal(class(eval), c("NINA", "eval"))

  expect_equal(nrow(eval$tab), length(levels(occ_data1$species)))

  expect_equal(dim(eval$threshold), c(length(levels(occ_data1$species)), 2))

  expect_equal(unique(eval$confusion$test), c("predicted", "random"))

  expect_equal(rowSums(eval$cases[,2:3]), rowSums(eval$n))

  eval = models_evaluation(g2_BC, plot = T, best.th = "similarity")

  expect_equal(class(eval), c("NINA", "eval"))

  expect_equal(nrow(eval$tab), length(levels(occ_data2$species)))

  expect_equal(dim(eval$threshold), c(length(levels(occ_data2$species)), 2))

  expect_equal(unique(eval$confusion$test), c("predicted", "random"))

  expect_equal(rowSums(eval$cases[,2:3]), rowSums(eval$n))

  Pseudo_abs <- sample_pseudoabsences(occ_data1, env_data, spsNames = NULL, plot = T)

  expect_equal(class(Pseudo_abs), "list")

  expect_equal(nrow(Pseudo_abs$tab), length(levels(occ_data2$species)))
})
