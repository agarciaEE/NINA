context("Testing")

test_that("Succes", {

  library(NINA)

  g1_EN = EN_model(env_data, occ_data1,  cluster = "env", n.clus = 5)

  eval = models_evaluation(g1_EN$pred.dis, g1_EN$obs, env_data, plot = T)

  expect_equal(class(eval), "NINA")

  expect_equal(nrow(eval$tab), length(levels(occ_data1$species)))

  expect_equal(dim(eval$threshold), c(length(levels(occ_data1$species)), 2))

  expect_equal(unique(eval$confusion$test), c("predicted", "random"))

  expect_equal(rowSums(eval$cases[,2:3]), rowSums(eval$n))

})
