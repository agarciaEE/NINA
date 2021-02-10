context("Testing")

test_that("Succes", {

  library(NINA)

  g1_EN = EN_model(env_data, occ_data1,  cluster = "env", n.clus = 5)
  g2_EN = EN_model(env_data, occ_data2,  cluster = g1_EN$clus)

  g2_BC <- BC_model(g2_EN, g1_EN, A.matrix = int_matrix, C.matrix = NULL, type = "region")

  g2_EC <- EC_model(g2_BC, type = "region")

  expect_equal(class(g2_EC), "NINA")

  expect_equal(raster::nlayers(g2_EC$maps), length(levels(occ_data2$species)))

  expect_equal(class(g2_EC$pca), c("pca", "dudi"))

  pred = g2_EC$pred.dis
  expect_equal(nrow(pred), nrow(na.exclude(raster::getValues(env_data))))

  expect_equal(length(g2_EC$z.mod), 5)

  expect_equal(length(g2_EC$w), 5)

  g1_EN = EN_model(env_data, occ_data1)
  g2_EN = EN_model(env_data, occ_data2)
  g2_BC <- BC_model(g2_EN, g1_EN, A.matrix = int_matrix, C.matrix = NULL, type = "global")
  g2_EC <- EC_model(g2_BC, type = "global")

  g2_NP <- niche_parameters(g2_BC, g2_EC, type = "global")

  expect_equal(class(g2_NP), "data.frame")

  expect_equal(all(names(g2_BC$maps) %in% g2_NP$species), TRUE)

  expect_equal(length(g2_EC$z.mod), 5)

  expect_equal(length(g2_EC$w), 5)
})
