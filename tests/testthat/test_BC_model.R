context("Testing")

test_that("Succes", {

  library(NINA)

  g1_EN = EN_model(env_data, occ_data1,  cluster = "env", n.clus = 5)
  g2_EN = EN_model(env_data, occ_data2,  cluster = g1_EN$clus)

  g2_BC <- BC_model(g2_EN, g1_EN, A.matrix = int_matrix, C.matrix = NULL, type = "region")

  expect_equal(class(g2_BC), c("NINA", "BCmodel"))

  expect_equal(raster::nlayers(g2_BC$maps), length(levels(occ_data2$species)))

  expect_equal(class(g2_BC$pca), c("pca", "dudi"))

  pred = g2_BC$pred.dis
  expect_equal(nrow(pred), nrow(na.exclude(raster::getValues(env_data))))

  expect_equal(length(g2_BC$z.mod), 5)

  expect_equal(length(g2_BC$w), 5)

  g1_EN = EN_model(env_data, occ_data1)
  g2_EN = EN_model(env_data, occ_data2)
  g2_BC <- BC_model(g2_EN, g1_EN, A.matrix = int_matrix, C.matrix = NULL, type = "global")

  expect_equal(length(g2_BC$z.mod), 5)

  expect_equal(length(g2_BC$w), 5)


})


