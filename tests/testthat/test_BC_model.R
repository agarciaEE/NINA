context("Testing")

test_that("Succes", {

  library(NINA)

  g1_EN = EN_model(env_data, occ_data1, res = NULL, path = "./", project.name	= "ENspG1_test", cluster = "env", n.clus = 5)
  g2_EN = EN_model(env_data, occ_data2, res = NULL, path = "./", project.name	= "ENspG2_test")

  g2_BC <- BC_model(g2_EN, g1_EN, A.matrix = int_matrix, C.matrix = NULL, type = "global")

  expect_equal(class(g2_BC), "NINA")

  pred = g2_BC$pred.dis
  expect_equal(nrow(pred), nrow(na.exclude(raster::getValues(env_data))))

  expect_equal(length(g2_BC$z.mod), length(levels(occ_data2$species)))
})


