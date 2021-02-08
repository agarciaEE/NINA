context("Testing")

test_that("Succes", {

  library(NINA)


  EN_sp1 = EN_model(env_data, occ_data1, cluster = "env", n.clus = 5, eval = T)

  expect_equal(class(EN_sp1), "NINA")

  pred = EN_sp1$pred.dis
  expect_equal(nrow(pred), nrow(na.exclude(raster::getValues(env_data))))

  expect_equal(length(EN_sp1$z.mod), 5)
  expect_equal(all(unlist(sapply(EN_sp1$z.mod, names)) %in% occ_data1$species), TRUE)

  EN_sp2 = EN_model(env_data, occ_data2, cluster = "env", n.clus = 5, bootstraps = 2, assemble.models = T, split.data = T, combine.clusters = T)

  expect_equal(ncol(EN_sp2$clus), 3)
  expect_equal(levels(EN_sp2$clus[,3]), LETTERS[1:5])

  pred = EN_sp2$pred.dis
  expect_equal(nrow(pred), nrow(na.exclude(raster::getValues(env_data))))

  expect_equal(length(EN_sp2$z.mod), 5)
  expect_equal(raster::nlayers(EN_sp2$maps), length(levels(occ_data2$species)))
  expect_equal(all(unlist(sapply(EN_sp1$z.mod, names)) %in% occ_data1$species), TRUE)
})


