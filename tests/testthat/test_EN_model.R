context("Testing")

test_that("Succes", {

  library(NINA)


  EN_sp1 = EN_model(env_data, occ_data1, cluster = "env", n.clus = 5, eval = T, combine.clusters = T)

  expect_equal(class(EN_sp1), "NINA")

  pred = EN_sp1$pred.dis
  expect_equal(nrow(pred), nrow(na.exclude(raster::getValues(env_data))))

  expect_equal(length(EN_sp1$z.mod), 5)
  expect_equal(all(unlist(sapply(EN_sp1$z.mod, names)) %in% occ_data1$species), TRUE)

  expect_equal(nrow(EN_sp1$eval$tab), length(levels(occ_data1$species)))

  expect_equal(dim(EN_sp1$eval$threshold), c(length(levels(occ_data1$species)), 2))

  expect_equal(unique(EN_sp1$eval$confusion$test), c("predicted", "random"))

  expect_equal(rowSums(EN_sp1$eval$cases[,2:3]), rowSums(EN_sp1$eval$n))

  EN_sp1 = EN_model(env_data, occ_data1, cluster = "obs", n.clus = 5)

  expect_equal(ncol(EN_sp1$clus), 3)

  EN_sp2 = EN_model(raster::as.data.frame(env_data, xy = T), occ_data2, bootstraps = 2, save.bootstraps = T, save.model = T, assemble.models = T, split.data = T)

  pred = EN_sp2$pred.dis
  expect_equal(nrow(pred), nrow(na.exclude(raster::getValues(env_data))))

  expect_equal(length(EN_sp2$z.mod), 5)
  expect_equal(raster::nlayers(EN_sp2$maps), length(levels(occ_data2$species)))

})


