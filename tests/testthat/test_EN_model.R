context("Testing")

test_that("Succes", {

  library(NINA)


  EN_sp1 = EN_model(env_data, occ_data1, cluster = "env", n.clus = 5)
  EN_sp2 = EN_model(env_data, occ_data2, cluster = "env", n.clus = 5, bootstraps = 2, assemble.models = T, split.data = T)
  plot(EN_sp2)
  EN_sp2$maps
  cluster_regions(EN_sp2$pred.dis, n.clus = 5)
  points(occ_data1[occ_data1$species == "spA",1:2], pch = 19)
  expect_equal(class(EN_sp1), "NINA")

  pred = EN_sp1$pred.dis
  expect_equal(nrow(pred), nrow(na.exclude(raster::getValues(env_data))))

  expect_equal(length(EN_sp1$z.mod), 5)
  expect_equal(all(unlist(sapply(EN_sp1$z.mod, names)) %in% occ_data1$species), TRUE)
})


