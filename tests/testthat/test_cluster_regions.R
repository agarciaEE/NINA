context("Testing")

test_that("Succes", {

  library(NINA)

  env_clus = cluster_regions(env_data, n.clus = 5, plot = T)

  expect_equal(class(env_clus), "data.frame")

  expect_equal(nrow(env_clus), nrow(na.exclude(raster::getValues(env_data))))

  expect_equal(length(unique(env_clus$cluster)), 5)
})

