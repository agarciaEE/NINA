context("Testing")

test_that("Succes", {

  set.seed(1)
  res = raster_projection(data.frame(x = 1:100, y = 1:100, z = rnorm(100), w =  rbinom(100, 1,.5)))
plot(res)
  expect_equal(class(res)[1], "RasterStack")

  cm = raster::as.matrix(res[[1]])
  expect_equal(dim(cm), c(100, 100))

  expect_equal(raster::nlayers(res), 2)
})
