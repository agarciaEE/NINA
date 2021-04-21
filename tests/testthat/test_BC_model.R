context("Testing")

test_that("Succes", {

  library(NINA)

  g1_EN = EN_model(env_data, occ_data1,  cluster = "env", n.clus = 5)
  g2_EN = EN_model(env_data, occ_data2,  cluster = g1_EN$clus)

  g2_BCcomp <- BC_model(g2_EN, g1_EN, A.matrix = int_matrix, C.matrix = NULL, method = "composition", cor  = F, type = "region")
  g2_BCdens <- BC_model(g2_EN, g1_EN, A.matrix = int_matrix, C.matrix = NULL, method = "denstities", cor  = T, type = "region")
dev.off()
plot(g2_EN$maps$sp2)
plot(g2_BCdens$maps$sp2)
plot(g2_BCcomp$maps$sp2)
plot(g2_BCcomp$z.mod$A$sp1$z.cor)
  plot(g2_BCdens$z.mod$A$sp1$z.cor)
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

  EN_sp1 = EN_model(raster::as.data.frame(env_data, xy = T), occ_data1, bootstraps = 5, save.bootstraps = F, save.model = F, assemble.models = T, split.data = T)$ensemble
  EN_sp2 = EN_model(raster::as.data.frame(env_data, xy = T), occ_data2, bootstraps = 5, save.bootstraps = F, save.model = F, assemble.models = F, split.data = T)

  expect_equal(class(EN_sp2), c("NINA", "modelsList"))

  BC_sp2 <- lapply(EN_sp2, function(i) BC_model(i, EN_sp1, D = 0, A.matrix = int_matrix, C.matrix = NULL, type = "global"))

  expect_equal(class(BC_sp2), c("list"))

  BC_sp2_ensemble = assemble_models(BC_sp2, method = "AUC")

  expect_equal(all(names(BC_sp2[[1]]$maps) %in% unique(occ_data2$species)), TRUE)

  expect_equal(length(BC_sp2), 5)

  expect_equal(class(BC_sp2_ensemble), c("NINA", "BCmodel"))

})


