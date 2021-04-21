context("Testing")

test_that("Succes", {

  library(NINA)

  g1_EN1 = EN_model(env_data, occ_data1,  cluster = "env", n.clus = 5, relative.niche = F, cor = T, eval = T)
  g1_EN2 = EN_model(env_data, occ_data1,  cluster = "env", n.clus = 5, relative.niche = F, cor = F, eval = T)
  g1_EN3 = EN_model(env_data, occ_data1,  cluster = "env", n.clus = 5, relative.niche = T, cor = T, eval = T)
  g1_EN4 = EN_model(env_data, occ_data1,  cluster = "env", n.clus = 5, relative.niche = T, cor = F, eval = T)

  plot(g1_EN1$eval)
  plot(g1_EN2$eval)
  plot(g1_EN3$eval)
  plot(g1_EN4$eval)

  g1_EN = EN_model(env_data, occ_data1,  cluster = "env", n.clus = 5, relative.niche = F, cor = F)
  g2_EN = EN_model(env_data, occ_data2,  cluster = g1_EN$clus, relative.niche = F, cor = F)
  summary(g1_EN)
  print(g1_EN)
  plot(g1_EN)
  g2_BC <- BC_model(g2_EN, g1_EN, A.matrix = int_matrix, C.matrix = NULL, type = "region")

  g2_EC <- EC_model(g2_BC, type = "region")

  expect_equal(class(g2_EC), c("NINA", "ECmodel"))

  expect_equal(raster::nlayers(g2_EC$maps), length(levels(occ_data2$species)))

  expect_equal(class(g2_EC$pca), c("pca", "dudi"))

  pred = g2_EC$pred.dis
  expect_equal(nrow(pred), nrow(na.exclude(raster::getValues(env_data))))

  expect_equal(length(g2_EC$z.mod), 5)

  expect_equal(length(g2_EC$g), 5)

  EN_sp1 = EN_model(raster::as.data.frame(env_data, xy = T), occ_data1, bootstraps = 5, save.bootstraps = F, save.model = F, assemble.models = T, split.data = T)$ensemble
  EN_sp2 = EN_model(raster::as.data.frame(env_data, xy = T), occ_data2, bootstraps = 5, save.bootstraps = F, save.model = F, assemble.models = F, split.data = T)

  expect_equal(class(EN_sp2), c("NINA", "modelsList"))

  EC_sp2 <- lapply(EN_sp2, function(i) EC_model(i, EN_sp1, D = 0, A.matrix = int_matrix, C.matrix = NULL, type = "global"))

  expect_equal(class(EC_sp2), c("list"))

  EC_sp2_ensemble = assemble_models(EC_sp2, method = "Jaccard Similarity")

  expect_equal(all(names(EC_sp2[[1]]$maps) %in% unique(occ_data2$species)), TRUE)

  expect_equal(length(EC_sp2), 5)

  expect_equal(class(EC_sp2_ensemble), c("NINA", "ECmodel"))
})
