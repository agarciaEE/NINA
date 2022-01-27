## code to prepare `DATASET` dataset goes here
library(raster)
library(sp)
library(virtualspecies)

env_data <- getData("worldclim",var="bio",res=10)
env_data <- env_data[[1:6]]
env_data <- crop(env_data, c(-12, 35, 35, 55))

set.seed(2021)
n_sp1 = 5
n_sp2 = 5

int_matrix <- matrix(sample(c(0,1), n_sp1*n_sp2, replace = T),n_sp2,n_sp1)
rownames(int_matrix) = sapply(1:n_sp2, function(i) paste0("sp", i))
colnames(int_matrix) = sapply(1:n_sp1, function(i) paste0("sp", LETTERS[i]))

occ_data1 <- data.frame()
sp1.l <- list()
for (i in 1:n_sp1){
  sp1 <- generateRandomSp(env_data, niche.breadth = sample(c("wide", "narrow"), 1))
  occ = sampleOccurrences(sp1, n = sample(10:100, 1))
  th1 <- as.numeric(sp1$PA.conversion["beta"])
  sp1$probability.of.occurrence[sp1$probability.of.occurrence < th1] = NA
  sp1.l[[i]] <- sp1$probability.of.occurrence
  occ_data1 = rbind(occ_data1, cbind(occ$sample.points[,1:2], species = as.factor(paste0("sp", LETTERS[i]))))
}

occ_data2 <- data.frame()
for (i in 1:n_sp2){
  Xvar <- which(int_matrix[i,] == 1)
  biomask <-  sum(stack(sp1.l[Xvar]), na.rm = T)
  biomask <- raster::mask(env_data, biomask)
  sp2 <- generateRandomSp(biomask, niche.breadth = sample(c("wide", "narrow"), 1))
  occ = sampleOccurrences(sp2, n = sample(10:100, 1))
  occ_data2 = rbind(occ_data2, cbind(occ$sample.points[,1:2], species = as.factor(paste0("sp", i))))
}

occ_data1 = occ_data1[!duplicated(occ_data1),]
occ_data2 = occ_data2[!duplicated(occ_data2),]

rm(n_sp1, n_sp2, sp1, sp2, sp1.l, th1,  biomask, Xvar, i, occ)

usethis::use_data(env_data, overwrite = TRUE, compress = "xz")
usethis::use_data(occ_data1, overwrite = TRUE, compress = "xz")
usethis::use_data(occ_data2, overwrite = TRUE, compress = "xz")
usethis::use_data(int_matrix, overwrite = TRUE, compress = "xz")

