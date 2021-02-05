## code to prepare `DATASET` dataset goes here
library(raster)
library(sp)
library(virtualspecies)

env_data <- getData("worldclim",var="bio",res=10, download = F)
env_data <- env_data[[1:6]]
env_data <- crop(env_data, c(-12, 35, 35, 55))

usethis::use_data(env_data, overwrite = TRUE, compress = "xz")

set.seed(2021)
n_sp = 5
occ_data1 <- data.frame()
occ_data2 <- data.frame()
int_matrix <- matrix(0,5,5)
for (i in 1:n_sp){
  sp <- generateRandomSp(env_data, niche.breadth = sample(c("wide", "narrow"), 1))
  occ = sampleOccurrences(sp, n = sample(1:100, 1))
  occ_data2 = rbind(occ_data2, cbind(occ$sample.points[,1:2], species = paste0("sp", i)))
  n_int = sample(1:5,1)
  for (j in 1:n_int) {
    int_matrix[i,j] = 1
    occ = sampleOccurrences(sp, n = sample(1:100, 1))
    occ_data1 = rbind(occ_data1, cbind(occ$sample.points[,1:2], species = paste0("sp", LETTERS[j])))
  }
}
occ_data1$species = as.factor(occ_data1$species)
occ_data2$species = as.factor(occ_data2$species)
rownames(int_matrix) = levels(occ_data2$species)
colnames(int_matrix) = sapply(1:5, function(ii) paste0("sp", LETTERS[ii]))

occ_data1 = occ_data1[!duplicated(occ_data1),]
occ_data2 = occ_data2[!duplicated(occ_data2),]

rm(n_sp, sp, i, j, occ, n_int)

usethis::use_data(occ_data1, overwrite = TRUE, compress = "xz")
usethis::use_data(occ_data2, overwrite = TRUE, compress = "xz")
usethis::use_data(int_matrix, overwrite = TRUE, compress = "xz")

