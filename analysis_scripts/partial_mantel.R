library(phytools)
library(tidyverse)

load("data/pairwise_matrices.Rdata")

# get cluster data to filter by region
locs = read.csv("data/population_data.csv") %>% dplyr::select(popcode, cluster)

north = locs %>% filter(cluster == "north") %>% dplyr::select(popcode)
north = lapply(north, as.character)
north = north$popcode

center = locs %>% filter(cluster == "central") %>% dplyr::select(popcode)
center = lapply(center, as.character)
center = center$popcode

fst.mat.north = fst.mat[north, north]
geodist.mat.north = geodist.mat[north, north]
tave_prism.mat.north = tave_prism.mat[north, north]
ppt_spring_prism.mat.north = ppt_spring_prism.mat[north, north]

fst.mat.center = fst.mat[center, center]
geodist.mat.center = geodist.mat[center, center]
tave_prism.mat.center = tave_prism.mat[center, center]
ppt_spring_prism.mat.center = ppt_spring_prism.mat[center, center]


# all pops
all.mod = multi.mantel(fst.mat, list(geodist.mat, tave_prism.mat, ppt_spring_prism.mat), nperm=1000)
all.mod$probt

north.mod = multi.mantel(fst.mat.north, list(geodist.mat.north, tave_prism.mat.north, ppt_spring_prism.mat.north), nperm=1000)
north.mod$probt

center.mod = multi.mantel(fst.mat.center, list(geodist.mat.center, tave_prism.mat.center, ppt_spring_prism.mat.center), nperm=1000)
center.mod$probt

# in all cases, geographic distance is the only predictor that matters.

