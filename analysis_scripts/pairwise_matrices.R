library(tidyverse)
library(geosphere)
library(cowplot)

# generate pairwise distances of:
# 1. FST
# 2. temperature differences
# 3. precip differences
# 4. geographic distance

# prepare the data in two ways: a set of pairwise matrices (for bedassle) and a tall dataframe (for plots)

# pairwise Fst ------------------------------------------------------------

# read in fst, this is in a matrix format already
fst = read_csv("data/fst_wc.csv")
fst.mat = data.matrix(column_to_rownames(fst, var = "X1"), rownames.force = TRUE)
fst.half = fst.mat
fst.half[upper.tri(fst.half)]  = NA
fst.half.df = data.frame(fst.half) %>% rownames_to_column()

# make this into a tall data frame 
fst.tall = fst.half.df %>% 
  gather("pop2", "fst", 2:33) %>% 
  select(rowname, pop2, fst) %>% 
  filter(!is.na(fst)) %>% 
  rename(pop1 = rowname)

# check out the fsts
hist(fst.tall$fst)


# pairwise geographic distance --------------------------------------------

# pull in spatial and climatic data
pop.info = read_csv("data/population_data.csv") %>% 
  arrange(popcode)

# make vector of ids
ids = pop.info$popcode

# make pairwise distance matrix in km
geodist.mat = distm(pop.info[,c('long','lat')], pop.info[,c('long','lat')], fun=distGeo)/1000

# append ids to matrix
rownames(geodist.mat) = ids
colnames(geodist.mat) = ids

geodist.half = geodist.mat
geodist.half[upper.tri(geodist.half, diag = TRUE)] = NA
geodist.half.df = data.frame(geodist.half) %>% rownames_to_column()

# now make a tall dataframe of pairwise geographic distances
geodist.tall = geodist.half.df %>% 
  gather("pop2", "distance_km", 2:33) %>% 
  select(rowname, pop2, distance_km) %>% 
  filter(!is.na(distance_km)) %>% 
  rename(pop1 = rowname)

# join this to the pairwise fsts. 
pairs1 = left_join(fst.tall, geodist.tall)
summary(pairs1)

# quick look
plot(pairs1$distance_km, pairs1$fst)
plot(geodist.half, fst.half)


# climate differences: tave_prism --------------------------------------------

clim = pop.info %>% 
  select(popcode, tave_5180_sep_jul, ppt_mm_5180_apr_jul)

tave_prism.mat = abs(outer(clim$tave_5180_sep_jul, clim$tave_5180_sep_jul, "-"))

# append ids
rownames(tave_prism.mat) = ids
colnames(tave_prism.mat) = ids

tave_prism.half = tave_prism.mat
tave_prism.half[upper.tri(tave_prism.half, diag = TRUE)] = NA
tave_prism.half.df = data.frame(tave_prism.half) %>% rownames_to_column()

# now make a tall dataframe of pairwise geographic distances
tave_prism.tall = tave_prism.half.df %>% 
  gather("pop2", "tave_prism_diff", 2:33) %>% 
  select(rowname, pop2, tave_prism_diff) %>% 
  filter(!is.na(tave_prism_diff)) %>% 
  rename(pop1 = rowname) %>% 
  mutate(tave_prism_diff = abs(tave_prism_diff))

# join this to the pairwise fsts. 
pairs2 = left_join(pairs1, tave_prism.tall)
summary(pairs2)

# quick look
plot(pairs2$distance_km, pairs2$tave_prism_diff)
plot(geodist.half, tave_prism.half)


# climate differences: ppt_spring_prism --------------------------------------------

ppt_spring_prism.mat = abs(outer(clim$ppt_mm_5180_apr_jul, clim$ppt_mm_5180_apr_jul, "-"))

# append ids
rownames(ppt_spring_prism.mat) = ids
colnames(ppt_spring_prism.mat) = ids

ppt_spring_prism.half = ppt_spring_prism.mat
ppt_spring_prism.half[upper.tri(ppt_spring_prism.half, diag = TRUE)] = NA
ppt_spring_prism.half.df = data.frame(ppt_spring_prism.half) %>% rownames_to_column()

# now make a tall dataframe of pairwise geographic distances
ppt_spring_prism.tall = ppt_spring_prism.half.df %>% 
  gather("pop2", "ppt_spring_prism_diff", 2:33) %>% 
  select(rowname, pop2, ppt_spring_prism_diff) %>% 
  filter(!is.na(ppt_spring_prism_diff)) %>% 
  rename(pop1 = rowname) %>% 
  mutate(ppt_spring_prism_diff = abs(ppt_spring_prism_diff))

# join this to the pairwise fsts. 
pairs3 = left_join(pairs2, ppt_spring_prism.tall)
summary(pairs3)

# quick look
plot(pairs3$distance_km, pairs3$ppt_spring_prism_diff)
plot(geodist.mat, ppt_spring_prism.mat)

summary(pairs3)


# save matrices and dataframes -------------------------------------------

write.csv(pairs3, "data/pairwise_diffs.csv", row.names = FALSE)

save(fst.mat, geodist.mat, tave_prism.mat, ppt_spring_prism.mat, file = "data/pairwise_matrices.Rdata")

# plot fst against geographic distance, colored by precip
ggplot(pairs3, aes(distance_km, fst, color = ppt_spring_prism_diff)) +
  geom_point(size = 2) +
  labs(x = "Geographic distance (km)", y = expression(F["st"])) +
  # scale_colour_gradientn(colors = c("blue", "red", "orange"), values = c(0, 0.5, 1)) +
  scale_colour_gradientn(colors = rainbow(10)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(), text = element_text(size = 14))

# plot fst against geographic distance, colored by temp
ggplot(pairs3, aes(distance_km, fst, color = tave_prism_diff)) +
  geom_point(size = 2) +
  labs(x = "Geographic distance (km)", y = expression(F["st"])) +
  # scale_colour_gradientn(colors = c("blue", "red", "orange"), values = c(0, 0.5, 1)) +
  scale_colour_gradientn(colors = rainbow(10)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(), text = element_text(size = 14))



