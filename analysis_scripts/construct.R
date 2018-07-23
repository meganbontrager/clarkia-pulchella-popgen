# running conStruct on clarkia pulchella data


# libraries and data ------------------------------------------------------

library(devtools)
install_github("gbradburd/conStruct/code/conStruct", build_vignettes=TRUE)
library(conStruct)
library(tidyverse)

vignette(topic = "format-data", package = "conStruct")
vignette(topic = "run-conStruct", package = "conStruct")

# load allele frequency matrix
load("popgen_tables/construct_matrix.Rdata")
# object is called construct.pop.matrix

# load lat/longs
locs = read.csv("data/population_data.csv") %>% select(popcode, long, lat) %>% column_to_rownames("popcode")
locs = as.matrix(locs)

# load geographic distance matrix
load("data/pairwise_matrices.Rdata")
# object is called geodist.mat



# run construct -----------------------------------------------------------

# run spatial and non-spatial versions with 1-5 clusters

conStruct(spatial = TRUE,
          K = 1,
          freqs = construct.pop.matrix,
          geoDist = geodist.mat,
          coords = locs,
          prefix = "spk1")

conStruct(spatial = TRUE,
          K = 2,
          freqs = construct.pop.matrix,
          geoDist = geodist.mat,
          coords = locs,
          prefix = "spk2")

conStruct(spatial = TRUE,
          K = 3,
          freqs = construct.pop.matrix,
          geoDist = geodist.mat,
          coords = locs,
          prefix = "spk3")

conStruct(spatial = TRUE,
          K = 4,
          freqs = construct.pop.matrix,
          geoDist = geodist.mat,
          coords = locs,
          prefix = "spk4")

conStruct(spatial = TRUE,
          K = 5,
          freqs = construct.pop.matrix,
          geoDist = geodist.mat,
          coords = locs,
          prefix = "spk5")

conStruct(spatial = FALSE,
          K = 1,
          freqs = construct.pop.matrix,
          geoDist = geodist.mat,
          coords = locs,
          prefix = "nsk1")

conStruct(spatial = FALSE,
          K = 2,
          freqs = construct.pop.matrix,
          geoDist = geodist.mat,
          coords = locs,
          prefix = "nsk2")

conStruct(spatial = FALSE,
          K = 3,
          freqs = construct.pop.matrix,
          geoDist = geodist.mat,
          coords = locs,
          prefix = "nsk3")

conStruct(spatial = FALSE,
          K = 4,
          freqs = construct.pop.matrix,
          geoDist = geodist.mat,
          coords = locs,
          prefix = "nsk4")

conStruct(spatial = FALSE,
          K = 5,
          freqs = construct.pop.matrix,
          geoDist = geodist.mat,
          coords = locs,
          prefix = "nsk5")



# run cross validation ----------------------------------------------------

x.validation(train.prop = 0.9, n.reps = 20, K = 1:5, 
             freqs = construct.pop.matrix,
             geoDist = geodist.mat,
             coords = locs,
             n.iter = 1e3,
             prefix = "xval_2_", make.figs = TRUE, save.files = TRUE)



# extract layer contributions ---------------------------------------------

lay_con = data.frame(sp1 = rep(NA, 5))

load("spk1_conStruct.results.Robj")
load("spk1_data.block.Robj")

lay_con$sp1 = c(calculate.layer.contribution(conStruct.results[[1]], data.block), rep(NA, 4))

load("conStruct_output/spk2_conStruct.results.Robj")
load("spk2_data.block.Robj")

lay_con$sp2 = c(calculate.layer.contribution(conStruct.results[[1]], data.block), rep(NA, 3))

load("spk3_conStruct.results.Robj")
load("spk3_data.block.Robj")

lay_con$sp3 = c(calculate.layer.contribution(conStruct.results[[1]], data.block), rep(NA, 2))

load("spk4_conStruct.results.Robj")
load("spk4_data.block.Robj")

lay_con$sp4 = c(calculate.layer.contribution(conStruct.results[[1]], data.block), rep(NA, 1))

load("spk5_conStruct.results.Robj")
load("spk5_data.block.Robj")

lay_con$sp5 = calculate.layer.contribution(conStruct.results[[1]], data.block)

load("nsk1_conStruct.results.Robj")
load("nsk1_data.block.Robj")

lay_con$nsp1 = c(calculate.layer.contribution(conStruct.results[[1]], data.block), rep(NA, 4))

load("nsk2_conStruct.results.Robj")
load("nsk2_data.block.Robj")

lay_con$nsp2 = c(calculate.layer.contribution(conStruct.results[[1]], data.block), rep(NA, 3))

load("nsk3_conStruct.results.Robj")
load("nsk3_data.block.Robj")

lay_con$nsp3 = c(calculate.layer.contribution(conStruct.results[[1]], data.block), rep(NA, 2))

load("nsk4_conStruct.results.Robj")
load("nsk4_data.block.Robj")

lay_con$nsp4 = c(calculate.layer.contribution(conStruct.results[[1]], data.block), rep(NA, 1))

load("nsk5_conStruct.results.Robj")
load("nsk5_data.block.Robj")

lay_con$nsp5 = calculate.layer.contribution(conStruct.results[[1]], data.block)

write.csv(lay_con, "data/construct_layer_contributions.csv", row.names = FALSE)



# construct xval results plot ---------------------------------------------

library(cowplot)

xval = read.csv("data/xval_summary.csv") %>%
  filter(sp_nsp == "sp")

ggplot(xval, aes(x = layers, y = value, group = layers)) +
  geom_boxplot() +
  ylab("Relative explanatory power") +
  xlab("Number of layers")
# ggsave("figs/xval_results.pdf", height = 4, width = 4)


# make map of populations -----------------------------------------------

library(raster)
library(sp)
library(rgdal)
library(stringr)

# load admixture proportions from best fitting model
load("conStruct_output/spk2_conStruct.results.Robj")
admix = conStruct.results$chain_1$MAP$admix.proportions
colnames(admix) = c("layer1", "layer2")

prj.aea = "+proj=aea +lat_1=41 +lat_2=47 +lat_0=44 +lon_0=-120 +x_0=0 +y_0=0 +units=km"

# load lat longs
locs = read_csv("data/population_data.csv") %>% dplyr::select(popcode, lat, long) %>% arrange(popcode)
locs = cbind(locs, admix)

# transform to aea projections
locs1 = locs
coordinates(locs1) = ~long + lat
proj4string(locs1) = CRS("+proj=longlat +ellps=WGS84")
locs_t = spTransform(locs1, CRS(prj.aea)) 
locs_p = data.frame(locs_t)
plot(locs$long, locs$lat)
plot(locs_t)
plot(locs1)

# load all localities
all = read.csv("data/all_localities_prism.csv", na.strings = c("", "NA")) %>%
  dplyr::filter(!is.na(id)) %>% 
  separate(date, into = c("month", "year2"), sep = "-") 
all1 = all %>% dplyr::select(long, lat) %>% distinct()
coordinates(all1) = ~long + lat
proj4string(all1) = CRS("+proj=longlat +ellps=WGS84")
all_t = spTransform(all1, CRS(prj.aea)) 
all_p = data.frame(all_t)

# load states
states = map_data("state")
states1 = states
coordinates(states1) =~long+lat
proj4string(states1) = CRS("+proj=longlat +ellps=WGS84")
states_t = spTransform(states1, CRS(prj.aea)) 
states_p = data.frame(states_t)

# load canada
canada = map_data("world", region = "canada*")
canada1 = canada
coordinates(canada1) =~long+lat
proj4string(canada1) = CRS("+proj=longlat +ellps=WGS84")
canada_t = spTransform(canada1, CRS(prj.aea)) 
canada_p = data.frame(canada_t)

# make an MCP of localities
hull = all_p[chull(all_p$long, all_p$lat),]

# labeled map
ggplot() +
  geom_polygon(data = canada_p, aes(x = long, y = lat, group = group), 
               color = "black", fill = "transparent", alpha = 0) +
  geom_polygon(data = states_p, aes(x = long, y = lat, group = group), 
               color = "black", fill = "transparent", alpha = 0) +
  # geom_text(data = locs_p, aes(x = long, y = lat, label = fig_name), size = 4) 
  geom_polygon(data = hull, aes(x = long, y = lat), color = "black", fill = "transparent", linetype = "dashed") +
  geom_point(data = locs_p, aes(x = long, y = lat), size = 2) +
  theme_bw() +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA))  + 
  coord_quickmap(xlim = c(min(hull$long) - 210,
                          max(hull$long) + 50),
                 ylim = c(min(hull$lat) - 5,
                          max(hull$lat) + 20)) 

# ggsave("figs/admix_base_map.pdf", height = 8, width = 8)

for_pie = locs %>% 
  dplyr::select(fig_name, layer1, layer2, lat) %>% 
  gather(key = "layer", value = "prop", 2:3) %>% 
  arrange(-lat) %>% 
  mutate(fig_name = factor(fig_name, unique(fig_name)))

ggplot(for_pie, aes(x = "", y = prop, fill = layer))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) + 
  facet_wrap(~fig_name) +
  scale_fill_manual(values=c("darkorchid4", "tomato1")) + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(), plot.background = element_blank(), panel.background = element_blank())

# ggsave("figs/pies_to_place.pdf", height = 12, width = 12)

ggplot(for_pie, aes(x = fig_name, y = prop, fill = layer)) +
  scale_fill_manual(values=c("darkorchid4", "tomato1")) + 
  geom_bar(stat = "identity") + 
  guides(fill = FALSE) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())

# ggsave("figs/admix_bars.pdf", height = 2, width = 6)

# manually place pies on map
# result is conStruct_main_fig.pdf
