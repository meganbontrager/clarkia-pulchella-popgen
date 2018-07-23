library(tidyverse)
library(cowplot)

# plot fst vs. temp, geo, dist

# read in pairwise distances
pwise = read.csv("data/pairwise_diffs.csv")

# read in localities for cluster assignments
locs = read_csv("data/population_data.csv")

north = locs$popcode[locs$cluster == "north"]
center = locs$popcode[locs$cluster == "central"]

# iterate, filtering to regional subsets if desired
# pwise = pwise %>% filter(pop1 %in% north & pop2 %in% north)
# pwise = pwise %>% filter(pop1 %in% center & pop2 %in% center)

# plot fst against geographic distance, with no color
fst_geo_ppt =
  ggplot(pwise, aes(distance_km, fst, color = ppt_spring_prism_diff)) +
  geom_point(size = 2) +
  labs(x = "Geographic distance (km)", y = expression(F["ST"])) +
  scale_colour_continuous(low = "greenyellow" , high =  "orangered2", name = "Spring/summer\nprecipitation\ndifference\n(mm)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 14), legend.title = element_text(size = 12), legend.text = element_text(size = 10)) +
  guides(color = guide_colorbar(ticks = FALSE))

fst_geo_temp =
  ggplot(pwise, aes(distance_km, fst, color = tave_prism_diff)) +
  geom_point(size = 2) +
  labs(x = "Geographic distance (km)", y = expression(F["ST"])) +
  scale_colour_continuous(low = "greenyellow" , high =  "orangered2", name = "Annual\ntemperature\ndifference\n(°C)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 14), legend.title = element_text(size = 12), legend.text = element_text(size = 10)) +
  guides(color = guide_colorbar(ticks = FALSE))

fst_ppt =
  ggplot(pwise, aes(ppt_spring_prism_diff, fst, color = distance_km)) +
  geom_point(size = 2) +
  labs(x = "Spring/summer precipitation\ndifference (mm)", y = expression(F["ST"])) +
  scale_colour_continuous(low = "palegreen" , high =  "darkblue", name = "Geographic\ndistance\n(km)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 14), legend.title = element_text(size = 12), legend.text = element_text(size = 10)) +
  guides(color = guide_colorbar(ticks = FALSE))

fst_temp =
  ggplot(pwise, aes(tave_prism_diff, fst, color = distance_km)) +
  geom_point(size = 2) +
  labs(x = "Annual temperature\ndifference (°C)", y = expression(F["ST"])) +
  scale_colour_continuous(low = "palegreen" , high =  "darkblue", name = "Geographic\ndistance\n(km)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 14), legend.title = element_text(size = 12), legend.text = element_text(size = 10)) +
  guides(color = guide_colorbar(ticks = FALSE))

plot_grid(fst_geo_temp, fst_geo_ppt, fst_temp, fst_ppt, align = "hv", labels = c("A", "B", "C", "D"))
# ggsave("figs_for_paper/fst_clim_dist.pdf", height = 7, width = 9.5)
# ggsave("figs_for_paper/fst_clim_dist_north.pdf", height = 7, width = 9.5)
# ggsave("figs_for_paper/fst_clim_dist_center.pdf", height = 7, width = 9.5)

plot(pwise$tave_prism_diff, pwise$distance_km)
plot(pwise$ppt_spring_prism_diff, pwise$distance_km)
# in the north, geodist and ppt are collinear

