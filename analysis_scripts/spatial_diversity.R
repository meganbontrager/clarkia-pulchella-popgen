library(tidyverse)
library(sp)
library(grid)
library(cowplot)

# spatial patterns in genetic diversity ----------------------------------------------------

# equal area projection for the pacific northwest
prj.aea = "+proj=aea +lat_1=41 +lat_2=47 +lat_0=44 +lon_0=-120 +x_0=0 +y_0=0 +units=km"

# read in locations, add equal area projection
locs = read_csv("data/population_data.csv")
locs1 = locs
coordinates(locs1) = ~long + lat
proj4string(locs1) = CRS("+proj=longlat +ellps=WGS84")
locs_t = spTransform(locs1, CRS(prj.aea)) 
locs_p = data.frame(locs_t)

# read in all locations
all = read.csv("data/all_localities_prism.csv", na.strings = c("", "NA")) %>%
  dplyr::filter(!is.na(id)) %>% 
  separate(date, into = c("month", "year2"), sep = "-") 
all1 = all %>% dplyr::select(long, lat) %>% distinct()
coordinates(all1) = ~long + lat
proj4string(all1) = CRS("+proj=longlat +ellps=WGS84")
all_t = spTransform(all1, CRS(prj.aea)) 
all_p = data.frame(all_t)

# get state outlines
states = map_data("state")
states1 = states
coordinates(states1) =~long+lat
proj4string(states1) = CRS("+proj=longlat +ellps=WGS84")
states_t = spTransform(states1, CRS(prj.aea)) 
states_p = data.frame(states_t)

# get canada outlines
canada = map_data("world", region = "canada*")
canada1 = canada
coordinates(canada1) =~long+lat
proj4string(canada1) = CRS("+proj=longlat +ellps=WGS84")
canada_t = spTransform(canada1, CRS(prj.aea)) 
canada_p = data.frame(canada_t)

# make range that is an MCP of all points
hull = all_p[chull(all_p$long, all_p$lat),]

# distance to the range edge
dist_to_edge = read.csv("data/pops_spatial_stats.csv")
gen_div = read.csv("data/pops_stats.csv")
sequenced = read.csv("data/representation.csv")

div = left_join(gen_div, dist_to_edge, by = c("pop" = "popcode"))
div = left_join(div, locs, by = c("pop" = "popcode"))
div = left_join(div, sequenced)
summary(div)

locs_he = left_join(locs_p, gen_div, by = c("popcode" = "pop"))

# does genetic diversity scale with sequence representation?
plot(div$exp_het ~ div$prop)
summary(lm(div$exp_het ~ div$prop))
cor.test(div$exp_het , div$prop)
# nope

# does geneitc diversity scale with latitude?
plot(div$exp_het ~ div$lat)
summary(lm(exp_het ~ lat, data = div))
# yep
# better described by quadratic?
summary(lm(exp_het ~ poly(lat,2), data = div))
# nope

# does genetic diversity scale with latitude?
plot(div$exp_het ~ div$lat)
summary(lm(exp_het ~ dist_to_edge, data = div))
# nope

# create manual adjustted lat longs so that colored points won't cover each other up so much.
# there must be a better way to do this! but I couldn't find one

locs_he$lat_dodge = c(553.57320, 299.68122, 266.76489, 251.07828, 194.69713, 549.89614, 530.10430, 525.37915, 505.43337, 396.07654, 376.50696, 347.52045, 326.49835, 561.11712, 568.17798, 550.79008, 548.05821, 547.13914, 541.91217, 534.52736, 494.93914, 493.02344, 486.94001, 395.83769, 369.20010, 340.72271, 331.46416, 309.55671, 252.70095,  52.28867,  42.81374, -74.39204)

locs_he$long_dodge = c(74.13265, 468.45918, 185.12892, 174.02597, 135.69023, 109.23145,  86.08754, 122.92373,  80.70117, 394.00691, 242.77959, 229.62006, 426.57386, 32.44419, 69.74939, 176.90362, 145.02152, 128.96014,  90.09049, 133.75600,  92.85738,  73.30780, 309.15411, 250.75523, 320.88869, 205.49929, 459.65691, 251.49386, 193.64248, -56.48335, -41.40437, 221.15653)

# make a map
het_map =
  ggplot() + 
  coord_quickmap(xlim = c(min(hull$long) - 200,
                          max(hull$long) + 100),
                 ylim = c(min(hull$lat) - 20,
                          max(hull$lat) + 20)) +
  geom_polygon(data = canada_p, aes(x = long, y = lat, group = group), 
               color = "black", fill = "transparent", alpha = 0) +
  geom_polygon(data = states_p, aes(x = long, y = lat, group = group), 
               color = "black", fill = "transparent", alpha = 0) +
  geom_polygon(data = hull, aes(x = long, y = lat), color = "black", fill = "transparent", linetype = "dashed") +
  geom_point(data = locs_he, aes(x = long_dodge, y = lat_dodge, fill = exp_het), size = 3.5, shape = 21, color = "black") +
  # geom_point(data = filter(locs_he, lat > 100), aes(x = long, y = lat, color = exp_het), size = 5) +
  # scale_color_continuous(low = "darkgoldenrod1", high = "mediumblue", name = "Expected\nheterozygosity") +
  scale_fill_continuous(low = "thistle", high = "deeppink4", name = "Expected\nheterozygosity") + 
  # geom_label_repel(data = locs_p, aes(x = long, y = lat, label = round(lat,3)), size = 3, segment.size  = 0.5, box.padding = 0.15, point.padding = 0.15,segment.color = "black", force = 0.5) +
  theme_bw() +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(0.9, 0.2), 
        legend.background = element_rect(fill = "white", color = "black"),
        legend.text = element_text(size = 10))+
  guides(fill = FALSE) +
  guides(fill = guide_colorbar(ticks = FALSE, title.position = "top"))
# ggsave("figs/het_map_only.pdf", height = 5)

het_lat = ggplot() +  
  geom_abline(intercept = -0.264745, slope = 0.010466 ) +
  geom_jitter(data = locs_he, aes(x = lat_old, y = exp_het, fill = exp_het), size = 3.5, shape = 21, color = "black") +
  # geom_point(data = filter(locs_he, lat > 100), aes(x = long, y = lat, color = exp_het), size = 5) +
  # scale_color_continuous(low = "darkgoldenrod1", high = "mediumblue", name = "Expected\nheterozygosity") +
  scale_fill_continuous(low = "thistle", high = "deeppink4", name = "Expected\nheterozygosity") + 
  scale_y_continuous(breaks = c(0.18, 0.20, 0.22, 0.24, 0.26)) +
  theme(plot.background  = element_rect(fill = "white", color = "black")) +
  guides(fill = FALSE) +
  xlab("Latitude (°)") +
  ylab("Expected heterozygosity")
# ggsave("figs/het_plot_only.pdf", height = 4.5, width = 4.5)

# put these two together and save
ggdraw() + 
  draw_plot(het_lat, 0.02, 0.04, 0.48, 0.92) +
  draw_plot(het_map, 0.52, 0.02, 0.46, 0.96) +
  # draw_grob(rectGrob(gp = gpar(col = "black", fill = "transparent")), 0.02, 0.025, 0.47, 0.96) +
  draw_plot_label(c("A", "B"), x = c(0, 0.49), y = c(1, 1))
# ggsave("figs_for_paper/het_map_combo.pdf", height = 5, width = 11.3)


