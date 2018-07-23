# map and climate plots

library(tidyverse)
library(sp)
library(cowplot)
library(ggrepel)

# make map of populations -----------------------------------------------
# equal area projection for the pacific northwest
prj.aea = "+proj=aea +lat_1=41 +lat_2=47 +lat_0=44 +lon_0=-120 +x_0=0 +y_0=0 +units=km"

# read in locations and convert to AEA
locs = read_csv("data/population_data.csv")
locs1 = locs
coordinates(locs1) = ~long + lat
proj4string(locs1) = CRS("+proj=longlat +ellps=WGS84")
locs_t = spTransform(locs1, CRS(prj.aea)) 
locs_p = data.frame(locs_t)
plot(locs$long, locs$lat)
plot(locs_t)
plot(locs1)

# plot elevation (not included on github, raster is a big file)
# load("data/elevation_8_res.Rdata")
# elevation8
# elev_aea = projectRaster(elevation8, crs = prj.aea)
# plot(elev_aea)
# elev_aea_land = elev_aea
# elev_aea_land[elev_aea_land < 0] = NA
# plot(elev_aea_land)
# elev_aea_land_df = rasterToPoints(elev_aea_land)
# elev_aea_land_df2 = data_frame(elev_aea_land_df)
# pdf("figs_for_paper/elev.pdf", height = 10, width = 10)
# plot(elev_aea_land, col = gray.colors(100, start = 0.9, end = 0))
# dev.off()

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

# map - can comment labels in or out
main =
  ggplot() + 
  coord_quickmap(xlim = c(min(hull$long) - 210,
                          max(hull$long) + 50),
                 ylim = c(min(hull$lat) - 5,
                          max(hull$lat) + 20)) +
  geom_polygon(data = canada_p, aes(x = long, y = lat, group = group), 
               color = "black", fill = "transparent", alpha = 0) +
  geom_polygon(data = states_p, aes(x = long, y = lat, group = group), 
               color = "black", fill = "transparent", alpha = 0) +
  # geom_point(data = all_p, aes(x = long, y = lat), size = 1, color = "black", shape = 21) +
  # geom_path(data = range_p, aes(x = long, y = lat), size = 0.5, linetype = "dashed") +
  # geom_text(data = locs_p, aes(x = long, y = lat, label = fig_name), size = 4) 
  geom_polygon(data = hull, aes(x = long, y = lat), color = "black", fill = "transparent", linetype = "dashed") +
  # geom_label_repel(data = locs_p, aes(x = long, y = lat, label = fig_name), size = 4, segment.size  = 0.5, box.padding = 0.15, point.padding = 0.15,
  # segment.color = "black", force = 0.5) +
  geom_point(data = locs_p, aes(x = long, y = lat), size = 2, color = "magenta3") +
  theme_bw() +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        
        panel.border = element_rect(colour = "black", fill=NA),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent")) 
# ggsave(filename = "figs_for_paper/map_nolabels.pdf", width = 6, height = 6)


# plot temp and ppt vs. latitude -----------------------------------------

# calculate averages from prism climate data

clim_temp = all %>% 
  filter(month != "Aug", tave > -1000) %>% 
  group_by(id) %>% 
  summarize(tave_5180_sep_jul = mean(tave))

clim_ppt = all %>% 
  filter(month %in% c("Apr", "May", "Jun", "Jul") , tave > -1000) %>% 
  group_by(id) %>% 
  summarize(ppt_mm_5180_apr_jul = mean(ppt_mm), lat = mean(lat), long = mean(long), elev = mean(elev))

clim_all = left_join(clim_temp, clim_ppt)

hist(clim_all$tsd_5180_sep_jul)
hist(clim_all$pcov_5180_apr_jul)

plot(clim_all$tsd_5180_sep_jul, clim_all$tave_5180_sep_jul)
plot(clim_all$pcov_5180_apr_jul, clim_all$ppt_mm_5180_apr_jul)

# linear models of temp, ppt, elev, lat, long
# precip is structured by longitude and elevation

ggplot(clim_all, aes(x = long, y = ppt_mm_5180_apr_jul)) +
  geom_point() +
  geom_smooth(method = "lm")
summary(lm(clim_all$ppt_mm_5180_apr_jul~clim_all$long))
# precip increases towards the east

ggplot(clim_all, aes(x = lat, y = ppt_mm_5180_apr_jul)) +
  geom_point() +
  geom_smooth(method = "lm")
summary(lm(clim_all$ppt_mm_5180_apr_jul~clim_all$lat))
# precip increases towards the north

ggplot(clim_all, aes(x = elev, y = ppt_mm_5180_apr_jul)) +
  geom_point() +
  geom_smooth(method = "lm")
summary(lm(clim_all$ppt_mm_5180_apr_jul~clim_all$elev))
# precip increases with elevation

# temperature is primarily structured by elevation

ggplot(clim_all, aes(x = long, y = tave_5180_sep_jul)) +
  geom_point() +
  geom_smooth(method = "lm")
summary(lm(clim_all$tave_5180_sep_jul~clim_all$long))
# temperatures decrease towards the east

ggplot(clim_all, aes(x = lat, y = tave_5180_sep_jul)) +
  geom_point() +
  geom_smooth(method = "lm")
summary(lm(clim_all$tave_5180_sep_jul~clim_all$lat))
# temperatures decrease towards the north

ggplot(clim_all, aes(x = elev, y = tave_5180_sep_jul)) +
  geom_point() +
  geom_smooth(method = "lm")
summary(lm(clim_all$tave_5180_sep_jul~clim_all$elev))
# temperatures decrease with elevation

# temperature vs. latitude subplots for main text
sub1 = ggplot() +
  geom_abline(slope = -0.24290, intercept = 17.63817) +
  geom_point(data = clim_all, aes(x = lat, y = tave_5180_sep_jul, color = elev), alpha = 0.8) +
  geom_point(data = locs, aes(x = lat, y = tave_5180_sep_jul, fill = elev), shape = 21, color = "black", size = 3, alpha = 1) +
  scale_color_continuous(low = "darkgoldenrod1", high = "mediumblue", name = "Elevation\n(m)") +
  scale_fill_continuous(low = "darkgoldenrod1", high = "mediumblue") +
  ylab("Mean annual temperature\n(°C, September-July)") +
  # guides(fill = FALSE, color = guide_colorbar(ticks = FALSE)) +
  guides(fill = FALSE, color = FALSE) +
  xlab("Latitude (°)")

# precip vs. latitude
sub2 = ggplot() +
  geom_abline(slope = 1.5057, intercept = -29.8506) +
  geom_point(data = clim_all, aes(x = lat, y = ppt_mm_5180_apr_jul, color = elev), alpha = 0.8) +
  geom_point(data = locs, aes(x = lat, y = ppt_mm_5180_apr_jul, fill = elev), shape = 21, color = "black", size = 3, alpha = 1) +
  scale_color_continuous(low = "darkgoldenrod1", high = "mediumblue", name = "Elevation\n(m)") +
  scale_fill_continuous(low = "darkgoldenrod1", high = "mediumblue") +
  ylab("Spring/summer precipitation\n(mm, April-July)") +
  xlab("Latitude (°)") +
  guides(fill = FALSE, color = guide_colorbar(ticks = FALSE))

sub = plot_grid(sub1, sub2, ncol = 2, labels = c("A", "B"), align = "h", rel_widths = c(0.5, 0.65))
# ggsave(filename = "figs_for_paper/latitude_climate.pdf", width = 9, height = 3.7)

sub3 = plot_grid(sub1, sub2, ncol = 1, align = "hv", axis = "tlrb")

ggdraw() +
  draw_plot(main, 0.02, 0, 0.55, 1) +
  draw_plot(sub3, 0.58, 0, 0.42, 1) +
  draw_plot_label(c("A", "B", "C"), x = c(0, 0.57, 0.57), y = c(1, 1, 0.5))
# ggsave(filename = "figs/map_climate.pdf", width = 12, height = 6.5)

# multipanel for supplement
lat.temp = ggplot() +
  geom_abline(slope = -0.24290, intercept = 17.63817) +
  geom_point(data = clim_all, aes(x = lat, y = tave_5180_sep_jul, color = elev), alpha = 0.8) +
  geom_point(data = locs, aes(x = lat, y = tave_5180_sep_jul, fill = elev), shape = 21, color = "black", size = 3, alpha = 1) +
  scale_color_continuous(low = "darkgoldenrod1", high = "mediumblue", name = "Elevation\n(m)") +
  scale_fill_continuous(low = "darkgoldenrod1", high = "mediumblue") +
  ylab("Mean annual temperature\n(°C, September-July)") +
  # guides(fill = FALSE, color = guide_colorbar(ticks = FALSE)) +
  guides(fill = FALSE, color = FALSE) +
  theme(axis.title.x = element_blank())

long.temp = ggplot() +
  geom_abline(slope = -0.23984, intercept = -21.77322) +
  geom_point(data = clim_all, aes(x = long, y = tave_5180_sep_jul, color = elev), alpha = 0.8) +
  geom_point(data = locs, aes(x = long, y = tave_5180_sep_jul, fill = elev), shape = 21, color = "black", size = 3, alpha = 1) +
  scale_color_continuous(low = "darkgoldenrod1", high = "mediumblue", name = "Elevation\n(m)") +
  scale_fill_continuous(low = "darkgoldenrod1", high = "mediumblue") +
  # ylab("Mean annual temperature\n(°C, September-July)") +
  guides(fill = FALSE, color = FALSE) +
  # guides(fill = FALSE, color = guide_colorbar(ticks = FALSE)) +
  # xlab("Longitude (°)") +
  theme(axis.title = element_blank())

elev.temp = ggplot() +
  geom_abline(slope = -0.0032009, intercept = 9.7381659) +
  geom_point(data = clim_all, aes(x = elev, y = tave_5180_sep_jul, color = elev), alpha = 0.8) +
  geom_point(data = locs, aes(x = elev, y = tave_5180_sep_jul, fill = elev), shape = 21, color = "black", size = 3, alpha = 1) +
  scale_color_continuous(low = "darkgoldenrod1", high = "mediumblue", name = "Elevation\n(m)") +
  scale_fill_continuous(low = "darkgoldenrod1", high = "mediumblue") +
  # ylab("Mean annual temperature\n(°C, September-July)") +
  guides(fill = FALSE, color = FALSE) +
  # guides(fill = FALSE, color = guide_colorbar(ticks = FALSE)) +
  # xlab("Elevation (m)") +
  theme(axis.title = element_blank())

lat.ppt = ggplot() +
  geom_abline(slope = 1.5057, intercept = -29.8506) +
  geom_point(data = clim_all, aes(x = lat, y = ppt_mm_5180_apr_jul, color = elev), alpha = 0.8) +
  geom_point(data = locs, aes(x = lat, y = ppt_mm_5180_apr_jul, fill = elev), shape = 21, color = "black", size = 3, alpha = 1) +
  scale_color_continuous(low = "darkgoldenrod1", high = "mediumblue", name = "Elevation\n(m)") +
  scale_fill_continuous(low = "darkgoldenrod1", high = "mediumblue") +
  ylab("Spring/summer precip.\n(mm, April-July)") +
  guides(fill = FALSE, color = FALSE) +
  # guides(fill = FALSE, color = guide_colorbar(ticks = FALSE)) +
  xlab("Latitude (°)")

long.ppt = ggplot() +
  geom_abline(slope = 4.441, intercept = 561.660) +
  geom_point(data = clim_all, aes(x = long, y = ppt_mm_5180_apr_jul, color = elev), alpha = 0.8) +
  geom_point(data = locs, aes(x = long, y = ppt_mm_5180_apr_jul, fill = elev), shape = 21, color = "black", size = 3, alpha = 1) +
  scale_color_continuous(low = "darkgoldenrod1", high = "mediumblue", name = "Elevation\n(m)") +
  scale_fill_continuous(low = "darkgoldenrod1", high = "mediumblue") +
  # ylab("Spring/summer precipitation\n(mm, April-July)") +
  guides(fill = FALSE, color = FALSE) +
  # guides(fill = FALSE, color = guide_colorbar(ticks = FALSE)) +
  xlab("Longitude (°)") +
  theme(axis.title.y = element_blank())

elev.ppt = ggplot() +
  geom_abline(slope = 0.018151, intercept = 20.876104) +
  geom_point(data = clim_all, aes(x = elev, y = ppt_mm_5180_apr_jul, color = elev), alpha = 0.8) +
  geom_point(data = locs, aes(x = elev, y = ppt_mm_5180_apr_jul, fill = elev), shape = 21, color = "black", size = 3, alpha = 1) +
  scale_color_continuous(low = "darkgoldenrod1", high = "mediumblue", name = "Elevation\n(m)") +
  scale_fill_continuous(low = "darkgoldenrod1", high = "mediumblue") +
  # ylab("Spring\summer precipitation\n(mm, April-July)") +
  guides(fill = FALSE, color = FALSE) +
  # guides(fill = FALSE, color = guide_colorbar(ticks = FALSE)) +
  xlab("Elevation (m)") +
  theme(axis.title.y = element_blank(), legend.box.background = element_rect(color = "black"), legend.background = element_blank())

legend = get_legend(elev.ppt)

p1 = plot_grid(lat.temp, long.temp, elev.temp, lat.ppt, long.ppt, elev.ppt, ncol = 3, align = "hv", labels = c("A", "B", "C", "D", "E", "F"))
ggdraw() +
  draw_plot(p1, 0, 0, 0.95, 1) + 
  draw_plot(legend, 0.92, 0.67, 0.4, 0.4)
# ggsave("figs_for_paper/climate_geo.pdf", width = 12, height = 7)
