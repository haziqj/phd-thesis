# Chapter 6
# Bayesian Variable Selection using I-priors simulation study
library(ggmap)
library(ggrepel)
chapter.no <- "06"

ozone.points <- read.table("data/ozone_map.txt", header = TRUE)
ozone.points <- cbind(ozone.points, labels = c(
  "El Monte, CA",
  "Sandberg, CA",
  "Upland, CA",
  "Vandenberg AFB",
  "LAX airport",
  "Dagget, CA"
))
ozone.box <- make_bbox(lat = ozone.points$lat, lon = ozone.points$lon, f = 0.2)
ozone.box["bottom"] <- ozone.box["bottom"] - 0.2
ozone.box["top"] <- ozone.box["top"] + 0.2
ozone.map <- get_map(location = ozone.box, maptype = "watercolor",
                     crop = TRUE, source = "stamen", zoom = 10)
ggmap(ozone.map) +
  geom_label_repel(data = ozone.points, aes(x = lon, y = lat, label = labels),
                   box.padding = 1, nudge_x = -0.05, col = "grey30",
                   segment.colour = "grey30") +
  geom_point(data = ozone.points, mapping = aes(x = lon, y = lat),
             fill = "grey30", colour = "grey30", size = 3,
             shape = 21) +
  theme_void() -> p

# Ignore warning messages

ggsave("figure/06-ozone_map.png", p, "png", width = 15, height = 6,
       units = "cm", dpi = 300)
move_fig_to_chapter()  # need to crop a little bit
