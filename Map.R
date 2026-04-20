# ============================================================
# Spatial application of the edge-ratio framework in Sabah
#
# Panel A: Sabah location inset + forest / non-forest map
# Panel B: Local edge-ratio map
# Panel C: Continuous model-derived relative spillover risk
#
# Notes:
# - Panel A uses binary forest / non-forest in Sabah
# - Edge detection is computed on a buffered regional raster,
#   then masked back to Sabah
# - Panel C applies a continuous hump-shaped risk function to
#   the local edge-ratio raster
# ============================================================

# ---- 0. Packages ----
library(terra)
library(geodata)
library(sf)
library(ggplot2)
library(viridis)
library(gridExtra)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)

# ---- 1. Settings ----
dir.create("geo_data", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

forest_threshold <- 0.50
agg_factor <- 1          # set to 2 if slow
window_n <- 21
buffer_width <- 50000    # 50 km buffer

# Continuous risk function settings
# Use the model peak around 0.45 and spread around 0.10
risk_peak <- 0.45
risk_sd   <- 0.10

# ---- 2. Download Sabah boundary ----
malaysia_adm1 <- gadm(country = "Malaysia", level = 1, path = "geo_data")

if (!"Sabah" %in% unique(malaysia_adm1$NAME_1)) {
  stop("Sabah not found in Malaysia level-1 boundaries.")
}

sabah <- malaysia_adm1[malaysia_adm1$NAME_1 == "Sabah", ]
writeVector(sabah, "geo_data/sabah_boundary.gpkg", overwrite = TRUE)

# ---- 3. Download tree-cover raster ----
trees <- landcover(var = "trees", path = "geo_data")
trees <- trees[[1]]  # ensure single layer

# Match CRS
sabah <- project(sabah, crs(trees))
sabah_sf <- st_as_sf(sabah)

# ---- 4. Crop and mask to Sabah ----
trees_sabah <- crop(trees, sabah)
trees_sabah <- mask(trees_sabah, sabah)
writeRaster(trees_sabah, "geo_data/trees_sabah.tif", overwrite = TRUE)

# ---- 5. Convert tree cover to binary forest / non-forest ----
forest_bin <- ifel(trees_sabah >= forest_threshold, 1, 0)

if (agg_factor > 1) {
  forest_bin <- aggregate(forest_bin, fact = agg_factor, fun = modal)
}

writeRaster(forest_bin, "geo_data/forest_bin_sabah.tif", overwrite = TRUE)

forest_df <- as.data.frame(forest_bin, xy = TRUE, na.rm = TRUE)
names(forest_df)[3] <- "forest_class"

# ---- 6. Better edge detection using a buffered regional raster ----
sabah_buffer <- buffer(sabah, width = buffer_width)

# Crop to buffered region, do not mask yet
trees_region <- crop(trees, sabah_buffer)

# Forest / non-forest in the larger region
forest_region <- ifel(trees_region >= forest_threshold, 1, 0)

# Treat no-data / ocean as non-forest
forest_region[is.na(forest_region)] <- 0

if (agg_factor > 1) {
  forest_region <- aggregate(forest_region, fact = agg_factor, fun = modal)
}

# Count forest neighbors
forest_neighbors <- focal(
  forest_region,
  w = matrix(1, 3, 3),
  fun = sum,
  na.rm = FALSE
)

# Force matching geometry if needed
if (!compareGeom(forest_region, forest_neighbors, stopOnError = FALSE)) {
  forest_neighbors <- resample(forest_neighbors, forest_region, method = "near")
}

# Forest edge cells
edge_region <- ifel((forest_region == 1) & (forest_neighbors < 9), 1, 0)

# Mask back to Sabah
edge_cells <- mask(edge_region, sabah)
writeRaster(edge_cells, "geo_data/edge_cells_sabah.tif", overwrite = TRUE)

# ---- 7. Local edge ratio ----
w <- matrix(1, window_n, window_n)

edge_ratio_local <- focal(
  edge_cells,
  w = w,
  fun = mean,
  na.rm = TRUE,
  expand = FALSE
)

# Strictly limit to Sabah
edge_ratio_local <- mask(edge_ratio_local, sabah)
writeRaster(edge_ratio_local, "geo_data/edge_ratio_local_sabah.tif", overwrite = TRUE)

edge_df <- as.data.frame(edge_ratio_local, xy = TRUE, na.rm = TRUE)
names(edge_df)[3] <- "edge_ratio"

# ---- 8. Continuous model-derived relative spillover risk ----
# Gaussian-shaped relative risk centered on model peak
risk_surface <- exp(-((edge_ratio_local - risk_peak)^2) / (2 * risk_sd^2))

# Normalize to max = 1
risk_surface <- risk_surface / global(risk_surface, "max", na.rm = TRUE)[1,1]

# Strictly limit to Sabah
risk_surface <- mask(risk_surface, sabah)
writeRaster(risk_surface, "geo_data/risk_surface_sabah.tif", overwrite = TRUE)

risk_df <- as.data.frame(risk_surface, xy = TRUE, na.rm = TRUE)
names(risk_df)[3] <- "risk"

# ---- 9. Larger inset map ----
world <- ne_countries(scale = "medium", returnclass = "sf")

inset_bbox <- st_bbox(c(
  xmin = 90, xmax = 135,
  ymin = -12, ymax = 25
), crs = st_crs(world))

inset_region <- suppressWarnings(st_crop(world, inset_bbox))
sabah_location_sf <- st_transform(sabah_sf, st_crs(world))

p_inset <- ggplot() +
  geom_sf(data = inset_region, fill = "grey90", color = "grey40", linewidth = 0.2) +
  geom_sf(data = sabah_location_sf, fill = "red", color = "grey40", linewidth = 0.3) +
  coord_sf(
    xlim = c(90, 135),
    ylim = c(-12, 25),
    expand = FALSE
  ) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.3)
  )

# ---- 10. Panel A: forest / non-forest map with inset ----
p_context_main <- ggplot() +
  geom_raster(data = forest_df, aes(x = x, y = y, fill = factor(forest_class))) +
  geom_sf(data = sabah_sf, fill = NA, color = "black", linewidth = 0.2) +
  coord_sf(expand = FALSE) +
  scale_fill_manual(
    values = c("0" = "#d9d9d9", "1" = "#238b45"),
    breaks = c("1", "0"),
    labels = c("1" = "Forest", "0" = "Non-forest"),
    name = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )

p_context <- p_context_main +
  inset_element(
    p_inset,
    left = 0.63, bottom = 0.58,
    right = 1, top = 1,
    align_to = "panel"
  )

# ---- 11. Panel B: local edge-ratio map ----
p_edge <- ggplot() +
  geom_raster(data = edge_df, aes(x = x, y = y, fill = edge_ratio)) +
  geom_sf(data = sabah_sf, fill = NA, color = "black", linewidth = 0.2) +
  coord_sf(expand = FALSE) +
  scale_fill_viridis(
    name = expression("Local edge ratio, " * E(x)),
    option = "cividis",
    direction = -1
  ) +
  labs(title = "(a)") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    legend.position = 'bottom'
  )

# ---- 12. Panel C: continuous relative risk map ----
p_risk <- ggplot() +
  geom_raster(data = risk_df, aes(x = x, y = y, fill = risk)) +
  geom_sf(data = sabah_sf, fill = NA, color = "black", linewidth = 0.2) +
  coord_sf(expand = FALSE) +
  scale_fill_viridis(
    name = expression("Relative spillover risk, " * R(x)),
    option = "magma",
    limits = c(0, 1),
    direction = -1
  ) +
  labs(title = "(b)") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    legend.position = 'bottom'
  )

# ---- 13. Combine panels ----
Fig_6 <- p_context
Fig_7 <- grid.arrange(p_edge, p_risk, nrow = 1) 

# ---- 14. Save outputs ----
ggsave("figures/Fig.6.png", Fig_6, width = 8, height = 6, dpi = 900)
ggsave("figures/Fig.7.png", Fig_7, width = 8, height = 6, dpi = 900)
