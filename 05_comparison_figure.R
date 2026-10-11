# Four-panel comparison figure for Section 6.4 / Figure 3
# Run after scripts 02--04 from the repository root.

library(sf)
library(ggplot2)
library(patchwork)
library(viridis)
library(lwgeom)

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Read data and fitted values
# -----------------------------------------------------------------------------
miho <- readRDS(file.path("data", "miho_tp.rds"))
obs_sf <- readRDS(file.path("output", "observed_tp.rds"))
metric_sf <- readRDS(file.path("output", "metricgraph_tp_pred.rds"))
ssn_sf <- readRDS(file.path("output", "ssn2_tp_pred.rds"))
smnet_sf <- readRDS(file.path("output", "smnet_tp_pred.rds"))

river_ll <- st_transform(miho$river, 4326)
river_m <- st_transform(miho$river, 5186)
metric_m <- st_transform(metric_sf, 5186)
ssn_m <- st_transform(ssn_sf, 5186)
smnet_m <- st_transform(smnet_sf, 5186)

# -----------------------------------------------------------------------------
# Common display network: approximately 500-m segments
# -----------------------------------------------------------------------------
river_lines <- st_cast(st_geometry(river_m), "LINESTRING", warn = FALSE)
river_lines <- river_lines[!st_is_empty(river_lines)]

make_short_segments <- function(line, h = 500) {
  L <- as.numeric(st_length(line))
  n_seg <- max(1, ceiling(L / h))
  out <- lapply(seq_len(n_seg), function(j) {
    lwgeom::st_linesubstring(
      line,
      from = (j - 1) / n_seg,
      to = j / n_seg
    )
  })
  do.call(c, out)
}

display_geom <- do.call(
  c,
  lapply(seq_along(river_lines), function(i) {
    make_short_segments(river_lines[i], h = 500)
  })
)

display_net <- st_sf(
  display_id = seq_along(display_geom),
  geometry = display_geom,
  crs = 5186
)

display_mid_geom <- st_line_sample(st_geometry(display_net), sample = 0.5)
display_mid_geom <- st_cast(display_mid_geom, "POINT")
display_mid <- st_sf(
  display_id = display_net$display_id,
  geometry = display_mid_geom,
  crs = 5186
)

# Map each fitted result to the common display segments.
idx_metric <- st_nearest_feature(display_mid, metric_m)
idx_ssn <- st_nearest_feature(display_mid, ssn_m)
idx_smnet <- st_nearest_feature(display_mid, smnet_m)

display_net$MetricGraph_TP <- metric_m$predicted_TP[idx_metric]
display_net$SSN2_TP <- ssn_m$predicted_TP[idx_ssn]
display_net$smnet_TP <- smnet_m$fitted_TP[idx_smnet]

stopifnot(
  all(is.finite(display_net$MetricGraph_TP)),
  all(is.finite(display_net$SSN2_TP)),
  all(is.finite(display_net$smnet_TP))
)

model_ranges <- rbind(
  Observed = range(obs_sf$mean_TP, na.rm = TRUE),
  MetricGraph = range(display_net$MetricGraph_TP, na.rm = TRUE),
  SSN2 = range(display_net$SSN2_TP, na.rm = TRUE),
  smnet = range(display_net$smnet_TP, na.rm = TRUE)
)
colnames(model_ranges) <- c("Min", "Max")
print(model_ranges)

display_ll <- st_transform(display_net, 4326)

# -----------------------------------------------------------------------------
# Common map extent and color scale
# -----------------------------------------------------------------------------
bb_river <- st_bbox(river_ll)
bb_obs <- st_bbox(obs_sf)
bb <- c(
  xmin = min(bb_river["xmin"], bb_obs["xmin"]),
  ymin = min(bb_river["ymin"], bb_obs["ymin"]),
  xmax = max(bb_river["xmax"], bb_obs["xmax"]),
  ymax = max(bb_river["ymax"], bb_obs["ymax"])
)

x_pad <- 0.015 * (bb["xmax"] - bb["xmin"])
y_pad <- 0.015 * (bb["ymax"] - bb["ymin"])
bb["xmin"] <- bb["xmin"] - x_pad
bb["xmax"] <- bb["xmax"] + x_pad
bb["ymin"] <- bb["ymin"] - y_pad
bb["ymax"] <- bb["ymax"] + y_pad

tp_lim <- range(
  c(
    obs_sf$mean_TP,
    display_net$MetricGraph_TP,
    display_net$SSN2_TP,
    display_net$smnet_TP
  ),
  na.rm = TRUE
)

tp_scale <- scale_color_viridis_c(
  limits = tp_lim,
  option = "viridis",
  name = "TP (mg/L)",
  guide = guide_colorbar(
    direction = "horizontal",
    title.position = "top",
    title.hjust = 0.5,
    barwidth = grid::unit(5.5, "cm"),
    barheight = grid::unit(0.35, "cm")
  )
)

common_coord <- coord_sf(
  xlim = c(bb["xmin"], bb["xmax"]),
  ylim = c(bb["ymin"], bb["ymax"]),
  expand = FALSE,
  datum = NA
)

common_theme <- theme_void(base_size = 10) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 10.5,
      margin = margin(b = 4)
    ),
    plot.margin = margin(4, 4, 4, 4)
  )

# -----------------------------------------------------------------------------
# Four panels
# -----------------------------------------------------------------------------
p_obs <- ggplot() +
  geom_sf(data = river_ll, color = "grey72", linewidth = 0.55, lineend = "round") +
  geom_sf(data = obs_sf, aes(color = mean_TP), size = 2.5) +
  tp_scale + common_coord + labs(title = "(a) Observed mean TP") + common_theme

p_metric <- ggplot() +
  geom_sf(data = display_ll, aes(color = MetricGraph_TP), linewidth = 1.05, lineend = "round") +
  tp_scale + common_coord + labs(title = "(b) MetricGraph") + common_theme

p_ssn <- ggplot() +
  geom_sf(data = display_ll, aes(color = SSN2_TP), linewidth = 1.05, lineend = "round") +
  tp_scale + common_coord + labs(title = "(c) SSN2") + common_theme

p_smnet <- ggplot() +
  geom_sf(data = display_ll, aes(color = smnet_TP), linewidth = 1.05, lineend = "round") +
  tp_scale + common_coord + labs(title = "(d) smnet") + common_theme

p_compare <-
  (p_obs | p_metric) /
  (p_ssn | p_smnet) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(p_compare)

ggsave(
  file.path("figures", "miho_four_panel_comparison.pdf"),
  p_compare,
  width = 7.2,
  height = 6.3,
  units = "in"
)

ggsave(
  file.path("figures", "miho_four_panel_comparison.png"),
  p_compare,
  width = 7.2,
  height = 6.3,
  units = "in",
  dpi = 300
)
