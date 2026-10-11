# MetricGraph illustration for Section 6.1
# Run from the repository root.

library(sf)
library(MetricGraph)

dir.create("output", showWarnings = FALSE, recursive = TRUE)

miho <- readRDS(file.path("data", "miho_tp.rds"))
river_m <- st_transform(miho$river, 5186)
sites_m <- st_transform(miho$sites, 5186)

# -----------------------------------------------------------------------------
# Construct metric graph and add observations
# -----------------------------------------------------------------------------
miho_graph <- metric_graph$new(lines = river_m, longlat = FALSE,
                               length_unit = "m",
                               vertex_unit = "m")

miho_graph$add_observations(
  data = sites_m,
  data_coords = "spatial",
  clear_obs = TRUE,
  tolerance = 100,
  include_distance_to_graph = TRUE,
  verbose = 1
)

# -----------------------------------------------------------------------------
# Whittle--Matern model and intercept-only comparison
# -----------------------------------------------------------------------------
fit_wm1 <- graph_lme(
  mean_TP ~ 1,
  graph = miho_graph,
  model = "WM1",
  possible_methods = c("L-BFGS-B", "Nelder-Mead", "BFGS"),
  improve_hessian = TRUE
)

fit_lm <- graph_lme(
  mean_TP ~ 1,
  graph = miho_graph,
  model = "lm"
)

model_comparison <- data.frame(
  model = c("Intercept only", "Whittle-Matern (WM1)"),
  logLik = c(as.numeric(logLik(fit_lm)), as.numeric(logLik(fit_wm1))),
  AIC = c(AIC(fit_lm), AIC(fit_wm1)),
  BIC = c(BIC(fit_lm), BIC(fit_wm1))
)

print(summary(fit_wm1))
print(model_comparison)

# -----------------------------------------------------------------------------
# Prediction on network mesh
# IMPORTANT: predict one location at a time to preserve location matching
# -----------------------------------------------------------------------------

miho_graph$build_mesh(h = 250) #time consuming
#miho_graph$build_mesh(h = 1000) #recommended (for the fast computation)

mesh_location <- data.frame(
  edge_number = miho_graph$mesh$VtE[, 1],
  distance_on_edge = miho_graph$mesh$VtE[, 2]
)

# Predict each mesh location separately
pred_list <- lapply(seq_len(nrow(mesh_location)), function(i) {
  
  loc_i <- mesh_location[i, , drop = FALSE]
  
  pr_i <- predict(
    fit_wm1,
    newdata = loc_i,
    normalized = TRUE,
    compute_variances = TRUE
  )
  
  data.frame(
    row_id = i,
    edge_number = loc_i$edge_number,
    distance_on_edge = loc_i$distance_on_edge,
    predicted_TP = as.numeric(pr_i$mean),
    prediction_variance = as.numeric(pr_i$variance)
  )
})

prediction_df <- do.call(rbind, pred_list)

# Coordinates for exactly the same locations
mesh_xy <- miho_graph$coordinates(
  PtE = as.matrix(
    prediction_df[, c("edge_number", "distance_on_edge")]
  ),
  normalized = TRUE
)

prediction_df$x <- mesh_xy[, 1]
prediction_df$y <- mesh_xy[, 2]

# Convert to sf
metricgraph_pred_sf <- st_as_sf(
  prediction_df,
  coords = c("x", "y"),
  crs = 5186
)

metricgraph_pred_sf <- st_transform(
  metricgraph_pred_sf,
  4326
)

observed_tp_sf <- st_transform(miho$sites, 4326)

saveRDS(
  metricgraph_pred_sf,
  file.path("output", "metricgraph_tp_pred.rds")
)

saveRDS(
  observed_tp_sf,
  file.path("output", "observed_tp.rds")
)