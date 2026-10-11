# SSN2 illustration for Section 6.2
# Run from the repository root.

library(sf)
library(SSNbler)
library(SSN2)

dir.create("output", showWarnings = FALSE, recursive = TRUE)
dir.create("work", showWarnings = FALSE, recursive = TRUE)

# Build/import the stream-network object if this script is run independently.
if (!exists("mf04p")) {
  source(file.path("scripts", "01_prepare_ssn.R"))
}

# -----------------------------------------------------------------------------
# Hydrologic distance matrices and candidate covariance specifications
# -----------------------------------------------------------------------------
ssn_create_distmat(
  ssn.object = mf04p,
  overwrite = TRUE
)

fit_tu <- ssn_lm(
  mean_TP ~ 1,
  ssn.object = mf04p,
  tailup_type = "exponential",
  additive = "afvArea"
)

fit_td <- ssn_lm(
  mean_TP ~ 1,
  ssn.object = mf04p,
  taildown_type = "spherical"
)

fit_tutd <- ssn_lm(
  mean_TP ~ 1,
  ssn.object = mf04p,
  tailup_type = "exponential",
  taildown_type = "spherical",
  additive = "afvArea"
)

fit_full <- ssn_lm(
  mean_TP ~ 1,
  ssn.object = mf04p,
  tailup_type = "exponential",
  taildown_type = "spherical",
  euclid_type = "gaussian",
  additive = "afvArea"
)

ssn_comparison <- data.frame(
  model = c(
    "Tail-up exponential",
    "Tail-down spherical",
    "Tail-up exponential + tail-down spherical",
    "Tail-up + tail-down + Euclidean Gaussian"
  ),
  AIC = c(AIC(fit_tu), AIC(fit_td), AIC(fit_tutd), AIC(fit_full)),
  BIC = c(BIC(fit_tu), BIC(fit_td), BIC(fit_tutd), BIC(fit_full))
)

print(ssn_comparison)

# -----------------------------------------------------------------------------
# Prediction sites: at least one site per segment, approximately 1 km apart
# -----------------------------------------------------------------------------
net_lines <- st_cast(st_geometry(net), "LINESTRING", warn = FALSE)
net_lines <- net_lines[!st_is_empty(net_lines)]
line_length <- as.numeric(st_length(net_lines))
n_each <- pmax(1, ceiling(line_length / 1000))

pred_list_geom <- lapply(seq_along(net_lines), function(i) {
  x <- st_line_sample(net_lines[i], n = n_each[i])
  if (st_is_empty(x)) return(NULL)
  st_cast(x, "POINT")
})

pred_list_geom <- pred_list_geom[!vapply(pred_list_geom, is.null, logical(1))]
pred_geom <- do.call(c, pred_list_geom)

pred_1km <- st_sf(
  pred_id = seq_along(pred_geom),
  geometry = pred_geom,
  crs = st_crs(net)
)

pred_net <- sites_to_lsn(
  sites = pred_1km,
  edges = net,
  lsn_path = lsn_path,
  file_name = "pred1km",
  snap_tolerance = 5,
  save_local = TRUE,
  overwrite = TRUE
)

pred_list <- list(pred1km = pred_net)
pred_list <- updist_sites(
  sites = pred_list,
  edges = net,
  lsn_path = lsn_path,
  length_col = "Shape_Leng"
)

pred_list <- afv_sites(
  sites = pred_list,
  edges = net,
  lsn_path = lsn_path,
  afv_col = "afvArea",
  save_local = TRUE
)

pred_net <- pred_list$pred1km

mf04p_pred <- ssn_assemble(
  edges = net,
  obs_sites = obs_net,
  preds_list = list(pred1km = pred_net),
  lsn_path = lsn_path,
  ssn_path = file.path("work", "Miho_Final_Pred.ssn"),
  import = TRUE,
  overwrite = TRUE
)

ssn_create_distmat(
  ssn.object = mf04p_pred,
  predpts = "pred1km",
  overwrite = TRUE
)

# Table 2 selects the spherical tail-down specification.
fit_td_pred <- ssn_lm(
  mean_TP ~ 1,
  ssn.object = mf04p_pred,
  taildown_type = "spherical"
)

ssn2_pred_sf <- augment(
  fit_td_pred,
  newdata = "pred1km"
)

ssn2_pred_sf$predicted_TP <- ssn2_pred_sf$.fitted
ssn2_pred_sf <- st_transform(ssn2_pred_sf, 4326)

saveRDS(
  ssn2_pred_sf,
  file.path("output", "ssn2_tp_pred.rds")
)
