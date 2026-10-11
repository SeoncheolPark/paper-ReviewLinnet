# smnet illustration for Section 6.3
# Run from the repository root.
#
# This analysis requires the compatibility-patched version of smnet used
# in the manuscript. In particular, smnet::network() must accept a `locs`
# argument while retaining the original penalized-smoothing routines.

library(sf)
library(smnet)

dir.create("output", showWarnings = FALSE, recursive = TRUE)

if (!("locs" %in% names(formals(smnet::network)))) {
  stop(
    "This script requires the compatibility-patched smnet version used in the manuscript."
  )
}

# Build/import the stream-network object if this script is run independently.
if (!exists("mf04p")) {
  source(file.path("scripts", "01_prepare_ssn.R"))
}

# -----------------------------------------------------------------------------
# Prepare adjacency, site locations, and Shreve-based network weights
# -----------------------------------------------------------------------------
edge_df <- st_drop_geometry(mf04p$edges)
obs_df <- st_drop_geometry(mf04p$obs)
obs_df <- obs_df[is.finite(obs_df$mean_TP), ]

adjacency <- smnet::get_adjacency(ssn_path, 1)
K <- nrow(adjacency$adjacency)
all_rid <- as.numeric(adjacency$rid_bid[, 1])

obs_df$segment_index <- smnet:::re_map_rid(
  as.numeric(as.character(obs_df$rid)),
  all_rid
)

edge_df$segment_index <- smnet:::re_map_rid(
  as.numeric(as.character(edge_df$rid)),
  all_rid
)
edge_df <- edge_df[order(edge_df$segment_index), ]

shreve_additive <- smnet:::adjacency_to_shreve(adjacency = adjacency)
weight_vec <- smnet:::get_shreve_weights(
  adjacency = adjacency$adjacency,
  shreve.order = as.matrix(shreve_additive)
)
weight_vec <- as.numeric(weight_vec)

stopifnot(
  length(obs_df$segment_index) == nrow(obs_df),
  length(weight_vec) == K,
  !anyNA(obs_df$segment_index),
  !anyNA(weight_vec),
  all(is.finite(weight_vec)),
  all(is.finite(obs_df$mean_TP))
)

# -----------------------------------------------------------------------------
# Penalized network smoother
# -----------------------------------------------------------------------------
fit_TP <- smnet::smnet(
  formula = mean_TP ~ network(
    adjacency = adjacency,
    weight = weight_vec,
    locs = obs_df$segment_index
  ),
  data.object = obs_df,
  method = "AICC"
)

print(summary(fit_TP))
cat("R2:", fit_TP$R2, "\n")
cat("Adjusted R2:", fit_TP$R2.adj, "\n")

# -----------------------------------------------------------------------------
# Fitted TP on all network segments
# -----------------------------------------------------------------------------
beta <- as.numeric(fit_TP$internals$beta_hat)
stopifnot(length(beta) == K + 1)

intercept <- beta[1]
network_effect <- beta[2:(K + 1)]
segment_fitted_TP <- intercept + network_effect

segment_pred <- data.frame(
  rid = as.numeric(adjacency$rid_bid[, 1]),
  fitted_TP = segment_fitted_TP
)

edges_pred <- mf04p$edges
idx <- match(as.character(edges_pred$rid), as.character(segment_pred$rid))
edges_pred$fitted_TP <- segment_pred$fitted_TP[idx]
stopifnot(!anyNA(edges_pred$fitted_TP))

smnet_pred_sf <- st_transform(edges_pred, 4326)

saveRDS(
  smnet_pred_sf,
  file.path("output", "smnet_tp_pred.rds")
)
