# Reproducibility code for Section 6 of "Statistics for Spatial Linear Networks"
# Run from the repository root.

library(sf)
library(SSNbler)
library(SSN2)

# -----------------------------------------------------------------------------
# Paths and data
# -----------------------------------------------------------------------------
#setwd("~/Library/CloudStorage/Dropbox/R files/ReviewLinnet/github")
dir.create("work", showWarnings = FALSE, recursive = TRUE)

data_file <- file.path("data", "miho_tp.rds")
miho <- readRDS(data_file)

stopifnot(all(c("river", "sites") %in% names(miho)))
stopifnot("mean_TP" %in% names(miho$sites))

# Project to a metric CRS for network-distance calculations.
river_m <- st_transform(miho$river, 5186)
sites_m <- st_transform(miho$sites, 5186)
river_m$Shape_Leng <- as.numeric(st_length(river_m))

# -----------------------------------------------------------------------------
# Build local spatial network (LSN) and .ssn object
# -----------------------------------------------------------------------------
lsn_path <- file.path("work", "Miho_LSN.ssn")
ssn_path <- file.path("work", "Miho_Final.ssn")

net <- lines_to_lsn(
  streams = river_m,
  lsn_path = lsn_path,
  check_topology = TRUE,
  overwrite = TRUE
)

obs_net <- sites_to_lsn(
  sites = sites_m,
  edges = net,
  lsn_path = lsn_path,
  file_name = "obs",
  snap_tolerance = 200,
  save_local = TRUE,
  overwrite = TRUE
)

net <- updist_edges(
  edges = net,
  lsn_path = lsn_path,
  length_col = "Shape_Leng"
)

site_list <- list(obs = obs_net)
site_list <- updist_sites(
  sites = site_list,
  edges = net,
  lsn_path = lsn_path,
  length_col = "Shape_Leng"
)

net <- afv_edges(
  edges = net,
  lsn_path = lsn_path,
  infl_col = "Shape_Leng",
  segpi_col = "areaPI",
  afv_col = "afvArea",
  save_local = TRUE
)

site_list <- afv_sites(
  sites = site_list,
  edges = net,
  lsn_path = lsn_path,
  afv_col = "afvArea",
  save_local = TRUE
)

obs_net <- site_list$obs

mf04p <- ssn_assemble(
  edges = net,
  obs_sites = obs_net,
  lsn_path = lsn_path,
  ssn_path = ssn_path,
  import = TRUE,
  overwrite = TRUE
)

cat("Prepared SSN object with", nrow(mf04p$obs), "monitoring sites.\n")
