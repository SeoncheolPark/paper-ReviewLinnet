# Miho River reproducibility code

This directory contains the code used for the real-data illustration in Section 6 of **Statistics for Spatial Linear Networks**.

## Data

`data/miho_tp.rds` should contain a list with two objects:

- `river`: Miho River network as an `sf` object in WGS84 (EPSG:4326)
- `sites`: 28 monitoring sites as an `sf` object with columns `name` and `mean_TP`

Only the 2023--2025 site-averaged total phosphorus (TP) values used in the manuscript are retained. 

## Run order

Run from the repository root:

1. `01_prepare_ssn.R`
2. `02_metricgraph.R`
3. `03_ssn2.R`
4. `04_smnet.R`
5. `05_comparison_figure.R`

Scripts 03 and 04 source script 01 automatically when needed.

## Important note on smnet

The manuscript analysis uses a compatibility-patched version of `smnet` that allows the processed stream-network location vector to be supplied through `network(..., locs = ...)`. The original penalized-smoothing routines are otherwise unchanged. For full reproducibility, the patched `smnet` source (or an explicit patch) is also included in the repository (`smnet_3.0.tar.gz`).

## Outputs

Intermediate fitted values are written to `output/`, generated SSN files to `work/`, and the final four-panel comparison figure to `figures/`.
