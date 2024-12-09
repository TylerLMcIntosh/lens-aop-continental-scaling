# This script conducts analyses on the representativeness of the NEON AOP footprints at a CONUS-wide scale
# Tyler L. McIntosh

# Setup ----

rm(list = ls())

### ### ### ### ### ### ##
### SET PARAMETER HERE ###
cyverse <- TRUE #Set to TRUE for cyverse processing or FALSE for local processing
### ### ### ### ### ### ##

options(scipen = 999)

if(!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}
library(here)

source(here::here("code", "functions.R"))

install_and_load_packages(
  package_list = c(
    "here",
    "terra",
    "sf",
    "tidyverse",
    "curl",
    "mapview",
    "gdalcubes",
    "tictoc",
    "tigris",
    "tmap",
    "stars",
    "furrr"),
  auto_install = "y"
)


if(cyverse) {
  dir_root <- "~/data-store/data/iplant/home/shared/earthlab/macrosystems/lens-aop-continental-scaling"
} else {
  dir_root <- here::here()
}


# Set up necessary data directories
dir_derived <- file.path(dir_root, "data", "derived")
dir_figs <- here::here("figs")
dir_ensure(dir_figs)



# Load context data for overlays
neon_region_polygons <- access_neon_domains_shp()
epa_region_polygons <- access_data_epa_l2_ecoregions_api() |>
  sf::st_transform(sf::st_crs(neon_region_polygons))




