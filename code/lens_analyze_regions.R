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

## Cyverse data store access if applicable ----
# Copy previous data over from data store if on cyverse and project has been run before
if(cyverse) {
  if(dir.exists("~/data-store/data/iplant/home/shared/earthlab/macrosystems/lens-aop-continental-scaling/data")) {
    system("cp -r ~/data-store/data/iplant/home/shared/earthlab/macrosystems/lens-aop-continental-scaling/data ~/lens-aop-continental-scaling/")
  }
  if(dir.exists("~/data-store/data/iplant/home/shared/earthlab/macrosystems/lens-aop-continental-scaling/figs")) {
    system("cp -r ~/data-store/data/iplant/home/shared/earthlab/macrosystems/lens-aop-continental-scaling/figs ~/lens-aop-continental-scaling/")
  }
}


# Set up necessary data directories
dir_ensure(here::here("data"))
dir_raw <- here::here("data/raw")
dir_derived <- here::here("data/derived")
dir_figs <- here::here("figs")
dir_ensure(dir_raw)
dir_ensure(dir_derived)
dir_ensure(dir_figs)

#For outputs

options(timeout = 1500)


## Download and stream data sources ----
tic()
raster <- access_landfire_evt_conus_2022(access = 'download',
                                         dir_path = dir_raw)
evt_conus_file <- here::here(dir_raw, "LF2022_EVT_230_CONUS/Tif/LC22_EVT_230.tif")
toc()

raster_cats <- access_landfire_evt_conus_2022_csv() |>
  dplyr::mutate(VALUE = as.integer(VALUE))
neon_region_polygons <- access_neon_domains_shp()
epa_region_polygons <- access_data_epa_l3_ecoregions_vsi() |>
  sf::st_transform(sf::st_crs(neon_region_polygons))
areas_of_interest <- access_neon_aop_flight_box_data() #note that flight boxes have the domain data as "D##" instead of just "##"


areas_of_interest <- areas_of_interest |>
  sf::st_transform(sf::st_crs(neon_region_polygons))


# Operate ----

## Prep region & AOI sets ----

# Get CONUS bounds to spatial subset
conus <- tigris::states() |>
  dplyr::filter(!STUSPS %in% c("AK", "HI", "PR", "VI", "MP", "GU", "AS")) |>
  dplyr::summarise(geometry = sf::st_union(geometry)) |>
  sf::st_transform(sf::st_crs(neon_region_polygons))

# Merge NEON regions and aois for analyses
neon_region_polygons_merged <- neon_region_polygons |>
  sf::st_cast("POLYGON") %>% #split out multipolygons
  dplyr::filter(sf::st_intersects(., conus, sparse = FALSE)) |>
  dplyr::group_by(DomainID, DomainName) |>
  dplyr::summarise(geometry = sf::st_union(geometry)) |>
  dplyr::ungroup() |>
  dplyr::arrange(DomainID)

areas_of_interest_merged <- areas_of_interest |>
  sf::st_make_valid() %>% #force validity, duplicate vertex error
  dplyr::filter(sf::st_intersects(., conus, sparse = FALSE)) |>
  dplyr::group_by(domain, domainName) |>
  dplyr::summarise(geometry = sf::st_union(geometry)) |>
  dplyr::ungroup() |>
  dplyr::rename(DomainName = domainName) |>
  dplyr::mutate(DomainID = as.integer(substr_right(domain, 2))) |>
  dplyr::arrange(DomainID)

#Ensure that the same sets are in both
common_domain_ids <- lubridate::intersect(neon_region_polygons_merged$DomainID, areas_of_interest_merged$DomainID)
neon_region_polygons_merged <- neon_region_polygons_merged %>%
  dplyr::filter(DomainID %in% common_domain_ids)
areas_of_interest_merged <- areas_of_interest_merged %>%
  dplyr::filter(DomainID %in% common_domain_ids)

neon_region_polygons_merged_file <- here::here(dir_derived, "neon_region_polygons_merged.gpkg")
sf::st_write(neon_region_polygons_merged, neon_region_polygons_merged_file)

areas_of_interest_merged_file <- here::here(dir_derived, "areas_of_interest_merged.gpkg")
sf::st_write(areas_of_interest_merged, areas_of_interest_merged_file)

#Set up to use EPA ecoregions as well
# epa_region_polygons_merged <- epa_region_polygons |>
#   sf::st_cast("POLYGON") %>% #split out multipolygons
#   dplyr::filter(sf::st_intersects(., conus, sparse = FALSE)) |>
#   # dplyr::group_by(DomainID, DomainName) |>
#   # dplyr::summarise(geometry = sf::st_union(geometry)) |>
#   # dplyr::ungroup() |>
#   # dplyr::arrange(DomainID)
#   
# epa_region_polygons_merged_file <- here::here(dir_derived, "epa_region_polygons_merged.gpkg")
# sf::st_write(epa_region_polygons_merged, epa_region_polygons_merged_file)


## Prep to run analyses ----
aoi_thresholds <- c(0.001, 0.01, 0.1, 1)

ag_dev_mine_evt_names <- raster_cats |> dplyr::filter(
  grepl(pattern = "Developed", x = EVT_PHYS, ignore.case = FALSE) |
    grepl(pattern = "Agricultural", x = EVT_PHYS, ignore.case = FALSE) |
    grepl(pattern = "Quarries", x = EVT_PHYS, ignore.case = FALSE)) |>
  dplyr::pull(EVT_NAME) |>
  as.vector()

#For conus analyses

# conus_neon <- neon_region_polygons |>
#   sf::st_cast("POLYGON") %>% #split out multipolygons
#   dplyr::filter(sf::st_intersects(., conus, sparse = FALSE)) |>
#   dplyr::summarise(geometry = sf::st_union(geometry))
#   
# conus_aoi <- areas_of_interest_merged |>
#   dplyr::summarise(geometry = sf::st_union(geometry))
# 
# 
# x <- conus_neon |> 
#   sf::st_transform(5070) |>
#   sf::st_buffer(2) |>
#   sf::st_transform(sf::st_crs(conus_neon))
# mapview(x)


## NEON domains analyses ----

# Run analysis with AOI threshold of 0.001% - 1% and all EVT classes
future::plan(multisession, workers = 4)
tic()
furrr::future_map2(.x = aoi_thresholds,
             .y = c("neon_domains_evt_raw_all_0001",
                    "neon_domains_evt_raw_all_001",
                    "neon_domains_evt_raw_all_01",
                    "neon_domains_evt_raw_all_1"),
             .f = function(x, y) conus_lens_analysis(region_polygons_merged = neon_region_polygons_merged_file,
                                                     areas_of_interest_merged = areas_of_interest_merged_file,
                                                     raster = evt_conus_file,
                                                     raster_cat_df = raster_cats,
                                                     run_name = y,
                                                     cat_base_column_name = "VALUE",
                                                     aoi_drop_perc = x,
                                                     drop_classes = NA,
                                                     drop_classes_column_name = NA,
                                                     out_rast_values = "BOTH",
                                                     out_rast_type = "BOTH",
                                                     out_dir = here::here('data/derived/')))
toc()


# Run analysis with AOI threshold of 0.001% - 1% and all EVT classes except ag_dev_mine
tic()
purrr::walk2(.x = aoi_thresholds,
             .y = c("neon_domains_evt_raw_no_ag_dev_mine_0001",
                    "neon_domains_evt_raw_no_ag_dev_mine_001",
                    "neon_domains_evt_raw_no_ag_dev_mine_01",
                    "neon_domains_evt_raw_no_ag_dev_mine_1"),
             .f = function(x, y) conus_lens_analysis(region_polygons_merged = region_polygons_merged,
                                                     areas_of_interest_merged = areas_of_interest_merged,
                                                     raster = raster,
                                                     raster_cat_df = raster_cats,
                                                     run_name = y,
                                                     cat_base_column_name = "VALUE",
                                                     aoi_drop_perc = x,
                                                     drop_classes = ag_dev_mine_evt_names,
                                                     drop_classes_column_name = "EVT_NAME",
                                                     out_rast_values = "BOTH",
                                                     out_rast_type = "BOTH",
                                                     out_dir = here::here('data/derived/')))
toc()


# Run analysis with AOI threshold of 0.001% - 1% and all EVT groups
tic()
purrr::walk2(.x = aoi_thresholds,
             .y = c("neon_domains_evt_groups_all_0001",
                    "neon_domains_evt_groups_all_001",
                    "neon_domains_evt_groups_all_01",
                    "neon_domains_evt_groups_all_1"),
             .f = function(x, y) conus_lens_analysis(region_polygons_merged = region_polygons_merged[1,],
                                                     areas_of_interest_merged = areas_of_interest_merged[1,],
                                                     raster = raster,
                                                     raster_cat_df = raster_cats,
                                                     run_name = y,
                                                     cat_base_column_name = "EVT_GP",
                                                     aoi_drop_perc = x,
                                                     drop_classes = NA,
                                                     drop_classes_column_name = NA,
                                                     out_rast_values = "BOTH",
                                                     out_rast_type = "BOTH",
                                                     out_dir = here::here('data/derived/')))
toc()



# Run analysis with AOI threshold of 0.001% - 1% and EVT groups except ag_dev_mine
tic()
purrr::walk2(.x = aoi_thresholds,
             .y = c("neon_domains_evt_groups_no_ag_dev_mine_0001",
                    "neon_domains_evt_groups_no_ag_dev_mine_001",
                    "neon_domains_evt_groups_no_ag_dev_mine_01",
                    "neon_domains_evt_groups_no_ag_dev_mine_1"),
             .f = function(x, y) conus_lens_analysis(region_polygons_merged = region_polygons_merged,
                                                     areas_of_interest_merged = areas_of_interest_merged,
                                                     raster = raster,
                                                     raster_cat_df = raster_cats,
                                                     run_name = y,
                                                     cat_base_column_name = "EVT_GP",
                                                     aoi_drop_perc = x,
                                                     drop_classes = ag_dev_mine_evt_names,
                                                     drop_classes_column_name = "EVT_NAME",
                                                     out_rast_values = "BOTH",
                                                     out_rast_type = "BOTH",
                                                     out_dir = here::here('data/derived/')))
toc()


## CONUS analyses ----

tic()
purrr::walk2(.x = aoi_thresholds,
             .y = c("conus_evt_raw_all_0001",
                    "conus_evt_raw_all_001",
                    "conus_evt_raw_all_01",
                    "conus_evt_raw_all_1"),
             .f = function(x, y) conus_lens_analysis(region_polygons_merged = region_polygons_merged,
                                                     areas_of_interest_merged = areas_of_interest_merged,
                                                     raster = raster,
                                                     raster_cat_df = raster_cats,
                                                     run_name = y,
                                                     cat_base_column_name = "VALUE",
                                                     aoi_drop_perc = x,
                                                     drop_classes = NA,
                                                     drop_classes_column_name = NA,
                                                     out_rast_values = "BOTH",
                                                     out_rast_type = "BOTH",
                                                     out_dir = here::here('data/derived/')))
toc()

tic()
purrr::walk2(.x = aoi_thresholds,
             .y = c("conus_evt_groups_all_0001",
                    "conus_evt_groups_all_001",
                    "conus_evt_groups_all_01",
                    "conus_evt_groups_all_1"),
             .f = function(x, y) conus_lens_analysis(region_polygons_merged = region_polygons_merged,
                                                     areas_of_interest_merged = areas_of_interest_merged,
                                                     raster = raster,
                                                     raster_cat_df = raster_cats,
                                                     run_name = y,
                                                     cat_base_column_name = "EVT_GP",
                                                     aoi_drop_perc = x,
                                                     drop_classes = NA,
                                                     drop_classes_column_name = NA,
                                                     out_rast_values = "BOTH",
                                                     out_rast_type = "BOTH",
                                                     out_dir = here::here('data/derived/')))
toc()








# Move data to cyverse data store if applicable ----
if(cyverse) {
  system("cp -r ~/lens-aop-continental-scaling/data ~/data-store/data/iplant/home/shared/earthlab/macrosystems/lens-aop-continental-scaling")
  system("cp -r ~/lens-aop-continental-scaling/figs ~/data-store/data/iplant/home/shared/earthlab/macrosystems/lens-aop-continental-scaling")
}




