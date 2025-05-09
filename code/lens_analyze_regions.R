# This script conducts analyses on the representativeness of the NEON AOP footprints at a CONUS-wide scale
# Tyler L. McIntosh

# Setup ----

rm(list = ls())

### ### ### ### ### ### ##
### SET PARAMETER HERE ###
cyverse <- FALSE #Set to TRUE for cyverse processing or FALSE for local processing
mem_tot <- 60 # Set to run on CyVerse since cyverse has problems reading RAM
### ### ### ### ### ### ##

options(scipen = 999)

if(!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}
library(here)

source(here::here("code", "functions.R"))

if(cyverse) {
  options(repos = c(CRAN = "https://packagemanager.posit.co/cran/2025-04-11"))
}
Sys.setenv(PAK_NONINTERACTIVE = "true")

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
    "furrr",
    "pals",
    "biscale",
    "RColorBrewer",
    "ggtext",
    "patchwork",
    "grid",
    "memuse"),
  auto_install = "y"
)


# Set up necessary data directories
dir_ensure(here::here("data"))
dir_raw <- here::here("data/raw")
dir_derived <- here::here("data/derived")
dir_figs <- here::here("figs")
dir_output <- here::here("data/output")

#If on CyVerse, copy data on CyVerse data store onto instance
if(cyverse) {
  if(dir.exists("~/data-store/data/iplant/home/shared/earthlab/macrosystems/lens-aop-continental-scaling/data/derived") &
     !dir.exists(dir_derived)) {
    system("cp -r ~/data-store/data/iplant/home/shared/earthlab/macrosystems/lens-aop-continental-scaling/data/derived ~/lens-aop-continental-scaling/data/")
  } 
  if(dir.exists("~/data-store/data/iplant/home/shared/earthlab/macrosystems/lens-aop-continental-scaling/data/raw") &
     !dir.exists(dir_raw)) {
    system("cp -r ~/data-store/data/iplant/home/shared/earthlab/macrosystems/lens-aop-continental-scaling/data/raw ~/lens-aop-continental-scaling/data/")
  }
}

#Ensure directory creation
dir_ensure(dir_raw)
dir_ensure(dir_derived)
dir_ensure(dir_figs)
dir_ensure(dir_output)


#Manage memory
if(cyverse) {
  terra::terraOptions(memmax = mem_tot)
  terra::terraOptions(memfrac = 0.8)
} else {
  mem_tot <- memuse::Sys.meminfo()[[1]]
}


## Download and stream data sources ----
tic()
options(timeout = 1500)
raster <- access_landfire_evt_conus_2022(access = 'download',
                                         dir_path = dir_raw)
evt_conus_file <- here::here(dir_raw, "LF2022_EVT_230_CONUS/Tif/LC22_EVT_230.tif")
toc()

raster_cats <- access_landfire_evt_conus_2022_csv() |>
  dplyr::mutate(VALUE = as.integer(VALUE))
neon_region_polygons <- access_neon_domains_shp()
epa_region_polygons <- access_data_epa_l2_ecoregions_api() |>
  sf::st_transform(sf::st_crs(neon_region_polygons))
areas_of_interest <- access_neon_aop_flight_box_data() #note that flight boxes have the domain data as "D##" instead of just "##"

# Create a grouped version of the raster
r_grp_fl <- here::here('data/derived/lfr_grouped.tif')
r_grp_cat_fl <- here::here('data/derived/lfr_grouped_cats.csv')
if(!file.exists(r_grp_fl)) {
  r_grp <- reclassify_raster_by_group(
    raster = raster,
    raster_cat_df = raster_cats,
    value_column = "VALUE",
    group_column = "EVT_GP_N",
    write_to_disk = TRUE,
    out_path = r_grp_fl
  )
  readr::write_csv(r_grp$cat_df,
                   r_grp_cat_fl)
}

raster_grouped <- terra::rast(r_grp_fl)
raster_cats_grouped <- readr::read_csv(r_grp_cat_fl) |>
  as.data.frame()


# Get CONUS bounds to spatial subset
conus <- tigris::states() |>
  dplyr::filter(!STUSPS %in% c("AK", "HI", "PR", "VI", "MP", "GU", "AS")) |>
  dplyr::summarise(geometry = sf::st_union(geometry)) |>
  sf::st_transform(sf::st_crs(neon_region_polygons))

areas_of_interest <- areas_of_interest |>
  sf::st_transform(sf::st_crs(neon_region_polygons)) |>
  sf::st_make_valid() %>% #force validity, duplicate vertex error
  dplyr::filter(sf::st_intersects(., conus, sparse = FALSE))


# MANUALLY CREATE TEST POLYGONS IF DESIRED
# install_and_load_packages(c("mapedit"))
# 
# test_rgn <- mapedit::drawFeatures()
# test_rgn <- test_rgn |>
#   sf::st_transform(sf::st_crs(neon_region_polygons))
# 
# test_aoi <- mapedit::drawFeatures()
# test_aoi <- test_aoi |>
#   sf::st_transform(sf::st_crs(neon_region_polygons))


# Operate ----

## Prep region & AOI sets ----

# Merge NEON regions and aois for analyses
neon_region_polygons_merged <- neon_region_polygons |>
  sf::st_cast("POLYGON") %>% #split out multipolygons
  dplyr::filter(sf::st_intersects(., conus, sparse = FALSE)) |>
  dplyr::group_by(DomainID, DomainName) |>
  dplyr::summarise(geometry = sf::st_union(geometry)) |>
  dplyr::ungroup() |>
  dplyr::arrange(DomainID)

neon_areas_of_interest_merged <- areas_of_interest |>
  dplyr::group_by(domain, domainName) |>
  dplyr::summarise(geometry = sf::st_union(geometry)) |>
  dplyr::ungroup() |>
  dplyr::rename(DomainName = domainName) |>
  dplyr::mutate(DomainID = as.integer(substr_right(domain, 2))) |>
  dplyr::arrange(DomainID)

#Ensure that the same sets are in both
neon_common_domain_ids <- lubridate::intersect(neon_region_polygons_merged$DomainID, neon_areas_of_interest_merged$DomainID)
neon_region_polygons_merged <- neon_region_polygons_merged %>%
  dplyr::filter(DomainID %in% neon_common_domain_ids)
neon_areas_of_interest_merged <- neon_areas_of_interest_merged %>%
  dplyr::filter(DomainID %in% neon_common_domain_ids)

neon_region_polygons_merged_file <- here::here(dir_derived, "neon_region_polygons_merged.gpkg")
sf::st_write(neon_region_polygons_merged, neon_region_polygons_merged_file, append = FALSE)

neon_areas_of_interest_merged_file <- here::here(dir_derived, "neon_areas_of_interest_merged.gpkg")
sf::st_write(neon_areas_of_interest_merged, neon_areas_of_interest_merged_file, append = FALSE)



# #Set up to use EPA ecoregions as well
# epa_region_polygons_merged <- epa_region_polygons |>
#   sf::st_cast("POLYGON") |> #split out multipolygons
#   sf::st_make_valid() %>%
#   dplyr::filter(sf::st_intersects(., conus, sparse = FALSE)) |>
#   dplyr::group_by(NA_L2CODE, NA_L2KEY) |>
#   dplyr::summarise(geometry = sf::st_union(geometry)) |>
#   dplyr::ungroup() |>
#   dplyr::arrange(NA_L2KEY) |>
#   dplyr::filter(NA_L2CODE != "0.0") |>
#   sf::st_intersection(conus)
# 
# epa_areas_of_interest_merged <- areas_of_interest |>
#   sf::st_intersection(epa_region_polygons_merged) |>
#   dplyr::group_by(NA_L2CODE, NA_L2KEY) |>
#   dplyr::summarise(geometry = sf::st_union(geometry)) |>
#   dplyr::ungroup() |>
#   dplyr::arrange(NA_L2KEY)
# 
# #Ensure that the same sets are in both
# common_epa_ids <- lubridate::intersect(epa_region_polygons_merged$NA_L2KEY, epa_areas_of_interest_merged$NA_L2KEY)
# epa_region_polygons_merged <- epa_region_polygons_merged %>%
#   dplyr::filter(NA_L2KEY %in% common_epa_ids)
# epa_areas_of_interest_merged <- epa_areas_of_interest_merged %>%
#   dplyr::filter(NA_L2KEY %in% common_epa_ids)
# 
# epa_region_polygons_merged_file <- here::here(dir_derived, "epa_region_polygons_merged.gpkg")
# sf::st_write(epa_region_polygons_merged, epa_region_polygons_merged_file, append = FALSE)
# 
# epa_areas_of_interest_merged_file <- here::here(dir_derived, "epa_areas_of_interest_merged.gpkg")
# sf::st_write(epa_areas_of_interest_merged, epa_areas_of_interest_merged_file, append = FALSE)


## Prep to run analyses ----
# aoi_thresholds <- c(0.001, 0.01, 0.1, 1)
# 
# ag_dev_mine_evt_names <- raster_cats |> dplyr::filter(
#   grepl(pattern = "Developed", x = EVT_PHYS, ignore.case = FALSE) |
#     grepl(pattern = "Agricultural", x = EVT_PHYS, ignore.case = FALSE) |
#     grepl(pattern = "Quarries", x = EVT_PHYS, ignore.case = FALSE)) |>
#   dplyr::pull(EVT_NAME) |>
#   as.vector()

#For conus analyses

conus_epa <- epa_region_polygons |>
  sf::st_cast("POLYGON") |> #split out multipolygons
  dplyr::filter(NA_L2CODE != "0.0") |>
  sf::st_make_valid() %>%
  dplyr::filter(sf::st_intersects(., conus, sparse = FALSE)) |>
  dplyr::summarise(geometry = sf::st_union(geometry)) |>
  sf::st_intersection(conus)

conus_epa_clean <- conus_epa |> 
  sf::st_transform(5070) |>
  sf::st_buffer(2) |> # clean up a few gaps in the polygon
  sf::st_transform(sf::st_crs(conus_epa)) |>
  dplyr::mutate(Region = "CONUS")
mapview(conus_epa_clean)


conus_aoi <- areas_of_interest |>
  sf::st_intersection(conus_epa_clean) |>
  dplyr::mutate(conus = "conus") |>
  dplyr::group_by(conus) |>
  dplyr::summarise(geometry = sf::st_union(geometry)) |>
  dplyr::ungroup()



# RUN ANALYSIS ----

# Full CONUS, grouped and ungrouped

tic()
test_full <- full_representative_categorical_analysis(raster = raster,
                                                      raster_cat_df = raster_cats,
                                                      region_shape = conus_epa_clean,
                                                      aoi_shape = conus_aoi,
                                                      run_name = "CONUS",
                                                      cat_base_column_name = "VALUE",
                                                      out_rast_values = c("PERC_COVER_AOI", "PERC_COVER_REGION"),
                                                      out_rast_type = c("FULL"),
                                                      out_dir = dir_output,
                                                      new_sub_dir = TRUE,
                                                      min_aoi_coverage = NA,
                                                      min_region_coverage = NA,
                                                      drop_classes = c("Open Water"),
                                                      drop_classes_column_name = "EVT_NAME",
                                                      perc_digits = 2,
                                                      raster_return = c("WRITE"))
toc()

tic() #91 minutes
test_full <- full_representative_categorical_analysis(raster = raster_grouped,
                                                      raster_cat_df = raster_cats_grouped,
                                                      region_shape = conus_epa_clean,
                                                      aoi_shape = conus_aoi,
                                                      run_name = "CONUS_grp",
                                                      cat_base_column_name = "VALUE",
                                                      out_rast_values = c("PERC_COVER_AOI", "PERC_COVER_REGION"),
                                                      out_rast_type = c("FULL"),
                                                      out_dir = dir_output,
                                                      new_sub_dir = TRUE,
                                                      min_aoi_coverage = NA,
                                                      min_region_coverage = NA,
                                                      drop_classes = c("Open Water"),
                                                      drop_classes_column_name = "EVT_GP_N",
                                                      perc_digits = 2,
                                                      raster_return = c("WRITE"))
toc()


# By NEON region, grouped and ungrouped

# 
# tic()
# test_full <- full_representative_categorical_analysis(raster = raster,
#                                                       raster_cat_df = raster_cats,
#                                                       region_shape = neon_region_polygons_merged[1,],
#                                                       aoi_shape = neon_areas_of_interest_merged[1,],
#                                                       run_name = "TEST_ALONE",
#                                                       cat_base_column_name = "VALUE",
#                                                       out_rast_values = c("PERC_COVER_AOI", "PERC_COVER_REGION"),
#                                                       out_rast_type = c("FULL"),
#                                                       out_dir = dir_output,
#                                                       new_sub_dir = TRUE,
#                                                       min_aoi_coverage = NA,
#                                                       min_region_coverage = NA,
#                                                       drop_classes = c("Open Water"),
#                                                       drop_classes_column_name = "EVT_NAME",
#                                                       perc_digits = 2,
#                                                       raster_return = c("WRITE"))
# toc()



tic()
test_full <- full_representative_categorical_analysis_set(full_run_nm = "full_run_2",
                                                          dir_out = dir_output,
                                                          region_polygons_merged = neon_region_polygons_merged[],
                                                          areas_of_interest_merged = neon_areas_of_interest_merged[],
                                                          region_name_col = "DomainName",
                                                          parallel = FALSE,
                                                          safe_parallel = FALSE,
                                                          n_workers = 16,
                                                          mem_tot = mem_tot * 0.7,
                                                          raster = raster,
                                                          raster_cat_df = raster_cats,
                                                          cat_base_column_name = "VALUE",
                                                          out_rast_values = c("PERC_COVER_AOI", "PERC_COVER_REGION"),
                                                          out_rast_type = c("FULL"),
                                                          new_sub_dir = TRUE,
                                                          min_aoi_coverage = NA,
                                                          min_region_coverage = NA,
                                                          drop_classes = c("Open Water"),
                                                          drop_classes_column_name = "EVT_NAME",
                                                          perc_digits = 2,
                                                          raster_return = c("WRITE"))
toc()



tic()
test_full <- full_representative_categorical_analysis_set(full_run_nm = "full_run_grp",
                                                          dir_out = dir_output,
                                                          region_polygons_merged = neon_region_polygons_merged[8:16,],
                                                          areas_of_interest_merged = neon_areas_of_interest_merged[8:16,],
                                                          region_name_col = "DomainName",
                                                          parallel = FALSE,
                                                          safe_parallel = FALSE,
                                                          n_workers = 16,
                                                          mem_tot = mem_tot * 0.7,
                                                          raster = raster_grouped,
                                                          raster_cat_df = raster_cats_grouped,
                                                          cat_base_column_name = "VALUE",
                                                          out_rast_values = c("PERC_COVER_AOI", "PERC_COVER_REGION"),
                                                          out_rast_type = c("FULL"),
                                                          new_sub_dir = TRUE,
                                                          min_aoi_coverage = NA,
                                                          min_region_coverage = NA,
                                                          drop_classes = c("Open Water"),
                                                          drop_classes_column_name = "EVT_GP_N",
                                                          perc_digits = 2,
                                                          raster_return = c("WRITE"))
toc()

system("cp -r ~/lens-aop-continental-scaling/data/output/full_run_grp ~/data-store/data/iplant/home/shared/earthlab/macrosystems/lens-aop-continental-scaling/data/output/")



if(cyverse) {
  system("cp -r ~/lens-aop-continental-scaling/data/output ~/data-store/data/iplant/home/shared/earthlab/macrosystems/lens-aop-continental-scaling/data")
  #system("cp -r ~/lens-aop-continental-scaling/figs ~/data-store/data/iplant/home/shared/earthlab/macrosystems/lens-aop-continental-scaling")
}

write_session_info(here::here(dir_output, 'session_info.txt'))






          