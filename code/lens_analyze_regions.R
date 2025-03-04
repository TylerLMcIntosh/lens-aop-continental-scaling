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

# ## Cyverse data store access if applicable ----
# # Copy previous data over from data store if on cyverse and project has been run before
# if(cyverse) {
#   if(dir.exists("~/data-store/data/iplant/home/shared/earthlab/macrosystems/lens-aop-continental-scaling/data")) {
#     system("cp -r ~/data-store/data/iplant/home/shared/earthlab/macrosystems/lens-aop-continental-scaling/data ~/lens-aop-continental-scaling/")
#   }
#   if(dir.exists("~/data-store/data/iplant/home/shared/earthlab/macrosystems/lens-aop-continental-scaling/figs")) {
#     system("cp -r ~/data-store/data/iplant/home/shared/earthlab/macrosystems/lens-aop-continental-scaling/figs ~/lens-aop-continental-scaling/")
#   }
# }


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
epa_region_polygons <- access_data_epa_l2_ecoregions_api() |>
  sf::st_transform(sf::st_crs(neon_region_polygons))
areas_of_interest <- access_neon_aop_flight_box_data() #note that flight boxes have the domain data as "D##" instead of just "##"

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
# test <- mapedit::drawFeatures()
# test <- test |>
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



#Set up to use EPA ecoregions as well
epa_region_polygons_merged <- epa_region_polygons |>
  sf::st_cast("POLYGON") |> #split out multipolygons
  sf::st_make_valid() %>%
  dplyr::filter(sf::st_intersects(., conus, sparse = FALSE)) |>
  dplyr::group_by(NA_L2CODE, NA_L2KEY) |>
  dplyr::summarise(geometry = sf::st_union(geometry)) |>
  dplyr::ungroup() |>
  dplyr::arrange(NA_L2KEY) |>
  dplyr::filter(NA_L2CODE != "0.0") |>
  sf::st_intersection(conus)

epa_areas_of_interest_merged <- areas_of_interest |>
  sf::st_intersection(epa_region_polygons_merged) |>
  dplyr::group_by(NA_L2CODE, NA_L2KEY) |>
  dplyr::summarise(geometry = sf::st_union(geometry)) |>
  dplyr::ungroup() |>
  dplyr::arrange(NA_L2KEY)

#Ensure that the same sets are in both
common_epa_ids <- lubridate::intersect(epa_region_polygons_merged$NA_L2KEY, epa_areas_of_interest_merged$NA_L2KEY)
epa_region_polygons_merged <- epa_region_polygons_merged %>%
  dplyr::filter(NA_L2KEY %in% common_epa_ids)
epa_areas_of_interest_merged <- epa_areas_of_interest_merged %>%
  dplyr::filter(NA_L2KEY %in% common_epa_ids)

epa_region_polygons_merged_file <- here::here(dir_derived, "epa_region_polygons_merged.gpkg")
sf::st_write(epa_region_polygons_merged, epa_region_polygons_merged_file, append = FALSE)

epa_areas_of_interest_merged_file <- here::here(dir_derived, "epa_areas_of_interest_merged.gpkg")
sf::st_write(epa_areas_of_interest_merged, epa_areas_of_interest_merged_file, append = FALSE)


## Prep to run analyses ----
aoi_thresholds <- c(0.001, 0.01, 0.1, 1)

ag_dev_mine_evt_names <- raster_cats |> dplyr::filter(
  grepl(pattern = "Developed", x = EVT_PHYS, ignore.case = FALSE) |
    grepl(pattern = "Agricultural", x = EVT_PHYS, ignore.case = FALSE) |
    grepl(pattern = "Quarries", x = EVT_PHYS, ignore.case = FALSE)) |>
  dplyr::pull(EVT_NAME) |>
  as.vector()

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
                                                     areas_of_interest_merged = neon_areas_of_interest_merged_file,
                                                     region_name_col = "DomainName",
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
             .f = function(x, y) conus_lens_analysis(region_polygons_merged = neon_region_polygons_merged,
                                                     areas_of_interest_merged = neon_areas_of_interest_merged,
                                                     region_name_col = "DomainName",
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
             .f = function(x, y) conus_lens_analysis(region_polygons_merged = neon_region_polygons_merged[1,],
                                                     areas_of_interest_merged = neon_areas_of_interest_merged[1,],
                                                     region_name_col = "DomainName",
                                                     raster = raster,
                                                     raster_cat_df = raster_cats,
                                                     run_name = y,
                                                     cat_base_column_name = "EVT_GP",
                                                     aoi_drop_perc = x,
                                                     drop_classes = NA,
                                                     drop_classes_column_name = NA,
                                                     out_rast_values = "PERC_COVER",
                                                     out_rast_type = "NOT_REP",
                                                     out_dir = here::here('data/derived/')))
toc()



# Run analysis with AOI threshold of 0.001% - 1% and EVT groups except ag_dev_mine
tic()
purrr::walk2(.x = aoi_thresholds,
             .y = c("neon_domains_evt_groups_no_ag_dev_mine_0001",
                    "neon_domains_evt_groups_no_ag_dev_mine_001",
                    "neon_domains_evt_groups_no_ag_dev_mine_01",
                    "neon_domains_evt_groups_no_ag_dev_mine_1"),
             .f = function(x, y) conus_lens_analysis(region_polygons_merged = neon_region_polygons_merged,
                                                     areas_of_interest_merged = neon_areas_of_interest_merged,
                                                     region_name_col = "DomainName",
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
             .f = function(x, y) conus_lens_analysis(region_polygons_merged = neon_region_polygons_merged,
                                                     areas_of_interest_merged = neon_areas_of_interest_merged,
                                                     region_name_col = "DomainName",
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
             .f = function(x, y) conus_lens_analysis(region_polygons_merged = neon_region_polygons_merged,
                                                     areas_of_interest_merged = neon_areas_of_interest_merged,
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






















# TESTING



#' Perform Representative Categorical Cover Analysis
#'
#' This function analyzes the categorical cover representation within a region and area of interest (AOI).
#' It computes the percentage of each class within the AOI and the larger region, then determines which 
#' classes are well represented or underrepresented based on a given threshold.
#'
#' @param raster A `SpatRaster` object containing categorical cover data.
#' @param raster_cat_df A `data.frame` mapping raster values to their respective categories.
#' @param region_shape A `sf` object representing the broader region of analysis.
#' @param aoi_shape A `sf` object representing the area of interest (AOI) to analyze.
#' @param run_name A `character` string specifying the name of the analysis run, which is used to name outputs and output directories.
#' @param cat_base_column_name A `character` string indicating the column in `raster_cat_df` that contains the categorical data values.
#' @param aoi_drop_perc A `numeric` value specifying the minimum percentage of a cover class within the AOI to be considered "represented." Defaults to `NA` (no threshold applied).
#' @param region_drop_perc A `numeric` value specifying the minimum percentage of a cover class within the larger region for inclusion. Defaults to `NA` (no filtering applied).
#' @param drop_classes A `vector` of category values to exclude from analysis. Defaults to `NA` (no classes dropped).
#' @param drop_classes_column_name A `character` string specifying the column in `raster_cat_df` used to filter `drop_classes`. Defaults to `NA`.
#' @param out_rast_values A `character` or `vector` of raster output types. Options:
#'   - `"RAW"`: Outputs the raw clipped categorical raster.
#'   - `"PERC_COVER_AOI"`: Outputs a raster where each pixel represents the percentage of the given class within the AOI.
#'   - `"PERC_COVER_REGION"`: Outputs a raster where each pixel represents the percentage of the given class within the region.
#' @param out_rast_type A `character` or `vector` specifying which raster types to output. Options:
#'   - `"REP"`: Raster of represented classes.
#'   - `"NOT_REP"`: Raster of underrepresented classes.
#'   - `"FULL"`: Raster including all classes, with no thresholding.
#'   - `"NONE"`: No raster output.
#' @param out_dir A `character` string specifying the directory where rasters should be saved.
#' @param new_sub_dir A `logical` value. If `TRUE`, creates a new subdirectory for output based on `run_name`. Defaults to `FALSE`.
#' @param perc_digits An `integer` specifying the number of decimal places for percentage values in raster outputs. If `NA`, values remain unrounded.
#' @param rasters_in_memory A `logical` value. If `TRUE`, returns raster outputs in memory; if `FALSE`, writes them to disk.
#'
#' @return A named `list` containing:
#'   - `df_raw`: The full categorical cover analysis results.
#'   - `df_included`: The filtered categorical cover analysis results.
#'   - `df_represented`: The subset of represented classes.
#'   - `df_not_represented`: The subset of underrepresented classes.
#'   - `perc_area_not_represented`: The total percentage of the region covered by underrepresented classes.
#'   - `rasters`: If `rasters_in_memory = TRUE`, a list of raster outputs.
#'
#' @export
representative_categorical_cover_analysis <- function(raster,
                                                      raster_cat_df,
                                                      region_shape,
                                                      aoi_shape,
                                                      run_name = "NotProvided",
                                                      cat_base_column_name, 
                                                      aoi_drop_perc = NA,
                                                      region_drop_perc = NA,
                                                      drop_classes = NA,
                                                      drop_classes_column_name = NA,
                                                      out_rast_values = c("RAW", "PERC_COVER_AOI", "PERC_COVER_REGION"),
                                                      out_rast_type = c("REP", "NOT_REP"),
                                                      out_dir = "",
                                                      new_sub_dir = FALSE,
                                                      perc_digits = NA,
                                                      rasters_in_memory = FALSE) {
  
  print(paste0("Operating on run: ", run_name))
  
  # Validate input parameters
  valid_rast_values <- c("RAW", "PERC_COVER_AOI", "PERC_COVER_REGION")
  valid_rast_types <- c("REP", "NOT_REP", "FULL", "NONE")
  
  if (!all(out_rast_values %in% valid_rast_values)) {
    stop("Invalid out_rast_values. Must be one or more of: 'RAW', 'PERC_COVER_AOI', 'PERC_COVER_REGION'.")
  }
  
  if (!all(out_rast_type %in% valid_rast_types)) {
    stop("Invalid out_rast_type. Must be one or more of: 'REP', 'NOT_REP', 'FULL', 'NONE'.")
  }
  
  if (!cat_base_column_name %in% names(raster_cat_df)) {
    stop("cat_base_column_name must be a column in raster_cat_df")
  }
  
  print(paste0("Operating on run: ", run_name))
  
  # Validate input parameters
  valid_rast_values <- c("RAW", "PERC_COVER_AOI", "PERC_COVER_REGION")
  valid_rast_types <- c("REP", "NOT_REP", "FULL", "NONE")
  
  if (!all(out_rast_values %in% valid_rast_values)) {
    stop("Invalid out_rast_values. Must be one or more of: 'RAW', 'PERC_COVER_AOI', 'PERC_COVER_REGION'.")
  }
  
  if (!all(out_rast_type %in% valid_rast_types)) {
    stop("Invalid out_rast_type. Must be one or more of: 'REP', 'NOT_REP', 'FULL', 'NONE'.")
  }
  
  if (!cat_base_column_name %in% names(raster_cat_df)) {
    stop("cat_base_column_name must be a column in raster_cat_df")
  }
  
  # Set up output directory
  clean_run_name <- gsub(" ", "", run_name)
  clean_aoi_dp <- gsub("\\.", "", as.character(aoi_drop_perc))
  clean_region_dp <- gsub("\\.", "", as.character(region_drop_perc))
  clean_run_name <- paste(clean_run_name, "_adp", clean_aoi_dp, "_rdp", clean_region_dp, sep = "")
  
  if (new_sub_dir) {
    out_dir <- here::here(out_dir, clean_run_name)
    dir_ensure(out_dir)
  }
  
  # Crop to region and AOI
  print("Cropping to region")
  larger_region_cover <- crop_careful_universal(raster, region_shape, mask = TRUE)
  print("Cropping to AOI")
  aoi_cover <- crop_careful_universal(larger_region_cover, aoi_shape, mask = TRUE)
  
  # Perform categorical cover analysis
  landcover_analysis_output_raw <- analyze_categorical_cover(aoi_cover, larger_region_cover, raster_cat_df, cat_base_column_name)
  
  # Process drop classes and thresholding
  landcover_analysis_output_included <- landcover_analysis_output_raw
  if (!is.na(region_drop_perc)) {
    landcover_analysis_output_included <- landcover_analysis_output_included |> dplyr::filter(region_perc > region_drop_perc)
  }
  if (length(drop_classes) > 0 && !all(is.na(drop_classes))) {
    landcover_analysis_output_included <- landcover_analysis_output_included |> 
      dplyr::filter(!(.data[[drop_classes_column_name]] %in% drop_classes))
  }
  
  # Split into represented and not represented classes
  df_represented <- landcover_analysis_output_included |> dplyr::filter(aoi_perc > aoi_drop_perc)
  df_not_represented <- landcover_analysis_output_included |> dplyr::filter(aoi_perc <= aoi_drop_perc)
  
  # Compute percent not represented
  perc_area_not_represented <- (sum(df_not_represented$region_count) / sum(landcover_analysis_output_included$region_count)) * 100
  
  # 🛠️ **Create filtered rasters for REP & NOT_REP cases**
  raster_represented <- terra::classify(larger_region_cover, 
                                        as.matrix(df_represented[, c(cat_base_column_name, cat_base_column_name)]), 
                                        others = NA)
  
  raster_not_represented <- terra::classify(larger_region_cover, 
                                            as.matrix(df_not_represented[, c(cat_base_column_name, cat_base_column_name)]), 
                                            others = NA)
  
  # Save rasters
  raster_outputs <- list()
  if ("FULL" %in% out_rast_type) {
    raster_outputs$full <- save_rasters(larger_region_cover, landcover_analysis_output_raw, "full",
                                        out_dir, clean_run_name, out_rast_values, perc_digits, rasters_in_memory, cat_base_column_name)
  }
  if ("REP" %in% out_rast_type) {
    raster_outputs$rep <- save_rasters(raster_represented, df_represented, "rep",
                                       out_dir, clean_run_name, out_rast_values, perc_digits, rasters_in_memory, cat_base_column_name)
  }
  if ("NOT_REP" %in% out_rast_type) {
    raster_outputs$not_rep <- save_rasters(raster_not_represented, df_not_represented, "not_rep",
                                           out_dir, clean_run_name, out_rast_values, perc_digits, rasters_in_memory, cat_base_column_name)
  }
  
  return(list(analysis_name = run_name,
              df_raw = landcover_analysis_output_raw,
              df_included = landcover_analysis_output_included,
              df_represented = df_represented,
              df_not_represented = df_not_represented,
              perc_area_not_represented = perc_area_not_represented,
              rasters = if (rasters_in_memory) raster_outputs else NULL))
}


#' Save Raster Outputs from Categorical Cover Analysis - A helper function for function representative_categorical_cover_analysis
#'
#' Saves raster outputs in various formats based on user selection.
#'
#' @param raster A `SpatRaster` object to be saved.
#' @param df A `data.frame` with categorical data used to generate raster values.
#' @param raster_type A `character` string indicating the type of raster (`"rep"`, `"not_rep"`, `"full"`).
#' @param out_dir A `character` string specifying the output directory.
#' @param run_name A `character` string specifying the run name.
#' @param out_rast_values A `character` or `vector` of output raster types:
#'   - `"RAW"`: Outputs the raw raster.
#'   - `"PERC_COVER_AOI"`: Outputs a raster with class percentages within the AOI.
#'   - `"PERC_COVER_REGION"`: Outputs a raster with class percentages within the region.
#' @param perc_digits An `integer` for rounding percentages in output rasters. `NA` means no rounding.
#' @param rasters_in_memory A `logical`. If `TRUE`, returns rasters in memory; if `FALSE`, writes to disk.
#' @param cat_base_column_name A `character` specifying the category column used in classification.
#'
#' @return A named list of raster outputs (if `rasters_in_memory = TRUE`), otherwise writes files and returns `NULL`.
#'
#' @export
save_rasters <- function(raster, df, type, out_dir, clean_run_name, 
                         out_rast_values, perc_digits, rasters_in_memory, 
                         cat_base_column_name) {  
  
  base_path <- here::here(out_dir, paste0(clean_run_name, "_", type))
  raster_list <- list()
  
  for (value_type in out_rast_values) {
    
    # If RAW, directly return/write the clipped raster without reclassification
    if (value_type == "RAW") {
      raster_raw <- raster  # Just retain the original values
      
      if (rasters_in_memory) {
        raster_list[[value_type]] <- raster_raw
      } else {
        output_file <- paste0(base_path, "_raw.tif")
        terra::writeRaster(raster_raw, output_file, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE"))
      }
      
      next  # Skip to the next raster type
    }
    
    # Determine the appropriate metric column
    metric_column <- switch(value_type,
                            "PERC_COVER_AOI" = "aoi_perc",
                            "PERC_COVER_REGION" = "region_perc",
                            stop(glue::glue("Invalid out_rast_values option: {value_type}")))
    
    # Ensure the necessary columns exist
    if (!cat_base_column_name %in% colnames(df)) {
      stop(glue::glue("Error: Column '{cat_base_column_name}' not found in the dataframe."))
    }
    
    if (!metric_column %in% colnames(df)) {
      stop(glue::glue("Error: Column '{metric_column}' not found in the dataframe."))
    }
    
    # Apply rounding if needed
    df <- df |>
      dplyr::mutate(!!metric_column := if (!is.na(perc_digits)) round(.data[[metric_column]], perc_digits) else .data[[metric_column]])
    
    # Create the classification matrix
    reclass_matrix <- df |>
      dplyr::select(all_of(cat_base_column_name), all_of(metric_column)) |> 
      as.matrix()
    
    # Classify the raster
    raster_reclassified <- terra::classify(raster, reclass_matrix)
    
    if (rasters_in_memory) {
      raster_list[[value_type]] <- raster_reclassified
    } else {
      output_file <- paste0(base_path, "_", tolower(value_type), ".tif")
      terra::writeRaster(raster_reclassified, output_file, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE"))
    }
  }
  
  if (rasters_in_memory) {
    return(raster_list)
  } else {
    return(invisible(NULL))
  }
}













t <- representative_categorical_cover_analysis(raster = raster,
                                               raster_cat_df = raster_cats,
                                               region_shape = test,
                                               aoi_shape = areas_of_interest |> dplyr::filter(siteID == "NIWO"),
                                               run_name = "TEST_D2_2",
                                               cat_base_column_name = "VALUE",
                                               out_rast_values = c("RAW", "PERC_COVER_AOI", "PERC_COVER_REGION"),
                                               out_rast_type = c("REP", "NOT_REP", "FULL"),
                                               out_dir = here::here("data/derived"),
                                               new_sub_dir = TRUE,
                                               aoi_drop_perc = 5,
                                               region_drop_perc = NA,
                                               drop_classes = NA,
                                               drop_classes_column_name = NA,
                                               perc_digits = 2,
                                               rasters_in_memory = TRUE)

