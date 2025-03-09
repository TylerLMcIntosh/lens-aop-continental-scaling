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
    "furrr",
    "pals",
    "biscale",
    "grid"),
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
install_and_load_packages(c("mapedit"))

test <- mapedit::drawFeatures()
test <- test |>
  sf::st_transform(sf::st_crs(neon_region_polygons))


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


conus_aoi <- areas_of_interest |>
  sf::st_intersection(conus_epa_clean) |>
  dplyr::mutate(conus = "conus") |>
  dplyr::group_by(conus) |>
  dplyr::summarise(geometry = sf::st_union(geometry)) |>
  dplyr::ungroup()


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









t <- representative_categorical_cover_analysis(raster = raster,
                                               raster_cat_df = raster_cats,
                                               region_shape = test,
                                               aoi_shape = areas_of_interest |> dplyr::filter(siteID == "NIWO"),
                                               run_name = "TUESDAY_CRUISDAY2",
                                               cat_base_column_name = "VALUE",
                                               out_rast_values = c("RAW", "PERC_COVER_AOI", "PERC_COVER_REGION"),
                                               out_rast_type = c("REP", "NOT_REP", "FULL"),
                                               out_dir = here::here("data/derived"),
                                               new_sub_dir = TRUE,
                                               min_aoi_coverage = NA,
                                               min_region_coverage = NA,
                                               drop_classes = NA,
                                               drop_classes_column_name = NA,
                                               perc_digits = 2,
                                               raster_return = c("WRITE", "MEMORY"))

t_gp <- representative_categorical_cover_analysis(raster = raster,
                                               raster_cat_df = raster_cats,
                                               region_shape = test,
                                               aoi_shape = areas_of_interest |> dplyr::filter(siteID == "NIWO"),
                                               run_name = "TUESDAY_CRUISDAY2",
                                               cat_base_column_name = "EVT_GP",
                                               out_rast_values = c("RAW", "PERC_COVER_AOI", "PERC_COVER_REGION"),
                                               out_rast_type = c("REP", "NOT_REP", "FULL"),
                                               out_dir = here::here("data/derived"),
                                               new_sub_dir = TRUE,
                                               min_aoi_coverage = NA,
                                               min_region_coverage = NA,
                                               drop_classes = NA,
                                               drop_classes_column_name = NA,
                                               perc_digits = 2,
                                               raster_return = c("WRITE", "MEMORY"))



x <- t$df_raw
gp <- t_gp$df_raw







t <- representative_categorical_cover_analysis(raster = raster,
                                               raster_cat_df = raster_cats,
                                               region_shape = neon_region_polygons_merged[11,],
                                               aoi_shape = neon_areas_of_interest_merged[11,],
                                               run_name = "TEST_NorthernRockies_GP",
                                               cat_base_column_name = "VALUE",
                                               out_rast_values = c("PERC_COVER_AOI", "PERC_COVER_REGION"),
                                               out_rast_type = c("FULL"),
                                               out_dir = here::here("data/derived"),
                                               new_sub_dir = TRUE,
                                               min_aoi_coverage = NA,
                                               min_region_coverage = NA,
                                               drop_classes = NA,
                                               drop_classes_column_name = NA,
                                               perc_digits = 2,
                                               raster_return = c("WRITE"))



# tic()
# t <- representative_categorical_cover_analysis(raster = raster,
#                                                raster_cat_df = raster_cats,
#                                                region_shape = conus_epa_clean,
#                                                aoi_shape = conus_aoi,
#                                                run_name = "TEST_CONUS",
#                                                cat_base_column_name = "VALUE",
#                                                out_rast_values = c("PERC_COVER_AOI", "PERC_COVER_REGION"),
#                                                out_rast_type = c("FULL"),
#                                                out_dir = here::here("data/derived"),
#                                                new_sub_dir = TRUE,
#                                                min_aoi_coverage = NA,
#                                                min_region_coverage = NA,
#                                                drop_classes = NA,
#                                                drop_classes_column_name = NA,
#                                                perc_digits = 2,
#                                                raster_return = c("WRITE"))
# toc()



# REPRESENT GRAPHS


results <- represent_stats(coverage_df = t$df_raw,
                           manual_lim = NA)


ggsave(plot = results$represent_plot_normal,
       width = 6,
       height = 6,
       units = "in",
       filename = here::here('data/derived/TEST_NorthernRockies_aoi0_region0/rep_plot_normal.jpg'))
ggsave(plot = results$represent_plot_log,
       width = 6,
       height = 6,
       units = "in",
      filename = here::here('data/derived/TEST_NorthernRockies_aoi0_region0/rep_plot_log.jpg'))
ggsave(plot = results$cumulative_plot,
       width = 6,
       height = 6,
       units = "in",
       filename = here::here('data/derived/TEST_NorthernRockies_aoi0_region0/cumulative_plot.jpg'))






# USE FUNCTIONS
region <- terra::rast(t$raster_file_names$full$PERC_COVER_REGION)
aoi <- terra::rast(t$raster_file_names$full$PERC_COVER_AOI)

# region <- t$rasters$full$PERC_COVER_REGION
# aoi <- t$rasters$full$PERC_COVER_AOI

tic()
bp_norm2_agg <- bivariate_raster_viz_3(x = terra::aggregate(region, fact = 30, fun = "median"),
                                       y = terra::aggregate(aoi, fact = 30, fun = "median"),
                                       bi_normal = TRUE,
                                       pals_pal = pals::brewer.seqseq2(n = 9),
                                       flip = FALSE,
                                       x_nm = "Region Coverage",
                                       y_nm = "AOI Coverage",
                                       title = "Aggregated bivariate coverage map, bi-normalized")
toc()


tic()
bp_norm1_agg <- bivariate_raster_viz_3(x = terra::aggregate(region, fact = 30, fun = "median"),
                             y = terra::aggregate(aoi, fact = 30, fun = "median"),
                             bi_normal = FALSE,
                             pals_pal = pals::brewer.seqseq2(n = 9),
                             flip = FALSE,
                             x_nm = "Region Coverage",
                             y_nm = "AOI Coverage",
                             title = "Aggregated bivariate coverage map, uni-normalized")
toc()


tic()
bp_norm2_log_agg <- bivariate_raster_viz_3(x = log1p(terra::aggregate(region, fact = 30, fun = "median")),
                             y = log1p(terra::aggregate(aoi, fact = 30, fun = "median")),
                             bi_normal = TRUE,
                             pals_pal = pals::brewer.seqseq2(n = 9),
                             flip = FALSE,
                             x_nm = "Log-transformed Region Coverage",
                             y_nm = "Log-transformed AOI Coverage",
                             title = "Aggregated and log-transformed bivariate coverage map, bi-normalized")
toc()



tic()
bp_norm1_log_agg <- bivariate_raster_viz_3(x = log1p(terra::aggregate(region, fact = 30, fun = "median")),
                             y = log1p(terra::aggregate(aoi, fact = 30, fun = "median")),
                             bi_normal = FALSE,
                             pals_pal = pals::brewer.seqseq2(n = 9),
                             flip = FALSE,
                             x_nm = "Log-transformed Region Coverage",
                             y_nm = "Log-transformed AOI Coverage",
                             title = "Aggregated and log-transformed bivariate coverage map, uni-normalized")
toc()



ggsave(plot = bp_norm1_agg$legend,
       filename = here::here("data/derived/TEST_NorthernRockies_aoi0_region0/biv_legend.jpg"))
tmap_save(tm = bp_norm2_agg$biv_plot,
          filename = here::here("data/derived/TEST_NorthernRockies_aoi0_region0/biv_norm2.jpg"))
tmap_save(tm = bp_norm1_agg$biv_plot,
          filename = here::here("data/derived/TEST_NorthernRockies_aoi0_region0/biv_norm1.jpg"))
tmap_save(tm = bp_norm2_log_agg$biv_plot,
          filename = here::here("data/derived/TEST_NorthernRockies_aoi0_region0/biv_norm2_log.jpg"))
tmap_save(tm = bp_norm1_log_agg$biv_plot,
          filename = here::here("data/derived/TEST_NorthernRockies_aoi0_region0/biv_norm1_log.jpg"))







