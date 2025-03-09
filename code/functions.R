# Tyler L. McIntosh

# LENS analyses ----




# A specific version of representative_categorical_cover_analysis that operates over a full set of matched
# region and AOI polygons, exports results, creates a csv summary, and a single CONUS-wide graphic
conus_lens_analysis <- function(region_polygons_merged,
                                region_name_col,
                                areas_of_interest_merged,
                                raster,
                                raster_cat_df,
                                run_name = "neon_domains",
                                cat_base_column_name,
                                aoi_drop_perc = NA,
                                drop_classes = NA,
                                drop_classes_column_name = NA,
                                out_rast_values = "PERC_COVER",
                                out_rast_type = "NOT_REP",
                                out_dir) {
  
  tic(paste("CONUS-lens-analysis run:", run_name, sep = " "))
  
  # Allow ease of parallel processing
  if(is.character(raster)) {
    raster <- terra::rast(raster)
  }
  
  if(is.character(region_polygons_merged)) {
    region_polygons_merged <- sf::st_read(region_polygons_merged)
  }
  
  if(is.character(areas_of_interest_merged)) {
    areas_of_interest_merged <- sf::st_read(areas_of_interest_merged)
  }
  
  # Setup output directory for rasters
  clean_aoi_dp <- gsub("\\.", "", as.character(aoi_drop_perc))
  dir_out <- here::here(out_dir, paste(run_name, clean_aoi_dp, sep = "_"))
  dir_ensure(dir_out)
  
  #Run analysis using representative_categorical_cover_analysis function
  all_region_results <- purrr::pmap(list(region_shape = split(region_polygons_merged, seq(nrow(region_polygons_merged))),
                                         aoi_shape = split(areas_of_interest_merged, seq(nrow(areas_of_interest_merged))),
                                         #run_name = paste(run_name, region_polygons_merged$DomainName, sep ="_")),
                                        run_name = paste(run_name, region_polygons_merged[[region_name_col]], sep = "_")),
                                    representative_categorical_cover_analysis,
                                    raster = raster,
                                    raster_cat_df = raster_cat_df,
                                    cat_base_column_name = cat_base_column_name,
                                    region_drop_perc = 0,
                                    aoi_drop_perc = aoi_drop_perc,
                                    drop_classes = drop_classes,
                                    drop_classes_column_name = drop_classes_column_name,
                                    out_rast_values = out_rast_values,
                                    out_rast_type = out_rast_type, #out_rast_type = "BOTH", "REP", "NOT_REP", or "NONE"
                                    out_dir = dir_out,
                                    new_sub_dir = FALSE)
  
  # Create a dataframe with all percentages and export
  result_df <- purrr::map_dfr(all_region_results, ~ tibble(
    region_name = .x$analysis_name,
    perc_area_not_represented = .x$perc_area_not_represented
  ))
  readr::write_csv(result_df, here::here(dir_out, paste0(run_name, "_results.csv")))
  
  toc()
}



conus_lens_figure <- function(dir_search,
                              pattern,
                              overlay_polygons,
                              name,
                              col_grad = scico::scico("grayC"),
                              downsample = TRUE) {
  
  # Read in the list of tif files to create the CONUS figure
  tif_files <- list.files(dir_search,
                          pattern = pattern,
                          full.names = TRUE,
                          recursive = TRUE,
                          breaks = c(0,15))
  
  # CREATE FIGURE AND SAVE
  tmap_options(max.raster = c(plot = 1e7, view = 1e5))
  conus <- tigris::states(cb = TRUE) |>  # `cb = TRUE` for a simplified "cartographic boundary" version
    dplyr::filter(!STUSPS %in% c("HI", "AK", "GU", "VI", "MP", "AS", "PR")) |>
    sf::st_transform(crs = terra::crs(terra::rast(tif_files[1])))
  
  # Loop through each raster file, simplify if needed, and add to the tmap object
  for (raster_path in tif_files) {
    
    # Load the raster
    r <- stars::read_stars(raster_path,
                           proxy = TRUE)
    
    if(raster_path == tif_files[1]) {
      tm_plot <- 
        tm_shape(overlay_polygons |>
                   sf::st_transform(terra::crs(terra::rast(tif_files[1]))),
                 bbox = sf::st_bbox(conus)) +
        tm_borders(col = "gray90", lwd = 1) +
        tm_fill(col = "gray90") +
        tmap::tm_shape(r,
                       bbox = sf::st_bbox(conus),
                       downsample = downsample) +
        tmap::tm_raster(palette = col_grad,
                        style = "cont",
                        breaks = breaks,
                        legend.show = TRUE,
                        title = "Unrepresented landscape\npercentage by class",
                        legend.reverse = FALSE,
                        legend.format = list(fun = function(x) paste0(x, "%")),
                        legend.is.portrait = FALSE)
    } else {
      # Add the raster to the tmap object
      tm_plot <- tm_plot +
        tmap::tm_shape(r,
                       bbox = sf::st_bbox(conus),
                       downsample = downsample) +
        tmap::tm_raster(palette = col_grad,
                        style = "cont",
                        breaks = breaks,
                        legend.show = FALSE)
    }
  }
  
  # Finalize the plot layout with the legend outside
  tm_plot <- tm_plot +
    tm_shape(overlay_polygons) +
    tm_borders(col = "gray20",
               lwd = 1) +
    tm_fill(col = NA, alpha = 0) +
    tm_shape(neon_areas_of_interest_merged) +
    tm_borders(col = "darkblue",
               lwd = 1) +
    tmap::tm_layout(legend.outside = FALSE,
                    legend.position = c("left", "bottom"),
                    title = name)
  
  # Save the plot
  tmap::tmap_save(tm_plot, here::here(dir_figs, paste0(name, ".jpeg")))
  
}



# Landscape representativeness ----

# Function to get presence of a raster value within NEON base plots
# PARAMETERS
# siteID :: the NEON ID of interest (e.g. "YELL")
# cat_raster :: a raster to pull values from (intended to be a landcover categorical raster), in SpatRaster format
get_base_plot_freqs <- function(site_id, cat_raster) {
  epsg <- terra::crs(cat_raster)
  plots <- access_neon_plots_shp() |>
    sf::st_transform(epsg)
  base_plots <- plots |>
    dplyr::filter(subtype == "basePlot" & site_id == site_id)
  crop_rast <- cat_raster |>
    terra::crop(y = base_plots, mask = TRUE, touches = FALSE)
  freqs <- terra::freq(crop_rast)  
  
  return(freqs)
}


#' Condense Frequency Table Groups
#'
#' This function takes a frequency table, typically returned by `terra::freq()`, and condenses
#' rows if there are duplicate values by summing their associated counts. If the table already 
#' contains unique values, no changes are made.
#'
#' @param freq_dats A data frame or tibble returned by `terra::freq()` that contains two columns: 
#'   `value` (the value names) and `count` (the corresponding frequency of each value).
#'
#' @return A condensed version of the frequency table where duplicate values have been merged 
#'   and their counts summed. If no duplicates are found, the original table is returned.
#' 
#' @importFrom dplyr group_by summarise ungroup
#'
#' @examples
#' # Example usage with a mock frequency table
#' freq_table <- data.frame(
#'   value = c(1, 1, 2, 3, 3, 3),
#'   count = c(10, 15, 5, 2, 8, 6)
#' )
#' 
#' condense_freq_groups(freq_table)
#'
#' @export
condense_freq_groups <- function(freq_dats) {
  if (length(unique(freq_dats$value)) == length(freq_dats$value)) {
    print('Frequency table is already correct')
    return(freq_dats)
  } else {
    print('Condensing frequency table')
    freq_dats <- freq_dats %>%
      group_by(value) %>%
      summarise(count = sum(count)) %>%
      ungroup()
    return(freq_dats)
  }
}



#' Perform Representative Categorical Cover Analysis
#'
#' This function analyzes the categorical cover representation within a region and area of interest (AOI).
#' It computes the percentage of each class within the AOI and the larger region, then determines which 
#' classes are well represented or underrepresented based on a given threshold.
#'
#' @param raster A SpatRaster object containing categorical cover data.
#' @param raster_cat_df A data.frame mapping raster values to their respective categories.
#' @param region_shape A sf object representing the broader region of analysis.
#' @param aoi_shape A sf object representing the area of interest (AOI) to analyze.
#' @param run_name A character string specifying the name of the analysis run, which is used to name outputs and output directories.
#' @param cat_base_column_name A character string indicating the column in raster_cat_df that contains the categorical data values.
#' @param min_aoi_coverage A numeric value specifying the minimum percentage of a cover class within the AOI to be considered "represented." Defaults to NA (no threshold applied).
#' @param min_region_coverage A numeric value specifying the minimum percentage of a cover class within the larger region for inclusion. Defaults to NA (no filtering applied).
#' @param drop_classes A vector of category values to exclude from analysis. Defaults to NA (no classes dropped).
#' @param drop_classes_column_name A character string specifying the column in raster_cat_df used to filter drop_classes. Defaults to NA.
#' @param out_rast_values A character or vector of raster output types. Options:
#'   - "RAW": Outputs the raw clipped categorical raster.
#'   - "PERC_COVER_AOI": Outputs a raster where each pixel represents the percentage of the given class within the AOI.
#'   - "PERC_COVER_REGION": Outputs a raster where each pixel represents the percentage of the given class within the region.
#' @param out_rast_type A character or vector specifying which raster types to output. Options:
#'   - "REP": Raster of represented classes.
#'   - "NOT_REP": Raster of underrepresented classes.
#'   - "FULL": Raster including all classes, with no thresholding.
#'   - "NONE": No raster output.
#' @param out_dir A character string specifying the directory where rasters should be saved.
#' @param new_sub_dir A logical value. If TRUE, creates a new subdirectory for output based on run_name. Defaults to FALSE.
#' @param perc_digits An integer specifying the number of decimal places for percentage values in raster outputs. If NA, values remain unrounded.
#' @param raster_return A character or vector of raster return options:
#'   - "MEMORY": Returns rasters in memory.
#'   - "WRITE": Writes to disk.
#'
#' @return A named list containing:
#'   - df_raw: The full categorical cover analysis results.
#'   - df_included: The filtered categorical cover analysis results.
#'   - df_represented: The subset of represented classes.
#'   - df_not_represented: The subset of underrepresented classes.
#'   - perc_area_not_represented: The total percentage of the region covered by underrepresented classes as defined by the thresholds.
#'   - rasters: Ifrasters_return contains "MEMORY", a list of raster outputs.
#'   - raster_file_names: If rasters_return contains "WRITE", a list of raster file names for written rasters.
#'
#' @export

representative_categorical_cover_analysis <- function(raster,
                                                      raster_cat_df,
                                                      region_shape,
                                                      aoi_shape,
                                                      run_name = "NotProvided",
                                                      cat_base_column_name, 
                                                      min_aoi_coverage = NA,
                                                      min_region_coverage = NA,
                                                      drop_classes = NA,
                                                      drop_classes_column_name = NA,
                                                      out_rast_values = c("RAW", "PERC_COVER_AOI", "PERC_COVER_REGION"),
                                                      out_rast_type = c("REP", "NOT_REP", "FULL"),
                                                      out_dir = "",
                                                      new_sub_dir = FALSE,
                                                      perc_digits = NA,
                                                      raster_return = "MEMORY") {
  
  print(paste0("Operating on run: ", run_name))
  
  # Ensure min_aoi_coverage and min_region_coverage are within valid ranges
  min_aoi_coverage <- ifelse(is.na(min_aoi_coverage) | min_aoi_coverage < 0, 0, min(min_aoi_coverage, 100))
  min_region_coverage <- ifelse(is.na(min_region_coverage) | min_region_coverage < 0, 0, min(min_region_coverage, 100))
  
  # Set up output directory
  clean_run_name <- gsub(" ", "", run_name)
  clean_aoi_dp <- gsub("\\.", "", as.character(min_aoi_coverage))
  clean_region_dp <- gsub("\\.", "", as.character(min_region_coverage))
  clean_run_name <- paste(clean_run_name, "_aoi", clean_aoi_dp, "_region", clean_region_dp, sep = "")
  
  if (new_sub_dir & "WRITE" %in% raster_return) {
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
  
  # Determine max observed percentages
  max_region_perc <- max(landcover_analysis_output_raw$region_perc, na.rm = TRUE)
  max_aoi_perc <- max(landcover_analysis_output_raw$aoi_perc, na.rm = TRUE)
  
  # Adjust extreme min_region_coverage and min_aoi_coverage values
  if (min_region_coverage > max_region_perc) {
    warning(glue::glue("min_region_coverage ({min_region_coverage}%) exceeds all observed region percentages. Adjusting to {max_region_perc}%."))
    min_region_coverage <- max_region_perc
  }
  
  if (min_aoi_coverage > max_aoi_perc) {
    warning(glue::glue("min_aoi_coverage ({min_aoi_coverage}%) exceeds all observed AOI percentages. Adjusting to {max_aoi_perc}%."))
    min_aoi_coverage <- max_aoi_perc
  }
  
  # Process drop classes and thresholding
  landcover_analysis_output_included <- landcover_analysis_output_raw
  if (!is.na(min_region_coverage)) {
    landcover_analysis_output_included <- landcover_analysis_output_included |> dplyr::filter(region_perc > min_region_coverage)
  }
  if (length(drop_classes) > 0 && !all(is.na(drop_classes))) {
    landcover_analysis_output_included <- landcover_analysis_output_included |> 
      dplyr::filter(!(.data[[drop_classes_column_name]] %in% drop_classes))
  }
  
  # Split into represented and not represented classes
  df_represented <- landcover_analysis_output_included |> dplyr::filter(aoi_perc > min_aoi_coverage)
  df_not_represented <- landcover_analysis_output_included |> dplyr::filter(aoi_perc <= min_aoi_coverage)
  
  # Compute percent not represented
  perc_area_not_represented <- ifelse(nrow(landcover_analysis_output_included) > 0,
                                      (sum(df_not_represented$region_count) / sum(landcover_analysis_output_included$region_count)) * 100,
                                      NA)
  
  # **Create filtered rasters for REP & NOT_REP cases**
  create_empty_raster <- function(raster) terra::ifel(!is.na(raster), NA, NA)
  
  if (nrow(df_represented) > 0) {
    raster_represented <- terra::classify(
      larger_region_cover,
      as.matrix(df_represented[, c(cat_base_column_name, cat_base_column_name)]),
      others = NA
    )
  } else {
    warning("No classes met the 'represented' threshold. Creating an empty raster.")
    raster_represented <- create_empty_raster(larger_region_cover)
  }
  
  if (nrow(df_not_represented) > 0) {
    raster_not_represented <- terra::classify(
      larger_region_cover,
      as.matrix(df_not_represented[, c(cat_base_column_name, cat_base_column_name)]),
      others = NA
    )
  } else {
    warning("No classes met the 'not represented' threshold. Creating an empty raster.")
    raster_not_represented <- create_empty_raster(larger_region_cover)
  }
  
  # Save rasters
  raster_outputs <- list()
  raster_files <- list()
  # 
  # if ("FULL" %in% out_rast_type) {
  #   raster_outputs$full <- save_rasters(larger_region_cover, landcover_analysis_output_raw, "full",
  #                                       out_dir, run_name, out_rast_values, perc_digits, raster_return,
  #                                       cat_base_column_name)
  #   raster_files$full <- raster_outputs$full$raster_files
  # }
  # if ("REP" %in% out_rast_type) {
  #   raster_outputs$rep <- save_rasters(raster_represented, df_represented, "rep",
  #                                      out_dir, run_name, out_rast_values, perc_digits, raster_return,
  #                                      cat_base_column_name)
  #   raster_files$rep <- raster_outputs$rep$raster_files
  # }
  # if ("NOT_REP" %in% out_rast_type) {
  #   raster_outputs$not_rep <- save_rasters(raster_not_represented, df_not_represented, "not_rep",
  #                                          out_dir, run_name, out_rast_values, perc_digits, raster_return,
  #                                          cat_base_column_name)
  #   raster_files$not_rep <- raster_outputs$not_rep$raster_files
  # }
  
  
  if ("FULL" %in% out_rast_type) {
    result <- save_rasters(larger_region_cover, landcover_analysis_output_raw, "full",
                           out_dir, run_name, out_rast_values, perc_digits, raster_return,
                           cat_base_column_name)
    raster_outputs$full <- result$raster_list
    raster_files$full <- result$raster_files
  }
  if ("REP" %in% out_rast_type) {
    result <- save_rasters(larger_region_cover, df_represented, "rep",
                           out_dir, run_name, out_rast_values, perc_digits, raster_return,
                           cat_base_column_name)
    raster_outputs$rep <- result$raster_list
    raster_files$rep <- result$raster_files
  }
  if ("NOT_REP" %in% out_rast_type) {
    result <- save_rasters(larger_region_cover, df_not_represented, "not_rep",
                           out_dir, run_name, out_rast_values, perc_digits, raster_return,
                           cat_base_column_name)
    raster_outputs$not_rep <- result$raster_list
    raster_files$not_rep <- result$raster_files
  }
  
  params <- list(analysis_name = run_name,
                 cat_base_column_name = cat_base_column_name,
                 region_shape = region_shape,
                 aoi_shape = aoi_shape,
                 min_aoi_coverage = min_aoi_coverage,
                 min_region_coverage = min_region_coverage,
                 drop_classes = drop_classes,
                 drop_classes_column_name = drop_classes_column_name,
                 perc_digits = perc_digits)
  
  
  return(list(params = params,
              df_raw = landcover_analysis_output_raw,
              df_included = landcover_analysis_output_included,
              df_represented = df_represented,
              df_not_represented = df_not_represented,
              perc_area_not_represented = perc_area_not_represented,
              rasters = if ("MEMORY" %in% raster_return) raster_outputs else NULL,
              raster_file_names = if ("WRITE" %in% raster_return) raster_files else NULL))
}




#' Save Raster Outputs from Categorical Cover Analysis - A helper function for function representative_categorical_cover_analysis
#'
#' Saves raster outputs in various formats based on user selection.
#'
#' @param raster A SpatRaster object to be saved.
#' @param df A data.frame with categorical data used to generate raster values.
#' @param raster_type A character string indicating the type of raster ("rep", "not_rep", "full").
#' @param out_dir A character string specifying the output directory.
#' @param run_name A character string specifying the run name.
#' @param out_rast_values A character or vector of output raster types:
#'   - "RAW": Outputs the raw raster.
#'   - "PERC_COVER_AOI": Outputs a raster with class percentages within the AOI.
#'   - "PERC_COVER_REGION": Outputs a raster with class percentages within the region.
#' @param perc_digits An integer for rounding percentages in output rasters. NA means no rounding.
#' @param raster_return A character or vector of raster return options:
#'   - "MEMORY": Returns rasters in memory.
#'   - "WRITE": Writes to disk.
#' @param cat_base_column_name A character specifying the category column used in classification.
#'
#' @return A named list of raster outputs (if rasters_in_memory = TRUE), otherwise writes files and returns NULL.
#'
#' @export

save_rasters <- function(raster, df, type, out_dir, run_name,
                         out_rast_values, perc_digits, raster_return,
                         cat_base_column_name) {
  
  base_path <- here::here(out_dir, paste0(run_name, "_", type))
  raster_list <- list()
  raster_files <- list()  # Store file paths of written rasters
  
  for (value_type in out_rast_values) {
    
    # If RAW, directly return/write the clipped raster without reclassification
    if (value_type == "RAW") {
      raster_raw <- raster  # Just retain the original values
      
      if ("MEMORY" %in% raster_return) {
        raster_list[[value_type]] <- raster_raw
      }
      if ("WRITE" %in% raster_return) {
        output_file <- paste0(base_path, "_raw.tif")
        terra::writeRaster(raster_raw, output_file, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE"))
        raster_files[[value_type]] <- output_file
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
      warning(glue::glue("Skipping {value_type} raster: Column '{cat_base_column_name}' not found in dataframe."))
      next
    }
    
    if (!metric_column %in% colnames(df)) {
      warning(glue::glue("Skipping {value_type} raster: Column '{metric_column}' not found in dataframe."))
      next
    }
    
    # Ensure the dataframe isn't empty before classification
    if (nrow(df) == 0) {
      warning(glue::glue("Skipping {value_type} raster: No data available after filtering. Creating an empty raster."))
      raster_reclassified <- terra::ifel(!is.na(raster), NA, NA)  # Raster filled with NA
    } else {
      # Apply rounding if needed
      df <- df |>
        dplyr::mutate(!!metric_column := if (!is.na(perc_digits)) round(.data[[metric_column]], perc_digits) else .data[[metric_column]])
      
      # Create the classification matrix
      reclass_matrix <- df |>
        dplyr::select(all_of(cat_base_column_name), all_of(metric_column)) |>
        as.matrix()
      
      # Ensure the classification matrix has at least one row before calling terra::classify()
      if (nrow(reclass_matrix) > 0) {
        raster_reclassified <- terra::classify(raster, reclass_matrix)
      } else {
        warning(glue::glue("Skipping {value_type} raster: Classification matrix is empty. Creating an empty raster."))
        raster_reclassified <- terra::ifel(!is.na(raster), NA, NA)
      }
    }
    
    if ("MEMORY" %in% raster_return) {
      raster_list[[value_type]] <- raster_reclassified
    }
    if ("WRITE" %in% raster_return) {
      output_file <- paste0(base_path, "_", tolower(value_type), ".tif")
      terra::writeRaster(raster_reclassified, output_file, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE"))
      raster_files[[value_type]] <- output_file
    }
  }
  
  result <- list()
  if ("MEMORY" %in% raster_return) result$raster_list <- raster_list
  if ("WRITE" %in% raster_return) result$raster_files <- raster_files
  
  return(result)
}






# save_rasters <- function(raster, df, type, out_dir, run_name, 
#                          out_rast_values, perc_digits, raster_return, 
#                          cat_base_column_name, new_sub_dir = FALSE) {  
#   
#   # Ensure output directory structure
#   if (new_sub_dir & "WRITE" %in% raster_return) {
#     out_dir <- file.path(out_dir, run_name)
#     dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#   }
#   
#   base_path <- file.path(out_dir, paste0(run_name, "_", type))
#   raster_list <- list()
#   
#   for (value_type in out_rast_values) {
#     
#     # If RAW, directly return/write the clipped raster without reclassification
#     if (value_type == "RAW") {
#       raster_raw <- raster  # Just retain the original values
#       
#       if ("MEMORY" %in% raster_return) {
#         raster_list[[value_type]] <- raster_raw
#       }
#       if ("WRITE" %in% raster_return) {
#         output_file <- paste0(base_path, "_raw.tif")
#         terra::writeRaster(raster_raw, output_file, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE"))
#       }
#       
#       next  # Skip to the next raster type
#     }
#     
#     # Determine the appropriate metric column
#     metric_column <- switch(value_type,
#                             "PERC_COVER_AOI" = "aoi_perc",
#                             "PERC_COVER_REGION" = "region_perc",
#                             stop(glue::glue("Invalid out_rast_values option: {value_type}")))
#     
#     # Ensure the necessary columns exist
#     if (!cat_base_column_name %in% colnames(df)) {
#       warning(glue::glue("Skipping {value_type} raster: Column '{cat_base_column_name}' not found in dataframe."))
#       next
#     }
#     
#     if (!metric_column %in% colnames(df)) {
#       warning(glue::glue("Skipping {value_type} raster: Column '{metric_column}' not found in dataframe."))
#       next
#     }
#     
#     # Ensure the dataframe isn't empty before classification
#     if (nrow(df) == 0) {
#       warning(glue::glue("Skipping {value_type} raster: No data available after filtering. Creating an empty raster."))
#       raster_reclassified <- terra::ifel(!is.na(raster), NA, NA)  # Raster filled with NA
#     } else {
#       # Apply rounding if needed
#       df <- df |>
#         dplyr::mutate(!!metric_column := if (!is.na(perc_digits)) round(.data[[metric_column]], perc_digits) else .data[[metric_column]])
#       
#       # Create the classification matrix
#       reclass_matrix <- df |>
#         dplyr::select(all_of(cat_base_column_name), all_of(metric_column)) |> 
#         as.matrix()
#       
#       # Ensure the classification matrix has at least one row before calling terra::classify()
#       if (nrow(reclass_matrix) > 0) {
#         raster_reclassified <- terra::classify(raster, reclass_matrix)
#       } else {
#         warning(glue::glue("Skipping {value_type} raster: Classification matrix is empty. Creating an empty raster."))
#         raster_reclassified <- terra::ifel(!is.na(raster), NA, NA)
#       }
#     }
#     
#     if ("MEMORY" %in% raster_return) {
#       raster_list[[value_type]] <- raster_reclassified
#     }
#     if ("WRITE" %in% raster_return) {
#       output_file <- paste0(base_path, "_", tolower(value_type), ".tif")
#       terra::writeRaster(raster_reclassified, output_file, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE"))
#     }
#   }
#   
#   if ("MEMORY" %in% raster_return) {
#     return(raster_list)
#   }
#   if ("WRITE" %in% raster_return) {
#     return(invisible(NULL))
#   }
# }






## Function to analyze landcover comparisons between AOI and a larger region.
# It will return a dataframe with columns for region_cover and aoi_cover
# both raw and in percentage of the total region and percentage of the aoi.
# The function will join the raster categorical data to the frequencies,
# and group the data and summarize it if desired based on a column of interest 
#aoiCover - land cover raster for the smaller region of interest
#regionCover - land cover raster for a region
# group: whether to group the data or not (TRUE, FALSE) - note that if grouped, all other columns from cats will be dropped, and the output will no longer have the associated raster values
analyze_categorical_cover <- function(aoi_raster,
                                      larger_region_raster,
                                      raster_cat_df,
                                      cat_base_column_name,
                                      group = FALSE,
                                      cat_group_column_name = NA) {
  
  # THIS IS NEW #
  #Ensure that raster categories and active category are set correctly
  levels(aoi_raster) <- raster_cat_df
  terra::activeCat(aoi_raster) <- cat_base_column_name
  levels(larger_region_raster) <- raster_cat_df
  terra::activeCat(larger_region_raster) <- cat_base_column_name
  # # #
  
  #Get frequencies & ensure that freq tables are clean
  print("Getting frequencies")
  aoi_freq <- terra::freq(aoi_raster) |>
    dplyr::select(-layer) |>
    condense_freq_groups()
  region_freq <- terra::freq(larger_region_raster) |>
    dplyr::select(-layer) |>
    condense_freq_groups()
  
  # Join together
  all_freqs <- region_freq |>
    dplyr::full_join(aoi_freq, by = c("value")) |>
    dplyr::mutate_if(is.numeric, coalesce, 0) |> #remove NAs
    dplyr::filter(value != 0) # remove background
  
  
  cat_base_column_sym <- rlang::sym(cat_base_column_name)
  
  names(all_freqs) <- c(cat_base_column_name,
                        "region_count",
                        "aoi_count")
  
  # Join raster cat codes
  
  #Ensure both columns are numeric
  raster_cat_df <- raster_cat_df |>
    dplyr::mutate({{cat_base_column_sym}} := as.numeric({{cat_base_column_sym}}))
  all_freqs <- all_freqs |>
    dplyr::mutate({{cat_base_column_sym}} := as.numeric({{cat_base_column_sym}}))
  
  all_freqs <- all_freqs |>
    dplyr::left_join(raster_cat_df, by = cat_base_column_name)
  
  
  # Perform grouping summarization on group of choice if desired
  if(group) {
    cat_group_column_sym <- rlang::sym(cat_group_column_name)
    all_freqs <- all_freqs |>
      dplyr::group_by({{cat_group_column_sym}}) |>
      dplyr::summarise(region_count = sum(region_count),
                       aoi_count = sum(aoi_count)) |>
      dplyr::ungroup()
  }
  
  
  # Add percentages and differences for analysis
  print("Calculating percentages")
  all_freqs <- all_freqs |>
    dplyr::mutate(region_perc = 100 * (region_count / sum(all_freqs$region_count)),
                  aoi_perc = 100 * (aoi_count / sum(all_freqs$aoi_count)),
                  diff_in_perc = region_perc - aoi_perc,
                  diff_in_num = region_count - aoi_perc)
  
  # Arrange and return
  all_freqs <- all_freqs |>
    dplyr::arrange(dplyr::desc(diff_in_perc))
  
  return(all_freqs)
}



# A function to create the new raster; all values in df will now be NA
# TIF must contain a column titled "VALUE"
keep_tif_values_in_df <- function(raster, df) {
  # values_not_in_df <- !terra::values(raster) %in% df$VALUE
  # new_raster <- raster
  # values(new_raster)[values_not_in_df] <- NA
  # 
  # rm(values_not_in_df)
  
  new_raster <- terra::ifel(raster %in% df$VALUE, raster, NA)
  
  #Maintain old cats
  c <- terra::cats(raster)[[1]]
  levels(new_raster) <- c
  
  return(new_raster)
}



#' Represent Graph of Coverage Comparison
#'
#' This function generates a scatter plot comparing regional coverage percentages
#' (`region_perc`) and area of interest (AOI) coverage percentages (`aoi_perc`).
#' It also calculates and displays absolute differences, including mean and median
#' values across different coverage thresholds.
#'
#' @param coverage_df A data frame containing columns `region_perc`, `aoi_perc`, 
#'   and `diff_in_perc`, where `region_perc` and `aoi_perc` represent coverage
#'   percentages, and `diff_in_perc` represents their difference.
#' @param manual_lim Optional numeric value to manually set the axis limits. If `NA` (default),
#'   the global maximum of `region_perc` and `aoi_perc` is used.
#' @param log Logical, whether to apply a pseudo-logarithmic scale to the axes. Default is `FALSE`.
#'
#' @return A `ggplot2` object representing the scatter plot.
#'
#' @details
#' - The plot includes a 45-degree reference line (`y = x`) in red to indicate perfect agreement.
#' - A caption summarizes mean and median absolute differences for all data points, as well
#'   as subsets with `region_perc` greater than 1% and 5%.
#' - If `log = TRUE`, a pseudo-logarithmic scale is applied to both axes.
#'
#' @importFrom dplyr filter pull
#' @importFrom ggplot2 ggplot geom_point geom_abline theme_minimal labs
#' @importFrom scales pseudo_log_trans
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' library(dplyr)
#' 
#' # Example data frame
#' coverage_df <- data.frame(
#'   region_perc = runif(100, 0, 100),
#'   aoi_perc = runif(100, 0, 100),
#'   diff_in_perc = runif(100, -10, 10)
#' )
#'
#' # Default plot
#' represent_graph(coverage_df)
#'
#' # Log scale plot
#' represent_graph(coverage_df, log = TRUE)
#'
#' # Custom axis limit
#' represent_graph(coverage_df, manual_lim = 50)
#' }
#'
#' @export
represent_graph <- function(coverage_df,
                            manual_lim = NA,
                            log = FALSE) {
  if (is.na(manual_lim)) {
    global_max <- max(coverage_df$region_perc, coverage_df$aoi_perc)
  } else {
    global_max <- manual_lim
  }
  
  # Calculate perpendicular residuals to the 45-degree line
  coverage_df$diff_in_perc_abs <- abs(coverage_df$diff_in_perc)
  mean_diff <- mean(coverage_df$diff_in_perc_abs, na.rm = TRUE)
  median_diff <- median(coverage_df$diff_in_perc_abs, na.rm = TRUE)
  
  mean_diff_over_1_perc <- mean(coverage_df |>
                                  dplyr::filter(region_perc > 1) |>
                                  dplyr::pull(diff_in_perc_abs), na.rm = TRUE)
  median_diff_over_1_perc <- median(coverage_df |>
                                      dplyr::filter(region_perc > 1) |>
                                      dplyr::pull(diff_in_perc_abs), na.rm = TRUE)
  
  mean_diff_over_5_perc <- mean(coverage_df |>
                                  dplyr::filter(region_perc > 5) |>
                                  dplyr::pull(diff_in_perc_abs), na.rm = TRUE)
  median_diff_over_5_perc <- median(coverage_df |>
                                      dplyr::filter(region_perc > 5) |>
                                      dplyr::pull(diff_in_perc_abs), na.rm = TRUE)
  
  # Construct the caption text
  caption_text <- sprintf("All classes: Mean Abs. Difference: %.2f | Median Abs. Difference: %.2f\nAll classes with >1%% regional coverage: Mean Abs. Difference: %.2f | Median Abs. Difference: %.2f\nAll classes with >5%% regional coverage: Mean Abs. Difference: %.2f | Median Abs. Difference: %.2f",
                          mean_diff, median_diff, 
                          mean_diff_over_1_perc, median_diff_over_1_perc,
                          mean_diff_over_5_perc, median_diff_over_5_perc)
  
  p <- ggplot(data = coverage_df) +
    geom_point(aes(x = region_perc, y = aoi_perc)) +
    geom_abline(slope = 1, intercept = 0, color = "red", lty = 2) +
    theme_minimal() +
    labs(caption = caption_text)
  
  if (log) {
    p <- p +
      scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                         limits = c(0, global_max)) +
      scale_x_continuous(trans = scales::pseudo_log_trans(base = 10),
                         limits = c(0, global_max)) +
      xlab("Regional Coverage Percentage (log scale)") +
      ylab("Area of Interest Coverage Percentage (log scale)")
  } else {
    p <- p +
      scale_y_continuous(limits = c(0, global_max)) +
      scale_x_continuous(limits = c(0, global_max)) +
      xlab("Regional Coverage Percentage") +
      ylab("Area of Interest Coverage Percentage")
  }
  
  # if(bivariate_background) {
  #
  #   p <- p +
  #     annotate("rect", xmin=10, xmax=15, ymin=0, ymax=Inf, alpha=0.1, fill="gold")
  #
  # }
  
  return(p)
}


# Bivariate raster mapping ----


#' Convert Hexadecimal Color to RGB
#'
#' This function converts a 6-character hexadecimal color code into its
#' corresponding RGB values.
#'
#' @param hex A string representing a 6-character hex color code prefixed with `#`
#'   (e.g., `"#FF5733"`).
#'
#' @return A named numeric vector with RGB values (`R`, `G`, `B`).
#'
#' @examples
#' hex_to_rgb("#FF5733")
#' hex_to_rgb("#00AABB")
#'
#' @export
hex_to_rgb <- function(hex) {
  # Ensure the input is a valid hex code with a leading "#"
  if (!grepl("^#([A-Fa-f0-9]{6})$", hex)) {
    stop("Invalid hex color format. Use a 6-character hex code, e.g., '#FF5733'.")
  }
  
  # Extract RGB components
  r <- strtoi(substr(hex, 2, 3), base = 16)
  g <- strtoi(substr(hex, 4, 5), base = 16)
  b <- strtoi(substr(hex, 6, 7), base = 16)
  
  # Check for NA values
  if (any(is.na(c(r, g, b)))) stop("Failed to convert hex to RGB.")
  
  return(c(R = r, G = g, B = b))
}

#' Create a Bivariate Palette Using a Color Vector
#'
#' Generates a biscale-compatible palette from a set of nine hex color codes.
#'
#' @param p_vec A character vector of exactly 9 hex color codes.
#' @param flip A logical value. If `TRUE`, the palette is reversed.
#'
#' @return A list containing:
#'   - `biscale_pal`: A named character vector of hex colors mapped to bivariate classes.
#'   - `rgb_df`: A data frame with class identifiers and corresponding RGB values.
#'
#' @details The function creates a bivariate color palette compatible with biscale mapping.
#'   The `pals` package provides various palettes that can be used to define `p_vec`.
#'   See available palettes here: \url{https://cran.r-project.org/web/packages/pals/vignettes/pals_examples.html}.
#'
#' @examples
#' p_vec <- c("#f7fbff", "#deebf7", "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#08519c", "#08306b")
#' pals_bi_3(p_vec)
#' pals_bi_3(p_vec, flip = TRUE)
#'
#' @export
pals_bi_3 <- function(p_vec, flip = FALSE) {
  if (length(p_vec) != 9) {
    stop("p_vec must contain exactly 9 hex color codes.")
  }
  
  if (flip) {
    p_vec <- rev(p_vec)
  }
  
  # Create biscale-compatible palette
  biscale_pal <- setNames(p_vec, c("1-1", "2-1", "3-1", "1-2", "2-2", "3-2", "1-3", "2-3", "3-3"))
  
  # Create dataframe
  class <- c(11, 12, 13, 21, 22, 23, 31, 32, 33)
  rgb_list <- lapply(p_vec, hex_to_rgb)
  rgb_df <- data.frame(class, do.call(rbind, rgb_list), stringsAsFactors = FALSE)
  
  return(list(biscale_pal = biscale_pal, rgb_df = rgb_df))
}

#' Normalize Values to a [0,1] Scale
#'
#' This function rescales numeric values to a 0-to-1 range.
#'
#' @param x A numeric vector or raster object to be normalized.
#'
#' @return A numeric vector or raster with values normalized between 0 and 1.
#'
#' @examples
#' normalize(c(10, 20, 30, 40))
#' normalize(c(5, 15, NA, 25, 35))
#'
#' @export
normalize <- function(x) (x - min(x[], na.rm = TRUE)) / (max(x[], na.rm = TRUE) - min(x[], na.rm = TRUE))

#' Normalize Raster Data Using Global Min-Max Scaling
#'
#' This function normalizes raster values from two datasets using a common
#' global minimum and maximum.
#'
#' @param x A `SpatRaster` object from the `terra` package.
#' @param y A second `SpatRaster` object to be used for determining the global range.
#'
#' @return A `SpatRaster` object where `x` has been normalized using the min/max
#'   values from both `x` and `y`.
#'
#' @details This function is useful when normalizing raster datasets that should
#'   be compared on the same scale. It ensures that both rasters share the same
#'   min-max normalization bounds.
#'
#' @examples
#' # Example assumes `x` and `y` are SpatRaster objects
#' \dontrun{
#' library(terra)
#' r1 <- rast(matrix(runif(100, 10, 50), 10, 10))
#' r2 <- rast(matrix(runif(100, 20, 60), 10, 10))
#' normalized_raster <- bi_dimensional_normalize_raster(r1, r2)
#' }
#'
#' @export
bi_dimensional_normalize_raster <- function(x, y) {
  global_min <- min(terra::minmax(x)[1], terra::minmax(y)[1])
  global_max <- max(terra::minmax(x)[2], terra::minmax(y)[2])
  
  x_norm <- (x - global_min) / (global_max - global_min)
  return(x_norm)                                          
}



#' Bivariate Raster Visualization
#'
#' This function generates a bivariate raster visualization by normalizing two raster datasets, classifying them
#' into categories, and mapping them to a bivariate color scale.
#'
#' @param x A `SpatRaster` object from the `terra` package representing the first variable.
#' @param y A `SpatRaster` object from the `terra` package representing the second variable.
#' @param bi_normal Logical; if `TRUE`, applies bivariate normalization using `bi_dimensional_normalize_raster()`;
#'   otherwise, applies standard normalization using `normalize()`.
#' @param pals_pal A character vector of exactly 9 hex color codes defining the bivariate color palette.
#' @param flip Logical; if `TRUE`, the color palette is reversed.
#' @param x_nm Character; label for the x-axis in the legend.
#' @param y_nm Character; label for the y-axis in the legend.
#'
#' @return A list containing:
#'   - `biv_plot`: A `tmap` object displaying the bivariate raster visualization.
#'   - `legend`: A `ggplot` object representing the legend for the bivariate color scale.
#'
#' @details The function:
#'   1. Normalizes the raster values.
#'   2. Categorizes them into three levels per variable (total 3x3 = 9 classes).
#'   3. Maps classified values to a predefined bivariate color palette.
#'   4. Generates an RGB raster for visualization.
#'   5. Uses `tmap::tm_rgb()` for rendering.
#'   6. Creates a legend using `biscale::bi_legend()`.
#'
#' @examples
#' \dontrun{
#' library(terra)
#' r1 <- rast(matrix(runif(100, 10, 50), 10, 10))
#' r2 <- rast(matrix(runif(100, 20, 60), 10, 10))
#' result <- bivariate_raster_viz_3(r1, r2, x_nm = "Variable X", y_nm = "Variable Y")
#' print(result$biv_plot)
#' print(result$legend)
#' }
#'
#' @export
bivariate_raster_viz_3 <- function(x,
                                   y,
                                   bi_normal = TRUE,
                                   pals_pal = pals::brewer.seqseq2(n = 9),
                                   flip = FALSE,
                                   x_nm,
                                   y_nm) {
  
  # ---- Parameter Checks ----
  
  if (!inherits(x, "SpatRaster")) stop("Error: 'x' must be a SpatRaster object from the terra package.")
  if (!inherits(y, "SpatRaster")) stop("Error: 'y' must be a SpatRaster object from the terra package.")
  if (!is.logical(bi_normal) || length(bi_normal) != 1) stop("Error: 'bi_normal' must be a logical value (TRUE/FALSE).")
  if (!is.character(pals_pal) || length(pals_pal) != 9) stop("Error: 'pals_pal' must be a character vector of exactly 9 hex color codes.")
  if (!is.logical(flip) || length(flip) != 1) stop("Error: 'flip' must be a logical value (TRUE/FALSE).")
  if (!is.character(x_nm) || length(x_nm) != 1) stop("Error: 'x_nm' must be a single character string.")
  if (!is.character(y_nm) || length(y_nm) != 1) stop("Error: 'y_nm' must be a single character string.")
  
  
  # ---- Spatial Consistency Check ----
  if (!terra::compareGeom(x, y)) stop("Error: 'x' and 'y' must have the same geometry. Check CRS, extent, resolution, and origin.")
  
  
  # ---- Normalize Data ----
  if (bi_normal) {
    x_norm <- bi_dimensional_normalize_raster(x, y)
    y_norm <- bi_dimensional_normalize_raster(y, x)
  } else {
    x_norm <- normalize(x)
    y_norm <- normalize(y)
  }
  
  # ---- Define Classification Breaks ----
  breaks <- seq(0, 1, length.out = 4)  # 3 categories (0-0.33, 0.33-0.67, 0.67-1)
  
  # ---- Classify Raster Data ----
  x_cat <- classify(x_norm, breaks, include.lowest = TRUE)
  y_cat <- classify(y_norm, breaks, include.lowest = TRUE)
  
  # ---- Create Bivariate Classification Raster ----
  bivariate_raster <- (y_cat + 1) * 10 + (x_cat + 1)  # Unique ID for each bivariate class
  
  # ---- Generate Bivariate Palette ----
  biv_p <- pals_bi_3(pals_pal, flip = flip)
  bivariate_palette <- biv_p$rgb_df
  
  # ---- Convert Classification Raster to RGB ----
  R <- classify(bivariate_raster, bivariate_palette[, c("class", "R")])
  G <- classify(bivariate_raster, bivariate_palette[, c("class", "G")])
  B <- classify(bivariate_raster, bivariate_palette[, c("class", "B")])
  
  # ---- Stack RGB Layers into a Single Raster ----
  bivariate_rgb <- c(R, G, B)
  names(bivariate_rgb) <- c("R", "G", "B")
  
  # ---- Plot the Bivariate Raster ----
  biv_plot <- tmap::tm_shape(bivariate_rgb) +
    tmap::tm_rgb(r = 1, g = 2, b = 3, tm_scale_rgb(max_color_value = 255))
  
  # ---- Create the Legend ----
  bi_leg <- biscale::bi_legend(
    pal = biv_p$biscale_pal,
    dim = 3,
    xlab = x_nm,
    ylab = y_nm,
    pad_width = 3
  )
  
  # ---- Return Plot and Legend ----
  return(list(
    biv_plot = biv_plot,
    legend = bi_leg
  ))
}


# Generating plots ----


# A function to randomly generate aop locations for Macrosystems
# PARAMETERS
# aop :: the aop area(s) of interest, as an sf object
# lc :: the landcover file to use for selections, as a SpatRaster
# numPlotsEach :: the number of plots to select for each landcover type present in the lc file
# name :: the human-readable name of the function run, which will be added to file paths (e.g. "Bend")
# lcName :: the human-readable name of the landcover tif (e.g. evt5window)
# outDir :: the directory to output data to
# sma :: the surface management agency. "USFS" or other string
generate_aop_plots <- function(aop,
                               lc,
                               numPlotsEach,
                               name,
                               lcName,
                               outDir,
                               sma = "USFS",
                               distFromRoadsMax) {
  
  epsg <- "EPSG:5070" #Albers equal area
  
  
  #Create operating area
  reasonableAreasOfAccess <- sf::st_read(here::here("data/manual/reasonable_areas_of_access.gpkg")) |>  #these are manually created areas of access. For YELL it is from create_yell_neon_aoi.R
    sf::st_transform(epsg) |>
    dplyr::mutate(grp = 1) |>
    dplyr::group_by(grp) |>
    dplyr::summarise(geom = st_union(geom)) |>
    dplyr::ungroup()
  
  genOpArea <- aop |>
    sf::st_transform(epsg) |>
    dplyr::mutate(group = 1) |>
    dplyr::group_by(group) |>
    dplyr::summarise(geometry = st_union(geometry)) |>
    dplyr::ungroup() |>
    sf::st_buffer(-50) |>
    sf::st_intersection(reasonableAreasOfAccess) |>
    dplyr::select(-group, -grp)
  
  
  #Access slope
  slope <- access_us_slope(slopeF = here::here('data/raw', 'usa_slope.tif')) #DEM Data
  # Clip slope
  thisSlope <- slope |>
    crop_careful_universal(aop, mask = TRUE, verbose = FALSE) |>
    terra::project(epsg)
  
  #Access roads
  roads <- access_osm_roads(aoi = genOpArea)
  
  
  #Buffer roads
  bigRoadBuff <- roads |> sf::st_buffer(distFromRoadsMax)
  
  #Surface Management Agency (SMA)
  
  if(sma == "USFS") {
    usfs <- access_us_sma(dir_raw,
                          layer = "SurfaceMgtAgy_USFS") |>
      sf::st_transform(epsg)
    
    thisUsfs <- usfs |>
      sf::st_filter(genOpArea) |>
      dplyr::mutate(group = 1) |>
      dplyr::group_by(group) |>
      dplyr::summarise(SHAPE = st_union(SHAPE)) |>
      dplyr::ungroup() |>
      sf::st_buffer(-400) |>
      sf::st_intersection(genOpArea)
  }
  
  
  
  #Clip land cover data
  evtHere <- lc |>
    crop_careful_universal(vector = genOpArea,
                           mask = TRUE,
                           verbose = FALSE) |>
    terra::project(epsg)
  
  
  #Slope mask
  thisSlopeResample <- thisSlope |>
    terra::crop(evtHere,
                snap = "near",
                extend = TRUE) |>
    terra::resample(evtHere,
                    method = "near") |>
    terra::project(epsg)
  slopeMask <- thisSlopeResample > 15
  
  
  #Apply all masks and generate sampling options raster
  evtHereMasked <- evtHere |>
    terra::mask(mask = slopeMask,
                maskvalues = TRUE) |>
    terra::mask(bigRoadBuff)
  
  if(sma == "USFS") {
    evtHereMasked <- evtHereMasked |>
      terra::mask(thisUsfs)
  }
  
  evtCsv <- access_landfire_evt_conus_2022_csv() |>
    dplyr::mutate(VALUE = as.integer(VALUE))
  
  # Sample raster
  set.seed(1)
  potentialPlots <- raster::sampleStratified(raster::raster(evtHereMasked), size = numPlotsEach, na.rm = TRUE, sp = TRUE) |>
    sf::st_as_sf() |>
    dplyr::mutate(VALUE = EVT_NAME) |>
    dplyr::select(-EVT_NAME) |>
    dplyr::left_join(evtCsv, by = join_by(VALUE == VALUE)) |>
    dplyr::filter(EVT_GP_N != "Open Water" & !grepl("developed", EVT_PHYS, ignore.case = TRUE) & EVT_GP_N != "Quarries-Strip Mines-Gravel Pits-Well and Wind Pads") |>
    dplyr::mutate(plotID = paste0(name, "_", row_number())) |>
    dplyr::mutate(sampleRaster = lcName)
  
  potentialPlotLocationsWGS <- potentialPlots |>
    sf::st_transform("EPSG:4326") |>
    sf::st_coordinates() %>%
    cbind(potentialPlots, .) |>
    dplyr::rename(lat = Y, long = X)
  
  potentialPlotLocationsClean <- potentialPlotLocationsWGS |>
    dplyr::select(plotID, EVT_NAME, lat, long, geometry) |>
    dplyr::rename(landcover = EVT_NAME)
  
  
  #sf::st_write(potentialPlotLocationsWGS, here::here(outDir, glue::glue('aop_{name}_points.gpkg')), append = FALSE)
  sf::st_write(potentialPlotLocationsClean, here::here(outDir, glue::glue('aop_{name}_{lcName}_points.gpkg')), append = FALSE)
  write_csv(potentialPlotLocationsClean |> sf::st_drop_geometry(), here::here(outDir, glue::glue('aop_{name}_{lcName}_points.csv')))
  st_write_shp(shp = potentialPlotLocationsClean,
               location = outDir,
               filename = glue::glue("aop_{name}_{lcName}_points"),
               zip_only = TRUE)
}




# A function to create the general operating area for plot selection
# PARAMETERS
# ranger :: USFS ranger district as an sf spatial object
# ecoRegion :: EPA ecoregion as an sf spatial object
create_operating_area <- function(ranger, ecoRegion) {
  
  epsg <- "EPSG:5070" #Albers equal area
  
  ranger <- ranger |>
    sf::st_transform(epsg)
  
  ecoRegion <- ecoRegion |>
    sf::st_transform(epsg)
  
  reasonableAreasOfAccess <- sf::st_read(here::here("data/manual/reasonable_areas_of_access.gpkg")) |>  #these are manually created areas of access. For YELL it is from create_yell_neon_aoi.R
    sf::st_transform(epsg) |>
    dplyr::mutate(grp = 1) |>
    dplyr::group_by(grp) |>
    dplyr::summarise(geom = st_union(geom)) |>
    dplyr::ungroup()
  
  #Wilderness & WSA
  wildernessAll <- access_us_wilderness(dest_path = here::here('data', 'raw', 'wild.gpkg')) |> 
    sf::st_transform(epsg) |>
    sf::st_buffer(400) #USFS standard is no drones w/in 300m of wilderness boundaries + 70, same reasoning as road buffer; add an extra 30m to be safe = 400
  
  wilderness <- wildernessAll |>
    sf::st_filter(ranger) |>
    dplyr::mutate(group = 1) |>
    dplyr::group_by(group) |>
    dplyr::summarise(geometry = st_union(geometry)) |>
    dplyr::ungroup()
  
  wsaAll <- access_us_wilderness_study_areas(dest_path = here::here('data', 'raw', 'wsa.gpkg')) |>
    sf::st_transform(epsg) |>
    sf::st_buffer(400) #USFS standard is no drones w/in 300m of wilderness boundaries + 70, same reasoning as road buffer; add an extra 30m to be safe = 400
  
  wsa <- wsaAll |>
    sf::st_filter(ranger) |>
    dplyr::mutate(group = 1) |>
    dplyr::group_by(group) |>
    dplyr::summarise(geometry = st_union(geometry)) |>
    dplyr::ungroup()
  
  genOpArea <- ecoRegion |>
    sf::st_buffer(-50) |>
    sf::st_intersection(ranger) |>
    dplyr::mutate(group = 1) |>
    dplyr::group_by(group) |>
    dplyr::summarise(geometry = st_union(geometry)) |>
    dplyr::ungroup() |>
    sf::st_intersection(reasonableAreasOfAccess)
  
  if(nrow(wilderness) > 0) {
    genOpArea <- genOpArea |>
      sf::st_difference(wilderness)
  }
  if(nrow(wsa) > 0) {
    genOpArea <- genOpArea |>
      sf::st_difference(wsa)
  }
  
  return(genOpArea)
}



# Generating AOP tile polygons ----

# Helper function to determine the appropriate NAD83 UTM zone EPSG code for a given polygon
find_nad83_utm_epsg <- function(polygon) {
  
  # Ensure the polygon is in WGS84 (EPSG:4326) for accurate longitude calculation
  polygonWgs84 <- sf::st_transform(polygon, crs = 4326)
  
  # Calculate the centroid
  centroid <- sf::st_centroid(polygonWgs84)
  centroidCoords <- sf::st_coordinates(centroid)
  
  # Get the longitude of the centroid
  longitude <- centroidCoords[1, "X"]
  
  # Calculate the UTM zone
  utmZone <- (floor((longitude + 180) / 6) %% 60) + 1
  
  # Determine the EPSG code for NAD83 UTM zone
  epsgCode <- 26900 + utmZone
  
  return(epsgCode)
}


# A function to generate a gridded set of polygons representing AOP image tiles for a NEON AOP footprint (the given aoi)
generate_aop_tile_polygons <- function(aoi) {
  
  # Find the NAD83 UTM EPSG code for the AOI
  epsgCode <- find_nad83_utm_epsg(aoi)
  
  # Transform the AOI to the appropriate UTM CRS
  aoiUtm <- sf::st_transform(aoi, crs = epsgCode)
  
  # Get the bounding box of the AOI
  bbox <- sf::st_bbox(aoiUtm)
  
  # Define the range of coordinates for the grid
  xCoords <- seq(floor(bbox["xmin"] / 1000) * 1000, ceiling(bbox["xmax"] / 1000) * 1000, by = 1000)
  yCoords <- seq(floor(bbox["ymin"] / 1000) * 1000, ceiling(bbox["ymax"] / 1000) * 1000, by = 1000)
  
  # Create the grid of polygons
  polygons <- list()
  idIndex <- 1
  
  for (x in xCoords) {
    for (y in yCoords) {
      # Create a 1km by 1km polygon
      polygon <- sf::st_polygon(list(matrix(c(x, y, 
                                          x + 1000, y, 
                                          x + 1000, y + 1000, 
                                          x, y + 1000, 
                                          x, y), 
                                        ncol = 2, byrow = TRUE)))
      
      # Create an sf object for the polygon with an id
      polygons[[idIndex]] <- sf::st_sf(id = paste0(x, "_", y), geometry = sf::st_sfc(polygon, crs = epsgCode))
      idIndex <- idIndex + 1
    }
  }
  
  # Combine all polygons into a single sf object using rbind
  grid <- do.call(rbind, polygons)
  
  # Clip the grid to the AOI
  gridClipped <- sf::st_intersection(grid, aoiUtm)
  
  return(gridClipped)
}


# Operational function to create and write out a polygon set with the names of all image tiles within an AOP footprint
# Example use
# 
# yellTest <- generate_aop_image_grid(site = "YELL",
#                                     year = 2023,
#                                     visit = 5)
generate_aop_image_grid <- function(site, year, visit) {
  
  aop_all <- access_neon_aop_flight_box_data()
  
  aop <- aop_all |>
    dplyr::filter(siteID == site)
  
  site_grid <- generate_aop_tile_polygons(aop) |>
    dplyr::mutate(tileName = paste0(year, "_", site, "_", visit, "_", id, "_image"))
  
  return(site_grid)
}




#' Find NAD83 UTM EPSG Code
#'
#' This helper function calculates the appropriate NAD83 UTM EPSG code for a given polygon based 
#' on its centroid's longitude. The polygon is first transformed to WGS84 (EPSG:4326) for accurate 
#' calculation.
#'
#' @param polygon An sf object representing a polygon. The coordinate reference system (CRS) is assumed to be set.
#'
#' @return An integer representing the EPSG code for the corresponding NAD83 UTM zone.
#'
#' @importFrom sf st_transform st_centroid st_coordinates
#' 
#' @export
find_nad83_utm_epsg <- function(polygon) {
  # Ensure the polygon is in WGS84 (EPSG:4326) for accurate longitude calculation
  polygonWgs84 <- sf::st_transform(polygon, crs = 4326)
  
  # Calculate the centroid
  centroid <- sf::st_centroid(polygonWgs84)
  centroidCoords <- sf::st_coordinates(centroid)
  
  # Get the longitude of the centroid
  longitude <- centroidCoords[1, "X"]
  
  # Calculate the UTM zone
  utmZone <- (floor((longitude + 180) / 6) %% 60) + 1
  
  # Determine the EPSG code for NAD83 UTM zone
  epsgCode <- 26900 + utmZone
  
  return(epsgCode)
}


#' Generate AOP Tile Polygons
#'
#' Generates a set of 1km x 1km polygons representing AOP image tiles within the bounding box 
#' of a given Area of Interest (AOI). The AOI is transformed to the appropriate UTM CRS based 
#' on its centroid before the grid of polygons is created.
#'
#' @param aoi An sf object representing the AOI (area of interest) polygon.
#'
#' @return An sf object containing the gridded polygons clipped to the AOI.
#'
#' @importFrom sf st_transform st_bbox st_polygon st_sfc st_sf st_intersection
#' @importFrom dplyr mutate
#' 
#' @export
generate_aop_tile_polygons <- function(aoi) {
  
  # Find the NAD83 UTM EPSG code for the AOI
  epsgCode <- find_nad83_utm_epsg(aoi)
  
  # Transform the AOI to the appropriate UTM CRS
  aoiUtm <- sf::st_transform(aoi, crs = epsgCode)
  
  # Get the bounding box of the AOI
  bbox <- sf::st_bbox(aoiUtm)
  
  # Define the range of coordinates for the grid
  xCoords <- seq(floor(bbox["xmin"] / 1000) * 1000, ceiling(bbox["xmax"] / 1000) * 1000, by = 1000)
  yCoords <- seq(floor(bbox["ymin"] / 1000) * 1000, ceiling(bbox["ymax"] / 1000) * 1000, by = 1000)
  
  # Create the grid of polygons
  polygons <- list()
  idIndex <- 1
  
  for (x in xCoords) {
    for (y in yCoords) {
      # Create a 1km by 1km polygon
      polygon <- sf::st_polygon(list(matrix(c(x, y, 
                                              x + 1000, y, 
                                              x + 1000, y + 1000, 
                                              x, y + 1000, 
                                              x, y), 
                                            ncol = 2, byrow = TRUE)))
      
      # Create an sf object for the polygon with an id
      polygons[[idIndex]] <- sf::st_sf(id = paste0(x, "_", y), geometry = sf::st_sfc(polygon, crs = epsgCode))
      idIndex <- idIndex + 1
    }
  }
  
  # Combine all polygons into a single sf object using rbind
  grid <- do.call(rbind, polygons)
  
  # Clip the grid to the AOI
  gridClipped <- sf::st_intersection(grid, aoiUtm)
  
  return(gridClipped)
}




#' Generate AOP Image Grid
#'
#' Generates and returns a gridded set of polygons representing AOP image tiles for a given NEON site, 
#' year, and visit number. The function retrieves the flight box data, filters it by site, and 
#' then calls \code{\link{generate_aop_tile_polygons}} to generate the grid. Each tile is named according 
#' to the input parameters.
#'
#' @param site A character string representing the NEON site identifier (e.g., "YELL").
#' @param year An integer specifying the year of the AOP data (e.g., 2023).
#' @param visit An integer representing the visit number (e.g., 5).
#'
#' @return An sf object representing the gridded AOP image tiles with names that include the year, site, 
#' visit number, and tile identifier.
#'
#' @importFrom dplyr filter mutate
#' 
#' @examples
#' \dontrun{
#' # Example: Generate AOP image grid for site "YELL", year 2023, and visit 5.
#' yell_image_grid <- generate_aop_image_grid(site = "YELL", year = 2023, visit = 5)
#' 
#' # View the result
#' print(yell_image_grid)
#' 
#' # Plot the generated grid
#' plot(sf::st_geometry(yell_image_grid))
#' }
#' @export
generate_aop_image_grid <- function(site, year, visit) {
  
  # Access the AOP flight box data
  aop_all <- access_neon_aop_flight_box_data()
  
  # Filter the data for the specified site
  aop <- aop_all |> 
    dplyr::filter(siteID == site)
  
  # Generate the tile polygons for the site
  site_grid <- generate_aop_tile_polygons(aop) |>
    dplyr::mutate(tileName = paste0(year, "_", site, "_", visit, "_", id, "_image"))
  
  return(site_grid)
}



# Generate plot boundaries ----

# Function to get the extent of pixels
get_pixel_extent <- function(row, col, raster) {
  cell_res <- terra::res(raster)  # Get cell resolution
  x_min <- terra::ext(raster)$xmin + cell_res[1] * (col - 1)
  y_max <- terra::ext(raster)$ymax - cell_res[2] * (row - 1)
  x_max <- x_min + cell_res[1]
  y_min <- y_max - cell_res[2]
  return(c(x_min, x_max, y_min, y_max))
}



# A function to create a 0&1 checkerboard SpatRaster that perfectly aligns with ARD pixels,
# providing coverage of the area over which plotPoints exists (buffered by an extra 200m)
create.landsat.checkerboard <- function(plotPoints) {
  
  #Get bounding boxes - both projected & in geographic coordinates
  bbox4326 <- sf::st_bbox(plotPoints |>
                            sf::st_buffer(200) |>
                            sf::st_transform("EPSG:4326"))
  
  bboxProj <- sf::st_bbox(plotPoints |>
                            sf::st_buffer(200))
  
  # Search the stac catalog
  stac("https://landsatlook.usgs.gov/stac-server") |>
    get_request()
  ## ###STACCatalog
  ## - id: earth-search-aws
  ## - description: A STAC API of public datasets on AWS
  ## - field(s): stac_version, type, id, title, description, links, conformsTo
  
  collection_formats()
  
  
  # Record start time
  a <- Sys.time()
  
  # Initialize STAC connection
  s = rstac::stac("https://landsatlook.usgs.gov/stac-server")
  
  
  # Search for landsat ARD images within specified bounding box and date range
  items = s |>
    rstac::stac_search(collections = "landsat-c2ard-sr",
                       bbox = c(bbox4326["xmin"],
                                bbox4326["ymin"],
                                bbox4326["xmax"],
                                bbox4326["ymax"]),
                       datetime = "2021-05-01T00:00:00Z/2021-06-01T00:00:00Z") |>
    post_request() |>
    items_fetch(progress = FALSE)
  
  
  
  
  
  # Print number of found items
  length(items$features)
  
  items
  
  
  # Prepare the assets for analysis
  library(gdalcubes)
  assets = c("qa_pixel")
  ard_collection = gdalcubes::stac_image_collection(items$features, asset_names = assets)
  
  b <- Sys.time()
  difftime(b, a)
  
  # Display the image collection
  ard_collection
  
  
  
  #Access the data
  
  # Record start time
  a <- Sys.time()
  
  # Define a specific view on the satellite image collection
  v = gdalcubes::cube_view(
    srs = epsg,
    dx = 30,
    dy = 30,
    dt = "P1M",
    aggregation = "median",
    resampling = "near",
    extent = list(
      t0 = "2021-05-01",
      t1 = "2021-06-01",
      left = bboxProj["xmin"],
      right = bboxProj["xmax"],
      top = bboxProj["ymax"],
      bottom = bboxProj["ymin"]
    )
  )
  
  b <- Sys.time()
  difftime(b, a)
  
  # Display the defined view
  v
  
  
  
  a <- Sys.time()
  
  x <- ard_collection |>
    raster_cube(v) |>
    write_tif() |>
    terra::rast()
  
  x
  
  b <- Sys.time()
  difftime(b, a)
  
  nrows <- nrow(x)
  ncols <- ncol(x)
  
  # Create a matrix with a checkerboard pattern
  checkerboard_matrix <- matrix(NA, nrow = nrows, ncol = ncols)
  for (i in 1:nrows) {
    for (j in 1:ncols) {
      checkerboard_matrix[i, j] <- ifelse((i + j) %% 2 == 0, 1, 0)
    }
  }
  
  # Convert the matrix to a raster
  filled_raster <- terra::rast(checkerboard_matrix, crs = terra::crs(x), extent = terra::ext(x))
  
  return(filled_raster)
  
}



# A function to create plot polygons from a raster template and the plot points
create.polygons.from.raster.and.points <- function(raster, plotPoints) {
  
  plotCoords <- sf::st_coordinates(plotPoints)
  
  cell_index <- terra::cellFromXY(raster, plotCoords)
  row <- terra::rowFromCell(raster, cell_index)
  col <- terra::colFromCell(raster, cell_index)
  
  
  # Find extents of neighboring pixels
  neighboring_extents <- list()
  for (i in seq_along(row)) {
    neighboring_extents[[i]] <- list()
    for (j in c(-1, 0, 1)) {
      for (k in c(-1, 0, 1)) {
        neighboring_row <- row[i] + j
        neighboring_col <- col[i] + k
        neighboring_extents[[i]][[length(neighboring_extents[[i]]) + 1]] <- get_pixel_extent(neighboring_row, neighboring_col, raster)
      }
    }
  }
  
  # Convert neighboring extents to a data frame
  neighbor_df <- data.frame()
  for (i in seq_along(neighboring_extents)) {
    for (extent in neighboring_extents[[i]]) {
      neighbor_df <- rbind(neighbor_df, c(i, extent))
    }
  }
  colnames(neighbor_df) <- c("ID", "xmin", "xmax", "ymin", "ymax")
  
  # Aggregate extents by ID and create plot polygons
  plot_extents <- neighbor_df %>%
    dplyr::group_by(ID) %>%
    dplyr::summarise(xmin = min(xmin), xmax = max(xmax), ymin = min(ymin), ymax = max(ymax)) %>%
    as.data.frame()
  
  # Create SpatialPolygons from the extents
  plot_polygons <- apply(plot_extents[2:5], 1, function(x) {
    terra::ext(x["xmin"], x["xmax"], x["ymin"], x["ymax"]) |> 
      terra::as.polygons() |> 
      st_as_sf()
  })
  
  # Combine the polygons into a single sf object
  plot_sf <- do.call(rbind, plot_polygons)
  st_crs(plot_sf) <- sf::st_crs(plotPoints)
  plotID <- plotPoints$plotID
  
  plot_sf <- cbind(plot_sf, plotID)
  plot_sf <- plot_sf |> dplyr::left_join(plotPoints |> sf::st_drop_geometry()) 
  
  return(plot_sf)
  
}


# Function to create a full set of macrosystems plots
## Main operational function ----
generate.and.write.macro.plots.and.dem <- function(plotPoints, areaName) {
  
  # Setup directory
  exportLocation = here::here('data', 'fieldwork', 'plots', areaName)
  if(!dir.exists(exportLocation)) {
    dir.create(exportLocation)
  }
  
  #Generate landsat ARD raster template
  ardTemplateFile <- here::here(exportLocation, paste0(areaName, '_ARDTemplate.tif'))
  if(!file.exists(ardTemplateFile)) {
    landsatARDTemplate <- create.landsat.checkerboard(plotPoints)
    #write ard template raster
    terra::writeRaster(landsatARDTemplate,
                       ardTemplateFile,
                       overwrite = TRUE,
                       gdal=c("COMPRESS=DEFLATE"))
  } else {
    landsatARDTemplate <- terra::rast(ardTemplateFile)
  }
  
  
  #Get DEM for area
  areaDEMFile <- here::here(exportLocation, paste0(areaName, '_DEM.tif'))
  
  if(!file.exists(areaDEMFile)) {
    areaDEM <- elevatr::get_elev_raster(landsatARDTemplate, z = 12) |> #https://cran.r-project.org/web/packages/elevatr/vignettes/introduction_to_elevatr.html#Key_information_about_version_0990_and_upcoming_versions_of_elevatr
      terra::rast() |> 
      terra::project(terra::crs('EPSG:4326'))
    
    #write DEM raster
    terra::writeRaster(areaDEM,
                       areaDEMFile,
                       overwrite = TRUE,
                       gdal=c("COMPRESS=DEFLATE"))
  } else {
    areaDEM <- terra::rast(areaDEMFile)
  }
  
  #Create plots that match with template
  plot_sf <- create.polygons.from.raster.and.points(raster = landsatARDTemplate, plotPoints)
  
  #Add functional plot buffer
  #This is area that can actually be used for data collection, and provides a guarantee of spatial overlap between the orthomosaics
  #and landsat pixels for testing unmixed product against classified imagery
  plot_sf <- plot_sf |> sf::st_buffer(10) #30 ft = 9.144 m. 
  
  # Print the final output
  print(plot_sf)
  
  
  #Export data
  
  #write plot packages
  sf::st_write(plot_sf, here::here(exportLocation, paste0(areaName, '.gpkg')), append = FALSE)
  sf::st_write(plot_sf, here::here(exportLocation, paste0(areaName, '.kml')), append = FALSE)
  
  #For each plot, export as .kml for UAS
  for(i in 1:nrow(plot_sf)) {
    
    #Write plot to shp
    plot <- plot_sf[i,]
    
    #Write plot to kml
    kmlnmP <- gsub(" ", "", paste(plot$plotID, ".kml", sep = ""))
    print(kmlnmP)
    
    #Write plot to kml for DJI use
    write.dji.kml(plot, kmlnmP, exportLocation)
    
    #Write plot to shapefile for use in metashape
    flNm <- gsub(" ", "", paste(plot$plotID, "_shp", sep = ""))
    st_write_shp(shp = plot,
                 location = exportLocation,
                 filename =  flNm,
                 zipOnly = FALSE,
                 overwrite = TRUE)
    #unlink(here::here(exportLocation, flNm, paste0(flNm, ".prj"))) #if need to remove .prj to make metashape accept it
    
    
    #Create pre-buffered plot for Map Pilot Pro
    plotBuffered <- plot |> sf::st_buffer(10) #30 ft = 9.144 m - this flight buffer is the ensure that the full area within the processed and clipped orthomosaic is usable and high quality
    kmlnmPBuffered <- gsub(" ", "", paste(plot$plotID, "_buffered.kml", sep = ""))
    sf::st_write(plotBuffered, here::here(exportLocation, kmlnmPBuffered), append = FALSE)
  }
  
}


# Create DJI-compatible KMLs


### ### ###
# A function to read in a KML file and turn it into a DJI-compatible KML, then export it
kml.to.dji.kml <- function(dir, kmlFile) {
  #Read KML file as raw txt
  kml <- readr::read_file(here::here(dir,kmlFile))
  
  #Get the coordinates from the original KML file
  coords <- kml %>% 
    stringr::str_match("<coordinates>\\s*(.*?)\\s*</coordinates>")
  coords <- coords[,2]
  
  #Add elevation to the original coordinates
  newCoords <- coords %>%
    gsub(" ", ",0 ", .) %>%
    paste0(., ",0", sep="")
  
  #Get name from original KML file
  nm <- kml %>%
    stringr::str_match("<name>\\s*(.*?)\\s*</name>")
  nm <- nm[,2]
  
  #Combine all text with the correctly formatted KML file text from a Google Earth Pro-created KML file
  newKMLtxt <- paste0("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<kml xmlns=\"http://www.opengis.net/kml/2.2\" xmlns:gx=\"http://www.google.com/kml/ext/2.2\" xmlns:kml=\"http://www.opengis.net/kml/2.2\" xmlns:atom=\"http://www.w3.org/2005/Atom\">\n<Document>\n\t<name>", nm, ".kml</name>\n\t<StyleMap id=\"m_ylw-pushpin\">\n\t\t<Pair>\n\t\t\t<key>normal</key>\n\t\t\t<styleUrl>#s_ylw-pushpin</styleUrl>\n\t\t</Pair>\n\t\t<Pair>\n\t\t\t<key>highlight</key>\n\t\t\t<styleUrl>#s_ylw-pushpin_hl</styleUrl>\n\t\t</Pair>\n\t</StyleMap>\n\t<Style id=\"s_ylw-pushpin\">\n\t\t<IconStyle>\n\t\t\t<scale>1.1</scale>\n\t\t\t<Icon>\n\t\t\t\t<href>http://maps.google.com/mapfiles/kml/pushpin/ylw-pushpin.png</href>\n\t\t\t</Icon>\n\t\t\t<hotSpot x=\"20\" y=\"2\" xunits=\"pixels\" yunits=\"pixels\"/>\n\t\t</IconStyle>\n\t</Style>\n\t<Style id=\"s_ylw-pushpin_hl\">\n\t\t<IconStyle>\n\t\t\t<scale>1.3</scale>\n\t\t\t<Icon>\n\t\t\t\t<href>http://maps.google.com/mapfiles/kml/pushpin/ylw-pushpin.png</href>\n\t\t\t</Icon>\n\t\t\t<hotSpot x=\"20\" y=\"2\" xunits=\"pixels\" yunits=\"pixels\"/>\n\t\t</IconStyle>\n\t</Style>\n\t<Placemark>\n\t\t<name>", nm, "</name>\n\t\t<styleUrl>#m_ylw-pushpin</styleUrl>\n\t\t<Polygon>\n\t\t\t<tessellate>1</tessellate>\n\t\t\t<outerBoundaryIs>\n\t\t\t\t<LinearRing>\n\t\t\t\t\t<coordinates>\n\t\t\t\t\t\t", newCoords, " \n\t\t\t\t\t</coordinates>\n\t\t\t\t</LinearRing>\n\t\t\t</outerBoundaryIs>\n\t\t</Polygon>\n\t</Placemark>\n</Document>\n</kml>\n", sep = "")
  
  #Sink to new kml file
  newFile <- paste(substr(kmlFile, 1, nchar(kmlFile)-4), "_for_dji.kml", sep="")
  sink(here::here(dir, newFile))
  cat(newKMLtxt)
  sink()
}

### ### ###
#The above, but as an EXPORT function
#This function takes in an sf object, exports it as a normal .kml, reads it back in, exports it as a dji.kml file, and deletes the original .kml file
#It is useful for keeping directories clean if you only need dji.kml files
#This function will download and load the 'sf' package if it is not already installed
#This function requires the function above (kml.to.dji.kml)
#sfObj is an sf object to write out, fileNm is the name of the file (e.g. mykml.kml), and outDir is the export directory (e.g.here::here('my', 'directory'))
write.dji.kml <- function(sfObj, fileNm, outDir) {
  #Check the required libraries and download if needed
  list.of.packages <- c("sf")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  library(sf)
  
  #Write as normal kml
  sf::st_write(sfObj, here(outDir, fileNm), append = FALSE)
  
  #Use kml.to.dji.kml function
  kml.to.dji.kml(outDir, fileNm)
  print(fileNm)
  
  #Remove original kml write-out
  if(file.exists(fileNm)) {
    unlink(here::here(outDir, fileNm))
  } else {
    print("File does not exist")
  }
}


# Data access ----


## Hydrosheds DEM data ----

# A function to mount and download slope data for the US
# Returns the slope data object
# PARAMETERS
# slopeF :: the path to which the slope data should be saved (if not already present).
#If already downloaded, data will be read without download
access_us_slope <- function(slopeF) {
  if(file.exists(slopeF)) {
    slope <- terra::rast(slopeF)
  } else {
    slope <- glue::glue(
      "/vsizip/vsicurl/", #magic remote connection 
      "https://data.hydrosheds.org/file/hydrosheds-v1-dem/hyd_na_dem_15s.zip", #copied link to download location
      "/hyd_na_dem_15s.tif") |> #path inside zip file
      terra::rast() |>
      terra::terrain("slope")
    terra::writeRaster(slope, slopeF,
                       gdal=c("COMPRESS=DEFLATE"))
  }
  return(slope)
}

## OSM Road Data ----


#A function to access road data from OSM
# PARAMETERS
# aoi :: an area of interest as an sf object - roads will be accessed within the area plus a 1km buffer
# Adapt function as necessary for filtering
access_osm_roads <- function(aoi) {
  roadsData <- osmdata::opq(bbox = sf::st_bbox(sf::st_transform(sf::st_buffer(aoi, 2000), 4326))) |>
    osmdata::add_osm_feature(key = "highway",
                             key_exact = FALSE,
                             value_exact = FALSE,
                             match_case = FALSE) |>
    osmdata::osmdata_sf()
  desiredColumns <- c("USFS", "highway", "access", "maintained", "motor_vehicle", "service", "smoothness", "surface", "tracktype")
  roads <- roadsData$osm_lines |>
    dplyr::select(dplyr::any_of(desiredColumns)) |>
    dplyr::filter(highway != "path" | is.na(highway)) |> 
    dplyr::filter(tracktype != "grade5" | is.na(tracktype)) |>
    dplyr::filter(access != "private" | is.na(access)) |>
    dplyr::mutate(group = 1) |>
    group_by(group) |>
    summarise(geometry = st_union(geometry)) |>
    ungroup() |>
    sf::st_transform(epsg) |>
    sf::st_intersection(sf::st_buffer(aoi, 1000)) #clip to district of interest + 1km
  
  return(roads)
}


## Administrative boundaries ----




access_ynp_bear_management_areas <- function() {
  # Write out the URL query
  base_url <- "https://services1.arcgis.com/fBc8EJBxQRMcHlei/arcgis/rest/services/YELL_BEAR_MANAGEMENT_AREAS_public_viewview/FeatureServer/0/query"
  query_params <- list(f = "json",
                       where = "1=1",
                       outFields = "*",
                       returnGeometry = "true")
  
  # Request data
  bma <- access_data_get_x_from_arcgis_rest_api_geojson(
    base_url = base_url, 
    query_params = query_params, 
    max_record = 1000, 
    n = "all", 
    timeout = 600
  )
  
  return(bma)
  
}

# The ranger districts file is quite small, so it is accessed via VSI
access_usfs_ranger_districts <- function() {
  usfs_rds <- paste0(
    "/vsizip/vsicurl/", #magic remote connection
    "https://data.fs.usda.gov/geodata/edw/edw_resources/shp/S_USA.RangerDistrict.zip", #copied link to download location
    "/S_USA.RangerDistrict.shp") |> #path inside zip file
    sf::st_read()
  return(usfs_rds)
}


# A function to access the US federal surface management agency polygon dataset
# The function downloads and unzips a geodatabase rather than accessing via VSI since this layer is 
# useful for visualization and field planning, as well as it being accessed multiple times
access_us_sma <- function(dir_path, layer) {
  
  loc <- here::here(dir_path, "SMA_WM.gdb")
  if(file.exists(loc)) {
    sma <- sf::st_read(loc, layer = layer)
  } else {
    download_unzip_file(url = "https://blm-egis.maps.arcgis.com/sharing/rest/content/items/6bf2e737c59d4111be92420ee5ab0b46/data",
                        extract_to = dir_path,
                        keep_zip = FALSE)
    sma <-  sf::st_read(loc, layer = layer)
  }
  return(sma)
}

access_us_sma_helper_show_layers <- function(dir_path) {
 sf::st_layers(here::here(dir_path, "SMA_WM.gdb"))
}

access_us_wilderness <- function(dest_path = NA) {
  if(is.na(dest_path)) {
    wild <- paste0(
      "/vsizip/vsicurl/", #magic remote connection
      "https://data.fs.usda.gov/geodata/edw/edw_resources/shp/S_USA.Wilderness.zip", #copied link to download location
      "/S_USA.Wilderness.shp") |> #path inside zip file
      sf::st_read()
  } else {
    if(file.exists(dest_path)) {
      wild <- sf::st_read(dest_path)
    } else {
      wild <- paste0(
        "/vsizip/vsicurl/", #magic remote connection
        "https://data.fs.usda.gov/geodata/edw/edw_resources/shp/S_USA.Wilderness.zip", #copied link to download location
        "/S_USA.Wilderness.shp") |> #path inside zip file
        sf::st_read()
      sf::st_write(wild, dest_path)
    }
  }
  return(wild)
}


# Queries the PADUS REST service for WSAs
# PADUS interactive online: https://usgs.maps.arcgis.com/home/item.html?id=98fce3fb0c8241ce8847e9f7d0d212e9
access_us_wilderness_study_areas <- function(dest_path = NA) {
  if(is.na(dest_path)) {
    query_params <- list(where = "DesTp_Desc='Wilderness Study Area'",
                         outFields = "*",
                         f = "json")
    base_url = "https://services.arcgis.com/v01gqwM5QqNysAAi/ArcGIS/rest/services/PADUS_Protection_Status_by_GAP_Status_Code/FeatureServer/0/QUERY"
    wsa <- access_data_get_x_from_arcgis_rest_api_geojson(base_url = base_url,
                                                          query_params = query_params,
                                                          max_record = 2000,
                                                          n = "all",
                                                          timeout = 500)
  } else {
    if(file.exists(dest_path)) {
      wsa <- sf::st_read(dest_path)
    } else {
      query_params <- list(where = "DesTp_Desc='Wilderness Study Area'",
                           outFields = "*",
                           f = "json")
      base_url = "https://services.arcgis.com/v01gqwM5QqNysAAi/ArcGIS/rest/services/PADUS_Protection_Status_by_GAP_Status_Code/FeatureServer/0/QUERY"
      wsa <- access_data_get_x_from_arcgis_rest_api_geojson(base_url = base_url,
                                                            query_params = query_params,
                                                            max_record = 2000,
                                                            n = "all",
                                                            timeout = 500)
      sf::st_write(wsa, dest_path)
    }
  }
  return(wsa)
}


## Landfire ----

#' Access LANDFIRE EVT Raster for CONUS (2023)
#'
#' This function remotely accesses and reads the 2023 LANDFIRE Existing Vegetation Type (EVT) raster data for the contiguous United States (CONUS). The data is accessed directly from a zipped online source using GDAL's VSI (Virtual File System) protocol.
#'
#' @details
#' The function utilizes GDAL's virtual file system (`/vsizip/vsicurl/`) to remotely access the LANDFIRE EVT raster file without needing to download or unzip it manually. The raster data is read into a `terra` raster object, suitable for geospatial analysis in R.
#'
#' @return
#' A `terra` raster object containing the 2023 LANDFIRE EVT data for CONUS.
#'
#' @examples
#' \dontrun{
#' lf_evt <- access_landfire_evt_conus_2023()
#' plot(lf_evt)
#' }
#' 
#' @importFrom terra rast
#' @export
access_landfire_evt_conus_2023 <- function() {
  lf_evt <- paste0(
    "/vsizip/vsicurl/", #magic remote connection
    "https://landfire.gov/data-downloads/US_240/LF2023_EVT_240_CONUS.zip", #copied link to download location
    "/LF2023_EVT_240_CONUS/Tif/LC23_EVT_240.tif") |> #path inside zip file
    terra::rast()
  return(lf_evt)
}

#' Access LANDFIRE EVT CSV for CONUS (2023)
#'
#' This function remotely accesses and reads the 2023 LANDFIRE Existing Vegetation Type (EVT) CSV data for the contiguous United States (CONUS). The data is accessed directly from a zipped online source using GDAL's VSI (Virtual File System) protocol.
#'
#' @details
#' This function uses GDAL's virtual file system (`/vsizip/vsicurl/`) to remotely access the LANDFIRE EVT CSV data without manual download or extraction. The CSV is read into an `sf` object using `sf::st_read()`, as GDAL's CSV handling is supported by spatial data functions. This method is necessary since standard R CSV readers do not natively support remote access via VSI.
#'
#' @return
#' An `sf` object containing the CSV data from the 2023 LANDFIRE EVT for CONUS.
#'
#' @examples
#' \dontrun{
#' lf_evt_csv <- access_landfire_evt_conus_2023_csv()
#' head(lf_evt_csv)
#' }
#' 
#' @importFrom sf st_read
#' @export
access_landfire_evt_conus_2023_csv <- function() {
  lf_evt_csv <- paste0(
    "/vsizip/vsicurl/", #magic remote connection
    "https://landfire.gov/data-downloads/US_240/LF2023_EVT_240_CONUS.zip", #copied link to download location
    "/LF2023_EVT_240_CONUS/CSV_Data/LF23_EVT_240.csv") |> #path inside zip file
    sf::st_read() #note that read_csv and other csv drivers in R don't talk to GDAL. Instead use st_read or terra::vect() to access CSV data in zip files
  return(lf_evt_csv)
}

#' Access LANDFIRE EVT Raster for CONUS (2022)
#'
#' This function remotely accesses and reads the 2022 LANDFIRE Existing Vegetation Type (EVT) raster data for the contiguous United States (CONUS). The data is accessed directly from a zipped online source using GDAL's VSI (Virtual File System) protocol.
#'
#' @details
#' The function utilizes GDAL's virtual file system (`/vsizip/vsicurl/`) to remotely access the LANDFIRE EVT raster file without needing to download or unzip it manually. The raster data is read into a `terra` raster object, suitable for geospatial analysis in R.
#'
#' @return
#' A `terra` raster object containing the 2022 LANDFIRE EVT data for CONUS.
#'
#' @examples
#' \dontrun{
#' lf_evt <- access_landfire_evt_conus_2022()
#' plot(lf_evt)
#' }
#' @param access Either "stream" or "download". "stream" will use VSI to access the dataset remotely.
#' "download" will check if the data is present within the directory given at dir_path, and, if not there, will download it before accessing it.
#'
#' @importFrom terra rast
#' @export
access_landfire_evt_conus_2022 <- function(access = "stream", dir_path = NA) {
  
  if(access == "stream") {
    lf_evt <- paste0(
      "/vsizip/vsicurl/", #magic remote connection
      "https://landfire.gov/data-downloads/US_230/LF2022_EVT_230_CONUS.zip", #copied link to download location
      "/LF2022_EVT_230_CONUS/Tif/LC22_EVT_230.tif") |> #path inside zip file
      terra::rast()
  }
  if(access == "download") {
    loc <- here::here(dir_path, "LF2022_EVT_230_CONUS/Tif/LC22_EVT_230.tif")
    if(file.exists(loc)) {
      lf_evt <- terra::rast(loc)
    } else {
      download_unzip_file(url = "https://landfire.gov/data-downloads/US_230/LF2022_EVT_230_CONUS.zip",
                          extract_to = dir_path,
                          keep_zip = FALSE)
      lf_evt <- terra::rast(loc)
    }
  }

  return(lf_evt)
}

#' Access LANDFIRE EVT CSV for CONUS (2022)
#'
#' This function remotely accesses and reads the 2022 LANDFIRE Existing Vegetation Type (EVT) CSV data for the contiguous United States (CONUS). The data is accessed directly from a zipped online source using GDAL's VSI (Virtual File System) protocol.
#'
#' @details
#' This function uses GDAL's virtual file system (`/vsizip/vsicurl/`) to remotely access the LANDFIRE EVT CSV data without manual download or extraction. The CSV is read into an `sf` object using `sf::st_read()`, as GDAL's CSV handling is supported by spatial data functions. This method is necessary since standard R CSV readers do not natively support remote access via VSI.
#'
#' @return
#' An `sf` object containing the CSV data from the 2022 LANDFIRE EVT for CONUS.
#'
#' @examples
#' \dontrun{
#' lf_evt_csv <- access_landfire_evt_conus_2022_csv()
#' head(lf_evt_csv)
#' }
#' 
#' @importFrom sf st_read
#' @export
access_landfire_evt_conus_2022_csv <- function() {
  lf_evt_csv <- paste0(
    "/vsizip/vsicurl/", #magic remote connection
    "https://landfire.gov/data-downloads/US_230/LF2022_EVT_230_CONUS.zip", #copied link to download location
    "/LF2022_EVT_230_CONUS/CSV_Data/LF22_EVT_230.csv") |> #path inside zip file
    sf::st_read() #note that read_csv and other csv drivers in R don't talk to GDAL. Instead use st_read or terra::vect() to access CSV data in zip files
  return(lf_evt_csv)
}

## NEON ----

#' Access AOP Flight Box Data
#'
#' Retrieves the flight box shapefile data for all NEON AOP sites by downloading and reading
#' it from the specified remote location.
#'
#' @return An sf object containing the flight box data for all NEON AOP sites.
#'
#' @importFrom sf st_read
#' 
#' @export
access_neon_aop_flight_box_data <- function() {
  aop_all <- paste0(
    "/vsizip/vsicurl/", # Magic remote connection
    "https://www.neonscience.org/sites/default/files/AOP_flightBoxes_0.zip", # Copied link to download location
    "/AOP_flightBoxes/AOP_flightboxesAllSites.shp") |> # Path inside zip file
    sf::st_read()
  return(aop_all)
}

#' Access NEON Plot Shapefiles
#'
#' This function remotely accesses and reads the NEON plot shapefiles directly from a zipped online source using GDAL's VSI (Virtual File System) protocol.
#' It downloads the shapefile for all NEON TOS (Tower Observation System) plots.
#'
#' @details
#' This function uses a magic remote connection through GDAL's virtual file system (`/vsizip/vsicurl/`) to access the shapefile directly from the zipped NEON data repository. The function does not require the user to download or unzip the file manually. The shapefile is read into an `sf` object for further spatial analysis in R.
#'
#' @return
#' An `sf` object containing the NEON TOS plot polygons.
#'
#' @examples
#' \dontrun{
#' neon_plots <- access_neon_plots_shp()
#' plot(neon_plots)
#' }
#' 
#' @importFrom sf st_read
#' @export
access_neon_plots_shp <- function() {
  neon_plots <- paste0(
    "/vsizip/vsicurl/", #magic remote connection
    "https://www.neonscience.org/sites/default/files/All_NEON_TOS_Plots_V10.zip", #copied link to download location
    "/All_NEON_TOS_Plots_V10/All_NEON_TOS_Plot_Polygons_V10.shp") |> #path inside zip file
    sf::st_read()
  return(neon_plots)
}

#' Access NEON Domain Shapefiles
#'
#' This function remotely accesses and reads the NEON domain shapefiles directly from a zipped online source using GDAL's VSI (Virtual File System) protocol.
#' It downloads the shapefile for NEON's defined geographic domains.
#'
#' @details
#' Similar to the `access_neon_plots_shp()` function, this function uses the GDAL virtual file system (`/vsizip/vsicurl/`) to access and read NEON domain shapefiles directly from a zipped source. The shapefile contains geographic domain boundaries for the NEON project, which is read into an `sf` object.
#'
#' @return
#' An `sf` object containing the NEON domain polygons.
#'
#' @examples
#' \dontrun{
#' neon_domains <- access_neon_domains_shp()
#' plot(neon_domains)
#' }
#' 
#' @importFrom sf st_read
#' @export
access_neon_domains_shp <- function() {
  neon_domains <- paste0(
    "/vsizip/vsicurl/", #magic remote connection
    "https://www.neonscience.org/sites/default/files/NEONDomains_0.zip", #copied link to download location
    "/NEON_Domains.shp") |> #path inside zip file
    sf::st_read()
  return(neon_domains)
}



## Ecoregions ----



#' Access EPA Level II Ecoregions Data via VSI
#'
#' This function retrieves the U.S. EPA Level II ecoregions shapefile from a remote server via VSI (Virtual Spatial Infrastructure).
#' The shapefile is stored in a ZIP file, and the function accesses it without downloading the file locally.
#'
#' @return A `sf` (simple features) object containing the EPA Level II ecoregions shapefile data.
#' 
#' @details
#' The function accesses the EPA Level II ecoregions shapefile directly from the EPA's data commons, utilizing the `/vsizip/` 
#' and `/vsicurl/` mechanisms to stream the shapefile from the zipped file. The file is accessed via a URL without the need to 
#' download it locally. This method allows efficient access to the shapefile data using the `sf` package.
#'
#' @source
#' U.S. EPA Ecoregions Data: \url{https://gaftp.epa.gov/EPADataCommons/ORD/Ecoregions/cec_na/}
#'
#' @references
#' U.S. EPA Ecoregions Information: \url{https://www.epa.gov/eco-research/ecoregions-north-america}
#'
#' @importFrom sf st_read
#' @export
#' @examples
#' # Example usage
#' epa_ecoregions <- access_data_epa_l2_ecoregions_vsi()
#'
access_data_epa_l2_ecoregions_vsi <- function() {
  epa_l2 <- paste0(
    "/vsizip/vsicurl/",
    "https://gaftp.epa.gov/EPADataCommons/ORD/Ecoregions/cec_na/na_cec_eco_l2.zip",
    "/NA_CEC_Eco_Level2.shp"
  ) |>
    sf::st_read()
  
  return(epa_l2)
}

#' Access EPA Level III Ecoregions Data via VSI
#'
#' This function retrieves the U.S. EPA Level III ecoregions shapefile from a remote server via VSI (Virtual Spatial Infrastructure).
#' The shapefile is stored in a ZIP file, and the function accesses it without downloading the file locally.
#'
#' @return A `sf` (simple features) object containing the EPA Level III ecoregions shapefile data.
#' 
#' @details
#' The function accesses the EPA Level III ecoregions shapefile directly from the EPA's data commons, utilizing the `/vsizip/` 
#' and `/vsicurl/` mechanisms to stream the shapefile from the zipped file. The file is accessed via a URL without the need to 
#' download it locally. This method allows efficient access to the shapefile data using the `sf` package.
#'
#' @source
#' U.S. EPA Ecoregions Data: \url{https://gaftp.epa.gov/EPADataCommons/ORD/Ecoregions/us/}
#' 
#' @references
#' U.S. EPA Ecoregions Information: \url{https://www.epa.gov/eco-research/ecoregions-north-america}
#'
#' @importFrom sf st_read
#' @export
#' @examples
#' # Example usage
#' epa_ecoregions <- access_data_epa_l3_ecoregions_vsi()
#'
access_data_epa_l3_ecoregions_vsi <- function() {
  epa_l3 <- paste0(
    "/vsizip/vsicurl/",
    "https://gaftp.epa.gov/EPADataCommons/ORD/Ecoregions/us/us_eco_l3.zip",
    "/us_eco_l3.shp"
  ) |>
    sf::st_read()
  
  return(epa_l3)
}



access_data_epa_l2_ecoregions_api <- function() {
  url <- "https://services2.arcgis.com/FiaPA4ga0iQKduv3/arcgis/rest/services/Ecological_Regions_of_North_America_Level_2/FeatureServer/0/QUERY"
  query_params <- list(where = "1=1",
                       outFields = "*",
                       f = "json")
  epa_l2 <- access_data_get_x_from_arcgis_rest_api_geojson(url,
                                                 query_params,
                                                 max_record = 2000,
                                                 n = "all",
                                                 timeout = 500)
  return(epa_l2)
}



# Utility ----


#' Install and Load Required Packages Using pak
#'
#' This function checks if the specified packages (both CRAN and GitHub) are installed and loads them. 
#' If any packages are missing, it installs them automatically.
#' It uses the `pak` package for faster and more efficient package installation.
#'
#' @param package_list A list of package names to check and install (non-string, e.g., `c(dplyr, here)`).
#' GitHub packages should be specified as `username/repo` in strings.
#' @param auto_install A character ("y" or "n", default is "n"). If "y", installs all required packages 
#' without asking for user permission. If "n", asks for permission from the user.
#' @return No return value. Installs and loads the specified packages as needed.
#' @examples
#' \dontrun{
#' install_and_load_packages(c(dplyr, here, "username/repo"))
#' }
#' @importFrom pak pkg_install
#' @export
install_and_load_packages <- function(package_list, auto_install = "n") {
  # Convert non-string package names to strings
  package_list <- lapply(package_list, function(pkg) {
    if (is.symbol(pkg)) {
      deparse(substitute(pkg))
    } else {
      pkg
    }
  })
  
  # # Check if 'renv' is installed; if not, skip the 'renv' check
  # if (requireNamespace("renv", quietly = TRUE) && renv::is_active()) {
  #   cat("renv is active. Only loading packages...\n")
  #   for (pkg in package_list) {
  #     package_name <- if (grepl("/", pkg)) unlist(strsplit(pkg, "/"))[2] else pkg
  #     if (!require(package_name, character.only = TRUE)) {
  #       cat("Failed to load package:", package_name, "\n")
  #     }
  #   }
  #   return(invisible())
  # }
  
  # Check if pak is installed; install if not
  if (!requireNamespace("pak", quietly = TRUE)) {
    cat("The 'pak' package is required for fast installation of packages, installing now.\n")
    install.packages("pak")
  }
  
  # Initialize lists to store missing CRAN and GitHub packages
  missing_cran_packages <- c()
  missing_github_packages <- c()
  
  # # Helper function to get user input
  # get_user_permission <- function(prompt_msg) {
  #   if (auto_install == "y") {
  #     return("y")
  #   } else {
  #     return(tolower(readline(prompt = prompt_msg)))
  #   }
  # }
  
  # Check for missing packages
  for (pkg in package_list) {
    if (grepl("/", pkg)) { # GitHub package
      package_name <- unlist(strsplit(pkg, "/"))[2]
      package_loaded <- require(package_name, character.only = TRUE, quietly = TRUE)
    } else { # CRAN package
      package_loaded <- require(pkg, character.only = TRUE, quietly = TRUE)
    }
    if (!package_loaded) {
      if (grepl("/", pkg)) {
        missing_github_packages <- c(missing_github_packages, pkg)
      } else {
        missing_cran_packages <- c(missing_cran_packages, pkg)
      }
    }
  }
  
  # Install missing CRAN packages using pak::pkg_install
  if (length(missing_cran_packages) > 0) {
    # cat("The following CRAN packages are missing: ", paste(missing_cran_packages, collapse = ", "), "\n")
    # response <- get_user_permission("\nDo you want to install the missing CRAN packages? (y/n): ")
    # if (response == "y") {
    pak::pkg_install(missing_cran_packages, upgrade = TRUE)
    # } else {
    #   cat("Skipping installation of missing CRAN packages.\n")
    # }
  }
  
  # Install missing GitHub packages using pak::pkg_install
  if (length(missing_github_packages) > 0) {
    # cat("The following GitHub packages are missing: ", paste(missing_github_packages, collapse = ", "), "\n")
    # response <- get_user_permission("\nDo you want to install the missing GitHub packages? (y/n): ")
    # if (response == "y") {
    pak::pkg_install(missing_github_packages, upgrade = TRUE)
    # } else {
    #   cat("Skipping installation of missing GitHub packages.\n")
    # }
  }
  
  # Load all packages after checking for installation
  for (pkg in package_list) {
    if (grepl("/", pkg)) { # GitHub package
      package_name <- unlist(strsplit(pkg, "/"))[2]
      if (!require(package_name, character.only = TRUE)) {
        cat("Failed to load GitHub package:", package_name, "\n")
      }
    } else { # CRAN package
      if (!require(pkg, character.only = TRUE)) {
        cat("Failed to load CRAN package:", pkg, "\n")
      }
    }
  }
  
  cat("All specified packages installed and loaded.\n")
}


#' Ensure Directory Exists
#'
#' This function checks if a directory exists at the specified path, and if not, creates a new directory.
#'
#' @param path A character string specifying the path to the new directory.
#' @return The function does not return any value. It creates a directory if it does not already exist.
#' @examples
#' # Ensure a directory named "data" exists
#' dir_ensure("data")
#'
#' @export
dir_ensure <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path)
    message("Directory created: ", path)
  } else {
    message("Directory already exists: ", path)
  }
}


# Function Description:
#   The function create.qgis.style.for.paletted.raster.from.csv generates a QGIS style file (.qml) for raster layers, specifically formatted for paletted rasters. It uses styling information provided in a data frame (styleData), allowing users to define colors and labels for different raster values. The function supports hexadecimal (hex) and RGB color schemes.
# 
# Parameters:
#   styleData: A data frame containing the styling information. The data frame should include the columns specified by valueColumn and labelColumn. For the hex color scheme, a color column is required. For RGB, columns R, G, and B are necessary.
# outputQmlPath: A string specifying the file path where the generated QGIS style file (.qml) will be saved.
# valueColumn: The name of the column in styleData that contains the raster values.
# labelColumn: The name of the column in styleData that contains the labels for each raster value.
# colorScheme: A string indicating the color scheme used in styleData. It can be "hex" for hexadecimal colors or "RGB" for separate red, green, and blue values. The default is "hex".
# Functionality:
#   The function iterates through each row of styleData, extracting the value, label, and color information to create palette entries in the QML file. For RGB color schemes, it converts the RGB values to hex using the rgb function. The function ensures proper XML formatting by escaping special characters in labels. After constructing the QML content, it is written to the specified output path.
# 
# Usage Example:
# # Assuming styleData is pre-defined with the appropriate columns
# create_qgis_style_for_paletted_raster_from_csv(styleData, "path/to/output.qml", "value", "label", "hex")
# Citation:
#   Function authored by R Code Stylist, GPT-4, OpenAI, in collaboration with the user.
create_qgis_style_for_paletted_raster_from_csv <- function(styleData, outputQmlPath, valueColumn, labelColumn, colorScheme = "hex") {
  
  # Check for the necessary columns in the CSV based on the color scheme
  if (!labelColumn %in% colnames(styleData)) {
    stop("CSV file must contain the specified label column.")
  }
  
  if (colorScheme == "hex" && !("color" %in% colnames(styleData))) {
    stop("CSV file must contain a 'color' column for hex color scheme.")
  } else if (colorScheme == "RGB" && !all(c("R", "G", "B") %in% colnames(styleData))) {
    stop("CSV file must contain 'R', 'G', 'B' columns for RGB color scheme.")
  }
  
  # Start creating the QML content
  qmlContent <- '<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE qgis PUBLIC \'http://mrcc.com/qgis.dtd\' \'SYSTEM\'>
<qgis hasScaleBasedVisibilityFlag="0" maxScale="0" version="3.22.12-Białowieża" styleCategories="AllStyleCategories" minScale="1e+08">
   <flags>
    <Identifiable>1</Identifiable>
    <Removable>1</Removable>
    <Searchable>1</Searchable>
    <Private>0</Private>
  </flags>
  <temporal enabled="0" fetchMode="0" mode="0">
    <fixedRange>
      <start></start>
      <end></end>
    </fixedRange>
  </temporal>
  <customproperties>
    <Option type="Map">
      <Option value="false" type="bool" name="WMSBackgroundLayer"/>
      <Option value="false" type="bool" name="WMSPublishDataSourceUrl"/>
      <Option value="0" type="int" name="embeddedWidgets/count"/>
      <Option value="Value" type="QString" name="identify/format"/>
    </Option>
  </customproperties>
  <pipe-data-defined-properties>
    <Option type="Map">
      <Option value="" type="QString" name="name"/>
      <Option name="properties"/>
      <Option value="collection" type="QString" name="type"/>
    </Option>
  </pipe-data-defined-properties>
  <pipe>
    <provider>
      <resampling enabled="false" zoomedInResamplingMethod="nearestNeighbour" zoomedOutResamplingMethod="nearestNeighbour" maxOversampling="2"/>
    </provider>
    <rasterrenderer opacity="1" nodataColor="" type="paletted" band="1" alphaBand="-1">
      <rasterTransparency/>
      <minMaxOrigin>
        <limits>None</limits>
        <extent>WholeRaster</extent>
        <statAccuracy>Estimated</statAccuracy>
        <cumulativeCutLower>0.02</cumulativeCutLower>
        <cumulativeCutUpper>0.98</cumulativeCutUpper>
        <stdDevFactor>2</stdDevFactor>
      </minMaxOrigin>
  <colorPalette>'  
  # Append palette entries from the CSV data
  for (i in 1:nrow(styleData)) {
    # Determine the color based on the scheme
    if (colorScheme == "hex") {
      color <- styleData$color[i]
    } else {
      color <- rgb(red = styleData$R[i], green = styleData$G[i], blue = styleData$B[i], maxColorValue = 255)
    }
    
    label <- styleData[[labelColumn]][i]
    label <- gsub("&", "and", label)
    label <- gsub('\\"', '', label)
    value <- styleData[[valueColumn]][i]
    
    
    qmlContent <- glue::glue('{qmlContent}
             <paletteEntry value="{value}" label="{label}" alpha="255" color="{color}"/>'
    )
  }
  
  # Finalize the QML content with closing tags
  qmlContent <- paste0(qmlContent, '\n      </colorPalette>
          <colorramp type="randomcolors" name="[source]">
        <Option/>
      </colorramp>
    </rasterrenderer>
    <brightnesscontrast gamma="1" brightness="0" contrast="0"/>
    <huesaturation grayscaleMode="0" colorizeOn="0" colorizeGreen="128" saturation="0" colorizeBlue="128" colorizeRed="255" colorizeStrength="100" invertColors="0"/>
    <rasterresampler maxOversampling="2"/>
    <resamplingStage>resamplingFilter</resamplingStage>
  </pipe>
  <blendMode>0</blendMode>
</qgis>')
  
  # Write the QML content to a file
  writeLines(qmlContent, outputQmlPath)
  
  print(qmlContent)
  
  return(paste("QGIS style file created at:", outputQmlPath))
}

#' Write Shapefile to a New Directory and Create a Zipped Version
#'
#' This function writes an `sf` object to a shapefile in a new, file-specific directory and optionally creates a zipped version of the shapefile.
#' It also allows for the removal of the original unzipped files and handles overwriting existing files.
#'
#' @param shp An `sf` object to write as a shapefile.
#' @param location A character string specifying the path of the directory to create the new file-specific subdirectory in.
#' @param filename A character string specifying the name of the file without the `.shp` extension.
#' @param zip_only A logical value indicating whether the original (unzipped) files should be removed after zipping. Defaults to `FALSE`.
#' @param overwrite A logical value indicating whether existing files should be overwritten. Defaults to `FALSE`.
#' @return No return value. The function writes a shapefile to a specified directory, optionally zips the files, and manages file cleanup based on user input.
#' @examples
#' \dontrun{
#' # Example usage
#' st_write_shp(shp = prepped_for_parks_etal,
#'              location = here::here("data/derived"),
#'              filename = "career_lba_for_parks_v1",
#'              zip_only = TRUE,
#'              overwrite = TRUE)
#' }
#' @importFrom sf st_write
#' @importFrom zip zip
#' @export
st_write_shp <- function(shp, location, filename, zip_only = FALSE, overwrite = FALSE) {
  
  # Define paths
  out_dir <- file.path(location, filename)
  zip_file <- file.path(out_dir, paste0(filename, ".zip"))
  zip_file_dest <- file.path(location, paste0(filename, ".zip"))
  
  # Manage overwriting and directory creation
  if (dir.exists(out_dir)) {
    if (overwrite) {
      unlink(out_dir, recursive = TRUE)
    } else {
      stop("Directory '", out_dir, "' already exists and overwrite is set to FALSE.")
    }
  }
  
  if (file.exists(zip_file_dest) && zip_only) {
    if (overwrite) {
      unlink(zip_file_dest)
    } else {
      stop("Zip file '", zip_file_dest, "' already exists and overwrite is set to FALSE.")
    }
  }
  
  # Create the directory if not there
  dir_ensure(out_dir)
  
  # Write the shapefile
  shapefile_path <- file.path(out_dir, paste0(filename, ".shp"))
  sf::st_write(shp, shapefile_path, append = FALSE)
  
  # Get all shapefile components
  all_shp_files <- list.files(out_dir, pattern = paste0(filename, ".*"), full.names = TRUE)
  
  # Create zip file
  zip::zip(zipfile = zip_file, files = all_shp_files, mode = "cherry-pick")
  
  # Remove raw files if zip_only is TRUE
  if (zip_only) {
    file.copy(zip_file, zip_file_dest)
    unlink(out_dir, recursive = TRUE)
  }
}


#Function to clip a raster to a vector, ensuring in same projection
#Returns raster in original projection, but clipped to vector
#Returns raster in the same form that it came in
# PARAMETERS
# raster : a SpatRaster, PackedSpatRaster, RasterLayer, RasterStack, or RasterBrick object
# vector : a SpatVector, PackedSpatVector or SF object
# mask : TRUE or FALSE; whether terra::clip should mask the raster as well
crop_careful_universal <- function(raster, vector, mask, verbose = FALSE) {
  pack <- FALSE
  
  #Unpack if parallelized inputs
  if(class(raster)[1] == "PackedSpatRaster") {
    raster <- terra::unwrap(raster)
    pack <- TRUE
  }
  if(class(vector)[1] == "PackedSpatVector") {
    vector <- sf::st_as_sf(terra::unwrap(vector))
  }
  
  #Handle unpacked spatVector
  if(class(vector)[1] == "SpatVector") {
    vector <- sf::st_as_sf(vector)
  }
  
  #If using raster package
  if(class(raster)[1] == "RasterLayer" | class(raster)[1] == "RasterStack" | class(raster)[1] == "RasterBrick") {
    
    #Perform operation
    if (raster::crs(vector) != raster::crs(raster)) { #if raster and vector aren't in same projection, change vector to match
      if(verbose) {print("Projecting vector")}
      vector <- sf::st_transform(vector, raster::crs(raster)) 
    } else {
      if(verbose) {print("Vector already in raster CRS")}
    }
    if(verbose) {print("Clipping")}
    r <- raster::crop(raster,
                      vector)
    if(mask) {
      r <- r |> raster::mask(vector)
    }
    
    return(r)
    
  } else { #terra package
    
    #Perform operation
    if (terra::crs(vector) != terra::crs(raster)) { #if raster and vector aren't in same projection, change vector to match
      if(verbose) {print("Projecting vector")}
      vector <- sf::st_transform(vector, terra::crs(raster)) 
    } else {
      if(verbose) {print("Vector already in raster CRS")}
    }
    if(verbose) {print("Clipping")}
    r <- terra::crop(raster,
                     vector,
                     mask = mask) #crop & mask
    
    #Repack if was packed coming in (i.e. parallelized)
    if(pack) {
      r <- terra::wrap(r)
    }
    return(r)
    
  }
}

#' Substring from the Right
#'
#' This function extracts a substring from the right side of a given string, retaining a specified number of characters.
#' It will work for either single strings or a vectorized input (e.g. when used on a data frame)
#'
#' @param str A character string from which to extract the substring.
#' @param n An integer specifying the number of characters to keep from the right of the string.
#' @return A character string containing the rightmost \code{n} characters of the input string.
#' @examples
#' # Extract the last 3 characters from a string
#' substr_right("Hello, World!", 3)
#'
#' @export
substr_right <- function(str, n) {
  # Check if input is a vector, and apply the function element-wise using sapply
  if (length(str) > 1) {
    return(sapply(str, function(x) {
      if (n > nchar(x)) {
        warning("n is greater than the length of the string. Returning the full string.")
        return(x)
      }
      return(substr(x, nchar(x) - n + 1, nchar(x)))
    }))
  }
  
  # If input is a single string, apply the logic directly
  if (n > nchar(str)) {
    warning("n is greater than the length of the string. Returning the full string.")
    return(str)
  }
  
  return(substr(str, nchar(str) - n + 1, nchar(str)))
}

#' Download and Unzip a File
#'
#' Downloads a ZIP file from a specified URL and extracts its contents to a specified directory.
#' Optionally, the ZIP file can be retained after extraction.
#'
#' @param url Character. The URL of the ZIP file to download.
#' @param extract_to Character. The directory where the contents should be extracted.
#' @param keep_zip Logical. If `TRUE`, retains the ZIP file after extraction. Defaults to `FALSE`.
#'
#' @return Invisible `NULL`. The function is used for its side effects of downloading and extracting files.
#'
#' @details The function downloads a ZIP file from a URL and extracts its contents to a specified directory.
#' If `keep_zip` is set to `FALSE`, the ZIP file will be deleted after extraction.
#'
#' @importFrom utils download.file unzip
#' @export
#'
#' @examples
#' \dontrun{
#' download_unzip_file("https://example.com/data.zip", "path/to/extract", keep_zip = TRUE)
#' }
download_unzip_file <- function(url, extract_to, keep_zip = FALSE) {
  # Validate URL and extraction path
  if (!is.character(url) || length(url) != 1) stop("`url` must be a single character string.")
  if (!is.character(extract_to) || length(extract_to) != 1) stop("`extract_to` must be a single character string.")
  if (!is.logical(keep_zip) || length(keep_zip) != 1) stop("`keep_zip` must be a single logical value.")
  
  # Ensure the extraction directory exists
  if (!dir.exists(extract_to)) dir.create(extract_to, recursive = TRUE)
  
  # Determine the path to save the ZIP file
  zip_path <- if (keep_zip) {
    # Save the ZIP file to the specified extraction directory
    file.path(extract_to, basename(url))
  } else {
    # Use a temporary file path for the ZIP file
    tempfile(fileext = ".zip")
  }
  
  # Ensure temporary file cleanup if there's an error and keep_zip is FALSE
  on.exit({
    if (!keep_zip && file.exists(zip_path)) {
      unlink(zip_path)
    }
  }, add = TRUE)
  
  # Attempt to download the ZIP file
  tryCatch({
    download.file(url, zip_path, mode = "wb")
  }, error = function(e) {
    stop("Failed to download the file from the specified URL: ", e$message)
  })
  
  # Attempt to unzip the file to the specified extraction directory
  tryCatch({
    unzip(zip_path, exdir = extract_to)
  }, error = function(e) {
    stop("Failed to unzip the file: ", e$message)
  })
  
  # Delete the ZIP file if 'keep_zip' is FALSE
  if (!keep_zip) {
    unlink(zip_path)
  }
  
  gc()
  
  invisible(NULL)
}

#' Merge a List of Raster Files and Optionally Write to Disk
#'
#' Merges a list of raster files into a single raster object. The merged raster can either
#' be saved to a specified file path or returned as an in-memory object.
#'
#' @param file_list Character vector. A list of file paths to the raster files to be merged.
#' @param file_final_path Character. The file path where the merged raster will be saved if `write = TRUE`.
#' @param datatype Character. The data type of the output raster. Defaults to `"INT2U"`.
#' @param compress Logical. If `TRUE`, compresses the output file with DEFLATE compression when writing to disk. Defaults to `TRUE`.
#' @param write Logical. If `TRUE`, writes the merged raster to `file_final_path`. If `FALSE`, returns the merged raster in memory. Defaults to `TRUE`.
#'
#' @return If `write = TRUE`, returns invisible `NULL` after writing to disk. If `write = FALSE`, returns the merged raster object.
#'
#' @details This function reads a list of raster files, merges them, and either writes the merged raster to a specified path
#' or returns it in memory. Compression is available when writing to disk to reduce file size.
#'
#' @importFrom purrr map
#' @importFrom terra rast sprc merge writeRaster
#' @export
#'
#' @examples
#' \dontrun{
#' file_paths <- c("path/to/raster1.tif", "path/to/raster2.tif")
#' # To write to disk
#' merge_list_of_rasters(file_paths, "path/to/final_raster.tif", datatype = "FLT4S", compress = TRUE, write = TRUE)
#' # To return in memory
#' merged_raster <- merge_list_of_rasters(file_paths, write = FALSE)
#' }
merge_list_of_rasters <- function(file_list, file_final_path = NULL, datatype = "INT2U", compress = TRUE, write = TRUE) {
  # Validate inputs
  if (!is.character(file_list) || length(file_list) < 1) stop("`file_list` must be a non-empty character vector.")
  if (write && (is.null(file_final_path) || !is.character(file_final_path) || length(file_final_path) != 1)) {
    stop("When `write = TRUE`, `file_final_path` must be a single, non-null character string.")
  }
  if (!is.logical(compress) || length(compress) != 1) stop("`compress` must be a single logical value.")
  if (!is.logical(write) || length(write) != 1) stop("`write` must be a single logical value.")
  
  # Load and merge the rasters
  combined_rasters <- file_list |>
    purrr::map(~ terra::rast(.x)) |>
    terra::sprc() |>
    terra::merge()
  
  # Write or return the merged raster
  if (write) {
    if (compress) {
      terra::writeRaster(combined_rasters,
                         file_final_path,
                         datatype = datatype,
                         gdal = c("COMPRESS=DEFLATE"))
    } else {
      terra::writeRaster(combined_rasters,
                         file_final_path,
                         datatype = datatype)
    }
    invisible(NULL)  # Return NULL after writing to disk
  } else {
    return(combined_rasters)  # Return the merged raster in memory
  }
}


## ArcGIS REST Data access


#' Fetch Data from an ArcGIS REST API Endpoint with Pagination
#'
#' This function retrieves geojson data from an ArcGIS REST API endpoint using pagination. It supports fetching a specified
#' number of entries or all available entries from the API endpoint. Written with ChatGPT 4o assistance.
#'
#' @param base_url A character string. The base URL of the ArcGIS REST API endpoint.
#' @param query_params A list. Query parameters to be used in the API request. The list should contain the necessary
#' parameters required by the API, such as `where`, `outFields`, and `f`.
#' @param max_record An integer. The maximum number of records that can be fetched in a single API request. This value is
#' usually defined by the ArcGIS REST API server limitations.
#' @param n An integer or character. Specifies the total number of entries to fetch. If `"all"`, the function fetches
#' all available records from the API. If an integer, it specifies the exact number of records to fetch.
#' @param timeout An integer. The time in seconds to wait before timing out the request.
#'
#' @return An `sf` object. A Simple Features (sf) object containing the fetched data.
#' @import httr sf
#' @examples
#' \dontrun{
#' base_url <- "https://example.com/arcgis/rest/services/your_service/FeatureServer/0/query"
#' query_params <- list(where = "1=1", outFields = "*", f = "geojson")
#' max_record <- 100
#' n <- 500  # Can also be "all"
#' result <- get.x.from.arcgis.rest.api(base_url, query_params, max_record, n)
#' print(result)
#' }
#' @importFrom httr GET status_code content timeout
#' @importFrom sf st_read
#' @export
access_data_get_x_from_arcgis_rest_api_geojson <- function(base_url, query_params, max_record, n, timeout) {
  # Input validation
  if (!is.character(base_url) || length(base_url) != 1) {
    stop("Parameter 'base_url' must be a single character string.")
  }
  if (!is.list(query_params)) {
    stop("Parameter 'query_params' must be a list.")
  }
  if (!is.numeric(max_record) || length(max_record) != 1 || max_record <= 0) {
    stop("Parameter 'max_record' must be a positive integer.")
  }
  if (!is.numeric(timeout) || length(timeout) != 1 || timeout <= 0) {
    stop("Parameter 'timeout' must be a positive integer.")
  }
  
  
  # Initialize variables
  total_features <- list()
  offset <- 0
  total_fetched <- 0  # Keep track of the total number of records fetched
  
  # Determine the limit for fetching records
  fetch_all <- FALSE
  if (n == "all") {
    fetch_all <- TRUE
  } else if (!is.numeric(n) || n <= 0) {
    stop("Parameter 'n' must be a positive integer or 'all'.")
  }
  
  repeat {
    # Update the resultOffset parameter in query_params
    query_params$resultOffset <- offset
    
    # Make the GET request using the base URL and query parameters
    response <- httr::GET(url = base_url, query = query_params, httr::timeout(timeout))
    
    # Check if the request was successful
    if (httr::status_code(response) == 200) {
      # Read the GeoJSON data directly into an sf object
      data <- sf::st_read(httr::content(response, as = "text"), quiet = TRUE)
      
      # Append the data to the list of features
      total_features <- append(total_features, list(data))
      
      # Update the total number of fetched records
      total_fetched <- total_fetched + nrow(data)
      
      # Provide user feedback for long-running processes
      cat(sprintf("Fetched %d records so far...\n", total_fetched))
      
      # Determine if we should stop fetching
      if ((nrow(data) < max_record) || (!fetch_all && total_fetched >= n)) {
        break
      }
      
      # Increment the offset by the maximum number of records for the next page
      offset <- offset + max_record
    } else {
      # Handle errors and provide meaningful messages
      error_message <- httr::content(response, "text", encoding = "UTF-8")
      stop("Failed to fetch data: ", httr::status_code(response), " - ", error_message)
    }
  }
  
  # Combine all pages into one sf object
  all_data_sf <- do.call(rbind, total_features)
  
  # If n is not "all", limit the output to the first n records
  if (!fetch_all) {
    all_data_sf <- all_data_sf[1:min(n, nrow(all_data_sf)), ]
  }
  
  return(all_data_sf)
}

