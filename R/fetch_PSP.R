#' Download and Process POLARIS Soil Data
#'
#' Downloads specified POLARIS soil property data from Google Earth Engine (GEE),
#' loads the rasters into memory, and optionally saves them to disk.
#'
#' @param aoi SpatRaster, SpatVector, or file path to a shapefile. The area of interest.
#' @param properties Character vector. List of soil properties to download (e.g., `"bd_mean"`, `"clay_mean"`, `"om_mean"`).
#' @param depths Character vector. Depth intervals to download (e.g., `"0_5"`, `"5_15"`).
#' @param output_dir Character. Directory where files will be saved. Set to `NULL` to keep only in memory.
#' @param suffix Character. Optional suffix for file names.
#' @param crs Character. Coordinate reference system for the output. Default is `"EPSG:4326"`.
#' @param scale Numeric. Resolution in meters. Default is `30`.
#' @param export Logical. If `TRUE`, saves downloaded rasters to `output_dir`. Default is `TRUE`.
#' @param tosoc Logical. If `TRUE`, converts Organic Matter (`om_mean`) to Soil Organic Carbon (SOC) using the Van Bemmelen factor (1.724).
#' @param convertOM Logical. If `TRUE`, converts Organic Matter (`om_mean`) from log10 to percentage. Default is `FALSE`.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{stack}}{A `SpatRaster` containing all downloaded layers.}
#'   \item{\code{file_paths}}{A character vector of file paths (if `export = TRUE`).}
#'   \item{\code{product}}{A character string `"PSP"` identifying the dataset.}
#' }
#'
#' @details
#' This function interfaces with the POLARIS soil dataset on Google Earth Engine (GEE),
#' allowing users to extract soil property maps for a specified area of interest (AOI).
#' Users can optionally apply transformations to Organic Matter (OM) and Soil Organic Carbon (SOC).
#'
#' @import terra rgeedim
#' @export
#'
#' @examples
#' \dontrun{
#' # Authenticate and initialize Google Earth Engine
#' gd_authenticate(auth_mode = "notebook")
#' gd_initialize()
#'
#' # Define AOI and output directory
#' output_directory <- "path/to/your/directory"
#' aoi_raster <- rast("path_to_aoi.tif")
#'
#' # Fetch POLARIS data (without OM to SOC conversion)
#' psp_data <- fetch_PSP(
#'   aoi = aoi_raster,
#'   properties = c("bd_mean", "clay_mean"),
#'   depths = c("0_5", "30_60"),
#'   output_dir = output_directory,
#'   suffix = "kit",
#'   crs = "EPSG:4326",
#'   scale = 30,
#'   export = TRUE,
#'   convertOM = FALSE,
#'   tosoc = FALSE
#' )
#'
#' # Fetch POLARIS data with OM converted to percentage
#' psp_data_om <- fetch_PSP(
#'   aoi = aoi_raster,
#'   properties = c("om_mean", "clay_mean"),
#'   depths = c("0_5", "30_60"),
#'   output_dir = output_directory,
#'   suffix = "kit",
#'   crs = "EPSG:4326",
#'   scale = 30,
#'   export = TRUE,
#'   convertOM = TRUE,
#'   tosoc = FALSE
#' )
#'
#' # Fetch POLARIS data with OM converted to SOC
#' psp_data_soc <- fetch_PSP(
#'   aoi = aoi_raster,
#'   properties = c("om_mean", "clay_mean"),
#'   depths = c("0_5", "30_60"),
#'   output_dir = output_directory,
#'   suffix = "kit",
#'   crs = "EPSG:4326",
#'   scale = 30,
#'   export = TRUE,
#'   convertOM = TRUE,
#'   tosoc = TRUE
#' )
#'
#' # Access the loaded SpatRaster
#' plot(psp_data$stack)
#' }

fetch_PSP <- function(
  aoi,
  properties,
  depths,
  output_dir = NULL,
  suffix = "",
  crs = "EPSG:4326",
  scale = 30,
  export = TRUE,
  tosoc = FALSE,
  convertOM = FALSE
) {
  # Load AOI as SpatVector if it's a shapefile
  if (is.character(aoi) && grepl("\\.shp$", aoi)) {
    aoi <- vect(aoi)
  }

  # Convert AOI to WGS84 if needed
  if (inherits(aoi, "SpatRaster") || inherits(aoi, "SpatVector")) {
    if (crs(aoi) != "EPSG:4326") {
      aoi <- project(aoi, "EPSG:4326")
    }
    bbox <- ext(aoi)
  } else {
    stop("AOI must be a SpatRaster, SpatVector, or a shapefile path.")
  }

  # Convert bounding box to GEE format
  bbox_gee <- gd_bbox(
    xmin = bbox$xmin,
    xmax = bbox$xmax,
    ymin = bbox$ymin,
    ymax = bbox$ymax
  )

  # POLARIS base path
  polaris_base <- "projects/sat-io/open-datasets/polaris"

  # Define available depth intervals and soil properties
  depth_intervals <- c("0_5", "5_15", "15_30", "30_60", "60_100", "100_200")
  polaris_properties <- c(
    "bd_mean",
    "clay_mean",
    "ksat_mean",
    "n_mean",
    "om_mean",
    "ph_mean",
    "sand_mean",
    "silt_mean",
    "theta_r_mean",
    "theta_s_mean",
    "lambda_mean",
    "hb_mean",
    "alpha_mean"
  )

  # Validate input properties and depths
  invalid_properties <- setdiff(properties, polaris_properties)
  if (length(invalid_properties) > 0) {
    stop(paste(
      "Invalid property selected:",
      paste(invalid_properties, collapse = ", ")
    ))
  }

  selected_depths <- intersect(depths, depth_intervals)
  if (length(selected_depths) == 0) {
    stop("Invalid depth(s) selected.")
  }

  raster_list <- list()
  file_paths <- list()

  for (prop in properties) {
    polaris_path <- paste0(polaris_base, "/", prop)

    for (depth in selected_depths) {
      depth_id <- paste0(prop, "_", depth)
      # Remove "_mean" from the property name when constructing img_path
      prop_clean <- gsub("_mean$", "", prop)
      img_path <- paste0(polaris_path, "/", prop_clean, "_", depth)

      img <- gd_image_from_id(img_path)

      # Adjust property name for output based on transformations
      final_prop_name <- prop
      if (prop == "om_mean") {
        if (tosoc) {
          final_prop_name <- "soc"
        } else if (convertOM) {
          final_prop_name <- "om_perc"
        }
      }

      # Define output file name
      output_file <- if (!is.null(output_dir) && export) {
        file.path(
          output_dir,
          paste0("PSP_", final_prop_name, "_", depth, suffix, ".tif")
        )
      } else {
        tempfile(fileext = ".tif")
      }

      # Download raster
      gd_download(
        img,
        filename = output_file,
        region = bbox_gee,
        scale = scale,
        crs = crs,
        resampling = "near",
        overwrite = TRUE,
        silent = FALSE
      )

      if (export) {
        file_paths <- append(file_paths, output_file)
      }

      # Load raster and remove extra bands
      full_raster <- rast(output_file)
      clean_raster <- full_raster[[1]]

      # Convert OM from log10 to percentage if selected
      if (prop == "om_mean" && convertOM) {
        clean_raster <- 10^clean_raster
        message("Converted Organic Matter from log10 to percentage.")
      }

      # Convert OM to SOC if selected
      if (prop == "om_mean" && tosoc) {
        clean_raster <- clean_raster / 1.724
        message(
          "Converted Organic Matter to Soil Organic Carbon using Van Bemmelen factor (1.724)."
        )
      }

      # Assign correct band name
      final_band_name <- paste0(final_prop_name, "_", depth)
      names(clean_raster) <- final_band_name

      # Save the cleaned raster
      temp_file <- paste0(output_file, "_tmp.tif")
      writeRaster(clean_raster, temp_file, overwrite = TRUE)
      file.rename(temp_file, output_file)

      raster_list[[final_band_name]] <- rast(output_file)

      message("Downloaded and cleaned: ", output_file)
    }
  }

  # Ensure at least one raster was downloaded
  if (length(raster_list) == 0) {
    stop("No rasters were successfully downloaded.")
  }

  # Create a raster stack from the list of downloaded rasters
  raster_stack <- rast()
  for (r in raster_list) {
    raster_stack <- c(raster_stack, r)
  }

  # Assign correct band names
  names(raster_stack) <- names(raster_list)

  if (nlyr(raster_stack) == 0) {
    stop("Failed to create a valid SpatRaster stack.")
  }

  # Assign "PSP" as an attribute to the SpatRaster
  attr(raster_stack, "product") <- "PSP"

  # Return structured output
  return(list(stack = raster_stack, file_paths = file_paths, product = "PSP"))
}
