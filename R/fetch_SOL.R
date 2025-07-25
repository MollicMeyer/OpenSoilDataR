#' Fetch SOLUS100 Soil Data
#'
#' Downloads and processes soil property data from SOLUS100, loading rasters into memory
#' while optionally saving them to disk.
#'
#' @param aoi SpatRaster or SpatVector. Area of interest.
#' @param properties Character vector. List of soil properties to download.
#' @param depths Character vector. Soil depth intervals to download.
#' @param measures Character vector. Measures to download ("l", "h", "p", "rpi").
#' @param output_dir Character. Directory where files will be saved. Set to NULL to keep in memory.
#' @param suffix Character. Optional suffix for filenames.
#' @param export Logical. If TRUE, saves downloaded rasters to `output_dir`.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{stack}}{A `SpatRaster` containing all downloaded layers.}
#'   \item{\code{file_paths}}{A character vector of file paths (if `export = TRUE`).}
#'   \item{\code{product}}{A character string `"SOL"` identifying the dataset.}
#' }
#'
#' @import terra httr utils
#' @export
#' @examples
#' \dontrun{
#' # Define AOI and fetch SOLUS100 data
#' aoi <- rast("path_to_aoi.tif")
#' output_directory <- "path/to/your/directory"
#' sol_data <- fetch_SOL(
#'   aoi = aoi,
#'   properties = c("claytotal", "soc"),
#'   depths = c("0_cm", "30_cm"),
#'   measures = c("p", "l", "h"), # rpi is data heavy
#'   output_dir = output_directory,
#'   suffix = "test",
#'   export = TRUE
#' )
#'
#' # Access the loaded SpatRaster
#' plot(sol_data$stack)
#' }

fetch_SOL <- function(
  aoi,
  properties,
  depths,
  measures,
  output_dir = NULL,
  suffix = "",
  export = TRUE
) {
  base_url <- "https://storage.googleapis.com/solus100pub/"

  # Load scalar lookup from package inst/extdata
  lookup_path <- system.file(
    "extdata",
    "scalar_lookup.csv",
    package = "OpenSoilDataR"
  )
  if (lookup_path == "") {
    stop("Could not find scalar_lookup.csv in inst/extdata/")
  }
  scalar_table <- read.csv(lookup_path)

  # Available SOLUS soil properties
  valid_properties <- c(
    "anylithicdpt",
    "caco3",
    "cec7",
    "claytotal",
    "dbovendry",
    "ec",
    "ecec",
    "fragvol",
    "gypsum",
    "ph1to1h2o",
    "resdept",
    "sandco",
    "sandfine",
    "sandmed",
    "sandtotal",
    "sandvc",
    "sandvf",
    "sar",
    "silttotal",
    "soc"
  )

  # Available depth intervals
  valid_depths <- c(
    "0_cm",
    "5_cm",
    "15_cm",
    "30_cm",
    "60_cm",
    "100_cm",
    "150_cm"
  )
  valid_measures <- c("l", "h", "p", "rpi")

  # Validate input parameters
  if (!all(properties %in% valid_properties)) {
    stop("Invalid property selected.")
  }
  if (!all(depths %in% valid_depths)) {
    stop("Invalid depth selected.")
  }
  if (!all(measures %in% valid_measures)) {
    stop("Invalid measure selected.")
  }

  # Convert AOI if needed
  if (is.character(aoi) && grepl("\\.shp$", aoi)) {
    aoi <- vect(aoi) # Load shapefile as SpatVector
  }

  if (!inherits(aoi, "SpatRaster") && !inherits(aoi, "SpatVector")) {
    stop("AOI must be a SpatRaster, SpatVector, or a shapefile path.")
  }

  raster_list <- list()
  file_paths <- list()

  for (prop in properties) {
    for (depth in depths) {
      for (measure in measures) {
        filename <- paste0(prop, "_", depth, "_", measure, ".tif")
        file_url <- paste0(base_url, filename)

        suffix_str <- ifelse(nchar(suffix) > 0, paste0("_", suffix), "")
        output_file <- if (!is.null(output_dir) && export) {
          file.path(
            output_dir,
            paste0("SOL_", prop, "_", depth, "_", measure, suffix_str, ".tif")
          )
        } else {
          tempfile(fileext = ".tif")
        }

        # Check if the file exists on the server before downloading
        response <- HEAD(file_url)
        if (status_code(response) != 200) {
          message("File does not exist on server: ", file_url, " Skipping...")
          next
        }

        # Download the raster
        tryCatch(
          {
            download.file(file_url, output_file, mode = "wb", timeout = 55)
          },
          error = function(e) {
            message("Download failed for: ", file_url, " | Trying httr::GET()")
            GET(
              file_url,
              write_disk(output_file, overwrite = TRUE),
              timeout(300)
            )
          }
        )

        # Store file path if exporting
        if (export) {
          file_paths <- append(file_paths, output_file)
        }

        # Load raster into R
        raster <- rast(output_file)

        # Reproject AOI to match raster CRS if necessary
        raster_crs <- crs(raster)
        if (crs(aoi) != raster_crs) {
          message("Reprojecting AOI to match raster CRS: ", raster_crs)
          aoi <- project(aoi, raster_crs)
        }

        # Crop the raster to AOI
        cropped_raster <- crop(raster, aoi)
        match_row <- scalar_table[
          scalar_table$property == prop &
            scalar_table$depth == depth &
            scalar_table$measure == measure,
        ]

        if (nrow(match_row) == 1) {
          scalar_value <- match_row$scalar[1]
          cropped_raster <- cropped_raster * scalar_value
        } else {
          warning(
            "No matching scalar found for ",
            prop,
            ", ",
            depth,
            ", ",
            measure
          )
        }

        # Create a temporary file to avoid overwrite error
        temp_file <- paste0(output_file, "_tmp.tif")

        # Save cropped raster to temp file first
        writeRaster(cropped_raster, temp_file, overwrite = TRUE)

        # Rename temp file to final output file
        file.rename(temp_file, output_file)

        # Store raster in list
        raster_list[[paste0(prop, "_", depth, "_", measure)]] <- rast(
          output_file
        )

        message("Processed and saved: ", output_file)
      }
    }
  }

  # Ensure there are downloaded rasters
  if (length(raster_list) == 0) {
    stop("No rasters were successfully downloaded.") # Handle empty case
  }

  # Initialize an empty SpatRaster
  raster_stack <- rast()

  # Append each raster to the stack
  for (r in raster_list) {
    raster_stack <- c(raster_stack, r) # Append raster using `c()`
  }

  # Rename bands based on property and depth
  names(raster_stack) <- names(raster_list)

  # Ensure raster_stack is a valid SpatRaster
  if (nlyr(raster_stack) == 0) {
    stop("Failed to create a valid SpatRaster stack.")
  }

  # Assign "SOL" as an attribute to the SpatRaster
  attr(raster_stack, "product") <- "SOL"

  # Return structured output
  return(list(stack = raster_stack, file_paths = file_paths, product = "SOL"))
}
