#' Fetch Soil Property Data from Soil Data Access SSURGO (SGO)
#'
#' Downloads soil property data from the USDA-NRCS Soil Data Access (SDA) service.
#' Supports only **horizon-level** (depth-dependent) properties. Uses `mukey.wcs()` to
#' retrieve mapunit raster and joins queried soil property values to raster via `mukey`.
#'
#' @param aoi A `SpatRaster`, `SpatVector`, or shapefile path representing the area of interest.
#' @param properties Character vector. Horizon-level property names (e.g. `"cec7_r"`, `"claytotal_r"`, `"om_r"`).
#' @param depths A list of numeric depth intervals (e.g., `list(c(0,5), c(5,15))`).
#' @param method Character. Method to use: `"Dominant Component (Numeric)"` (default), `"Weighted Average"`, or `"Min/Max"`.
#' @param tosoc Logical. If `TRUE`, converts `"om_r"` to `"oc_r"` using the Van Bemmelen factor (1.724).
#' @param crs Character. CRS for raster output. Default is `"EPSG:4326"`.
#' @param res Numeric. Raster resolution in meters. Default is `30`.
#' @param db Character. SDA raster database to query: `"gssurgo"` (default), `"gnatsgo"`, `"RSS"`, or `"statsgo"`.
#' @param export Logical. If `TRUE`, saves GeoTIFFs to disk.
#' @param output_dir Character. Directory to save exported files (required if `export = TRUE`).
#' @param suffix Optional character string to append to output filenames.
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{stack}: SpatRaster of queried layers
#'     \item \code{file_paths}: Vector of written file paths (if `export = TRUE`)
#'     \item \code{mukey_raster}: The original mukey raster (SpatRaster)
#'   }
#'
#' @import terra soilDB
#' @export
#'
#' @section Horizon-Level Properties (Require Depth):
#' \tabular{ll}{
#' \strong{Property} \tab \strong{Column} \cr
#' Available Water Capacity \tab awc_r \cr
#' Bulk Density (0.1 bar) \tab dbtenthbar_r \cr
#' Bulk Density (Oven Dry) \tab dbovendry_r \cr
#' Calcium Carbonate \tab caco3_r \cr
#' Cation Exchange Capacity \tab cec7_r \cr
#' Coarse Sand \tab sandco_r \cr
#' Electrical Conductivity \tab ec_r \cr
#' Gypsum \tab gypsum_r \cr
#' Organic Matter \tab om_r \cr
#' pH (1:1 H2O) \tab ph1to1h2o_r \cr
#' Saturated Hydraulic Conductivity \tab ksat_r \cr
#' Total Clay \tab claytotal_r \cr
#' Total Sand \tab sandtotal_r \cr
#' Total Silt \tab silttotal_r \cr
#' And more (see SDA documentation) \tab ...
#' }
#'
#' @examples
#' \dontrun{
#' aoi <- rast("your_dem.tif")
#' sgo <- fetch_SGO(
#'   aoi = aoi,
#'   properties = c("ph1to1h2o_r", "cec7_r", "om_r"),
#'   depths = list(c(0,5), c(5,15)),
#'   method = "Dominant Component (Numeric)",
#'   tosoc = TRUE,
#'   crs = "EPSG:4326",
#'   res = 30,
#'   db = "gssurgo",
#'   export = TRUE,
#'   output_dir = "SGO_outputs",
#'   suffix = "example"
#' )
#' }

fetch_SGO <- function(
  aoi,
  properties,
  depths = list(
    c(0, 5),
    c(5, 15),
    c(15, 30),
    c(30, 60),
    c(60, 100),
    c(100, 200)
  ),
  method = "Dominant Component (Numeric)",
  tosoc = FALSE,
  crs = "EPSG:4326",
  res = 30,
  db = "gssurgo",
  export = FALSE,
  output_dir = NULL,
  suffix = ""
) {
  stopifnot(requireNamespace("terra"))
  stopifnot(requireNamespace("soilDB"))

  # Supported numeric methods
  valid_methods <- c(
    "Dominant Component (Numeric)",
    "Weighted Average",
    "Min/Max"
  )
  if (!method %in% valid_methods) {
    stop("Invalid method: choose from ", paste(valid_methods, collapse = ", "))
  }

  # Valid horizon-level properties
  valid_props <- c(
    "wtenthbar_r",
    "wthirdbar_r",
    "wfifteenbar_r",
    "awc_r",
    "pbray1_r",
    "dbtenthbar_r",
    "dbthirdbar_r",
    "dbfifteenbar_r",
    "dbovendry_r",
    "claysizedcarb_r",
    "caco3_r",
    "cec7_r",
    "sandco_r",
    "siltco_r",
    "ecec_r",
    "ec15_r",
    "ec_r",
    "esp_r",
    "extral_r",
    "extracid_r",
    "sandfine_r",
    "siltfine_r",
    "freeiron_r",
    "gypsum_r",
    "kffact",
    "kifact",
    "krfact",
    "kwfact",
    "lep_r",
    "ll_r",
    "sandmed_r",
    "om_r",
    "aloxalate_r",
    "feoxalate_r",
    "poxalate_r",
    "pi_r",
    "frag3to10_r",
    "fraggt10_r",
    "fiberrubbedpct_r",
    "wsatiated_r",
    "ksat_r",
    "sar_r",
    "sumbases_r",
    "claytotal_r",
    "ptotal_r",
    "sandtotal_r",
    "silttotal_r",
    "fiberunrubbedpct_r",
    "sandvc_r",
    "sandvf_r",
    "ph2osoluble_r",
    "sieveno10_r",
    "sieveno200_r",
    "sieveno4_r",
    "sieveno40_r",
    "ph01mcacl2_r",
    "ph1to1h2o_r",
    "phoxidized_r"
  )

  if (any(!properties %in% valid_props)) {
    stop(
      "Invalid property(ies): ",
      paste(properties[!properties %in% valid_props], collapse = ", ")
    )
  }

  # Get mukey raster
  cat("Querying mukey raster from:", db, "\n")
  mu_raster <- soilDB::mukey.wcs(aoi = aoi, db = db, res = res)
  mukeys <- unique(terra::values(mu_raster))
  mukeys <- mukeys[!is.na(mukeys)]

  raster_stack <- list()
  file_paths <- c()

  for (prop in properties) {
    for (depth in depths) {
      top <- depth[1]
      bottom <- depth[2]

      query_prop <- prop
      out_prop <- prop

      if (tosoc && prop == "om_r") {
        query_prop <- "om_r" # Always query om_r
        out_prop <- "oc_r" # Rename output to oc_r
      }

      cat("Querying", query_prop, "|", top, "-", bottom, "cm...\n")
      df <- get_SDA_property(
        property = query_prop,
        method = method,
        mukeys = mukeys,
        top_depth = top,
        bottom_depth = bottom,
        include_minors = FALSE
      )[, c("mukey", query_prop)]

      lyr <- terra::classify(mu_raster, df, others = NA)
      nm <- paste0(out_prop, "_", top, "_", bottom)
      names(lyr) <- nm

      # Convert OM to SOC
      if (tosoc && prop == "om_r") {
        lyr <- lyr / 1.724
        names(lyr) <- nm
        message("Converted 'om_r' to 'oc_r' using Van Bemmelen factor (1.724).")
      }

      # Export if enabled
      if (export && !is.null(output_dir)) {
        suffix_str <- ifelse(nzchar(suffix), paste0("_", suffix), "")
        out_file <- file.path(
          output_dir,
          paste0("SGO_", nm, suffix_str, ".tif")
        )
        terra::writeRaster(lyr, out_file, overwrite = TRUE)
        file_paths[nm] <- out_file
      }

      raster_stack[[nm]] <- lyr
    }
  }

  final_stack <- terra::rast(raster_stack)
  attr(final_stack, "product") <- "SGO"

  # Reproject to user-specified CRS
  if (!is.null(crs)) {
    target_crs <- terra::crs(crs, proj = TRUE)
    final_stack <- terra::project(final_stack, target_crs, method = "bilinear")
  }

  return(list(stack = final_stack, file_paths = file_paths, product = "SGO"))
}
