#' Compute Zonal Statistics for Multiple Soil Products
#'
#' This function calculates zonal statistics for soil property rasters from POLARIS (PSP),
#' SoilGrids v2 (SG2), SOLUS100 (SOL), CSRL Soil props 800, and SSURGO/gNATSGO/RSS (SGO). It ensures that depth intervals from different
#' datasets are consistently translated, and supports depth-weighted averaging.
#'
#' @param soil_data A `SpatRaster` or list of `SpatRaster` objects from `fetch_PSP`, `fetch_SG2`, `fetch_SOL` or `fetch_SGO`.
#' @param tdepth Numeric. Top depth of interest (cm).
#' @param bdepth Numeric. Bottom depth of interest (cm).
#' @param props Character vector. Standardized soil properties (e.g., `"clay"`, `"sand"`, `"bd"`).
#' @param shapes `sf` or `SpatVector`. Polygon boundaries for statistics. **The "Name" field must match plot names.**
#' @param plots Character vector. Names of plots to process (must match `shapes$Name`).
#' @param stats Character vector. Statistics to calculate (`"mean"`, `"min"`, `"max"`, `"sd"`).
#' @param wtd.mean Logical. If `TRUE`, computes depth-weighted mean across specified depth range.
#' @param MakePlot Logical. If `TRUE`, generates raster plots for each statistic.
#' @param output_dir Character. Directory to save output files. Set to `NULL` to disable file saving.
#'
#' #' @section Property Lookup:
#' **Standardized soil properties across datasets.**
#'
#' \preformatted{
#' property_lookup <- list(
#' # Water Retention
#' "wtenthbar"      = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "wtenthbar_r"),
#' "wthirdbar"      = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "wthirdbar_r"),
#' "wfifteenbar"    = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "wfifteenbar_r"),
#' "awc"            = c(PSP = NA, SG2 = NA, SOL = "awc_r", CSRL = NA, SGO = "awc_r"),
#'
#' # Phosphorus
#' "pbray1"         = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "pbray1_r"),
#' "ptotal"         = c(PSP = NA, SG2 = NA, SOL = "ptotal_r", CSRL = NA, SGO = "ptotal_r"),
#' "ph2osoluble"    = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "ph2osoluble_r"),
#' "poxalate"       = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "poxalate_r"),
#'
#' # Bulk Density
#' "bd"             = c(PSP = "bd_mean", SG2 = "bdod", SOL = "dbovendry_r", CSRL = "bulk_density", SGO = "dbovendry_r"),
#' "bd_tenthbar"    = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "dbtenthbar_r"),
#' "bd_thirdbar"    = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "dbthirdbar_r"),
#' "bd_fifteenbar"  = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "dbfifteenbar_r"),
#'
#' # Texture
#' "clay"           = c(PSP = "clay_mean", SG2 = "clay", SOL = "claytotal", CSRL = "clay_profile", SGO = "claytotal_r"),
#' "sand"           = c(PSP = "sand_mean", SG2 = "sand", SOL = "sandtotal", CSRL = "sand_profile", SGO = "sandtotal_r"),
#' "silt"           = c(PSP = "silt_mean", SG2 = "silt", SOL = "silttotal", CSRL = "silt_profile", SGO = "silttotal_r"),
#' "sandco"         = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "sandco_r"),
#' "sandfine"       = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "sandfine_r"),
#' "sandmed"        = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "sandmed_r"),
#' "sandvc"         = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "sandvc_r"),
#' "sandvf"         = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "sandvf_r"),
#' "siltco"         = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "siltco_r"),
#' "siltfine"       = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "siltfine_r"),
#'
#' # Organic Carbon
#' "som"            = c(PSP = "om_mean", SG2 = NA, SOL = NA, CSRL = "som_max", SGO = "om_r"),
#' "soc"            = c(PSP = "soc", SG2 = "soc", SOL = "soc", CSRL = "soc_max", SGO = "oc_r"),
#' "soc_stock"      = c(PSP = NA, SG2 = "ocs", SOL = NA, CSRL = "soc_stock", SGO = NA),
#'
#' # Chemistry
#' "ph"             = c(PSP = "ph_mean", SG2 = "phh2o", SOL = "ph1to1h2o", CSRL = "ph_profile", SGO = "ph1to1h2o_r"),
#' "ph_cacl2"       = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "ph01mcacl2_r"),
#' "ph_oxidized"    = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "phoxidized_r"),
#' "cec"            = c(PSP = NA, SG2 = "cec", SOL = "cec7", CSRL = "cec_profile", SGO = "cec7_r"),
#' "ecec"           = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "ecec_r"),
#' "ec"             = c(PSP = NA, SG2 = NA, SOL = "ec", CSRL = "ec_profile", SGO = "ec_r"),
#' "ec15"           = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "ec15_r"),
#' "esp"            = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "esp_r"),
#' "sar"            = c(PSP = NA, SG2 = NA, SOL = "sar", CSRL = "sar", SGO = "sar_r"),
#' "gypsum"         = c(PSP = NA, SG2 = NA, SOL = "gypsum", CSRL = NA, SGO = "gypsum_r"),
#' "caco3"          = c(PSP = NA, SG2 = NA, SOL = "caco3", CSRL = "caco3", SGO = "caco3_r"),
#'
#' # Extractables
#' "extracid"       = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "extracid_r"),
#' "extral"         = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "extral_r"),
#' "aloxalate"      = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "aloxalate_r"),
#' "feoxalate"      = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "feoxalate_r"),
#'
#' # Mechanics
#' "lep"            = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "lep_r"),
#' "ll"             = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "ll_r"),
#' "pi"             = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "pi_r"),
#'
#' # Rock Fragments
#' "frag3to10"      = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "frag3to10_r"),
#' "fraggt10"       = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "fraggt10_r"),
#'
#' # Fiber
#' "fiberrubbedpct" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "fiberrubbedpct_r"),
#' "fiberunrubbedpct" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "fiberunrubbedpct_r"),
#'
#' # Sieve
#' "sieveno10"      = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "sieveno10_r"),
#' "sieveno200"     = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "sieveno200_r"),
#' "sieveno4"       = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "sieveno4_r"),
#' "sieveno40"      = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "sieveno40_r"),
#'
#' # Hydrology
#' "ksat"           = c(PSP = "ksat_mean", SG2 = NA, SOL = NA, CSRL = "ksat_max", SGO = "ksat_r"),
#' "wsatiated"      = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "wsatiated_r"),

#' # Erosion
#' "kw"             = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "kwfact"),
#' "kffact"         = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "kffact"),
#' "kifact"         = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "kifact"),
#' "krfact"         = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "krfact")
#' )
#' }
#'
#' @return A list containing:
#'   - `$Unweighted`: Data frame of zonal statistics per depth interval.
#'   - `$Weighted`: Data frame of depth-weighted means (if `wtd.mean = TRUE`).
#'
#' @import terra sf dplyr tidyr
#' @export
#'
#' @examples
#' \dontrun{
#' # Load required libraries
#' library(terra)
#' library(sf)
#'
#' # Load polygon shapefile
#' plot_shapes <- st_read("path/to/plots.shp")
#'
#' # Fetch PSP soil data
#' psp_data <- fetch_PSP(aoi = plot_shapes, properties = c("clay_mean", "sand_mean"), depths = c("0_5", "5_15", "15_30"))
#'
#' # Compute zonal statistics
#' results <- s.zonalstats(
#'   soil_data = psp_data,
#'   tdepth = 0, bdepth = 20,
#'   props = c("clay", "sand"),
#'   shapes = plot_shapes,
#'   plots = c("NW_plot", "SW_plot"),
#'   stats = c("mean", "min", "max", "sd"),
#'   wtd.mean = TRUE,
#'   MakePlot = TRUE,
#'   output_dir = "path/to/output"
#' )
#'
#' # View results
#' head(results$Unweighted)
#' head(results$Weighted)
#' }

s.zonalstats <- function(
  soil_data,
  tdepth = 0,
  bdepth = 20,
  props = c("sand", "clay"),
  shapes = plots,
  plots = c("NW_plot", "SW_plot"),
  stats = c("mean", "min", "max", "sd"),
  wtd.mean = TRUE,
  MakePlot = FALSE,
  output_dir = NULL
) {
  if (MakePlot) {
    cat(
      "You've selected the MakePlot = TRUE, this will slow processing a tad bit. Set to MakePlot = FALSE if you do not need plots. \n"
    )
  }

  # Property Lookup Table for use in s.zonalstats()
  property_lookup <- list(
    # Water Retention
    "wtenthbar" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "wtenthbar_r"
    ),
    "wthirdbar" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "wthirdbar_r"
    ),
    "wfifteenbar" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "wfifteenbar_r"
    ),
    "awc" = c(PSP = NA, SG2 = NA, SOL = "awc_r", CSRL = NA, SGO = "awc_r"),

    # Phosphorus
    "pbray1" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "pbray1_r"),
    "ptotal" = c(
      PSP = NA,
      SG2 = NA,
      SOL = "ptotal_r",
      CSRL = NA,
      SGO = "ptotal_r"
    ),
    "ph2osoluble" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "ph2osoluble_r"
    ),
    "poxalate" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "poxalate_r"),

    # Bulk Density
    "bd" = c(
      PSP = "bd_mean",
      SG2 = "bdod",
      SOL = "dbovendry_r",
      CSRL = "bulk_density",
      SGO = "dbovendry_r"
    ),
    "bd_tenthbar" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "dbtenthbar_r"
    ),
    "bd_thirdbar" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "dbthirdbar_r"
    ),
    "bd_fifteenbar" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "dbfifteenbar_r"
    ),

    # Texture
    "clay" = c(
      PSP = "clay_mean",
      SG2 = "clay",
      SOL = "claytotal",
      CSRL = "clay_profile",
      SGO = "claytotal_r"
    ),
    "sand" = c(
      PSP = "sand_mean",
      SG2 = "sand",
      SOL = "sandtotal",
      CSRL = "sand_profile",
      SGO = "sandtotal_r"
    ),
    "silt" = c(
      PSP = "silt_mean",
      SG2 = "silt",
      SOL = "silttotal",
      CSRL = "silt_profile",
      SGO = "silttotal_r"
    ),
    "sandco" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "sandco_r"),
    "sandfine" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "sandfine_r"),
    "sandmed" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "sandmed_r"),
    "sandvc" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "sandvc_r"),
    "sandvf" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "sandvf_r"),
    "siltco" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "siltco_r"),
    "siltfine" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "siltfine_r"),

    # Organic Carbon
    "som" = c(
      PSP = "om_mean",
      SG2 = NA,
      SOL = NA,
      CSRL = "som_max",
      SGO = "om_r"
    ),
    "soc" = c(
      PSP = "soc",
      SG2 = "soc",
      SOL = "soc",
      CSRL = "soc_max",
      SGO = "oc_r"
    ),
    "soc_stock" = c(
      PSP = NA,
      SG2 = "ocs",
      SOL = NA,
      CSRL = "soc_stock",
      SGO = NA
    ),

    # Chemistry
    "ph" = c(
      PSP = "ph_mean",
      SG2 = "phh2o",
      SOL = "ph1to1h2o",
      CSRL = "ph_profile",
      SGO = "ph1to1h2o_r"
    ),
    "ph_cacl2" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "ph01mcacl2_r"
    ),
    "ph_oxidized" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "phoxidized_r"
    ),
    "cec" = c(
      PSP = NA,
      SG2 = "cec",
      SOL = "cec7",
      CSRL = "cec_profile",
      SGO = "cec7_r"
    ),
    "ecec" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "ecec_r"),
    "ec" = c(PSP = NA, SG2 = NA, SOL = "ec", CSRL = "ec_profile", SGO = "ec_r"),
    "ec15" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "ec15_r"),
    "esp" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "esp_r"),
    "sar" = c(PSP = NA, SG2 = NA, SOL = "sar", CSRL = "sar", SGO = "sar_r"),
    "gypsum" = c(
      PSP = NA,
      SG2 = NA,
      SOL = "gypsum",
      CSRL = NA,
      SGO = "gypsum_r"
    ),
    "caco3" = c(
      PSP = NA,
      SG2 = NA,
      SOL = "caco3",
      CSRL = "caco3",
      SGO = "caco3_r"
    ),

    # Extractables
    "extracid" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "extracid_r"),
    "extral" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "extral_r"),
    "aloxalate" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "aloxalate_r"
    ),
    "feoxalate" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "feoxalate_r"
    ),

    # Mechanics
    "lep" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "lep_r"),
    "ll" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "ll_r"),
    "pi" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "pi_r"),

    # Rock Fragments
    "frag3to10" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "frag3to10_r"
    ),
    "fraggt10" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "fraggt10_r"),

    # Fiber
    "fiberrubbedpct" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "fiberrubbedpct_r"
    ),
    "fiberunrubbedpct" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "fiberunrubbedpct_r"
    ),

    # Sieve
    "sieveno10" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "sieveno10_r"
    ),
    "sieveno200" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "sieveno200_r"
    ),
    "sieveno4" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "sieveno4_r"),
    "sieveno40" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "sieveno40_r"
    ),

    # Hydrology
    "ksat" = c(
      PSP = "ksat_mean",
      SG2 = NA,
      SOL = NA,
      CSRL = "ksat_max",
      SGO = "ksat_r"
    ),
    "wsatiated" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "wsatiated_r"
    ),

    # Erosion
    "kw" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "kwfact"),
    "kffact" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "kffact"),
    "kifact" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "kifact"),
    "krfact" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "krfact")
  )

  # **Depth Interval Lookup Table** - Standardizing depth labels across datasets
  depth_interval_lookup <- list(
    "0_5" = c("0_5", "0-5cm", "0_cm", "0-5"),
    "5_15" = c("5_15", "5-15cm", "5_cm", "0-25"),
    "15_30" = c("15_30", "15-30cm", "15_cm", "0-25", "25-50"),
    "30_60" = c("30_60", "30-60cm", "30_cm", "30-60", "25-50"),
    "60_100" = c("60_100", "60-100cm", "60_cm", "30-60"),
    "100_200" = c("100_200", "100-200cm", "100_cm", "150_cm")
  )

  # **Depth Range Lookup Table** - Assigning numeric depths
  depth_range_lookup <- list(
    "0_5" = c(0, 5),
    "5_15" = c(5, 15),
    "15_30" = c(15, 30),
    "30_60" = c(30, 60),
    "60_100" = c(60, 100),
    "100_200" = c(100, 200)
  )

  # **Detect Dataset Types (PSP, SG2, SOL)**
  datasets <- list()
  if (inherits(soil_data, "SpatRaster")) {
    dataset_name <- attr(soil_data, "product")
    if (is.null(dataset_name)) {
      stop(
        "Dataset type unknown. Ensure data comes from fetch_PSP, fetch_SG2, fetch_SOL, fetch_CSRL, or fetch_SGO."
      )
    }
    datasets[[dataset_name]] <- soil_data
  } else if (is.list(soil_data)) {
    for (raster_stack in soil_data) {
      dataset_name <- attr(raster_stack, "product")
      if (is.null(dataset_name)) {
        stop(
          "Dataset type unknown. Ensure data comes from fetch_PSP, fetch_SG2, or fetch_SOL."
        )
      }
      datasets[[dataset_name]] <- raster_stack
    }
  } else {
    stop(
      "Invalid input: Provide a SpatRaster object or a list of SpatRasters from fetch_PSP, fetch_SG2, or fetch_SOL."
    )
  }

  # **Reproject rasters to match zones CRS**
  zones_crs <- crs(shapes)
  for (dataset_name in names(datasets)) {
    if (crs(datasets[[dataset_name]]) != zones_crs) {
      datasets[[dataset_name]] <- project(datasets[[dataset_name]], zones_crs)
    }
  }

  # **Initialize Results DataFrame**
  results_df <- data.frame()
  weighted_mean_df <- data.frame()

  for (dataset_name in names(datasets)) {
    cat("\nProcessing dataset:", dataset_name, "\n")

    soil_rasters <- datasets[[dataset_name]]

    for (plot in plots) {
      if (!plot %in% shapes$Name) {
        warning(paste("Plot", plot, "not found in shapes$Name"))
        next
      }
      poly <- shapes[shapes$Name == plot, ]
      if (nrow(poly) == 0 || any(is.na(sf::st_geometry(poly)))) {
        warning(paste("No valid geometry found for plot:", plot))
        next
      }
      if (!all(sf::st_geometry_type(poly) %in% c("POLYGON", "MULTIPOLYGON"))) {
        warning(paste("Non-polygon geometry detected for plot:", plot))
      }
      poly <- terra::vect(poly)

      for (prop in props) {
        dataset_prop <- property_lookup[[prop]][[dataset_name]]
        if (is.na(dataset_prop)) {
          next
        }

        # **Match raster layers based on property & depth**
        matched_layers <- names(soil_rasters)[grepl(
          dataset_prop,
          names(soil_rasters)
        )]

        if (length(matched_layers) == 0) {
          warning(paste(
            "No raster files found for",
            dataset_prop,
            "in",
            dataset_name,
            ". Skipping..."
          ))
          next
        }

        selected_layers <- list()
        depth_weights <- list()

        for (raster_name in matched_layers) {
          matched_depth <- NA

          # **Extract depth from filename using interval lookup**
          for (depth_name in names(depth_interval_lookup)) {
            if (
              any(grepl(
                paste(depth_interval_lookup[[depth_name]], collapse = "|"),
                raster_name
              ))
            ) {
              matched_depth <- depth_name
              break
            }
          }
          if (is.na(matched_depth)) {
            next
          }

          layer_depths <- depth_range_lookup[[matched_depth]]
          overlap <- pmin(bdepth, layer_depths[2]) -
            pmax(tdepth, layer_depths[1])

          if (overlap > 0) {
            selected_layers[[raster_name]] <- soil_rasters[[raster_name]]
            depth_weights[[raster_name]] <- overlap /
              (layer_depths[2] - layer_depths[1])
          }

          raster_layer <- soil_rasters[[raster_name]]

          for (stat in stats) {
            zonal_stats <- zonal(raster_layer, poly, stat)

            # Store per-depth statistics
            results_df <- rbind(
              results_df,
              data.frame(
                Product = dataset_name,
                SoilProperty = prop,
                Depth = matched_depth,
                Statistic = stat,
                Plot = plot,
                ZonalStats = round(zonal_stats[[1]], 2)
              )
            )
          }
        }

        # **Compute Weighted Mean Across Depths**
        if (wtd.mean && length(selected_layers) > 0) {
          selected_rasters <- rast(selected_layers)
          weights <- unlist(depth_weights)
          weighted_raster <- sum(selected_rasters * weights) / sum(weights)

          for (stat in stats) {
            zonal_stats <- zonal(weighted_raster, poly, stat)

            # Store weighted mean statistics separately
            weighted_mean_df <- rbind(
              weighted_mean_df,
              data.frame(
                Product = dataset_name,
                SoilProperty = prop,
                Depth = paste0(tdepth, "_", bdepth), # Weighted depth interval
                Statistic = stat,
                Plot = plot,
                ZonalStats = round(zonal_stats[[1]], 2)
              )
            )

            # ✅ Plot raster if enabled
            if (MakePlot) {
              options(terra.pal = NULL)
              buffer_size = 50
              bbox <- st_bbox(shapes) +
                c(-buffer_size, -buffer_size, buffer_size, buffer_size)

              plot(
                raster_layer,
                col = terrain.colors(10),
                main = paste(
                  "Map:",
                  dataset_name,
                  "| Plot:",
                  plot,
                  "| Property:",
                  prop,
                  "|",
                  stat,
                  ":",
                  round(zonal_stats[[1]], 2)
                ),
                ext = bbox
              )
              plot(st_geometry(shapes), add = TRUE, border = "black", lwd = 2)
              plot(poly, add = TRUE, border = "red", lwd = 2.5)
            }
          }
        }
      }
    }
  }

  # **Save Results**
  if (!is.null(output_dir)) {
    write.csv(
      results_df,
      file = file.path(output_dir, "DSM_zonal_stats_merged.csv"),
      row.names = FALSE
    )
    if (wtd.mean) {
      write.csv(
        weighted_mean_df,
        file = file.path(output_dir, "DSM_zonal_stats_weighted.csv"),
        row.names = FALSE
      )
    }
  }

  # **Return both datasets**
  return(list(
    Unweighted = results_df,
    Weighted = if (wtd.mean) weighted_mean_df else NULL
  ))
}
