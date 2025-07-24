#' Convert zonal summaries of raster soil data into a SoilProfileCollection
#'
#' @param rstack A `SpatRaster` of soil property layers
#' @param zones A polygon `SpatVector` or `sf` object
#' @param stat One or more zonal statistics to compute, e.g. `"mean"`, `"max"`
#' @param id_column Unique identifier column in `zones`
#'
#' @return A `SoilProfileCollection` object with horizons populated from zonal summaries
#' @export
zones_to_SPC <- function(rstack, zones, stat = "mean", id_column = "Name") {
  stopifnot(
    requireNamespace("terra"),
    requireNamespace("dplyr"),
    requireNamespace("tidyr"),
    requireNamespace("stringr"),
    requireNamespace("purrr"),
    requireNamespace("aqp")
  )

  # Lookup tables (user-supplied)
  depth_interval_lookup <- list(
    "0_5" = c("0_5", "0-5cm", "0_cm", "0-5"),
    "5_15" = c("5_15", "5-15cm", "5_cm", "0-25"),
    "15_30" = c("15_30", "15-30cm", "15_cm", "0-25", "25-50"),
    "30_60" = c("30_60", "30-60cm", "30_cm", "30-60", "25-50"),
    "60_100" = c("60_100", "60-100cm", "60_cm", "30-60"),
    "100_200" = c("100_200", "100-200cm", "100_cm", "150_cm")
  )

  depth_range_lookup <- list(
    "0_5" = c(0, 5),
    "5_15" = c(5, 15),
    "15_30" = c(15, 30),
    "30_60" = c(30, 60),
    "60_100" = c(60, 100),
    "100_200" = c(100, 200)
  )

  # Reproject zones if needed
  if (inherits(zones, "sf")) {
    zones <- terra::vect(zones)
  }
  if (!terra::same.crs(rstack, zones)) {
    zones <- terra::project(zones, terra::crs(rstack))
  }

  # Zonal summary
  stat_list <- lapply(stat, function(s) terra::zonal(rstack, zones, fun = s))
  df <- Reduce(function(x, y) dplyr::full_join(x, y, by = id_column), stat_list)

  # Add unique ID
  df$peiid <- df[[id_column]]

  # Convert wide → long
  long_df <- df %>%
    dplyr::select(-any_of(id_column)) %>%
    tidyr::pivot_longer(
      cols = -peiid,
      names_to = "layer",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      # extract depth label using regex
      depth_label = stringr::str_extract(
        layer,
        "[0-9]+[_-][0-9]+(cm)?|[0-9]+_cm|[0-9]+-[0-9]+"
      ),
      # match against lookup
      matched_interval = purrr::map_chr(depth_label, function(lbl) {
        hit <- purrr::keep(
          depth_interval_lookup,
          ~ any(stringr::str_detect(
            lbl,
            paste0("^(", paste(.x, collapse = "|"), ")$")
          ))
        )
        if (length(hit) > 0) names(hit)[1] else NA
      }),
      hzdept = purrr::map_dbl(matched_interval, ~ depth_range_lookup[[.x]][1]),
      hzdepb = purrr::map_dbl(matched_interval, ~ depth_range_lookup[[.x]][2]),
      property = stringr::str_remove(layer, "_[0-9]+.*$") # strip depth suffix
    ) %>%
    dplyr::filter(!is.na(hzdept) & !is.na(property))

  horizon_data <- long_df %>%
    dplyr::select(peiid, hzdept, hzdepb, property, value) %>%
    tidyr::pivot_wider(names_from = property, values_from = value)

  # Construct SoilProfileCollection
  depths(horizon_data) <- peiid ~ hzdept + hzdepb
  site(horizon_data) <- df %>% dplyr::select(peiid)

  return(horizon_data)
}
