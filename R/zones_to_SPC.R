#' Convert zonal raster statistics into a SoilProfileCollection
#'
#' This function extracts zonal statistics from a raster stack using polygons (zones),
#' reshapes the data into horizon format, and converts it into an `aqp::SoilProfileCollection`.
#'
#' @param rstack A `SpatRaster` containing stacked soil property rasters.
#' @param zones A `SpatVector` or `sf` object of polygons used for zonal summaries.
#' @param stat A summary statistic to apply per zone (e.g., "mean", "median"). Default: "mean".
#' @param id_column Character. Name of the column in `zones` to use as profile IDs.
#'
#' @return A `SoilProfileCollection` with horizon data extracted from the raster stack.
#' @export
zones_to_SPC <- function(rstack, zones, stat = "mean", id_column = "Name") {
  require(terra)
  require(dplyr)
  require(tidyr)
  require(stringr)
  require(aqp)

  # Depth lookup tables
  depth_interval_lookup <- list(
    "0_5" = c("0_5", "0-5cm", "_0_cm_p", "0-5"),
    "5_15" = c("5_15", "5-15cm", "_5_cm_p", "0-25"),
    "15_30" = c("15_30", "15-30cm", "_15_cm_p", "0-25", "25-50"),
    "30_60" = c("30_60", "30-60cm", "_30_cm_p", "30-60", "25-50"),
    "60_100" = c("60_100", "60-100cm", "_60_cm_p", "30-60"),
    "100_200" = c("100_200", "100-200cm", "_100_cm_p", "_150_cm_p")
  )

  depth_range_lookup <- list(
    "0_5" = c(0, 5),
    "5_15" = c(5, 15),
    "15_30" = c(15, 30),
    "30_60" = c(30, 60),
    "60_100" = c(60, 100),
    "100_200" = c(100, 200)
  )

  # CRS alignment
  if (!terra::same.crs(rstack, zones)) {
    zones <- if (inherits(zones, "SpatVector")) {
      terra::project(zones, terra::crs(rstack))
    } else {
      sf::st_transform(zones, terra::crs(rstack))
    }
  }

  zones_vect <- if (inherits(zones, "SpatVector")) zones else terra::vect(zones)

  # Extract zonal stats
  zstats <- terra::extract(rstack, zones_vect, fun = stat, na.rm = TRUE)

  # Add unique profile ID
  zone_ids <- as.data.frame(zones)[[id_column]]
  zstats$peiid <- zone_ids

  if ("ID" %in% names(zstats)) {
    zstats <- zstats %>% select(-ID)
  }

  df <- zstats %>% relocate(peiid)

  long_df <- df %>%
    pivot_longer(cols = -peiid, names_to = "layer", values_to = "value") %>%
    rowwise() %>%
    mutate(
      matched_label = {
        matched <- NA_character_
        for (key in names(depth_interval_lookup)) {
          for (val in depth_interval_lookup[[key]]) {
            if (str_detect(layer, fixed(val))) {
              matched <- key
              break
            }
          }
          if (!is.na(matched)) break
        }
        matched
      },
      hzdept = if (!is.na(matched_label)) {
        depth_range_lookup[[matched_label]][1]
      } else {
        NA_real_
      },
      hzdepb = if (!is.na(matched_label)) {
        depth_range_lookup[[matched_label]][2]
      } else {
        NA_real_
      },
      matched_string = if (!is.na(matched_label)) {
        match_found <- NA_character_
        for (val in depth_interval_lookup[[matched_label]]) {
          if (str_detect(layer, fixed(val))) {
            match_found <- val
            break
          }
        }
        match_found
      } else {
        NA_character_
      },
      variable = if (!is.na(matched_string)) {
        var <- str_replace(layer, fixed(matched_string), "")
        var <- str_replace_all(var, "__+", "_")
        var <- str_remove_all(var, "^_|_$")
        var
      } else {
        NA_character_
      }
    ) %>%
    ungroup() %>%
    filter(!is.na(hzdept), !is.na(hzdepb), !is.na(variable)) %>%
    mutate(value = as.numeric(value)) %>%
    group_by(peiid, hzdept, hzdepb, variable) %>%
    summarize(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = variable, values_from = value)

  site_meta <- as.data.frame(zones)[, id_column, drop = FALSE]
  colnames(site_meta)[1] <- "peiid"
  site(long_df) <- left_join(site(long_df), site_meta, by = "peiid")

  return(long_df)
}
