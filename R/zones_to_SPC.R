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
  library(terra)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(aqp)

  # Lookup tables
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
  if (!terra::same.crs(rstack, zones)) {
    zones <- sf::st_transform(zones, terra::crs(rstack))
  }

  zones_vect <- terra::vect(zones)
  zstats <- terra::extract(rstack, zones_vect, fun = stat, na.rm = TRUE)

  # Add ID from zones to extracted
  zstats[[id_column]] <- zones[[id_column]]

  # Rename ID column to peiid
  df <- zstats %>% dplyr::rename(peiid = !!sym(id_column))

  # Ensure ID column doesn't conflict
  if ("ID" %in% names(df)) {
    df <- df %>% dplyr::select(-ID)
  }

  # Pivot longer to depth layers
  long_df <- df %>%
    pivot_longer(cols = -peiid, names_to = "layer", values_to = "value") %>%
    mutate(
      matched_label = sapply(layer, function(lab) {
        match <- NA_character_
        for (label in names(depth_interval_lookup)) {
          if (
            any(str_detect(
              lab,
              paste0(
                "(",
                paste(depth_interval_lookup[[label]], collapse = "|"),
                ")$"
              )
            ))
          ) {
            match <- label
            break
          }
        }
        match
      }),
      hzdept = ifelse(
        !is.na(matched_label),
        sapply(matched_label, function(m) depth_range_lookup[[m]][1]),
        NA_real_
      ),
      hzdepb = ifelse(
        !is.na(matched_label),
        sapply(matched_label, function(m) depth_range_lookup[[m]][2]),
        NA_real_
      ),
      variable = ifelse(
        !is.na(matched_label),
        mapply(
          function(lab, label) {
            str_remove(
              lab,
              paste0(
                "_(",
                paste(depth_interval_lookup[[label]], collapse = "|"),
                ")$"
              )
            )
          },
          layer,
          matched_label
        ),
        NA_character_
      )
    ) %>%
    filter(!is.na(matched_label) & !is.na(variable))

  hz_data <- long_df %>%
    select(peiid, hzdept, hzdepb, variable, value) %>%
    pivot_wider(names_from = variable, values_from = value)

  depths(hz_data) <- peiid ~ hzdept + hzdepb

  # Add site metadata from original zones
  site_meta <- as.data.frame(zones)[, id_column, drop = FALSE]
  colnames(site_meta)[1] <- "peiid"
  site(hz_data) <- left_join(site(hz_data), site_meta, by = "peiid")

  return(hz_data)
}
