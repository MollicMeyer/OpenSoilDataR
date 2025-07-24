#' Convert zonal raster statistics into a SoilProfileCollection
#'
#' Takes a raster stack of soil properties and zonal vector features,
#' computes summary statistics per zone, parses depth intervals, and
#' returns a `SoilProfileCollection` with horizon-level data.
#'
#' @param rstack A `SpatRaster` with soil property layers named using depth suffixes (e.g. `_0_5`, `_5_15`).
#' @param zones An `sf` or `SpatVector` object representing zones/polygons for which zonal stats are calculated.
#' @param stat Summary function to use in `terra::extract()` (e.g., `"mean"`, `"median"`).
#' @param id_column Column name in `zones` used to assign unique profile IDs.
#'
#' @return A `SoilProfileCollection` object with horizon-level attributes for each zone.
#' @import terra
#' @importFrom dplyr select rename filter
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom stringr str_detect str_remove
#' @importFrom aqp depths site
#' @export
zones_to_SPC <- function(rstack, zones, stat = "mean", id_column = "ID") {
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

  if (!inherits(zones, "SpatVector")) {
    zones <- terra::vect(zones)
  }

  if (!terra::same.crs(rstack, zones)) {
    zones <- terra::project(zones, terra::crs(rstack))
  }

  zstats <- terra::extract(rstack, zones, fun = stat, na.rm = TRUE)

  zstats[[id_column]] <- zones[[id_column]]
  df <- zstats %>% dplyr::rename(peiid = !!sym(id_column))

  if ("ID" %in% colnames(df)) {
    colnames(df)[colnames(df) == "ID"] <- "peiid"
  }

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
  site(hz_data) <- df[, "peiid", drop = FALSE]

  return(hz_data)
}
