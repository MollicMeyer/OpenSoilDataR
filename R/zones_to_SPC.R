#' Convert Zonal Summary Raster Data to SoilProfileCollection
#'
#' @param rstack A `SpatRaster` stack from a fetch_*() function
#' @param zones A polygon `sf` or `SpatVector` object defining zones
#' @param stat Character. One or more statistics passed to `terra::zonal()` (e.g., "mean", "min")
#' @param id_column Character. Name of the unique identifier column in zones (default: "Name")
#' @return A SoilProfileCollection with one profile per zone
#' @export
zones_to_SPC <- function(rstack, zones, stat = "mean", id_column = "Name") {
  require(terra)
  require(dplyr)
  require(tidyr)
  require(aqp)
  require(stringr)
  require(sf)

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

  parse_layers <- function(df, idcol = "peiid") {
    long_df <- df %>%
      pivot_longer(
        cols = -all_of(idcol),
        names_to = "layer",
        values_to = "value"
      ) %>%
      separate(
        layer,
        into = c("property", "top", "bottom"),
        sep = "_",
        remove = FALSE
      ) %>%
      mutate(
        depth_label = paste0(top, "_", bottom),
        matched_label = names(depth_interval_lookup)[
          sapply(depth_interval_lookup, function(p) {
            any(str_detect(
              depth_label,
              paste0("^(", paste(p, collapse = "|"), ")$")
            ))
          })
        ],
        hzdept = depth_range_lookup[[matched_label]][1],
        hzdepb = depth_range_lookup[[matched_label]][2]
      ) %>%
      ungroup()

    long_df %>%
      select(all_of(idcol), hzdept, hzdepb, property, value) %>%
      pivot_wider(names_from = property, values_from = value)
  }

  if (!inherits(zones, "SpatVector")) {
    zones <- terra::vect(zones)
  }
  if (!id_column %in% names(zones)) {
    stop("id_column not found in zones.")
  }
  if (!same.crs(rstack, zones)) {
    zones <- project(zones, crs(rstack))
  }

  stat_list <- lapply(stat, function(s) terra::zonal(rstack, zones, fun = s))
  df <- Reduce(function(x, y) full_join(x, y, by = id_column), stat_list)
  df$peiid <- df[[id_column]]

  hz <- parse_layers(df, idcol = "peiid")
  depths(hz) <- peiid ~ hzdept + hzdepb
  site(hz) <- df %>% select(peiid)

  return(hz)
}
