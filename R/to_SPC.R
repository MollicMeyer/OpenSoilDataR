#' Convert Raster, Point, or Zonal Soil Data to SoilProfileCollection
#'
#' @param rstack A `SpatRaster` (or list of) from `fetch_*()` functions
#' @param mode One of "raster", "points", "zonal"
#' @param locations Optional `sf` or `SpatVector` object containing points or polygons. Must include a unique identifier column (e.g., "Name").
#' @param stat For "zonal" mode, one or more of "mean", "min", "max", "median"
#' @param id_column ID column in `locations` (default: "Name")
#' @return A `SoilProfileCollection` object
#' @export
to_spc <- function(
  rstack,
  mode = "raster",
  locations = NULL,
  stat = "mean",
  id_column = "Name"
) {
  require(terra)
  require(dplyr)
  require(tidyr)
  require(aqp)
  require(sf)
  require(stringr)

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
        into = c("property", "mid", "top", "bottom"),
        sep = "_",
        remove = FALSE
      ) %>%
      mutate(
        depth_label = paste0(top, "_", bottom),
        variable = paste(property, mid, sep = "_")
      ) %>%
      rowwise() %>%
      mutate(
        matched_label = {
          matches <- sapply(depth_interval_lookup, function(patterns) {
            any(str_detect(
              depth_label,
              paste0("^(", paste(patterns, collapse = "|"), ")$")
            ))
          })
          names(depth_interval_lookup)[which(matches)[1]]
        },
        hzdept = depth_range_lookup[[matched_label]][1],
        hzdepb = depth_range_lookup[[matched_label]][2]
      ) %>%
      ungroup()

    long_df %>%
      select(all_of(idcol), hzdept, hzdepb, variable, value) %>%
      pivot_wider(names_from = variable, values_from = value)
  }

  if (mode == "raster") {
    df <- as.data.frame(rstack, xy = TRUE, cells = TRUE, na.rm = TRUE)
    df$peiid <- paste0("cell_", df$cell)
    site_data <- df %>% select(peiid, x, y)
    horizon_data <- df %>% select(-x, -y, -cell)
    hz <- parse_layers(horizon_data)
    depths(hz) <- peiid ~ hzdept + hzdepb
    site(hz) <- site_data
    return(hz)
  }

  if (mode == "points") {
    stopifnot(!is.null(locations))
    if (!inherits(locations, "SpatVector")) {
      locations <- terra::vect(locations)
    }
    values <- terra::extract(rstack, locations)
    values[[id_column]] <- locations[[id_column]]
    names(values)[1] <- "peiid"
    hz <- parse_layers(values, idcol = "peiid")
    depths(hz) <- peiid ~ hzdept + hzdepb
    site(hz) <- data.frame(peiid = values$peiid)
    return(hz)
  }

  if (mode == "zonal") {
    stopifnot(!is.null(locations))
    # Subset only polygons with an ID present
    if (!inherits(locations, "SpatVector")) {
      locations <- terra::vect(locations)
    }
    if (!id_column %in% names(locations)) {
      stop("ID column not found in 'locations'.")
    }
    locations <- locations[!is.na(locations[[id_column]]), ]
    stopifnot(!is.null(locations))
    # Already a SpatVector at this point
    crs_rstack <- crs(rstack)
    if (!terra::same.crs(rstack, locations)) {
      locations <- terra::project(locations, crs(rstack))
    }

    stat_list <- lapply(stat, function(s) {
      terra::zonal(rstack, locations, fun = s)
    })
    df <- Reduce(function(x, y) full_join(x, y, by = id_column), stat_list)
    df$peiid <- locations[[id_column]]
    hz <- parse_layers(df, idcol = "peiid")
    depths(hz) <- peiid ~ hzdept + hzdepb
    site(hz) <- locations |>
      as.data.frame() |>
      select(peiid = all_of(id_column))
    return(hz)
  }

  stop("Unsupported mode. Choose from 'raster', 'points', or 'zonal'")
}
