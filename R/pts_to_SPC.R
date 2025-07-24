#' Convert Extracted Point Soil Data to SoilProfileCollection
#'
#' @param rstack A `SpatRaster` object
#' @param locations An `sf` or `SpatVector` object of point locations
#' @param id_column A column name in `locations` that uniquely identifies each point (e.g., "Name")
#' @return A `SoilProfileCollection` object
#' @export
pts_to_SPC <- function(rstack, locations, id_column = "Name") {
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

  if (!inherits(locations, "SpatVector")) {
    locations <- terra::vect(locations)
  }

  extracted <- terra::extract(rstack, locations)
  extracted[[id_column]] <- locations[[id_column]]
  names(extracted)[1] <- "peiid"

  long_df <- extracted %>%
    pivot_longer(
      cols = -c(peiid, all_of(id_column)),
      names_to = "layer",
      values_to = "value"
    ) %>%
    rowwise() %>%
    mutate(
      matched_label = {
        matches <- sapply(depth_interval_lookup, function(pats) {
          any(str_detect(
            layer,
            paste0("_(", paste(pats, collapse = "|"), ")$")
          ))
        })
        if (any(matches)) {
          names(depth_interval_lookup)[which(matches)[1]]
        } else {
          NA_character_
        }
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
      variable = if (!is.na(matched_label)) {
        str_remove(
          layer,
          paste0(
            "_(",
            paste(depth_interval_lookup[[matched_label]], collapse = "|"),
            ")$"
          )
        )
      } else {
        NA_character_
      }
    ) %>%
    ungroup()

  hz_data <- long_df %>%
    select(peiid, hzdept, hzdepb, variable, value) %>%
    pivot_wider(names_from = variable, values_from = value)

  depths(hz_data) <- peiid ~ hzdept + hzdepb
  site(hz_data) <- extracted %>% select(peiid)

  return(hz_data)
}
