#' Convert Extracted Point Soil Data to SoilProfileCollection
#'
#' This function extracts raster values at point locations, reshapes the data into horizon format,
#' and converts it into an `aqp::SoilProfileCollection`.
#'
#' @param rstack A `SpatRaster` object with depth-structured layers.
#' @param locations An `sf` or `SpatVector` object of point locations.
#' @param id_column A column name in `locations` that uniquely identifies each point (e.g., "Name").
#'
#' @return A `SoilProfileCollection` object with extracted and structured soil horizon data.
#' @export
pts_to_SPC <- function(rstack, locations, id_column = "Name") {
  require(terra)
  require(dplyr)
  require(tidyr)
  require(aqp)
  require(stringr)

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
        str_remove_all(var, "^_|_$")
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

  spc <- long_df
  depths(spc) <- peiid ~ hzdept + hzdepb
  site(spc) <- tibble(peiid = unique(long_df$peiid)) # Add minimal site data

  return(spc)
}
