#' Convert RasterStack to SoilProfileCollection
#'
#' @param rstack A `SpatRaster` with layer names indicating depth intervals (e.g., "claytotal_r_0_5")
#' @return A `SoilProfileCollection` object
#' @export
ras_to_SPC <- function(rstack, source = "R") {
  require(terra)
  require(dplyr)
  require(tidyr)
  require(aqp)
  require(stringr)

  # Lookup tables
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

  # Convert raster to data.frame
  df <- as.data.frame(rstack, xy = TRUE, cells = TRUE, na.rm = TRUE)
  df$peiid <- paste0(source, "_cell_", df$cell)
  site_data <- df %>% select(peiid, x, y)

  # Long format with depth & variable parsing
  long_df <- df %>%
    select(-x, -y, -cell) %>%
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
        var <- if (
          str_detect(layer, paste0("_?", fixed(matched_string), "_?$"))
        ) {
          str_replace(layer, paste0("_?", fixed(matched_string), "_?$"), "")
        } else {
          str_replace(layer, fixed(matched_string), "")
        }
        var <- str_replace_all(var, "__+", "_") # Replace double underscore
        var <- str_remove_all(var, "^_|_$") # Trim leading/trailing _
        var
      } else {
        NA_character_
      }
    ) %>%
    ungroup() %>%
    filter(!is.na(hzdept), !is.na(hzdepb), !is.na(variable)) %>%
    mutate(value = as.numeric(value)) %>%
    pivot_wider(names_from = variable, values_from = value)

  # Construct SoilProfileCollection
  depths(long_df) <- peiid ~ hzdept + hzdepb
  site(long_df) <- site_data

  return(long_df)
}
