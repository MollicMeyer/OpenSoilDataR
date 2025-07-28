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

  # Raster to dataframe
  df <- as.data.frame(rstack, xy = TRUE, cells = TRUE, na.rm = TRUE)
  df$peiid <- paste0(source, "_cell_", df$cell)
  site_data <- df %>% select(peiid, x, y)

  long_df <- df %>%
    select(-x, -y, -cell) %>%
    pivot_longer(cols = -peiid, names_to = "layer", values_to = "value") %>%
    rowwise() %>%
    mutate(
      # Match the key in depth_interval_lookup whose values match the layer name
      matched_label = {
        match_found <- NA_character_
        for (k in names(depth_interval_lookup)) {
          pats <- depth_interval_lookup[[k]]
          if (
            any(str_detect(
              layer,
              paste0("(", paste(pats, collapse = "|"), ")")
            ))
          ) {
            match_found <- k
            break
          }
        }
        match_found
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

      # Remove matched depth label + suffix from variable name
      variable = if (!is.na(matched_label)) {
        str_remove(
          layer,
          paste0(
            "(_)?(",
            paste(depth_interval_lookup[[matched_label]], collapse = "|"),
            ").*$"
          )
        )
      } else {
        NA_character_
      }
    ) %>%
    ungroup() %>%
    mutate(value = as.numeric(value)) %>%
    filter(!is.na(hzdept), !is.na(hzdepb)) %>%
    # NO summarise — keep all depths for each variable
    pivot_wider(names_from = variable, values_from = value)

  # Construct SoilProfileCollection
  depths(long_df) <- peiid ~ hzdept + hzdepb
  site(long_df) <- site_data

  return(long_df)
}
