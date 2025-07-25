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

  # Extract raster to data.frame
  df <- as.data.frame(rstack, xy = TRUE, cells = TRUE, na.rm = TRUE)
  df$peiid <- paste0(source, "_cell_", df$cell)
  site_data <- df %>% select(peiid, x, y)

  long_df <- df %>%
    select(-x, -y, -cell) %>%
    pivot_longer(cols = -peiid, names_to = "layer", values_to = "value") %>%
    rowwise() %>%
    mutate(
      matched_label = {
        matches <- sapply(depth_interval_lookup, function(pats) {
          any(str_detect(layer, paste0("(", paste(pats, collapse = "|"), ")")))
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

      variable = case_when(
        str_detect(
          layer,
          paste0(
            "_\\d+_cm_(",
            paste(c("p", "r", "l", "h", "rpi", "mean"), collapse = "|"),
            ")$"
          )
        ) ~
          str_remove(
            layer,
            paste0(
              "_\\d+_cm_(",
              paste(c("p", "r", "l", "h", "rpi", "mean"), collapse = "|"),
              ")$"
            )
          ),
        str_detect(
          layer,
          paste0(
            "_(\\d+(-|_)\\d+|\\d+)(cm)?(_(",
            paste(c("p", "r", "l", "h", "rpi", "mean"), collapse = "|"),
            "))?$"
          )
        ) ~
          str_remove(
            layer,
            paste0(
              "_(\\d+(-|_)\\d+|\\d+)(cm)?(_(",
              paste(c("p", "r", "l", "h", "rpi", "mean"), collapse = "|"),
              "))?$"
            )
          ),
        TRUE ~ layer
      )
    ) %>%
    ungroup()

  # Pivot into horizon format
  hz_data <- long_df %>%
    select(peiid, hzdept, hzdepb, variable, value) %>%
    pivot_wider(names_from = variable, values_from = value)

  # Construct SPC
  depths(hz_data) <- peiid ~ hzdept + hzdepb
  site(hz_data) <- site_data

  return(hz_data)
}
