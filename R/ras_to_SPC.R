#' Convert Raster Stack to SoilProfileCollection
#'
#' Converts a horizon-layered `SpatRaster` into a `SoilProfileCollection` (SPC).
#' Assumes layer names follow the format: `property_TOP_BOTTOM` (e.g., `claytotal_0_5`)
#'
#' @param rstack A `SpatRaster` with depth-layered soil properties.
#'
#' @return A `SoilProfileCollection` object with spatial site metadata.
#' @export
ras_to_SPC <- function(rstack) {
  require(terra)
  require(dplyr)
  require(tidyr)
  require(aqp)
  require(stringr)

  depth_interval_lookup <- list(
    "0_5" = c("0_5", "0-5cm", "0_cm", "0-5"),
    "5_15" = c("5_15", "5-15cm", "5_cm", "0-25"),
    "15_30" = c("15_30", "15-30cm", "15_cm", "25-50"),
    "30_60" = c("30_60", "30-60cm", "30_cm", "25-50"),
    "60_100" = c("60_100", "60-100cm", "60_cm"),
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

  df <- as.data.frame(rstack, xy = TRUE, cells = TRUE, na.rm = TRUE)
  df$peiid <- paste0("cell_", df$cell)
  site_data <- df %>% select(peiid, x, y)

  long_df <- df %>%
    select(-x, -y, -cell) %>%
    pivot_longer(cols = -peiid, names_to = "layer", values_to = "value") %>%
    separate(
      layer,
      into = c("property", "top", "bottom"),
      sep = "_",
      remove = FALSE
    ) %>%
    mutate(depth_label = paste0(top, "_", bottom)) %>%
    rowwise() %>%
    mutate(
      matched_label = {
        matches <- sapply(depth_interval_lookup, function(p) {
          any(str_detect(
            depth_label,
            paste0("^(", paste(p, collapse = "|"), ")$")
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
      }
    ) %>%
    ungroup()

  hz <- long_df %>%
    select(peiid, hzdept, hzdepb, property, value) %>%
    pivot_wider(names_from = property, values_from = value)

  depths(hz) <- peiid ~ hzdept + hzdepb
  site(hz) <- site_data
  return(hz)
}
