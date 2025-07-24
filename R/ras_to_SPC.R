#' Convert Raster Stack to SoilProfileCollection
#'
#' Converts a raster stack from a fetch_*() function into a SoilProfileCollection.
#'
#' @param rstack A `SpatRaster` of soil property layers, e.g., from `fetch_SGO()`.
#'
#' @return A SoilProfileCollection with site info from raster cells.
#' @export
ras_to_SPC <- function(rstack) {
  require(terra); require(dplyr); require(tidyr); require(aqp); require(stringr)

  df <- as.data.frame(rstack, xy = TRUE, cells = TRUE, na.rm = TRUE)
  df$peiid <- paste0("cell_", df$cell)
  site_data <- df %>% select(peiid, x, y)

  long_df <- df %>%
    select(-x, -y, -cell) %>%
    pivot_longer(cols = -peiid, names_to = "layer", values_to = "value") %>%
    separate(layer, into = c("property", "top", "bottom"), sep = "_", remove = FALSE) %>%
    mutate(
      depth_label = paste0(top, "_", bottom),
      matched_label = names(depth_interval_lookup)[
        sapply(depth_interval_lookup, function(p) any(str_detect(depth_label, paste0("^(", paste(p, collapse = "|"), ")$")))
      ],
      hzdept = depth_range_lookup[[matched_label]][1],
      hzdepb = depth_range_lookup[[matched_label]][2]
    ) %>%
    mutate(variable = property) %>%
    ungroup()

  hz <- long_df %>%
    select(peiid, hzdept, hzdepb, variable, value) %>%
    pivot_wider(names_from = variable, values_from = value)

  depths(hz) <- peiid ~ hzdept + hzdepb
  site(hz) <- site_data
  return(hz)
}
