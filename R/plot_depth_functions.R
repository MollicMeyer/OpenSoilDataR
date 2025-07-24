#' Plot Depth Functions from Multiple SoilProfileCollections
#'
#' @param spc_list A list of SoilProfileCollection objects to compare.
#' @param source_labels A character vector of labels (e.g., "PSP", "SGO") matching `spc_list`.
#' @param variables Character vector of horizon-level variables to plot.
#' @param slab_structure Numeric vector of depth breaks (e.g., c(0, 5, 15, 30, 60, 100)).
#'
#' @return A lattice depth function plot (mean ± SD) comparing groups.
#' @export
plot_depth_functions <- function(
  spc_list,
  source_labels,
  variables,
  slab_structure
) {
  stopifnot(length(spc_list) == length(source_labels))

  # Load required packages
  library(aqp)
  library(dplyr)
  library(lattice)
  library(latticeExtra)

  # Assign source label to each SPC
  for (i in seq_along(spc_list)) {
    source_name <- source_labels[i]
    profile_id(spc_list[[i]]) <- paste0(
      source_name,
      "_",
      profile_id(spc_list[[i]])
    )
    site(spc_list[[i]])$source <- source_name
  }

  # Combine SPCs using aqp::combine
  combined_spc <- do.call(aqp::combine, spc_list)

  # Remove any residual `source` from horizons
  if ("source" %in% names(horizons(combined_spc))) {
    horizons(combined_spc)$source <- NULL
  }

  # Define slab summary function
  mean.and.sd <- function(values) {
    m <- mean(values, na.rm = TRUE)
    s <- sd(values, na.rm = TRUE)
    c(mean = m, lower = m - s, upper = m + s)
  }

  # Compute slabbed summary
  slab_df <- slab(
    combined_spc,
    fm = source ~ .,
    slab.structure = slab_structure,
    slab.fun = mean.and.sd,
    variables = variables
  )

  # Okabe-Ito colorblind-friendly palette (max 8 groups)
  okabe_ito <- c(
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7",
    "#999999"
  )

  group_levels <- unique(slab_df$source)
  n_groups <- length(group_levels)

  # Set color palette
  palette_colors <- okabe_ito[seq_len(n_groups)]

  # Build lattice plot
  xyplot(
    bottom ~ mean | variable,
    data = slab_df,
    groups = source,
    lower = slab_df$lower,
    upper = slab_df$upper,
    sync.colors = TRUE,
    alpha = 0.4,
    ylab = "Depth (cm)",
    xlab = "",
    ylim = c(max(slab_structure), 0),
    layout = c(length(variables), 1),
    scales = list(
      x = list(relation = 'free', alternating = 2),
      y = list(alternating = 3)
    ),
    par.settings = list(
      superpose.line = list(
        col = palette_colors,
        lwd = 2
      )
    ),
    panel = panel.depth_function,
    prepanel = prepanel.depth_function,
    strip = strip.custom(bg = grey(0.8)),
    auto.key = list(columns = 2, lines = TRUE, points = FALSE)
  )
}
