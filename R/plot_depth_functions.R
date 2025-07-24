#' Plot Depth Functions from Multiple SoilProfileCollections
#'
#' @param spc_list A list of SoilProfileCollection objects to compare.
#' @param source_labels A character vector of labels (e.g., "PSP", "SGO") matching `spc_list`.
#' @param variables Character vector of standardized property names to plot (e.g., "sand", "clay").
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

  # Property Lookup Table for use in s.zonalstats()
  # Property Lookup Table for use in s.zonalstats()
  property_lookup <- list(
    # Water Retention
    "wtenthbar" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "wtenthbar_r"
    ),
    "wthirdbar" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "wthirdbar_r"
    ),
    "wfifteenbar" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "wfifteenbar_r"
    ),
    "awc" = c(PSP = NA, SG2 = NA, SOL = "awc_r", CSRL = NA, SGO = "awc_r"),

    # Phosphorus
    "pbray1" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "pbray1_r"),
    "ptotal" = c(
      PSP = NA,
      SG2 = NA,
      SOL = "ptotal_r",
      CSRL = NA,
      SGO = "ptotal_r"
    ),
    "ph2osoluble" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "ph2osoluble_r"
    ),
    "poxalate" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "poxalate_r"),

    # Bulk Density
    "bd" = c(
      PSP = "bd_mean",
      SG2 = "bdod",
      SOL = "dbovendry_r",
      CSRL = "bulk_density",
      SGO = "dbovendry_r"
    ),
    "bd_tenthbar" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "dbtenthbar_r"
    ),
    "bd_thirdbar" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "dbthirdbar_r"
    ),
    "bd_fifteenbar" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "dbfifteenbar_r"
    ),

    # Texture
    "clay" = c(
      PSP = "clay_mean",
      SG2 = "clay",
      SOL = "claytotal",
      CSRL = "clay_profile",
      SGO = "claytotal_r"
    ),
    "sand" = c(
      PSP = "sand_mean",
      SG2 = "sand",
      SOL = "sandtotal",
      CSRL = "sand_profile",
      SGO = "sandtotal_r"
    ),
    "silt" = c(
      PSP = "silt_mean",
      SG2 = "silt",
      SOL = "silttotal",
      CSRL = "silt_profile",
      SGO = "silttotal_r"
    ),
    "sandco" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "sandco_r"),
    "sandfine" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "sandfine_r"),
    "sandmed" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "sandmed_r"),
    "sandvc" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "sandvc_r"),
    "sandvf" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "sandvf_r"),
    "siltco" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "siltco_r"),
    "siltfine" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "siltfine_r"),

    # Organic Carbon
    "som" = c(
      PSP = "om_mean",
      SG2 = NA,
      SOL = NA,
      CSRL = "som_max",
      SGO = "om_r"
    ),
    "soc" = c(
      PSP = "soc",
      SG2 = "soc",
      SOL = "soc",
      CSRL = "soc_max",
      SGO = "oc_r"
    ),
    "soc_stock" = c(
      PSP = NA,
      SG2 = "ocs",
      SOL = NA,
      CSRL = "soc_stock",
      SGO = NA
    ),

    # Chemistry
    "ph" = c(
      PSP = "ph_mean",
      SG2 = "phh2o",
      SOL = "ph1to1h2o",
      CSRL = "ph_profile",
      SGO = "ph1to1h2o_r"
    ),
    "ph_cacl2" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "ph01mcacl2_r"
    ),
    "ph_oxidized" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "phoxidized_r"
    ),
    "cec" = c(
      PSP = NA,
      SG2 = "cec",
      SOL = "cec7",
      CSRL = "cec_profile",
      SGO = "cec7_r"
    ),
    "ecec" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "ecec_r"),
    "ec" = c(PSP = NA, SG2 = NA, SOL = "ec", CSRL = "ec_profile", SGO = "ec_r"),
    "ec15" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "ec15_r"),
    "esp" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "esp_r"),
    "sar" = c(PSP = NA, SG2 = NA, SOL = "sar", CSRL = "sar", SGO = "sar_r"),
    "gypsum" = c(
      PSP = NA,
      SG2 = NA,
      SOL = "gypsum",
      CSRL = NA,
      SGO = "gypsum_r"
    ),
    "caco3" = c(
      PSP = NA,
      SG2 = NA,
      SOL = "caco3",
      CSRL = "caco3",
      SGO = "caco3_r"
    ),

    # Extractables
    "extracid" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "extracid_r"),
    "extral" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "extral_r"),
    "aloxalate" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "aloxalate_r"
    ),
    "feoxalate" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "feoxalate_r"
    ),

    # Mechanics
    "lep" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "lep_r"),
    "ll" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "ll_r"),
    "pi" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "pi_r"),

    # Rock Fragments
    "frag3to10" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "frag3to10_r"
    ),
    "fraggt10" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "fraggt10_r"),

    # Fiber
    "fiberrubbedpct" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "fiberrubbedpct_r"
    ),
    "fiberunrubbedpct" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "fiberunrubbedpct_r"
    ),

    # Sieve
    "sieveno10" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "sieveno10_r"
    ),
    "sieveno200" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "sieveno200_r"
    ),
    "sieveno4" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "sieveno4_r"),
    "sieveno40" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "sieveno40_r"
    ),

    # Hydrology
    "ksat" = c(
      PSP = "ksat_mean",
      SG2 = NA,
      SOL = NA,
      CSRL = "ksat_max",
      SGO = "ksat_r"
    ),
    "wsatiated" = c(
      PSP = NA,
      SG2 = NA,
      SOL = NA,
      CSRL = NA,
      SGO = "wsatiated_r"
    ),

    # Erosion
    "kw" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "kwfact"),
    "kffact" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "kffact"),
    "kifact" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "kifact"),
    "krfact" = c(PSP = NA, SG2 = NA, SOL = NA, CSRL = NA, SGO = "krfact")
  )

  # Rename variables in horizons
  for (i in seq_along(spc_list)) {
    src <- source_labels[i]
    spc <- spc_list[[i]]

    for (v in variables) {
      actual_name <- property_lookup[[v]][[src]]
      if (!is.na(actual_name)) {
        if (actual_name %in% horizonNames(spc)) {
          horizons(spc)[[v]] <- horizons(spc)[[actual_name]]
        } else {
          stop(sprintf(
            "Variable '%s' (mapped as '%s') not found in horizons of SPC for source '%s'",
            v,
            actual_name,
            src
          ))
        }
      } else {
        stop(sprintf(
          "Variable '%s' is not defined for source '%s' in the lookup table",
          v,
          src
        ))
      }
    }

    profile_id(spc) <- paste0(src, "_", profile_id(spc))
    site(spc)$source <- src
    spc_list[[i]] <- spc
  }

  # Combine SPCs
  combined_spc <- do.call(aqp::combine, spc_list)

  if ("source" %in% names(horizons(combined_spc))) {
    horizons(combined_spc)$source <- NULL
  }

  mean.and.sd <- function(values) {
    m <- mean(values, na.rm = TRUE)
    s <- sd(values, na.rm = TRUE)
    c(mean = m, lower = m - s, upper = m + s)
  }

  slab_list <- lapply(variables, function(var) {
    slab(
      combined_spc,
      fm = stats::as.formula(paste(var, "~ source")),
      slab.structure = slab_structure,
      slab.fun = mean.and.sd
    )
  })
  slab_df <- bind_rows(slab_list)

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
  palette_colors <- okabe_ito[seq_len(length(unique(slab_df$source)))]

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
      superpose.line = list(col = palette_colors, lwd = 2)
    ),
    panel = panel.depth_function,
    prepanel = prepanel.depth_function,
    strip = strip.custom(bg = grey(0.8)),
    auto.key = list(columns = 2, lines = TRUE, points = FALSE)
  )
}
