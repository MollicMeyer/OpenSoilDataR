#' Lookup scalar conversion factor for SOLUS properties
#'
#' @param property Character. Soil property name (e.g., "soc").
#' @param depth Character. Depth string (e.g., "30_cm").
#' @param measure Character. Measure ("p", "l", "h", or "rpi").
#'
#' @return Numeric scalar factor.
#' @export
lookup_scalar <- function(property, depth, measure) {
  ...
}

lookup_scalar <- function(property, depth, measure) {
  filetype <- switch(
    measure,
    "l" = "95% low prediction interval",
    "h" = "95% high prediction interval",
    "p" = "prediction",
    "rpi" = "relative prediction interval"
  )

  scalar_table <- structure(
    list(
      property = c(
        "anylithicdpt",
        "anylithicdpt",
        "anylithicdpt",
        "anylithicdpt",
        "anylithicdpt",
        "caco3",
        "caco3",
        "caco3",
        "caco3",
        "cec7",
        "cec7",
        "cec7",
        "cec7",
        "claytotal",
        "claytotal",
        "claytotal",
        "claytotal",
        "dbovendry",
        "dbovendry",
        "dbovendry",
        "dbovendry",
        "ec",
        "ec",
        "ec",
        "ec",
        "ecec",
        "ecec",
        "ecec",
        "ecec",
        "fragvol",
        "fragvol",
        "fragvol",
        "fragvol",
        "gypsum",
        "gypsum",
        "gypsum",
        "gypsum",
        "ph1to1h2o",
        "ph1to1h2o",
        "ph1to1h2o",
        "ph1to1h2o",
        "resdept",
        "resdept",
        "resdept",
        "resdept",
        "sandco",
        "sandco",
        "sandco",
        "sandco",
        "sandfine",
        "sandfine",
        "sandfine",
        "sandfine",
        "sandmed",
        "sandmed",
        "sandmed",
        "sandmed",
        "sandtotal",
        "sandtotal",
        "sandtotal",
        "sandtotal",
        "sandvc",
        "sandvc",
        "sandvc",
        "sandvc",
        "sandvf",
        "sandvf",
        "sandvf",
        "sandvf",
        "sar",
        "sar",
        "sar",
        "sar",
        "silttotal",
        "silttotal",
        "silttotal",
        "silttotal",
        "soc",
        "soc",
        "soc",
        "soc"
      ),
      depth = rep(c("0_cm", "5_cm", "15_cm", "30_cm"), each = 20),
      filetype = rep(c("prediction"), 80),
      scalar = c(
        1,
        1,
        1,
        1,
        0.001, # anylithicdpt
        1,
        1,
        1,
        0.01, # caco3
        0.1,
        0.1,
        0.1,
        0.01, # cec7
        1,
        1,
        1,
        0.01, # claytotal
        0.01,
        0.01,
        0.01,
        0.01, # dbovendry
        0.1,
        0.1,
        0.1,
        0.01, # ec
        0.1,
        0.1,
        0.1,
        0.01, # ecec
        1,
        1,
        1,
        0.01, # fragvol
        0.1,
        0.1,
        0.1,
        0.01, # gypsum
        0.01,
        0.01,
        0.01,
        0.01, # ph1to1h2o
        1,
        1,
        1,
        0.01, # resdept
        1,
        1,
        1,
        0.01, # sandco
        1,
        1,
        1,
        0.01, # sandfine
        1,
        1,
        1,
        0.01, # sandmed
        1,
        1,
        1,
        0.01, # sandtotal
        1,
        1,
        1,
        0.01, # sandvc
        1,
        1,
        1,
        0.01, # sandvf
        1,
        1,
        1,
        0.01, # sar
        1,
        1,
        1,
        0.1, # silttotal
        0.001,
        0.001,
        0.001,
        0.01
      ) # soc
    ),
    class = "data.frame"
  )

  match_row <- subset(
    scalar_table,
    property == property &
      depth == depth &
      filetype == filetype
  )

  if (nrow(match_row) == 1) {
    return(as.numeric(match_row$scalar))
  } else {
    warning(
      "No matching scalar found for ",
      property,
      ", ",
      depth,
      ", ",
      measure
    )
    return(1)
  }
}
