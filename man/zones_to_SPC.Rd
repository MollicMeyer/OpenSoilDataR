% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zones_to_SPC.R
\name{zones_to_SPC}
\alias{zones_to_SPC}
\title{Convert zonal raster statistics into a SoilProfileCollection}
\usage{
zones_to_SPC(rstack, zones, stat = "mean", id_column = "Name")
}
\arguments{
\item{rstack}{A \code{SpatRaster} containing stacked soil property rasters.}

\item{zones}{A \code{SpatVector} or \code{sf} object of polygons used for zonal summaries.}

\item{stat}{A summary statistic to apply per zone (e.g., "mean", "median"). Default: "mean".}

\item{id_column}{Character. Name of the column in \code{zones} to use as profile IDs.}
}
\value{
A \code{SoilProfileCollection} with horizon data extracted from the raster stack.
}
\description{
This function extracts zonal statistics from a raster stack using polygons (zones),
reshapes the data into horizon format, and converts it into an \code{aqp::SoilProfileCollection}.
}
