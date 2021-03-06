% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interpolateImageValues.R
\name{interpolateImageValues}
\alias{interpolateImageValues}
\title{interpolateImageValues}
\usage{
interpolateImageValues(img, points, type = "point",
  interpolation = "linear")
}
\arguments{
\item{img}{an antsImage}

\item{points}{the locations of interest}

\item{type}{'point' or 'index'}

\item{interpolation}{options are: 'linear'}
}
\description{
return image values at points or indices
}
