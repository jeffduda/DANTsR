% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/antsrMesh_class.R
\name{write.antsrMesh}
\alias{write.antsrMesh}
\title{write.antsrMesh}
\usage{
write.antsrMesh(mesh, filename, image = NULL, seeds = NULL,
  cells.as = "polygon", binary = TRUE)
}
\arguments{
\item{mesh}{antsrMesh to write}

\item{filename}{name of the file to read}

\item{image}{reference image (for TrackVis files)}

\item{seeds}{tract seed indices (for Camino files)}

\item{cells.as}{what VTK-type of cell should cell data be written as ("NA", "polygon", "line")}

\item{binary}{write data as binary (default is TRUE)}
}
\description{
write an antsrMesh to a file
}
