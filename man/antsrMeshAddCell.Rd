% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/antsrMesh_class.R
\name{antsrMeshAddCell}
\alias{antsrMeshAddCell}
\title{antsrMeshAddCell}
\usage{
antsrMeshAddCell(mesh, points, type, identifier = NA)
}
\arguments{
\item{mesh}{an 'antsrMesh'}

\item{points}{array of point indices}

\item{type}{of cell to add}

\item{identifier}{index of cell to add}
}
\description{
add a cell to the mesh
}
\examples{
x =  antsrMeshCreate( 3, "float", reserve=128 )
antsrMeshAddPoint( x, c(0,0,0), 1 )
antsrMeshAddPoint( x, c(1,0,0), 2 )
antsrMeshAddPoint( x, c(1,1,0), 3 )
antsrMeshAddCell( x, c(0,1,2), "polyline", 1)
}
