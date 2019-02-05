#' @title cellImageValueSummary
#' @description summary of image value sover all points in a cell (mean, max, etc)
#' @param mesh 'antsrMesh defining a point set and cells
#' @param image an 'antsImage' to average over
#' @param measure summary measure to obtain (default="mean")
#' @param index the cell/s to examine

cellImageValueSummary = function( mesh, image, measure="mean", index=NULL ) {

  if ( is.null(index) ) {
    index = c(1:antsrMeshGetNumberOfCells(mesh))
  }

  vImage = image
  if ( image@precision != mesh@precision ) {
    vImage = antsImageTypeCast(image, mesh@precision)
  }

  return( .Call( "cellImageValueSummary", mesh, vImage, measure, index, PACKAGE="DANTsR") )

}
