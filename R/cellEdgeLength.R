#' @title cellEdgeLength
#' @description sum of Lengths of cell edges
#' @param mesh 'antsrMesh defining a point set and cells
#' @param index the cell to examine

cellEdgeLength = function( mesh, index=NULL ) {

  if ( is.null(index) ) {
    index = c(1:antsrMeshGetNumberOfCells(mesh))
  }

  return( .Call( "cellEdgeLength", mesh, index, PACKAGE="DANTsR") )

}
