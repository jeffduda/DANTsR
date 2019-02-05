#' @title cellPointInMask
#' @description indicate if a cell has a point inside an image mask
#' @param mesh 'antsrMesh defining a point set and cells
#' @param mask an 'antsImage'
#' @param index the cell/s to examine

cellPointInMask = function( mesh, mask, index=NULL ) {

  if ( is.null(index) ) {
    index = c(1:antsrMeshGetNumberOfCells(mesh))
  }

  iMask = mask
  if ( mask@precision != mesh@precision ) {
    iMask = antsImageTypeCast(mask, mesh@precision)
  }

  return( .Call( "cellPointInMask", mesh, iMask, index, PACKAGE="DANTsR") )

}
