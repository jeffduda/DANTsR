#' @title streamlineTargetsFromSeed
#' @description from seed, find labeled targets
#' @param mesh 'antsrMesh defining a point set and cells (streamlines)
#' @param seeds a vector indicating the seed index for each streamline
#' @param image an 'antsImage' of labeled regions
#' @param index the cell/s to examine

streamlineTargetsFromSeed = function( mesh, seeds, image, index ) {

  if ( is.null(index) ) {
    index = c(1:antsrMeshGetNumberOfCells(mesh))
  }

  iRegions = image
  if ( image@precision != mesh@precision ) {
    iRegions = antsImageTypeCast(mask, mesh@precision)
  }

  return( .Call( "streamlineTargetsFromSeed", mesh, seeds, iRegions, index, PACKAGE="DANTsR") )

}
