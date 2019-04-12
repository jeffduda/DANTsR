#' @title cellCountImage
#' @description return image that counts the number of unique cells that contain a point in each voxel
#' @param mesh 'antsrMesh' of polylines
#' @param reference antsImage that defines the physical space to consider (points outside of image space are ignored)
#' @export
cellCountImage = function(mesh, reference)
{

  if ( !class(mesh)=="antsrMesh") {
    stop("mesh must be an 'antsrMesh'")
  }
  if ( !class(reference)=="antsImage") {
    stop("reference must be an 'antsImage'")
  }

  return( .Call( "cellCountImage", mesh, reference, PACKAGE="DANTsR") )
}
