#' @title polylineCountImage
#' @description return point coordinates of all voxels with specified label
#' @param mesh 'antsrMesh' of polylines
#' @param reference antsImage that defines the physical space to consider
#' @export

polylineCountImage = function(mesh, reference)
{

  if ( !class(mesh)=="antsrMesh") {
    stop("mesh must be an 'antsrMesh'")
  }
  if ( !class(reference)=="antsImage") {
    stop("reference must be an 'antsImage'")
  }

  return(.Call("polylineCountImage", mesh, reference, package="DANTsR"))

}
