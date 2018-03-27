#' @title niftiToItk
#' @description convert world coordinates from NIFTI to ITK
#' @param niftiImage an niftiImage
#' @param antsImage an antsImage
#' @param filename name of image file
#' @param points the locations of interest (matrix or antsrMesh)
#' @param in.place transform points "in-place", for antsrMesh only
#' @export

#library(RNifti)
niftiToItk = function(points, niftiImage=NA, antsImage=NA, filename=NA, in.place=FALSE)
{

  if ( !is.na(filename) ) {
    if ( is.na(niftiImage) ) {
      niftiImage = readNifti(filename)
    }
    if ( is.na(antsImage) ) {
      antsImage = antsImageRead(filename)
    }
  }

  if ( class(niftiImage)[1] != "niftiImage") {
    stop("Invalid niftiImage")
  }
  if ( class(antsImage)[1] != "antsImage") {
    stop("Invalid antsImage")
  }



  if ( class(points)=="matrix" ) {

    if ( dim(points)[2] != antsImage@dimension ) {
      stop("Incompatible dimensions")
    }

    return(antsTransformIndexToPhysicalPoint(antsImage, RNifti::worldToVoxel(points, niftiImage)))
  }
  else if ( class(points)=="antsrMesh" ) {

    if ( points@dimension != antsImage@dimension ) {
      stop("Incompatible dimensions")
    }

    ptMatrix = antsTransformIndexToPhysicalPoint(antsImage, RNifti::worldToVoxel(antsrMeshGetPoints(points), niftiImage))

    if ( in.place ) {
      for ( i in 1:antsrMeshGetNumberOfPoints(points) ) {
        antsrMeshSetPoint(points, ptMatrix[i,], i)
      }
      invisible(NA)
    }
    else {
      return( antsrMeshCreate(dimension=dim(ptMatrix)[2], points=ptMatrix))
    }

  }
  else {
    stop("Invalid points argument, must be 'matrix' or 'antsrMesh'")
  }


}

#' @title itkToNifti
#' @description convert world coordinates from ITK to NIFTI
#' @param niftiImage an niftiImage
#' @param antsImage an antsImage
#' @param filename name of image file
#' @param points the locations of interest
#' @export
itkToNifti = function(points, niftiImage=NA, antsImage=NA, filename=NA)
{

  if ( !is.na(filename) ) {
    if ( is.na(niftiImage) ) {
      niftiImage = readNifti(filename)
    }
    if ( is.na(antsImage) ) {
      antsImage = antsImageRead(filename)
    }
  }

  if ( class(niftiImage)[1] != "niftiImage") {
    stop("Invalid niftiImage")
  }
  if ( class(antsImage)[1] != "antsImage") {
    stop("Invalid antsImage")
  }

  if ( dim(points)[2] == antsImage@dimension ) {
    stop("Incompatible dimensions")
  }

  return(RNifti::voxelToWorld(antsTransformIndexToPhysicalPoint(antsImage, points), niftiImage))

}
