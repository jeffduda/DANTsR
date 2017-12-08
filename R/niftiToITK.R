#' @title niftiToITK
#' @description convert world coordinates from NIFTI to ITK
#' @param niftiImage an antsImage
#' @param points the locations of interest
#' @param type 'point' or 'index'
#' @param interpolation options are: 'linear'
#' @export

library(RNifti)

niftiToItk = function(points, niftiImage=NA, antsImage=NA, filename=NA)
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

  if ( dim(points)[2] != antsImage@dimension ) {
    print(dim(points))
    print(antsImage)
    print(niftiImage)
    stop("Incompatible dimensions")

  }

  return(antsTransformIndexToPhysicalPoint(antsImage, worldToVoxel(points, niftiImage)))

}

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

  return(voxelToWorld(antsTransformIndexToPhysicalPoint(antsImage, points), niftiImage))

}
