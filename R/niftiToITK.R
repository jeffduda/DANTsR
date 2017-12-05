#' @title niftiToITK
#' @description return image values at points or indices
#' @param img an antsImage
#' @param points the locations of interest
#' @param type 'point' or 'index'
#' @param interpolation options are: 'linear'
#' @export
niftiToItk = function(points, direction, dimension=3)
{

  if ( ! dim(points)[2] == dimension ) {
    print(dim(points))
    stop("Incompatible dimensions")

  }

  niftimat = matrix(0,dimension,dimension)
  niftimat[1,1] = -1
  niftimat[2,2] = -1
  if ( dimension > 2 ) {
    niftimat[3,3] = 1
  }

  mat = t(niftimat %*% direction) %*% direction

  return (points %*% mat)

}
