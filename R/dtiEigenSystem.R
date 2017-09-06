#' @title dtiEigenSystem
#' @description calculate eigen values and vectors
#' @param img a 6-channel antsImage
#' @export
dtiEigenSystem = function(img)
{
  if ( !class(img)=="antsImage") {
    stop("Must input an 'antsImage'")
  }

  if ( img@components != 6 ) {
    stop("Input image must have 6 components")
  }

  .Call( "dtiFilters", img, "eigensystem", package="DANTsR")

}
