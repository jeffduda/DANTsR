#' @title dtiExp
#' @description calculate exp euclidean tensor
#' @param img a 6-channel antsImage
#' @export
dtiExp = function(img)
{
  if ( !class(img)=="antsImage") {
    stop("Must input an 'antsImage'")
  }

  if ( img@components != 6 ) {
    stop("Input image must have 6 components")
  }

  return(.Call( "dtiFilters", img, "exp", PACKAGE="DANTsR"))

}
