#' @title dtiLog
#' @description calculate log euclidean tensor
#' @param img a 6-channel antsImage
#' @export
dtiLog= function(img)
{
  if ( !class(img)=="antsImage") {
    stop("Must input an 'antsImage'")
  }

  if ( img@components != 6 ) {
    stop("Input image must have 6 components")
  }

  return(.Call( "dtiFilters", img, "log", PACKAGE="DANTsR"))

}
