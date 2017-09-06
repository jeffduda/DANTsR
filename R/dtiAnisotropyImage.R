#' @title dtiAnisotropyImage
#' @description calculate an anistropy image from DTI
#' @param img a 6-channel antsImage
#' @param method which type of anisotropy index to calculate
#' @export
dtiAnisotropyImage = function(img, method="FractionalAnisotropy")
{
  if ( !class(img)=="antsImage") {
    stop("Must input an 'antsImage'")
  }

  if ( img@components != 6 ) {
    stop("Input image must have 6 components")
  }

  method = tolower(method)

  .Call( "dtiFilters", img, method, package="DANTsR")

}
