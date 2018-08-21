#' @title dtiAnisotropy
#' @description calculate an anistropy image from DTI
#' @param img a 6-channel antsImage
#' @param method the measure to return
#' \itemize{
#'   \item{FractionalAnisotropy}{}
#'   \item{RelativeAnisotropy}{}
#' }
#' @export
dtiAnisotropy= function(img, method="FractionalAnisotropy")
{
  if ( !class(img)=="antsImage") {
    stop("Must input an 'antsImage'")
  }

  if ( img@components != 6 ) {
    stop("Input image must have 6 components")
  }

  method = tolower(method)

  .Call( "dtiFilters", img, method, PACKAGE="DANTsR")

}
