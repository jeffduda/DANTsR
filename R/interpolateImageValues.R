#' @title interpolateImageValues
#' @description return image values at points or indices
#' @param img an antsImage
#' @param points the locations of interest
#' @param type 'point' or 'index'
#' @param interpolation options are: 'linear'
#' @export
interpolateImageValues = function(img, points, type="point", interpolation="linear")
{
  if ( !class(img)=="antsImage") {
    stop("Must input an 'antsImage'")
  }

  if ( !class(points)=="matrix") {
    stop("points must be a 'matrix'")
  }

  .Call( "interpolateImageValues", img, points, type, interpolation, PACKAGE="DANTsR")

}
