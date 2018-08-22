#' @title labelsToPoints
#' @description return point coordinates of all voxels with specified label
#' @param img an antsImage
#' @param label the value of interest
#' @export
labelsToPoints = function(img, label=1)
{
  if ( !class(img)=="antsImage") {
    stop("Must input an 'antsImage'")
  }

  .Call( "labelsToPoints", img, label, package="DANTsR")

}
