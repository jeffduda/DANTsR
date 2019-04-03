#' @title physicalVectorsToIndexVectors
#' @description convert image of physical vectors to index vectors
#' @param image antsImage of pixels to transform
#' @param in.place modify the input image instead of returning a new image
#' @return antsImage
#' @examples
#' @export
physicalVectorsToIndexVectors <- function(image, in.place=FALSE) {

  check_ants(image)

  if ( image@dimension != image@components ) {
    stop("Image dimension must equal number of components")
  }

  returnImg = image
  if ( !in.place ) {
    returnImg = antsImageClone(image)
  }

  .Call("physicalVectorsToIndexVectors", returnImage, TRUE, package="DANTsR")
  return(returnImage)

}

#' @title indexVectorsToPhysicalVectors
#' @description convert image of index vectors to physical vectors
#' @param image antsImage of pixels to transform
#' @param in.place modify the input image instead of returning a new image
#' @return antsImage
#' @examples
#' @export
indexVectorsToPhysicalVectors <- function(image, in.place=FALSE) {

  check_ants(image)

  if ( image@dimension != image@components ) {
    stop("Image dimension must equal number of components")
  }

  returnImg = image
  if ( !in.place ) {
    returnImg = antsImageClone(image)
  }

  .Call("physicalVectorsToIndexVectors", returnImage, FALSE, package="DANTsR")
  return(returnImage)

}
