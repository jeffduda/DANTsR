#' @title physicalTensorsToIndexTensors
#' @description convert image of physical Tensors to index Tensors
#' @param image antsImage of pixels to transform
#' @param in.place modify the input image instead of returning a new image
#' @return antsImage
#' @examples
#' @export
physicalTensorsToIndexTensors <- function(image, in.place=FALSE) {

  check_ants(image)

  if ( image@components != 6 ) {
    stop("Image components must be 6")
  }

  returnImg = image
  if ( !in.place ) {
    returnImg = antsImageClone(image)
  }

  dir = antsGetDirection(image)
  tx = createAntsrTransform(dim=3, type="AffineTransform", precision="double", matrix=dir)

  .Call("antsrTransform_TransformPixels", tx, returnImg, "diffusiontensor", PACKAGE="DANTsR")

  if ( in.place ) {
    return(NA)
  }
  return(returnImg)

}

#' @title indexTensorsToPhysicalTensors
#' @description convert image of index Tensors to physical Tensors
#' @param image antsImage of pixels to transform
#' @param in.place modify the input image instead of returning a new image
#' @return antsImage
#' @examples
#' @export
indexTensorsToPhysicalTensors <- function(image, in.place=FALSE) {

  check_ants(image)

  if ( image@components != 6 ) {
    stop("Image components must be 6")
  }

  returnImg = image
  if ( !in.place ) {
    returnImg = antsImageClone(image)
  }

  dir = antsGetDirection(image)
  tx = createAntsrTransform(dim=3, type="AffineTransform", precision="double", matrix=dir)
  tx = invertAntsrTransform(tx)

  .Call("antsrTransform_TransformPixels", tx, returnImg, "diffusiontensor", PACKAGE="DANTsR")

  if ( in.place ) {
    return(NA)
  }
  return(returnImg)

}
