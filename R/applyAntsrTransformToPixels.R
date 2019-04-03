#' @title applyAntsrTransformToPixels
#' @description Apply transform to pixel data (vectors, tensors, etc)
#' @param transform antsrTransform
#' @param image antsImage of pixels to transform
#' @param pixels.as type of data stored as pixels ("vector", "diffusiontensor")
#' @param in.place modify the input image instead of returning a new image
#' @return antsImage
#' @examples
#' img <- antsImageRead(getANTsRData("r16"))
#' tx = new("antsrTransform", precision="float", type="AffineTransform", dimension=2 )
#' setAntsrTransformParameters(tx, c(0.9,0,0,1.1,10,11) )
#' img2 = applyAntsrTransformToImage(tx, img, img)
#' # plot(img,img2)
#' @export
applyAntsrTransformToPixels <- function(transform, image, pixels.as="vector", in.place=FALSE) {
  if ( typeof(transform) == "list")
  {
    transform <- composeAntsrTransforms(transform)
  }

  if ( pixels.as=="vector" ) {
    if ( image@components != transform@dimension ) {
      stop( "Dimensionality of transform must match number of components in image for pixels.as='vector' " )
    }
  }
  else if ( pixels.as=="diffusiontensor") {
    if ( image@components != 6 ) {
      stop("diffusiontensor pixels must have 6 componenets")
    }
    if ( tranform@dimension != 3 ) {
      stop("diffusion tensorpixels must have a 3D transform")
    }
  }
  else {
    stop("Invalid 'pixels.as' type")
  }

  outImg = NA
  if ( in.place ) {
    outImg = image
  }
  else {
    outImg = antsImageClone(image)
  }

  .Call("antsrTransform_TransformPixels", transform, outImg, pixels.as, PACKAGE = "DANTsR")
  return(outImg)
}
