#' @title applyAntsrTransformToPixels
#' @description Apply transform to pixel data (vectors, tensors, etc)
#' @param transform antsrTransform
#' @param image antsImage of pixels to transform
#' @param pixel.as type of data stored as pixels ("vector", "diffusiontensor")
#' @param in.place modify the input image instead of returning a new image
#' @return antsImage
#' @examples
#' img <- antsImageRead(getANTsRData("r16"))
#' tx = new("antsrTransform", precision="float", type="AffineTransform", dimension=2 )
#' setAntsrTransformParameters(tx, c(0.9,0,0,1.1,10,11) )
#' img2 = applyAntsrTransformToImage(tx, img, img)
#' # plot(img,img2)
#' @export
applyAntsrTransformToPixels <- function(transform, image, pixel.as="vector", in.place=FALSE) {

  if ( class(transform) == "antsrTransform")
  {
    transform <- list(transform)
  }

  pixel.as = tolower(pixel.as)

  for ( tx in transform ) {
    if ( pixel.as=="vector" | pixel.as=="covariantvector") {
      if ( image@components != tx@dimension ) {
        stop( "Dimensionality of transform must match number of components in image for pixel vector types" )
      }
    }
    else if ( pixel.as=="diffusiontensor") {
      if ( image@components != 6 ) {
        stop("diffusiontensor pixels must have 6 componenets")
      }
      if ( tx@dimension != 3 ) {
        stop("diffusion tensor pixels must have a 3D transform")
      }
    }
    else {
      stop("Invalid 'pixel.as' type")
    }

    if ( tx@type == "DisplacementFieldTransform" ) {
      if ( tx@dimension != image@dimension ) {
        stop("for displacement field transforms, field and image must have the same dimension")
      }
    }
  }


  outImg = NA
  if ( in.place ) {
    outImg = image
  }
  else {
    outImg = antsImageClone(image)
  }

  for ( t in rev(transform) ) {
    .Call("antsrTransform_TransformPixels", t, outImg, pixel.as, PACKAGE="DANTsR")
  }


  if ( in.place ) {
    return(NA)
  }

  return(outImg)

}
