#' @title dimensionToChannel
#' @description create multichannel image by reducing dimension by one
#' @param img an antsImage
#' @export
dimensionToChannel = function(img, dimension=length(dim(img)))
{
  if ( !class(img)=="antsImage") {
    stop("Must input an 'antsImage'")
  }

  if ( img@isVector ) {
    stop("Input image must be a scalar image")
  }

  dim = 1:length(dim(img))
  newdim = dim[which(dim!=dimension)]

  mimg = as.antsImage(aperm(as.array(img), c(dimension, newdim)), components=T)
  antsSetSpacing( mimg, antsGetSpacing(img)[newdim])
  antsSetOrigin( mimg, antsGetOrigin(img)[newdim])
  antsSetDirection( mimg, antsGetDirection(img)[newdim, newdim] )

  return(mimg)
}
