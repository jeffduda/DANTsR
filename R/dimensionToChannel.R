#' @title dimensionToChannel
#' @description create multichannel image by reducing dimension by one
#' @param img an antsImage
#' @param dimension the dimension that is converted to a channel
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

#' @title channelToDimension
#' @description convert channels to a dimension
#' @param img a multi component antsImage
#' @param dimension the channels are converted to
#' @param spacing spacing for new dimension
#' @param origin origin for new dimension
#' @param dircol column values for direction matrix in new dimension
#' @param dirrow row values for direction matrix in new dimension
#' @export
channelToDimension = function(img, dimension=(img@dimension+1), spacing=1, origin=0, dircol=NULL, dirrow=NULL)
{
  if ( !class(img)=="antsImage") {
    stop("Must input an 'antsImage'")
  }

  if ( !img@isVector ) {
    stop("Input image must be a vector image")
  }

  #newdims = append(1:img@dimension, img@dimension+1, dimension-1)
  sp = append(antsGetSpacing(img), spacing, dimension-1)
  og = append(antsGetOrigin(img), origin, dimension-1)

  if (is.null(dircol)) {
    dircol = rep(0, img@dimension+1)
    dircol[dimension] = 1
  }
  if (is.null(dirrow)) {
    dirrow = rep(0, img@dimension+1)
    dirrow[dimension] = 1
  }

  dat = as.array(img)
  ap = 1:(img@dimension+1)
  if ( dimension != 1 ) {
    ap = 2:(img@dimension+1)
    ap = append(ap, 1, dimension-1)
    dat = aperm(dat, ap)
  }

  dir = antsGetDirection(img)
  dir = rbind(rep(0,img@dimension), dir )
  dir = cbind(rep(0,img@dimension+1), dir)
  dir = dir[ap,ap]
  dir[dimension,] = dirrow
  dir[,dimension] = dircol

  mimg = as.antsImage(dat, components=F)
  antsSetSpacing( mimg, sp )
  antsSetOrigin( mimg, og )
  antsSetDirection( mimg, dir )

  return(mimg)
}
