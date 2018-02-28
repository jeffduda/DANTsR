#' @title dtiColorMap
#' @description calculate an RGB map of primary direction weighted by anisotropy
#' @param img a 6-channel antsImage
#' @param weight method used to weight color values'
#' @export
dtiColorMap = function(img, weight="FractionalAnisotropy")
{
  if ( !class(img)=="antsImage") {
    stop("Must input an 'antsImage'")
  }

  if ( img@components != 6 ) {
    stop("Input image must have 6 components")
  }

  weight = tolower(weight)

  fa = dtiAnisotropy(img, method=weight)
  fa[fa<0]=0
  fa[fa>1]=1

  eig = dtiEigenSystem(img)

  rgb = 255*mergeChannels(list(fa,fa,fa))*abs(eig$v1)

}
