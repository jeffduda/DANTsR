#' @title pointCountImage
#' @description return image values at points or indices
#' @param img an antsImage
#' @param points the locations of interest
#' @param type 'point' or 'index'
#' @param interpolation options are: 'linear'
#' @export
pointCountImage = function(points, img)
{
  countImg = img*0

  idx = round(antsTransformPhysicalPointToIndex(img, points))

  for ( i in 1:dim(points)[1] ) {
    if ( prod(idx[i,]>0) & (prod(idx[i,]<=dim(img))) ) {
      countImg[ idx[i,1], idx[i,2], idx[i,3] ] = countImg[idx[i,1], idx[i,2], idx[i,3]][1] + 1
    }
  }
  return(countImg)
}