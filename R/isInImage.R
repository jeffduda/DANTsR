
#' @title pointIsInImage
#' @description is this point in the image space
#' @param image 'antsImage' defining a physical space
#' @param point spatial point
pointIsInImage = function( image, point ) {
  return( isInImage(image, point, type="point") )
}

#' @title indexIsInImage
#' @description is this index in the image space
#' @param image 'antsImage' defining a physical space
#' @param index location to test
indexIsInImage = function( image, index ) {
  return( isInImage(image, index, type="index") )
}

#' @title isInImage
#' @description is this point or index in the image space
#' @param image 'antsImage' defining a physical space
#' @param coordinate point or index to test
#' @param type either "point" or "index"

isInImage = function( image, coordinate, type="point") {
  return (.Call( "isInImage", image, coordinate, type, PACKAGE="DANTsR"))
}
