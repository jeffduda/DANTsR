#' @title dtiConvert
#' @description convert between package formats
#' @param img a DT image to convert
#' @param input the format of the input image
#' \itemize{
#'   \item{ANTsR}{}
#'   \item{MRtrix}{}
#' }
#' @param output the format of the returned image
#' \itemize{
#'   \item{ANTsR}{}
#'   \item{MRTrix}{}
#' }
#' @export
dtiConvert = function(img, input, output="ANTsR")
{
  if ( input==output ) {
    return(input)
  }

  dt = NA
  if ( output=="ANTsR" ) {
    if ( input=="MRTrix") {
      x = splitChannels(dimensionToChannel(img))
      dt = mergeChannels(list( x[[1]], x[[4]], x[[5]], x[[2]], x[[6]], x[[3]]))
    }
    else {
      stop( "Unsupported input image format")
    }

  }
  else if ( output=="MRTrix") {
    if ( input=="ANTsR" ) {
      x = splitChannels(img)
      dt = channelToDimension(mergeChannels(list(x[[1]], x[[4]], x[[6]], x[[2]], x[[3]], x[[5]])))
    }
    else {
      stop("Unsupported input image format")
    }

  }
  else {
    stop("Unsupported output image format")
  }

  return(dt)

}
