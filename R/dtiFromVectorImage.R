  #' @title dtiFromVectorImage
  #' @description create a diffusion image from a 6-channel antsImage
  #' @param vImg a 6-channel antsImage
  #' @export
dtiFromVectorImage = function( vImg ) {

  if ( !(class(vImg)=='antsImage') ) {
    stop("Input must be an 'antsImage'")
  }
  if ( vImg@components != 6 ) {
    stop("Input must have 6 channels")
  }

  x = .Call( "dtiFromVectorImage", vImg, PACKAGE="DANTsR" )

  return(0)

}
