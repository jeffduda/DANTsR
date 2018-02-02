# this file defines the S4 classes related to 'antsrRegistration' and the associated
# methods



#' @rdname antsrRegistration
#' @title An S4 class for registration
#'
#' @description C++ type used to represent an ITK image registration filter
#'
#' @param object input object to convert
#' @param .Object input object to convert
#' @param precision string e.g. "float" or "double"
#' @param dimension dimensionality of the transform (2,3,or 4)
#' @slot dimension usually 2 or 3 but can be 4
#' @slot precision math precision is float or double'
#' @slot pointer to the memory location of the itk object
setClass(Class = "antsrRegistration",
         representation(precision= "character", dimension = "integer",
           components = "integer", pixeltype="character",
           smoothing = "numeric", transform="character",
           optimizer="character", interpolator="character",
           sampling="numeric", metric="character",
           iterations="numeric", shrink="numeric",
           samplingStrategy="character", convergence="numeric",
           smoothingPhysical="integer", nLevels="integer" ),
         prototype( precision="float", dimension=as.integer(2),
           components=as.integer(1), pixeltype="float",
           smoothing=c(1), transform="translation",
           optimizer="regularStepGradientDescent", interpolator="linear",
           sampling=0.20, metric="meanSquares",
           iterations=c(20), shrink=c(1),
           samplingStrategy="regular", convergence=c(1e-6),
           smoothingPhysical=as.integer(0), nLevels=as.integer(1) ) )

#' @rdname antsrRegistration
#' @aliases show,antsrRegistration-method
setMethod(f = "show", "antsrRegistration", function(object){
    cat("antsrRegistration\n")
    cat("  Dimension    :", object@dimension, "\n")
    cat("  Precision    :", object@precision, "\n")
    cat("  Iterations:  :", object@iterations, "\n")
    cat("  Components   :", object@components, "\n")
    cat("  PixelType    :", object@pixeltype, "\n")
    cat("  Smoothing    :", object@smoothing, "\n")
    cat("  Transform    :", object@transform, "\n")
    cat("  Optimizer    :", object@optimizer, "\n")
    cat("  Interpolator :", object@interpolator, "\n")
    cat("  Sampling     :", object@sampling, "\n")
    cat("\n")
})

#' @rdname antsrRegsitration
#' @aliases initialize,antsrRegistration-method
#setMethod(f = "initialize", signature(.Object = "antsrRegistration"), definition = function(.Object,
#  dimension = 3, precision = "float") {
#  reg = .Call("antsrRegistration", precision, dimension, PACKAGE = "DANTsR")
#  return( reg )
#})


#' @title antsrRegistrationCreate
#' @description create object that registers two images
#' @param dimension number of dimensions
#' @param percision use 'float' or 'double' for values
#' @param reserve number of points to allocate on creation
#' @return antsrRegistration
#' @examples
#' x =  antsrRegistrationCreate( 3, "float", reserve=128 )
#' @export
#antsrRegsitrationCreate <- function(dimension=3, precision="float", reserve=0)
#{

  # Check for valid dimension
#  if ( (dimension < 2) | (dimension > 4) )
#  {
#    stop(paste("Unsupported dimension:", dimension))
#  }

#  if ( (precision != "float") & (precision != "double"))
#  {
#    stop(paste("Unsupported precision:", precision))
#  }

#  if ( reserve < 0 ) {
#    stop(paste("Unsupported reserve number:", reserve))
#  }

#  reg = .Call("antsrRegistration", precision, dimension, reserve, PACKAGE = "DANTsR")

#  return(reg)
#  }


antsrRegistrationRun <- function( stages, fixed, moving ) {

  # Check parameters for validity


  .Call( "antsrRegistrationRun", stages, fixed, moving, package="DANTsR")
  #return(NA)

}
