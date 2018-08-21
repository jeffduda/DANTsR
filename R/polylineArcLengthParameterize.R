#' @title polylineArcLengthParameterize
#' @description parameterize polyline by arc-length, coordinates save as point data
#' @param img an antsrMesh
#' @export
polylineArcLengthParameterize = function(mesh)
{

  if ( class(mesh)!="antsrMesh") {
    stop("Input must be an 'anstsrMesh'")
    }

  return(.Call("polylineArcLengthParameterize", mesh, package="DANTsR"))
}
