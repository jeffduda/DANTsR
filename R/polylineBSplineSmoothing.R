#' @title polylineBSplineSmoothing
#' @description parameterize polyline by arc-length, coordinates save as point data
#' @param img an antsrMesh
#' @export
polylineBSplineSmoothing = function(mesh)
{

  if ( class(mesh)!="antsrMesh") {
    stop("Input must be an 'anstsrMesh'")
    }

  return(.Call("polylineBSplineSmoothing", mesh, package="DANTsR"))
}
