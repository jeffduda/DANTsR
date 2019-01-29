#' @title isInMesh
#' @description is this point or in the mesh
#' @param mesh 'antsrMesh defining a point set
#' @param coordinate point or index to test
#' @param tolerance any point with a euclidan distance less than this is considered equal

isInMesh = function( mesh, coordinate, tolerance=1e-5) {
  return (.Call( "isInMesh", mesh, coordinate, tolerance, PACKAGE="DANTsR"))
}
