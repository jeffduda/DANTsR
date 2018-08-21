#' @title extractMesh
#' @description extract a submesh
#' @param img an antsrMesh
#' @param points vector of point ids to extract
#' @param cells vector of cell ids to extract
#' @export
extractMesh = function(mesh, points=NULL, cells=NULL)
{

  if ( class(mesh)!="antsrMesh") {
    stop("Input must be an 'anstsrMesh'")
    }

  return(.Call("extractMesh", mesh, points, cells, PACKAGE="DANTsR"))
}
