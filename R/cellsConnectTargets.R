#' @title cellsConnectTargets
#' @description return indices of cells that hit the target ROI
#' @param mesh an anstsrMesh
#' @param target antsImage with target ROI
#' @param value1 value of first target
#' @param value2 value of second target
#' @export
cellsConnectTargets = function(mesh, target, value1, value2)
{
  hits = .Call("cellsConnectTargets", mesh, target, value1, value2, PACKAGE="DANTsR")
  return(hits)
}
