#' @title cellsHitTarget
#' @description return indices of cells that hit the target ROI
#' @param mesh an anstsrMesh
#' @param target antsImage with target ROI
#' @param value value of the labeled ROI
#' @export
cellsHitTarget = function(mesh, target, value=1)
{

  hits = .Call("cellsHitTarget", mesh, target, value, PACKAGE="DANTsR")
  return(hits)
}
