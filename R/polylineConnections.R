#' @title polylineConnections
#' @description return point coordinates of all voxels with specified label
#' @param mesh 'antsrMesh' of polylines
#' @param regions 'antsImage' of labeled regions
#' @param seeds vector of seed indices
#' @param method method for determining connections
#' @param max.regions maximum number of regions that a single polyline can connect
#' @export
polylineConnections = function(mesh, regions, seeds=NULL, method="forward", max.regions=2)
{

  if ( !class(mesh)=="antsrMesh") {
    stop("mesh must be an 'antsrMesh'")
  }
  if ( !class(regions)=="antsImage") {
    stop("regions must be an 'antsImage'")
  }

  if ( !is.null(seeds) ) {
    if ( length(seeds) != antsrMeshGetNumberOfCells(mesh) ) {
      stop("must have a seed for each cell in the mesh")
    }
    method="seed"
  }

  connections = list()
  for ( i in 1:antsrMeshGetNumberOfCells(mesh) ) {
    seed = -1
    if ( !is.null(seeds) ) {
      seed = seeds[i]
    }

    lineLabels = interpolateImageValues( regions, antsrMeshGetCellPoints(mesh, i), interpolation="nearestneighbor")

    if (is.null(seeds)) {
      lineLabels = lineLabels[lineLabels!=0]
      if ( length(lineLabels) > max.regions ) {
        if ( method=="forward" ) {
          lineLabels = lineLabels[1:max.regions]
        }
        else if (method=="backward") {
          lineLabels = lineLabels[ (length(lineLabels)-max.regions+1):length(lineLabels) ]
        }
      }
    }
    else {
       labels1 = lineLabels[1:seed[i]]
       labels2 = lineLabels[seed[i]:length(lineLabels)]
       labels1 = labels1[labels1!=0]
       labels2 = labels2[labels2!=0]

       lim = max.regions %/% 2
       extra = max.regions %% 2

       if ( length(labels1>lim) ) {
         labels1 = labels1[ (length(labels1)-lim+1):length(labels1) ]
       }
       if ( length(labels2>lim) ) {
         labels2 = labels2[1:lim]
       }

       lineLabels = c(labels1, labels2)
    }

    connections[[i]] = lineLabels
  }

  return(connections)

}
