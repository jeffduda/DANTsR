#' @title polylineVoxelCount
#' @description return point coordinates of all voxels with specified label
#' @param img an antsImage
#' @param reference antsImage that defines the physical space to consider
#' @export

isInImage = function( image, index ) {
  d = dim(image)
  for ( i in 1:length(index) ) {
    if ( index[i] < 1 ) {
      return(FALSE)
    }
    else if ( index[i] > (d[i]+1) ) {
      return(FALSE)
    }
  }
  return(TRUE)
}

polylineVoxelCount = function(mesh, reference)
{

  if ( !class(mesh)=="antsrMesh") {
    stop("mesh must be an 'antsrMesh'")
  }
  if ( !class(reference)=="antsImage") {
    stop("reference must be an 'antsImage'")
  }

  count = reference*0

  connections = list()
  for ( i in 1:antsrMeshGetNumberOfCells(mesh) ) {
    pts = antsrMeshGetCellPoints(mesh,i)
    idx = antsTransformPhysicalPointToIndex(count, pts)
    for ( j in 1:dim(idx)[1] ) {
      skip = FALSE
      if ( j > 1 ) {
        skip = sum(idx[j,]==idx[j-1,])==dim(idx)[2]
      }
      if ( !skip ) {
        if ( isInImage(count, idx) ) {
          count[ idx[j,1], idx[j,2], idx[j,3] ] = count[ idx[j,1], idx[j,2], idx[j,3] ][1] + 1
        }
        else {
          print(paste("Invalid point/index: ", point, " / ", idx[j,]))
        }
      }
    }

  }

  return(count)

}
