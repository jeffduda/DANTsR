getConnectivityMatrix = function( mesh, labelImage, labels=NULL, measure=NULL ) {

  if (is.null(labels) ) {
    labels = unique(as.vector(as.array(labelImage)))
    labels = labels[ labels != 0 ]
    labels = labels[order(labels)]
  }

  values=NULL
  if ( !is.null(measure) ) {
    values = cellImageValueSummary( mesh, measure, "mean" )
  }

  nLabels = length(labels)
  cmat = matrix(0, nLabels, nLabels)
  for ( i in 1:nLabels ) {
    for ( j in (i+1):nLabels ) {
      l1 = labels[i]
      l2 = labels[j]

      hits = cellsConnectTargets( mesh, labelImage, l1, l2 )
      if ( is.null(values) ) {
        cmat[i,j] = length(hits)
      } else {
        cmat[i,j] = mean( values[hits] )
      }
      cmat[j,i] = cmat[i,j]
      print(paste(i,j,cmat[i,j]))
    }
  }

  return(cmat)

}
