sparseDecomR.InitializeVector = function( size, seed ) {
  return( rnorm(size) )
}


sparseDecomR.OrthogonalizeVector = function( vector, reference ) {

  ipv = sum(reference*reference)
  if (ipv == 0) {
    return(vector)
  }

  ratio = sum(vector*reference)/ipv
  return( vector - (reference*ratio) )

}

# FIXME - add surface image smoothing option
spareDecomR.SpatiallySmoothVector = function( vector, mask, smoothing ) {
  img = mask*0
  img[mask > 0] = vector
  img = smoothImage(img, smoothing)
  return( img[mask>0] )
}

sparseDecomR = function( matrix, nvecs, sparseness, mask ) {

  smoothing = 1.0

  sparseParam = rep(sparseness, nvecs)

  # Data driven initialization
  variates = matrix(0, ncol(matrix), nvecs)
  for ( i in 1:nvecs ) {
    vec = sparseDecomOne.InitializeVector(ncol(matrix))
    if ( i > 1) for ( j in 1:(i-1) ) {
      vec = sparseDecomOne.OrthogonalizeVector( vec, variates[,j])
    }
    vec = sparseDeconOne.SpatiallySmoothVector( vec, mask, smoothing)
    vec = as.vector(scale(vec))
    variates[,j] = vec
  }




}
