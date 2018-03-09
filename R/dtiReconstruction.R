dtiReconstruction.r.svd <- function(x, basis, bMat) {
  if ( dim(bMat)[2] > 6 ) {
    x = bMat %*% x
  }
  solve(basis, x)
}

#' @title dtiReconstruction
#' @description calculate an anistropy image from DTI
#' @param dwi an N-channel antsImage
#' @param gradients Nx4 matrix of gradient directions and b-values
#' @param method the reconstruction algorithm to use (default = "itk-svd")
#' \itemize{
#'   \item{itk-svd}{uses the itk::DiffusionTensor3DReconstructionImageFilter filter}
#' }
#' @export
dtiReconstruction = function(dwi, gradients, method, mask=NULL)
{

  method = tolower(method)

  if ( method=="r-svd" ) {
    dat = channelToDimension(dwi)

    if ( is.null(mask) ) {
      mask = extractSlice(dat, 1, 4)*0+1
    }

    mat = timeseries2matrix(dat, mask)

    bmat = gradients[which(gradients[,4]!=0),1:3]
    bmat = cbind( bmat[,1]*bmat[,1],
                  2*bmat[,1]*bmat[,2],
                  2*bmat[,1]*bmat[,3],
                  bmat[,2]*bmat[,2],
                  2*bmat[,2]*bmat[,3],
                  bmat[,3]*bmat[,3])

    bvalMat = gradients[which(gradients[,4]!=0)]
    bValue = bvalMat[1]

    tensorBasis=NA
    if (dim(bmat)[1] > 6) {
      tensorBasis = t(bmat) %*% bmat
    }
    else {
      tensorBasis = bmat
    }
    bmat = t(bmat)

    bvals = gradients[,4]
    b0ids = which(bvals==0)

    b0=mat[b0ids,]
    if ( length(b0ids) > 1 ) {
      b0 = colMeans(b0)
    }
    mat = mat[-b0ids,]

    invalidB0 = which(b0==0)
    b0[invalidB0] = 1

    mat = apply(mat, 1, function(x) ( -log(x/b0)/bValue ) )
    mat = apply(mat, 1, dtiReconstruction.r.svd, basis=tensorBasis, bMat=bmat)

    mat

  }
  else {
    .Call( "dtiReconstruction", dwi, gradients, method, package="DANTsR")
  }

}
