#' @title dtiEigenSystem
#' @description calculate eigen values and vectors
#' @param img a 6-channel antsImage
#' @param method method used to calculate the Eigen system, default is 'itk-svd'
#' \itemize{
#'   \item{itk-svd}{Uses ITK based call to netlib/vnl SVD. Is faster.}
#'   \item{eigen}{Uses R's 'eigen' function}
#' }
#' @export
dtiEigenSystem = function(img, method="itk-svd")
{
  if ( !class(img)=="antsImage") {
    stop("Must input an 'antsImage'")
  }

  if ( img@components != 6 ) {
    stop("Input image must have 6 components")
  }

  if ( method=="eigen") {
    dta = as.array(img)
    d = dim(dta)
    dim(dta) = c(d[1],d[2]*d[3]*d[4])

    all_eigs = apply(dta, 2, function(x) {
      tens = diag(3)
      tens[lower.tri(tens, diag=TRUE)] = x
      tens = tens + t(tens) - diag(diag(tens))
      #tens = matrix(0,3,3)
      #tens[1,1] =
      eigen( tens )
      } )

    eigValues = lapply(all_eigs, function(x) {
      x$values
    })

    # Get eigenvalues into antsImage format
    eigValues = lapply(eigValues, unlist)
    eigValues = do.call(rbind, eigValues)
    e1 = eigValues[,1]
    e2 = eigValues[,2]
    e3 = eigValues[,3]
    dim(e1) = dim(img)
    dim(e2) = dim(img)
    dim(e3) = dim(img)
    # FIXME - passing a reference image fails for multichannel
    e1 = antsCopyImageInfo(img, as.antsImage(e1))
    e2 = antsCopyImageInfo(img, as.antsImage(e2))
    e3 = antsCopyImageInfo(img, as.antsImage(e3))

    # Get eigenvectors into antsImage format
    eig_vecs = lapply(all_eigs, function(x) {
      x$vectors
    })
    eig_vecs = lapply(eig_vecs, c)
    eig_vecs = do.call(rbind, eig_vecs)

    v1 = aperm(eig_vecs[,1:3], c(2,1))
    v2 = aperm(eig_vecs[,4:6], c(2,1))
    v3 = aperm(eig_vecs[,7:9], c(2,1))
    dim(v1) = c(3,dim(img))
    dim(v2) = c(3,dim(img))
    dim(v3) = c(3,dim(img))

    # FIXME - passing a reference image fails for multichannel
    v1 = antsCopyImageInfo(img,as.antsImage(v1, components=T))
    v2 = antsCopyImageInfo(img,as.antsImage(v2, components=T))
    v3 = antsCopyImageInfo(img,as.antsImage(v3, components=T))

    dt = list(e1=e1, e2=e2, e3=e3, v1=v1, v2=v2, v3=v3)

  }
  else {
   dt = .Call( "dtiFilters", img, "eigensystem", package="DANTsR")
 }

 return(dt)

}
