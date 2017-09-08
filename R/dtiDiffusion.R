#' @title dtiDiffusion
#' @description calculate an anistropy image from DTI
#' @param img a 6-channel antsImage or and eigen sytem list
#' @param method the summary meaure to return
#' \itemize{
#'   \time{Trace}{Sum of eigenvalues}
#'   \item{MeanDiffusion}{Average eigenvalue}
#'   \item{AxialDiffusion}{Largest eigenvalue}
#'   \item{RadialDiffusion}{Mean of two smallest eigenvalues}
#' }
#' @param method which type of anisotropy index to calculate
#' @export
dtiDiffusion = function(img, method="MeanDiffusion")
{

  method = tolower(method)
  if ( class(img) == "antsImage" ) {
    img = dtiEigenSystem(img)
  }

  if ( class(img)=="list" ) {
    if (method=="trace") {
      return (img$e1 + img$e2 + img$e3)
    }
    else if (method=="meandiffusion") {
      return (img$e1 + img$e2 + img$e3)/3
    }
    else if (method=="axialdiffusion") {
      return (img$e1)
    }
    else if (method=="radialdiffusion") {
      return (img$e2 + img$e3)/2
    }
    else {
      stop("Unknown method passed")
    }
  }
  else {
    stop("Input must be 'antsImage' or Eigen system list")
  }

  #.Call( "dtiFilters", img, method, package="DANTsR")

}
