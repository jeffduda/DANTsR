#' @title dtiReconstruction
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
dtiReconstruction = function(dwi, gradients, method)
{

  method = tolower(method)
  .Call( "dtiReconstruction", dwi, gradients, method, package="DANTsR")

}
