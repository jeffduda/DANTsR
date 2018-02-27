#' @title dtiReconstruction
#' @description calculate an anistropy image from DTI
#' @param img a 6-channel antsImage or and eigen sytem list
#' @param method the summary meaure to return
#' \itemize{
#'   \item{Trace}{Sum of eigenvalues}
#'   \item{MeanDiffusion}{Average eigenvalue}
#'   \item{AxialDiffusion}{Largest eigenvalue}
#'   \item{RadialDiffusion}{Mean of two smallest eigenvalues}
#' }
#' @param method the reconstruction algorithm to use (default = "itk-svd")
#' \itemize{
#'   \item{itk-svd}{uses the itk::DiffusionTensor3DReconstructionImageFilter filter}
#' }
#' @export
dtiReconstruction = function(dwi, gradients, method)
{

  method = tolower(method)
  .Call( "dtiReconstruction", dwi, gradients, method, package="DANTsR")

}
