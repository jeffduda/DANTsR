#' @title dtiReconstruction
#' @description calculate an anistropy image from DTI
#' @param dwi an N-channel antsImage
#' @param gradients Nx4 matrix of gradient directions and b-values
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
