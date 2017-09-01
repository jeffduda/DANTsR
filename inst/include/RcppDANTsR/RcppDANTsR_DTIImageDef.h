#ifndef __RCPPANTSR_IMAGE_H
#define __RCPPANTSR_IMAGE_H

#include "itkMacro.h"
#include "itkImage.h"
#include "itkDiffusionTensor3D.h"
#include <RcppCommon.h>

namespace Rcpp {

typedef itk::DiffusionTensor3D<float> DTFloat;
typedef itk::DiffusionTensor3D<double> DTDouble;

// itk::Image to antsImage
template <> inline SEXP wrap( const itk::Image<DTDouble,3>::Pointer &image );
template <> inline SEXP wrap( const itk::Image<DTFloat,3>::Pointer &image );

// antsImage to itk::Image
template <> inline itk::Image<DTDouble,3>::Pointer as( SEXP itkImageR );
template <> inline itk::Image<DTFloat,3>::Pointer as( SEXP itkImageR );

}

#endif
