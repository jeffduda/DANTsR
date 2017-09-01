#ifndef __RCPPDANTSR_H
#define __RCPPDANTSR_H

#include "itkMacro.h"
#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkVector.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "vnl/vnl_vector_ref.h"
#include "itkTransform.h"
#include "itkAffineTransform.h"

#include <RcppCommon.h>

// itk::Image to antsImage
#include <RcppANTsR/RcppANTsR_VectorImageDef.h>
#include <RcppDANTsR/RcppDANTsR_DTIImageDef.h>

// This needs to go after wrap declarations and before implementations
#include <Rcpp.h>

#include <RcppANTsR/RcppANTsR_VectorImageImp.h>
#include <RcppDANTsR/RcppDANTsR_DTIImageImp.h>

#endif
