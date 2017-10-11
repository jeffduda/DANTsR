#ifndef __RCPPDANTSR_H
#define __RCPPDANTSR_H

#include "itkMacro.h"
#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkVector.h"
#include "itkMesh.h"

#include <RcppCommon.h>

// itk::Image to antsImage
#include <RcppANTsR/RcppANTsR_ImageDef.h>
#include <RcppANTsR/RcppANTsR_VectorImageDef.h>
#include <RcppDANTsR/RcppDANTsR_MeshDef.h>

// This needs to go after wrap declarations and before implementations
#include <Rcpp.h>

#include <RcppANTsR/RcppANTsR_ImageImp.h>
#include <RcppANTsR/RcppANTsR_VectorImageImp.h>
#include <RcppDANTsR/RcppDANTsR_MeshImp.h>

#endif
