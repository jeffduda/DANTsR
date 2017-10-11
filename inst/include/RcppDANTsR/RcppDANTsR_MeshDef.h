#ifndef __RCPPANTSR_MESH_H
#define __RCPPANTSR_MESH_H

#include "itkMacro.h"
#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkVector.h"
#include "itkMesh.h"
#include <RcppCommon.h>

namespace Rcpp {

// itk to antsr
template <> inline SEXP wrap( const itk::Mesh< double, 2 >::Pointer & mesh );
template <> inline SEXP wrap( const itk::Mesh< double, 3 >::Pointer & mesh );
template <> inline SEXP wrap( const itk::Mesh< double, 4 >::Pointer & mesh );

template <> inline SEXP wrap( const itk::Mesh< float, 2 >::Pointer & mesh );
template <> inline SEXP wrap( const itk::Mesh< float, 3 >::Pointer & mesh );
template <> inline SEXP wrap( const itk::Mesh< float, 4 >::Pointer & mesh );

// antsr to itk
template <> inline itk::Mesh< double, 2 >::Pointer as( SEXP itkMesh );
template <> inline itk::Mesh< double, 3 >::Pointer as( SEXP itkMesh );
template <> inline itk::Mesh< double, 4 >::Pointer as( SEXP itkMesh );

template <> inline itk::Mesh< float, 2 >::Pointer as( SEXP itkMesh );
template <> inline itk::Mesh< float, 3 >::Pointer as( SEXP itkMesh );
template <> inline itk::Mesh< float, 4 >::Pointer as( SEXP itkMesh );

}


#endif
