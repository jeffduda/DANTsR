
#include <RcppDANTsR.h>

#include "itkImage.h"
#include "itkPointCountImageFilter.h"

template< class InputMeshType, class OutputMeshType >
SEXP polylineABSplineSmoothing( SEXP r_mesh );

RcppExport SEXP polylineArBSplineSmoothing( SEXP r_mesh );
