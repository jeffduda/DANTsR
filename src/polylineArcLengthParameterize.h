
#include <RcppDANTsR.h>

#include "itkImage.h"
#include "itkPointCountImageFilter.h"

template< class InputMeshType, class OutputMeshType >
SEXP polylineArcLengthParameterize( SEXP r_mesh );

RcppExport SEXP polylineArcLengthParameterize( SEXP r_mesh );
