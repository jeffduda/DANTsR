
#include <RcppDANTsR.h>

#include "itkImage.h"
#include "itkMeshExtractionFilter.h"

template< class InputMeshType, class OutputMeshType >
SEXP extractMesh( SEXP r_mesh, SEXP r_points, SEXP r_cells );

RcppExport SEXP extractMesh( SEXP r_mesh, SEXP r_points, SEXP r_cells );
