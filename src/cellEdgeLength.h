
#include <RcppDANTsR.h>

#include "itkMesh.h"
#include "itkEdgeLengthCellFunction.h"

template< class MeshType >
SEXP cellEdgeLength( SEXP r_mesh, SEXP r_index );

RcppExport SEXP cellEdgeLength( SEXP r_image, SEXP r_index );
