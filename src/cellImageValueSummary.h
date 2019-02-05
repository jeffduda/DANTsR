
#include <RcppDANTsR.h>

#include "itkMesh.h"
#include "itkImage.h"

template< class MeshType, class ImageType, class TOutput >
SEXP cellImageValueSummary( SEXP r_mesh, SEXP r_mask, SEXP r_measure, SEXP r_index );

RcppExport SEXP cellImageValueSummary( SEXP r_image, SEXP r_mask, SEXP r_measure, SEXP r_index );
