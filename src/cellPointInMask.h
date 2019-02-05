
#include <RcppDANTsR.h>

#include "itkMesh.h"
#include "itkImage.h"

template< class MeshType, class ImageType >
SEXP cellPointInMask( SEXP r_mesh, SEXP r_mask, SEXP r_index );

RcppExport SEXP cellPointInMask( SEXP r_image, SEXP r_mask, SEXP r_index );
