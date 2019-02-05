
#include <RcppDANTsR.h>

#include "itkMesh.h"
#include "itkImage.h"

template< class MeshType, class ImageType >
SEXP streamlineTargetsFromSeed( SEXP r_mesh, SEXP r_image, SEXP r_seeds, SEXP r_index );

RcppExport SEXP streamlineTargetsFromSeed( SEXP r_mesh, SEXP r_image, SEXP r_seeds, SEXP r_index );
