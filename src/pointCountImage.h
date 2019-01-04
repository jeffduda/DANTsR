#include <RcppDANTsR.h>

template< class MeshType, class ImageType >
SEXP pointCountImageFunction( SEXP r_mesh, SEXP r_image );

RcppExport SEXP pointCountImage( SEXP r_mesh, SEXP r_image );
