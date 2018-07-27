
#include <RcppDANTsR.h>

#include "itkImage.h"
#include "itkPointCountImageFilter.h"

template< class MeshType, class ImageType >
SEXP pointCountImage( SEXP r_mesh, SEXP r_image );

RcppExport SEXP pointCountImage( SEXP r_mesh, SEXP r_image );
