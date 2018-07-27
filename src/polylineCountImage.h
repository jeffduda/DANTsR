
#include <RcppDANTsR.h>

#include "itkImage.h"
#include "itkPolyLineMeshToCellCountImageFilter.h"

template< class MeshType, class ImageType >
SEXP polylineCountImage( SEXP r_mesh, SEXP r_image );

RcppExport SEXP polylineCountImage( SEXP r_mesh, SEXP r_image );
