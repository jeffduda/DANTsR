
#include <RcppDANTsR.h>

#include "itkImage.h"

template< class ImageType >
SEXP isInImage( SEXP r_image, SEXP r_coordinate, SEXP r_type );

RcppExport SEXP isInImage( SEXP r_image, SEXP r_coordinate, SEXP r_type );
