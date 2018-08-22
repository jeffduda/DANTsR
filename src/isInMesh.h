#include <RcppDANTsR.h>

template< class MeshType >
SEXP isInMesh( SEXP r_mesh, SEXP r_coordinate );

RcppExport SEXP isInMesh( SEXP r_mesh, SEXP r_coordinate, SEXP r_tolerance );
