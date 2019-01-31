
#include <algorithm>
#include <vector>
#include <string>
#include <RcppDANTsR.h>

#include "antsUtilities.h"
#include "itkMesh.h"
#include "itkByteSwapper.h"
#include "itkCaminoStreamlineFileReader.h"
#include "itkCaminoStreamlineFileWriter.h"
#include "itkMRTrixStreamlineFileReader.h"
#include "itkTrackVisStreamlineFileReader.h"
#include "itkTrackVisStreamlineFileWriter.h"
#include "itkVtkPolyDataFileReader.h"
#include "itkVtkPolyDataFileWriter.h"

#include "itkLineCell.h"
#include "itkPolygonCell.h"
#include "itkPolyLineCell.h"

template< class MeshType >
typename MeshType::Pointer antsrMesh( itk::IdentifierType reserve, SEXP r_points )
{

  typename MeshType::Pointer mesh = MeshType::New();
  typename MeshType::PointType itkPoint;

  Rcpp::NumericMatrix points(r_points);
  //Rcpp::NumericMatrix cells(r_cells);

  if ( reserve > 0 ) {
    mesh->GetPoints()->Reserve(reserve);
  }

  if ( points.ncol() > 1 ) {

    for ( unsigned int i=0; i<points.nrow(); i++ ) {

      for ( unsigned int j=0; j<points.ncol(); j++ ) {
        itkPoint[j] = points(i,j);
      }

      mesh->SetPoint(i, itkPoint);
    }
  }

  // options here are:
  // CellsAllocationMethodUndefined
  // CellsAllocatedAsStaticArray
  // CellsAllocatedAsADynamicArray
  // CellsAllocatedDynamicallyCellByCell
  mesh->SetCellsAllocationMethod( MeshType::CellsAllocatedDynamicallyCellByCell );

  //Rcpp::Rcout << "Created mesh: " << mesh->GetNumberOfPoints() << std::endl;

  return mesh;

}

//pixeltype, precision, dimension, type, isVector
RcppExport SEXP antsrMesh( SEXP r_precision, SEXP r_dimension, SEXP r_reserve, SEXP r_points)
{
try
{
  std::string precision = Rcpp::as< std::string >( r_precision );
  unsigned int dimension = Rcpp::as< int >( r_dimension );
  itk::IdentifierType reserve = Rcpp::as< itk::IdentifierType >( r_reserve );

  if ( (dimension < 2) || (dimension > 4) ) {
    Rcpp::stop("Unsupported dimension type - must be 2,3, or 4");
  }

  if ( precision=="double") {
    using PrecisionType = double;
    if ( dimension == 2 ) {
      using MeshType = itk::Mesh<PrecisionType,2>;
      return Rcpp::wrap( antsrMesh<MeshType>(reserve, r_points) );
    }
    else if ( dimension == 3 ) {
      using MeshType = itk::Mesh<PrecisionType,3>;
      return Rcpp::wrap( antsrMesh<MeshType>(reserve, r_points) );
    }
    else if ( dimension == 4 ) {
      using MeshType = itk::Mesh<PrecisionType,4>;
      return Rcpp::wrap( antsrMesh<MeshType>(reserve, r_points) );
    }
  }
  else if (precision=="float") {
    using PrecisionType = float;
    if ( dimension == 2 ) {
      using MeshType = itk::Mesh<PrecisionType,2>;
      return Rcpp::wrap( antsrMesh<MeshType>(reserve, r_points) );
    }
    else if ( dimension == 3 ) {
      using MeshType = itk::Mesh<PrecisionType,3>;
      return Rcpp::wrap( antsrMesh<MeshType>(reserve, r_points) );
        }
    else if ( dimension == 4 ) {
      using MeshType = itk::Mesh<PrecisionType,4>;
      return Rcpp::wrap( antsrMesh<MeshType>(reserve, r_points) );
    }
  }
  else {
    Rcpp::stop( "Unsupported precision type - must be 'float' or 'double'");
  }

  // Never reached
  return( Rcpp::wrap(NA_REAL) );

}
catch( itk::ExceptionObject & err )
  {
  Rcpp::Rcout << "ITK ExceptionObject caught !" << std::endl;
  Rcpp::Rcout << err << std::endl;
  Rcpp::stop("ITK exception caught");
  }
catch( const std::exception& exc )
  {
  forward_exception_to_r( exc ) ;
  }
catch(...)
  {
	Rcpp::stop("c++ exception (unknown reason)");
  }
return Rcpp::wrap(NA_REAL); //not reached
}

void antsrMesh_Valid( SEXP r_mesh )
{
  Rcpp::S4 rMesh( r_mesh );
  std::string precision = rMesh.slot("precision");
  unsigned int dimension = rMesh.slot("dimension");

  if ( (dimension < 2) || (dimension > 4) ) {
    Rcpp::stop("Unsupported dimension type - must be 2,3, or 4");
  }
  if ( (precision != "float") && (precision != "double")) {
    Rcpp::stop("Unsupported precision type - must be 'float' or 'double'");
  }
}

std::string antsrMesh_GetPrecision( SEXP r_mesh )
{
  Rcpp::S4 rMesh(r_mesh);
  return rMesh.slot("precision");
}

unsigned int antsrMesh_GetDimension( SEXP r_mesh )
{
  Rcpp::S4 rMesh(r_mesh);
  return rMesh.slot("dimension");
}

template< class MeshType >
SEXP
antsrMesh_GetNumberOfPoints( SEXP rMesh )
{
  using MeshPointerType = typename MeshType::Pointer;
  MeshPointerType mesh = Rcpp::as<MeshPointerType>( rMesh );
  return( Rcpp::wrap(mesh->GetNumberOfPoints()) );
}

//pixeltype, precision, dimension, type, isVector
RcppExport SEXP antsrMesh_GetNumberOfPoints( SEXP r_mesh )
{
try
{
  antsrMesh_Valid(r_mesh);
  std::string precision = antsrMesh_GetPrecision(r_mesh);
  unsigned int dimension = antsrMesh_GetDimension(r_mesh);

  if ( precision=="double") {
    using PrecisionType = double;
    if ( dimension == 2 ) {
      using MeshType = itk::Mesh<PrecisionType,2>;
      return antsrMesh_GetNumberOfPoints<MeshType>(r_mesh);
    }
    else if ( dimension == 3 ) {
      using MeshType = itk::Mesh<PrecisionType,3>;
      return antsrMesh_GetNumberOfPoints<MeshType>(r_mesh);
    }
    else if ( dimension == 4 ) {
      using MeshType = itk::Mesh<PrecisionType,4>;
      return antsrMesh_GetNumberOfPoints<MeshType>(r_mesh);
    }
  }
  else if ( precision=="float") {
    using PrecisionType = float;
    if ( dimension == 2 ) {
      using MeshType = itk::Mesh<PrecisionType,2>;
      return antsrMesh_GetNumberOfPoints<MeshType>(r_mesh);
    }
    else if ( dimension == 3 ) {
      using MeshType = itk::Mesh<PrecisionType,3>;
      return antsrMesh_GetNumberOfPoints<MeshType>(r_mesh);
        }
    else if ( dimension == 4 ) {
      using MeshType = itk::Mesh<PrecisionType,4>;
      return antsrMesh_GetNumberOfPoints<MeshType>(r_mesh);
    }
  }
  else {
    Rcpp::stop( "Unsupported precision type - must be 'float' or 'double'");
  }

  // Never reached
  return( Rcpp::wrap(NA_REAL) );

}
catch( itk::ExceptionObject & err )
  {
  Rcpp::Rcout << "ITK ExceptionObject caught !" << std::endl;
  Rcpp::Rcout << err << std::endl;
  Rcpp::stop("ITK exception caught");
  }
catch( const std::exception& exc )
  {
  forward_exception_to_r( exc ) ;
  }
catch(...)
  {
	Rcpp::stop("c++ exception (unknown reason)");
  }
return Rcpp::wrap(NA_REAL); //not reached
}

template< class MeshType >
SEXP
antsrMesh_GetNumberOfCells( SEXP rMesh )
{
  using MeshPointerType = typename MeshType::Pointer;
  MeshPointerType mesh = Rcpp::as<MeshPointerType>( rMesh );
  return( Rcpp::wrap(mesh->GetNumberOfCells()) );
}

//pixeltype, precision, dimension, type, isVector
RcppExport SEXP antsrMesh_GetNumberOfCells( SEXP r_mesh )
{
try
{
  Rcpp::S4 rMesh( r_mesh );
  std::string precision = rMesh.slot("precision");
  unsigned int dimension = rMesh.slot("dimension");

  if ( (dimension < 2) || (dimension > 4) ) {
    Rcpp::stop("Unsupported dimension type - must be 2,3, or 4");
  }

  if ( precision=="double") {
    using PrecisionType = double;
    if ( dimension == 2 ) {
      using MeshType = itk::Mesh<PrecisionType,2>;
      return antsrMesh_GetNumberOfCells<MeshType>(r_mesh);
    }
    else if ( dimension == 3 ) {
      using MeshType = itk::Mesh<PrecisionType,3>;
      return antsrMesh_GetNumberOfCells<MeshType>(r_mesh);
    }
    else if ( dimension == 4 ) {
      using MeshType = itk::Mesh<PrecisionType,4>;
      return antsrMesh_GetNumberOfCells<MeshType>(r_mesh);
    }
  }
  else if (precision=="float") {
    using PrecisionType = float;
    if ( dimension == 2 ) {
      using MeshType = itk::Mesh<PrecisionType,2>;
      return antsrMesh_GetNumberOfCells<MeshType>(r_mesh);
    }
    else if ( dimension == 3 ) {
      using MeshType = itk::Mesh<PrecisionType,3>;
      return antsrMesh_GetNumberOfCells<MeshType>(r_mesh);
        }
    else if ( dimension == 4 ) {
      using MeshType = itk::Mesh<PrecisionType,4>;
      return antsrMesh_GetNumberOfCells<MeshType>(r_mesh);
    }
  }
  else {
    Rcpp::stop( "Unsupported precision type - must be 'float' or 'double'");
  }

  // Never reached
  return( Rcpp::wrap(NA_REAL) );

}
catch( itk::ExceptionObject & err )
  {
  Rcpp::Rcout << "ITK ExceptionObject caught !" << std::endl;
  Rcpp::Rcout << err << std::endl;
  Rcpp::stop("ITK exception caught");
  }
catch( const std::exception& exc )
  {
  forward_exception_to_r( exc ) ;
  }
catch(...)
  {
	Rcpp::stop("c++ exception (unknown reason)");
  }
return Rcpp::wrap(NA_REAL); //not reached
}

template< class MeshType >
SEXP
antsrMesh_GetPoint( SEXP r_mesh, SEXP r_identifier )
{
  Rcpp::S4 rMesh( r_mesh );
  unsigned int dimension = rMesh.slot("dimension");

  using MeshPointerType = typename MeshType::Pointer;
  MeshPointerType mesh = Rcpp::as<MeshPointerType>( rMesh );
  itk::IdentifierType id = Rcpp::as<itk::IdentifierType>( r_identifier )-1;

  Rcpp::NumericVector pt( dimension );
  typename MeshType::PointType itkPoint;
  mesh->GetPoint(id, &itkPoint);

  for ( unsigned int i=0; i<dimension; i++ ) {
    pt[i] = itkPoint[i];
  }

  return( Rcpp::wrap(pt) );
}

//pixeltype, precision, dimension, type, isVector
RcppExport SEXP antsrMesh_GetPoint( SEXP r_mesh, SEXP r_identifier )
{
try
{
  Rcpp::S4 rMesh( r_mesh );
  std::string precision = rMesh.slot("precision");
  unsigned int dimension = rMesh.slot("dimension");

  if ( (dimension < 2) || (dimension > 4) ) {
    Rcpp::stop("Unsupported dimension type - must be 2,3, or 4");
  }

  if ( precision=="double") {
    using PrecisionType = double;
    if ( dimension == 2 ) {
      using MeshType = itk::Mesh<PrecisionType,2>;
      return antsrMesh_GetPoint<MeshType>(r_mesh, r_identifier);
    }
    else if ( dimension == 3 ) {
      using MeshType = itk::Mesh<PrecisionType,3>;
      return antsrMesh_GetPoint<MeshType>(r_mesh, r_identifier);
    }
    else if ( dimension == 4 ) {
      using MeshType = itk::Mesh<PrecisionType,4>;
      return antsrMesh_GetPoint<MeshType>(r_mesh, r_identifier);
    }
  }
  else if (precision=="float") {
    using PrecisionType = float;
    if ( dimension == 2 ) {
      using MeshType = itk::Mesh<PrecisionType,2>;
      return antsrMesh_GetPoint<MeshType>(r_mesh, r_identifier);
    }
    else if ( dimension == 3 ) {
      using MeshType = itk::Mesh<PrecisionType,3>;
      return antsrMesh_GetPoint<MeshType>(r_mesh, r_identifier);
        }
    else if ( dimension == 4 ) {
      using MeshType = itk::Mesh<PrecisionType,4>;
      return antsrMesh_GetPoint<MeshType>(r_mesh, r_identifier);
    }
  }
  else {
    Rcpp::stop( "Unsupported precision type - must be 'float' or 'double'");
  }

  // Never reached
  return( Rcpp::wrap(NA_REAL) );

}
catch( itk::ExceptionObject & err )
  {
  Rcpp::Rcout << "ITK ExceptionObject caught !" << std::endl;
  Rcpp::Rcout << err << std::endl;
  Rcpp::stop("ITK exception caught");
  }
catch( const std::exception& exc )
  {
  forward_exception_to_r( exc ) ;
  }
catch(...)
  {
	Rcpp::stop("c++ exception (unknown reason)");
  }
return Rcpp::wrap(NA_REAL); //not reached
}

template< class MeshType >
SEXP
antsrMesh_GetCell( SEXP r_mesh, SEXP r_identifier )
{
  Rcpp::S4 rMesh( r_mesh );
  //unsigned int dimension = rMesh.slot("dimension");

  typedef typename MeshType::Pointer            MeshPointerType;
  typedef typename MeshType::CellType           CellType;
  typedef typename CellType::CellAutoPointer    CellAutoPointer;

  MeshPointerType mesh = Rcpp::as<MeshPointerType>( rMesh );
  itk::IdentifierType id = Rcpp::as<itk::IdentifierType>( r_identifier )-1;

  CellAutoPointer cell;
  if ( !mesh->GetCell(id, cell) ) {
    Rcpp::stop("Cell not read");
  }

  Rcpp::NumericVector ids( cell->GetNumberOfPoints() );
  for (unsigned int i=0; i<cell->GetNumberOfPoints(); i++ ) {
    ids[i] = cell->GetPointIds()[i] + 1; // 1-based indexing in R
  }

  return( Rcpp::wrap(ids) );
}

//pixeltype, precision, dimension, type, isVector
RcppExport SEXP antsrMesh_GetCell( SEXP r_mesh, SEXP r_identifier )
{
try
{
  Rcpp::S4 rMesh( r_mesh );
  std::string precision = rMesh.slot("precision");
  unsigned int dimension = rMesh.slot("dimension");

  if ( (dimension < 2) || (dimension > 4) ) {
    Rcpp::stop("Unsupported dimension type - must be 2,3, or 4");
  }

  if ( precision=="double") {
    using PrecisionType = double;
    if ( dimension == 2 ) {
      using MeshType = itk::Mesh<PrecisionType,2>;
      return antsrMesh_GetCell<MeshType>(r_mesh, r_identifier);
    }
    else if ( dimension == 3 ) {
      using MeshType = itk::Mesh<PrecisionType,3>;
      return antsrMesh_GetCell<MeshType>(r_mesh, r_identifier);
    }
    else if ( dimension == 4 ) {
      using MeshType = itk::Mesh<PrecisionType,4>;
      return antsrMesh_GetCell<MeshType>(r_mesh, r_identifier);
    }
  }
  else if (precision=="float") {
    using PrecisionType = float;
    if ( dimension == 2 ) {
      using MeshType = itk::Mesh<PrecisionType,2>;
      return antsrMesh_GetCell<MeshType>(r_mesh, r_identifier);
    }
    else if ( dimension == 3 ) {
      using MeshType = itk::Mesh<PrecisionType,3>;
      return antsrMesh_GetCell<MeshType>(r_mesh, r_identifier);
        }
    else if ( dimension == 4 ) {
      using MeshType = itk::Mesh<PrecisionType,4>;
      return antsrMesh_GetCell<MeshType>(r_mesh, r_identifier);
    }
  }
  else {
    Rcpp::stop( "Unsupported precision type - must be 'float' or 'double'");
  }

  // Never reached
  return( Rcpp::wrap(NA_REAL) );

}
catch( itk::ExceptionObject & err )
  {
  Rcpp::Rcout << "ITK ExceptionObject caught !" << std::endl;
  Rcpp::Rcout << err << std::endl;
  Rcpp::stop("ITK exception caught");
  }
catch( const std::exception& exc )
  {
  forward_exception_to_r( exc ) ;
  }
catch(...)
  {
	Rcpp::stop("c++ exception (unknown reason)");
  }
return Rcpp::wrap(NA_REAL); //not reached
}

template< class MeshType >
SEXP
antsrMesh_GetCellPoints( SEXP r_mesh, SEXP r_identifier )
{
  Rcpp::S4 rMesh( r_mesh );
  unsigned int dimension = rMesh.slot("dimension");

  typedef typename MeshType::Pointer            MeshPointerType;
  typedef typename MeshType::CellType           CellType;
  typedef typename CellType::CellAutoPointer    CellAutoPointer;

  MeshPointerType mesh = Rcpp::as<MeshPointerType>( rMesh );
  itk::IdentifierType id = Rcpp::as<itk::IdentifierType>( r_identifier )-1;

  CellAutoPointer cell;
  if ( !mesh->GetCell(id, cell) ) {
    Rcpp::stop("Cell not read");
  }

  Rcpp::NumericMatrix pts( cell->GetNumberOfPoints(), dimension );
  for (unsigned int i=0; i<cell->GetNumberOfPoints(); i++ ) {
    typename MeshType::PointType pt = mesh->GetPoint( cell->GetPointIds()[i] );
    for (unsigned int j=0; j<dimension; j++) {
      pts(i,j) = pt[j];
    }
  }

  return( Rcpp::wrap(pts) );
}

//pixeltype, precision, dimension, type, isVector
RcppExport SEXP antsrMesh_GetCellPoints( SEXP r_mesh, SEXP r_identifier )
{
try
{
  Rcpp::S4 rMesh( r_mesh );
  std::string precision = rMesh.slot("precision");
  unsigned int dimension = rMesh.slot("dimension");

  if ( (dimension < 2) || (dimension > 4) ) {
    Rcpp::stop("Unsupported dimension type - must be 2,3, or 4");
  }

  if ( precision=="double") {
    using PrecisionType = double;
    if ( dimension == 2 ) {
      using MeshType = itk::Mesh<PrecisionType,2>;
      return antsrMesh_GetCellPoints<MeshType>(r_mesh, r_identifier);
    }
    else if ( dimension == 3 ) {
      using MeshType = itk::Mesh<PrecisionType,3>;
      return antsrMesh_GetCellPoints<MeshType>(r_mesh, r_identifier);
    }
    else if ( dimension == 4 ) {
      using MeshType = itk::Mesh<PrecisionType,4>;
      return antsrMesh_GetCellPoints<MeshType>(r_mesh, r_identifier);
    }
  }
  else if (precision=="float") {
    using PrecisionType = float;
    if ( dimension == 2 ) {
      using MeshType = itk::Mesh<PrecisionType,2>;
      return antsrMesh_GetCellPoints<MeshType>(r_mesh, r_identifier);
    }
    else if ( dimension == 3 ) {
      using MeshType = itk::Mesh<PrecisionType,3>;
      return antsrMesh_GetCellPoints<MeshType>(r_mesh, r_identifier);
        }
    else if ( dimension == 4 ) {
      using MeshType = itk::Mesh<PrecisionType,4>;
      return antsrMesh_GetCellPoints<MeshType>(r_mesh, r_identifier);
    }
  }
  else {
    Rcpp::stop( "Unsupported precision type - must be 'float' or 'double'");
  }

  // Never reached
  return( Rcpp::wrap(NA_REAL) );

}
catch( itk::ExceptionObject & err )
  {
  Rcpp::Rcout << "ITK ExceptionObject caught !" << std::endl;
  Rcpp::Rcout << err << std::endl;
  Rcpp::stop("ITK exception caught");
  }
catch( const std::exception& exc )
  {
  forward_exception_to_r( exc ) ;
  }
catch(...)
  {
	Rcpp::stop("c++ exception (unknown reason)");
  }
return Rcpp::wrap(NA_REAL); //not reached
}

template< class MeshType >
SEXP
antsrMesh_GetPoints( SEXP r_mesh, SEXP r_identifiers )
{
  Rcpp::S4 rMesh( r_mesh );
  unsigned int dimension = rMesh.slot("dimension");

  using MeshPointerType = typename MeshType::Pointer;
  MeshPointerType mesh = Rcpp::as<MeshPointerType>( rMesh );
  Rcpp::NumericVector ids = Rcpp::as<Rcpp::NumericVector>(r_identifiers);

  typename MeshType::PointType itkPoint;
  itk::IdentifierType nPoints = mesh->GetNumberOfPoints();
  if ( ids.size() != 0 ) {
    nPoints = ids.size();
  }
  Rcpp::NumericMatrix pts( nPoints, dimension);

  for ( itk::IdentifierType i=0; i<nPoints; i++ ) {

    itk::IdentifierType id = i;
    if (ids.size() != 0 ) {
      id = ids[i]-1;
    }

    mesh->GetPoint(id, &itkPoint);
    for ( unsigned int j=0; j<dimension; j++) {
      pts(i,j) = itkPoint[j];
    }
  }

  return( Rcpp::wrap(pts) );
}

//pixeltype, precision, dimension, type, isVector
RcppExport SEXP antsrMesh_GetPoints( SEXP r_mesh, SEXP r_identifiers )
{
try
{
  Rcpp::S4 rMesh( r_mesh );
  std::string precision = rMesh.slot("precision");
  unsigned int dimension = rMesh.slot("dimension");

  if ( (dimension < 2) || (dimension > 4) ) {
    Rcpp::stop("Unsupported dimension type - must be 2,3, or 4");
  }

  if ( precision=="double") {
    using PrecisionType = double;
    if ( dimension == 2 ) {
      using MeshType = itk::Mesh<PrecisionType,2>;
      return antsrMesh_GetPoints<MeshType>(r_mesh, r_identifiers);
    }
    else if ( dimension == 3 ) {
      using MeshType = itk::Mesh<PrecisionType,3>;
      return antsrMesh_GetPoints<MeshType>(r_mesh, r_identifiers);
    }
    else if ( dimension == 4 ) {
      using MeshType = itk::Mesh<PrecisionType,4>;
      return antsrMesh_GetPoints<MeshType>(r_mesh, r_identifiers);
    }
  }
  else if (precision=="float") {
    using PrecisionType = float;
    if ( dimension == 2 ) {
      using MeshType = itk::Mesh<PrecisionType,2>;
      return antsrMesh_GetPoints<MeshType>(r_mesh, r_identifiers);
    }
    else if ( dimension == 3 ) {
      using MeshType = itk::Mesh<PrecisionType,3>;
      return antsrMesh_GetPoints<MeshType>(r_mesh, r_identifiers);
        }
    else if ( dimension == 4 ) {
      using MeshType = itk::Mesh<PrecisionType,4>;
      return antsrMesh_GetPoints<MeshType>(r_mesh, r_identifiers);
    }
  }
  else {
    Rcpp::stop( "Unsupported precision type - must be 'float' or 'double'");
  }

  // Never reached
  return( Rcpp::wrap(NA_REAL) );

}
catch( itk::ExceptionObject & err )
  {
  Rcpp::Rcout << "ITK ExceptionObject caught !" << std::endl;
  Rcpp::Rcout << err << std::endl;
  Rcpp::stop("ITK exception caught");
  }
catch( const std::exception& exc )
  {
  forward_exception_to_r( exc ) ;
  }
catch(...)
  {
	Rcpp::stop("c++ exception (unknown reason)");
  }
return Rcpp::wrap(NA_REAL); //not reached
}

template< class MeshType >
SEXP
antsrMesh_AddPoint( SEXP r_mesh, SEXP r_identifier, SEXP r_point )
{
  Rcpp::S4 rMesh( r_mesh );
  unsigned int dimension = rMesh.slot("dimension");

  using MeshPointerType = typename MeshType::Pointer;
  MeshPointerType mesh = Rcpp::as<MeshPointerType>( rMesh );
  itk::IdentifierType id = Rcpp::as<itk::IdentifierType>( r_identifier )-1;

  Rcpp::NumericVector pt = Rcpp::as<Rcpp::NumericVector>( r_point );
  if ( pt.size() != dimension ) {
    Rcpp::stop("Point has incorrect size");
  }

  typename MeshType::PointType itkPoint;

  for ( unsigned int i=0; i<dimension; i++ ) {
    itkPoint[i] = pt(i);
  }

  mesh->GetPoints()->InsertElement(id, itkPoint);

  return( Rcpp::wrap(NA_REAL) );
}

//pixeltype, precision, dimension, type, isVector
RcppExport SEXP antsrMesh_AddPoint( SEXP r_mesh, SEXP r_identifier, SEXP r_point )
{
try
{
  Rcpp::S4 rMesh( r_mesh );
  std::string precision = rMesh.slot("precision");
  unsigned int dimension = rMesh.slot("dimension");

  if ( (dimension < 2) || (dimension > 4) ) {
    Rcpp::stop("Unsupported dimension type - must be 2,3, or 4");
  }

  if ( precision=="double") {
    using PrecisionType = double;
    if ( dimension == 2 ) {
      using MeshType = itk::Mesh<PrecisionType,2>;
      return antsrMesh_AddPoint<MeshType>(r_mesh, r_identifier, r_point);
    }
    else if ( dimension == 3 ) {
      using MeshType = itk::Mesh<PrecisionType,3>;
      return antsrMesh_AddPoint<MeshType>(r_mesh, r_identifier, r_point);
    }
    else if ( dimension == 4 ) {
      using MeshType = itk::Mesh<PrecisionType,4>;
      return antsrMesh_AddPoint<MeshType>(r_mesh, r_identifier, r_point);
    }
  }
  else if (precision=="float") {
    using PrecisionType = float;
    if ( dimension == 2 ) {
      using MeshType = itk::Mesh<PrecisionType,2>;
      return antsrMesh_AddPoint<MeshType>(r_mesh, r_identifier, r_point);
    }
    else if ( dimension == 3 ) {
      using MeshType = itk::Mesh<PrecisionType,3>;
      return antsrMesh_AddPoint<MeshType>(r_mesh, r_identifier, r_point);
        }
    else if ( dimension == 4 ) {
      using MeshType = itk::Mesh<PrecisionType,4>;
      return antsrMesh_AddPoint<MeshType>(r_mesh, r_identifier, r_point);
    }
  }
  else {
    Rcpp::stop( "Unsupported precision type - must be 'float' or 'double'");
  }

  // Never reached
  return( Rcpp::wrap(NA_REAL) );

}
catch( itk::ExceptionObject & err )
  {
  Rcpp::Rcout << "ITK ExceptionObject caught !" << std::endl;
  Rcpp::Rcout << err << std::endl;
  Rcpp::stop("ITK exception caught");
  }
catch( const std::exception& exc )
  {
  forward_exception_to_r( exc ) ;
  }
catch(...)
  {
	Rcpp::stop("c++ exception (unknown reason)");
  }
return Rcpp::wrap(NA_REAL); //not reached
}

template< class MeshType >
SEXP
antsrMesh_SetPoint( SEXP r_mesh, SEXP r_identifier, SEXP r_point )
{
  Rcpp::S4 rMesh( r_mesh );
  unsigned int dimension = rMesh.slot("dimension");

  using MeshPointerType = typename MeshType::Pointer;
  MeshPointerType mesh = Rcpp::as<MeshPointerType>( rMesh );
  itk::IdentifierType id = Rcpp::as<itk::IdentifierType>( r_identifier )-1;

  Rcpp::NumericVector pt = Rcpp::as<Rcpp::NumericVector>( r_point );
  if ( pt.size() != dimension ) {
    Rcpp::stop("Point has incorrect size");
  }

  typename MeshType::PointType itkPoint;

  for ( unsigned int i=0; i<dimension; i++ ) {
    itkPoint[i] = pt(i);
  }

  mesh->SetPoint(id, itkPoint);

  return( Rcpp::wrap(NA_REAL) );
}

//pixeltype, precision, dimension, type, isVector
RcppExport SEXP antsrMesh_SetPoint( SEXP r_mesh, SEXP r_identifier, SEXP r_point )
{
try
{
  Rcpp::S4 rMesh( r_mesh );
  std::string precision = rMesh.slot("precision");
  unsigned int dimension = rMesh.slot("dimension");

  if ( (dimension < 2) || (dimension > 4) ) {
    Rcpp::stop("Unsupported dimension type - must be 2,3, or 4");
  }

  if ( precision=="double") {
    using PrecisionType = double;
    if ( dimension == 2 ) {
      using MeshType = itk::Mesh<PrecisionType,2>;
      return antsrMesh_SetPoint<MeshType>(r_mesh, r_identifier, r_point);
    }
    else if ( dimension == 3 ) {
      using MeshType = itk::Mesh<PrecisionType,3>;
      return antsrMesh_SetPoint<MeshType>(r_mesh, r_identifier, r_point);
    }
    else if ( dimension == 4 ) {
      using MeshType = itk::Mesh<PrecisionType,4>;
      return antsrMesh_SetPoint<MeshType>(r_mesh, r_identifier, r_point);
    }
  }
  else if (precision=="float") {
    using PrecisionType = float;
    if ( dimension == 2 ) {
      using MeshType = itk::Mesh<PrecisionType,2>;
      return antsrMesh_SetPoint<MeshType>(r_mesh, r_identifier, r_point);
    }
    else if ( dimension == 3 ) {
      using MeshType = itk::Mesh<PrecisionType,3>;
      return antsrMesh_SetPoint<MeshType>(r_mesh, r_identifier, r_point);
        }
    else if ( dimension == 4 ) {
      using MeshType = itk::Mesh<PrecisionType,4>;
      return antsrMesh_SetPoint<MeshType>(r_mesh, r_identifier, r_point);
    }
  }
  else {
    Rcpp::stop( "Unsupported precision type - must be 'float' or 'double'");
  }

  // Never reached
  return( Rcpp::wrap(NA_REAL) );

}
catch( itk::ExceptionObject & err )
  {
  Rcpp::Rcout << "ITK ExceptionObject caught !" << std::endl;
  Rcpp::Rcout << err << std::endl;
  Rcpp::stop("ITK exception caught");
  }
catch( const std::exception& exc )
  {
  forward_exception_to_r( exc ) ;
  }
catch(...)
  {
	Rcpp::stop("c++ exception (unknown reason)");
  }
return Rcpp::wrap(NA_REAL); //not reached
}

template< class MeshType >
SEXP
antsrMesh_AddPolyline( SEXP r_mesh, SEXP r_identifier, SEXP r_points )
{
  Rcpp::S4 rMesh( r_mesh );
  //unsigned int dimension = rMesh.slot("dimension");

  using MeshPointerType = typename MeshType::Pointer;
  MeshPointerType mesh = Rcpp::as<MeshPointerType>( rMesh );
  itk::IdentifierType id = Rcpp::as<itk::IdentifierType>( r_identifier )-1;

  typedef typename MeshType::CellType       CellType;
  typedef typename CellType::CellAutoPointer      CellAutoPointer;
  typedef typename itk::PolyLineCell< CellType >  PolyLineCellType;

  Rcpp::NumericVector pts = Rcpp::as<Rcpp::NumericVector>( r_points );
  unsigned int nPoints = pts.size();

  typename MeshType::PointIdentifier polyPoints[ nPoints ];

  for (unsigned long i=0; i<nPoints; i++) {
    polyPoints[i] = pts[i]-1;
  }

  PolyLineCellType * polyline = new PolyLineCellType;
  polyline->SetPointIds( 0, nPoints, polyPoints );
  CellAutoPointer streamline;
  streamline.TakeOwnership( polyline );
  mesh->SetCell(id, streamline);


  return( Rcpp::wrap(NA_REAL) );
}

//pixeltype, precision, dimension, type, isVector
RcppExport SEXP antsrMesh_AddPolyline( SEXP r_mesh, SEXP r_identifier, SEXP r_points )
{
try
{
  Rcpp::S4 rMesh( r_mesh );
  std::string precision = rMesh.slot("precision");
  unsigned int dimension = rMesh.slot("dimension");

  if ( (dimension < 2) || (dimension > 4) ) {
    Rcpp::stop("Unsupported dimension type - must be 2,3, or 4");
  }

  if ( precision=="double") {
    using PrecisionType = double;
    if ( dimension == 2 ) {
      using MeshType = itk::Mesh<PrecisionType,2>;
      return antsrMesh_AddPolyline<MeshType>(r_mesh, r_identifier, r_points);
    }
    else if ( dimension == 3 ) {
      using MeshType = itk::Mesh<PrecisionType,3>;
      return antsrMesh_AddPolyline<MeshType>(r_mesh, r_identifier, r_points);
    }
    else if ( dimension == 4 ) {
      using MeshType = itk::Mesh<PrecisionType,4>;
      return antsrMesh_AddPolyline<MeshType>(r_mesh, r_identifier, r_points);
    }
  }
  else if (precision=="float") {
    using PrecisionType = float;
    if ( dimension == 2 ) {
      using MeshType = itk::Mesh<PrecisionType,2>;
      return antsrMesh_AddPolyline<MeshType>(r_mesh, r_identifier, r_points);
    }
    else if ( dimension == 3 ) {
      using MeshType = itk::Mesh<PrecisionType,3>;
      return antsrMesh_AddPolyline<MeshType>(r_mesh, r_identifier, r_points);
        }
    else if ( dimension == 4 ) {
      using MeshType = itk::Mesh<PrecisionType,4>;
      return antsrMesh_AddPolyline<MeshType>(r_mesh, r_identifier, r_points);
    }
  }
  else {
    Rcpp::stop( "Unsupported precision type - must be 'float' or 'double'");
  }

  // Never reached
  return( Rcpp::wrap(NA_REAL) );

}
catch( itk::ExceptionObject & err )
  {
  Rcpp::Rcout << "ITK ExceptionObject caught !" << std::endl;
  Rcpp::Rcout << err << std::endl;
  Rcpp::stop("ITK exception caught");
  }
catch( const std::exception& exc )
  {
  forward_exception_to_r( exc ) ;
  }
catch(...)
  {
	Rcpp::stop("c++ exception (unknown reason)");
  }
return Rcpp::wrap(NA_REAL); //not reached
}


template< class MeshType >
SEXP
antsrMesh_ReadVTK( SEXP r_filename )
{
  typedef typename MeshType::Pointer      MeshPointerType;
  typedef itk::VtkPolyDataFileReader<MeshType>  ReaderType;

  std::string filename = Rcpp::as<std::string>( r_filename );

  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( filename.c_str() );
  reader->Update();
  MeshPointerType mesh = reader->GetOutput();

  //typename ReaderType::LineSetType * lines = reader->GetLines();

  Rcpp::List pointScalarList( reader->GetPointScalars()->Size() );
  Rcpp::CharacterVector pointScalarNames( reader->GetPointScalars()->Size() );
  for (unsigned int i=0; i<reader->GetPointScalars()->Size(); i++ ) {
   pointScalarList[i] = reader->GetPointScalars()->ElementAt(i);
   pointScalarNames[i] = reader->GetPointScalarsNames()->ElementAt(i);
  }
  pointScalarList.attr("names") = pointScalarNames;

  Rcpp::List cellScalarList( reader->GetCellScalars()->Size() );
  Rcpp::CharacterVector cellScalarNames( reader->GetCellScalars()->Size() );
  for (unsigned int i=0; i<reader->GetCellScalars()->Size(); i++ ) {
   cellScalarList[i] = reader->GetCellScalars()->ElementAt(i);
   cellScalarNames[i] = reader->GetCellScalarsNames()->ElementAt(i);
  }
  cellScalarList.attr("names") = cellScalarNames;

  Rcpp::List list = Rcpp::List::create(Rcpp::Named("Mesh")=Rcpp::wrap(mesh),
                                       Rcpp::Named("PointScalars")=pointScalarList,
                                       Rcpp::Named("CellScalars")=cellScalarList);
  return Rcpp::wrap(list);




}

//pixeltype, precision, dimension, type, isVector
RcppExport SEXP antsrMesh_ReadVTK( SEXP r_filename, SEXP r_dimension, SEXP r_pixeltype )
{
try
{
  std::string precision = Rcpp::as<std::string>(r_pixeltype);
  unsigned int dimension = Rcpp::as<unsigned int>(r_dimension);

  if ( (dimension < 2) || (dimension > 4) ) {
    Rcpp::stop("Unsupported dimension type - must be 2,3, or 4");
  }

  if ( precision=="double") {
    using PrecisionType = double;
    if ( dimension == 2 ) {
      using MeshType = itk::Mesh<PrecisionType,2>;
      return antsrMesh_ReadVTK<MeshType>(r_filename);
    }
    else if ( dimension == 3 ) {
      using MeshType = itk::Mesh<PrecisionType,3>;
      return antsrMesh_ReadVTK<MeshType>(r_filename);
    }
    else if ( dimension == 4 ) {
      using MeshType = itk::Mesh<PrecisionType,4>;
      return antsrMesh_ReadVTK<MeshType>(r_filename);
    }
  }
  else if (precision=="float") {
    using PrecisionType = float;
    if ( dimension == 2 ) {
      using MeshType = itk::Mesh<PrecisionType,2>;
      return antsrMesh_ReadVTK<MeshType>(r_filename);
    }
    else if ( dimension == 3 ) {
      using MeshType = itk::Mesh<PrecisionType,3>;
      return antsrMesh_ReadVTK<MeshType>(r_filename);
      }
    else if ( dimension == 4 ) {
      using MeshType = itk::Mesh<PrecisionType,4>;
      return antsrMesh_ReadVTK<MeshType>(r_filename);
    }
  }
  else {
    Rcpp::stop( "Unsupported precision type - must be 'float' or 'double'");
  }

  // Never reached
  return( Rcpp::wrap(NA_REAL) );

}
catch( itk::ExceptionObject & err )
  {
  Rcpp::Rcout << "ITK ExceptionObject caught !" << std::endl;
  Rcpp::Rcout << err << std::endl;
  Rcpp::stop("ITK exception caught");
  }
catch( const std::exception& exc )
  {
  forward_exception_to_r( exc ) ;
  }
catch(...)
  {
	Rcpp::stop("c++ exception (unknown reason)");
  }
return Rcpp::wrap(NA_REAL); //not reached
}

template< class MeshType >
SEXP
antsrMesh_ReadCamino( SEXP r_filename )
{
  typedef typename MeshType::Pointer                 MeshPointerType;
  typedef itk::CaminoStreamlineFileReader<MeshType>  ReaderType;

  std::string filename = Rcpp::as<std::string>( r_filename );

  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( filename.c_str() );
  reader->Update();
  MeshPointerType mesh = reader->GetOutput();

  Rcpp::NumericVector seeds(reader->GetOutput()->GetNumberOfCells());

  for ( unsigned long i=0; i<reader->GetOutput()->GetNumberOfCells(); i++ ) {
    seeds[i] = static_cast<unsigned long>(reader->GetSeeds()->GetElement(i));
  }

  Rcpp::List list = Rcpp::List::create(Rcpp::Named("Mesh")=Rcpp::wrap(mesh),
                                       Rcpp::Named("Seeds")=seeds);
  return Rcpp::wrap(list);
}

//pixeltype, precision, dimension, type, isVector
RcppExport SEXP antsrMesh_ReadCamino( SEXP r_filename, SEXP r_pixeltype )
{
try
{
  std::string precision = Rcpp::as<std::string>(r_pixeltype);

  if ( precision=="double") {
    using PrecisionType = double;
    using MeshType = itk::Mesh<PrecisionType,3>;
    return antsrMesh_ReadCamino<MeshType>(r_filename);
  }
  else if (precision=="float") {
    using PrecisionType = float;
    using MeshType = itk::Mesh<PrecisionType,3>;
    return antsrMesh_ReadCamino<MeshType>(r_filename);
  }
  else {
    Rcpp::stop( "Unsupported precision type - must be 'float' or 'double'");
  }

  // Never reached
  return( Rcpp::wrap(NA_REAL) );

}
catch( itk::ExceptionObject & err )
  {
  Rcpp::Rcout << "ITK ExceptionObject caught !" << std::endl;
  Rcpp::Rcout << err << std::endl;
  Rcpp::stop("ITK exception caught");
  }
catch( const std::exception& exc )
  {
  forward_exception_to_r( exc ) ;
  }
catch(...)
  {
	Rcpp::stop("c++ exception (unknown reason)");
  }
return Rcpp::wrap(NA_REAL); //not reached
}

template< class MeshType >
SEXP
antsrMesh_ReadTrk( SEXP r_filename )
{
  typedef typename MeshType::Pointer                   MeshPointerType;
  typedef itk::TrackVisStreamlineFileReader<MeshType>  ReaderType;

  std::string filename = Rcpp::as<std::string>( r_filename );

  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( filename.c_str() );
  reader->Update();
  MeshPointerType mesh = reader->GetOutput();

  typename ReaderType::ImagePointerType img = reader->GetReferenceImage();
  Rcpp::List list = Rcpp::List::create(Rcpp::Named("Mesh")=Rcpp::wrap(mesh),
                                       Rcpp::Named("Image")=Rcpp::wrap(img));

  return Rcpp::wrap(list);
}

//pixeltype, precision, dimension, type, isVector
RcppExport SEXP antsrMesh_ReadTrk( SEXP r_filename, SEXP r_pixeltype )
{
try
{
  std::string precision = Rcpp::as<std::string>(r_pixeltype);

  if ( precision=="double") {
    using PrecisionType = double;
    using MeshType = itk::Mesh<PrecisionType,3>;
    return antsrMesh_ReadTrk<MeshType>(r_filename);
  }
  else if (precision=="float") {
    using PrecisionType = float;
    using MeshType = itk::Mesh<PrecisionType,3>;
    return antsrMesh_ReadTrk<MeshType>(r_filename);
  }
  else {
    Rcpp::stop( "Unsupported precision type - must be 'float' or 'double'");
  }

  // Never reached
  return( Rcpp::wrap(NA_REAL) );

}
catch( itk::ExceptionObject & err )
  {
  Rcpp::Rcout << "ITK ExceptionObject caught !" << std::endl;
  Rcpp::Rcout << err << std::endl;
  Rcpp::stop("ITK exception caught");
  }
catch( const std::exception& exc )
  {
  forward_exception_to_r( exc ) ;
  }
catch(...)
  {
	Rcpp::stop("c++ exception (unknown reason)");
  }
return Rcpp::wrap(NA_REAL); //not reached
}

template< class MeshType >
SEXP
antsrMesh_ReadTck( SEXP r_filename )
{
  typedef typename MeshType::Pointer                 MeshPointerType;
  typedef itk::MRTrixStreamlineFileReader<MeshType>  ReaderType;

  std::string filename = Rcpp::as<std::string>( r_filename );

  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( filename.c_str() );
  reader->Update();
  MeshPointerType mesh = reader->GetOutput();

  Rcpp::List list = Rcpp::List::create(Rcpp::Named("Mesh")=Rcpp::wrap(mesh));
  return Rcpp::wrap(list);
}

//pixeltype, precision, dimension, type, isVector
RcppExport SEXP antsrMesh_ReadTck( SEXP r_filename, SEXP r_pixeltype )
{
try
{
  std::string precision = Rcpp::as<std::string>(r_pixeltype);

  if ( precision=="double") {
    using PrecisionType = double;
    using MeshType = itk::Mesh<PrecisionType,3>;
    return antsrMesh_ReadTck<MeshType>(r_filename);
  }
  else if (precision=="float") {
    using PrecisionType = float;
    using MeshType = itk::Mesh<PrecisionType,3>;
    return antsrMesh_ReadTck<MeshType>(r_filename);
  }
  else {
    Rcpp::stop( "Unsupported precision type - must be 'float' or 'double'");
  }

  // Never reached
  return( Rcpp::wrap(NA_REAL) );

}
catch( itk::ExceptionObject & err )
  {
  Rcpp::Rcout << "ITK ExceptionObject caught !" << std::endl;
  Rcpp::Rcout << err << std::endl;
  Rcpp::stop("ITK exception caught");
  }
catch( const std::exception& exc )
  {
  forward_exception_to_r( exc ) ;
  }
catch(...)
  {
	Rcpp::stop("c++ exception (unknown reason)");
  }
return Rcpp::wrap(NA_REAL); //not reached
}

template< class MeshType >
SEXP
antsrMesh_WriteCamino( SEXP r_mesh, SEXP r_filename, SEXP r_seeds )
{
  using MeshPointerType = typename MeshType::Pointer;
  using WriterType = itk::CaminoStreamlineFileWriter< MeshType >;
  using WriterPointerType = typename WriterType::Pointer;

  std::string filename = Rcpp::as<std::string>( r_filename );
  MeshPointerType mesh = Rcpp::as<MeshPointerType>( r_mesh );
  Rcpp::NumericVector seeds( r_seeds );

  WriterPointerType writer = WriterType::New();
  writer->SetFileName( filename );
  writer->SetInput( mesh );
  for (long i=0; i<seeds.size(); i++ )
  {
    writer->SetSeed(i, seeds[i]);
  }
  writer->Update();


  return Rcpp::wrap(1);
}

//pixeltype, precision, dimension, type, isVector
RcppExport SEXP antsrMesh_WriteCamino( SEXP r_mesh, SEXP r_filename, SEXP r_seeds )
{
try
{
  Rcpp::S4 rMesh( r_mesh );
  std::string precision = rMesh.slot("precision");
  //unsigned int dimension = rMesh.slot("dimension");

  if ( precision=="double") {
    using PrecisionType = double;
    using MeshType = itk::Mesh<PrecisionType,3>;
    return antsrMesh_WriteCamino<MeshType>(r_mesh, r_filename, r_seeds);
  }
  else if (precision=="float") {
    using PrecisionType = float;
    using MeshType = itk::Mesh<PrecisionType,3>;
    return antsrMesh_WriteCamino<MeshType>(r_mesh, r_filename, r_seeds);
  }
  else {
    Rcpp::stop( "Unsupported precision type - must be 'float' or 'double'");
  }

  // Never reached
  return( Rcpp::wrap(NA_REAL) );

}
catch( itk::ExceptionObject & err )
  {
  Rcpp::Rcout << "ITK ExceptionObject caught !" << std::endl;
  Rcpp::Rcout << err << std::endl;
  Rcpp::stop("ITK exception caught");
  }
catch( const std::exception& exc )
  {
  forward_exception_to_r( exc ) ;
  }
catch(...)
  {
	Rcpp::stop("c++ exception (unknown reason)");
  }
return Rcpp::wrap(NA_REAL); //not reached
}



template< class MeshType >
SEXP
antsrMesh_WriteVTK( SEXP r_mesh, SEXP r_filename, SEXP r_cellsAs, SEXP r_binary )
{

  //Rcpp::Rcout << "antsrMesh_WriteVTK<MeshType>()" << std::endl;


  using MeshPointerType = typename MeshType::Pointer;
  //using CellType = typename MeshType::CellType;
  //using CellAutoPointer = typename CellType::CellAutoPointer;
  using WriterType = itk::VtkPolyDataFileWriter< MeshType >;
  using WriterPointerType = typename WriterType::Pointer;

  std::string filename = Rcpp::as<std::string>( r_filename );
  MeshPointerType mesh = Rcpp::as<MeshPointerType>( r_mesh );
  std::string cellsAs = Rcpp::as<std::string>( r_cellsAs );
  bool isbinary = Rcpp::as<bool>( r_binary );

  WriterPointerType writer = WriterType::New();
  writer->SetFileName( filename );
  writer->SetInput( mesh );

  if ( cellsAs == "lines" ) {
    writer->SetCellsAsLines(true);
  }
  else if ( cellsAs == "polygons" ) {
    writer->SetCellsAsPolygons(true);
  }

  writer->SetWriteBinary(isbinary);
  writer->Update();

  return Rcpp::wrap(1);
}

//pixeltype, precision, dimension, type, isVector
RcppExport SEXP antsrMesh_WriteVTK( SEXP r_mesh, SEXP r_filename, SEXP r_cellsAs, SEXP r_binary )
{
  //Rcpp::Rcout << "antsrMesh_WriteVTK()" << std::endl;

try
{
  Rcpp::S4 rMesh( r_mesh );
  std::string precision = rMesh.slot("precision");
  //unsigned int dimension = rMesh.slot("dimension");

  if ( precision=="double") {
    using PrecisionType = double;
    using MeshType = itk::Mesh<PrecisionType,3>;
    return antsrMesh_WriteVTK<MeshType>(r_mesh, r_filename, r_cellsAs, r_binary);
  }
  else if (precision=="float") {
    using PrecisionType = float;
    using MeshType = itk::Mesh<PrecisionType,3>;
    return antsrMesh_WriteVTK<MeshType>(r_mesh, r_filename, r_cellsAs, r_binary);
  }
  else {
    Rcpp::stop( "Unsupported precision type - must be 'float' or 'double'");
  }

  // Never reached
  return( Rcpp::wrap(NA_REAL) );

}
catch( itk::ExceptionObject & err )
  {
  Rcpp::Rcout << "ITK ExceptionObject caught !" << std::endl;
  Rcpp::Rcout << err << std::endl;
  Rcpp::stop("ITK exception caught");
  }
catch( const std::exception& exc )
  {
  forward_exception_to_r( exc ) ;
  }
catch(...)
  {
	Rcpp::stop("c++ exception (unknown reason)");
  }
return Rcpp::wrap(NA_REAL); //not reached
}

template< class MeshType, class ImageType >
SEXP
antsrMesh_WriteTrk( SEXP r_mesh, SEXP r_filename, SEXP r_Image)
{
  //Rcpp::Rcout << "antsrMesh_WriteTrk<MeshType>()" << std::endl;

  using MeshPointerType = typename MeshType::Pointer;
  using WriterType = itk::TrackVisStreamlineFileWriter< MeshType, ImageType >;
  using WriterPointerType = typename WriterType::Pointer;
  using ImagePointerType = typename ImageType::Pointer;

  std::string filename = Rcpp::as<std::string>( r_filename );
  MeshPointerType mesh = Rcpp::as<MeshPointerType>( r_mesh );
  ImagePointerType image = Rcpp::as<ImagePointerType>( r_Image );

  WriterPointerType writer = WriterType::New();
  writer->SetFileName( filename );
  writer->SetInput( mesh );
  writer->SetReferenceImage( image );
  writer->Update();

  return Rcpp::wrap(1);
}

//pixeltype, precision, dimension, type, isVector
RcppExport SEXP antsrMesh_WriteTrk( SEXP r_mesh, SEXP r_filename, SEXP r_image)
{
//Rcpp::Rcout << "antsrMesh_WriteTrk()" << std::endl;

try
{
  Rcpp::S4 rMesh( r_mesh );
  std::string precision = rMesh.slot("precision");
  //unsigned int dimension = rMesh.slot("dimension");

  if ( precision=="double") {
    using PrecisionType = double;
    using MeshType = itk::Mesh<PrecisionType,3>;
    using ImageType = itk::Image<PrecisionType,3>;
    return antsrMesh_WriteTrk<MeshType, ImageType>(r_mesh, r_filename, r_image);
  }
  else if (precision=="float") {
    using PrecisionType = float;
    using MeshType = itk::Mesh<PrecisionType,3>;
    using ImageType = itk::Image<PrecisionType,3>;
    return antsrMesh_WriteTrk<MeshType, ImageType>(r_mesh, r_filename, r_image);
  }
  else {
    Rcpp::stop( "Unsupported precision type - must be 'float' or 'double'");
  }

  // Never reached
  return( Rcpp::wrap(NA_REAL) );

}
catch( itk::ExceptionObject & err )
  {
  Rcpp::Rcout << "ITK ExceptionObject caught !" << std::endl;
  Rcpp::Rcout << err << std::endl;
  Rcpp::stop("ITK exception caught");
  }
catch( const std::exception& exc )
  {
  forward_exception_to_r( exc ) ;
  }
catch(...)
  {
	Rcpp::stop("c++ exception (unknown reason)");
  }
return Rcpp::wrap(NA_REAL); //not reached
}

// Apply transform to image
template< class MeshType, class TransformType >
SEXP antsrMesh_TransformMesh( SEXP r_transform, SEXP r_mesh, SEXP r_inplace )
{
  typedef typename MeshType::Pointer                 MeshPointerType;
  typedef typename TransformType::Pointer            TransformPointerType;

  const unsigned int Dimension = TransformType::InputSpaceDimension;

  bool inplace = Rcpp::as<bool>(r_inplace);

  MeshPointerType mesh = Rcpp::as<MeshPointerType>( r_mesh );
  TransformPointerType transform = Rcpp::as<TransformPointerType>( r_transform );

  typedef typename TransformType::ParametersValueType                   PrecisionType;
  typedef typename MeshType::PixelType PixelType;

  typedef typename MeshType::PointType      MeshPointType;

  typename MeshType::Pointer outMesh = nullptr;
  if ( !inplace ) {
    outMesh = MeshType::New();
    outMesh->GetPoints()->Reserve(mesh->GetNumberOfPoints());
  }
  else {
    outMesh = mesh;
  }

  typename TransformType::InputPointType inPoint;
  typename TransformType::OutputPointType outPoint;
  for (unsigned int i=0; i<mesh->GetNumberOfPoints(); i++ ) {
    MeshPointType mPoint = mesh->GetPoint(i);
    for (unsigned int j=0; j<Dimension; j++) {
      inPoint[j] = static_cast<PrecisionType>(mPoint[j]);
    }

    outPoint = transform->TransformPoint( inPoint );
    for (unsigned int j=0; j<Dimension; j++) {
      mPoint[j] = static_cast<PixelType>(outPoint[j]);
    }

    outMesh->SetPoint(i, mPoint);
  }

  return Rcpp::wrap<MeshPointerType>( outMesh );
}

template< class MeshType >
SEXP antsrMesh_TransformMesh( SEXP r_transform, SEXP r_mesh, SEXP r_inplace )
{
  Rcpp::S4 transform( r_transform );
  std::string type = Rcpp::as<std::string>( transform.slot("type") );
  std::string precision = Rcpp::as<std::string>( transform.slot("precision") );

  Rcpp::S4 rMesh( r_mesh );
  std::string pixeltype = rMesh.slot("precision");

  const unsigned int Dimension = MeshType::PointDimension;

  if ( precision == "double" )
  {
    typedef itk::Transform<double, Dimension, Dimension> TransformType;
    return antsrMesh_TransformMesh<MeshType, TransformType>( r_transform, r_mesh, r_inplace );
  }
  else if ( precision == "float" )
  {
    typedef itk::Transform<float, Dimension, Dimension> TransformType;
    return antsrMesh_TransformMesh<MeshType, TransformType>( r_transform, r_mesh, r_inplace );
  }
  else
  {
    Rcpp::stop("Unsupported precision in antsrTransform");
  }

  return Rcpp::wrap(NA_REAL);

}

RcppExport SEXP antsrMesh_TransformMesh( SEXP r_transform, SEXP r_mesh, SEXP r_inplace )
{
try
{
  Rcpp::S4 rMesh( r_mesh );
  std::string pixeltype = rMesh.slot("precision");
  unsigned int idimension = rMesh.slot("dimension");

  Rcpp::S4 transform( r_transform );
  std::string precision = Rcpp::as<std::string>( transform.slot("precision") );
  unsigned int dimension = Rcpp::as<int>( transform.slot("dimension") );

  if ( (dimension < 2) || (dimension > 4) ) {
    Rcpp::stop("Unsupported dimension type - must be 2,3, or 4");
  }
  if ( dimension != idimension ) {
    Rcpp::stop("Image and transform must have same dimension");
  }


  if ( pixeltype=="double") {
    typedef double PixelType;
    if ( dimension == 2 ) {
      typedef itk::Mesh<PixelType,2> MeshType;
      return antsrMesh_TransformMesh<MeshType>(r_transform, r_mesh, r_inplace);
    }
    else if ( dimension == 3 ) {
      typedef itk::Mesh<PixelType,3> MeshType;
      return antsrMesh_TransformMesh<MeshType>(r_transform, r_mesh, r_inplace);
    }
    else if ( dimension == 4 ) {
      typedef itk::Mesh<PixelType,4> MeshType;
      return antsrMesh_TransformMesh<MeshType>(r_transform, r_mesh, r_inplace);
    }
  }
  else if (pixeltype=="float") {
    typedef float PixelType;
    if ( dimension == 2 ) {
      typedef itk::Mesh<PixelType,2> MeshType;
      return antsrMesh_TransformMesh<MeshType>(r_transform, r_mesh, r_inplace);
    }
    else if ( dimension == 3 ) {
      typedef itk::Mesh<PixelType,3> MeshType;
      return antsrMesh_TransformMesh<MeshType>(r_transform, r_mesh, r_inplace);
      }
    else if ( dimension == 4 ) {
      typedef itk::Mesh<PixelType,4> MeshType;
      return antsrMesh_TransformMesh<MeshType>(r_transform, r_mesh, r_inplace);
    }
  }
  else {
    Rcpp::stop( "Unsupported pixel type in mesh - must be 'float' or 'double'");
  }

  // Never reached
  return( Rcpp::wrap(NA_REAL) );

}
catch( itk::ExceptionObject & err )
  {
  Rcpp::Rcout << "ITK ExceptionObject caught !" << std::endl;
  Rcpp::Rcout << err << std::endl;
  Rcpp::stop("ITK exception caught");
  }
catch( const std::exception& exc )
  {
  forward_exception_to_r( exc ) ;
  }
catch(...)
  {
	Rcpp::stop("c++ exception (unknown reason)");
  }
return Rcpp::wrap(NA_REAL); //not reached
}
