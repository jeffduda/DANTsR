#include <string>
#include "isInMesh.h"
#include "itkMesh.h"

template< class MeshType >
SEXP isInMesh( SEXP r_mesh, SEXP r_coordinate, SEXP r_tolerance )
{
  using MeshPointerType = typename MeshType::Pointer;
  MeshPointerType mesh = Rcpp::as<MeshPointerType>(r_mesh);
  if ( ! mesh.IsNotNull() )
    {
    Rcpp::stop("Mesh not yet allocated");
    }
  Rcpp::NumericVector coord( r_coordinate );
  double tolerance = Rcpp::as<double>( r_tolerance );

  typename MeshType::PointType point;
  for (unsigned int i=0; i<MeshType::PointDimension; i++) {
    point[i] = coord[i];
  }

  unsigned long idx = 0;

  for ( unsigned int i=0; i<mesh->GetNumberOfPoints(); i++ ) {
    double distance = point.EuclideanDistanceTo( mesh->GetPoint(i) );
    if ( distance < tolerance ) {
      idx = i+1;
    }
  }

  return Rcpp::wrap( idx );
}

RcppExport SEXP isInMesh( SEXP r_mesh, SEXP r_coordinate, SEXP r_tolerance )
{
try
{
  if( r_mesh == nullptr )
    {
    Rcpp::Rcout << "Unspecified Argument" << std::endl ;
    return Rcpp::wrap( 1 ) ;
    }

  Rcpp::S4 antsrmesh( r_mesh );
  std::string pixeltype = Rcpp::as< std::string >( antsrmesh.slot( "precision" ) );
  unsigned int dimension = Rcpp::as< int >( antsrmesh.slot( "dimension" ) );

  if( pixeltype == "double" )
    {
    typedef double PixelType;
    if ( dimension == 2 )
      {
      const unsigned int Dimension = 2;
      typedef itk::Mesh< PixelType, Dimension >   MeshType;
      return  isInMesh<MeshType>( r_mesh, r_coordinate, r_tolerance );
      }
    else if ( dimension == 3 )
      {
      const unsigned int Dimension = 3;
      typedef itk::Mesh< PixelType, Dimension >   MeshType;
      return  isInMesh<MeshType>( r_mesh, r_coordinate, r_tolerance );
      }
    else if ( dimension == 4 )
      {
      const unsigned int Dimension = 4;
      typedef itk::Mesh< PixelType, Dimension >   MeshType;
      return  isInMesh<MeshType>( r_mesh, r_coordinate, r_tolerance );
      }
    }
  else if( pixeltype == "float" )
    {
    typedef float PixelType;
    if ( dimension == 2 )
      {
      const unsigned int Dimension = 2;
      typedef itk::Mesh< PixelType, Dimension >   MeshType;
      return  isInMesh<MeshType>( r_mesh, r_coordinate, r_tolerance );
      }
    else if ( dimension == 3 )
      {
      const unsigned int Dimension = 3;
      typedef itk::Mesh< PixelType, Dimension >   MeshType;
      return  isInMesh<MeshType>( r_mesh, r_coordinate, r_tolerance );
      }
    else if ( dimension == 4 )
      {
      const unsigned int Dimension = 4;
      typedef itk::Mesh< PixelType, Dimension >   MeshType;
      return  isInMesh<MeshType>( r_mesh, r_coordinate, r_tolerance );
      }
    }
  else
    {
    Rcpp::Rcout << "Unsupported antsrMesh@Precision" << std::endl ;
    return Rcpp::wrap( NA_REAL ) ;
    }

  return Rcpp::wrap(NA_REAL);

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
