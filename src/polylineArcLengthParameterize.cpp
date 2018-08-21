#include <string>
#include "itkSpatialPolyLineMeshToArcLengthMeshFilter.h"
#include "polylineArcLengthParameterize.h"


template< class InputMeshType, class OutputMeshType >
SEXP polylineArcLengthParameterize( SEXP r_mesh )
{
  Rcpp::Rcout << "polylineArcLengthParameterize<MeshType,MeshType>()" << std::endl;

  using InputMeshPointerType = typename InputMeshType::Pointer;
  using InputPixelType = typename InputMeshType::PixelType;

  using OutputMeshPointerType = typename OutputMeshType::Pointer;

  using FilterType = itk::SpatialPolyLineMeshToArcLengthMeshFilter<InputMeshType, OutputMeshType>;
  using FilterPointerType = typename FilterType::Pointer;

  InputMeshPointerType mesh = Rcpp::as<InputMeshPointerType>(r_mesh);
  if ( ! mesh.IsNotNull() )
    {
    Rcpp::stop("Mesh is not allocated");
    }

  FilterPointerType filter = FilterType::New();
  filter->SetInput( mesh );
  filter->Update();

  OutputMeshPointerType arcMesh = filter->GetOutput();
  //return Rcpp::wrap( 1 );
  return Rcpp::wrap( arcMesh );
}

RcppExport SEXP polylineArcLengthParameterize( SEXP r_mesh )
{
try
{
  Rcpp::Rcout << "polylineArcLengthParameterize()" << std::endl;

  if (r_mesh == nullptr)
    {
    Rcpp::Rcout << "Unspecified Argument" << std::endl ;
    return Rcpp::wrap( 1 ) ;
    }

  Rcpp::S4 antsrmesh( r_mesh );
  std::string pixeltype = Rcpp::as< std::string >( antsrmesh.slot( "precision" ) );
  unsigned int dimension = Rcpp::as< int >( antsrmesh.slot( "dimension" ) );

  if ( ( dimension > 4) || (dimension < 1) ) {
    Rcpp::Rcout << "Invalid Dimension: must be 1,2,3, or 4" << std::endl;
    return Rcpp::wrap(NA_REAL);
  }

  if( pixeltype == "double" )
    {
    using InputPixelType = double;
    if ( dimension == 1 )
      {
      const unsigned int Dimension = 1;
      using InputMeshType = itk::Mesh< InputPixelType, Dimension >;
      using OutputPixelType = typename InputMeshType::PointType;
      using OutputMeshType = InputMeshType;
      return  polylineArcLengthParameterize< InputMeshType, OutputMeshType >( r_mesh);
      }
    else if ( dimension == 2 )
      {
      const unsigned int Dimension = 2;
      using InputMeshType = itk::Mesh< InputPixelType, Dimension >;
      using OutputPixelType = typename InputMeshType::PointType;
      using OutputMeshType = InputMeshType;
      return  polylineArcLengthParameterize< InputMeshType, OutputMeshType >( r_mesh);
      }
    else if ( dimension == 3 )
      {
      const unsigned int Dimension = 3;
      using InputMeshType = itk::Mesh< InputPixelType, Dimension >;
      using OutputPixelType = typename InputMeshType::PointType;
      using OutputMeshType = InputMeshType;
      return  polylineArcLengthParameterize< InputMeshType, OutputMeshType >( r_mesh);
      }
    else if ( dimension == 4 )
      {
      const unsigned int Dimension = 4;
      using InputMeshType = itk::Mesh< InputPixelType, Dimension >;
      using OutputPixelType = typename InputMeshType::PointType;
      using OutputMeshType = InputMeshType;
      return  polylineArcLengthParameterize< InputMeshType, OutputMeshType >( r_mesh);
      }
    }
  else if( pixeltype == "float" )
    {
    using InputPixelType = float;
    if ( dimension == 1 )
      {
      const unsigned int Dimension = 1;
      using InputMeshType = itk::Mesh< InputPixelType, Dimension >;
      using OutputPixelType = typename InputMeshType::PointType;
      using OutputMeshType = InputMeshType;
      return  polylineArcLengthParameterize< InputMeshType, OutputMeshType >( r_mesh);
      }
    else if ( dimension == 2 )
      {
      const unsigned int Dimension = 2;
      using InputMeshType = itk::Mesh< InputPixelType, Dimension >;
      using OutputPixelType = typename InputMeshType::PointType;
      using OutputMeshType = InputMeshType;
      return  polylineArcLengthParameterize< InputMeshType, OutputMeshType >( r_mesh);
      }
    else if ( dimension == 3 )
      {
      const unsigned int Dimension = 3;
      using InputMeshType = itk::Mesh< InputPixelType, Dimension >;
      using OutputPixelType = typename InputMeshType::PointType;
      using OutputMeshType = InputMeshType;
      return  polylineArcLengthParameterize< InputMeshType, OutputMeshType >( r_mesh);
      }
    else if ( dimension == 4 )
      {
      const unsigned int Dimension = 4;
      using InputMeshType = itk::Mesh< InputPixelType, Dimension >;
      using OutputPixelType = typename InputMeshType::PointType;
      using OutputMeshType = InputMeshType;
      return  polylineArcLengthParameterize< InputMeshType, OutputMeshType >( r_mesh);
      }
    }
  else
    {
    Rcpp::Rcout << "Unsupported antsrMesh@precision" << std::endl ;
    return Rcpp::wrap( NA_REAL ) ;
    }

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
