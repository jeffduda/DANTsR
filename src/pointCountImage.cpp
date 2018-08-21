#include <RcppDANTsR.h>
#include <string>
#include "itkImage.h"
#include "itkPointCountImageFilter.h"

//#include "pointCountImage.h"

template< class MeshType, class ImageType >
SEXP pointCountImageFunction( SEXP r_mesh, SEXP r_image )
{
  //Rcpp::Rcout << "pointCountImage<MeshType,ImageType>()" << std::endl;

  using ImagePointerType = typename ImageType::Pointer;
  using MeshPointerType = typename MeshType::Pointer;
  using PixelType = typename MeshType::PixelType;
  using OutputImageType = itk::Image< PixelType, ImageType::ImageDimension >;
  using OutputImagePointer = typename OutputImageType::Pointer;

  using FilterType = itk::PointCountImageFilter<MeshType, OutputImageType>;
  using FilterPointerType = typename FilterType::Pointer;


  ImagePointerType image = Rcpp::as<ImagePointerType>(r_image);
  if ( ! image.IsNotNull() )
    {
    Rcpp::stop("Image not yet allocated");
    }

  MeshPointerType mesh = Rcpp::as<MeshPointerType>(r_mesh);
  if ( ! image.IsNotNull() )
    {
    Rcpp::stop("Mesh is not allocated");
    }

  FilterPointerType filter = FilterType::New();
  filter->SetInput( mesh );
  filter->SetInfoImage( image );
  filter->Update();

  OutputImagePointer countImage = filter->GetOutput();
  return Rcpp::wrap( countImage );
}

RcppExport SEXP pointCountImageCall( SEXP r_mesh, SEXP r_image )
{
try
{
  //Rcpp::Rcout << "pointCountImage()" << std::endl;

  if( (r_image == nullptr) || (r_mesh == nullptr) )
    {
    Rcpp::Rcout << "Unspecified Argument" << std::endl ;
    return Rcpp::wrap( 1 ) ;
    }

  Rcpp::S4 antsrmesh( r_mesh );
  std::string pixeltype = Rcpp::as< std::string >( antsrmesh.slot( "precision" ) );
  unsigned int dimension = Rcpp::as< int >( antsrmesh.slot( "dimension" ) );

  Rcpp::S4 antsimage( r_image ) ;
  unsigned int idimension = Rcpp::as< int >( antsimage.slot( "dimension" ) );

  if ( dimension != idimension ) {
    Rcpp::Rcout << "Mesh and image must have same dimension" << std::endl;
    return Rcpp::wrap(NA_REAL);
  }

  if ( ( dimension > 4) || (dimension < 2) ) {
    Rcpp::Rcout << "Invalid Dimension: must be 2,3, or 4" << std::endl;
    return Rcpp::wrap(NA_REAL);
  }

  if( pixeltype == "double" )
    {
    typedef double PixelType;
    if ( dimension == 2 )
      {
      const unsigned int Dimension = 2;
      typedef itk::ImageBase< Dimension >  ImageType;
      typedef itk::Mesh< PixelType, Dimension >   MeshType;
      return  pointCountImageFunction< MeshType, ImageType >( r_mesh, r_image );
      }
    else if ( dimension == 3 )
      {
      const unsigned int Dimension = 3;
      typedef itk::ImageBase< Dimension >  ImageType;
      typedef itk::Mesh< PixelType, Dimension >   MeshType;
      return  pointCountImageFunction< MeshType, ImageType >( r_mesh, r_image );
      }
    else if ( dimension == 4 )
      {
      const unsigned int Dimension = 4;
      typedef itk::ImageBase< Dimension >  ImageType;
      typedef itk::Mesh< PixelType, Dimension >   MeshType;
      return  pointCountImageFunction< MeshType, ImageType >( r_mesh, r_image );
      }
    }
  else if( pixeltype == "float" )
    {
    typedef float PixelType;
    if ( dimension == 2 )
      {
      const unsigned int Dimension = 2;
      typedef itk::Image< PixelType, Dimension >  ImageType;
      typedef itk::Mesh< PixelType, Dimension >   MeshType;
      return  pointCountImageFunction< MeshType, ImageType >( r_mesh, r_image );
      }
    else if ( dimension == 3 )
      {
      const unsigned int Dimension = 3;
      typedef itk::Image< PixelType, Dimension >  ImageType;
      typedef itk::Mesh< PixelType, Dimension >   MeshType;
      return  pointCountImageFunction< MeshType, ImageType >( r_mesh, r_image );
      }
    else if ( dimension == 4 )
      {
      const unsigned int Dimension = 4;
      typedef itk::Image< PixelType, Dimension >  ImageType;
      typedef itk::Mesh< PixelType, Dimension >   MeshType;
      return  pointCountImageFunction< MeshType, ImageType >( r_mesh, r_image );
      }
    }
  else
    {
    Rcpp::Rcout << "Unsupported antsrMesh@Precision" << std::endl ;
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
