#include <RcppDANTsR.h>
#include <string>
#include "itkImage.h"
#include "itkCellCountImageFilter.h"

template< class MeshType, class ImageType >
SEXP cellCountImageFunction( SEXP r_mesh, SEXP r_image, SEXP r_target, SEXP r_subset )
{

  using ImagePointerType = typename ImageType::Pointer;
  using MeshPointerType = typename MeshType::Pointer;
  using PixelType = typename MeshType::PixelType;
  using OutputImageType = itk::Image< PixelType, ImageType::ImageDimension >;
  using OutputImagePointer = typename OutputImageType::Pointer;

  using FilterType = itk::CellCountImageFilter<MeshType, OutputImageType>;
  using FilterPointerType = typename FilterType::Pointer;

  ImagePointerType image = Rcpp::as<ImagePointerType>(r_image);
  if ( ! image.IsNotNull() )
    {
    Rcpp::stop("Image not yet allocated");
    }

  MeshPointerType mesh = Rcpp::as<MeshPointerType>(r_mesh);
  if ( ! mesh.IsNotNull() )
    {
    Rcpp::stop("Mesh is not allocated");
    }

  Rcpp::NumericVector subset( r_subset );
  std::vector<double> stdsubset = Rcpp::as< std::vector<double> >( subset );

  bool target = Rcpp::as<bool>(r_target);

  FilterPointerType filter = FilterType::New();
  filter->SetInput( mesh );
  filter->SetInfoImage( image );
  filter->SetTarget( target );
  filter->SetSubset( stdsubset );
  //if ( r_target != nullptr ) {
  //  ImagePointerType target = Rcpp::as<ImagePointerType>(r_target);
  //  filter->SetTarget( target );
  //}
  filter->Update();

  OutputImagePointer countImage = filter->GetOutput();
  countImage->DisconnectPipeline();

  return Rcpp::wrap( countImage );
}

RcppExport SEXP cellCountImage( SEXP r_mesh, SEXP r_image, SEXP r_target, SEXP r_subset )
{
try
{
  if( (r_image == nullptr) || (r_mesh == nullptr) )
    {
    Rcpp::stop("Unspecified Argument");
    }

  Rcpp::S4 antsrmesh( r_mesh );
  std::string pixeltype = Rcpp::antsrMeshPrecision(r_mesh);
  unsigned int dimension = Rcpp::antsrMeshDimension(r_mesh);

  Rcpp::S4 antsimage( r_image ) ;
  unsigned int idimension = Rcpp::as< int >( antsimage.slot( "dimension" ) );

  if ( dimension != idimension ) {
    Rcpp::stop("Mesh and image must have same dimension");
  }

  if ( ( dimension > 4) || (dimension < 2) ) {
    Rcpp::stop("Invalid Dimension: must be 2,3, or 4");
  }

  if( pixeltype == "double" )
    {
    typedef double PixelType;
    if ( dimension == 2 )
      {
      const unsigned int Dimension = 2;
      typedef itk::Image< PixelType,Dimension >  ImageType;
      typedef itk::Mesh< PixelType, Dimension >   MeshType;
      return  cellCountImageFunction< MeshType, ImageType >( r_mesh, r_image, r_target, r_subset );
      }
    else if ( dimension == 3 )
      {
      const unsigned int Dimension = 3;
      typedef itk::Image< PixelType,Dimension >  ImageType;
      typedef itk::Mesh< PixelType, Dimension >   MeshType;
      return  cellCountImageFunction< MeshType, ImageType >( r_mesh, r_image, r_target, r_subset  );
      }
    else if ( dimension == 4 )
      {
      const unsigned int Dimension = 4;
      typedef itk::Image< PixelType,Dimension >  ImageType;
      typedef itk::Mesh< PixelType, Dimension >   MeshType;
      return  cellCountImageFunction< MeshType, ImageType >( r_mesh, r_image, r_target, r_subset  );
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
      return  cellCountImageFunction< MeshType, ImageType >( r_mesh, r_image, r_target, r_subset  );
      }
    else if ( dimension == 3 )
      {
      const unsigned int Dimension = 3;
      typedef itk::Image< PixelType, Dimension >  ImageType;
      typedef itk::Mesh< PixelType, Dimension >   MeshType;
      return  cellCountImageFunction< MeshType, ImageType >( r_mesh, r_image, r_target, r_subset  );
      }
    else if ( dimension == 4 )
      {
      const unsigned int Dimension = 4;
      typedef itk::Image< PixelType, Dimension >  ImageType;
      typedef itk::Mesh< PixelType, Dimension >   MeshType;
      return  cellCountImageFunction< MeshType, ImageType >( r_mesh, r_image, r_target, r_subset  );
      }
    }
  else
    {
    Rcpp::stop("Unsupported antsrMesh@Precision");
    return Rcpp::wrap( NA_REAL ) ;
    }

  // never reached
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
