#include <RcppDANTsR.h>
#include <string>
#include "itkImage.h"
#include "itkOrientImageFilter.h"

itk::SpatialOrientation::ValidCoordinateOrientationFlags getOrientationCode( std::string code )
{
  if ( code=="RIP" )  {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RIP;
  }
  else if (code=="LIP") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LIP;
  }
  else if (code=="RSP") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RSP;
  }
  else if (code=="LSP") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LSP;
  }
  else if (code=="RIA") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RIA;
  }
  else if (code=="LIA") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LIA;
  }
  else if (code=="RSA") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RSA;
  }
  else if (code=="LSA") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LSA;
  }
  else if (code=="IRP") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IRP;
  }
  else if (code=="ILP") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ILP;
  }
  else if (code=="SRP") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SRP;
  }
  else if (code=="SLP") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SLP;
  }
  else if (code=="IRA") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IRA;
  }
  else if (code=="ILA") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ILA;
  }
  else if (code=="SRA") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SRA;
  }
  else if (code=="SLA") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SLA;
  }
  else if (code=="RPI") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RPI;
  }
  else if (code=="LPI") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LPI;
  }
  else if (code=="RAI") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI;
  }
  else if (code=="LAI") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LAI;
  }
  else if (code=="RPS") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RPS;
  }
  else if (code=="LPS") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LPS;
  }
  else if (code=="RAS") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAS;
  }
  else if (code=="LAS") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LAS;
  }
  else if (code=="PRI") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PRI;
  }
  else if (code=="PLI") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PLI;
  }
  else if (code=="ARI") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ARI;
  }
  else if (code=="ALI") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ALI;
  }
  else if (code=="PRS") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PRS;
  }
  else if (code=="PLS") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PLS;
  }
  else if (code=="ARS") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ARS;
  }
  else if (code=="ALS") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ALS;
  }
  else if (code=="IPR") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IPR;
  }
  else if (code=="SPR") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SPR;
  }
  else if (code=="IAR") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IAR;
  }
  else if (code=="SAR") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SAR;
  }
  else if (code=="IPL") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IPL;
  }
  else if (code=="SPL") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SPL;
  }
  else if (code=="IAL") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IAL;
  }
  else if (code=="SAL") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SAL;
  }
  else if (code=="PIR") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PIR;
  }
  else if (code=="PSR") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PSR;
  }
  else if (code=="AIR") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_AIR;
  }
  else if (code=="ASR") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ASR;
  }
  else if (code=="PIL") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PIL;
  }
  else if (code=="PSL") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PSL;
  }
  else if (code=="AIL") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_AIL;
  }
  else if (code=="ASL") {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ASL;
  }
  else {
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_INVALID;
  }

}

template< class ImageType >
SEXP orientImageFunction( SEXP r_image, SEXP r_given, SEXP r_desired )
{

  using ImagePointerType = typename ImageType::Pointer;

  using FilterType = itk::OrientImageFilter<ImageType, ImageType>;
  using FilterPointerType = typename FilterType::Pointer;

  ImagePointerType image = Rcpp::as<ImagePointerType>(r_image);
  if ( ! image.IsNotNull() )
    {
    Rcpp::stop("Image not yet allocated");
    }

  std::string given = Rcpp::as<std::string>(r_given);
  std::string desired = Rcpp::as<std::string>(r_desired);

  std::cout << given << " -> " << desired << std::endl;

  FilterPointerType filter = FilterType::New();
  filter->SetInput( image );
  if ( given=="DIR" ) {
    filter->UseImageDirectionOn();
  }
  else {
    filter->SetGivenCoordinateOrientation( getOrientationCode(given) );
  }

  filter->SetDesiredCoordinateOrientation( getOrientationCode(desired) );
  filter->Update();

  ImagePointerType orientedImage = filter->GetOutput();
  orientedImage->DisconnectPipeline();

  return Rcpp::wrap( orientedImage );
}

RcppExport SEXP orientImage( SEXP r_image, SEXP r_given, SEXP r_desired )
{
try
{

  if( r_image == nullptr )
    {
    Rcpp::Rcout << "Unspecified Argument" << std::endl ;
    return Rcpp::wrap( 1 ) ;
    }

  Rcpp::S4 antsimage( r_image ) ;
  unsigned int dimension = Rcpp::as< int >( antsimage.slot( "dimension" ) );
  std::string pixeltype = Rcpp::as<std::string>( antsimage.slot("pixeltype"));

  if ( dimension != 3  ) {
    Rcpp::Rcout << "Invalid Dimension: must 3" << std::endl;
    return Rcpp::wrap(NA_REAL);
  }

  //typedef itk::Image<float, 3> ImageType;
  //return orientImageFunction< ImageType >( r_image, r_given, r_desired );

  const unsigned int Dimension = 3;
  if( pixeltype == "double" )
    {
    typedef double PixelType;
    typedef itk::Image< PixelType, Dimension >  ImageType;
    return  orientImageFunction< ImageType >( r_image, r_given, r_desired );
    }
  else if( pixeltype == "float" )
    {
    typedef float PixelType;
    typedef itk::Image< PixelType, Dimension >  ImageType;
    return  orientImageFunction< ImageType >( r_image, r_given, r_desired );
    }
  else
    {
    Rcpp::Rcout << "Unsupported antsImage@pixeltype" << std::endl ;
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
