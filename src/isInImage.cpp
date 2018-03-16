
#include <algorithm>
#include <vector>
#include <string>
#include <RcppDANTsR.h>

#include "itkImage.h"

template< class ImageType >
SEXP isInImage( SEXP r_image, SEXP r_coordinate, SEXP r_type )
{
  typedef typename ImageType::Pointer ImagePointerType;
  ImagePointerType image = Rcpp::as<ImagePointerType>(r_image);
  if ( ! image.IsNotNull() )
    {
    Rcpp::stop("Image not yet allocated");
    }
  Rcpp::NumericVector coord( r_coordinate );
  std::string type = Rcpp::as<std::string>( r_type );

  typename itk::ContinuousIndex<float, ImageType::ImageDimension> index;

  if ( type=="point" ) {
    typename ImageType::PointType point;
    for (unsigned int i=0; i<ImageType::ImageDimension; i++) {
      point[i] = coord[i];
    }
    image->TransformPhysicalPointToContinuousIndex(point, index);
  }
  else if ( type=="index" ) {
    for (unsigned int i=0; i<ImageType::ImageDimension; i++) {
      index[i] = coord[i];
    }
  }
  else {
    Rcpp::stop("Only point and index types are allowed");
  }

  return Rcpp::wrap( image->GetLargestPossibleRegion().IsInside(index));
}

RcppExport SEXP isInImage( SEXP r_image, SEXP r_coordinate, SEXP r_type )
{
try
{
  if( r_image == NULL )
    {
    Rcpp::Rcout << "Unspecified Argument" << std::endl ;
    return Rcpp::wrap( 1 ) ;
    }

  Rcpp::S4 image( r_image );
  unsigned int dimension = Rcpp::as< int >( image.slot( "dimension" ) );

  if (dimension == 2)
    {
    typedef itk::ImageBase<2>      ImageType;
    return isInImage<ImageType>( r_image, r_coordinate, r_type );
    }
  else if (dimension == 3)
    {
    typedef itk::ImageBase<3>      ImageType;
    return isInImage<ImageType>( r_image, r_coordinate, r_type );
    }
  else if (dimension == 4)
    {
    typedef itk::ImageBase<4>      ImageType;
    return isInImage<ImageType>( r_image, r_coordinate, r_type );
    }
  else
    {
    Rcpp::stop( "Invalid image dimension");
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
