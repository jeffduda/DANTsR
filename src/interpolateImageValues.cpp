
#include <algorithm>
#include <vector>
#include <string>
#include <RcppDANTsR.h>

#include "itkImage.h"
#include "itkInterpolateImageFunction.h"

template< class ImageType >
SEXP interpolateImageValues( SEXP r_img, SEXP r_points, SEXP r_type, SEXP r_interp )
{

  //Rcpp::Rcout << "interpolateImageValues<ImageType>" << std::endl;

  typedef typename ImageType::PixelType PixelType;
  typedef typename ImageType::Pointer ImagePointerType;

  bool isPoint = true;

  ImagePointerType image = Rcpp::as<ImagePointerType>(r_img);

  Rcpp::NumericMatrix points( r_points );
  Rcpp::NumericVector values( points.nrow() );

  typedef itk::InterpolateImageFunction<ImageType> FunctionType;
  typedef typename FunctionType::Pointer           FunctionPointerType;

  typedef itk::LinearInterpolateImageFunction<ImageType> LinearInterpolatorType;

  FunctionPointerType function = dynamic_cast<FunctionType *>(LinearInterpolatorType::New().GetPointer());

  function->SetInputImage( image );
  typename ImageType::PointType pt;
  typename FunctionType::ContinuousIndexType idx;

  for ( unsigned long i=0; i<points.nrow(); i++ ) {
    for (unsigned int j=0; j<ImageType::ImageDimension; j++ ) {
      if ( isPoint )
        pt[j] = points(i,j);
      else
        idx[j] = points(i,j);
    }

    if ( isPoint ) {
      values[i] = function->Evaluate( pt );
      //Rcpp::Rcout << pt << " - " << values[i] << std::endl;
    }
    else {
      values[i] = function->EvaluateAtContinuousIndex( idx );
      //Rcpp::Rcout << idx << " - " << values[i] << std::endl;
    }
  }

  return( Rcpp::wrap(values) );

}



RcppExport SEXP interpolateImageValues( SEXP r_img, SEXP r_points, SEXP r_type, SEXP r_interp ) {
try
  {
  //Rcpp::Rcout << "interpolateImageValues()" << std::endl;
  if( r_img == NULL || r_points == NULL )
    {
      Rcpp::Rcout << "Unspecified Arguments" << std::endl;
      return Rcpp::wrap( NA_REAL );
    }

  Rcpp::S4 antsimage( r_img ) ;
  std::string pixeltype = Rcpp::as< std::string >( antsimage.slot( "pixeltype" ) );
  unsigned int dimension = Rcpp::as< int >( antsimage.slot( "dimension" ) );
  unsigned int components = Rcpp::as< int >( antsimage.slot( "components" ) );
  //bool isVector = Rcpp::as<bool>( antsimage.slot("isVector") );

  Rcpp::NumericMatrix mat( r_points );
  //Rcpp::Rcout << "Matrix is of size " << mat.nrow() << " x " << mat.ncol() << std::endl;

  if ( components != 1 ) {
    Rcpp::Rcout << "Input must have 1 components" << std::endl;
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
      typedef itk::Image< PixelType, Dimension >  ImageType;
      return interpolateImageValues< ImageType >( r_img, r_points, r_type, r_interp );
      }
    else if ( dimension == 3 )
      {
      const unsigned int Dimension = 3;
      typedef itk::Image< PixelType, Dimension >  ImageType;
      return interpolateImageValues< ImageType >( r_img, r_points, r_type, r_interp );
      }
    else if ( dimension == 4 )
      {
      const unsigned int Dimension = 4;
      typedef itk::Image< PixelType, Dimension >  ImageType;
      return interpolateImageValues< ImageType >( r_img, r_points, r_type, r_interp );
      }
    }
  else if( pixeltype == "float" )
    {
    typedef float PixelType;
    if ( dimension == 2 )
      {
      const unsigned int Dimension = 2;
      typedef itk::Image< PixelType, Dimension >  ImageType;
      return interpolateImageValues< ImageType >( r_img, r_points, r_type, r_interp );
      }
    else if ( dimension == 3 )
      {
      const unsigned int Dimension = 3;
      typedef itk::Image< PixelType, Dimension >  ImageType;
      return interpolateImageValues< ImageType >( r_img, r_points, r_type, r_interp );
      }
    else if ( dimension == 4 )
      {
      const unsigned int Dimension = 4;
      typedef itk::Image< PixelType, Dimension >  ImageType;
      return interpolateImageValues< ImageType >( r_img, r_points, r_type, r_interp );
      }
    }
  else if( pixeltype == "unsigned int" )
    {
    typedef unsigned int PixelType;
    if ( dimension == 2 )
      {
      const unsigned int Dimension = 2;
      typedef itk::Image< PixelType, Dimension >  ImageType;
      return interpolateImageValues< ImageType >( r_img, r_points, r_type, r_interp );
      }
    else if ( dimension == 3 )
      {
      const unsigned int Dimension = 3;
      typedef itk::Image< PixelType, Dimension >  ImageType;
      return interpolateImageValues< ImageType >( r_img, r_points, r_type, r_interp );
      }
    else if ( dimension == 4 )
      {
      const unsigned int Dimension = 4;
      typedef itk::Image< PixelType, Dimension >  ImageType;
      return interpolateImageValues< ImageType >( r_img, r_points, r_type, r_interp );
      }
    }
  else if( pixeltype == "unsigned char" )
    {
    typedef unsigned char PixelType;
    if ( dimension == 2 )
      {
      const unsigned int Dimension = 2;
      typedef itk::Image< PixelType, Dimension >  ImageType;
      return interpolateImageValues< ImageType >( r_img, r_points, r_type, r_interp );
      }
    else if ( dimension == 3 )
      {
      const unsigned int Dimension = 3;
      typedef itk::Image< PixelType, Dimension >  ImageType;
      return interpolateImageValues< ImageType >( r_img, r_points, r_type, r_interp );
      }
    else if ( dimension == 4 )
      {
      const unsigned int Dimension = 4;
      typedef itk::Image< PixelType, Dimension >  ImageType;
      return interpolateImageValues< ImageType >( r_img, r_points, r_type, r_interp );
      }
    }
  else
    {
    Rcpp::Rcout << "Unsupported PixelType" << std::endl ;
    return Rcpp::wrap( NA_REAL ) ;
    }

  }
catch( const std::exception& exc )
  {
    Rcpp::Rcout<< exc.what() << std::endl;
    return Rcpp::wrap( NA_REAL );
  }

// Never reached
return Rcpp::wrap( NA_REAL );

}
