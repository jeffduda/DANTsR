#ifndef __RCPPANTSR_IMAGE_HPP
#define __RCPPANTSR_IMAGE_HPP

#include "itkMacro.h"
#include "itkImage.h"
#include "itkCastImageFilter.h"
#include "itkDiffusionTensor3D.h"
#include <Rcpp.h>

namespace Rcpp {

/* Example code for implementing 'wrap' for an itk image class
template <> inline
SEXP wrap( const itk::IMAGETYPE<PIXELTYPE,DIMENSION>::Pointer &image )
{
  typedef itk::IMAGETYPE<PIXELTYPE,DIMENSION> ImageType;
  typedef ImageType::Pointer                  ImagePointerType;

  ImagePointerType* rawPointer = new ImagePointerType( image );
  Rcpp::XPtr< ImagePointerType > xptr( rawPointer , true ) ;

  Rcpp::S4 itkImage( std::string( "antsImage" ) );
  itkImage.slot( "pixeltype" ) = "PIXELTYPE";
  itkImage.slot( "dimension" ) = DIMENSION;
  itkImage.slot( "components" ) = image->GetNumberOfComponentsPerPixel();
  itkImage.slot( "pointer") = xptr;

  return(wrap(itkImage));
}
*/

typedef itk::DiffusionTensor3D<double> DTDouble;
typedef itk::DiffusionTensor3D<float>  DTFloat;

template <> inline
SEXP wrap( const itk::Image<DTDouble,3>::Pointer &image )
{
  Rcpp::Rcout << "wrap( itk::Image<DiffusionTensor3D<double>( image ) )" << std::endl;
  /*
  typedef itk::Image<DTDouble,3>        ImageType;
  typedef ImageType::Pointer            ImagePointerType;
  //typedef itkImageFinalizer<double,2> FinalizerType;
  //typedef Rcpp::XPtr<ImagePointerType,PreserveStorage,FinalizerType::Finalize> XPtrType;
  typedef Rcpp::XPtr<ImagePointerType> XPtrType;

  ImagePointerType* rawPointer = new ImagePointerType( image );
  XPtrType xptr( rawPointer , true );

  Rcpp::S4 itkImage( std::string( "antsImage" ) );
  itkImage.slot( "pixeltype" ) = "DTdouble";
  itkImage.slot( "dimension" ) = 3;
  itkImage.slot( "components" ) = 6;
  itkImage.slot( "isVector" ) = true;
  itkImage.slot( "pointer") = xptr;

  return(wrap(itkImage));
  */
  return(NULL)
}

template <> inline
SEXP wrap( const itk::Image<DTFloat,3>::Pointer &image )
{
  Rcpp::Rcout << "wrap( itk::Image<DiffusionTensor3D<double>( image ) )" << std::endl;
  wrap(NULL);

  /*
  typedef itk::Image<DTFloat,3>         ImageType;
  typedef ImageType::Pointer            ImagePointerType;
  //typedef itkImageFinalizer<double,2> FinalizerType;
  //typedef Rcpp::XPtr<ImagePointerType,PreserveStorage,FinalizerType::Finalize> XPtrType;
  typedef Rcpp::XPtr<ImagePointerType> XPtrType;

  ImagePointerType* rawPointer = new ImagePointerType( image );
  XPtrType xptr( rawPointer , true );

  Rcpp::S4 itkImage( std::string( "antsImage" ) );
  itkImage.slot( "pixeltype" ) = "DTfloat";
  itkImage.slot( "dimension" ) = 3;
  itkImage.slot( "components" ) = 6;
  itkImage.slot( "isVector" ) = true;
  itkImage.slot( "pointer") = xptr;

  return(wrap(itkImage));
  */
}

template <> inline
itk::Image<DTDouble,3>::Pointer as( SEXP itkImageR )
{
  const unsigned int Dim = 3;
  typedef itk::Image<DTDouble,Dim>        ImageType;

  Rcpp::S4 itkImageObject( itkImageR );

  if (!itkImageObject.is( "antsImage") ||
      (Rcpp::as<int>(itkImageObject.slot("dimension")) != Dim) ||
      (Rcpp::as<int>(itkImageObject.slot("components")) != 6) )
    {
    Rcpp::stop( "Invalid S4 object type in itk::Image<DTdouble,3>::Pointer Rcpp::as()");
    }

  if ( Rcpp::as<std::string>(itkImageObject.slot("pixeltype")) == "DTdouble" ) {
    XPtr<typename ImageType::Pointer> xptr( static_cast<SEXP>( itkImageObject.slot("pointer") ));
    return *xptr;
  }

  // never reached
  Rcpp::Rcout << "Unsupported pixeltype: " <<
    Rcpp::as<std::string>(itkImageObject.slot("pixeltype")) << std::endl;
  return NULL;

}

template <> inline
itk::Image<DTFloat,3>::Pointer as( SEXP itkImageR )
{
  const unsigned int Dim = 3;
  typedef itk::Image<DTFloat,Dim>        ImageType;

  Rcpp::S4 itkImageObject( itkImageR );

  if (!itkImageObject.is( "antsImage") ||
      (Rcpp::as<int>(itkImageObject.slot("dimension")) != Dim) ||
      (Rcpp::as<int>(itkImageObject.slot("components")) != 6) )
    {
    Rcpp::stop( "Invalid S4 object type in itk::Image<DTfloat,3>::Pointer Rcpp::as()");
    }

  if ( Rcpp::as<std::string>(itkImageObject.slot("pixeltype")) == "DTfloat" ) {
    XPtr<typename ImageType::Pointer> xptr( static_cast<SEXP>( itkImageObject.slot("pointer") ));
    return *xptr;
  }

  // never reached
  Rcpp::Rcout << "Unsupported pixeltype: " <<
    Rcpp::as<std::string>(itkImageObject.slot("pixeltype")) << std::endl;
  return NULL;

}

}

#endif
