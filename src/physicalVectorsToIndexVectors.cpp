
#include <algorithm>
#include <vector>
#include <string>
#include <RcppANTsR.h>

#include "antsUtilities.h"
#include "itkTransform.h"
#include "itkImage.h"
#include "itkDiffusionTensor3D.h"
#include "itkFixedArray.h"



template<class ImageType>
void physicalVectorsToIndexVectors( SEXP r_image, SEXP r_inverse )
{
  using ImagePointerType = typename ImageType::Pointer;
  using IteratorType = typename itk::ImageRegionIteratorWithIndex<ImageType>;
  using PixelType = typename ImageType::PixelType;
  using VectorType = typename itk::FixedArray<double, ImageType::ImageDimension>;

  ImagePointerType image = Rcpp::as<ImagePointerType>(r_image);
  bool inverse = Rcpp::as<bool>( r_inverse );

  IteratorType it( image, image->GetLargestPossibleRegion() );
  while (!it.IsAtEnd()) {
    PixelType pix = it.Get();

    VectorType pVec;
    for ( unsigned int i=0; i<ImageType::ImageDimension; i++ ) {
      pVec[i] = pix[i];
    }

    VectorType lVec;
    if ( inverse ) {
      image->TransformLocalVectorToPhysicalVector(pVec, lVec);
    }
    else {
      image->TransformPhysicalVectorToLocalVector(pVec, lVec);
    }

    for ( unsigned int i=0; i<ImageType::ImageDimension; i++ ) {
      pix[i] = lVec[i];
    }

    it.Set(pix);

    ++it;
  }

}



RcppExport SEXP physicalVectorsToIndexVectors( SEXP r_image, SEXP r_inverse )
{
try
{
  Rcpp::S4 image(r_image);
  std::string pixeltype = Rcpp::as<std::string>( image.slot("precision") );
  int dimension = Rcpp::as<int>( image.slot("dimension"));

  if ( pixeltype=="double") {
    using PixelType = double;
    if (dimension==2) {
      using ImageType = itk::VectorImage<PixelType, 2>;
      physicalVectorsToIndexVectors<ImageType>(r_image, r_inverse);
    }
    if (dimension==3) {
      using ImageType = itk::VectorImage<PixelType, 3>;
      physicalVectorsToIndexVectors<ImageType>(r_image, r_inverse);
    }
    if (dimension==4) {
      using ImageType = itk::VectorImage<PixelType, 4>;
      physicalVectorsToIndexVectors<ImageType>(r_image, r_inverse);
    }
  }
  else if (pixeltype=="float") {
    using PixelType = float;
    if (dimension==2) {
      using ImageType = itk::VectorImage<PixelType, 2>;
      physicalVectorsToIndexVectors<ImageType>(r_image, r_inverse);
    }
    if (dimension==3) {
      using ImageType = itk::VectorImage<PixelType, 3>;
      physicalVectorsToIndexVectors<ImageType>(r_image, r_inverse);
    }
    if (dimension==4) {
      using ImageType = itk::VectorImage<PixelType, 4>;
      physicalVectorsToIndexVectors<ImageType>(r_image, r_inverse);
    }
  }
  else if (pixeltype=="unsigned int") {
    using PixelType = unsigned int;
    if (dimension==2) {
      using ImageType = itk::VectorImage<PixelType, 2>;
      physicalVectorsToIndexVectors<ImageType>(r_image, r_inverse);
    }
    if (dimension==3) {
      using ImageType = itk::VectorImage<PixelType, 3>;
      physicalVectorsToIndexVectors<ImageType>(r_image, r_inverse);
    }
    if (dimension==4) {
      using ImageType = itk::VectorImage<PixelType, 4>;
      physicalVectorsToIndexVectors<ImageType>(r_image, r_inverse);
    }
  }
  else if (pixeltype=="unsigned char") {
    using PixelType = unsigned char;
    if (dimension==2) {
      using ImageType = itk::VectorImage<PixelType, 2>;
      physicalVectorsToIndexVectors<ImageType>(r_image, r_inverse);
    }
    if (dimension==3) {
      using ImageType = itk::VectorImage<PixelType, 3>;
      physicalVectorsToIndexVectors<ImageType>(r_image, r_inverse);
    }
    if (dimension==4) {
      using ImageType = itk::VectorImage<PixelType, 4>;
      physicalVectorsToIndexVectors<ImageType>(r_image, r_inverse);
    }
  }

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
