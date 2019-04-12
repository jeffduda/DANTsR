
#include <algorithm>
#include <vector>
#include <string>
#include <RcppANTsR.h>

#include "antsUtilities.h"
#include "itkTransform.h"
#include "itkImage.h"
#include "itkDiffusionTensor3D.h"



// Apply transform to vector pixel data
template< class TransformType, class ImageType >
void antsrTransform_WarpVectorPixels( SEXP r_transform, SEXP r_image )
{
  using TransformPointerType = typename TransformType::Pointer;
  using ImagePointerType = typename ImageType::Pointer;
  using IteratorType = typename itk::ImageRegionIteratorWithIndex<ImageType>;
  using PixelType = typename ImageType::PixelType;
  using InputVectorType = typename TransformType::InputVectorType;
  using OutputVectorType = typename TransformType::OutputVectorType;
  using PointType = typename ImageType::PointType;
  using WarpType = itk::DisplacementFieldTransform< typename TransformType::ParametersValueType, TransformType::InputSpaceDimension >;
  using WarpPointerType = typename WarpType::Pointer;

  TransformPointerType transform2 = Rcpp::as<TransformPointerType>( r_transform );
  WarpPointerType transform = dynamic_cast<WarpType *>( transform2.GetPointer() );
  ImagePointerType image = Rcpp::as<ImagePointerType>( r_image );

  IteratorType it( image, image->GetLargestPossibleRegion() );
  while ( !it.IsAtEnd() )
  {
    PixelType pix = it.Get();
    InputVectorType inVector;
    for (unsigned int i=0; i<image->GetNumberOfComponentsPerPixel(); i++)
    {
      inVector[i] = pix[i];
    }

    PointType pt;
    image->TransformIndexToPhysicalPoint( it.GetIndex(), pt );

    OutputVectorType outVec = transform->TransformVector(inVector, pt);

    for (unsigned int i=0; i<image->GetNumberOfComponentsPerPixel(); i++)
    {
      pix[i] = outVec[i];
    }

    it.Set(pix);

    ++it;
  }

}

// Apply transform to vector pixel data
template< class TransformType>
void antsrTransform_WarpVectorPixels( SEXP r_transform, SEXP r_image )
{
  Rcpp::S4 image( r_image );
  std::string pixeltype = Rcpp::as<std::string>( image.slot("pixeltype") );

  if ( pixeltype=="double" ) {
    using PixelType = double;
    using ImageType = itk::VectorImage<PixelType,TransformType::InputSpaceDimension>;
    antsrTransform_WarpVectorPixels<TransformType, ImageType>( r_transform, r_image );
  }
  else if ( pixeltype=="float" ) {
    using PixelType = float;
    using ImageType = itk::VectorImage<PixelType,TransformType::InputSpaceDimension>;
    antsrTransform_WarpVectorPixels<TransformType, ImageType>( r_transform, r_image );
  }
  else if ( pixeltype=="unsigned int" ) {
    using PixelType = unsigned int;
    using ImageType = itk::VectorImage<PixelType,TransformType::InputSpaceDimension>;
    antsrTransform_WarpVectorPixels<TransformType, ImageType>( r_transform, r_image );
  }
  else if ( pixeltype=="unsigned char" ) {
    using PixelType = unsigned char;
    using ImageType = itk::VectorImage<PixelType,TransformType::InputSpaceDimension>;
    antsrTransform_WarpVectorPixels<TransformType, ImageType>( r_transform, r_image );
  }
}

// Apply transform to vector pixel data
template< class TransformType, class ImageType >
void antsrTransform_WarpCovariantVectorPixels( SEXP r_transform, SEXP r_image )
{
  using TransformPointerType = typename TransformType::Pointer;
  using ImagePointerType = typename ImageType::Pointer;
  using IteratorType = typename itk::ImageRegionIteratorWithIndex<ImageType>;
  using PixelType = typename ImageType::PixelType;
  using InputVectorType = typename TransformType::InputCovariantVectorType;
  using OutputVectorType = typename TransformType::OutputCovariantVectorType;
  using PointType = typename ImageType::PointType;
  using WarpType = itk::DisplacementFieldTransform< typename TransformType::ParametersValueType, TransformType::InputSpaceDimension >;
  using WarpPointerType = typename WarpType::Pointer;

  TransformPointerType transform2 = Rcpp::as<TransformPointerType>( r_transform );
  WarpPointerType transform = dynamic_cast<WarpType *>( transform2.GetPointer() );
  ImagePointerType image = Rcpp::as<ImagePointerType>( r_image );

  IteratorType it( image, image->GetLargestPossibleRegion() );
  while ( !it.IsAtEnd() )
  {
    PixelType pix = it.Get();
    InputVectorType inVector;
    for (unsigned int i=0; i<image->GetNumberOfComponentsPerPixel(); i++)
    {
      inVector[i] = pix[i];
    }

    PointType pt;
    image->TransformIndexToPhysicalPoint( it.GetIndex(), pt );

    OutputVectorType outVec = transform->TransformCovariantVector(inVector, pt);

    for (unsigned int i=0; i<image->GetNumberOfComponentsPerPixel(); i++)
    {
      pix[i] = outVec[i];
    }

    it.Set(pix);

    ++it;
  }

}

// Apply transform to vector pixel data
template< class TransformType>
void antsrTransform_WarpCovariantVectorPixels( SEXP r_transform, SEXP r_image )
{
  Rcpp::S4 image( r_image );
  std::string pixeltype = Rcpp::as<std::string>( image.slot("pixeltype") );

  if ( pixeltype=="double" ) {
    using PixelType = double;
    using ImageType = itk::VectorImage<PixelType,TransformType::InputSpaceDimension>;
    antsrTransform_WarpCovariantVectorPixels<TransformType, ImageType>( r_transform, r_image );
    }
  else if ( pixeltype=="float" ) {
    using PixelType = float;
    using ImageType = itk::VectorImage<PixelType,TransformType::InputSpaceDimension>;
    antsrTransform_WarpCovariantVectorPixels<TransformType, ImageType>( r_transform, r_image );
  }
  else if ( pixeltype=="unsigned int" ) {
    using PixelType = unsigned int;
    using ImageType = itk::VectorImage<PixelType,TransformType::InputSpaceDimension>;
    antsrTransform_WarpCovariantVectorPixels<TransformType, ImageType>( r_transform, r_image );
  }
  else if ( pixeltype=="unsigned char" ) {
    using PixelType = unsigned char;
    using ImageType = itk::VectorImage<PixelType,TransformType::InputSpaceDimension>;
    antsrTransform_WarpCovariantVectorPixels<TransformType, ImageType>( r_transform, r_image );
  }
}

// Apply transform to dt pixel data
template< class TransformType, class ImageType >
void antsrTransform_WarpDiffusionTensorPixels( SEXP r_transform, SEXP r_image )
{
  using TransformPointerType = typename TransformType::Pointer;
  using ImagePointerType = typename ImageType::Pointer;
  using IteratorType = typename itk::ImageRegionIteratorWithIndex<ImageType>;
  using PixelType = typename ImageType::PixelType;
  using InputTensorType = typename TransformType::InputDiffusionTensor3DType;
  using OutputTensorType = typename TransformType::OutputDiffusionTensor3DType;
  using PointType = typename ImageType::PointType;
  using WarpType = itk::DisplacementFieldTransform< typename TransformType::ParametersValueType, TransformType::InputSpaceDimension >;
  using WarpPointerType = typename WarpType::Pointer;

  TransformPointerType transform2 = Rcpp::as<TransformPointerType>( r_transform );
  WarpPointerType transform = dynamic_cast<WarpType *>( transform2.GetPointer() );

  ImagePointerType image = Rcpp::as<ImagePointerType>( r_image );


  IteratorType it( image, image->GetLargestPossibleRegion() );
  while ( !it.IsAtEnd() )
  {
    PixelType pix = it.Get();
    InputTensorType inDT;
    for (unsigned int i=0; i<image->GetNumberOfComponentsPerPixel(); i++)
    {
      inDT[i] = pix[i];
    }

    PointType pt;
    image->TransformIndexToPhysicalPoint( it.GetIndex(), pt );

    OutputTensorType outDT = transform->TransformDiffusionTensor3D(inDT,pt);

    for (unsigned int i=0; i<image->GetNumberOfComponentsPerPixel(); i++)
    {
      pix[i] = outDT[i];
    }

    it.Set(pix);

    ++it;
  }

}

// Apply transform to vector pixel data
template< class TransformType>
void antsrTransform_WarpDiffusionTensorPixels( SEXP r_transform, SEXP r_image )
{
  Rcpp::S4 image( r_image );
  std::string pixeltype = Rcpp::as<std::string>( image.slot("pixeltype") );

  if ( pixeltype=="double" ) {
    using PixelType = double;
    using ImageType = itk::VectorImage<PixelType,TransformType::InputSpaceDimension>;
    antsrTransform_WarpDiffusionTensorPixels<TransformType, ImageType>( r_transform, r_image );
  }
  else if ( pixeltype=="float" ) {
    using PixelType = float;
    using ImageType = itk::VectorImage<PixelType,TransformType::InputSpaceDimension>;
    antsrTransform_WarpDiffusionTensorPixels<TransformType, ImageType>( r_transform, r_image );
  }
  else if ( pixeltype=="unsigned int" ) {
    using PixelType = unsigned int;
    using ImageType = itk::VectorImage<PixelType,TransformType::InputSpaceDimension>;
    antsrTransform_WarpDiffusionTensorPixels<TransformType, ImageType>( r_transform, r_image );
  }
  else if ( pixeltype=="unsigned char" ) {
    using PixelType = unsigned char;
    using ImageType = itk::VectorImage<PixelType,TransformType::InputSpaceDimension>;
    antsrTransform_WarpDiffusionTensorPixels<TransformType, ImageType>( r_transform, r_image );
  }

}

// Apply transform to pixel data
template< class PrecisionType, unsigned int Dimension >
void antsrTransform_WarpPixels( SEXP r_transform, SEXP r_image, SEXP r_pixelAs )
{
  //using TransformType = itk::DisplacementFieldTransform<PrecisionType,Dimension>;
  using TransformType = itk::Transform<PrecisionType,Dimension,Dimension>;
  std::string pixelAs = Rcpp::as<std::string>( r_pixelAs );

  if ( pixelAs=="vector" ) {
    antsrTransform_WarpVectorPixels<TransformType>( r_transform, r_image );
  }
  else if ( pixelAs=="covariantvector" ) {
    antsrTransform_WarpCovariantVectorPixels<TransformType>( r_transform, r_image );
  }
  else if ( pixelAs=="diffusiontensor" ) {
    antsrTransform_WarpDiffusionTensorPixels<TransformType>( r_transform, r_image );
  }
  else {
    Rcpp::stop("Invalid pixel.as selection");
  }

}


// Apply transform to vector pixel data
template< class TransformType, class ImageType >
void antsrTransform_TransformVectorPixels( SEXP r_transform, SEXP r_image )
{
  using TransformPointerType = typename TransformType::Pointer;
  using ImagePointerType = typename ImageType::Pointer;
  using IteratorType = typename itk::ImageRegionIteratorWithIndex<ImageType>;
  using PixelType = typename ImageType::PixelType;
  using InputVectorType = typename TransformType::InputVectorType;
  using OutputVectorType = typename TransformType::OutputVectorType;

  TransformPointerType transform = Rcpp::as<TransformPointerType>( r_transform );
  ImagePointerType image = Rcpp::as<ImagePointerType>( r_image );

  IteratorType it( image, image->GetLargestPossibleRegion() );
  while ( !it.IsAtEnd() )
  {
    PixelType pix = it.Get();
    InputVectorType inVector;
    for (unsigned int i=0; i<image->GetNumberOfComponentsPerPixel(); i++)
    {
      inVector[i] = pix[i];
    }

    OutputVectorType outVec = transform->TransformVector(inVector);

    for (unsigned int i=0; i<image->GetNumberOfComponentsPerPixel(); i++)
    {
      pix[i] = outVec[i];
    }

    it.Set(pix);

    ++it;
  }

}

// Apply transform to vector pixel data
template< class TransformType>
void antsrTransform_TransformVectorPixels( SEXP r_transform, SEXP r_image )
{
  Rcpp::S4 image( r_image );
  std::string pixeltype = Rcpp::as<std::string>( image.slot("pixeltype") );
  int dimension = Rcpp::as<int>( image.slot("dimension"));

  if ( pixeltype=="double" ) {
    using PixelType = double;
    if ( dimension==2  ) {
      using ImageType = itk::VectorImage<PixelType,2>;
      antsrTransform_TransformVectorPixels<TransformType, ImageType>( r_transform, r_image );
    }
    else if ( dimension==3  ) {
      using ImageType = itk::VectorImage<PixelType,3>;
      antsrTransform_TransformVectorPixels<TransformType, ImageType>( r_transform, r_image );
    }
    else if ( dimension==4  ) {
      using ImageType = itk::VectorImage<PixelType,4>;
      antsrTransform_TransformVectorPixels<TransformType, ImageType>( r_transform, r_image );
    }
  }
  else if ( pixeltype=="float" ) {
    using PixelType = float;
    if ( dimension==2  ) {
      using ImageType = itk::VectorImage<PixelType,2>;
      antsrTransform_TransformVectorPixels<TransformType, ImageType>( r_transform, r_image );
    }
    else if ( dimension==3  ) {
      using ImageType = itk::VectorImage<PixelType,3>;
      antsrTransform_TransformVectorPixels<TransformType, ImageType>( r_transform, r_image );
    }
    else if ( dimension==4  ) {
      using ImageType = itk::VectorImage<PixelType,4>;
      antsrTransform_TransformVectorPixels<TransformType, ImageType>( r_transform, r_image );
    }
  }
  else if ( pixeltype=="unsigned int" ) {
    using PixelType = unsigned int;
    if ( dimension==2  ) {
      using ImageType = itk::VectorImage<PixelType,2>;
      antsrTransform_TransformVectorPixels<TransformType, ImageType>( r_transform, r_image );
    }
    else if ( dimension==3  ) {
      using ImageType = itk::VectorImage<PixelType,3>;
      antsrTransform_TransformVectorPixels<TransformType, ImageType>( r_transform, r_image );
    }
    else if ( dimension==4  ) {
      using ImageType = itk::VectorImage<PixelType,4>;
      antsrTransform_TransformVectorPixels<TransformType, ImageType>( r_transform, r_image );
    }
  }
  else if ( pixeltype=="unsigned char" ) {
    using PixelType = unsigned char;
    if ( dimension==2  ) {
      using ImageType = itk::VectorImage<PixelType,2>;
      antsrTransform_TransformVectorPixels<TransformType, ImageType>( r_transform, r_image );
    }
    else if ( dimension==3  ) {
      using ImageType = itk::VectorImage<PixelType,3>;
      antsrTransform_TransformVectorPixels<TransformType, ImageType>( r_transform, r_image );
    }
    else if ( dimension==4  ) {
      using ImageType = itk::VectorImage<PixelType,4>;
      antsrTransform_TransformVectorPixels<TransformType, ImageType>( r_transform, r_image );
    }
  }
}

// Apply transform to vector pixel data
template< class TransformType, class ImageType >
void antsrTransform_TransformCovariantVectorPixels( SEXP r_transform, SEXP r_image )
{
  using TransformPointerType = typename TransformType::Pointer;
  using ImagePointerType = typename ImageType::Pointer;
  using IteratorType = typename itk::ImageRegionIteratorWithIndex<ImageType>;
  using PixelType = typename ImageType::PixelType;
  using InputVectorType = typename TransformType::InputCovariantVectorType;
  using OutputVectorType = typename TransformType::OutputCovariantVectorType;

  TransformPointerType transform = Rcpp::as<TransformPointerType>( r_transform );
  ImagePointerType image = Rcpp::as<ImagePointerType>( r_image );

  IteratorType it( image, image->GetLargestPossibleRegion() );
  while ( !it.IsAtEnd() )
  {
    PixelType pix = it.Get();
    InputVectorType inVector;
    for (unsigned int i=0; i<image->GetNumberOfComponentsPerPixel(); i++)
    {
      inVector[i] = pix[i];
    }

    OutputVectorType outVec = transform->TransformCovariantVector(inVector);

    for (unsigned int i=0; i<image->GetNumberOfComponentsPerPixel(); i++)
    {
      pix[i] = outVec[i];
    }

    it.Set(pix);

    ++it;
  }

}

// Apply transform to vector pixel data
template< class TransformType>
void antsrTransform_TransformCovariantVectorPixels( SEXP r_transform, SEXP r_image )
{
  Rcpp::S4 image( r_image );
  std::string pixeltype = Rcpp::as<std::string>( image.slot("pixeltype") );
  int dimension = Rcpp::as<int>( image.slot("dimension"));

  if ( pixeltype=="double" ) {
    using PixelType = double;
    if ( dimension==2  ) {
      using ImageType = itk::VectorImage<PixelType,2>;
      antsrTransform_TransformCovariantVectorPixels<TransformType, ImageType>( r_transform, r_image );
    }
    else if ( dimension==3  ) {
      using ImageType = itk::VectorImage<PixelType,3>;
      antsrTransform_TransformCovariantVectorPixels<TransformType, ImageType>( r_transform, r_image );
    }
    else if ( dimension==4  ) {
      using ImageType = itk::VectorImage<PixelType,4>;
      antsrTransform_TransformCovariantVectorPixels<TransformType, ImageType>( r_transform, r_image );
    }
  }
  else if ( pixeltype=="float" ) {
    using PixelType = float;
    if ( dimension==2  ) {
      using ImageType = itk::VectorImage<PixelType,2>;
      antsrTransform_TransformCovariantVectorPixels<TransformType, ImageType>( r_transform, r_image );
    }
    else if ( dimension==3  ) {
      using ImageType = itk::VectorImage<PixelType,3>;
      antsrTransform_TransformCovariantVectorPixels<TransformType, ImageType>( r_transform, r_image );
    }
    else if ( dimension==4  ) {
      using ImageType = itk::VectorImage<PixelType,4>;
      antsrTransform_TransformCovariantVectorPixels<TransformType, ImageType>( r_transform, r_image );
    }
  }
  else if ( pixeltype=="unsigned int" ) {
    using PixelType = unsigned int;
    if ( dimension==2  ) {
      using ImageType = itk::VectorImage<PixelType,2>;
      antsrTransform_TransformCovariantVectorPixels<TransformType, ImageType>( r_transform, r_image );
    }
    else if ( dimension==3  ) {
      using ImageType = itk::VectorImage<PixelType,3>;
      antsrTransform_TransformCovariantVectorPixels<TransformType, ImageType>( r_transform, r_image );
    }
    else if ( dimension==4  ) {
      using ImageType = itk::VectorImage<PixelType,4>;
      antsrTransform_TransformCovariantVectorPixels<TransformType, ImageType>( r_transform, r_image );
    }
  }
  else if ( pixeltype=="unsigned char" ) {
    using PixelType = unsigned char;
    if ( dimension==2  ) {
      using ImageType = itk::VectorImage<PixelType,2>;
      antsrTransform_TransformCovariantVectorPixels<TransformType, ImageType>( r_transform, r_image );
    }
    else if ( dimension==3  ) {
      using ImageType = itk::VectorImage<PixelType,3>;
      antsrTransform_TransformCovariantVectorPixels<TransformType, ImageType>( r_transform, r_image );
    }
    else if ( dimension==4  ) {
      using ImageType = itk::VectorImage<PixelType,4>;
      antsrTransform_TransformCovariantVectorPixels<TransformType, ImageType>( r_transform, r_image );
    }
  }
}

// Apply transform to dt pixel data
template< class TransformType, class ImageType >
void antsrTransform_TransformDiffusionTensorPixels( SEXP r_transform, SEXP r_image )
{
  using TransformPointerType = typename TransformType::Pointer;
  using ImagePointerType = typename ImageType::Pointer;
  using IteratorType = typename itk::ImageRegionIteratorWithIndex<ImageType>;
  using PixelType = typename ImageType::PixelType;
  using InputTensorType = typename TransformType::InputDiffusionTensor3DType;
  using OutputTensorType = typename TransformType::OutputDiffusionTensor3DType;

  TransformPointerType transform = Rcpp::as<TransformPointerType>( r_transform );
  ImagePointerType image = Rcpp::as<ImagePointerType>( r_image );

  IteratorType it( image, image->GetLargestPossibleRegion() );
  while ( !it.IsAtEnd() )
  {
    PixelType pix = it.Get();
    InputTensorType inDT;
    for (unsigned int i=0; i<image->GetNumberOfComponentsPerPixel(); i++)
    {
      inDT[i] = pix[i];
    }

    OutputTensorType outDT = transform->TransformDiffusionTensor3D(inDT);

    for (unsigned int i=0; i<image->GetNumberOfComponentsPerPixel(); i++)
    {
      pix[i] = outDT[i];
    }

    it.Set(pix);

    ++it;
  }

}

// Apply transform to vector pixel data
template< class TransformType>
void antsrTransform_TransformDiffusionTensorPixels( SEXP r_transform, SEXP r_image )
{
  Rcpp::S4 image( r_image );
  std::string pixeltype = Rcpp::as<std::string>( image.slot("pixeltype") );
  int dimension = Rcpp::as<int>( image.slot("dimension"));

  if ( pixeltype=="double" ) {
    using PixelType = double;
    if ( dimension==2  ) {
      using ImageType = itk::VectorImage<PixelType,2>;
      antsrTransform_TransformDiffusionTensorPixels<TransformType, ImageType>( r_transform, r_image );
    }
    else if ( dimension==3  ) {
      using ImageType = itk::VectorImage<PixelType,3>;
      antsrTransform_TransformDiffusionTensorPixels<TransformType, ImageType>( r_transform, r_image );
    }
    else if ( dimension==4  ) {
      using ImageType = itk::VectorImage<PixelType,4>;
      antsrTransform_TransformDiffusionTensorPixels<TransformType, ImageType>( r_transform, r_image );
    }
  }
  else if ( pixeltype=="float" ) {
    using PixelType = float;
    if ( dimension==2  ) {
      using ImageType = itk::VectorImage<PixelType,2>;
      antsrTransform_TransformDiffusionTensorPixels<TransformType, ImageType>( r_transform, r_image );
    }
    else if ( dimension==3  ) {
      using ImageType = itk::VectorImage<PixelType,3>;
      antsrTransform_TransformDiffusionTensorPixels<TransformType, ImageType>( r_transform, r_image );
    }
    else if ( dimension==4  ) {
      using ImageType = itk::VectorImage<PixelType,4>;
      antsrTransform_TransformDiffusionTensorPixels<TransformType, ImageType>( r_transform, r_image );
    }
  }
  else if ( pixeltype=="unsigned int" ) {
    using PixelType = unsigned int;
    if ( dimension==2  ) {
      using ImageType = itk::VectorImage<PixelType,2>;
      antsrTransform_TransformDiffusionTensorPixels<TransformType, ImageType>( r_transform, r_image );
    }
    else if ( dimension==3  ) {
      using ImageType = itk::VectorImage<PixelType,3>;
      antsrTransform_TransformDiffusionTensorPixels<TransformType, ImageType>( r_transform, r_image );
    }
    else if ( dimension==4  ) {
      using ImageType = itk::VectorImage<PixelType,4>;
      antsrTransform_TransformDiffusionTensorPixels<TransformType, ImageType>( r_transform, r_image );
    }
  }
  else if ( pixeltype=="unsigned char" ) {
    using PixelType = unsigned char;
    if ( dimension==2  ) {
      using ImageType = itk::VectorImage<PixelType,2>;
      antsrTransform_TransformDiffusionTensorPixels<TransformType, ImageType>( r_transform, r_image );
    }
    else if ( dimension==3  ) {
      using ImageType = itk::VectorImage<PixelType,3>;
      antsrTransform_TransformDiffusionTensorPixels<TransformType, ImageType>( r_transform, r_image );
    }
    else if ( dimension==4  ) {
      using ImageType = itk::VectorImage<PixelType,4>;
      antsrTransform_TransformDiffusionTensorPixels<TransformType, ImageType>( r_transform, r_image );
    }
  }

}


// Apply transform to pixel data
template< class PrecisionType, unsigned int Dimension >
void antsrTransform_TransformPixels( SEXP r_transform, SEXP r_image, SEXP r_pixelAs )
{
  using TransformType = itk::Transform<PrecisionType,Dimension,Dimension>;
  std::string pixelAs = Rcpp::as<std::string>( r_pixelAs );

  if ( pixelAs=="vector" ) {
    antsrTransform_TransformVectorPixels<TransformType>( r_transform, r_image );
  }
  else if ( pixelAs=="covariantvector" ) {
    antsrTransform_TransformCovariantVectorPixels<TransformType>( r_transform, r_image );
  }
  else if ( pixelAs=="diffusiontensor" ) {
    antsrTransform_TransformDiffusionTensorPixels<TransformType>( r_transform, r_image );
  }
  else {
    Rcpp::stop("Invalid pixel.as selection");
  }

}


RcppExport SEXP antsrTransform_TransformPixels( SEXP r_transform, SEXP r_image, SEXP r_pixelAs )
{
try
{
  Rcpp::S4 transform( r_transform );

  std::string precision = Rcpp::as<std::string>( transform.slot("precision") );
  unsigned int dimension = Rcpp::as<int>( transform.slot("dimension") );
  std::string txType = Rcpp::as<std::string>( transform.slot("type") );

  bool warp = ( txType == "DisplacementFieldTransform" );

  if ( (dimension < 1) || (dimension > 4) )
    {
    Rcpp::stop("Unsupported image dimension");
    }

  if ( (precision != "float") && (precision != "double"))
    {
    Rcpp::stop( "Precision must be 'float' or 'double'");
    }

  if( precision == "double" )
    {
    using PrecisionType = double;
    if( dimension == 4 )
	    {
      warp ?
        antsrTransform_WarpPixels<PrecisionType,4>( r_transform, r_image, r_pixelAs )  :
        antsrTransform_TransformPixels<PrecisionType,4>( r_transform, r_image, r_pixelAs );
      }
    else if( dimension == 3 )
	    {
      warp ?
        antsrTransform_WarpPixels<PrecisionType,3>( r_transform, r_image, r_pixelAs )  :
        antsrTransform_TransformPixels<PrecisionType,3>( r_transform, r_image, r_pixelAs );
	    }
    else if( dimension == 2 )
	    {
      warp ?
        antsrTransform_WarpPixels<PrecisionType,2>( r_transform, r_image, r_pixelAs )  :
        antsrTransform_TransformPixels<PrecisionType,2>( r_transform, r_image, r_pixelAs );
	    }
	  }
  else if( precision == "float" )
    {
    using PrecisionType = float;
    if( dimension == 4 )
	    {
      warp ?
        antsrTransform_WarpPixels<PrecisionType,4>( r_transform, r_image, r_pixelAs )  :
        antsrTransform_TransformPixels<PrecisionType,4>( r_transform, r_image, r_pixelAs );
      }
    else if( dimension == 3 )
	    {
      warp ?
        antsrTransform_WarpPixels<PrecisionType,3>( r_transform, r_image, r_pixelAs )  :
        antsrTransform_TransformPixels<PrecisionType,3>( r_transform, r_image, r_pixelAs );
	    }
    else if( dimension == 2 )
	    {
      warp ?
        antsrTransform_WarpPixels<PrecisionType,2>( r_transform, r_image, r_pixelAs )  :
        antsrTransform_TransformPixels<PrecisionType,2>( r_transform, r_image, r_pixelAs );
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
