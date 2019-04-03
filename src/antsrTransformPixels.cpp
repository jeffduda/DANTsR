
#include <algorithm>
#include <vector>
#include <string>
#include <RcppANTsR.h>

#include "antsUtilities.h"
#include "itkAffineTransform.h"
#include "itkCenteredAffineTransform.h"
#include "itkEuler2DTransform.h"
#include "itkEuler3DTransform.h"
#include "itkRigid2DTransform.h"
#include "itkRigid3DTransform.h"
#include "itkCenteredRigid2DTransform.h"
#include "itkCenteredEuler3DTransform.h"
#include "itkSimilarity2DTransform.h"
#include "itkCenteredSimilarity2DTransform.h"
#include "itkSimilarity3DTransform.h"
#include "itkQuaternionRigidTransform.h"
#include "itkTranslationTransform.h"
#include "itkResampleImageFilter.h"
#include "itkTransformFileReader.h"
#include "itkCompositeTransform.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkDisplacementFieldTransform.h"
#include "itkConstantBoundaryCondition.h"

#include "itkBSplineInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkGaussianInterpolateImageFunction.h"
#include "itkInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkWindowedSincInterpolateImageFunction.h"
#include "itkLabelImageGaussianInterpolateImageFunction.h"
#include "itkTransformFileWriter.h"
#include "itkDiffusionTensor3D.h"


/*
template< class TransformType >
Rcpp::XPtr<typename TransformType::Pointer> antsrTransformGetXPtr()
{
  using TransformPointerType = typename TransformType::Pointer;
  TransformPointerType transformPtr = TransformType::New();

  TransformPointerType * rawPointer = new TransformPointerType( transformPtr );
  Rcpp::XPtr<TransformPointerType> xptr( rawPointer, true );
  return xptr;
}*/

/*
// Apply transform to point
template< class PrecisionType, unsigned int Dimension >
SEXP antsrTransform_TransformPoint( SEXP r_transform, SEXP r_point )
{
  using TransformType = itk::Transform<PrecisionType,Dimension,Dimension>;
  using TransformPointerType = typename TransformType::Pointer;
  using InputPointType = typename TransformType::InputPointType;
  using OutputPointType = typename TransformType::OutputPointType;

  Rcpp::S4 transform( r_transform );
  std::string type = Rcpp::as<std::string>( transform.slot("type") );

  TransformPointerType itkTransform = Rcpp::as<TransformPointerType>( r_transform );

  Rcpp::NumericMatrix inPoints( r_point );
  Rcpp::NumericMatrix outPoints( inPoints.nrow(), inPoints.ncol() );

  for (unsigned int n=0; n<inPoints.nrow(); n++) {

    InputPointType inItkPoint;
    for (unsigned int i=0; i<InputPointType::PointDimension; i++)
      {
      inItkPoint[i] = inPoints(n,i);
      }

    OutputPointType outItkPoint = itkTransform->TransformPoint( inItkPoint );

    for (unsigned int i=0; i<OutputPointType::PointDimension; i++)
      {
      outPoints(n,i) = outItkPoint[i];
      }
    }

  return Rcpp::wrap(outPoints);
}


RcppExport SEXP antsrTransform_TransformPoint( SEXP r_transform, SEXP r_point )
{
try
{
  Rcpp::S4 transform( r_transform );

  std::string precision = Rcpp::as<std::string>( transform.slot("precision") );
  unsigned int dimension = Rcpp::as<int>( transform.slot("dimension") );

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
      return antsrTransform_TransformPoint<PrecisionType,4>( r_transform, r_point  );
      }
    else if( dimension == 3 )
	    {
      return antsrTransform_TransformPoint<PrecisionType,3>( r_transform, r_point  );
	    }
    else if( dimension == 2 )
	    {
      return antsrTransform_TransformPoint<PrecisionType,2>( r_transform, r_point );
	    }
	  }
  else if( precision == "float" )
    {
    using PrecisionType = float;
    if( dimension == 4 )
	    {
      return antsrTransform_TransformPoint<PrecisionType,4>( r_transform, r_point );
      }
    else if( dimension == 3 )
	    {
      return antsrTransform_TransformPoint<PrecisionType,3>( r_transform, r_point );
	    }
    else if( dimension == 2 )
	    {
      return antsrTransform_TransformPoint<PrecisionType,2>( r_transform, r_point );
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

// Apply transform to vector
template< class PrecisionType, unsigned int Dimension >
SEXP antsrTransform_TransformVector( SEXP r_transform, SEXP r_vector )
{
  using TransformType = itk::Transform<PrecisionType,Dimension,Dimension>;
  using TransformPointerType = typename TransformType::Pointer;
  using InputVectorType = typename TransformType::InputVectorType;
  using OutputVectorType = typename TransformType::OutputVectorType;

  Rcpp::S4 transform( r_transform );
  std::string type = Rcpp::as<std::string>( transform.slot("type") );

  TransformPointerType itkTransform = Rcpp::as<TransformPointerType>( r_transform );
  Rcpp::NumericMatrix inVectors( r_vector );
  Rcpp::NumericMatrix outVectors( inVectors.nrow(), inVectors.ncol() );

  for (unsigned int n=0; n<inVectors.nrow(); n++ ) {
    InputVectorType inItkVector;
    for (unsigned int i=0; i<InputVectorType::Dimension; i++)
      {
      inItkVector[i] = inVectors(n,i);
      }

    OutputVectorType outItkVector = itkTransform->TransformVector( inItkVector );

    for (unsigned int i=0; i<OutputVectorType::Dimension; i++)
      {
      outVectors(n,i) = outItkVector[i];
      }
    }

  return Rcpp::wrap(outVectors);
}


RcppExport SEXP antsrTransform_TransformVector( SEXP r_transform, SEXP r_vector )
{
try
{
  Rcpp::S4 transform( r_transform );

  std::string precision = Rcpp::as<std::string>( transform.slot("precision") );
  unsigned int dimension = Rcpp::as<int>( transform.slot("dimension") );

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
      return antsrTransform_TransformVector<PrecisionType,4>( r_transform, r_vector  );
      }
    else if( dimension == 3 )
	    {
      return antsrTransform_TransformVector<PrecisionType,3>( r_transform, r_vector  );
	    }
    else if( dimension == 2 )
	    {
      return antsrTransform_TransformVector<PrecisionType,2>( r_transform, r_vector );
	    }
	  }
  else if( precision == "float" )
    {
    using PrecisionType = float;
    if( dimension == 4 )
	    {
      return antsrTransform_TransformVector<PrecisionType,4>( r_transform, r_vector );
      }
    else if( dimension == 3 )
	    {
      return antsrTransform_TransformVector<PrecisionType,3>( r_transform, r_vector );
	    }
    else if( dimension == 2 )
	    {
      return antsrTransform_TransformVector<PrecisionType,2>( r_transform, r_vector );
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
*/


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
  std::string pixeltype = Rcpp::as<std::string>( image.slot("precision") );
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
  std::string pixeltype = Rcpp::as<std::string>( image.slot("precision") );
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
      antsrTransform_TransformPixels<PrecisionType,4>( r_transform, r_image, r_pixelAs );
      }
    else if( dimension == 3 )
	    {
      antsrTransform_TransformPixels<PrecisionType,3>( r_transform, r_image, r_pixelAs );
	    }
    else if( dimension == 2 )
	    {
      antsrTransform_TransformPixels<PrecisionType,2>( r_transform, r_image, r_pixelAs );
	    }
	  }
  else if( precision == "float" )
    {
    using PrecisionType = float;
    if( dimension == 4 )
	    {
      antsrTransform_TransformPixels<PrecisionType,4>( r_transform, r_image, r_pixelAs );
      }
    else if( dimension == 3 )
	    {
      antsrTransform_TransformPixels<PrecisionType,3>( r_transform, r_image, r_pixelAs );
	    }
    else if( dimension == 2 )
	    {
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

template<class PrecisionType>
unsigned int antsrTransform_GetDimensionFromFile( SEXP r_filename )
{
  using TransformReaderType = itk::TransformFileReaderTemplate<PrecisionType>;
  using TransformListType = typename TransformReaderType::TransformListType;
  using TransformIteratorType = typename TransformListType::const_iterator;

  std::string filename = Rcpp::as<std::string>( r_filename );

  typename TransformReaderType::Pointer reader = TransformReaderType::New();
  reader->SetFileName( filename );
  reader->Update();
  const typename TransformReaderType::TransformListType * transformList = reader->GetTransformList();

  unsigned int dim = 0;
  unsigned int count = 0;

  for (TransformIteratorType i = transformList->begin(); i != transformList->end(); ++i)
  {
    unsigned int inDim = (*i)->GetInputSpaceDimension();
    unsigned int outDim = (*i)->GetOutputSpaceDimension();
    if (count == 0 ) {
      dim = inDim;
    }
    else {
      if (inDim != dim ) {
        Rcpp::stop( "All transforms must have the same dimension");
      }
    }

    if ( inDim != outDim ) {
      Rcpp::stop( "Must have same input and output dimensions");
    }

    ++count;
  }

  return(dim);
}

template< class PrecisionType, unsigned int Dimension >
SEXP antsrTransform_Read( SEXP r_filename, SEXP r_precision )
{
  using TransformBaseType = itk::Transform<PrecisionType,Dimension,Dimension>;
  using TransformBasePointerType = typename TransformBaseType::Pointer;
  using CompositeTransformType = typename itk::CompositeTransform<PrecisionType, Dimension>;
  using TransformReaderType = itk::TransformFileReaderTemplate<PrecisionType>;
  using TransformListType = typename TransformReaderType::TransformListType;
  using TransformIteratorType = typename TransformListType::const_iterator;

  std::string filename = Rcpp::as<std::string>( r_filename );

  typename TransformReaderType::Pointer reader = TransformReaderType::New();
  reader->SetFileName( filename );
  reader->Update();

  const typename TransformReaderType::TransformListType * transformList = reader->GetTransformList();
  for (TransformIteratorType i = transformList->begin(); i != transformList->end(); ++i)
  {
    unsigned int inDim = (*i)->GetInputSpaceDimension();
    unsigned int outDim = (*i)->GetOutputSpaceDimension();
    if ( inDim != Dimension ) {
      Rcpp::stop( "Invalid input space dimension");
    }
    if ( outDim != Dimension ) {
      Rcpp::stop( "Invalid output space dimension");
    }
  }

  Rcpp::S4 antsrTransform( "antsrTransform" );
  antsrTransform.slot("dimension") = Dimension;
  antsrTransform.slot("precision") = Rcpp::as<std::string>( r_precision );

  TransformBasePointerType transform;

  if ( transformList->size() > 1 )
  {
    typename CompositeTransformType::Pointer comp_transform = CompositeTransformType::New();
    for (TransformIteratorType i = transformList->begin(); i != transformList->end(); ++i)
    {
      comp_transform->AddTransform( dynamic_cast<TransformBaseType *>( i->GetPointer()) );
    }
    transform = dynamic_cast<TransformBaseType *>(comp_transform.GetPointer());
  }
  else
  {
    transform = dynamic_cast<TransformBaseType *>( transformList->front().GetPointer() );
  }

  std::string type = transform->GetNameOfClass();
  antsrTransform.slot("type") = type;

  TransformBasePointerType * rawPointer = new TransformBasePointerType( transform );
  Rcpp::XPtr<TransformBasePointerType> xptr( rawPointer, true );
  antsrTransform.slot("pointer") = xptr;

  return antsrTransform;
}

RcppExport SEXP antsrTransform_Read( SEXP r_filename, SEXP r_dimension, SEXP r_precision )
{
try
{
  unsigned int dimension = Rcpp::as<int>( r_dimension );
  std::string precision = Rcpp::as<std::string>( r_precision );

  if ( (precision != "float") && (precision != "double"))
    {
    Rcpp::stop( "Precision must be 'float' or 'double'");
    }

  if (dimension == NA_INTEGER) {
    if ( precision=="float" ) {
      dimension = antsrTransform_GetDimensionFromFile<float>( r_filename );
    }
    else if (precision == "double" ) {
      dimension = antsrTransform_GetDimensionFromFile<double>( r_filename );
    }
  }

  if ( precision == "float")
  {
    using PrecisionType = float;
    if ( dimension == 4 )
    {
      return antsrTransform_Read<PrecisionType,4>( r_filename, r_precision );
    }
    else if ( dimension == 3)
    {
      return antsrTransform_Read<PrecisionType,3>( r_filename, r_precision );
    }
    else if ( dimension == 2 )
    {
      return antsrTransform_Read<PrecisionType,2>( r_filename, r_precision );
    }
    else
    {
      Rcpp::stop( "Unsupported dimension" );
    }
  }
  else if ( precision == "double")
  {
    using PrecisionType = double;
    if ( dimension == 4 )
    {
      return antsrTransform_Read<PrecisionType,4>( r_filename, r_precision );
    }
    else if ( dimension == 3)
    {
      return antsrTransform_Read<PrecisionType,3>( r_filename, r_precision );
    }
    else if ( dimension == 2 )
    {
      return antsrTransform_Read<PrecisionType,2>( r_filename, r_precision );
    }
    else
    {
      Rcpp::stop( "Unsupported dimension" );
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


template< class PrecisionType, unsigned int Dimension >
SEXP antsrTransform_Compose( SEXP r_list, SEXP r_precision )
{
  using TransformBaseType = itk::Transform<PrecisionType,Dimension,Dimension>;
  using TransformBasePointerType = typename TransformBaseType::Pointer;
  using CompositeTransformType = typename itk::CompositeTransform<PrecisionType, Dimension>;

  Rcpp::List transforms( r_list );
  Rcpp::S4 antsrTransform( "antsrTransform" );
  antsrTransform.slot("dimension") = Dimension;
  antsrTransform.slot("precision") = Rcpp::as<std::string>( r_precision );

  typename CompositeTransformType::Pointer comp_transform = CompositeTransformType::New();
  for ( unsigned int i=0; i<transforms.size(); i++ )
    {
    TransformBasePointerType t = Rcpp::as<TransformBasePointerType>( transforms[i] );
    comp_transform->AddTransform( t );
    Rcpp::S4 tran(transforms[i]);
    }

  TransformBasePointerType transform = dynamic_cast<TransformBaseType *>(comp_transform.GetPointer());

  std::string type = comp_transform->GetNameOfClass();
  antsrTransform.slot("type") = type;

  TransformBasePointerType * rawPointer = new TransformBasePointerType( transform );
  Rcpp::XPtr<TransformBasePointerType> xptr( rawPointer, true );
  antsrTransform.slot("pointer") = xptr;

  return antsrTransform;
}

RcppExport SEXP antsrTransform_Compose( SEXP r_list, SEXP r_dimension, SEXP r_precision )
{
try
{
  unsigned int dimension = Rcpp::as<int>( r_dimension );
  std::string precision = Rcpp::as<std::string>( r_precision );

  if ( (precision != "float") && (precision != "double"))
    {
    Rcpp::stop( "Precision must be 'float' or 'double'");
    }

  if ( precision == "float")
  {
    using PrecisionType = float;
    if ( dimension == 4 )
    {
      return antsrTransform_Compose<PrecisionType,4>( r_list, r_precision );
    }
    else if ( dimension == 3)
    {
      return antsrTransform_Compose<PrecisionType,3>( r_list, r_precision );
    }
    else if ( dimension == 2 )
    {
      return antsrTransform_Compose<PrecisionType,2>( r_list, r_precision );
    }
    else
    {
      Rcpp::stop( "Unsupported dimension" );
    }
  }
  else if ( precision == "double")
  {
    using PrecisionType = double;
    if ( dimension == 4 )
    {
      return antsrTransform_Compose<PrecisionType,4>( r_list, r_precision );
    }
    else if ( dimension == 3)
    {
      return antsrTransform_Compose<PrecisionType,3>( r_list, r_precision );
    }
    else if ( dimension == 2 )
    {
      return antsrTransform_Compose<PrecisionType,2>( r_list, r_precision );
    }
    else
    {
      Rcpp::stop( "Unsupported dimension" );
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


template< class PrecisionType, unsigned int Dimension >
SEXP antsrTransform_FromDisplacementField( SEXP r_field, std::string precision )
{
  using TransformType = itk::Transform<PrecisionType,Dimension,Dimension>;
  using TransformPointerType = typename TransformType::Pointer;
  using DisplacementFieldTransformType = typename itk::DisplacementFieldTransform<PrecisionType, Dimension>;
  using DisplacementFieldType = typename DisplacementFieldTransformType::DisplacementFieldType;
  using VectorType = typename DisplacementFieldType::PixelType;
  using IteratorType = itk::ImageRegionIteratorWithIndex<DisplacementFieldType>;

  // Displacement field is an itk::Image with vector pixels, while in ANTsR we use the
  // itk::VectorImage class for multichannel data. So we must copy the field
  // and pass it to the transform
  using AntsrFieldType = itk::VectorImage<PrecisionType, Dimension>;
  using AntsrFieldPointerType = typename AntsrFieldType::Pointer;

  AntsrFieldPointerType antsrField = Rcpp::as<AntsrFieldPointerType>( r_field );
  typename DisplacementFieldType::Pointer itkField = DisplacementFieldType::New();
  itkField->SetRegions( antsrField->GetLargestPossibleRegion() );
  itkField->SetSpacing( antsrField->GetSpacing() );
  itkField->SetOrigin( antsrField->GetOrigin() );
  itkField->SetDirection( antsrField->GetDirection() );
  itkField->Allocate();

  IteratorType it( itkField, itkField->GetLargestPossibleRegion() );
  while ( !it.IsAtEnd() )
  {
    typename AntsrFieldType::PixelType vec = antsrField->GetPixel( it.GetIndex() );
    VectorType dvec;
    for ( unsigned int i=0; i<Dimension; i++)
      {
      dvec[i] = vec[i];
      }
    itkField->SetPixel(it.GetIndex(), dvec);
    ++it;
  }

  typename DisplacementFieldTransformType::Pointer displacementTransform =
    DisplacementFieldTransformType::New();
  displacementTransform->SetDisplacementField( itkField );

  TransformPointerType transform = dynamic_cast<TransformType *>( displacementTransform.GetPointer() );

  Rcpp::S4 antsrTransform( "antsrTransform" );
  antsrTransform.slot("dimension") = Dimension;
  antsrTransform.slot("precision") = precision;
  std::string type = displacementTransform->GetNameOfClass();
  antsrTransform.slot("type") = type;
  TransformPointerType * rawPointer = new TransformPointerType( transform );
  Rcpp::XPtr<TransformPointerType> xptr( rawPointer, true );
  antsrTransform.slot("pointer") = xptr;

  return antsrTransform;
}

RcppExport SEXP antsrTransform_FromDisplacementField( SEXP r_field )
{
try
{
  Rcpp::S4 antsImage( r_field );
  unsigned int dimension = Rcpp::as<int>( antsImage.slot("dimension"));
  unsigned int components = Rcpp::as<int>( antsImage.slot("components"));
  std::string precision = Rcpp::as<std::string>( antsImage.slot("pixeltype"));

  if ( (precision != "float") && (precision != "double") )
  {
    Rcpp::stop("Field must have pixeltype of either float or double");
  }

  if ( components != dimension )
  {
    Rcpp::stop("Field must have number of pixel compenents equal to image dimension");
  }

  if ( precision == "float")
  {
    using PrecisionType = float;
    if ( dimension == 4 )
    {
      return antsrTransform_FromDisplacementField<PrecisionType,4>( r_field, precision );
    }
    else if ( dimension == 3)
    {
      return antsrTransform_FromDisplacementField<PrecisionType,3>( r_field, precision );
    }
    else if ( dimension == 2 )
    {
      return antsrTransform_FromDisplacementField<PrecisionType,2>( r_field, precision );
    }
    else
    {
      Rcpp::stop( "Unsupported dimension" );
    }
  }
  else if ( precision == "double")
  {
    using PrecisionType = float;
    if ( dimension == 4 )
    {
      return antsrTransform_FromDisplacementField<PrecisionType,4>( r_field, precision );
    }
    else if ( dimension == 3)
    {
      return antsrTransform_FromDisplacementField<PrecisionType,3>( r_field, precision );
    }
    else if ( dimension == 2 )
    {
      return antsrTransform_FromDisplacementField<PrecisionType,2>( r_field, precision );
    }
    else
    {
      Rcpp::stop( "Unsupported dimension" );
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


// Apply transform to point
template< class PrecisionType, unsigned int Dimension >
SEXP antsrTransform_Inverse( SEXP r_transform )
{
  using TransformType = itk::Transform<PrecisionType,Dimension,Dimension>;
  using TransformPointerType = typename TransformType::Pointer;

  Rcpp::S4 transform( r_transform );
  std::string type = Rcpp::as<std::string>( transform.slot("type") );

  TransformPointerType itkTransform = Rcpp::as<TransformPointerType>( r_transform );
  if ( !itkTransform->IsLinear() )
  {
    Rcpp::stop("Only linear transforms may be inverted with this method");
  }

  TransformPointerType inverse = itkTransform->GetInverseTransform();
  return Rcpp::wrap(inverse);
}


RcppExport SEXP antsrTransform_Inverse( SEXP r_transform )
{
try
{
  Rcpp::S4 transform( r_transform );

  std::string precision = Rcpp::as<std::string>( transform.slot("precision") );
  unsigned int dimension = Rcpp::as<int>( transform.slot("dimension") );

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
      return antsrTransform_Inverse<PrecisionType,4>( r_transform );
      }
    else if( dimension == 3 )
	    {
      return antsrTransform_Inverse<PrecisionType,3>( r_transform );
	    }
    else if( dimension == 2 )
	    {
      return antsrTransform_Inverse<PrecisionType,2>( r_transform );
	    }
	  }
  else if( precision == "float" )
    {
    using PrecisionType = float;
    if( dimension == 4 )
	    {
      return antsrTransform_Inverse<PrecisionType,4>( r_transform );
      }
    else if( dimension == 3 )
	    {
      return antsrTransform_Inverse<PrecisionType,3>( r_transform );
	    }
    else if( dimension == 2 )
	    {
      return antsrTransform_Inverse<PrecisionType,2>( r_transform );
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

template< class PrecisionType, unsigned int Dimension >
SEXP antsrTransform_Write( SEXP r_transform ,SEXP filename_ )
{
  using TransformType = itk::Transform<PrecisionType,Dimension,Dimension>;
  using TransformPointerType = typename TransformType::Pointer;
  using TransformWriterType = itk::TransformFileWriter;

  std::string filename = Rcpp::as<std::string>(filename_);
  Rcpp::S4 transform( r_transform );
  std::string type = Rcpp::as<std::string>( transform.slot("type") );

  TransformPointerType itkTransform = Rcpp::as<TransformPointerType>( r_transform );
  typename TransformWriterType::Pointer transformWriter = TransformWriterType::New();
  transformWriter->SetInput( itkTransform );
  transformWriter->SetFileName( filename.c_str() );
  transformWriter->Update();
  return Rcpp::wrap(true);

}
RcppExport SEXP antsrTransform_Write( SEXP r_transform ,SEXP filename_)
{
  try
  {
  Rcpp::S4 transform( r_transform );

  std::string precision = Rcpp::as<std::string>( transform.slot("precision") );
  unsigned int dimension = Rcpp::as<int>( transform.slot("dimension") );
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
	      return antsrTransform_Write<PrecisionType,4>( r_transform ,filename_);
      }
    else if( dimension == 3 )
	    {
      return antsrTransform_Write<PrecisionType,3>( r_transform ,filename_);
	    }
    else if( dimension == 2 )
	    {
      return antsrTransform_Write<PrecisionType,2>( r_transform ,filename_);
	    }
	  }
  else if( precision == "float" )
    {
    using PrecisionType = float;
    if( dimension == 4 )
	    {
      return antsrTransform_Write<PrecisionType,4>( r_transform ,filename_);
      }
    else if( dimension == 3 )
	    {
      return antsrTransform_Write<PrecisionType,3>( r_transform ,filename_);
	    }
    else if( dimension == 2 )
	    {
      return antsrTransform_Write<PrecisionType,2>( r_transform ,filename_);
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
