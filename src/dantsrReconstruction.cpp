
#include <algorithm>
#include <vector>
#include <string>
#include <RcppANTsR.h>

#include "itkDiffusionTensor3D.h"
#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkDiffusionTensor3DReconstructionImageFilter.h"

typedef itk::DiffusionTensor3D<double> DTDouble;
typedef itk::DiffusionTensor3D<float>  DTFloat;


template< class VectorImageType >
SEXP dtiReconstruction( SEXP r_dwi, SEXP r_gradients, SEXP r_method )
{
  typedef typename VectorImageType::PixelType VectorType;
  typedef typename VectorImageType::Pointer VectorImagePointerType;
  typedef typename VectorType::ValueType ValueType;

  VectorImagePointerType image = Rcpp::as<VectorImagePointerType>(r_dwi);

  std::string method = Rcpp::as< std::string >( r_method );
  Rcpp::NumericMatrix gradients( r_gradients );

  if ( gradients.ncol() != 4 ) {
    Rcpp::stop( "Gradient matrix must have 4 columns: xdir, ydir, zdir, bvalue");
  }
  if ( gradients.nrow() != image->GetVectorLength()  ) {
    Rcpp::stop( "Number of matrix rows must equals components in image");
  }

  if ( method == "itk-svd" )
  {
    typedef itk::DiffusionTensor3DReconstructionImageFilter<ValueType,ValueType,ValueType>
      FilterType;

    typename FilterType::GradientDirectionType vect3d;
    typename FilterType::GradientDirectionContainerType::Pointer diffusionVectors =
      FilterType::GradientDirectionContainerType::New();

    ValueType bValue = 0;
    for ( unsigned int i=0; i<gradients.nrow(); i++) {
      ValueType thisB = gradients(i,3);
      if ( thisB != 0 ) {
        if (bValue==0) {
          bValue = thisB;
        }
        if ( thisB != bValue ) {
          Rcpp::stop("This method only supports one non-zero b-value");
        }
      }
      vect3d[0] = gradients(i,0);
      vect3d[1] = gradients(i,1);
      vect3d[2] = gradients(i,2);
      diffusionVectors->InsertElement( i, vect3d );
      //Rcpp::Rcout << "Set gradient: " << vect3d << " with bvalue=" << thisB << std::endl;
    }

    typename FilterType::Pointer filter = FilterType::New();
    filter->SetNumberOfThreads(1); // due to  netlib/dsvdc.c
    filter->SetGradientImage( diffusionVectors, image );
    filter->SetBValue( bValue );

    //filter->SetThreshold( lowValue ); add this as option?
    filter->Update();
    typename FilterType::OutputImageType::Pointer dtiImage = filter->GetOutput();

    //Copy to ANTsR compatible image
    VectorImagePointerType outImage = VectorImageType::New();
    outImage->SetRegions( dtiImage->GetBufferedRegion() );
    outImage->SetSpacing( dtiImage->GetSpacing() );
    outImage->SetOrigin( dtiImage->GetOrigin() );
    outImage->SetDirection( dtiImage->GetDirection() );
    outImage->SetVectorLength(6);
    outImage->Allocate();

    itk::ImageRegionIteratorWithIndex<typename FilterType::OutputImageType> it(dtiImage, dtiImage->GetBufferedRegion() );
    itk::ImageRegionIteratorWithIndex<VectorImageType> itOut(outImage, outImage->GetBufferedRegion() );

    it.GoToBegin();
    while( !it.IsAtEnd() ) {
      VectorType dt(6);
      dt.SetData( it.Get().GetDataPointer(), 6, false);
      itOut.Set(dt);
      ++it;
      ++itOut;
    }

    return Rcpp::wrap(outImage);
  }

  Rcpp::Rcout << "Unknown DTI reconstruction method passed: " << method << std::endl;
  return(Rcpp::wrap(NA_REAL));

}



RcppExport SEXP dtiReconstruction( SEXP r_dwi, SEXP r_gradients, SEXP r_method ) {
try
  {
  if( r_dwi == nullptr || r_method == nullptr )
    {
      Rcpp::Rcout << "Unspecified Arguments" << std::endl;
      return Rcpp::wrap( NA_REAL );
    }

  Rcpp::S4 antsimage( r_dwi ) ;
  std::string pixeltype = Rcpp::as< std::string >( antsimage.slot( "pixeltype" ) );
  unsigned int dimension = Rcpp::as< int >( antsimage.slot( "dimension" ) );
  //unsigned int components = Rcpp::as< int >( antsimage.slot( "components" ) );
  //bool isVector = Rcpp::as<bool>( antsimage.slot("isVector") );

  if ( dimension != 3 ) {
    Rcpp::Rcout << "Only 3D multichannel images are supported" << std::endl;
  }
  const int ImageDimension = 3;

  if( pixeltype == "double" )
    {
    typedef double PixelType;
    typedef itk::VectorImage< PixelType, ImageDimension > VectorImageType;
    return  dtiReconstruction< VectorImageType >( r_dwi, r_gradients, r_method );
    }
  else if( pixeltype == "float" )
    {
    typedef float PixelType;
    typedef itk::VectorImage< PixelType, ImageDimension > VectorImageType;
    return dtiReconstruction< VectorImageType >( r_dwi, r_gradients, r_method );
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
