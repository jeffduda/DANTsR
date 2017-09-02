
#include <algorithm>
#include <vector>
#include <string>
#include <RcppDANTsR.h>

#include "itkDiffusionTensor3D.h"
#include "itkCastImageFilter.h"
//#include "itkVectorToDiffusionTensor3DAccessor.h"
#include "itkImportImageFilter.h"

#include "itkAddImageFilter.h"
#include "itkDefaultConvertPixelTraits.h"
#include "itkMultiplyImageFilter.h"
#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkImageBase.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNeighborhoodIterator.h"
#include "itkPermuteAxesImageFilter.h"
#include "itkCentralDifferenceImageFunction.h"
#include "itkContinuousIndex.h"

typedef itk::DiffusionTensor3D<double> DTDouble;
typedef itk::DiffusionTensor3D<float>  DTFloat;

template< class VectorImageType, class TensorImageType >
SEXP dtiFromVectorImage( SEXP r_vImg, bool copyData=true )
{
  Rcpp::Rcout << "dtiFromVectorImage<VectorImageType,TensorImageType>" << std::endl;
  typedef typename VectorImageType::Pointer VectorImagePointer;
  VectorImagePointer vImg = Rcpp::as<VectorImagePointer>( r_vImg );

  typename TensorImageType::Pointer dti;


  typename VectorImageType::PixelType v;
  typename TensorImageType::PixelType d;

  if (copyData) {
    dti = TensorImageType::New();
    dti->SetRegions( vImg->GetRequestedRegion() );
    dti->SetSpacing( vImg->GetSpacing() );
    dti->SetOrigin( vImg->GetOrigin() );
    dti->SetDirection( vImg->GetDirection() );
    dti->Allocate();

    typename itk::ImageRegionIteratorWithIndex<VectorImageType>
      it( vImg, vImg->GetRequestedRegion() );
    for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
      typename TensorImageType::PixelType dt;
      typename VectorImageType::PixelType vec = it.Get();
      for (unsigned int i=0; i<6; i++) {
        dt[i] = vec[i];
      }
      dti->SetPixel(it.GetIndex(), dt);
    }
  }
  else {
    typedef typename VectorImageType::InternalPixelType ValueType;
    ValueType * dat = vImg->GetBufferPointer();
    //typedef itk::itkImportImageFilter
  }
  //typedef itk::CastImageFilter<VectorImageType,TensorImageType> CastType;
  //typename CastType::Pointer cast = CastType::New();
  //cast->SetInput( vImg );

  //typedef typename TensorImageType::Pointer TensorImagePointer;
  //TensorImagePointer dti = cast->GetOutput();



  return Rcpp::wrap( NULL );


}

RcppExport SEXP dtiFromVectorImage( SEXP r_vImg )
{
Rcpp::Rcout << "dtiFromVectorImage(vImg)" << std::endl;
try
{
  if( r_vImg == NULL )
    {
    Rcpp::stop("Unspecified Arguments");
    }

  Rcpp::S4 antsImage(r_vImg);
  std::string pixeltype = Rcpp::as<std::string>( antsImage.slot("pixeltype") );
  unsigned int dimension = Rcpp::as<int>( antsImage.slot("dimension") );

  if ( dimension != 3 )
    {
    Rcpp::stop("Unsupported image dimension");
    }

  if( pixeltype == "double" )
    {
    typedef itk::Image<DTDouble,3>        TensorImageType;
    typedef itk::VectorImage<double,3>    VectorImageType;
    dtiFromVectorImage<VectorImageType,TensorImageType>( r_vImg );
	  }
  else if( pixeltype == "float" )
    {
    typedef itk::Image<DTFloat,3>        TensorImageType;
    typedef itk::VectorImage<float,3>    VectorImageType;
    dtiFromVectorImage<VectorImageType,TensorImageType>( r_vImg );
    }
  else
    {
    Rcpp::stop("Unsupported PixelType");
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
  //forward_exception_to_r( exc ) ;
  }
catch(...)
  {
	Rcpp::stop("c++ exception (unknown reason)");
  }
return Rcpp::wrap(NA_REAL); //not reached
}
