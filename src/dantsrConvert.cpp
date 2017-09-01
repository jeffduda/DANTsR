
#include <algorithm>
#include <vector>
#include <string>
#include <RcppDANTsR.h>

#include "itkDiffusionTensor3D.h"
#include "itkCastImageFilter.h"
#include "itkVectorToDiffusionTensor3DAccessor.h"

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
SEXP dtiFromVectorImage( SEXP r_vImg )
{
  typedef typename VectorImageType::Pointer VectorImagePointer;
  VectorImagePointer vImg = Rcpp::as<VectorImagePointer>( r_vImg );


  typename VectorImageType::PixelType v;
  typename TensorImageType::PixelType d;



  //typedef itk::CastImageFilter<VectorImageType,TensorImageType> CastType;
  //typename CastType::Pointer cast = CastType::New();
  //cast->SetInput( vImg );

  //typedef typename TensorImageType::Pointer TensorImagePointer;
  //TensorImagePointer dti = cast->GetOutput();



  return Rcpp::wrap(NULL);


}

RcppExport SEXP dtiFromVectorImage( SEXP r_vImg )
{
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
