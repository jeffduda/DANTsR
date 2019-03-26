
#include <algorithm>
#include <vector>
#include <string>
#include <RcppANTsR.h>

#include "itkDiffusionTensor3D.h"
#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkEigenAnalysis3DImageFilter.h"

typedef itk::DiffusionTensor3D<double> DTDouble;
typedef itk::DiffusionTensor3D<float>  DTFloat;

namespace itk
  {

   /* Class to return the nth element for data types with more than one component
    * and the first value for one element data type
    */
   template< typename DataType >
   class DANTsR_GetDataValue
   {
   public:
     DANTsR_GetDataValue() {}
     ~DANTsR_GetDataValue() {}

     typedef itk::DefaultConvertPixelTraits<DataType> ConvertType;
     typedef itk::NumericTraits<DataType>             TraitsType;
     typedef typename TraitsType::ValueType           ValueType;

     bool operator!=(const DANTsR_GetDataValue &) const
     {
       return false;
     }

     bool operator==(const DANTsR_GetDataValue & other) const
     {
       return !( *this != other );
     }

     static inline ValueType NthValue (unsigned long i, const DataType & A ) {
       if ( TraitsType::GetLength(A)==1 ) {
         return ConvertType::GetNthComponent(0,A);
       }
       else {
         return ConvertType::GetNthComponent(i,A);
       }

     }
   };

   namespace Functor
   {
     template< typename VectorType, typename DataType >
     class DANTsR_FractionalAnisotropy
     {
     public:
       DANTsR_FractionalAnisotropy() {}
       ~DANTsR_FractionalAnisotropy() {}

       typedef itk::DefaultConvertPixelTraits<VectorType> ConvertType;
       typedef itk::NumericTraits<VectorType>             TraitsType;
       typedef typename TraitsType::ValueType             ValueType;
       typedef DANTsR_GetDataValue<VectorType>            GetValueType;

       bool operator!=(const DANTsR_FractionalAnisotropy &) const
       {
         return false;
       }

      bool operator==(const DANTsR_FractionalAnisotropy & other) const
       {
         return !( *this != other );
       }

      inline DataType operator()(const VectorType & A ) const
       {
         typedef typename itk::DiffusionTensor3D<DataType> TensorType;
         TensorType dt( A.GetDataPointer() );
         return dt.GetFractionalAnisotropy();
       }
     };

     template< typename VectorType, typename DataType >
     class DANTsR_RelativeAnisotropy
     {
     public:
       DANTsR_RelativeAnisotropy() {}
       ~DANTsR_RelativeAnisotropy() {}

       typedef itk::DefaultConvertPixelTraits<VectorType> ConvertType;
       typedef itk::NumericTraits<VectorType>             TraitsType;
       typedef typename TraitsType::ValueType             ValueType;
       typedef DANTsR_GetDataValue<VectorType>            GetValueType;

       bool operator!=(const DANTsR_RelativeAnisotropy &) const
       {
         return false;
       }

      bool operator==(const DANTsR_RelativeAnisotropy & other) const
       {
         return !( *this != other );
       }

      inline DataType operator()(const VectorType & A ) const
       {
         typedef typename itk::DiffusionTensor3D<DataType> TensorType;
         TensorType dt( A.GetDataPointer() );
         return dt.GetRelativeAnisotropy();
       }
     };

     template< typename VectorType, typename DataType >
     class DANTsR_InnerScalarProduct
     {
     public:
       DANTsR_InnerScalarProduct() {}
       ~DANTsR_InnerScalarProduct() {}

       typedef itk::DefaultConvertPixelTraits<VectorType> ConvertType;
       typedef itk::NumericTraits<VectorType>             TraitsType;
       typedef typename TraitsType::ValueType             ValueType;
       typedef DANTsR_GetDataValue<VectorType>            GetValueType;

       bool operator!=(const DANTsR_InnerScalarProduct &) const
       {
         return false;
       }

      bool operator==(const DANTsR_InnerScalarProduct & other) const
       {
         return !( *this != other );
       }

      inline DataType operator()(const VectorType & A ) const
       {
         typedef typename itk::DiffusionTensor3D<DataType> TensorType;
         TensorType dt( A.GetDataPointer() );
         return dt.GetInnerScalarProduct();
       }
     };

   }
 }



template< class VectorImageType, class ImageType >
SEXP dtiFilters( SEXP r_dti, SEXP r_filter )
{

  typedef typename ImageType::PixelType PixelType;
  typedef typename ImageType::Pointer ImagePointerType;

  typedef typename VectorImageType::PixelType VectorType;
  typedef typename VectorImageType::Pointer VectorImagePointerType;

  VectorImagePointerType image = Rcpp::as<VectorImagePointerType>(r_dti);

  std::string filter = Rcpp::as< std::string >( r_filter );

  if ( filter == "fractionalanisotropy" )
  {
    typedef itk::UnaryFunctorImageFilter<VectorImageType, ImageType,
      itk::Functor::DANTsR_FractionalAnisotropy<VectorType,PixelType>  > FilterType;

    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput( image );
    filter->Update();
    typename ImageType::Pointer outImage = filter->GetOutput();
    outImage->DisconnectPipeline();

    return Rcpp::wrap(outImage);
  }
  else if ( filter == "relativeanisotropy" )
  {
    typedef itk::UnaryFunctorImageFilter<VectorImageType, ImageType,
      itk::Functor::DANTsR_RelativeAnisotropy<VectorType,PixelType>  > FilterType;

    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput( image );
    filter->Update();
    typename ImageType::Pointer outImage = filter->GetOutput();
    outImage->DisconnectPipeline();

    return Rcpp::wrap(outImage);
  }
  else if ( filter == "innerscalarproduct" )
  {
    typedef itk::UnaryFunctorImageFilter<VectorImageType, ImageType,
      itk::Functor::DANTsR_InnerScalarProduct<VectorType,PixelType>  > FilterType;

    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput( image );
    filter->Update();
    typename ImageType::Pointer outImage = filter->GetOutput();
    outImage->DisconnectPipeline();

    return Rcpp::wrap(outImage);
  }
  else if ( filter == "eigensystem")
  {
    typedef itk::EigenAnalysis3DImageFilter<VectorImageType, ImageType, VectorImageType> EigenFilterType;
    typename EigenFilterType::Pointer eigen = EigenFilterType::New();
    eigen->SetInput( image );
    eigen->Update();
    ImagePointerType e1 = eigen->GetEigenValue(1);
    ImagePointerType e2 = eigen->GetEigenValue(2);
    ImagePointerType e3 = eigen->GetEigenValue(3);
    VectorImagePointerType v1 = eigen->GetEigenVector(1);
    VectorImagePointerType v2 = eigen->GetEigenVector(2);
    VectorImagePointerType v3 = eigen->GetEigenVector(3);

    return Rcpp::wrap( Rcpp::List::create(
      Rcpp::Named("e1")=e1, Rcpp::Named("e2")=e2, Rcpp::Named("e3")=e3,
      Rcpp::Named("v1")=v1, Rcpp::Named("v2")=v2, Rcpp::Named("v3")=v3 ) );
  }

  Rcpp::Rcout << "Unknown DTI filter passed: " << filter << std::endl;
  return(Rcpp::wrap(NA_REAL));

}



RcppExport SEXP dtiFilters( SEXP r_dti, SEXP r_filter ) {
try
  {
  if( r_dti == nullptr || r_filter == nullptr )
    {
      Rcpp::Rcout << "Unspecified Arguments" << std::endl;
      return Rcpp::wrap( NA_REAL );
    }

  Rcpp::S4 antsimage( r_dti ) ;
  std::string pixeltype = Rcpp::as< std::string >( antsimage.slot( "pixeltype" ) );
  unsigned int dimension = Rcpp::as< int >( antsimage.slot( "dimension" ) );
  unsigned int components = Rcpp::as< int >( antsimage.slot( "components" ) );
  //bool isVector = Rcpp::as<bool>( antsimage.slot("isVector") );

  if ( components != 6 ) {
    Rcpp::Rcout << "Input must have 6 components" << std::endl;
    return Rcpp::wrap(NA_REAL);
  }

  if ( dimension != 3 ) {
    Rcpp::Rcout << "Only 3D images are supported" << std::endl;
  }
  const int ImageDimension = 3;

  if( pixeltype == "double" )
    {
    typedef double PixelType;
    typedef itk::Image< PixelType , ImageDimension >      ImageType;
    typedef itk::VectorImage< PixelType, ImageDimension > VectorImageType;
    return  dtiFilters< VectorImageType, ImageType >( r_dti, r_filter );
    }
  else if( pixeltype == "float" )
    {
    typedef float PixelType;
    typedef itk::Image< PixelType , ImageDimension >      ImageType;
    typedef itk::VectorImage< PixelType, ImageDimension > VectorImageType;
    return dtiFilters< VectorImageType, ImageType >( r_dti, r_filter );
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
