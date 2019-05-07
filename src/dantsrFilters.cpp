
#include <algorithm>
#include <vector>
#include <string>
#include <RcppANTsR.h>

#include "itkDiffusionTensor3D.h"
#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkEigenAnalysis3DImageFilter.h"
#include "itkSymmetricEigenAnalysis.h"
#include "itkRGBPixel.h"
#include "itkRGBAPixel.h"
#include "TensorFunctions.h"

typedef itk::DiffusionTensor3D<double> DTDouble;
typedef itk::DiffusionTensor3D<float>  DTFloat;

/*
template< class VectorType >
MatrixType tensorAsDiffusionTensor3D VectorType dtVec )
{
  using TensorType = typename itk::DiffusionTensor3D< typename VectorType::ValueType >;

  TensorType dt()  ;

}
*/


template< class VectorType, class MatrixType >
MatrixType tensorAsMatrix( VectorType dtVec )
{
  MatrixType dtMat;
  dtMat[0][0] = dtVec[0];
  dtMat[0][1] = dtMat[1][0] = dtVec[1];
  dtMat[0][2] = dtMat[2][0] = dtVec[2];
  dtMat[1][1] = dtVec[3];
  dtMat[1][2] = dtVec[4];
  dtMat[2][2] = dtVec[5];
  return dtMat;
}

template< class MatrixType, class VectorType >
VectorType tensorAsVector( MatrixType dtMat )
{
  VectorType dtVec;
  itk::NumericTraits<VectorType>::SetLength(dtVec, 6);

  dtVec[0] = dtMat[0][0];
  dtVec[1] = dtMat[0][1];
  dtVec[2] = dtMat[0][2];
  dtVec[3] = dtMat[1][1];
  dtVec[4] = dtMat[1][2];
  dtVec[5] = dtMat[2][2];
  return dtVec;
}

template< class VectorType >
VectorType zeroNonFiniteTensors( VectorType dt )
{
  bool valid = true;
  VectorType dtOut;
  itk::NumericTraits<VectorType>::SetLength(dtOut,6);

  for (unsigned int i=0; i<6; i++)
    {
    valid = valid && (std::isfinite(dt[i]));
    }

  for (unsigned int i=0; i<6; i++)
    {
    if (valid)
      {
      dtOut[i] = dt[i];
      }
    else
      {
      dtOut[i] = 0;
      }
    }

  return(dtOut);

}

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

     template< typename VectorType, typename DataType >
     class DANTsR_Log
     {
     public:
       DANTsR_Log() {}
       ~DANTsR_Log() {}

       typedef itk::DefaultConvertPixelTraits<VectorType> ConvertType;
       typedef itk::NumericTraits<VectorType>             TraitsType;
       typedef typename TraitsType::ValueType             ValueType;
       typedef DANTsR_GetDataValue<VectorType>            GetValueType;

       bool operator!=(const DANTsR_Log &) const
       {
         return false;
       }

      bool operator==(const DANTsR_Log & other) const
       {
         return !( *this != other );
       }

      inline DataType operator()(const VectorType & A ) const
       {
        using TensorType = typename itk::DiffusionTensor3D< typename VectorType::ValueType>;
        TensorType dt( A.GetDataPointer() );
        TensorType dt2 = TensorLog<TensorType>( dt );
        VectorType logDt( dt2.GetDataPointer(),6 );
        return( logDt );

        /*
        using MatrixType = typename itk::Matrix<typename DataType::ValueType, 3,3>;
        using ArrayType = typename itk::Vector<typename DataType::ValueType, 3>;
        using SolverType = itk::SymmetricEigenAnalysis<MatrixType,ArrayType,MatrixType>;
        SolverType solver;
        ArrayType eigenValues;
        MatrixType eigenVectors;

        MatrixType dtMat = tensorAsMatrix<VectorType, MatrixType>(A);

        solver.SetDimension(3);
        solver.ComputeEigenValuesAndVectors(dtMat, eigenValues, eigenVectors);

        MatrixType eValues;
        eValues.Fill(0);

        for (unsigned int i=0; i<3; i++)
          {
          eValues(i,i) = std::log(eigenValues[i]);
          }

        MatrixType DTrec = eigenVectors * eValues * eigenVectors.GetTranspose();
        DataType dtv = tensorAsVector<MatrixType,DataType>(DTrec);
        return(zeroNonFiniteTensors<DataType>(dtv));
        */
       }

     };

     template< typename VectorType, typename DataType >
     class DANTsR_Exp
     {
     public:
       DANTsR_Exp() {}
       ~DANTsR_Exp() {}

       typedef itk::DefaultConvertPixelTraits<VectorType> ConvertType;
       typedef itk::NumericTraits<VectorType>             TraitsType;
       typedef typename TraitsType::ValueType             ValueType;
       typedef DANTsR_GetDataValue<VectorType>            GetValueType;

       bool operator!=(const DANTsR_Exp &) const
       {
         return false;
       }

      bool operator==(const DANTsR_Exp & other) const
       {
         return !( *this != other );
       }

       inline DataType operator()(const VectorType & A ) const
        {
          using TensorType = typename itk::DiffusionTensor3D< typename VectorType::ValueType>;
          TensorType dt( A.GetDataPointer() );
          bool junk=false;
          TensorType dt2 = TensorExp<TensorType>( dt, junk );
          VectorType logDt( dt2.GetDataPointer(),6 );
          return( logDt );

         /*
         using MatrixType = typename itk::Matrix<typename DataType::ValueType, 3,3>;
         using ArrayType = typename itk::Vector<typename DataType::ValueType, 3>;
         using SolverType = itk::SymmetricEigenAnalysis<MatrixType,ArrayType,MatrixType>;
         SolverType solver;
         ArrayType eigenValues;
         MatrixType eigenVectors;

         MatrixType dtMat = tensorAsMatrix<VectorType, MatrixType>(A);

         solver.SetDimension(3);
         solver.ComputeEigenValuesAndVectors(dtMat, eigenValues, eigenVectors);

         MatrixType eValues;
         eValues.Fill(0);

         for (unsigned int i=0; i<3; i++)
           {
           eValues(i,i) = std::exp(eigenValues[i]);
           }

         MatrixType DTrec = eigenVectors * eValues * eigenVectors.GetTranspose();
         DataType dtv = tensorAsVector<MatrixType,DataType>(DTrec);
         return(zeroNonFiniteTensors<DataType>(dtv));
         */
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
  else if ( filter == "log" )
  {
    typedef itk::UnaryFunctorImageFilter<VectorImageType, VectorImageType,
      itk::Functor::DANTsR_Log<VectorType,VectorType>  > FilterType;

    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput( image );
    filter->DebugOn();
    filter->Update();
    typename VectorImageType::Pointer outImage = filter->GetOutput();
    outImage->DisconnectPipeline();

    return Rcpp::wrap(outImage);
  }
  else if ( filter == "exp" )
  {
    typedef itk::UnaryFunctorImageFilter<VectorImageType, VectorImageType,
      itk::Functor::DANTsR_Exp<VectorType,VectorType>  > FilterType;

    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput( image );
    filter->Update();
    typename VectorImageType::Pointer outImage = filter->GetOutput();
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
