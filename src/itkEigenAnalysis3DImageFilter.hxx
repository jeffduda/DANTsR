/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkEigenAnalysis3DImageFilter_hxx
#define itkEigenAnalysis3DImageFilter_hxx

#include "itkEigenAnalysis3DImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkProgressReporter.h"

namespace itk
{
/**
 * Constructor
 */
template< typename TInputImage, typename TEigenValueImage, typename TEigenVectorImage >
EigenAnalysis3DImageFilter< TInputImage, TEigenValueImage, TEigenVectorImage >
::EigenAnalysis3DImageFilter()
{
  this->SetNumberOfRequiredInputs(1);
  this->SetNumberOfRequiredOutputs(6);
  this->SetNthOutput( 0, this->MakeOutput(0) );
  this->SetNthOutput( 1, this->MakeOutput(1) );
  this->SetNthOutput( 2, this->MakeOutput(2) );
  this->SetNthOutput( 3, this->MakeOutput(3) );
  this->SetNthOutput( 4, this->MakeOutput(4) );
  this->SetNthOutput( 5, this->MakeOutput(5) );
}

/**
 * Get the largest eigenvalue considering the sign
 */
template< typename TInputImage, typename TEigenValueImage, typename TEigenVectorImage >
typename EigenAnalysis3DImageFilter< TInputImage, TEigenValueImage, TEigenVectorImage >::EigenValueImageType *
EigenAnalysis3DImageFilter< TInputImage, TEigenValueImage, TEigenVectorImage >
::GetMaxEigenValue(void)
{
  return dynamic_cast< EigenValueImageType * >(
           this->ProcessObject::GetOutput(2) );
}

/**
 * Get the eigenvalue considering the sign
 */
template< typename TInputImage, typename TEigenValueImage, typename TEigenVectorImage >
typename EigenAnalysis3DImageFilter< TInputImage, TEigenValueImage, TEigenVectorImage >::EigenValueImageType *
EigenAnalysis3DImageFilter< TInputImage, TEigenValueImage, TEigenVectorImage >
::GetEigenValue(unsigned int i)
{
  if ( (i>0) && (i<4) ) {
    return dynamic_cast< EigenValueImageType * >(
             this->ProcessObject::GetOutput(i-1) );
  }
  return NULL;

}

/**
 * Get the eigenvector
 */
template< typename TInputImage, typename TEigenValueImage, typename TEigenVectorImage >
typename EigenAnalysis3DImageFilter< TInputImage, TEigenValueImage, TEigenVectorImage >::EigenVectorImageType *
EigenAnalysis3DImageFilter< TInputImage, TEigenValueImage, TEigenVectorImage >
::GetEigenVector(unsigned int i)
{
  EigenVectorImageType *eigenVector = dynamic_cast< EigenVectorImageType * >(
    this->ProcessObject::GetOutput(i+2) );

  if ( eigenVector )
    {
    return eigenVector;
    }
  else
    {
    itkWarningMacro(

      <<
      "EigenAnalysis3DImageFilter::GetMaxEigenVector(): dynamic_cast has failed. A reinterpret_cast is being attempted."
      << std::endl << "Type name is: "
      << typeid( *this->GetOutput(i+2) ).name() );
    return reinterpret_cast< EigenVectorImageType * >(
             this->ProcessObject::GetOutput(i+2) );
    }
}


/**
 * Get the smallest eigenvalue considering the sign
 */
template< typename TInputImage, typename TEigenValueImage, typename TEigenVectorImage >
typename EigenAnalysis3DImageFilter< TInputImage, TEigenValueImage, TEigenVectorImage >::EigenValueImageType *
EigenAnalysis3DImageFilter< TInputImage, TEigenValueImage, TEigenVectorImage >
::GetMinEigenValue(void)
{
  return dynamic_cast< EigenValueImageType * >(
           this->ProcessObject::GetOutput(0) );
}

/**
 * Get the eigenvector corresponding to the largest eigenvalue (considering the sign)
 */
template< typename TInputImage, typename TEigenValueImage, typename TEigenVectorImage >
typename EigenAnalysis3DImageFilter< TInputImage, TEigenValueImage, TEigenVectorImage >::EigenVectorImageType *
EigenAnalysis3DImageFilter< TInputImage, TEigenValueImage, TEigenVectorImage >
::GetMaxEigenVector(void)
{
  EigenVectorImageType *eigenVector = dynamic_cast< EigenVectorImageType * >(
    this->ProcessObject::GetOutput(5) );

  if ( eigenVector )
    {
    return eigenVector;
    }
  else
    {
    itkWarningMacro(

      <<
      "EigenAnalysis3DImageFilter::GetMaxEigenVector(): dynamic_cast has failed. A reinterpret_cast is being attempted."
      << std::endl << "Type name is: "
      << typeid( *this->GetOutput(5) ).name() );
    return reinterpret_cast< EigenVectorImageType * >(
             this->ProcessObject::GetOutput(5) );
    }
}

/**
 *   Make Ouput
 * \todo Verify that MakeOutput is createing the right type of objects
 *  this could be the cause of the reinterpret_cast bug in this class
 */
template< typename TInputImage, typename TEigenValueImage, typename TEigenVectorImage >
DataObject::Pointer
EigenAnalysis3DImageFilter< TInputImage, TEigenValueImage, TEigenVectorImage >
::MakeOutput(DataObjectPointerArraySizeType idx)
{
  DataObject::Pointer output;

  switch ( idx )
    {
    case 0:
      output = ( EigenValueImageType::New() ).GetPointer();
      break;
    case 1:
      output = ( EigenValueImageType::New() ).GetPointer();
      break;
    case 2:
      output = ( EigenValueImageType::New() ).GetPointer();
      break;
    case 3:
      output = ( EigenVectorImageType::New() ).GetPointer();
      break;
    case 4:
      output = ( EigenVectorImageType::New() ).GetPointer();
      break;
    case 5:
      output = ( EigenVectorImageType::New() ).GetPointer();
      break;
    }
  return output.GetPointer();
}

/* Eigen analysis not threadsafe? */
template< typename TInputImage, typename TEigenValueImage, typename TEigenVectorImage >
void
EigenAnalysis3DImageFilter< TInputImage, TEigenValueImage, TEigenVectorImage >
::GenerateData()
{
  typename TInputImage::ConstPointer inputPtr(
    dynamic_cast< const TInputImage * >(
      ( ProcessObject::GetInput(0) ) ) );

  EigenValueImagePointer outputPtr1 = dynamic_cast< EigenValueImageType * >(
           this->ProcessObject::GetOutput(0) );
  EigenValueImagePointer outputPtr2 = dynamic_cast< EigenValueImageType * >(
           this->ProcessObject::GetOutput(1) );
  EigenValueImagePointer outputPtr3 = dynamic_cast< EigenValueImageType * >(
           this->ProcessObject::GetOutput(2) );
  EigenVectorImagePointer outputPtr4 = dynamic_cast< EigenVectorImageType * >(
           this->ProcessObject::GetOutput(3) );
  EigenVectorImagePointer outputPtr5 = dynamic_cast< EigenVectorImageType * >(
           this->ProcessObject::GetOutput(4) );
  EigenVectorImagePointer outputPtr6 = dynamic_cast< EigenVectorImageType * >(
           this->ProcessObject::GetOutput(5) );

  // FIXME - need this?
  if ( !outputPtr4 )
    {
    itkWarningMacro(

      <<
      "EigenAnalysis3DImageFilter::GetMaxEigenVector(): dynamic_cast has failed. A reinterpret_cast is being attempted."
      << std::endl << "Type name is: "
      << typeid( *this->GetOutput(3) ).name() );
    outputPtr4 =  reinterpret_cast< EigenVectorImageType * >(
             this->ProcessObject::GetOutput(3) );
    }
  if ( !outputPtr5 )
    {
    itkWarningMacro(

      <<
      "EigenAnalysis3DImageFilter::GetMaxEigenVector(): dynamic_cast has failed. A reinterpret_cast is being attempted."
      << std::endl << "Type name is: "
      << typeid( *this->GetOutput(4) ).name() );
    outputPtr5 =  reinterpret_cast< EigenVectorImageType * >(
             this->ProcessObject::GetOutput(4) );
    }
  if ( !outputPtr6 )
    {
    itkWarningMacro(

      <<
      "EigenAnalysis3DImageFilter::GetMaxEigenVector(): dynamic_cast has failed. A reinterpret_cast is being attempted."
      << std::endl << "Type name is: "
      << typeid( *this->GetOutput(5) ).name() );
    outputPtr6 =  reinterpret_cast< EigenVectorImageType * >(
             this->ProcessObject::GetOutput(5) );
    }

  outputPtr4->SetVectorLength(3);
  outputPtr5->SetVectorLength(3);
  outputPtr6->SetVectorLength(3);

  outputPtr1->SetBufferedRegion( inputPtr->GetBufferedRegion() );
  outputPtr2->SetBufferedRegion( inputPtr->GetBufferedRegion() );
  outputPtr3->SetBufferedRegion( inputPtr->GetBufferedRegion() );
  outputPtr4->SetBufferedRegion( inputPtr->GetBufferedRegion() );
  outputPtr5->SetBufferedRegion( inputPtr->GetBufferedRegion() );
  outputPtr6->SetBufferedRegion( inputPtr->GetBufferedRegion() );

  outputPtr1->Allocate();
  outputPtr2->Allocate();
  outputPtr3->Allocate();
  outputPtr4->Allocate();
  outputPtr5->Allocate();
  outputPtr6->Allocate();

  EigenValueImageRegionType region = outputPtr1->GetRequestedRegion();

  ImageRegionConstIteratorWithIndex< TInputImage > inputIt(inputPtr, region);

  ImageRegionIteratorWithIndex< EigenValueImageType >  outputIt1(outputPtr1, region);
  ImageRegionIteratorWithIndex< EigenValueImageType >  outputIt2(outputPtr2, region);
  ImageRegionIteratorWithIndex< EigenValueImageType >  outputIt3(outputPtr3, region);
  ImageRegionIteratorWithIndex< EigenVectorImageType > outputIt4(outputPtr4, region);
  ImageRegionIteratorWithIndex< EigenVectorImageType > outputIt5(outputPtr5, region);
  ImageRegionIteratorWithIndex< EigenVectorImageType > outputIt6(outputPtr6, region);

  EigenVectorType nullVector;
  nullVector.Fill(0.0);

  // support progress methods/callbacks
  ProgressReporter progress( this, 0, region.GetNumberOfPixels() );

  inputIt.GoToBegin();

  outputIt1.GoToBegin();
  outputIt2.GoToBegin();
  outputIt3.GoToBegin();
  outputIt4.GoToBegin();
  outputIt5.GoToBegin();
  outputIt6.GoToBegin();

  EigenVectorType eigenVector;

  while ( !inputIt.IsAtEnd() )
    {
    TensorType dt( inputIt.Get().GetDataPointer() );

    EigenVectorImagePixelType v1;
    EigenVectorImagePixelType v2;
    EigenVectorImagePixelType v3;
    v1.SetSize(3);
    v2.SetSize(3);
    v3.SetSize(3);

    typename TensorType::EigenValuesArrayType eigenValues;
    typename TensorType::EigenVectorsMatrixType eigenVectors;

    dt.ComputeEigenAnalysis(eigenValues, eigenVectors);

    EigenValueImagePixelType e1 = eigenValues[2];
    EigenValueImagePixelType e2 = eigenValues[1];
    EigenValueImagePixelType e3 = eigenValues[0];

    for (unsigned int i=0; i<3; i++) {
      v1[i] = eigenVectors[2][i];
      v2[i] = eigenVectors[1][i];
      v3[i] = eigenVectors[0][i];
    }

    outputIt1.Set(e1);
    outputIt2.Set(e2);
    outputIt3.Set(e3);

    outputIt4.Set(v1);
    outputIt5.Set(v2);
    outputIt6.Set(v3);

    ++inputIt;

    ++outputIt1;
    ++outputIt2;
    ++outputIt3;
    ++outputIt4;
    ++outputIt5;
    ++outputIt6;

    progress.CompletedPixel();
    }
}
} // end namespace itk

#endif
