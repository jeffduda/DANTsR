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
#ifndef itkLabeledImageToPointSetFilter_hxx
#define itkLabeledImageToPointSetFilter_hxx

#include "itkLabeledImageToPointSetFilter.h"
#include "itkNumericTraits.h"
#include "itkProgressReporter.h"

namespace itk
{
/**
 *
 */
template< typename TInputImage, typename TOutputMesh >
LabeledImageToPointSetFilter< TInputImage, TOutputMesh >
::LabeledImageToPointSetFilter()
{
  // Modify superclass default values, can be overridden by subclasses
  this->SetNumberOfRequiredInputs(1);

  m_Label = 1;
  OutputMeshPointer mesh = this->GetOutput();

}

/**
 *
 */
template< typename TInputImage, typename TOutputMesh >
LabeledImageToPointSetFilter< TInputImage, TOutputMesh >
::~LabeledImageToPointSetFilter()
{}

/**
 *
 */
template< typename TInputImage, typename TOutputMesh >
void
LabeledImageToPointSetFilter< TInputImage, TOutputMesh >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Label: " << m_Label << std::endl;
}

/**
 *
 */
template< typename TInputImage, typename TOutputMesh >
void
LabeledImageToPointSetFilter< TInputImage, TOutputMesh >
::GenerateOutputInformation()
{}

/**
 *
 */
template< typename TInputImage, typename TOutputMesh >
void
LabeledImageToPointSetFilter< TInputImage, TOutputMesh >
::SetInput(const InputImageType *inputImage)
{
  // This const_cast is needed due to the lack of
  // const-correctness in the ProcessObject.
  this->SetNthInput( 0,
                     const_cast< InputImageType * >( inputImage ) );
}

/**
 *
 */
template< typename TInputImage, typename TOutputMesh >
void
LabeledImageToPointSetFilter< TInputImage, TOutputMesh >
::GenerateData(void)
{

  OutputMeshPointer      mesh      = this->GetOutput();
  InputImageConstPointer image     = this->GetInput(0);
  PointsContainerPointer points    = PointsContainer::New();

  InputImageIterator it( image, image->GetRequestedRegion() );
  PointType point;

  while ( !it.IsAtEnd() )
    {
    InputImagePixelType val = it.Value();

    if ( val ==  m_Label )
      {
      image->TransformIndexToPhysicalPoint(it.GetIndex(), point);
      points->push_back(point);
      }
    ++it;
    //progress.CompletedPixel();
    }

  mesh->SetPoints(points);

  // This indicates that the current BufferedRegion is equal to the
  // requested region. This action prevents useless rexecutions of
  // the pipeline.
  mesh->SetBufferedRegion( mesh->GetRequestedRegion() );
}
} // end namespace itk

#endif
