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
#ifndef itkDeterministicDTITractography_hxx
#define itkDeterministicDTITractography_hxx

#include "itkDeterministicDTITractography.h"
#include "itkContinuousIndex.h"
#include "itkNumericTraits.h"
#include "itkMath.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

namespace itk
{
template< typename TInputImage, typename TOutputMesh >
DeterministicDTITractography< TInputImage, TOutputMesh >
::DeterministicDTITractography()
{
  // Modify superclass default values, can be overridden by subclasses
  this->SetNumberOfRequiredInputs(2);

  //this->GetOutput()->GetPoints()->Reserve(3000);
  m_Mask = nullptr;

  m_MinimumNumberOfPoints = 2;
  m_MaximumNumberOfPoints = itk::NumericTraits<unsigned long>::max();
  m_StepSize = 0.2;

  //this->GetOutput()->GetCells()->Reserve(m_CellLimit);
}

template< typename TInputImage, typename TOutputMesh >
DeterministicDTITractography< TInputImage, TOutputMesh >
::~DeterministicDTITractography()
{

}

template< typename TInputImage, typename TOutputMesh >
void
DeterministicDTITractography< TInputImage, TOutputMesh >
::SetField(const InputImageType *image)
{
  this->SetNthInput( 0, const_cast< TInputImage * >( image ) );
}

template< typename TInputImage, typename TOutputMesh >
typename TInputImage::ConstPointer
DeterministicDTITractography< TInputImage, TOutputMesh >
::GetField()
{
  return static_cast< const TInputImage * >
         ( this->ProcessObject::GetInput(0) );
}

template< typename TInputImage, typename TOutputMesh >
void
DeterministicDTITractography< TInputImage, TOutputMesh >
::SetSeeds(const TOutputMesh * seeds)
{
  this->ProcessObject::SetNthInput( 1, const_cast< TOutputMesh * >( seeds ) );
}

template< typename TInputImage, typename TOutputMesh >
typename TOutputMesh::ConstPointer
DeterministicDTITractography< TInputImage, TOutputMesh >
::GetSeeds()
{
  return static_cast< const TOutputMesh * >
         ( this->ProcessObject::GetInput(1) );
}



template< typename TInputImage, typename TOutputMesh >
typename DeterministicDTITractography< TInputImage, TOutputMesh >::PointsContainerPointer
DeterministicDTITractography< TInputImage, TOutputMesh >
::TrackFiber( PointType point, bool forward )
{
  using InterpolatorType = itk::LinearInterpolateImageFunction< InputImageType, ValueType>;
  using IterpolatorPointer = typename InterpolatorType::Pointer;
  using IndexType = typename itk::ContinuousIndex<ValueType, InputImageType::ImageDimension>;

  using NNInterpolatorType = itk::NearestNeighborInterpolateImageFunction< MaskType, ValueType>;
  using NNInterpolatorPointer = typename NNInterpolatorType::Pointer;
  NNInterpolatorPointer maskinterp = NNInterpolatorType::New();
  if ( m_Mask != nullptr ) {
    maskinterp->SetInputImage( m_Mask );
  }

  PointsContainerPointer points = PointsContainer::New();

  IterpolatorPointer interp = InterpolatorType::New();
  interp->SetInputImage( this->GetField() );

  IndexType idx;
  this->GetField()->TransformPhysicalPointToContinuousIndex( point, idx );
  InputPixelType iDirection = interp->EvaluateAtContinuousIndex(idx);
  if ( !forward ) {
    iDirection *= -1;
  }

  InputPixelType previousDirection = iDirection;
  PointType currentPoint = point;
  PointType previousPoint = point;

  points->InsertElement(0, point);

  ValueType stepSize = 0.2;

  bool continueTracking = true;

  unsigned long id = 1;

  while( continueTracking ) {

    this->GetField()->TransformPhysicalPointToContinuousIndex( currentPoint, idx );
    iDirection = interp->EvaluateAtContinuousIndex(idx);

    ValueType dirCheck = iDirection[0]*previousDirection[0] + iDirection[1]*previousDirection[1] + iDirection[2]*previousDirection[2];

    if ( dirCheck < 0 ) {
      iDirection *= -1;
    }
    VectorType vec;
    for (unsigned i=0; i<InputImageType::ImageDimension; i++) {
      vec[i] = iDirection[i];
    }
    vec /= vec.GetNorm();
    vec *= m_StepSize;


    //std::cout << id << " : " << currentPoint[0] << "," << currentPoint[1] << "," << currentPoint[2] << std::endl;

    previousPoint = currentPoint;
    currentPoint = previousPoint + vec;
    previousDirection = iDirection;

    bool addPoint = true;

    if ( m_Mask != nullptr ) {
      if ( maskinterp->EvaluateAtContinuousIndex(idx) < 1 ) {
        addPoint = false;
      }
    }

    if ( addPoint ) {
      points->InsertElement(id, currentPoint);
      id++;
    }
    else {
      continueTracking = false;
    }

    if ( id > 1000 ) {
      continueTracking = false;
    }

  }


  return(points);
}



/** Generate the data */
template< typename TInputImage, typename TOutputMesh >
void
DeterministicDTITractography< TInputImage, TOutputMesh >
::GenerateData()
{

  // FIXME - make this parallel as each seed is independent

  std::cout << "GenerateData()" << std::endl;
  // Initialize variables


  m_OutputMesh = this->GetOutput();

  m_InputImage = this->GetField();
  m_SeedMesh = this->GetSeeds();

  std::cout << "Tracking from " << m_SeedMesh->GetNumberOfPoints() << " seed points" << std::endl;
  //this->m_OutputMesh->GetPoints()->Reserve(1);

  IdentifierType nPoints = 0;


  for ( IdentifierType i=0; i<m_SeedMesh->GetNumberOfPoints(); i++ )
  //for ( IdentifierType i=0; i<1000; i++ )
  {

    //if ( !m_OutputMesh->GetPoints()->IndexExists( nPoints ) )
    //{
    //  m_OutputMesh->GetPoints()->CreateIndex( nPoints ); // FIXME - add memory in chunks?
    //}

    PointsContainerPointer forwardPts = this->TrackFiber( m_SeedMesh->GetPoints()->GetElement(i), true );
    PointsContainerPointer backwardPts = this->TrackFiber( m_SeedMesh->GetPoints()->GetElement(i), false );

    unsigned long tractPoints = forwardPts->Size() + backwardPts->Size() - 1;
    std::cout << "Tract " << i << " has " << tractPoints << " points" << std::endl;

    bool keepTract = true;
    if ( tractPoints < m_MinimumNumberOfPoints ) {
      keepTract = false;
      std::cout << "Too few points" << std::endl;
    }
    else if ( tractPoints > m_MaximumNumberOfPoints ) {
      keepTract = false;
      std::cout << "too many points" << std::endl;
    }

    if ( keepTract ) {
      std::cout << "Add tract" << std::endl;
      typename OutputMeshType::PointIdentifier polyPoints[ tractPoints ];

      unsigned long tractId=0;
      for ( unsigned long i=(backwardPts->Size()-1); i>0; i-- ) {
        m_OutputMesh->GetPoints()->InsertElement( nPoints, backwardPts->GetElement(i) );
        polyPoints[tractId] = nPoints;
        ++tractId;
        ++nPoints;
      }
      for ( unsigned long i=0; i<forwardPts->Size(); i++ )  {
        m_OutputMesh->GetPoints()->InsertElement( nPoints, forwardPts->GetElement(i) );
        polyPoints[tractId] = nPoints;
        ++tractId;
        ++nPoints;
      }

      PolyLineCellType * polyline = new PolyLineCellType;
      polyline->SetPointIds( 0, tractPoints, polyPoints );
      CellAutoPointer streamline;
      streamline.TakeOwnership( polyline );
      m_OutputMesh->SetCell(i, streamline);
    }

  }


  // This indicates that the current BufferedRegion is equal to the
  // requested region. This action prevents useless rexecutions of
  // the pipeline.
  //this->m_OutputMesh->SetBufferedRegion( this->GetOutput()->GetRequestedRegion() );
  std::cout << "End GenerateData()" << std::endl;
}

/** PrintSelf */
template< typename TInputImage, typename TOutputMesh >
void
DeterministicDTITractography< TInputImage, TOutputMesh >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent
     << "NumberOfNodes: "
     << m_NumberOfNodes
     << std::endl;

  os << indent
     << "NumberOfCells: "
     << m_NumberOfCells
     << std::endl;

}
} /** end namespace itk. */

#endif
