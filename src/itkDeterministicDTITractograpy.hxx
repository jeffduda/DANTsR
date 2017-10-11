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

namespace itk
{
template< typename TInputImage, typename TOutputMesh, typename TMaskImage >
DeterministicDTITractography< TInputImage, TOutputMesh, TMaskImage >
::DeterministicDTITractography() :
  m_SeedMesh(ITK_NULLPTR),
  m_OutputMesh(ITK_NULLPTR),
  m_InputImage(ITK_NULLPTR),
  m_ValidRegion(ITK_NULLPTR)
{
  // Modify superclass default values, can be overridden by subclasses
  this->SetNumberOfRequiredInputs(1);

  this->GetOutput()->GetPoints()->Reserve(3000);
  //this->GetOutput()->GetCells()->Reserve(m_CellLimit);
}

template< typename TInputImage, typename TOutputMesh, typename TMaskImage >
DeterministicDTITractography< TInputImage, TOutputMesh, TMaskImage >
::~DeterministicDTITractography()
{

}

template< typename TInputImage, typename TOutputMesh, typename TMaskImage >
void
DeterministicDTITractography< TInputImage, TOutputMesh, TMaskImage >
::SetInput(const InputImageType *image)
{
  this->ProcessObject::SetNthInput( 0,
                                    const_cast< InputImageType * >( image ) );
}

template< typename TInputImage, typename TOutputMesh, typename TMaskImage >
void
DeterministicDTITractography< TInputImage, TOutputMesh, TMaskImage >
::SetSeeds(const OutputMeshType *seeds)
{
  this->ProcessObject::SetNthInput( 1,
                                    const_cast< OutputMeshType * >( seeds ) );
}

/** Generate the data */
template< typename TInputImage, typename TOutputMesh, typename TMaskImage >
void
DeterministicDTITractography< TInputImage, TOutputMesh, TMaskImage >
::GenerateData()
{

  // FIXME - make this parallel as each seed is independent

  std::cout << "GenerateData()" << std::endl;
  // Initialize variables
  m_OutputMesh = this->GetOutput();
  m_InputImage =
    static_cast< const InputImageType * >( this->ProcessObject::GetInput(0) );
  m_SeedMesh =
    static_cast< const OutputMeshType * >( this->ProcessObject::GetInput(1) );



  std::cout << "Tracking from " << m_SeedMesh->GetNumberOfPoints() << " seed points" << std::endl;
  this->m_OutputMesh->GetPoints()->Reserve()

  IdentifierType nPoints = 0;

  typename OutputMeshType::PointsContainer points = this->m_OutputMesh->GetPoints();

  for ( IdentifierType i=0; i<m_SeedMesh->GetNumberOfPoints(); i++ )
  {
    if ( !points->IndexExists( nPoints ) )
    {
      points->CreateIndex( nPoints ); // FIXME - add memory in chunks?
    }
    points->SetElement( nPoints, m_SeedMesh->GetPoints()->GetElement(i) );
    ++nPoints;
  }


  // This indicates that the current BufferedRegion is equal to the
  // requested region. This action prevents useless rexecutions of
  // the pipeline.
  this->m_OutputMesh->SetBufferedRegion( this->GetOutput()->GetRequestedRegion() );

}

/** PrintSelf */
template< typename TInputImage, typename TOutputMesh, typename TMaskImage >
void
DeterministicDTITractography< TInputImage, TOutputMesh, TMaskImage >
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

  os << indent
     << "SeedMeshProvidedByUser: "
     << m_SeedMeshProvidedByUser
     << std::endl;
}
} /** end namespace itk. */

#endif
