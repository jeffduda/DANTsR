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
#ifndef itkSpatialPolyLineMeshToArcLengthMeshFilter_hxx
#define itkSpatialPolyLineMeshToArcLengthMeshFilter_hxx

#include "itkSpatialPolyLineMeshToArcLengthMeshFilter.h"
#include "itkMacro.h"

namespace itk
{
/**
 *
 */
template< typename TInputMesh, typename TOutputMesh >
SpatialPolyLineMeshToArcLengthMeshFilter< TInputMesh, TOutputMesh >
::SpatialPolyLineMeshToArcLengthMeshFilter()
{
}

/**
 *
 */
template< typename TInputMesh, typename TOutputMesh >
void
SpatialPolyLineMeshToArcLengthMeshFilter< TInputMesh, TOutputMesh >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

/**
 * This method causes the filter to generate its output.
 */
template< typename TInputMesh, typename TOutputMesh >
void
SpatialPolyLineMeshToArcLengthMeshFilter< TInputMesh, TOutputMesh >
::GenerateData(void)
{
  using InputPointsContainer = typename TInputMesh::PointsContainer;
  using OutputPointsContainer = typename TOutputMesh::PointsContainer;

  using InputPointsContainerConstPointer = typename TInputMesh::PointsContainerConstPointer;
  using OutputPointsContainerPointer = typename TOutputMesh::PointsContainerPointer;

  const InputMeshType *inputMesh   =  this->GetInput();
  OutputMeshPointer    outputMesh   =  this->GetOutput();

  if ( !inputMesh )
    {
    itkExceptionMacro(<< "Missing Input Mesh");
    }

  if ( !outputMesh )
    {
    itkExceptionMacro(<< "Missing Output Mesh");
    }


  outputMesh->SetBufferedRegion( outputMesh->GetRequestedRegion() );

  InputPointsContainerConstPointer inPoints  = inputMesh->GetPoints();
  OutputPointsContainerPointer     outPoints = outputMesh->GetPoints();

  outPoints->Reserve( inputMesh->GetNumberOfPoints() );
  outPoints->Squeeze();  // in case the previous mesh had
                         // allocated a larger memory

  unsigned long nCells = inputMesh->GetNumberOfCells();

  for (unsigned long id=0; id<nCells; id++)
    {
    CellAutoPointer cell;
    if ( this->GetInput()->GetCell(id, cell) )
      {
      CoordRepType length=0.0;
      for (unsigned int i=0; i<cell->GetNumberOfPoints(); i++ )
        {
        InputPointType pt = this->GetInput()->GetPoint( cell->GetPointIds()[i] );
        if ( i > 0 )
          {
          InputPointType lastPt = this->GetInput()->GetPoint( cell->GetPointIds()[i-1] );
          length += lastPt.EuclideanDistanceTo(pt);
          }

        OutputPixelType outPix = length;

        OutputPointType outPt;
        for ( unsigned int d=0; d<InputMeshType::PointDimension; d++)
          {
          outPt[d] = pt[d];
          }

        outputMesh->SetPoint(cell->GetPointIds()[i], outPt);
        outputMesh->SetPointData(cell->GetPointIds()[i], outPix);

        }
      }
    }

  // Create duplicate references to the rest of data on the mesh
  //this->CopyInputMeshToOutputMeshPointData();
  this->CopyInputMeshToOutputMeshCellLinks();
  this->CopyInputMeshToOutputMeshCells();
  this->CopyInputMeshToOutputMeshCellData();

  //unsigned int maxDimension = TInputMesh::MaxTopologicalDimension;
  //for ( unsigned int dim = 0; dim < maxDimension; dim++ )
  //  {
  //  outputMesh->SetBoundaryAssignments( dim,
  //                                      inputMesh->GetBoundaryAssignments(dim) );
  //  }

}
} // end namespace itk

#endif
