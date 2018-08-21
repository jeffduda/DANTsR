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
/*=========================================================================
 *
 *  Portions of this file are subject to the VTK Toolkit Version 3 copyright.
 *
 *  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
 *
 *  For complete copyright, license and disclaimer of warranty information
 *  please refer to the NOTICE file at the top of the ITK source tree.
 *
 *=========================================================================*/
#ifndef itkMeshExtractionFilter_hxx
#define itkMeshExtractionFilter_hxx

#include "itkMeshExtractionFilter.h"
#include "itkNumericTraits.h"
#include "itkObjectFactory.h"

#include <set>

namespace itk
{
/**
 * ------------------------------------------------
 */
template< typename TInputMesh, typename TOutputMesh >
MeshExtractionFilter< TInputMesh, TOutputMesh >
::MeshExtractionFilter() :
  m_NumberOfCellsInRegion(NumericTraits< SizeValueType >::ZeroValue()),
  m_RegionNumber(NumericTraits< IdentifierType >::ZeroValue())
{
}



/**
 * ------------------------------------------------
 */
template< typename TInputMesh, typename TOutputMesh >
void
MeshExtractionFilter< TInputMesh, TOutputMesh >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Extraction Mode: ";
  if ( m_ExtractionMode == Self::PointRegions )
    {
    os << "Point Seeded Regions" << std::endl;
    }
  else if ( m_ExtractionMode == Self::CellRegions )
    {
    os << "Cell Seeded Regions" << std::endl;
    }

}

/**
 *
 */
template< typename TInputMesh, typename TOutputMesh >
void
MeshExtractionFilter< TInputMesh, TOutputMesh >
::GenerateData()
{
  std::cout << "GenerateData()" << std::endl;

  InputMeshConstPointer                  input = this->GetInput();
  OutputMeshPointer                      output = this->GetOutput();
  InputMeshPointsContainerConstPointer   inPts = input->GetPoints();
  InputMeshCellsContainerConstPointer    inCells = input->GetCells();
  InputMeshCellDataContainerConstPointer inCellData = input->GetCellData();

  itkDebugMacro(<< "Executing connectivity");

  //  Check input/allocate storage
  IdentifierType numCells = input->GetNumberOfCells();
  IdentifierType numPts = input->GetNumberOfPoints();

  std::cout << "Input mesh has " << numPts << " points and " << numCells << " cells " << std::endl;

  std::cout << "Number of points for extraction: " << m_PointList.size() << std::endl;
  std::cout << "Number of cells for extraction: " << m_CellList.size() << std::endl;

  // Initialize.  Keep track of points and cells visited.
  //input->BuildCellLinks(); //needed to get neighbors

  std::vector< IdentifierType > cellVisited(numCells);

  // Add all cells that were passed in
  for (unsigned long i=0; i<m_CellList.size(); i++ )
    {
    cellVisited[ m_CellList[i] ] = 1;
    }

  std::cout << "Found cells passed in" << std::endl;

  // Add all cells that contain a passed point
  InputMeshCellAutoPointer cell;
  for ( unsigned long i=0; i < m_PointList.size(); i++ )
    {
    std::cout << "checking point " << i << std::endl;
    for ( unsigned long j=0; j < numCells; j++ )
      {
        this->GetInput()->GetCell(j, cell);
        for ( unsigned long k=0; k < cell->GetNumberOfPoints(); k++ )
          {
          if ( cell->GetPointIds()[k] == m_PointList[i] )
            {
            cellVisited[j] = 1;
            std::cout << "point " << i << " contained in cell " << j << std::endl;
            }
          }
      }
    }

  std::cout << "Found cells that include passed points" << std::endl;
  for (unsigned long i=0; i<numCells; i++ )
    {
    if ( cellVisited[i] )
      {
      std::cout << "cell " << i << " included" << std::endl;
      }
    }



  std::vector< IdentifierType > pointVisited(numPts);

  // points that were passed in
  for ( unsigned long i=0; i<m_PointList.size(); i++)
    {
    std::cout << "Passed point " << m_PointList[i] << std::endl;
    pointVisited[ m_PointList[i] ] = 1;
    }

  for ( unsigned long i=0; i<numPts; i++)
    {
    std::cout << i << ":" << pointVisited[i] << std::endl;
    }

  std::cout << "Found points passed in" << std::endl;

  // add all points in included cells
  for (unsigned long i=0; i<numCells; i++ )
    {
    if ( cellVisited[i] )
      {
      std::cout << "add points for cell " << i << std::endl;
      this->GetInput()->GetCell(i, cell);
      for ( unsigned long j=0; j<cell->GetNumberOfPoints(); j++)
        {
        pointVisited[ cell->GetPointIds()[j] ] = 1;
        std::cout << "cell " << i << " adds point " << cell->GetPointIds()[j] << std::endl;
        }
      }
    }

  std::cout << "Found points in all included cells" << std::endl;
  std::vector< IdentifierType > pointMap;
  for ( unsigned long i=0; i<pointVisited.size(); i++ )
    {
    if ( pointVisited[i] )
      {
      std::cout << "Retaining point " << i << std::endl;
      pointMap.push_back(i);
      }
    }

  std::vector< IdentifierType > cellMap;
  for ( unsigned int i=0; i<cellVisited.size(); i++ )
    {
    if ( cellVisited[i] )
      {
      cellMap.push_back(i);
      }
    }

  // Copy points to output
  std::cout << "Output mesh will have " << pointMap.size() << " points" << std::endl;
  output->GetPoints()->Reserve( pointMap.size() );
  for (unsigned int i=0; i<pointMap.size(); i++)
    {
    output->GetPoints()->InsertElement( i, input->GetPoint( pointMap[i] ) );
    pointVisited[ pointMap[i] ] = i; // reverse of pointMap
    }

  std::cout << "Set output mesh points" << std::endl;

  std::cout << "Output mesh will have " << cellMap.size() << " cells" << std::endl;

  // Copy cells to output mesh
  OutputMeshCellAutoPointer ocell;
  for (unsigned long i=0; i<cellMap.size(); i++)
    {
    // get the cell
    //this->GetInput()->GetCell(cellMap[i], cell);
    //InputMeshCellType c = this->GetInput()->GetCells()[ cellMap[i] ];

    this->GetInput()->GetCell( cellMap[i], cell );
    std::cout << "Add cell " << cellMap[i] << " with " << cell->GetNumberOfPoints() << " points " << std::endl;

    if ( cell->GetType() == InputMeshCellType::LINE_CELL )
      {
      ocell.TakeOwnership( new LineType );
      }
    else if ( cell->GetType() == InputMeshCellType::POLYGON_CELL )
      {
      ocell.TakeOwnership( new PolyType );
      }

    for ( unsigned long j=0; j<cell->GetNumberOfPoints(); j++)
      {
      ocell->SetPointId( j, pointVisited[ cell->GetPointIds()[j] ] );
      }

    output->SetCell(i, ocell);

    }

  // Copy point data

  // Copy cell data



  // Report some statistics
  if ( this->GetDebug() )
    {
    itkDebugMacro (<< "Extracted " << output->GetNumberOfPoints() << " points");
    itkDebugMacro (<< "Extracted " << output->GetNumberOfCells() << " cells");
    }

  // This prevents unnecessary re-executions of the pipeline.
  output->SetBufferedRegion( output->GetRequestedRegion() );
}

} // end namespace itk

#endif
