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
#ifndef itkPolyLineMeshBSplineSmoothingFilter_hxx
#define itkPolyLineMeshBSplineSmoothingFilter_hxx

#include "itkPolyLineMeshBSplineSmoothingFilter.h"
#include "itkMacro.h"

namespace itk
{
/**
 *
 */
template< typename TInputMesh, typename TOutputMesh >
PolyLineMeshBSplineSmoothingFilter< TInputMesh, TOutputMesh >
::PolyLineMeshBSplineSmoothingFilter()
{
  this->m_UsePointDataAsParameter = false;
}

/**
 *
 */
template< typename TInputMesh, typename TOutputMesh >
void
PolyLineMeshBSplineSmoothingFilter< TInputMesh, TOutputMesh >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

/**
 * This method causes the filter to generate its output.
 */
template< typename TInputMesh, typename TOutputMesh >
void
PolyLineMeshBSplineSmoothingFilter< TInputMesh, TOutputMesh >
::GenerateData(void)
{
  // Hard code for now as these parameters tend to work well for this application (i.e. fiber tractography)
  unsigned int SplineOrder = 8;
  unsigned int NumberOfControlPoints = SplineOrder+2;
  unsigned int NumberOfLevels = 4;
  unsigned int nPoints = 100;
  unsigned int nIterations = 1;
  double convergence_thresh = 0.001;


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
    std::cout  << "Smoothing cell " << id << " of " << nCells << std::endl;

    CellAutoPointer cell;
    if ( this->GetInput()->GetCell(id, cell) )
      {

      ParametricMeshPointer splineMesh = ParametricMeshType::New();
      splineMesh->GetPoints()->Initialize();
      splineMesh->GetPointData()->Initialize();

      // Get full length of line if needed
      CoordRepType cellLength=0.0;
      if ( ! this->m_UsePointDataAsParameter )
        {
        for (unsigned int i=0; i<cell->GetNumberOfPoints(); i++ )
          {
            InputPointType pt = this->GetInput()->GetPoint( cell->GetPointIds()[i] );
            if ( i > 0 )
              {
              InputPointType lastPt = this->GetInput()->GetPoint( cell->GetPointIds()[i-1] );
              cellLength += lastPt.EuclideanDistanceTo(pt);
              }
          }
        }

        // Set param value for each point
        CoordRepType length = 0;
        for (unsigned int i=0; i<cell->GetNumberOfPoints(); i++ )
          {

          InputPointType pt = this->GetInput()->GetPoint( cell->GetPointIds()[i] );
          if ( ! this->m_UsePointDataAsParameter )
            {
            if ( i > 0 )
              {
              InputPointType lastPt = this->GetInput()->GetPoint( cell->GetPointIds()[i-1] );
              length += lastPt.EuclideanDistanceTo(pt);
              }
            splineMesh->SetPoint(i, length/cellLength);
            }
          else
            {
            splineMesh->SetPoint(i, this->GetInput()->GetPointData()->GetElement( cell->GetPointIds()[i]  ));
            }

          ParametricPixelType splineCoord;
          for ( unsigned int d=0; d<InputMeshType::PointDimension; d++)
            {
            splineCoord[d] = pt[d];
            }

          //splineMesh->SetPoint(i, splineParam);
          splineMesh->SetPointData(i, splineCoord);
          }

        typename ParametricImageType::SizeType psize;
        typename ParametricImageType::PointType porigin;
        typename ParametricImageType::SpacingType pspacing;
        typename BSplineFilterType::ArrayType ncps;


        psize.Fill( cell->GetNumberOfPoints() );
        porigin.Fill( 0.0 );
        pspacing.Fill( 1.0 / (cell->GetNumberOfPoints()-1) );
        ncps[0] = NumberOfControlPoints;

        typename BSplineFilterType::Pointer smoother = BSplineFilterType::New();
        smoother->SetSize( psize );
        smoother->SetSpacing( pspacing );
        smoother->SetOrigin( porigin );
        smoother->SetSplineOrder( SplineOrder );
        smoother->SetNumberOfControlPoints(ncps);
        smoother->SetNumberOfLevels(NumberOfLevels);

        smoother->SetInput( splineMesh );
        smoother->Update();

        ParametricImagePointer img = smoother->GetOutput();
        typename ParametricImageType::IndexType idx;
        for (unsigned int i=0; i<cell->GetNumberOfPoints(); i++)
          {
          idx[0] = i;
          ParametricPixelType coord = img->GetPixel(idx);
          OutputPointType pt;
          for ( unsigned int d=0; d<OutputMeshType::PointDimension; d++)
            {
            pt[d] = coord[d];
            }

          outPoints->InsertElement( cell->GetPointIds()[i], pt);
          }

      }
    }

  outputMesh->SetPoints( outPoints );

  // Create duplicate references to the rest of data on the mesh
  this->CopyInputMeshToOutputMeshCellLinks();
  this->CopyInputMeshToOutputMeshCells();
  this->CopyInputMeshToOutputMeshCellData();

  unsigned int maxDimension = TInputMesh::MaxTopologicalDimension;
  for ( unsigned int dim = 0; dim < maxDimension; dim++ )
    {
    outputMesh->SetBoundaryAssignments( dim,
                                        inputMesh->GetBoundaryAssignments(dim) );
    }

}
} // end namespace itk

#endif
