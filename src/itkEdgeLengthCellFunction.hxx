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
#ifndef itkEdgeLengthCellFunction_hxx
#define itkEdgeLengthCellFunction_hxx

#include "itkEdgeLengthCellFunction.h"
#include "itkMath.h"

namespace itk
{
template< typename TInputMesh, typename TOutput >
EdgeLengthCellFunction< TInputMesh, TOutput >
::EdgeLengthCellFunction()
{
}

template< typename TInputMesh, typename TOutput >
TOutput
EdgeLengthCellFunction< TInputMesh, TOutput >
::Evaluate(const IndexType & index) const
{
  CellAutoPointer cell;
  CellAutoPointer edge;

  if ( !this->GetInputMesh()->GetCell(index, cell) ) {
    return std::numeric_limits<TOutput>::quiet_NaN();
  }

  TOutput length = 0;
  unsigned int nEdges = cell->GetNumberOfBoundaryFeatures(1);

  for (unsigned int i=0; i<nEdges; i++) {
    cell->GetBoundaryFeature(1, i, edge);
    CellPointIterator it = edge->PointIdsBegin();
    PointType p1 = this->GetInputMesh()->GetPoint(*it);
    ++it;
    PointType p2 = this->GetInputMesh()->GetPoint(*it);

    typename PointType::RealType dist = p1.EuclideanDistanceTo(p2);
    length += static_cast<TOutput>(dist);
  }

  return length;
}

template< typename TInputMesh, typename TOutput >
void
EdgeLengthCellFunction< TInputMesh, TOutput >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}
} // end namespace itk

#endif
