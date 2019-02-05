#ifndef itkCellPointInMaskCellFunction_hxx
#define itkCellPointInMaskCellFunction_hxx

#include "itkCellPointInMaskCellFunction.h"
#include "itkMath.h"

namespace itk
{
template< typename TInputMesh, typename TInputImage >
CellPointInMaskCellFunction< TInputMesh, TInputImage >
::CellPointInMaskCellFunction()
{
  m_Mask = nullptr;
}

template< typename TInputMesh, typename TInputImage >
bool
CellPointInMaskCellFunction< TInputMesh, TInputImage >
::Evaluate(const IndexType & index) const
{
  if ( m_Mask == nullptr ) {
    itkExceptionMacro("Mask image not defined");
  }

  CellAutoPointer cell;
  CellAutoPointer edge;

  if ( !this->GetInputMesh()->GetCell(index, cell) ) {
    return false;
  }

  bool hit=false;

  CellPointIterator it = cell->PointIdsBegin();
  ImagePointType iPt;
  ImageIndexType idx;

  while ( (it != cell->PointIdsEnd()) && !hit ) {

    PointType pt = this->GetInputMesh()->GetPoint(*it);

    // FIXME - possible no need to copy here
    for (unsigned int i=0; i<InputImageType::ImageDimension; i++) {
      iPt[i] = pt[i];
    }

    m_Mask->TransformPhysicalPointToIndex( iPt, idx );
    if (m_Mask->GetPixel(idx) > 0) {
      hit = true;
    }

    ++it;
  }

  return hit;
}

template< typename TInputMesh, typename TInputImage >
void
CellPointInMaskCellFunction< TInputMesh, TInputImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}
} // end namespace itk

#endif
