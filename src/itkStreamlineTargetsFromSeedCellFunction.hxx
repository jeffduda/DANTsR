#ifndef itkStreamlineTargetsFromSeedCellFunction_hxx
#define itkStreamlineTargetsFromSeedCellFunction_hxx

#include "itkStreamlineTargetsFromSeedCellFunction.h"
#include "itkMath.h"

namespace itk
{
template< typename TInputMesh, typename TInputImage >
StreamlineTargetsFromSeedCellFunction< TInputMesh, TInputImage >
::StreamlineTargetsFromSeedCellFunction()
{
  m_Image = nullptr;
}

template< typename TInputMesh, typename TInputImage >
itk::FixedArray<unsigned long, 2>
StreamlineTargetsFromSeedCellFunction< TInputMesh, TInputImage >
::Evaluate(const IndexType & index) const
{
  if ( m_Image == nullptr ) {
    itkExceptionMacro("Mask image not defined");
  }

  CellAutoPointer cell;
  itk::FixedArray<unsigned long, 2> hits;
  hits[0]=0;
  hits[1]=0;

  if ( !this->GetInputMesh()->GetCell(index, cell) ) {
    hits[0] = std::numeric_limits<unsigned long>::quiet_NaN();
    hits[1] = hits[0];
    return(hits);
  }

  unsigned long seed = m_Seeds[index];
  unsigned long position = 0;

  CellPointIterator it = cell->PointIdsBegin();
  ImagePointType iPt;
  ImageIndexType idx;

  while ( (it != cell->PointIdsEnd()) && !hits[0] && !hits[1] ) {

    PointType pt = this->GetInputMesh()->GetPoint(*it);

    // FIXME - possible no need to copy here
    for (unsigned int i=0; i<InputImageType::ImageDimension; i++) {
      iPt[i] = pt[i];
    }

    m_Image->TransformPhysicalPointToIndex( iPt, idx );
    if (m_Image->GetPixel(idx) > 0) {

      if ( position <= seed ) {
        hits[0] = m_Image->GetPixel(idx);
      }
      else if ( !hits[1] ){
        hits[1] = m_Image->GetPixel(idx);
      }

    }

    ++it;
    ++position;
  }

  return hits;
}

template< typename TInputMesh, typename TInputImage >
void
StreamlineTargetsFromSeedCellFunction< TInputMesh, TInputImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}
} // end namespace itk

#endif
