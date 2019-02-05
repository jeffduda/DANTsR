#ifndef itkCellImageValueSummaryCellFunction_hxx
#define itkCellImageValueSummaryCellFunction_hxx

#include "itkCellImageValueSummaryCellFunction.h"
#include "itkMath.h"

namespace itk
{
template< typename TInputMesh, typename TInputImage, typename TOutput >
CellImageValueSummaryCellFunction< TInputMesh, TInputImage, TOutput >
::CellImageValueSummaryCellFunction()
{
  m_Image = nullptr;
  m_Measure = CellImageValueSummaryCellFunction< TInputMesh, TInputImage, TOutput >::MEAN;
}

template< typename TInputMesh, typename TInputImage, typename TOutput >
TOutput
CellImageValueSummaryCellFunction< TInputMesh, TInputImage, TOutput >
::Evaluate(const IndexType & index) const
{
  if ( m_Image == nullptr ) {
    itkExceptionMacro("Image not defined");
  }

  CellAutoPointer cell;

  if ( !this->GetInputMesh()->GetCell(index, cell) ) {
    return std::numeric_limits<TOutput>::quiet_NaN();
  }

  std::vector< TOutput > values;

  CellPointIterator it = cell->PointIdsBegin();
  ImagePointType iPt;
  ImageIndexType idx;

  while ( it != cell->PointIdsEnd() ) {

    PointType pt = this->GetInputMesh()->GetPoint(*it);

    // FIXME - possible no need to copy here
    for (unsigned int i=0; i<InputImageType::ImageDimension; i++) {
      iPt[i] = pt[i];
    }

    m_Image->TransformPhysicalPointToIndex( iPt, idx );
    values.push_back( static_cast<TOutput>(m_Image->GetPixel(idx)) );

    ++it;
  }

  TOutput summary=0;


  if ( m_Measure == MEAN ) {
    for (TOutput v : values) {
      summary += v;
    }
    summary /= values.size();
  }
  else if ( m_Measure == MEDIAN ) {
    size_t size = values.size();
    sort( values.begin(), values.end() );
    if (size % 2 == 0)
    {
      summary = (values[size / 2 - 1] + values[size / 2]) / 2;
    }
    else
    {
      summary = values[size / 2];
    }
  }
  else if ( m_Measure == MAX ) {
    for (TOutput v : values) {
      if ( v > summary ) {summary=v;}
    }
  }
  else if ( m_Measure == MIN ) {
    summary = std::numeric_limits<TOutput>::max();
    for (TOutput v : values) {
      if ( v < summary ) {summary=v;}
    }
  }

  return summary;
}

template< typename TInputMesh, typename TInputImage, typename TOutput >
void
CellImageValueSummaryCellFunction< TInputMesh, TInputImage, TOutput >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}
} // end namespace itk

#endif
