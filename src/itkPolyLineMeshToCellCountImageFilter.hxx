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
#ifndef itkPolyLineMeshToCellCountImageFilter_hxx
#define itkPolyLineMeshToCellCountImageFilter_hxx

#include "itkPolyLineMeshToCellCountImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkContinuousIndex.h"
#include "itkNumericTraits.h"
#include <cstdlib>

namespace itk
{
/** Constructor */
template< typename TInputMesh, typename TOutputImage >
PolyLineMeshToCellCountImageFilter< TInputMesh, TOutputImage >
::PolyLineMeshToCellCountImageFilter()
{
  this->SetNumberOfRequiredInputs(1);

  m_Size.Fill(0);
  m_Index.Fill(0);

  for ( unsigned int i = 0; i < 3; i++ )
    {
    m_Spacing[i] = 1.0;
    m_Origin[i] = 0;
    }

  m_InsideValue = NumericTraits< ValueType >::OneValue();
  m_OutsideValue = NumericTraits< ValueType >::ZeroValue();
  m_Direction.GetVnlMatrix().set_identity();

  m_Tolerance = 1e-5;
  m_InfoImage = nullptr;
}

/** Destructor */
template< typename TInputMesh, typename TOutputImage >
PolyLineMeshToCellCountImageFilter< TInputMesh, TOutputImage >
::~PolyLineMeshToCellCountImageFilter()
{}

/** Set the Input Mesh */
template< typename TInputMesh, typename TOutputImage >
void
PolyLineMeshToCellCountImageFilter< TInputMesh, TOutputImage >
::SetInput(TInputMesh *input)
{
  this->ProcessObject::SetNthInput(0, input);
}

/** Get the input Mesh */
template< typename TInputMesh, typename TOutputImage >
typename PolyLineMeshToCellCountImageFilter< TInputMesh, TOutputImage >::InputMeshType *
PolyLineMeshToCellCountImageFilter< TInputMesh, TOutputImage >
::GetInput(void)
{
  return static_cast< TInputMesh * >( this->GetInput() );
}

/** Get the input Mesh */
template< typename TInputMesh, typename TOutputImage >
typename PolyLineMeshToCellCountImageFilter< TInputMesh, TOutputImage >::InputMeshType *
PolyLineMeshToCellCountImageFilter< TInputMesh, TOutputImage >
::GetInput(unsigned int idx)
{
  return itkDynamicCastInDebugMode< TInputMesh * >
         ( this->ProcessObject::GetInput(idx) );
}

//----------------------------------------------------------------------------
template< typename TInputMesh, typename TOutputImage >
void
PolyLineMeshToCellCountImageFilter< TInputMesh, TOutputImage >
::SetSpacing(const double spacing[3])
{
  SpacingType s;
  for(unsigned int i = 0; i < TOutputImage::ImageDimension; ++i)
    {
    s[i] = static_cast<SpacePrecisionType>(spacing[i]);
    }
  this->SetSpacing(s);
}

template< typename TInputMesh, typename TOutputImage >
void
PolyLineMeshToCellCountImageFilter< TInputMesh, TOutputImage >
::SetSpacing(const float spacing[3])
{
  //Vector< float, 3 > sf(spacing);
  //SpacingType        s;
  //s.CastFrom(sf);
  //this->SetSpacing(s);
}

//----------------------------------------------------------------------------
template< typename TInputMesh, typename TOutputImage >
void
PolyLineMeshToCellCountImageFilter< TInputMesh, TOutputImage >
::SetOrigin(const double origin[3])
{
  PointType p(origin);
  this->SetOrigin(p);
}

template< typename TInputMesh, typename TOutputImage >
void
PolyLineMeshToCellCountImageFilter< TInputMesh, TOutputImage >
::SetOrigin(const float origin[3])
{
  //Point< float, 3 > of(origin);
  //PointType         p;
  //p.CastFrom(of);
  //this->SetOrigin(p);
}





//----------------------------------------------------------------------------

/** Update */
template< typename TInputMesh, typename TOutputImage >
void
PolyLineMeshToCellCountImageFilter< TInputMesh, TOutputImage >
::GenerateData(void)
{

  itkDebugMacro(<< "PolyLineMeshToCellCountImageFilter::Update() called");

  // Get the input and output pointers
  OutputImagePointer OutputImage = this->GetOutput();
  InputMeshPointer InputMesh = this->GetInput(0);

  if ( m_InfoImage != nullptr)
    {
    itkDebugMacro(<< "Using info image");
    m_Size = m_InfoImage->GetLargestPossibleRegion().GetSize();
    m_Index = m_InfoImage->GetLargestPossibleRegion().GetIndex();
    m_Spacing = m_InfoImage->GetSpacing();
    m_Origin = m_InfoImage->GetOrigin();
    m_Direction = m_InfoImage->GetDirection();
    }

  if ( m_Size[0] == 0 ||  m_Size[1] == 0 ||  m_Size[2] == 0 )
    {
    itkExceptionMacro(<< "Must Set Image Size");
    }

  typename OutputImageType::RegionType region;

  region.SetSize (m_Size);
  region.SetIndex(m_Index);

  OutputImage->SetLargestPossibleRegion(region); //
  OutputImage->SetBufferedRegion(region);        // set the region
  OutputImage->SetRequestedRegion(region);       //
  OutputImage->SetSpacing(m_Spacing);            // set spacing
  OutputImage->SetOrigin(m_Origin);              //   and origin
  OutputImage->SetDirection(m_Direction);        // direction cosines
  OutputImage->Allocate();
  OutputImage->FillBuffer(0);

  m_CellMask = OutputImageType::New();
  m_CellMask->CopyInformation( OutputImage );
  m_CellMask->SetLargestPossibleRegion(region); //
  m_CellMask->SetBufferedRegion(region);        // set the region
  m_CellMask->SetRequestedRegion(region);
  m_CellMask->Allocate();
  m_CellMask->FillBuffer(0);

  for (unsigned long i=0; i<this->GetInput(0)->GetNumberOfCells(); i++ ) {
    CountCellVoxelHits( i );
  }

  itkDebugMacro(<< "PolyLineMeshToCellCountImageFilter::Update() finished");

} // end update function

template< typename TInputMesh, typename TOutputImage >
void
PolyLineMeshToCellCountImageFilter< TInputMesh, TOutputImage >::
CountCellVoxelHits( unsigned long id ) {

  // Binary mask for a given cell
  m_CellMask->FillBuffer(0);
  CellAutoPointer cell;
  if ( this->GetInput(0)->GetCell(id, cell) ) {
    ContinuousIndexType cidx;
    IndexType idx;
    for (unsigned int i=0; i<cell->GetNumberOfPoints(); i++ ) {
      PointType pt = this->GetInput(0)->GetPoint( cell->GetPointIds()[i] );
      this->GetOutput()->TransformPhysicalPointToContinuousIndex(pt, cidx);
      idx.CopyWithRound(cidx);
      m_CellMask->SetPixel(idx,1);
    }

    // Update output by adding the binary mask
    itk::ImageRegionIteratorWithIndex< OutputImageType > it( this->GetOutput(),
      this->GetOutput()->GetLargestPossibleRegion() );
    while ( !it.IsAtEnd() ) {
      if ( m_CellMask->GetPixel( it.GetIndex() )) {
        it.Set( it.Value() + 1 );
      }
      ++it;
    }

  }


}


template< typename TInputMesh, typename TOutputImage >
void
PolyLineMeshToCellCountImageFilter< TInputMesh, TOutputImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Size : " << m_Size << std::endl;
  os << indent << "Inside Value : "
     << static_cast< typename NumericTraits< ValueType >::PrintType >( m_InsideValue ) << std::endl;
  os << indent << "Outside Value : "
     << static_cast< typename NumericTraits< ValueType >::PrintType >( m_OutsideValue ) << std::endl;
  os << indent << "Tolerance: " << m_Tolerance << std::endl;
  os << indent << "Origin: " << m_Origin << std::endl;
  os << indent << "Spacing: " << m_Spacing << std::endl;
  os << indent << "Direction: " << std::endl << m_Direction << std::endl;
  os << indent << "Index: " << m_Index << std::endl;
}
} // end namespace itk

#endif
