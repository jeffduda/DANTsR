/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkCaminoStreamlineFileReader.hxx,v $
  Language:  C++
  Date:      $Date: 2009/03/13 19:48:16 $
  Version:   $Revision: 1.23 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkCaminoStreamlineFileReader_hxx
#define __itkCaminoStreamlineFileReader_hxx

#include "itkCaminoStreamlineFileReader.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLabelContourImageFilter.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkByteSwapper.h"

#include <fstream>
#include <stdio.h>
#include <string>

namespace itk {

//
// Constructor
//
template<class TOutputMesh>
CaminoStreamlineFileReader<TOutputMesh>
::CaminoStreamlineFileReader()
{
  //
  // Create the output
  //
  typename TOutputMesh::Pointer output = TOutputMesh::New();
  this->ProcessObject::SetNumberOfRequiredOutputs( 1 );
  this->ProcessObject::SetNthOutput( 0, output.GetPointer() );

  this->m_ReadFloat = true;
  this->m_Seeds = SeedSetType::New();
}

template<class TOutputMesh>
template<typename PrecisionType>
void
CaminoStreamlineFileReader<TOutputMesh>
::ReadCaminoTract()
{
  //std::cout << "itk::CaminoStreamlineFileReader::ReadCaminoTract<PrecisionType>()" << std::endl;

  typename OutputMeshType::Pointer outputMesh = this->GetOutput();
  long lineId = outputMesh->GetNumberOfCells();
  CellAutoPointer streamline;

  PrecisionType nTractPoints;
  this->m_InputFile.read( reinterpret_cast< char * >( &nTractPoints ), sizeof(PrecisionType) );
  ByteSwapper<PrecisionType>::SwapRangeFromSystemToBigEndian(&nTractPoints,1);

  //std::cout << "  nPoints=" << nTractPoints << std::endl;

  if ( !this->m_InputFile.eof() ) {

    PrecisionType seed;
    this->m_InputFile.read( reinterpret_cast< char * >( &seed ), sizeof(PrecisionType) );
    ByteSwapper<PrecisionType>::SwapRangeFromSystemToBigEndian(&seed,1);

    PrecisionType * pts = new PrecisionType[static_cast<unsigned int>(3*nTractPoints)];
    this->m_InputFile.read( reinterpret_cast< char * >( pts ), (3*nTractPoints)*sizeof(PrecisionType) );
    ByteSwapper<PrecisionType>::SwapRangeFromSystemToBigEndian(pts, (3*nTractPoints));

    this->m_Seeds->InsertElement(lineId, static_cast<unsigned long>(seed) );

    unsigned int nPoints = static_cast<unsigned int>( nTractPoints );

    LineType polyLine;
    polyLine.SetSize( nPoints );

    typename OutputMeshType::PointIdentifier polyPoints[ nPoints ];

    for (unsigned long i=0; i<nPoints; i++) {
      typename OutputMeshType::PointType point;
      point[0] = static_cast<PixelType>(pts[3*i]);
      point[1] = static_cast<PixelType>(pts[3*i + 1]);
      point[2] = static_cast<PixelType>(pts[3*i + 2]);
      long pointId = outputMesh->GetNumberOfPoints();
      outputMesh->SetPoint(pointId,point);
      polyPoints[i] = pointId;

      polyLine[i] = pointId;
      ++pointId;
    }

    PolyLineCellType * polygon = new PolyLineCellType;
    polygon->SetPointIds( 0, nPoints, polyPoints );
    streamline.TakeOwnership( polygon );
    outputMesh->SetCell(lineId, streamline);

    delete[] pts;

  }

}


template<class TOutputMesh>
void
CaminoStreamlineFileReader<TOutputMesh>
::GenerateData()
{
  if( this->m_FileName == "" )
    {
    itkExceptionMacro( "No input FileName" );
    return;
    }

    /**
     * Get filename extension
     */
    std::string::size_type pos = this->m_FileName.rfind( "." );
    std::string extension( this->m_FileName, pos+1, this->m_FileName.length()-1 );
    if( extension == "Bfloat" )
      {
      m_ReadFloat = true;
      }
    else if( extension == "Bdouble" )
      {
      m_ReadFloat = false;
      }
    else
      {
      itkExceptionMacro( "Unknown extension: " << extension );
      return;
      }

  this->m_InputFile.open( this->m_FileName.c_str(), std::ifstream::binary );
  if (!this->m_InputFile )
    {
    itkExceptionMacro("Unable to open input file: ");
    return;
    }

  typename OutputMeshType::Pointer outputMesh = this->GetOutput();
  outputMesh->GetPoints()->Initialize();

  outputMesh->SetCellsAllocationMethod( OutputMeshType::CellsAllocatedDynamicallyCellByCell );

  this->m_Seeds->Initialize();

  CellAutoPointer streamline;

  while (!this->m_InputFile.eof()) {

    if ( m_ReadFloat )
    {
      this->ReadCaminoTract<float>();
    }
    else
    {
      this->ReadCaminoTract<double>();
    }
  }
}


template<class TOutputMesh>
void
CaminoStreamlineFileReader<TOutputMesh>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "FileName: "
     << this->m_FileName << std::endl;
}

} //end of namespace itk


#endif
