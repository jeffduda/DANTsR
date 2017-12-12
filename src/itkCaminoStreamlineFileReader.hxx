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

  this->m_Lines = LineSetType::New();
  this->m_Seeds = SeedSetType::New();
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

  this->m_InputFile.open( this->m_FileName.c_str(), std::ifstream::binary );
  if (!this->m_InputFile )
    {
    itkExceptionMacro("Unable to open input file: ");
    return;
    }

  unsigned long nPoints=0;
  unsigned long nTracts=0;
  float n;

  while (!this->m_InputFile.eof()) {
    float * nTractPoints = new float;
    this->m_InputFile.read( reinterpret_cast< char * >( nTractPoints ), sizeof(n) );
    ByteSwapper<float>::SwapRangeFromSystemToBigEndian(nTractPoints,1);

    float skip[(int)(3*nTractPoints[0])+1];
    this->m_InputFile.read( reinterpret_cast< char * >( skip ), (3*nTractPoints[0]+1)*sizeof(n) );
    //ByteSwapper<float>::SwapRangeFromSystemToBigEndian(skip, (3*nTractPoints[0])+1);

    nPoints += nTractPoints[0];
    nTracts++;
  }

  this->m_InputFile.close();
  this->m_InputFile.open( this->m_FileName.c_str(), std::ifstream::binary );

  typename OutputMeshType::Pointer outputMesh = this->GetOutput();
  outputMesh->GetPoints()->Initialize();
  outputMesh->GetPoints()->Reserve(nPoints);

  outputMesh->SetCellsAllocationMethod( OutputMeshType::CellsAllocatedDynamicallyCellByCell );

  long pointId = 0;
  long cellId = 0;
  long lineId = 0;
  this->m_Lines->Initialize();
  this->m_Seeds->Initialize();

  CellAutoPointer streamline;

  while (!this->m_InputFile.eof()) {
    float * nTractPoints = new float;
    this->m_InputFile.read( reinterpret_cast< char * >( nTractPoints ), sizeof(n) );
    ByteSwapper<float>::SwapRangeFromSystemToBigEndian(nTractPoints,1);

    if ( nTractPoints[0] > 0 ) {

      float pts[(int)(3*nTractPoints[0])+1];
      this->m_InputFile.read( reinterpret_cast< char * >( pts ), (3*nTractPoints[0]+1)*sizeof(n) );
      ByteSwapper<float>::SwapRangeFromSystemToBigEndian(pts, (3*nTractPoints[0])+1);

      this->m_Seeds->InsertElement(lineId,(unsigned long)pts[0]);

      unsigned int nPoints = static_cast<unsigned int>( nTractPoints[0] );

      LineType polyLine;
      polyLine.SetSize( nPoints );

      typename OutputMeshType::PointIdentifier polyPoints[ nPoints ];

      for (unsigned long i=0; i<nPoints; i++) {
        typename OutputMeshType::PointType point;
        point[0] = pts[3*i + 1];
        point[1] = pts[3*i + 2];
        point[2] = pts[3*i + 3];
        outputMesh->SetPoint(pointId,point);
        polyPoints[i] = pointId;

        polyLine[i] = pointId;
        ++pointId;
      }

      PolygonCellType * polygon = new PolygonCellType;
      polygon->SetPointIds( 0, nPoints, polyPoints );
      streamline.TakeOwnership( polygon );
      outputMesh->SetCell(lineId, streamline);

      this->m_Lines->InsertElement(lineId, polyLine);

      if ( lineId == 0 ) {
        Rcpp::Rcout << nPoints << std::endl;
        Rcpp::Rcout << polygon->GetNumberOfPoints() << std::endl;
        Rcpp::Rcout << streamline->GetNumberOfPoints() << std::endl;
        CellAutoPointer getCell;
        outputMesh->GetCell(0, getCell);
        Rcpp::Rcout << getCell->GetNumberOfPoints() << std::endl;
      }

      ++lineId;

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
