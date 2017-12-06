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
    std::cerr << "Unable to open file: " << this->m_FileName << std::endl;
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

  Rcpp::Rcout << "Reading " << nTracts << " streamlines with " << nPoints << " points" << std::endl;

  this->m_InputFile.close();
  this->m_InputFile.open( this->m_FileName.c_str(), std::ifstream::binary );

  typename OutputMeshType::Pointer outputMesh = this->GetOutput();
  outputMesh->GetPoints()->Initialize();
  outputMesh->GetPoints()->Reserve(nPoints);

  long pointId = 0;
  long cellId = 0;
  while (!this->m_InputFile.eof()) {
    float * nTractPoints = new float;
    this->m_InputFile.read( reinterpret_cast< char * >( nTractPoints ), sizeof(n) );
    ByteSwapper<float>::SwapRangeFromSystemToBigEndian(nTractPoints,1);

    float pts[(int)(3*nTractPoints[0])+1];
    this->m_InputFile.read( reinterpret_cast< char * >( pts ), (3*nTractPoints[0]+1)*sizeof(n) );
    ByteSwapper<float>::SwapRangeFromSystemToBigEndian(pts, (3*nTractPoints[0])+1);

    for (unsigned long i=0; i<nTractPoints[0]; i++) {
      typename OutputMeshType::PointType point;
      point[0] = pts[3*i + 1];
      point[1] = pts[3*i + 2];
      point[2] = pts[3*i + 3];
      outputMesh->GetPoints()->InsertElement(i,point);
    }

    Rcpp::Rcout << "Done reading, testing stuff" << std::endl;

    typename OutputMeshType::Pointer m2 = OutputMeshType::New();
    m2->Initialize();
    PointType p0;
    PointType p1;
    PointType p2;
    p0[0] = -1.0; p0[1] = 0.0; p0[2] = 0.0;
    p1[0] =  1.0; p1[1] = 0.0; p1[2] = 0.0;
    p2[0] =  1.0; p2[1] = 1.0; p2[2] = 0.0;
    m2->SetPoint( 0, p0 );
    m2->SetPoint( 1, p1 );
    m2->SetPoint( 2, p2 );


    CellAutoPointer line0;
    line0.TakeOwnership( new LineCellType );
    line0->SetPointId( 0, 0 ); // line between points 0 and 1
    line0->SetPointId( 1, 1 );
    m2->SetCell( 0, line0 );
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
