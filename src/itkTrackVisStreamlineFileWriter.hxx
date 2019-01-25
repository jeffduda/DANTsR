/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTrackVisStreamlineFileWriter.hxx,v $
  Language:  C++
  Date:      $Date: 2009/03/04 23:10:58 $
  Version:   $Revision: 1.18 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTrackVisStreamlineFileWriter_hxx
#define __itkTrackVisStreamlineFileWriter_hxx

#include "itkTrackVisStreamlineFileWriter.h"
#include "itkImageFileWriter.h"
#include "itkByteSwapper.h"

#include <fstream>

namespace itk {

//
// Constructor
//
template<class TInputMesh,class TInputImage>
TrackVisStreamlineFileWriter<TInputMesh,TInputImage>
::TrackVisStreamlineFileWriter()
{
  this->m_Input = NULL;
  this->m_FileName = "";
  this->m_MultiComponentScalarSets = NULL;
  this->m_ImageSize.Fill( 0 );
  this->m_ReferenceImage = nullptr;

  this->m_NScalars = 0;
  this->m_NProperties = 0;
}

//
// Destructor
//
template<class TInputMesh,class TInputImage>
TrackVisStreamlineFileWriter<TInputMesh,TInputImage>
::~TrackVisStreamlineFileWriter()
{
}

//
//
//
template<class TInputMesh,class TInputImage>
void
TrackVisStreamlineFileWriter<TInputMesh,TInputImage>
::SetInput(InputMeshType * input)
{
  this->m_Input = input;
}

//
// Write the input mesh to the output file
//
template<class TInputMesh,class TInputImage>
void TrackVisStreamlineFileWriter<TInputMesh,TInputImage>
::Update()
{
  this->GenerateData();
}

//
// Write the input mesh to the output file
//
template<class TInputMesh,class TInputImage>
void TrackVisStreamlineFileWriter<TInputMesh,TInputImage>
::Write()
{
  this->GenerateData();
}

template<class TInputMesh,class TInputImage>
void
TrackVisStreamlineFileWriter<TInputMesh,TInputImage>
::GenerateData()
{
  if( this->m_FileName == "" )
    {
    itkExceptionMacro( "No FileName" );
    return;
    }

  if (this->m_ReferenceImage == nullptr) {
    itkExceptionMacro( "No Reference Image" );
    return;
    }

  if (this->m_NScalars > 0 ) {
    itkExceptionMacro( "Scalars not yet implemented" );
    return;
  }

  if (this->m_NProperties > 0 ) {
    itkExceptionMacro( "Properties not yet implemented" );
    return;
  }

  //
  // Read output file
  //
  std::ofstream outputFile( m_FileName.c_str() );

  if( !outputFile.is_open() )
    {
    itkExceptionMacro( "Unable to open file\n"
        "outputFilename= " << m_FileName );
    return;
    }
  else
    {
    outputFile.close();
    }

  /**
   * Get filename extension
   */
  std::string::size_type pos = this->m_FileName.rfind( "." );
  std::string extension( this->m_FileName, pos+1, this->m_FileName.length()-1 );

  if( extension == "trk" )
    {
    this->WriteTrkFile();
    }
  else
    {
    itkExceptionMacro( "Unknown extension: " << extension );
    return;
    }

}
template<class TInputMesh,class TInputImage>
void
TrackVisStreamlineFileWriter<TInputMesh,TInputImage>
::WriteTrkFile()
{
  this->WriteTrkHeader();
  this->WriteTrkTracts();
}

template<class TInputMesh,class TInputImage>
void
TrackVisStreamlineFileWriter<TInputMesh,TInputImage>
::WriteTrkHeader( ) {
  std::cout << "TrackVisStreamlineFileWriter<TInputMesh,TInputImage>::WriteTrkHeader()"  << std::endl;
  std::ofstream outputFile( this->m_FileName.c_str(), std::ofstream::out );
  outputFile << "TRACK ";

  // Write reference image info (needed why?)
  short int dim[3];
  float voxelSize[3];
  float origin[3];
  for ( unsigned int i=0; i<3; i++ ) {
    dim[i] = this->m_ReferenceImage->GetLargestPossibleRegion().GetSize()[i];
    voxelSize[i] = this->m_ReferenceImage->GetSpacing()[i];
    origin[i] = this->m_ReferenceImage->GetOrigin()[i];
  }

  outputFile.write( reinterpret_cast<char *>( dim ), 3 * sizeof(short int) );
  outputFile.write( reinterpret_cast<char *>( voxelSize ), 3 * sizeof(float) );
  outputFile.write( reinterpret_cast<char *>( origin ), 3 * sizeof(float) );

  short int n_scalars = this->m_NScalars;
  outputFile.write( reinterpret_cast<char *>( &n_scalars ), sizeof(short int) );

  // scalar names char[10][20]
  char scalar_names[200];
  for ( unsigned x=0; x<200; x++ ) {
    scalar_names[x] = 0;
  }
  outputFile.write( reinterpret_cast<char *>( scalar_names ), 200 );

  short int n_properties = this->m_NProperties;
  outputFile.write( reinterpret_cast<char *>( &n_properties ), sizeof(short int) );

  // property_name char[10][20]
  char property_names[200];
  for ( unsigned x=0; x<200; x++ ) {
    property_names[x] = 0;
  }
  outputFile.write( reinterpret_cast<char *>( property_names ), 200 );

  // float vox_to_ras[4][4]
  float vox_to_ras[16];
  for (unsigned int x=0; x<16; x++ ) {
    vox_to_ras[x] = 0.0;
  }
  outputFile.write( reinterpret_cast<char *>( vox_to_ras ), 16 * sizeof(float) );

  char reserved[444];
  for ( unsigned x=0; x<444; x++ ) {
    reserved[x] = 0;
  }
  outputFile.write( reserved, 444 );

  char voxel_order[4]; // Empty is "LPS"
  for ( unsigned x=0; x<4; x++ ) {
    voxel_order[x] = 0;
  }
  outputFile.write( voxel_order, 4 );

  char pad2[4]; // paddings?
  for ( unsigned x=0; x<4; x++ ) {
    pad2[x] = 0;
  }
  outputFile.write( pad2, 4 );

  float image_orientation_patient[6];
  image_orientation_patient[0] = 1;
  image_orientation_patient[1] = 0;
  image_orientation_patient[2] = 0;
  image_orientation_patient[3] = 0;
  image_orientation_patient[4] = 1;
  image_orientation_patient[5] = 0;
  outputFile.write( reinterpret_cast<char *>( image_orientation_patient ), 6 * sizeof(float) );

  char pad1[2];
  pad1[0] = 0;
  pad1[1] = 0;
  outputFile.write( pad1, 2 );

  unsigned char invert_x = 0;
  outputFile.write( reinterpret_cast<char *>(&invert_x), 1 );

  unsigned char invert_y = 0;
  outputFile.write( reinterpret_cast<char *>(&invert_y), 1 );

  unsigned char invert_z = 0;
  outputFile.write( reinterpret_cast<char *>(&invert_z), 1 );

  unsigned char swap_xy = 0;
  outputFile.write( reinterpret_cast<char *>(&swap_xy), 1 );

  unsigned char swap_yz = 0;
  outputFile.write( reinterpret_cast<char *>(&swap_yz), 1 );

  unsigned char swap_zx = 0;
  outputFile.write( reinterpret_cast<char *>(&swap_zx), 1 );

  int n_count = this->m_Input->GetNumberOfCells();
  n_count = 1;
  outputFile.write( reinterpret_cast<char *>(&n_count), sizeof(int) );

  int version=2;
  outputFile.write( reinterpret_cast<char *>( &version ), sizeof(int) );

  int hdr_size=1000;
  outputFile.write( reinterpret_cast<char *>( &hdr_size ), sizeof(int) );
  //ByteSwapper<float>::SwapRangeFromSystemToBigEndian(ptData,numberOfPoints*3);

  outputFile.close();
}

template<class TInputMesh,class TInputImage>
void
TrackVisStreamlineFileWriter<TInputMesh,TInputImage>
::WriteTrkTracts()
{
  //
  // Write to output file
  //
  std::cout << "TrackVisStreamlineFileWriter<TInputMesh,TInputImage>::WriteTrkTracts()"  << std::endl;
  std::ofstream outputFile( this->m_FileName.c_str(), std::ofstream::out | std::ofstream::app );

  unsigned long numberOfCells = this->m_Input->GetNumberOfCells();
  numberOfCells = 1;
  for (unsigned int i=0; i<numberOfCells; i++) {
    CellAutoPointer cell;
    if ( !this->m_Input->GetCell(i, cell) ) {
      itkExceptionMacro( "Cell could not be read" );
    }

    long nPoints = cell->GetNumberOfPoints();
    unsigned long dataSize = nPoints*(3+this->m_NScalars) + this->m_NProperties;
    float cellData[dataSize];

    std::cout << "Writing tract " << i << " with " << nPoints << " points " << std::endl;

    // Write point coordinates and associated scalars
    unsigned long idx = 0;
    for ( unsigned long j=0; j<nPoints; j++ ) {
      int pointId = cell->GetPointIds()[j];
      PointType point;
      bool pointExists =  this->m_Input->GetPoint( pointId, &point );
      if ( !pointExists ) {
        itkExceptionMacro( "Unknown point id found: " << j );
      }

      cellData[idx++] = point[0];
      cellData[idx++] = point[1];
      cellData[idx++] = point[2];

      for ( unsigned int k=0; k<this->m_NScalars; k++ ) {
        //cellData[idx++] = SCALAR VALUE GOES HERE
        }
      }

    for ( unsigned int k=0; k<this->m_NProperties; k++) {
      //cellData[idx++] = PROPERTTY VALUE GOES HERE
      }

    outputFile.write( reinterpret_cast<char *>( cellData ), dataSize * sizeof(float) );

    }

  outputFile.close();
}




template<class TInputMesh,class TInputImage>
void
TrackVisStreamlineFileWriter<TInputMesh,TInputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "FileName: " << this->m_FileName << std::endl;
}

} //end of namespace itk

#endif
