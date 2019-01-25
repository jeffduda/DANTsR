/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTrackVisStreamlineFileReader.hxx,v $
  Language:  C++
  Date:      $Date: 2009/03/04 23:10:58 $
  Version:   $Revision: 1.18 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTrackVisStreamlineFileReader_hxx
#define __itkTrackVisStreamlineFileReader_hxx

#include "itkTrackVisStreamlineFileReader.h"
#include "itkImageFileReader.h"
#include "itkByteSwapper.h"

#include <fstream>

namespace itk {

//
// Constructor
//
template<class TOutputMesh>
TrackVisStreamlineFileReader<TOutputMesh>
::TrackVisStreamlineFileReader()
{
  this->m_FileName = "";
  this->m_MultiComponentScalarSets = NULL;

  this->m_NScalars = 0;
  this->m_NProperties = 0;
}

//
// Destructor
//
template<class TOutputMesh>
TrackVisStreamlineFileReader<TOutputMesh>
::~TrackVisStreamlineFileReader()
{
}

//
// Write the input mesh to the output file
//
template<class TOutputMesh>
void TrackVisStreamlineFileReader<TOutputMesh>
::Update()
{
  this->GenerateData();
}

//
// Write the input mesh to the output file
//
template<class TOutputMesh>
void TrackVisStreamlineFileReader<TOutputMesh>
::Read()
{
  this->GenerateData();
}

template<class TOutputMesh>
void
TrackVisStreamlineFileReader<TOutputMesh>
::GenerateData()
{
  if( this->m_FileName == "" )
    {
    itkExceptionMacro( "No FileName" );
    return;
    }

  /**
   * Get filename extension
   */
  std::string::size_type pos = this->m_FileName.rfind( "." );
  std::string extension( this->m_FileName, pos+1, this->m_FileName.length()-1 );

  if( extension != "trk" )
    {
    itkExceptionMacro( "Unknown extension: " << extension );
    return;
    }

  //
  // Read input file
  //
  std::ifstream inputFile( m_FileName.c_str() );
  m_InputFile.open ( m_FileName.c_str(), std::ifstream::in);

  if( !m_InputFile.is_open() )
    {
    itkExceptionMacro( "Unable to open file\n"
        "inputFilename= " << m_FileName );
    return;
    }

  this->ReadTrkFile();
  m_InputFile.close();

}
template<class TOutputMesh>
void
TrackVisStreamlineFileReader<TOutputMesh>
::ReadTrkFile()
{
  TRACKVIS_HEADER_V2 hdr = this->ReadTrkHeader();

  typename OutputMeshType::Pointer outputMesh = this->GetOutput();
  outputMesh->SetCellsAllocationMethod( OutputMeshType::CellsAllocatedDynamicallyCellByCell );
  outputMesh->GetPoints()->Initialize();

  m_NScalars = hdr.n_scalars;
  m_NProperties = hdr.n_properties;

  m_ReferenceImage = ImageType::New();
  m_ReferenceImage->SetSpacing( hdr.voxel_size );
  m_ReferenceImage->SetOrigin( hdr.origin );
  typename ImageType::RegionType region;
  typename ImageType::DirectionType dir;

  for (unsigned int i=0; i<3; i++ ) {
    region.SetSize( i, hdr.dim[i] );
    dir(0,i) = hdr.image_orientation_patient[i];
    dir(1,i) = hdr.image_orientation_patient[i+3];
  }
  dir(2,0) = dir(0,1)*dir(1,2) - dir(0,2)*dir(1,1);
  dir(2,1) = dir(0,2)*dir(1,0) - dir(0,0)*dir(1,2);
  dir(2,2) = dir(0,0)*dir(1,1) - dir(0,1)*dir(1,0);

  m_ReferenceImage->SetRegions( region );
  m_ReferenceImage->SetDirection( dir );

  for ( int i=0; i<hdr.n_count; i++ ) {
    this->ReadTrkTract();
  }
}

template<class TOutputMesh>
void
TrackVisStreamlineFileReader<TOutputMesh>
::PrintTrkHeader( TRACKVIS_HEADER_V2 hdr ) {

  std::cout << "---------- " << this->m_FileName << " ----------" << std::endl;
  std::cout << "id_string: " << hdr.id_string << std::endl;
  std::cout << "dim: " << hdr.dim[0] << " x " << hdr.dim[1] << " x " << hdr.dim[2] << std::endl;
  std::cout << "voxel_size: " << hdr.voxel_size[0] << " x " << hdr.voxel_size[1] << " x " << hdr.voxel_size[2] << std::endl;
  std::cout << "origin: " << hdr.origin[0] << " x " << hdr.origin[1] << " x " << hdr.origin[2] << std::endl;
  std::cout << "n_scalars: " << hdr.n_scalars << std::endl;
  std::cout << "n_properties: " << hdr.n_properties << std::endl;
  std::cout << "voxel_order: " << hdr.voxel_order << std::endl;
  std::cout << "image_orientation_patient: (" << hdr.image_orientation_patient[0] << "," <<
    hdr.image_orientation_patient[1] << "," <<
    hdr.image_orientation_patient[2] << "), (" <<
    hdr.image_orientation_patient[3] << "," <<
    hdr.image_orientation_patient[4] << "," <<
    hdr.image_orientation_patient[5] << ")" << std::endl;
  std::cout << "invert_x: " << (int) hdr.invert_x << std::endl;
  std::cout << "invert_y: " << (int)hdr.invert_y << std::endl;
  std::cout << "invert_z: " << (int)hdr.invert_z << std::endl;
  std::cout << "swap_xy: " << (int)hdr.swap_xy << std::endl;
  std::cout << "swap_yz: " << (int)hdr.swap_yz << std::endl;
  std::cout << "swap_zx: " << (int) hdr.swap_zx << std::endl;
  std::cout << "n_count: " << hdr.n_count << std::endl;
  std::cout << "version: " << hdr.version << std::endl;
  std::cout << "hdr_size: " << hdr.hdr_size << std::endl;
  std::cout << "---------- " << this->m_FileName << " ----------" << std::endl;

}

template<class TOutputMesh>
typename TrackVisStreamlineFileReader<TOutputMesh>::TRACKVIS_HEADER_V2
TrackVisStreamlineFileReader<TOutputMesh>
::ReadTrkHeader( ) {
  std::cout << "TrackVisStreamlineFileReader<TOutputMesh>::ReadTrkHeader()"  << std::endl;

  TRACKVIS_HEADER_V2 hdr;
  m_InputFile.read( reinterpret_cast<char *>(&hdr), sizeof(hdr));

  this->PrintTrkHeader( hdr );

  return(hdr);

}

template<class TOutputMesh>
void
TrackVisStreamlineFileReader<TOutputMesh>
::ReadTrkTract()
{
  typename OutputMeshType::Pointer outputMesh = this->GetOutput();
  int nPoints;
  m_InputFile.read( reinterpret_cast<char *>( &nPoints), sizeof(int) );
  int bufferSize = nPoints*(3+m_NScalars) + m_NProperties;
  float *buffer = new float[bufferSize];
  m_InputFile.read( reinterpret_cast<char *>( buffer), bufferSize*sizeof(float));

  CellAutoPointer streamline;
  typename OutputMeshType::PointIdentifier polyPoints[ nPoints ];

  int idx = 0;
  for ( int i=0; i<nPoints; i++ ) {

    typename OutputMeshType::PointType point;
    point[0] = buffer[idx++];
    point[1] = buffer[idx++];
    point[2] = buffer[idx++];
    polyPoints[i] = outputMesh->GetNumberOfPoints();
    outputMesh->GetPoints()->InsertElement(outputMesh->GetNumberOfPoints(),point);

    // FIXME - read in point scalars
    for (unsigned int j=0; j<m_NScalars; j++) {
      idx++;
    }

  }
  // FIXME - read in tract properties
  for (int i=0; i<m_NProperties; i++) {
    idx++;
  }

  PolyLineCellType * polyline = new PolyLineCellType;
  polyline->SetPointIds( 0, nPoints, polyPoints );
  streamline.TakeOwnership( polyline );
  outputMesh->SetCell(outputMesh->GetNumberOfCells(), streamline);

  delete[] buffer;
}




template<class TOutputMesh>
void
TrackVisStreamlineFileReader<TOutputMesh>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "FileName: " << this->m_FileName << std::endl;
}

} //end of namespace itk

#endif
