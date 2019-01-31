#ifndef __itkMRTrixStreamlineFileWriter_hxx
#define __itkMRTrixStreamlineFileWriter_hxx

#include "itkMRTrixStreamlineFileWriter.h"
#include "itkImageFileWriter.h"
#include "itkByteSwapper.h"

#include <fstream>

namespace itk {

//
// Constructor
//
template<class TInputMesh,class TInputImage>
MRTrixStreamlineFileWriter<TInputMesh,TInputImage>
::MRTrixStreamlineFileWriter()
{
  this->m_Input = nullptr;
  this->m_FileName = "";
  this->m_MultiComponentScalarSets = nullptr;
  this->m_ImageSize.Fill( 0 );
  this->m_ReferenceImage = nullptr;

  this->m_NScalars = 0;
  this->m_NProperties = 0;
}

//
// Destructor
//
template<class TInputMesh,class TInputImage>
MRTrixStreamlineFileWriter<TInputMesh,TInputImage>
::~MRTrixStreamlineFileWriter()
{
}

//
//
//
template<class TInputMesh,class TInputImage>
void
MRTrixStreamlineFileWriter<TInputMesh,TInputImage>
::SetInput(InputMeshType * input)
{
  this->m_Input = input;
}

//
// Write the input mesh to the output file
//
template<class TInputMesh,class TInputImage>
void MRTrixStreamlineFileWriter<TInputMesh,TInputImage>
::Update()
{
  this->GenerateData();
}

//
// Write the input mesh to the output file
//
template<class TInputMesh,class TInputImage>
void MRTrixStreamlineFileWriter<TInputMesh,TInputImage>
::Write()
{
  this->GenerateData();
}

template<class TInputMesh,class TInputImage>
void
MRTrixStreamlineFileWriter<TInputMesh,TInputImage>
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
  // Open output file
  //
  m_OutputFile.open( m_FileName.c_str(), std::ofstream::out );

  if( !m_OutputFile.is_open() )
    {
    itkExceptionMacro( "Unable to open file\n"
        "outputFilename= " << m_FileName );
    return;
    }

  this->WriteTrkFile();

  m_OutputFile.close();

}
template<class TInputMesh,class TInputImage>
void
MRTrixStreamlineFileWriter<TInputMesh,TInputImage>
::WriteTrkFile()
{
  this->WriteTrkHeader();

  for( unsigned int i=0; i<this->m_Input->GetNumberOfCells(); i++) {
    this->WriteTrkTract(i);
  }
}

template<class TInputMesh,class TInputImage>
void
MRTrixStreamlineFileWriter<TInputMesh,TInputImage>
::WriteTrkHeader( ) {
  //std::cout << "MRTrixStreamlineFileWriter<TInputMesh,TInputImage>::WriteTrkHeader()"  << std::endl;

  MRTrix_HEADER_V2 hdr;
  hdr.id_string[0] = 'T';
  hdr.id_string[1] = 'R';
  hdr.id_string[2] = 'A';
  hdr.id_string[3] = 'C';
  hdr.id_string[4] = 'K';
  hdr.id_string[5] = '\0';

  for ( unsigned int i=0; i<3; i++ ) {
    hdr.dim[i] = static_cast<short int>( this->m_ReferenceImage->GetLargestPossibleRegion().GetSize()[i] );
    hdr.voxel_size[i] = static_cast<float>( this->m_ReferenceImage->GetSpacing()[i] );
    hdr.origin[i] = static_cast<float>( this->m_ReferenceImage->GetOrigin()[i] );
  }

  hdr.n_scalars = static_cast<short int>(this->m_NScalars);
  hdr.n_properties = static_cast<short int>(this->m_NProperties);

  for ( unsigned int x=0; x<10; x++ ) {
    for (unsigned int y=0; y<20; y++) {
      hdr.scalar_names[x][y] = 0;
      hdr.property_names[x][y] = 0;
    }
  }

  for (unsigned int x=0; x<4; x++ ) {
    for (unsigned int y=0; y<4; y++) {
      hdr.vox_to_ras[x][y] = 0.0;
    }
  }

  for ( unsigned x=0; x<444; x++ ) {
    hdr.reserved[x] = 0;
  }

  hdr.voxel_order[0] = 'L';
  hdr.voxel_order[1] = 'P';
  hdr.voxel_order[2] = 'S';
  hdr.voxel_order[3] = '\0';

  for ( unsigned x=0; x<4; x++ ) {
    hdr.pad2[x] = 0;
  }

  hdr.image_orientation_patient[0] = 1;
  hdr.image_orientation_patient[1] = 0;
  hdr.image_orientation_patient[2] = 0;
  hdr.image_orientation_patient[3] = 0;
  hdr.image_orientation_patient[4] = 1;
  hdr.image_orientation_patient[5] = 0;

  hdr.pad1[0] = 0;
  hdr.pad1[1] = 0;

  hdr.invert_x = 0;
  hdr.invert_y = 0;
  hdr.invert_z = 0;
  hdr.swap_xy = 0;
  hdr.swap_yz = 0;
  hdr.swap_zx = 0;

  hdr.n_count = static_cast<int>( this->m_Input->GetNumberOfCells() );

  hdr.version = 2;

  hdr.hdr_size=1000;

  m_OutputFile.write( reinterpret_cast<char *>(&hdr), sizeof(hdr) );

}

template<class TInputMesh,class TInputImage>
void
MRTrixStreamlineFileWriter<TInputMesh,TInputImage>
::WriteTrkTract(unsigned int id)
{
  //
  // Write to output file
  //
  //std::cout << "MRTrixStreamlineFileWriter<TInputMesh,TInputImage>::WriteTrkTracts()"  << std::endl;

  CellAutoPointer cell;
  if ( !this->m_Input->GetCell(id, cell) ) {
    itkExceptionMacro( "Cell could not be read" );
  }

  int nPoints = cell->GetNumberOfPoints();
  m_OutputFile.write( reinterpret_cast<char *>(&nPoints), sizeof(int) );

  unsigned long dataSize = nPoints*(3+this->m_NScalars) + this->m_NProperties;
  float *cellData = new float[dataSize];

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

  m_OutputFile.write( reinterpret_cast<char *>( cellData ), dataSize * sizeof(float) );

  delete[] cellData;

}




template<class TInputMesh,class TInputImage>
void
MRTrixStreamlineFileWriter<TInputMesh,TInputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "FileName: " << this->m_FileName << std::endl;
}

} //end of namespace itk

#endif
