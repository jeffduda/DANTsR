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
template<class TInputMesh>
MRTrixStreamlineFileWriter<TInputMesh>
::MRTrixStreamlineFileWriter()
{
  this->m_Input = nullptr;
  this->m_FileName = "";
  this->m_SwapFlag = 2;
}

//
// Destructor
//
template<class TInputMesh>
MRTrixStreamlineFileWriter<TInputMesh>
::~MRTrixStreamlineFileWriter()
{
}

//
//
//
template<class TInputMesh>
void
MRTrixStreamlineFileWriter<TInputMesh>
::SetInput(InputMeshType * input)
{
  this->m_Input = input;
}

//
// Write the input mesh to the output file
//
template<class TInputMesh>
void MRTrixStreamlineFileWriter<TInputMesh>
::Update()
{
  this->GenerateData();
}

//
// Write the input mesh to the output file
//
template<class TInputMesh>
void MRTrixStreamlineFileWriter<TInputMesh>
::Write()
{
  this->GenerateData();
}

template<class TInputMesh>
void
MRTrixStreamlineFileWriter<TInputMesh>
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
  if( extension != "tck" )
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

  this->WriteTckFile();

  m_OutputFile.close();

}
template<class TInputMesh>
void
MRTrixStreamlineFileWriter<TInputMesh>
::WriteTckFile()
{
  this->WriteTckHeader();

  for( unsigned int i=0; i<this->m_Input->GetNumberOfCells(); i++) {
    this->WriteTckTract(i);
  }

  // write triplit of inf to finish
  float * finisher = new float[3];
  finisher[0] = finisher[1] = finisher[2] = std::numeric_limits<float>::infinity();

  if ( m_SwapFlag == 1 ) {
    ByteSwapper<float>::SwapRangeFromSystemToBigEndian(finisher, 3);
  }
  else if ( m_SwapFlag == 2 ) {
    ByteSwapper<float>::SwapRangeFromSystemToLittleEndian(finisher, 3);
  }

  m_OutputFile.write( reinterpret_cast<char *>( finisher ), 3 * sizeof(float) );
  delete[] finisher;

}

template<class TInputMesh>
void
MRTrixStreamlineFileWriter<TInputMesh>
::WriteTckHeader( ) {
  //std::cout << "MRTrixStreamlineFileWriter<TInputMesh>::WriteTckHeader()"  << std::endl;

  m_OutputFile << "mrtrix tracks" << std::endl;
  m_OutputFile << "count: " << std::to_string(this->m_Input->GetNumberOfCells()) << std::endl;

  if ( m_SwapFlag == 0 ) {
    m_OutputFile << "datatype: Float32" << std::endl;
  }
  else if ( m_SwapFlag == 1 ) {
    m_OutputFile << "datatype: Float32BE" << std::endl;
  }
  else if ( m_SwapFlag == 2 ) {
    m_OutputFile << "datatype: Float32LE" << std::endl;
  }

  /* other options - required?
  downsample_factor: 3
  fod_power: 0.25
  init_threshold: 0.1
  lmax: 8
  max_angle: 45
  max_dist: 250
  max_num_attempts: 50000
  max_num_tracks: 500
  max_seed_attempts: 1
  max_trials: 1000
  method: iFOD2
  min_dist: 4
  mrtrix_version: 0.3.12-325-gc203eda9
  output_step_size: 1.25
  rk4: 0
  samples_per_step: 4
  sh_precomputed: 1
  source: dwi2fod/out.mif
  step_size: 1.25
  stop_on_all_include: 0
  threshold: 0.1
  timestamp: 1443190374.822715044
  unidirectional: 0
  roi: seed mask.mif
  roi: mask mask.mif
  roi: mask mask.mif
  total_count: 1748
  */

  // Estimate size and add buffer of zeros
  std::streampos hsize = m_OutputFile.tellp();
  hsize += 30;
  long extra = 4 - ( static_cast<long>(hsize) % 4 );
  hsize += extra; // make divis by 4 for easier hex editor reading

  m_OutputFile << "file: . " + std::to_string( hsize ) << std::endl;
  m_OutputFile << "END" << std::endl;


  while ( m_OutputFile.tellp() < hsize ) {
    m_OutputFile << (char) 0;
  }

}

template<class TInputMesh>
void
MRTrixStreamlineFileWriter<TInputMesh>
::WriteTckTract(unsigned int id)
{
  //
  // Write to output file
  //
  //std::cout << "MRTrixStreamlineFileWriter<TInputMesh>::WriteTckTracts()"  << std::endl;

  CellAutoPointer cell;
  if ( !this->m_Input->GetCell(id, cell) ) {
    itkExceptionMacro( "Cell could not be read" );
  }

  int nPoints = cell->GetNumberOfPoints();
  //m_OutputFile.write( reinterpret_cast<char *>(&nPoints), sizeof(int) );

  unsigned long dataSize = (nPoints+1)*3;
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

    }

  cellData[idx++] = std::numeric_limits<float>::quiet_NaN();
  cellData[idx++] = std::numeric_limits<float>::quiet_NaN();
  cellData[idx++] = std::numeric_limits<float>::quiet_NaN();

  if ( m_SwapFlag == 1 ) {
    ByteSwapper<float>::SwapRangeFromSystemToBigEndian(cellData, dataSize);
  }
  else if ( m_SwapFlag == 2 ) {
    ByteSwapper<float>::SwapRangeFromSystemToLittleEndian(cellData, dataSize);
  }

  m_OutputFile.write( reinterpret_cast<char *>( cellData ), dataSize * sizeof(float) );

  delete[] cellData;
}




template<class TInputMesh>
void
MRTrixStreamlineFileWriter<TInputMesh>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "FileName: " << this->m_FileName << std::endl;
}

} //end of namespace itk

#endif
