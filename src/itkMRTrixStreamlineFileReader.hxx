#ifndef __itkMRTrixStreamlineFileReader_hxx
#define __itkMRTrixStreamlineFileReader_hxx

#include "itkMRTrixStreamlineFileReader.h"
#include "itkImageFileReader.h"
#include "itkByteSwapper.h"

#include <fstream>

namespace itk {

//
// Constructor
//
template<class TOutputMesh>
MRTrixStreamlineFileReader<TOutputMesh>
::MRTrixStreamlineFileReader()
{
  this->m_FileName = "";
  this->m_SwapFlag = 0;
}

//
// Destructor
//
template<class TOutputMesh>
MRTrixStreamlineFileReader<TOutputMesh>
::~MRTrixStreamlineFileReader()
{
}

//
// Write the input mesh to the output file
//
template<class TOutputMesh>
void MRTrixStreamlineFileReader<TOutputMesh>
::Update()
{
  this->GenerateData();
}

//
// Write the input mesh to the output file
//
template<class TOutputMesh>
void MRTrixStreamlineFileReader<TOutputMesh>
::Read()
{
  this->GenerateData();
}

template<class TOutputMesh>
void
MRTrixStreamlineFileReader<TOutputMesh>
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

  this->ReadTckFile();
  m_InputFile.close();

}
template<class TOutputMesh>
void
MRTrixStreamlineFileReader<TOutputMesh>
::ReadTckFile()
{
  m_Header.clear();
  this->ReadTckHeader();
  //this->PrintTckHeader();

  std::string countStr = this->GetHeaderValue( "count", true );
  std::string fileStr  = this->GetHeaderValue( "file", true );
  std::string datatype = this->GetHeaderValue( "datatype", true );

  if ( datatype.find("Float32") > datatype.size() ) {
    m_InputFile.close();
    itkExceptionMacro( "Unsupported data type" << datatype )
    return;
  }

  if ( !(datatype.find("BE") > datatype.size()) ) {
    m_SwapFlag = 1;
  }
  else if ( !(datatype.find("LE") > datatype.size()) ) {
    m_SwapFlag = 2;
  }

  unsigned long count = std::stoi( countStr );
  unsigned long offset = std::stoi( fileStr.substr( fileStr.find(".")+1, fileStr.size() - 1 - fileStr.find(".")) );

  m_InputFile.seekg( 0 );
  std::streampos fsize = m_InputFile.tellg();
  m_InputFile.seekg( 0, std::ios::end );
  fsize = m_InputFile.tellg() - fsize;

  m_InputFile.seekg(offset);

  //unsigned long nValues = fsize / sizeof(float);
  //unsigned long nPoints = std::floor( (nValues - 3*count)/3 );
  //std::cout << "file size: " << fsize << std::endl;
  //std::cout << "n values: " << nValues << std::endl;
  //std::cout << "n tracks: " << count << std::endl;
  //std::cout << "n points: " << nPoints << std::endl;

  typename OutputMeshType::Pointer outputMesh = this->GetOutput();
  outputMesh->SetCellsAllocationMethod( OutputMeshType::CellsAllocatedDynamicallyCellByCell );
  outputMesh->GetPoints()->Initialize();

  long dataSize = (unsigned long)fsize - offset;
  this->ReadTckTracts(dataSize, count);

}

template<class TOutputMesh>
std::string
MRTrixStreamlineFileReader<TOutputMesh>
::GetHeaderValue( std::string key, bool required ) {
  typename HeaderType::const_iterator pos = m_Header.find(key);

  if (pos == m_Header.end()) {
    if ( required ) {
      itkExceptionMacro( "Required key not found in header:" << key );
    }
    return("");
  }

  return(pos->second);
}

template<class TOutputMesh>
void
MRTrixStreamlineFileReader<TOutputMesh>
::PrintTckHeader( void ) {
  std::cout << "-------" << m_FileName << "-------" << std::endl;
  typename HeaderType::iterator it = m_Header.begin();
  while(it != m_Header.end())
  {
    std::cout<<it->first<<" :: "<<it->second<<std::endl;
    it++;
  }
  std::cout << "-------" << m_FileName << "-------" << std::endl;
}

template<class TOutputMesh>
void
MRTrixStreamlineFileReader<TOutputMesh>
::ReadTckHeader( ) {
  //std::cout << "MRTrixStreamlineFileReader<TOutputMesh>::ReadTckHeader()"  << std::endl;

  std::string line;
  std::getline(m_InputFile, line);
  if ( line.find("mrtrix mrtracks") < 0 ) {
    m_InputFile.close();
    itkExceptionMacro( "Does not appear to be an 'mrtrix tracks' file: '" << line << "'" );
    return;
  }

  while ( (!m_InputFile.eof()) && (line.find("END") != 0) ) {
    std::getline(m_InputFile, line);
    std::istringstream valueLine(line);
    std::vector<std::string> pair;
    std::string tmp;
    while ( std::getline(valueLine,tmp,':') ) {
      pair.push_back(tmp);
    }
    if ( pair.size() == 2) {
      m_Header.insert( std::make_pair(pair[0],pair[1]));
    }
  }

}

template<class TOutputMesh>
void
MRTrixStreamlineFileReader<TOutputMesh>
::ReadTckTracts(long dataSize, unsigned long count)
{
  typename OutputMeshType::Pointer outputMesh = this->GetOutput();
  CellAutoPointer streamline;

  unsigned long nFloats = dataSize / sizeof(float);
  float * data = new float[nFloats];

  this->m_InputFile.read( reinterpret_cast< char * >( data ), dataSize );

  if ( m_SwapFlag == 1) {
    ByteSwapper<float>::SwapRangeFromSystemToBigEndian(data, nFloats);
  }
  else if ( m_SwapFlag == 2 ) {
    ByteSwapper<float>::SwapRangeFromSystemToLittleEndian(data, nFloats);
  }

  typename OutputMeshType::PointIdentifier * cellBuffer = new typename OutputMeshType::PointIdentifier[nFloats];

  float value=0;
  unsigned long index = 0;
  unsigned long pointId = 0;
  unsigned long cellId = 0;
  unsigned long pointIndex = 0;

  while ( !std::isinf(value) && (index < nFloats) ) {

    value = data[index];

    if ( !( std::isnan(value) || std::isinf(value) ) ) {
      PointType pt;
      pt[0] = value;
      pt[1] = data[index+1];
      pt[2] = data[index+2];
      cellBuffer[pointIndex++] = pointId;
      outputMesh->SetPoint(pointId++,pt);
    }
    else if ( !std::isinf(value) ) {
        PolyLineCellType * polyline = new PolyLineCellType;
        polyline->SetPointIds( 0, pointIndex, cellBuffer );
        streamline.TakeOwnership( polyline );
        outputMesh->SetCell(cellId++, streamline);
        pointIndex = 0;
      }

    index += 3;
  }

  delete[] cellBuffer;

  if ( outputMesh->GetNumberOfCells() != count ) {
    itkExceptionMacro( "Number of streamlines read does not match header count" );
    return;
  }

}




template<class TOutputMesh>
void
MRTrixStreamlineFileReader<TOutputMesh>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "FileName: " << this->m_FileName << std::endl;
}

} //end of namespace itk

#endif
