#ifndef __itkCaminoStreamlineFileWriter_hxx
#define __itkCaminoStreamlineFileWriter_hxx

#include "itkCaminoStreamlineFileWriter.h"
#include "itkImageFileWriter.h"
#include "itkByteSwapper.h"

#include <fstream>

namespace itk {

//
// Constructor
//
template<class TInputMesh>
CaminoStreamlineFileWriter<TInputMesh>
::CaminoStreamlineFileWriter()
{
  this->m_Input = nullptr;
  this->m_FileName = "";

  this->m_WriteFloat = true;
  this->m_Seeds = SeedSetType::New();
  this->m_Seeds->Initialize();

}

//
// Destructor
//
template<class TInputMesh>
CaminoStreamlineFileWriter<TInputMesh>
::~CaminoStreamlineFileWriter()
{
}

//
//
//
template<class TInputMesh>
void
CaminoStreamlineFileWriter<TInputMesh>
::SetInput(InputMeshType * input)
{
  this->m_Input = input;

}

template<class TInputMesh>
void
CaminoStreamlineFileWriter<TInputMesh>
::SetSeed(long i, unsigned long seed)
{
  //std::cout << "itk::CaminoStreamlineFileWrite::SetSeed()" << std::endl;
  this->m_Seeds->CreateIndex(i);
  this->m_Seeds->InsertElement(i, seed);
}

//
// Write the input mesh to the output file
//
template<class TInputMesh>
void CaminoStreamlineFileWriter<TInputMesh>
::Update()
{
  this->GenerateData();
}

//
// Write the input mesh to the output file
//
template<class TInputMesh>
void CaminoStreamlineFileWriter<TInputMesh>
::Write()
{
  this->GenerateData();
}

template<class TInputMesh>
void
CaminoStreamlineFileWriter<TInputMesh>
::GenerateData()
{
  //std::cout << "CaminoStreamlineFileWriter<TInputMesh>::GenerateData()"  << std::endl;

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
  if( extension == "Bfloat" )
    {
    m_WriteFloat = true;
    }
  else if( extension == "Bdouble" )
    {
    m_WriteFloat = false;
    }
  else
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

  this->WriteCaminoFile();
  m_OutputFile.close();

}
template<class TInputMesh>
void
CaminoStreamlineFileWriter<TInputMesh>
::WriteCaminoFile()
{
  //std::cout << "CaminoStreamlineFileWriter<TInputMesh>::WriteCaminoFile()"  << std::endl;
  for( unsigned int i=0; i<this->m_Input->GetNumberOfCells(); i++) {

    if ( m_WriteFloat ) {
      this->WriteCaminoTract<float>(i);
    }
    else {
      this->WriteCaminoTract<double>(i);
    }

  }
}

template<class TInputMesh>
template<typename PrecisionType>
void
CaminoStreamlineFileWriter<TInputMesh>
::WriteCaminoTract(unsigned int id)
{
  //
  // Write to output file
  //
  //std::cout << "CaminoStreamlineFileWriter<TInputMesh>::WriteCaminoTract<PrecisonType>()"  << std::endl;

  CellAutoPointer cell;
  if ( !this->m_Input->GetCell(id, cell) ) {
    itkExceptionMacro( "Cell could not be read" );
  }

  int nPoints = cell->GetNumberOfPoints();

  unsigned long dataSize = nPoints*3+2;
  PrecisionType * cellData = new PrecisionType[dataSize];

  unsigned long idx = 0;
  cellData[idx++] = static_cast<PrecisionType>(nPoints);
  cellData[idx++] = static_cast<PrecisionType>(m_Seeds->GetElement(id)); // FIXME: seed goes here


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

  itk::ByteSwapper<PrecisionType>::SwapRangeFromSystemToBigEndian(cellData, dataSize );
  m_OutputFile.write( reinterpret_cast<char *>( cellData ), dataSize * sizeof(PrecisionType) );

  delete[] cellData;

}

template<class TInputMesh>
void
CaminoStreamlineFileWriter<TInputMesh>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "FileName: " << this->m_FileName << std::endl;
}

} //end of namespace itk

#endif
