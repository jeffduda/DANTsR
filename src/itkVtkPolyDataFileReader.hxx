/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkVtkPolyDataFileReader.hxx,v $
  Language:  C++
  Date:      $Date: 2009/03/13 19:48:16 $
  Version:   $Revision: 1.23 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkVtkPolyDataFileReader_hxx
#define __itkVtkPolyDataFileReader_hxx

#include "itkVtkPolyDataFileReader.h"

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
VtkPolyDataFileReader<TOutputMesh>
::VtkPolyDataFileReader()
{
  this->m_RandomPercentage = 1.0;
  this->m_ExtractBoundaryPoints = false;

  this->m_MultiComponentScalarSets = NULL;
  this->m_CellMultiComponentScalarSets = NULL;
  this->m_Lines = LineSetType::New();
  this->m_Polygons = LineSetType::New();
  this->m_Vertices = LineSetType::New();
  this->m_Strips = LineSetType::New();

  this->m_PointScalars = MatrixSetType::New();
  this->m_CellScalars = MatrixSetType::New();
  this->m_PointNormals = MatrixSetType::New();
  this->m_CellNormals = MatrixSetType::New();
  this->m_PointVectors = MatrixSetType::New();
  this->m_CellVectors = MatrixSetType::New();
  this->m_PointTextureCoordinates = MatrixSetType::New();
  this->m_CellTextureCoordinates = MatrixSetType::New();
  this->m_PointTensors = MatrixSetType::New();
  this->m_CellTensors = MatrixSetType::New();

  this->m_PointScalarsNames = DataNameSetType::New();
  this->m_CellScalarsNames = DataNameSetType::New();
  this->m_PointNormalsNames = DataNameSetType::New();
  this->m_CellNormalsNames = DataNameSetType::New();
  this->m_PointVectorsNames = DataNameSetType::New();
  this->m_CellVectorsNames = DataNameSetType::New();
  this->m_PointTextureCoordinatesNames = DataNameSetType::New();
  this->m_CellTextureCoordinatesNames = DataNameSetType::New();
  this->m_PointTensorsNames = DataNameSetType::New();
  this->m_CellTensorsNames = DataNameSetType::New();

  //
  // Create the output
  //
  typename TOutputMesh::Pointer output = TOutputMesh::New();
  this->ProcessObject::SetNumberOfRequiredOutputs( 1 );
  this->ProcessObject::SetNthOutput( 0, output.GetPointer() );
}


template<class TOutputMesh>
void
VtkPolyDataFileReader<TOutputMesh>
::GenerateData()
{
  if( this->m_FileName == "" )
    {
    itkExceptionMacro( "No input FileName" );
    return;
    }

  //
  // Read input file
  //
  std::ifstream inputFile( m_FileName.c_str() );

  if( !inputFile.is_open() )
    {
    itkExceptionMacro("Unable to open file\n"
        "inputFilename= " << m_FileName );
    return;
    }
  else
    {
    inputFile.close();
    }

  /**
   * Get filename extension
   */

  std::string::size_type pos = this->m_FileName.rfind( "." );
  std::string extension( this->m_FileName, pos+1, this->m_FileName.length()-1 );

  if( extension == "vtk" )
    {
    this->ReadVTKFile();
    }
  else if ( extension == "vtp" )
    {
    std::cerr << "VTK XML files not yet supported" << std::endl;
    return;
    }
  else // try reading the file as an image
    {
    std::cerr << "Invalid file extension" << std::endl;
    return;
    }


  // FIXME Debugging output
  if ( false ) {
    Rcpp::Rcout << "# Points: " << this->GetOutput()->GetNumberOfPoints() << std::endl;
    Rcpp::Rcout << "# Cells: " << this->GetOutput()->GetNumberOfCells() << std::endl;

    Rcpp::Rcout << "# PointScalarSets: " << this->m_PointScalars->Size() << std::endl;
    for ( unsigned int i=0; i<this->m_PointScalars->Size(); i++ ) {
      Rcpp::Rcout << "  " << this->m_PointScalarsNames->ElementAt((long)i) << std::endl;
    }
    Rcpp::Rcout << "# CellScalarSets: " << this->m_CellScalars->Size() << std::endl;
    for ( unsigned int i=0; i<this->m_CellScalars->Size(); i++ ) {
      Rcpp::Rcout << "  " << this->m_CellScalarsNames->ElementAt((long)i) << std::endl;
    }
    Rcpp::Rcout << "# PointVectorSets: " << this->m_PointVectors->Size() << std::endl;
    for ( unsigned int i=0; i<this->m_PointVectors->Size(); i++ ) {
      Rcpp::Rcout << "  " << this->m_PointVectorsNames->ElementAt((long)i) << std::endl;
    }
    Rcpp::Rcout << "# CellVectorSets: " << this->m_CellVectors->Size() << std::endl;
    for ( unsigned int i=0; i<this->m_CellVectors->Size(); i++ ) {
      Rcpp::Rcout << "  " << this->m_CellVectorsNames->ElementAt((long)i) << std::endl;
    }
    Rcpp::Rcout << "# PointNormalsSets: " << this->m_PointNormals->Size() << std::endl;
    for ( unsigned int i=0; i<this->m_PointNormals->Size(); i++ ) {
      Rcpp::Rcout << "  " << this->m_PointNormalsNames->ElementAt((long)i) << std::endl;
    }
    Rcpp::Rcout << "# CellNormalsSets: " << this->m_CellNormals->Size() << std::endl;
    for ( unsigned int i=0; i<this->m_CellNormals->Size(); i++ ) {
      Rcpp::Rcout << "  " << this->m_CellNormalsNames->ElementAt((long)i) << std::endl;
    }
    Rcpp::Rcout << "# PointTextureCoordinatesSets: " << this->m_PointTextureCoordinates->Size() << std::endl;
    for ( unsigned int i=0; i<this->m_PointTextureCoordinates->Size(); i++ ) {
      Rcpp::Rcout << "  " << this->m_PointTextureCoordinatesNames->ElementAt((long)i) << std::endl;
    }
    Rcpp::Rcout << "# CellNormalsSets: " << this->m_CellTextureCoordinates->Size() << std::endl;
    for ( unsigned int i=0; i<this->m_CellTextureCoordinates->Size(); i++ ) {
      Rcpp::Rcout << "  " << this->m_CellTextureCoordinatesNames->ElementAt((long)i) << std::endl;
    }
    Rcpp::Rcout << "# PointTensorsSets: " << this->m_PointTensors->Size() << std::endl;
    for ( unsigned int i=0; i<this->m_PointTensors->Size(); i++ ) {
      Rcpp::Rcout << "  " << this->m_PointTensorsNames->ElementAt((long)i) << std::endl;
    }
    Rcpp::Rcout << "# CellTensorsSets: " << this->m_CellTensors->Size() << std::endl;
    for ( unsigned int i=0; i<this->m_CellTensors->Size(); i++ ) {
      Rcpp::Rcout << "  " << this->m_CellTensorsNames->ElementAt((long)i) << std::endl;
    }  }

}

template<class TOutputMesh>
void
VtkPolyDataFileReader<TOutputMesh>
::ReadVTKFile()
{
  this->ReadVTKPolyDataFile();
}


template<class TOutputMesh>
typename VtkPolyDataFileReader<TOutputMesh>::MultiComponentScalarSetType *
VtkPolyDataFileReader<TOutputMesh>
::GetMultiComponentScalarSetByName( const char * name )
{

  for (unsigned int i=0; i<this->m_MultiComponentScalarSetNames->Size(); i++)
  {
    if (!(this->m_MultiComponentScalarSetNames->GetElement(i).compare(name)))
    {
      return this->m_MultiComponentScalarSets->GetElement(i).GetPointer();
    }
  }

  //std::cout << "WARNING: No point data found with name: " << name << std::endl;

  return NULL;

}

template<class TOutputMesh>
bool
VtkPolyDataFileReader<TOutputMesh>
::ReadVTKPolyDataFile( )
{

  this->m_InputFile.open( this->m_FileName.c_str(), std::ifstream::in );
  if (!this->m_InputFile )
    {
    std::cerr << "Unable to open file: " << this->m_FileName << std::endl;
    return false;
    }

  std::string line;
  std::getline(this->m_InputFile, line);

  if (line.find( "# vtk DataFile Version" ) == std::string::npos)
    {
    std::cerr << "Invalid header: " << line << std::endl;
    return false;
    }

  // Get header info line (ignore for now)
  std::getline(this->m_InputFile, line);

  // Get ASCII || BINARY flag
  std::getline(this->m_InputFile, line);
  if (line.find("BINARY") != std::string::npos)
    {
    this->m_BinaryData = true;
    }
  else if (line.find("ASCII") != std::string::npos)
    {
    this->m_BinaryData = false;
    }
  else
    {
    std::cerr << "Unknown data format: " << line << std::endl;
    return false;
    }

  // Get Dataset type
  std::getline(this->m_InputFile, line);
  if (line.find("DATASET POLYDATA") != std::string::npos)
    {
    return this->ReadVTKPolyData( );
    }
  else
    {
    std::cerr << "Only POLYDATA currently supported" << std::endl;
    return false;
    }

  return false;
}

template<class TOutputMesh>
bool
VtkPolyDataFileReader<TOutputMesh>
::ReadVTKPolyData( )
{
  std::string line;

  while ( ! this->m_InputFile.eof() )
    {
    std::getline( this->m_InputFile, line );

    if ( line.find( "POINTS" ) != (std::string::npos) )
      {

      std::string::size_type pos = line.rfind( " " );
      std::string dataType = std::string( line, pos+1, line.length()-1 );

      pos = line.find( " " );
      std::string::size_type pos2 = line.find( " ", pos+1 );
      std::string temp = std::string( line, pos+1, (pos2-pos) );
      unsigned long nPoints = std::atoi( temp.c_str() );

      if (! this->ReadVTKPoints(nPoints, dataType))
        {
        return false;
        }
      }
    else if (line.find("LINES") != std::string::npos)
      {
      std::string::size_type sp1 = line.find( " " );
      std::string::size_type sp2 = line.find( " ", sp1+1);

      unsigned long nLines = std::atoi( std::string( line, sp1+1, (sp2-sp1) ).c_str() );
      unsigned long nValues = std::atoi( std::string( line, sp2+1, line.length()-sp2 ).c_str()  );

      if (!this->ReadVTKLines( nLines, nValues ))
        {
        return false;
        }
      }
    else if (line.find("POLYGONS") != std::string::npos)
      {
      std::string::size_type sp1 = line.find( " " );
      std::string::size_type sp2 = line.find( " ", sp1+1);

      unsigned long nPolys = std::atoi( std::string( line, sp1+1, (sp2-sp1) ).c_str() );
      unsigned long nValues = std::atoi( std::string( line, sp2+1, line.length()-sp2 ).c_str()  );

      if (!this->ReadVTKPolygons( nPolys, nValues ))
        {
        return false;
        }
      }
    else if (line.find("VERTICES") != std::string::npos)
      {
      std::string::size_type sp1 = line.find( " " );
      std::string::size_type sp2 = line.find( " ", sp1+1);

      unsigned long nVertices = std::atoi( std::string( line, sp1+1, (sp2-sp1) ).c_str() );
      unsigned long nValues = std::atoi( std::string( line, sp2+1, line.length()-sp2 ).c_str()  );

      if (!this->ReadVTKVertices( nVertices, nValues ))
        {
        return false;
        }
      }
    else if (line.find("TRIANGLE_STRIPS") != std::string::npos)
      {
      std::string::size_type sp1 = line.find( " " );
      std::string::size_type sp2 = line.find( " ", sp1+1);

      unsigned long nStrips = std::atoi( std::string( line, sp1+1, (sp2-sp1) ).c_str() );
      unsigned long nValues = std::atoi( std::string( line, sp2+1, line.length()-sp2 ).c_str()  );

      if (!this->ReadVTKStrips( nStrips, nValues ))
        {
        return false;
        }
      }
    else if (line.find("POINT_DATA") != std::string::npos)
      {
      //Rcpp::Rcout << "Found POINT_DATA" << std::endl;
      std::string::size_type sp1 = line.find( " " );

      unsigned long nPoints = std::atoi( std::string( line, sp1, line.length()-1 ).c_str() );

      if (!this->ReadVTKData( nPoints, true ))
        {
        return false;
        }
      }
    else if (line.find("CELL_DATA") != std::string::npos)
      {
      std::string::size_type sp1 = line.find( " " );
      unsigned long nCells = std::atoi( std::string( line, sp1, line.length()-1 ).c_str() );

      //Rcpp::Rcout << "Found CELL_DATA " << nCells << std::endl;

      if (!this->ReadVTKData( nCells, false ))
        {
        return false;
        }

      }
    else {
      //Rcpp::Rcout << "Skipping: " << line << std::endl;
    }

    }

  this->m_InputFile.close();
  return true;

}


template<class TOutputMesh>
bool
VtkPolyDataFileReader<TOutputMesh>
::ReadVTKPoints( unsigned long nPoints, std::string itkNotUsed(dataType) )
{
  //Rcpp::Rcout << "Reading " << nPoints << " points" << std::endl;

  typename OutputMeshType::Pointer outputMesh = this->GetOutput();
  outputMesh->GetPoints()->Initialize();
  outputMesh->GetPoints()->Reserve(nPoints);

  typename OutputMeshType::PointDataContainer::Pointer pointData = OutputMeshType::PointDataContainer::New();
  pointData->Initialize();
  outputMesh->SetPointData( pointData );

  //std::cout << "Reading " << nPoints << " points of type " << dataType << std::endl;

  if (!this->m_BinaryData)
  {
    for (unsigned long i=0; i<nPoints; i++)
    {
      float x,y,z;
      this->m_InputFile >> x >> y >> z;
      typename OutputMeshType::PointType point;
      point[0] = x;
      if (point.Size() > 1)
        point[1] = y;
      if (point.Size() > 2)
        point[2] = z;

      outputMesh->GetPoints()->InsertElement(i,point);

    }
  }
  else
  {
    itkDebugMacro( "Data is binary" );

    float p;
    float * ptData = new float [ nPoints*3 ];
    this->m_InputFile.read( reinterpret_cast< char * >( ptData ), 3 * nPoints * sizeof(p) );
    ByteSwapper<float>::SwapRangeFromSystemToBigEndian(ptData,nPoints*3);

    for (unsigned long i = 0; i < nPoints; i++ )
    {
      typename OutputMeshType::PointType point;
      for (long j = 0; j < 3; j++ )
      {
        if (j < OutputMeshType::PointType::Dimension)
          point[j] = ptData[i*3+j];
      }

      outputMesh->SetPoint( i, point );
    }

    delete [] ptData;

  }

  return true;
}

template<class TOutputMesh>
bool
VtkPolyDataFileReader<TOutputMesh>
::ReadVTKLines( unsigned long nLines, unsigned long nValues )
{
  //std::cout << "Reading " << nLines << " lines with " << nValues << " values" << std::endl;
  this->m_Lines->Initialize();
  typename OutputMeshType::Pointer outputMesh = this->GetOutput();

  CellAutoPointer streamline;
  //unsigned long cellId = 0;
  unsigned long index;
  unsigned long lineSize;


  if (!this->m_BinaryData)
  {
    for (unsigned long i=0; i<nLines; i++)
    {
      streamline.TakeOwnership( new PolyLineCellType );
      this->m_InputFile >> lineSize;
      for (unsigned long j=0; j<lineSize; j++)
      {
        this->m_InputFile >> index;
        streamline->SetPointId(j, index);
      }
      outputMesh->SetCell(i, streamline);
    }
  }
  else
  {
    int p;
    int * lineData = new int [ nValues ];
    this->m_InputFile.read( reinterpret_cast< char * >( lineData ), nValues * sizeof(p) );
    ByteSwapper<int>::SwapRangeFromSystemToBigEndian(lineData,nValues);

    unsigned long valueId = 0;
    unsigned long lineId = 0;

    CellAutoPointer streamline;

    while (valueId < nValues)
    {
      unsigned int lineLength = lineData[valueId];
      ++valueId;
      streamline.TakeOwnership( new PolyLineCellType );

      for (unsigned long i = 0; i < lineLength; i++)
        {
        streamline->SetPointId(i, lineData[valueId]);
        ++valueId;
        }

      outputMesh->SetCell(lineId, streamline);
      ++lineId;
    }
    delete [] lineData;
  }

  return true;
}

template<class TOutputMesh>
bool
VtkPolyDataFileReader<TOutputMesh>
::ReadVTKPolygons( unsigned long nPolygons, unsigned long nValues )
{
  //Rcpp::Rcout << "Reading " << nPolygons << " polygons" << std::endl;
  typename OutputMeshType::Pointer outputMesh = this->GetOutput();
  this->m_Polygons->Initialize();

  CellAutoPointer polyCell;

  unsigned long polySize;
  unsigned long index;

  if (!this->m_BinaryData)
  {
    for (unsigned long i=0; i<nPolygons; i++)
    {
      polyCell.TakeOwnership( new PolygonCellType );

      this->m_InputFile >> polySize;
      for (unsigned long j=0; j<polySize; j++)
      {
        this->m_InputFile >> index;
        polyCell->SetPointId(j, index);
      }
      outputMesh->SetCell(i, polyCell);
    }
  }
  else
  {
    int p;
    int * polyData = new int [ nValues ];
    this->m_InputFile.read( reinterpret_cast< char * >( polyData ), nValues * sizeof(p) );
    ByteSwapper<int>::SwapRangeFromSystemToBigEndian(polyData,nValues);

    unsigned long valueId = 0;
    unsigned long polyId = 0;

    CellAutoPointer polyCell;

    while (valueId < nValues)
    {
      unsigned int polySize = polyData[valueId];
      ++valueId;

      polyCell.TakeOwnership( new PolygonCellType );
      for (unsigned long i = 0; i < polySize; i++)
        {
        polyCell->SetPointId(i, polyData[valueId]);
        ++valueId;
        }
      outputMesh->SetCell(polyId, polyCell);
      ++polyId;
    }

  }

  return true;
}

template<class TOutputMesh>
bool
VtkPolyDataFileReader<TOutputMesh>
::ReadVTKVertices( unsigned long nVertices, unsigned long nValues )
{
  //Rcpp::Rcout << "Reading " << nVertices << " vertices" << std::endl;
  typename OutputMeshType::Pointer outputMesh = this->GetOutput();
  this->m_Vertices->Initialize();

  CellAutoPointer vertexCell;

  unsigned long vertexSize;
  unsigned long index;


  if (!this->m_BinaryData)
    {
    for (unsigned long i=0; i<nVertices; i++)
      {
      this->m_InputFile >> vertexSize;

      //typename OutputMeshType::PointIdentifier polyPoints[ vertexSize ];
      VertexCellType * vertex = new VertexCellType;

      if ( vertexSize != 1 ) {
        itkExceptionMacro("Only single-point vertices are supported");
        return false;
      }

      this->m_InputFile >> index;
      vertex->SetPointId(0, index);

      vertexCell.TakeOwnership( vertex );
      outputMesh->SetCell(i, vertexCell);
      }
    }
  else
    {

    int p;
    int * vertexData = new int [ nValues ];
    this->m_InputFile.read( reinterpret_cast< char * >( vertexData ), nValues * sizeof(p) );
    ByteSwapper<int>::SwapRangeFromSystemToBigEndian(vertexData,nValues);

    unsigned long valueId = 0;
    unsigned long lineId = 0;

    CellAutoPointer vertexCell;

    unsigned int vertexLength;

    while (valueId < nValues)
      {
      vertexLength = vertexData[valueId];
      ++valueId;

      if ( vertexLength != 1 ) {
        itkExceptionMacro("Only single-point vertices are supported");
        return false;
      }
      VertexCellType * vertex = new VertexCellType;
      vertex->SetPointId(0, vertexData[valueId]);
      ++valueId;
      vertexCell.TakeOwnership( vertex );

      outputMesh->SetCell(lineId, vertexCell);
      ++lineId;
      }
    delete [] vertexData;
  }

  return true;
}

template<class TOutputMesh>
bool
VtkPolyDataFileReader<TOutputMesh>
::ReadVTKStrips( unsigned long nStrips, unsigned long nValues )
{
  //Rcpp::Rcout << "Reading " << nStrips << " strips" << std::endl;
  typename OutputMeshType::Pointer outputMesh = this->GetOutput();
  this->m_Vertices->Initialize();

  if (!this->m_BinaryData)
    {
    for (unsigned long i=0; i<nStrips; i++)
      {
      unsigned long stripSize;
      this->m_InputFile >> stripSize;

      LineType strip(stripSize);

      for (unsigned long j=0; j<stripSize; j++)
        {
        unsigned long index;
        this->m_InputFile >> index;
        strip[j] = index;
        }
      this->m_Strips->InsertElement(i,strip);
      }
    }
  else
    {

    int p;
    int * stripData = new int [ nValues ];
    this->m_InputFile.read( reinterpret_cast< char * >( stripData ), nValues * sizeof(p) );
    ByteSwapper<int>::SwapRangeFromSystemToBigEndian(stripData,nValues);

    unsigned long valueId = 0;
    unsigned long lineId = 0;
    while (valueId < nValues)
      {
      unsigned int stripLength = stripData[valueId];
      ++valueId;

      LineType polyStrip;
      polyStrip.SetSize( stripLength );

      for (unsigned long i=0; i < stripLength; i++)
        {
        polyStrip[i] = stripData[valueId];
        ++valueId;
        }

      this->m_Strips->InsertElement( lineId, polyStrip );
      ++lineId;
      }
    delete [] stripData;
  }

  return true;
}

template<class TOutputMesh>
bool
VtkPolyDataFileReader<TOutputMesh>
::ReadVTKData( unsigned long n, bool isPointData )
{

  std::string line;
  std::getline(this->m_InputFile, line);

  // When a POINT_DATA or CELl_DATA line is included with no data after it
  if ( line.find("POINT_DATA") != std::string::npos) {
    std::string::size_type sp1 = line.find( " " );
    unsigned long nPoints = std::atoi( std::string( line, sp1, line.length()-1 ).c_str() );
    //Rcpp::Rcout << "Found POINT_DATA " << nPoints << std::endl;
    return this->ReadVTKData(nPoints, true);

  }
  else if (line.find("CELL_DATA") != std::string::npos) {
    std::string::size_type sp1 = line.find( " " );
    unsigned long nCells = std::atoi( std::string( line, sp1, line.length()-1 ).c_str() );
    //Rcpp::Rcout << "Found CELL_DATA " << nCells << std::endl;
    return this->ReadVTKData(nCells, false);
  }

  while (line.length() > 0)
  {
    if (line.find("SCALARS") != std::string::npos)
    {
      unsigned long nComponents = 1;
      std::string::size_type sp1 = line.find( " " );
      std::string::size_type sp2 = line.find( " ", sp1+1);
      std::string::size_type sp3 = line.find( " ", sp2+1);
      if (sp3 == std::string::npos)
      {
        sp3 = line.length()-1;
      }
      else
      {
        nComponents = std::atoi( std::string( line, sp3+1, line.length()-sp2-1 ).c_str()  );
      }

      std::string dataName = std::string( line, sp1+1, (sp2-sp1)-1 );
      std::string dataType = std::string( line, sp2+1, (sp3-sp2)-1 );

      std::getline(this->m_InputFile, line);
      sp1 = line.find( " " );
      std::string tableName = std::string( line, sp1+1, line.length()-1 );

      MatrixType dat = this->ReadVTKDataMatrix(dataType, n, nComponents);

      if ( tableName.compare("default") ) {
        std::getline(this->m_InputFile, line);
        sp1 = line.find( " " );
        sp2 = line.find( " ", sp1+1);
        sp3 = line.find( " ", sp2+1);
        std::string tableValuesName = std::string( line, sp1+1, line.length()-1 );
        unsigned long nValues = std::atoi( std::string( line, sp3+1, line.length()-sp2-1 ).c_str()  );

        // unsigned char for binary, float for ASCII
        if ( this->m_BinaryData ) {
          dataType = "unsigned char";
        }
        MatrixType dat2 = this->ReadVTKDataMatrix(dataType, nValues, 4);

      }
      if (isPointData) {
        long idx = this->m_PointScalars->Size();
        this->m_PointScalars->CreateIndex(idx);
        this->m_PointScalars->InsertElement(idx, dat);
        this->m_PointScalarsNames->CreateIndex(idx);
        this->m_PointScalarsNames->InsertElement(idx, dataName);
      }
      else {
        long idx = this->m_CellScalars->Size();
        this->m_CellScalars->CreateIndex(idx);
        this->m_CellScalars->InsertElement(idx, dat);
        this->m_CellScalarsNames->CreateIndex(idx);
        this->m_CellScalarsNames->InsertElement(idx, dataName);
      }
      //Rcpp::Rcout << "Read " << n << " SCALARS named " << dataName << " isPointData=" << (int)isPointData << std::endl;
    }
    else if (line.find("COLOR_SCALARS") != std::string::npos)
    {
      Rcpp::Rcout << "COLOR_SCALARS not yet supported" << std::endl;
    }
    else if (line.find("LOOKUP_TABLE") != std::string::npos)
    {
      Rcpp::Rcout << "LOOKUP_TABLE not yet supported" << std::endl;
    }
    else if (line.find("VECTORS") != std::string::npos)
    {
      std::string::size_type sp1 = line.find( " " );
      std::string::size_type sp2 = line.find( " ", sp1+1);
      std::string dataName = std::string( line, sp1+1, (sp2-sp1)-1 );
      std::string dataType = std::string( line, sp2+1, line.length()-sp1-1 );
      MatrixType dat = this->ReadVTKDataMatrix(dataType, n, 3);
      if ( isPointData) {
        this->m_PointVectors->InsertElement(this->m_PointVectors->Size(), dat);
        this->m_PointVectorsNames->InsertElement(this->m_PointVectorsNames->Size(), dataName);
      }
      else {
        this->m_CellVectors->InsertElement(this->m_CellVectors->Size(), dat);
        this->m_CellVectorsNames->InsertElement(this->m_CellVectorsNames->Size(), dataName);
      }
      //Rcpp::Rcout << "Read " << n << " VECTORS: isPointData=" << (int)isPointData << std::endl;
    }
    else if (line.find("NORMALS") != std::string::npos)
    {
      std::string::size_type sp1 = line.find( " " );
      std::string::size_type sp2 = line.find( " ", sp1+1);
      std::string dataName = std::string( line, sp1+1, (sp2-sp1)-1 );
      std::string dataType = std::string( line, sp2+1, line.length()-sp1-1 );
      MatrixType dat = this->ReadVTKDataMatrix(dataType, n, 3);
      if ( isPointData) {
        this->m_PointNormals->InsertElement(this->m_PointNormals->Size(), dat);
        this->m_PointNormalsNames->InsertElement(this->m_PointNormalsNames->Size(), dataName);
      }
      else {
        this->m_CellNormals->InsertElement(this->m_CellNormals->Size(), dat);
        this->m_CellNormalsNames->InsertElement(this->m_CellNormalsNames->Size(), dataName);
      }
      //Rcpp::Rcout << "Read " << n << " NORMALS: isPointData=" << (int)isPointData << std::endl;
    }
    else if (line.find("TEXTURE_COORDINATES") != std::string::npos)
    {
      unsigned int dim = 1;
      std::string::size_type sp1 = line.find( " " );
      std::string::size_type sp2 = line.find( " ", sp1+1);
      std::string::size_type sp3 = line.find( " ", sp2+1);
      dim = std::atoi( std::string( line, sp2+1, (sp3-sp2)-1 ).c_str() );
      std::string dataName = std::string( line, sp1+1, (sp2-sp1)-1 );
      std::string dataType = std::string( line, sp3+1, line.length()-sp2-1 );
      MatrixType dat = this->ReadVTKDataMatrix(dataType, n, dim);
      if ( isPointData) {
        this->m_PointTextureCoordinates->InsertElement(this->m_PointTextureCoordinates->Size(), dat);
        this->m_PointTextureCoordinatesNames->InsertElement(this->m_PointTextureCoordinatesNames->Size(), dataName);
      }
      else {
        this->m_CellTextureCoordinates->InsertElement(this->m_CellTextureCoordinates->Size(), dat);
        this->m_CellTextureCoordinatesNames->InsertElement(this->m_CellTextureCoordinatesNames->Size(), dataName);
      }
      //Rcpp::Rcout << "Read " << n << " TEXTURE_COORDINATES: isPointData=" << (int)isPointData << std::endl;
    }
    else if (line.find("TENSORS") != std::string::npos)
    {
      std::string::size_type sp1 = line.find( " " );
      std::string::size_type sp2 = line.find( " ", sp1+1);
      std::string dataName = std::string( line, sp1+1, (sp2-sp1)-1 );
      std::string dataType = std::string( line, sp2+1, line.length()-sp1-1 );
      MatrixType dat = this->ReadVTKDataMatrix(dataType, n, 9);
      if ( isPointData) {
        this->m_PointTensors->InsertElement(this->m_PointTensors->Size(), dat);
        this->m_PointTensorsNames->InsertElement(this->m_PointTensorsNames->Size(), dataName);
      }
      else {
        this->m_CellTensors->InsertElement(this->m_CellTensors->Size(), dat);
        this->m_CellTensorsNames->InsertElement(this->m_CellTensorsNames->Size(), dataName);
      }
      //Rcpp::Rcout << "Read " << n << " TENSORS: isPointData=" << (int)isPointData << std::endl;
    }
    else if (line.find("FIELD") != std::string::npos)
    {
      //Rcpp::Rcout << "FIELD not yet fully supported" << std::endl;
      std::string::size_type sp1 = line.find( " " );
      std::string::size_type sp2 = line.find( " ", sp1+1);
      std::string dataName = std::string( line, sp1+1, (sp2-sp1)-1 );
      unsigned int nFields = std::atoi( std::string( line, sp2+1, line.length()-sp1-1 ).c_str() );
      //Rcpp::Rcout << "FIELD named " << dataName << " has " << nFields << " data arrays" << std::endl;
      for (unsigned int i=0; i<nFields; i++) {
        std::getline(this->m_InputFile, line);
        sp1 = line.find( " " );
        sp2 = line.find( " ", sp1+1);
        std::string::size_type sp3 = line.find( " ", sp2+1);

        std::string dataName = std::string( line, 0, sp1-1 );
        unsigned int nComponents = std::atoi( std::string( line, sp1+1, (sp2-sp1)-1 ).c_str() );
        unsigned int nTuples = std::atoi( std::string( line, sp2+1, (sp3-sp2)-1 ).c_str() );
        std::string dataType = std::string( line, sp3+1, line.length()-sp3-1 );

        MatrixType dat = this->ReadVTKDataMatrix(dataType, nComponents, nTuples);
        //Rcpp::Rcout << "Read array " << dataName << " : " << nComponents << "x" << nTuples << " values of type " << dataType << std::endl;
      }
    }

    std::getline( this->m_InputFile, line );
    //Rcpp::Rcout << "Next line=" << line << std::endl;
  }

  return true;
}


template<class TOutputMesh>
bool
VtkPolyDataFileReader<TOutputMesh>
::ReadVTKScalars( std::string dataName, std::string itkNotUsed(dataType), unsigned long nPoints, unsigned long nComponents, bool isPointData )
{

  //std::cout << "Reading " << nPoints << " scalar data points with " << nComponents << " components of type " << dataType << " named " << dataName << std::endl;

  // get LOOKUP_TABLE name
  std::string line;
  std::getline(this->m_InputFile, line);

  bool labels = false;
  if (dataName.find("pointLabels") != std::string::npos)
    {
    labels = true;
    }

  typename MultiComponentScalarSetType::Pointer set = MultiComponentScalarSetType::New();
  set->Initialize();
  set->Reserve( nPoints );

  if (!this->m_BinaryData)
  {
    for (unsigned long i=0; i<nPoints; i++)
    {
      MultiComponentScalarType value(nComponents);
      for (unsigned long j=0; j<nComponents; j++)
      {
        float component;
        this->m_InputFile >> component;
        value[j] = component;
      }
      if (labels)
      {
        typename OutputMeshType::PixelType label = static_cast<typename OutputMeshType::PixelType>( value[0] );
        this->GetOutput()->GetPointData()->InsertElement(i, label);
      }
      else
      {
        set->InsertElement(i, value);
      }
    }
  }
  else
  {

  }

  if (!labels)
  {
    if (isPointData)
      {
      long id = this->m_MultiComponentScalarSets->Size();
      this->m_MultiComponentScalarSets->InsertElement( id, set );
      this->m_MultiComponentScalarSetNames->InsertElement( id, dataName );
    }
    else
      {
      long id = this->m_CellMultiComponentScalarSets->Size();
      this->m_CellMultiComponentScalarSets->InsertElement( id, set );
      this->m_CellMultiComponentScalarSetNames->InsertElement( id, dataName );
      }
  }
  return true;
}

template<class TOutputMesh>
typename VtkPolyDataFileReader<TOutputMesh>::MatrixType
VtkPolyDataFileReader<TOutputMesh>
::ReadVTKDataMatrix( std::string dataType, unsigned long nRows, unsigned long nCols )
{
  //Rcpp::Rcout << "ReadVTKDataMatrix(" << nRows << "," << nCols << ") - ";
  MatrixType matrix(nRows, nCols);

  if (!this->m_BinaryData)
  {
    this->ReadVTKASCIIMatrix( matrix );
  }
  else
  {
    if ( !dataType.compare("unsigned char") ) {
      this->ReadVTKBinaryMatrix<unsigned char>( matrix );
    }
    else if ( !dataType.compare("char") ) {
      this->ReadVTKBinaryMatrix<char>( matrix );
    }
    else if ( !dataType.compare("int") ) {
      this->ReadVTKBinaryMatrix<int>( matrix );
    }
    else if ( !dataType.compare("long") ) {
      this->ReadVTKBinaryMatrix<long>( matrix );
    }
    else if ( !dataType.compare("float") ) {
      this->ReadVTKBinaryMatrix<float>( matrix );
    }
    else if ( !dataType.compare("double") ) {
      this->ReadVTKBinaryMatrix<double>( matrix );
    }
    else {
      Rcpp::Rcout << "Datatype not currently supported: " << dataType << std::endl;
    }
  }

  return matrix;
}






template<class TOutputMesh>
void
VtkPolyDataFileReader<TOutputMesh>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "FileName: "
     << this->m_FileName << std::endl;
  os << indent << "RandomPercentage: "
     << this->m_RandomPercentage << std::endl;

}

} //end of namespace itk


#endif
