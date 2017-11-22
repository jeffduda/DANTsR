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

  this->m_PointNormals = MatrixSetType::New();
  this->m_CellNormals = MatrixSetType::New();
  this->m_PointVectors = MatrixSetType::New();
  this->m_CellVectors = MatrixSetType::New();
  this->m_PointTextureCoordinates = MatrixSetType::New();
  this->m_CellTextureCoordinates = MatrixSetType::New();
  this->m_PointTensors = MatrixSetType::New();
  this->m_CellTensors = MatrixSetType::New();

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

  if( this->m_RandomPercentage < 1.0 )
    {
    typedef Statistics::MersenneTwisterRandomVariateGenerator GeneratorType;
    typename GeneratorType::Pointer generator = GeneratorType::New();
    generator->SetSeed();

    typename OutputMeshType::Pointer output = OutputMeshType::New();
    output->Initialize();

    if( this->GetOutput()->GetNumberOfPoints() > 0 )
      {
      typename OutputMeshType::PointsContainerIterator It =
        this->GetOutput()->GetPoints()->Begin();

      unsigned long count = 0;
      while( It != this->GetOutput()->GetPoints()->End() )
        {
        if( generator->GetVariateWithClosedRange() <= this->m_RandomPercentage )
          {
          output->SetPoint( count, It.Value() );
          PixelType label;
          bool elementExists = this->GetOutput()->GetPointData()
            ->GetElementIfIndexExists( It.Index(), &label );
          if( elementExists )
            {
            output->SetPointData( count, label );
            }
          count++;
          }
        ++It;
        }
      }
    this->GraftOutput( output );
    }
  this->m_LabelSet.clear();

  /**
   * If the number of points does not match the number of
   * point data, fill the point data with zeros
   */
  if( !this->GetOutput()->GetPointData() ||
    this->GetOutput()->GetPointData()->Size() !=
    this->GetOutput()->GetPoints()->Size() )
    {
    //itkWarningMacro( "Number of points does not match number of labels. "
    //  << "Filling point data with label zero." );
     typename OutputMeshType::PointsContainerIterator It =
       this->GetOutput()->GetPoints()->Begin();

     while( It != this->GetOutput()->GetPoints()->End() )
       {
       this->GetOutput()->SetPointData( It.Index(),
         NumericTraits<PixelType>::Zero );
       ++It;
       }
    }

  if( this->GetOutput()->GetNumberOfPoints() > 0 )
    {
    typename OutputMeshType::PointDataContainerIterator ItD =
      this->GetOutput()->GetPointData()->Begin();
    while( ItD != this->GetOutput()->GetPointData()->End() )
      {
      if( find( this->m_LabelSet.begin(),
        this->m_LabelSet.end(), ItD.Value() ) == this->m_LabelSet.end() )
        {
        this->m_LabelSet.push_back( ItD.Value() );
        }
      ++ItD;
      }
    }
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

  std::cout << "WARNING: No point data found with name: " << name << std::endl;

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
      Rcpp::Rcout << "Found POINT_DATA" << std::endl;
      std::string::size_type sp1 = line.find( " " );

      unsigned long nPoints = std::atoi( std::string( line, sp1, line.length()-1 ).c_str() );

      if (!this->ReadVTKData( nPoints, true ))
        {
        return false;
        }
      }
    else if (line.find("CELL_DATA") != std::string::npos)
      {
      Rcpp::Rcout << "Found CELL_DATA" << std::endl;
      std::string::size_type sp1 = line.find( " " );

      unsigned long nCells = std::atoi( std::string( line, sp1, line.length()-1 ).c_str() );

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
  Rcpp::Rcout << "Reading " << nPoints << " points" << std::endl;

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

  Rcpp::Rcout << "Reading " << nLines << " lines" << std::endl;
  typename OutputMeshType::Pointer outputMesh = this->GetOutput();
  this->m_Lines->Initialize();

  //std::cout << "Reading " << nLines << " lines with " << nValues << " values" << std::endl;


  if (!this->m_BinaryData)
  {
    for (unsigned long i=0; i<nLines; i++)
    {


      unsigned long lineSize;
      this->m_InputFile >> lineSize;

      LineType line(lineSize);

      for (unsigned long j=0; j<lineSize; j++)
      {
        unsigned long index;
        this->m_InputFile >> index;
        line[j] = index;
      }

      this->m_Lines->InsertElement(i,line);

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
    while (valueId < nValues)
    {
      unsigned int lineLength = lineData[valueId];
      ++valueId;

      LineType polyLine;
      polyLine.SetSize( lineLength );

      for (unsigned long i = 0; i < lineLength; i++)
        {
        polyLine[i] = lineData[valueId];
        ++valueId;
        }

      this->m_Lines->InsertElement( lineId, polyLine );
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
  Rcpp::Rcout << "Reading " << nPolygons << " polygons" << std::endl;
  typename OutputMeshType::Pointer outputMesh = this->GetOutput();
  this->m_Polygons->Initialize();

  if (!this->m_BinaryData)
  {
    for (unsigned long i=0; i<nPolygons; i++)
    {
      unsigned long lineSize;
      this->m_InputFile >> lineSize;

      LineType line(lineSize);

      for (unsigned long j=0; j<lineSize; j++)
      {
        unsigned long index;
        this->m_InputFile >> index;
        line[j] = index;
      }

      this->m_Polygons->InsertElement(i,line);

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
    while (valueId < nValues)
    {
      unsigned int lineLength = lineData[valueId];
      ++valueId;

      LineType polyLine;
      polyLine.SetSize( lineLength );

      for (unsigned long i = 0; i < lineLength; i++)
        {
        polyLine[i] = lineData[valueId];
        ++valueId;
        }

      this->m_Polygons->InsertElement( lineId, polyLine );
      ++lineId;
    }

    delete [] lineData;

  }

  return true;
}

template<class TOutputMesh>
bool
VtkPolyDataFileReader<TOutputMesh>
::ReadVTKVertices( unsigned long nVertices, unsigned long nValues )
{
  Rcpp::Rcout << "Reading " << nVertices << " vertices" << std::endl;
  typename OutputMeshType::Pointer outputMesh = this->GetOutput();
  this->m_Vertices->Initialize();

  if (!this->m_BinaryData)
    {
    for (unsigned long i=0; i<nVertices; i++)
      {
      unsigned long vertexSize;
      this->m_InputFile >> vertexSize;

      LineType vertex(vertexSize);

      for (unsigned long j=0; j<vertexSize; j++)
        {
        unsigned long index;
        this->m_InputFile >> index;
        vertex[j] = index;
        }
      this->m_Vertices->InsertElement(i,vertex);
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
    while (valueId < nValues)
      {
      unsigned int vertexLength = vertexData[valueId];
      ++valueId;

      LineType polyVertex;
      polyVertex.SetSize( vertexLength );

      for (unsigned long i = 0; i < vertexLength; i++)
        {
        polyVertex[i] = vertexData[valueId];
        ++valueId;
        }

      this->m_Vertices->InsertElement( lineId, polyVertex );
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
  Rcpp::Rcout << "Reading " << nStrips << " strips" << std::endl;
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
        if ( this->m_BinaryData) {
          dataType = "unsigned char";
        }
        MatrixType dat2 = this->ReadVTKDataMatrix(dataType, nValues, 4);

      }

      Rcpp::Rcout << "Read " << n << " SCALARS: isPointData=" << (int)isPointData << std::endl;
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
      Rcpp::Rcout << "Read " << n << " VECTORS: isPointData=" << (int)isPointData << std::endl;
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
      Rcpp::Rcout << "Read " << n << " NORMALS: isPointData=" << (int)isPointData << std::endl;
    }
    else if (line.find("TEXTURE_COORDINATES") != std::string::npos)
    {
      unsigned int dim = 1;
      std::string::size_type sp1 = line.find( " " );
      std::string::size_type sp2 = line.find( " ", sp1+1);
      std::string::size_type sp3 = line.find( " ", sp2+1);
      dim = std::atoi( std::string( line, sp2+1, (sp3-sp2)-1 ).c_str()  );
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
      Rcpp::Rcout << "Read " << n << " TEXTURE_COORDINATES: isPointData=" << (int)isPointData << std::endl;
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
      Rcpp::Rcout << "Read " << n << " TENSORS: isPointData=" << (int)isPointData << std::endl;
    }
    else if (line.find("FIELD") != std::string::npos)
    {
      Rcpp::Rcout << "FIELD not yet supported" << std::endl;
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
  Rcpp::Rcout << "ReadVTKDataMatrix() - ";
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
    else if ( !dataType.compare("int") ) {
      this->ReadVTKBinaryMatrix<int>( matrix );
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
