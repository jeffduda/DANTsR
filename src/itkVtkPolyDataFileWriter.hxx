/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkVtkPolyDataFileWriter.hxx,v $
  Language:  C++
  Date:      $Date: 2009/03/04 23:10:58 $
  Version:   $Revision: 1.18 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkVtkPolyDataFileWriter_hxx
#define __itkVtkPolyDataFileWriter_hxx

#include "itkVtkPolyDataFileWriter.h"
#include "itkImageFileWriter.h"
#include "itkByteSwapper.h"

#include <fstream>

namespace itk {

//
// Constructor
//
template<class TInputMesh>
VtkPolyDataFileWriter<TInputMesh>
::VtkPolyDataFileWriter()
{
  this->m_WriteBinary = true;
  this->m_Input = NULL;
  this->m_FileName = "";
  this->m_MultiComponentScalarSets = NULL;
  this->m_Lines = NULL;
  this->m_ImageSize.Fill( 0 );
  this->m_CellsAsPolygons = false;
  this->m_CellsAsLines = false;
}

//
// Destructor
//
template<class TInputMesh>
VtkPolyDataFileWriter<TInputMesh>
::~VtkPolyDataFileWriter()
{
}

//
// Set the input mesh
//
template<class TInputMesh>
void
VtkPolyDataFileWriter<TInputMesh>
::SetCellsAsPolygons(bool flag)
{
  this->m_CellsAsPolygons = flag;
  if ( flag ) {
    this->m_CellsAsLines = false;
  }
}

//
template<class TInputMesh>
void
VtkPolyDataFileWriter<TInputMesh>
::SetCellsAsLines(bool flag)
{
  this->m_CellsAsLines = flag;
  if ( flag ) {
    this->m_CellsAsPolygons = false;
  }
}

//
//
//
template<class TInputMesh>
void
VtkPolyDataFileWriter<TInputMesh>
::SetInput(InputMeshType * input)
{
  this->m_Input = input;
}

//
// Write the input mesh to the output file
//
template<class TInputMesh>
void VtkPolyDataFileWriter<TInputMesh>
::Update()
{
  this->GenerateData();
}

//
// Write the input mesh to the output file
//
template<class TInputMesh>
void VtkPolyDataFileWriter<TInputMesh>
::Write()
{
  this->GenerateData();
}

template<class TInputMesh>
void
VtkPolyDataFileWriter<TInputMesh>
::GenerateData()
{
  if( this->m_FileName == "" )
    {
    itkExceptionMacro( "No FileName" );
    return;
    }

  if( this->m_ImageSize[0] == 0 )
    {
    this->m_ImageSize.Fill( 100 );
    this->m_ImageOrigin.CastFrom(
      this->m_Input->GetBoundingBox()->GetMinimum() );
    for( unsigned int d = 0; d < Dimension; d++ )
      {
      this->m_ImageSpacing[d] = (
        this->m_Input->GetBoundingBox()->GetMaximum()[d] -
        this->m_Input->GetBoundingBox()->GetMinimum()[d] )
        / static_cast<double>( this->m_ImageSize[d] + 1 );
      }
    this->m_ImageDirection.SetIdentity();
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

  if( extension == "txt" )
    {
    this->WritePointsToAvantsFile();
    }
  else if( extension == "vtk" )
    {
    this->WriteVTKFile();
    }
  else
    {
    try
      {
      this->WritePointsToImageFile();
      }
    catch(...)
      {
      itkExceptionMacro( "Unknown extension: " << extension );
      }
    }

}
template<class TInputMesh>
void
VtkPolyDataFileWriter<TInputMesh>
::WriteVTKFile()
{
  this->WritePointsToVTKFile();
  this->WriteLinesToVTKFile();
  this->WriteVTKPolygons();
  this->WritePointDataToVTKFile();
}

template<class TInputMesh>
void
VtkPolyDataFileWriter<TInputMesh>
::WritePointsToVTKFile()
{
  //
  // Write to output file
  //
  std::ofstream outputFile( this->m_FileName.c_str() );

  outputFile << "# vtk DataFile Version 2.0" << std::endl;
  outputFile << "File written by itkVtkPolyDataFileWriter" << std::endl;
  if (this->m_WriteBinary)
    {
    outputFile << "BINARY" << std::endl;
    }
  else
    {
    outputFile << "ASCII" << std::endl;
    }
  outputFile << "DATASET POLYDATA" << std::endl;

  // POINTS go first

  unsigned int numberOfPoints = this->m_Input->GetNumberOfPoints();
  outputFile << "POINTS " << numberOfPoints << " float" << std::endl;

  typename InputMeshType::PointsContainerIterator pointIterator
    = this->m_Input->GetPoints()->Begin();
  typename InputMeshType::PointsContainerIterator pointEnd
    = this->m_Input->GetPoints()->End();

  if ( this->m_WriteBinary )
    {

    itkDebugMacro( "Writing data as binary" );
    float p;
    float * ptData = new float [ numberOfPoints*3 ];
    //inputFile.read( reinterpret_cast< char * >( ptData ), 3 * numberOfPoints * sizeof(p) );
    //ByteSwapper<float>::SwapRangeFromSystemToBigEndian(ptData,numberOfPoints*3);

    unsigned long idx = 0;
    while ( pointIterator !=  pointEnd )
      {
      ptData[idx] = pointIterator.Value()[0];
      ++idx;
      ptData[idx] = pointIterator.Value()[1];
      ++idx;
      ptData[idx] = pointIterator.Value()[2];
      ++idx;
      pointIterator++;
      }
    ByteSwapper<float>::SwapRangeFromSystemToBigEndian(ptData,numberOfPoints*3);
    outputFile.write( reinterpret_cast<char *>( ptData ), 3 * numberOfPoints * sizeof(p) );

    delete [] ptData;

    }
  else
    {
    while( pointIterator != pointEnd )
      {
      PointType point = pointIterator.Value();
      outputFile << point[0] << " " << point[1];
      if( Dimension == 2 )
        {
        outputFile << " 0 " << std::endl;
        }
      else if( Dimension == 3 )
        {
        outputFile << " " << point[2] << " " << std::endl;
        }
      pointIterator++;
      }
    }
  outputFile.close();
}

template<class TInputMesh>
void
VtkPolyDataFileWriter<TInputMesh>
::WritePointDataToVTKFile()
{
  //
  // Write to output file
  //
  std::ofstream outputFile( this->m_FileName.c_str(), std::ios::app );

  // No point data conditions
  //if (!this->m_Input->GetPointData()) return;
  //if (this->m_Input->GetPointData()->Size() == 0) return;

  unsigned int numberOfPoints = this->m_Input->GetNumberOfPoints();

  if( !this->m_MultiComponentScalarSets )
    {

    if (!this->m_Input->GetPointData())
    {
      outputFile.close();
      return;
    }

    if (this->m_Input->GetPointData()->Size() == 0)
    {
      outputFile.close();
      return;
    }

    outputFile << std::endl;
    outputFile << "POINT_DATA " << numberOfPoints << std::endl;
    std::string type = std::string( "float" );

    outputFile << "SCALARS pointLabels " << type
      << " 1" << std::endl;
    outputFile << "LOOKUP_TABLE default" << std::endl;

    typename InputMeshType::PointDataContainerIterator pointDataIterator
      = this->m_Input->GetPointData()->Begin();
    typename InputMeshType::PointDataContainerIterator pointDataEnd
      = this->m_Input->GetPointData()->End();

    while( pointDataIterator != pointDataEnd )
      {
      outputFile << pointDataIterator.Value() << " ";
      pointDataIterator++;
      }
    outputFile << std::endl;
    }
  else
    {
    outputFile << std::endl;
    outputFile << "POINT_DATA " << numberOfPoints << std::endl;
    std::string type = std::string( "float" );

    for (unsigned int nSet = 0; nSet < this->m_MultiComponentScalarSets->Size(); nSet++)
      {
      MultiComponentScalarType scalar
        = this->m_MultiComponentScalarSets->GetElement(nSet)->GetElement( 0 );
      unsigned int numberOfComponents = scalar.GetSize();

      std::ostringstream defaultName;
      defaultName << "scalars" << nSet;
      std::string dataName = defaultName.str();
      if (this->m_MultiComponentScalarSetNames)
      {
        if (nSet < this->m_MultiComponentScalarSetNames->Size())
        {
          dataName = this->m_MultiComponentScalarSetNames->GetElement(nSet);
        }
      }


      outputFile << "SCALARS " << dataName << " " << type
        << " " << numberOfComponents << std::endl;
      outputFile << "LOOKUP_TABLE default" << std::endl;

      typename MultiComponentScalarSetType::Iterator It
        = this->m_MultiComponentScalarSets->GetElement(nSet)->Begin();
      typename MultiComponentScalarSetType::Iterator ItEnd
        = this->m_MultiComponentScalarSets->GetElement(nSet)->End();

      while( It != ItEnd )
      {
        for (unsigned int i=0; i<It.Value().Size(); i++)
        {
          outputFile << It.Value()[i] << " ";
        }
        It++;
      }

      outputFile << std::endl;
      }

    }
  outputFile.close();
}


template<class TInputMesh>
void
VtkPolyDataFileWriter<TInputMesh>
::WriteCellDataToVTKFile()
{
  //
  // Write to output file
  //
  std::ofstream outputFile( this->m_FileName.c_str(), std::ios::app );

  // No point data conditions
  //if (!this->m_Input->GetPointData()) return;
  //if (this->m_Input->GetPointData()->Size() == 0) return;

  unsigned int numberOfPoints = this->m_Input->GetNumberOfPoints();

  if( !this->m_MultiComponentScalarSets )
    {

    if (!this->m_Input->GetPointData())
    {
      outputFile.close();
      return;
    }

    if (this->m_Input->GetPointData()->Size() == 0)
    {
      outputFile.close();
      return;
    }

    outputFile << std::endl;
    outputFile << "CELL_DATA " << numberOfPoints << std::endl;
    std::string type = std::string( "float" );

    outputFile << "SCALARS pointLabels " << type
      << " 1" << std::endl;
    outputFile << "LOOKUP_TABLE default" << std::endl;

    typename InputMeshType::PointDataContainerIterator pointDataIterator
      = this->m_Input->GetPointData()->Begin();
    typename InputMeshType::PointDataContainerIterator pointDataEnd
      = this->m_Input->GetPointData()->End();

    while( pointDataIterator != pointDataEnd )
      {
      outputFile << pointDataIterator.Value() << " ";
      pointDataIterator++;
      }
    outputFile << std::endl;
    }
  else
    {
    outputFile << std::endl;
    outputFile << "CELL_DATA " << numberOfPoints << std::endl;
    std::string type = std::string( "float" );

    for (unsigned int nSet = 0; nSet < this->m_MultiComponentScalarSets->Size(); nSet++)
      {
      MultiComponentScalarType scalar
        = this->m_MultiComponentScalarSets->GetElement(nSet)->GetElement( 0 );
      unsigned int numberOfComponents = scalar.GetSize();

      std::ostringstream defaultName;
      defaultName << "scalars" << nSet;
      std::string dataName = defaultName.str();
      if (this->m_MultiComponentScalarSetNames)
      {
        if (nSet < this->m_MultiComponentScalarSetNames->Size())
        {
          dataName = this->m_MultiComponentScalarSetNames->GetElement(nSet);
        }
      }


      outputFile << "SCALARS " << dataName << " " << type
        << " " << numberOfComponents << std::endl;
      outputFile << "LOOKUP_TABLE default" << std::endl;

      typename MultiComponentScalarSetType::Iterator It
        = this->m_MultiComponentScalarSets->GetElement(nSet)->Begin();
      typename MultiComponentScalarSetType::Iterator ItEnd
        = this->m_MultiComponentScalarSets->GetElement(nSet)->End();

      while( It != ItEnd )
      {
        for (unsigned int i=0; i<It.Value().Size(); i++)
        {
          outputFile << It.Value()[i] << " ";
        }
        It++;
      }

      outputFile << std::endl;
      }

    }
  outputFile.close();
}


template<class TInputMesh>
void
VtkPolyDataFileWriter<TInputMesh>
::WriteLinesToVTKFile()
{


  if ( this->m_CellsAsLines ) {
    std::ofstream outputFile( this->m_FileName.c_str(), std::ios::app );

    unsigned long numberOfLines = this->m_Input->GetNumberOfCells();
    unsigned long totalSize = 0;

    CellAutoPointer cell;
    for ( unsigned long i=0; i<numberOfLines; i++ ) {
      if ( !this->m_Input->GetCell(i, cell) ) {
        Rcpp::stop("Cell could not be read");
      }
      totalSize += cell->GetNumberOfPoints() + 1;
    }

    // Debugging
    for ( unsigned long i=0; i<numberOfLines; i++)  {
      this->m_Input->GetCell(i, cell);

      Rcpp::Rcout << "Cell " << i << ": " << cell->GetNumberOfPoints() << " ";
      for ( unsigned long j=0; j<cell->GetNumberOfPoints(); j++) {
        Rcpp::Rcout << cell->GetPointIds()[j] << " ";
      }
      Rcpp::Rcout << std::endl;
    }

    outputFile << "LINES " <<
      numberOfLines << " " << totalSize << std::endl;

    if ( this->m_WriteBinary )
      {
      int p;
      int * ptData = new int [ totalSize ];

      unsigned long idx = 0;
      for ( unsigned long i=0; i<numberOfLines; i++)  {
        this->m_Input->GetCell(i, cell);
        ptData[idx] = cell->GetNumberOfPoints();
        ++idx;

        for ( unsigned long j=0; j<cell->GetNumberOfPoints(); j++) {
          ptData[idx] = cell->GetPointIds()[j];
          ++idx;
        }
      }

      ByteSwapper<int>::SwapRangeFromSystemToBigEndian(ptData,totalSize);
      outputFile.write( reinterpret_cast<char *>( ptData ), totalSize * sizeof(p) );
      delete [] ptData;
      }
    else
      {
      for ( unsigned long i=0; i<numberOfLines; i++)  {
        this->m_Input->GetCell(i, cell);

        outputFile << cell->GetNumberOfPoints() << " ";
        for ( unsigned long j=0; j<cell->GetNumberOfPoints(); j++) {
          outputFile << cell->GetPointIds()[j] << " ";
        }

        outputFile << std::endl;
      }

    outputFile << std::endl;
    outputFile.close();
    }

  }
  else if( this->m_Lines )
    {

    std::ofstream outputFile( this->m_FileName.c_str(), std::ios::app );

    unsigned int numberOfLines = this->m_Lines->Size();
    unsigned int totalSize = 0;

    typename LineSetType::Iterator It
      = this->m_Lines->Begin();
    typename LineSetType::Iterator ItEnd = this->m_Lines->End();

    while( It != ItEnd )
      {
      totalSize += ( It.Value() ).Size();
      totalSize++;
      It++;
      }

    outputFile << "LINES " <<
      numberOfLines << " " << totalSize << std::endl;

    It = this->m_Lines->Begin();

    if ( this->m_WriteBinary )
      {
      int p;
      int * ptData = new int [ totalSize ];

      unsigned long idx = 0;
      while( It != ItEnd )
        {
        ptData[idx] = It.Value().Size();
        ++idx;
        for (unsigned int i=0; i<It.Value().Size(); i++)
          {
          ptData[idx] = It.Value()[i];
          ++idx;
          }
        ++It;
        }
      ByteSwapper<int>::SwapRangeFromSystemToBigEndian(ptData,totalSize);
      outputFile.write( reinterpret_cast<char *>( ptData ), totalSize * sizeof(p) );
      delete [] ptData;
      }
    else
      {
      while( It != ItEnd )
        {
        unsigned int numberOfPoints = ( It.Value() ).Size();
        outputFile << numberOfPoints << " ";
        for( unsigned int d = 0; d < numberOfPoints; d++ )
          {
          outputFile << ( It.Value() )[d] << " ";
          }
        outputFile << std::endl;
        ++It;
        }
      }
    outputFile << std::endl;
    outputFile.close();
    }

}

template<class TInputMesh>
void
VtkPolyDataFileWriter<TInputMesh>
::WriteVTKPolygons()
{
  if( this->m_Polygons )
    {

    std::ofstream outputFile( this->m_FileName.c_str(), std::ios::app );

    unsigned int numberOfPolys = this->m_Polygons->Size();
    unsigned int totalSize = 0;

    typename LineSetType::Iterator It
      = this->m_Polygons->Begin();
    typename LineSetType::Iterator ItEnd = this->m_Polygons->End();

    while( It != ItEnd )
      {
      totalSize += ( It.Value() ).Size();
      totalSize++;
      It++;
      }

    outputFile << "POLYGONS " <<
      numberOfPolys << " " << totalSize << std::endl;

    It = this->m_Polygons->Begin();
    while( It != ItEnd )
      {
      unsigned int numberOfPoints = ( It.Value() ).Size();
      outputFile << numberOfPoints << " ";
      for( unsigned int d = 0; d < numberOfPoints; d++ )
        {
        outputFile << ( It.Value() )[d] << " ";
        }
      outputFile << std::endl;
      ++It;
      }
    outputFile << std::endl;
    outputFile.close();
    }
}

template<class TInputMesh>
void
VtkPolyDataFileWriter<TInputMesh>
::WritePointsToAvantsFile()
{
  //
  // Write to output file
  //
  std::ofstream outputFile( this->m_FileName.c_str() );

  outputFile << "0 0 0 0" << std::endl;

  if( this->m_Input->GetNumberOfPoints() > 0 )
    {
    typename InputMeshType::PointsContainerIterator pointIterator
      = this->m_Input->GetPoints()->Begin();
    typename InputMeshType::PointsContainerIterator pointEnd
      = this->m_Input->GetPoints()->End();

    typename InputMeshType::PointDataContainerIterator pointDataIterator
      = this->m_Input->GetPointData()->Begin();

    while( pointIterator != pointEnd )
      {
      PointType point = pointIterator.Value();
      outputFile << point[0] << " " << point[1];
      if( Dimension == 2 )
        {
        outputFile << " 0 ";
        }
      else if( Dimension == 3 )
        {
        outputFile << " " << point[2] << " ";
        }
      outputFile << pointDataIterator.Value() << std::endl;
      pointIterator++;
      pointDataIterator++;
      }
    }

  outputFile << "0 0 0 0" << std::endl;

  outputFile.close();
}

template<class TInputMesh>
void
VtkPolyDataFileWriter<TInputMesh>
::WritePointsToImageFile()
{
  typename VtkPolyDataImageType::Pointer outputImage
    = VtkPolyDataImageType::New();
  outputImage->SetDirection( this->m_ImageDirection );
  outputImage->SetRegions( this->m_ImageSize );
  outputImage->SetOrigin( this->m_ImageOrigin );
  outputImage->SetSpacing( this->m_ImageSpacing );
  outputImage->Allocate();
  outputImage->FillBuffer( NumericTraits<PixelType>::Zero );

  if( this->m_Input->GetNumberOfPoints() > 0 )
    {
    typename InputMeshType::PointsContainerIterator pointIterator
      = this->m_Input->GetPoints()->Begin();
    typename InputMeshType::PointsContainerIterator pointEnd
      = this->m_Input->GetPoints()->End();

    typename InputMeshType::PointDataContainerIterator pointDataIterator
      = this->m_Input->GetPointData()->Begin();

    while( pointIterator != pointEnd )
      {
      PointType point = pointIterator.Value();
      PixelType label = pointDataIterator.Value();

      typename VtkPolyDataImageType::IndexType index;
      typename VtkPolyDataImageType::PointType ipoint;
      ipoint.CastFrom( point );
      if( outputImage->TransformPhysicalPointToIndex( ipoint, index ) )
        {
        outputImage->SetPixel( index, label );
        }
      pointIterator++;
      pointDataIterator++;
      }
    }

  typedef ImageFileWriter<VtkPolyDataImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( this->m_FileName.c_str() );
  writer->SetInput( outputImage );
  writer->Update();
}

template<class TInputMesh>
void
VtkPolyDataFileWriter<TInputMesh>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "FileName: " << this->m_FileName << std::endl;
}

} //end of namespace itk

#endif
