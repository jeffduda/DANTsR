/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkVtkPolyDataFileReader.h,v $
  Language:  C++
  Date:      $Date: 2009/03/04 23:10:58 $
  Version:   $Revision: 1.17 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkVtkPolyDataFileReader_h
#define __itkVtkPolyDataFileReader_h

#include "itkMesh.h"
#include "itkMeshSource.h"

#include "itkArray.h"
#include "itkByteSwapper.h"
#include "itkImage.h"
#include "itkVectorContainer.h"

#include <vector>
#include <fstream>

#include <RcppDANTsR.h>

namespace itk {

/** \class VtkPolyDataFileReader
 * \brief
 * Reads a file and creates an itkMesh.
 *
 */
template <class TOutputMesh>
class VtkPolyDataFileReader
: public MeshSource<TOutputMesh>
{
public:
  /** Standard "Self" typedef. */
  typedef VtkPolyDataFileReader               Self;
  typedef MeshSource<TOutputMesh>             Superclass;
  typedef SmartPointer<Self>                  Pointer;
  typedef SmartPointer<const Self>            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Extract dimension from the output mesh. */
  itkStaticConstMacro( Dimension, unsigned int,
                       TOutputMesh::PointType::Dimension );


  /** Run-time type information (and related methods). */
  itkTypeMacro( VtkPolyDataFileReader, MeshSource );

  /** Hold on to the type information specified by the template parameters. */
  typedef TOutputMesh                             OutputMeshType;
  typedef typename OutputMeshType::MeshTraits     MeshTraits;
  typedef typename OutputMeshType::Superclass     PointSetType;
  typedef typename OutputMeshType::PointType      PointType;
  typedef typename MeshTraits::PixelType          PixelType;
  typedef Array<float>                            MultiComponentScalarType;

  typedef typename Rcpp::NumericMatrix            MatrixType;
  typedef VectorContainer<long, MatrixType>       MatrixSetType;

  typedef Array<unsigned long>                    LineType;
  typedef VectorContainer<long,
    MultiComponentScalarType>                     MultiComponentScalarSetType;

  typedef VectorContainer<long, std::string>      MultiComponentScalarSetNamesType;
  typedef VectorContainer<long,
    SmartPointer<MultiComponentScalarSetType> >   MultiComponentScalarMultiSetType;

  typedef VectorContainer<long, LineType>         LineSetType;

  typedef Image<PixelType,
    itkGetStaticConstMacro( Dimension )>          VtkPolyDataImageType;

  typedef std::vector<PixelType>                  LabelSetType;

  /** Set/Get the name of the file to be read. */
  itkSetStringMacro( FileName );
  itkGetStringMacro( FileName );

  itkSetMacro( ExtractBoundaryPoints, bool );
  itkGetMacro( ExtractBoundaryPoints, bool );
  itkBooleanMacro( ExtractBoundaryPoints );

  /**
   * Percentage of points selected randomnly
   */
  itkSetClampMacro( RandomPercentage, double, 0.0, 1.0 );
  itkGetConstMacro( RandomPercentage, double );

  LabelSetType* GetLabelSet() { return &this->m_LabelSet; }
  unsigned int GetNumberOfLabels() { return this->m_LabelSet.size(); }

  MultiComponentScalarMultiSetType* GetMultiComponentScalarSets()
    { return this->m_MultiComponentScalarSets.GetPointer(); }

  MultiComponentScalarMultiSetType* GetCellMultiComponentScalarSets()
    { return this->m_CellMultiComponentScalarSets.GetPointer(); }

  MultiComponentScalarSetNamesType* GetMultiComponentScalarSetNames()
    { return this->m_MultiComponentScalarSetNames.GetPointer(); }

  MultiComponentScalarSetType* GetMultiComponentScalarSetByName( const char * );

  MultiComponentScalarSetType* GetMultiComponentScalarSet( long i )
  { return this->m_MultiComponentScalarSets->GetElement( i ).GetPointer(); }

  MultiComponentScalarSetType* GetCellMultiComponentScalarSet( long i )
  { return this->m_CellMultiComponentScalarSets->GetElement( i ).GetPointer(); }

  //MultiComponentScalarSetType* GetMultiComponentScalarsByName()
  //  { return this->m_MultiComponentScalars.GetPointer(); }

  LineSetType* GetLines()
    { return this->m_Lines.GetPointer(); }

  LineSetType* GetPolygons()
    { return this->m_Polygons.GetPointer(); }

  LineSetType* GetVertices()
    { return this->m_Vertices.GetPointer(); }

  LineSetType* GetTriangularStrips()
    { return this->m_Strips.GetPointer(); }

  template <typename ValueType>
  ValueType * ReadVTKBinaryData( unsigned long nValues ) {
    ValueType p;
    ValueType * data = new ValueType [ nValues ];
    this->m_InputFile.read( reinterpret_cast< char * >( data ), nValues * sizeof(p) );
    ByteSwapper<ValueType>::SwapRangeFromSystemToBigEndian(data,nValues);
    return data;
  }

  template <typename ValueType>
  bool ReadVTKBinaryMatrix(  MatrixType normMatrix ) {
    ValueType * vec = ReadVTKBinaryData<ValueType>(normMatrix.rows()*normMatrix.cols());
    for ( unsigned long i=0; i<normMatrix.rows()*normMatrix.cols(); i++ ) {
      normMatrix[i] = vec[i];
      }
    delete [] vec;
    return true;
  }

protected:
  VtkPolyDataFileReader();
  ~VtkPolyDataFileReader() {}
  void PrintSelf( std::ostream& os, Indent indent ) const override;

  /** Reads the file */
  void GenerateData() override;

  bool                                            m_ExtractBoundaryPoints;

  std::string                                     m_FileName;
  double                                          m_RandomPercentage;
  LabelSetType                                    m_LabelSet;

  typename MultiComponentScalarMultiSetType::Pointer   m_MultiComponentScalarSets;
  typename MultiComponentScalarMultiSetType::Pointer   m_CellMultiComponentScalarSets;
  typename LineSetType::Pointer                        m_Lines;
  typename LineSetType::Pointer                        m_Polygons;
  typename LineSetType::Pointer                        m_Vertices;
  typename LineSetType::Pointer                        m_Strips;

  typedef enum { Scalars, ColorScalars, LookupTable, Vectors, Normals, TextureCoordinates, Tensors, Field } VTKDataSetType;

  typename MatrixSetType::Pointer                 m_PointNormals;
  typename MatrixSetType::Pointer                 m_CellNormals;

  typename MatrixSetType::Pointer                 m_PointVectors;
  typename MatrixSetType::Pointer                 m_CellVectors;

  typename MatrixSetType::Pointer                 m_PointTextureCoordinates;
  typename MatrixSetType::Pointer                 m_CellTextureCoordinates;

  bool                                            m_BinaryData;
  std::ifstream                                   m_InputFile;

  typename MultiComponentScalarSetNamesType::Pointer m_MultiComponentScalarSetNames;
  typename MultiComponentScalarSetNamesType::Pointer m_CellMultiComponentScalarSetNames;

private:
  VtkPolyDataFileReader( const Self& ); // purposely not implemented
  void operator=( const Self& ); // purposely not implemented

  void ReadVTKFile();
  void ReadPointsFromVTKFile();
  void ReadScalarsFromVTKFile();
  void ReadLinesFromVTKFile();

  bool ReadVTKPolyDataFile();
  bool ReadVTKPolyData();
  bool ReadVTKDataAttributes();
  bool ReadVTKPoints(unsigned long nPoints, std::string dataType);
  bool ReadVTKLines(unsigned long nLines, unsigned long nValues);
  bool ReadVTKPolygons(unsigned long nLines, unsigned long nValues);
  bool ReadVTKVertices(unsigned long nLines, unsigned long nValues);
  bool ReadVTKStrips(unsigned long nLines, unsigned long nValues);
  bool ReadVTKPointData(unsigned long nPoints);
  bool ReadVTKCellData(unsigned long nCells);

  bool ReadVTKDataSet(std::string dataName, std::string dataType, VTKDataSetType type, unsigned long nRows, unsigned long nCols, bool isPointData);

  bool ReadVTKScalars(std::string dataName, std::string dataType, unsigned long nPoints, unsigned long nComponents, bool isPointData);
  bool ReadVTKVectors(std::string dataName, std::string dataType, unsigned long nPoints, bool isPointData, bool isNormal);
  bool ReadVTKTextureCoordinates(std::string dataName, std::string dataType, unsigned long nRows, unsigned long dim, bool isPointData);


};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkVtkPolyDataFileReader.hxx"
#endif

#endif //_itkVtkPolyDataFileReader_h
