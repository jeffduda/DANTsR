/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkCaminoStreamlineFileReader.h,v $
  Language:  C++
  Date:      $Date: 2009/03/04 23:10:58 $
  Version:   $Revision: 1.17 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkCaminoStreamlineFileReader_h
#define __itkCaminoStreamlineFileReader_h

#include "itkMesh.h"
#include "itkMeshSource.h"
#include "itkLineCell.h"
#include "itkPolygonCell.h"

#include "itkArray.h"
#include "itkByteSwapper.h"
#include "itkImage.h"
#include "itkVectorContainer.h"

#include <vector>
#include <fstream>

#include <RcppDANTsR.h>

namespace itk {

/** \class CaminoStreamlineFileReader
 * \brief
 * Reads a file and creates an itkMesh.
 *
 */
template <class TOutputMesh>
class CaminoStreamlineFileReader
: public MeshSource<TOutputMesh>
{
public:
  /** Standard "Self" typedef. */
  typedef CaminoStreamlineFileReader               Self;
  typedef MeshSource<TOutputMesh>             Superclass;
  typedef SmartPointer<Self>                  Pointer;
  typedef SmartPointer<const Self>            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Extract dimension from the output mesh. */
  itkStaticConstMacro( Dimension, unsigned int,
                       TOutputMesh::PointType::Dimension );


  /** Run-time type information (and related methods). */
  itkTypeMacro( CamionStreamlineFileReader, MeshSource );

  /** Hold on to the type information specified by the template parameters. */
  typedef TOutputMesh                             OutputMeshType;
  typedef typename OutputMeshType::CellType       CellType;
  typedef typename CellType::CellAutoPointer      CellAutoPointer;

  typedef LineCell<CellType>                      LineCellType;
  typedef PolygonCell<CellType>                   PolygonCellType;

  typedef typename OutputMeshType::MeshTraits     MeshTraits;
  typedef typename OutputMeshType::Superclass     PointSetType;
  typedef typename OutputMeshType::PointType      PointType;
  typedef typename MeshTraits::PixelType          PixelType;
  typedef Array<float>                            MultiComponentScalarType;

  typedef typename Rcpp::NumericMatrix            MatrixType;
  typedef VectorContainer<long, MatrixType>       MatrixSetType;

  typedef VectorContainer<long, std::string>      DataNameSetType;
  typedef VectorContainer<long, unsigned long>    SeedSetType;

  typedef Array<unsigned long>                    LineType;
  typedef VectorContainer<long,
    MultiComponentScalarType>                     MultiComponentScalarSetType;

  typedef VectorContainer<long, std::string>      MultiComponentScalarSetNamesType;
  typedef VectorContainer<long,
    SmartPointer<MultiComponentScalarSetType> >   MultiComponentScalarMultiSetType;

  typedef VectorContainer<long, LineType>         LineSetType;

  typedef Image<PixelType,
    itkGetStaticConstMacro( Dimension )>          CaminoStreamlineImageType;

  typedef std::vector<PixelType>                  LabelSetType;

  /** Set/Get the name of the file to be read. */
  itkSetStringMacro( FileName );
  itkGetStringMacro( FileName );

  itkSetMacro( ExtractBoundaryPoints, bool );
  itkGetMacro( ExtractBoundaryPoints, bool );
  itkBooleanMacro( ExtractBoundaryPoints );

  itkGetObjectMacro( Lines, LineSetType);

  itkGetObjectMacro( Seeds, SeedSetType);

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

  bool ReadVTKASCIIMatrix( MatrixType matrix ) {
    double value=0;
    for (unsigned long i=0; i<matrix.rows()*matrix.cols(); i++) {
      this->m_InputFile >> value;
      matrix[i] = value;
    }
    std::string line;
    std::getline( this->m_InputFile, line );
    return true;
  }

protected:
  CaminoStreamlineFileReader();
  ~CaminoStreamlineFileReader() {}
  void PrintSelf( std::ostream& os, Indent indent ) const override;

  /** Reads the file */
  void GenerateData() override;

  bool                                            m_ExtractBoundaryPoints;

  std::string                                     m_FileName;

  typename MultiComponentScalarMultiSetType::Pointer   m_MultiComponentScalarSets;
  typename MultiComponentScalarMultiSetType::Pointer   m_CellMultiComponentScalarSets;
  typename SeedSetType::Pointer                        m_Seeds;
  typename LineSetType::Pointer                        m_Lines;
  typename LineSetType::Pointer                        m_Polygons;
  typename LineSetType::Pointer                        m_Vertices;
  typename LineSetType::Pointer                        m_Strips;

  typedef enum { Scalars, ColorScalars, LookupTable, Vectors, Normals, TextureCoordinates, Tensors, Field } VTKDataSetType;

  typename MatrixSetType::Pointer                 m_PointNormals;
  typename MatrixSetType::Pointer                 m_CellNormals;

  typename DataNameSetType::Pointer               m_PointNormalsNames;
  typename DataNameSetType::Pointer               m_CellNormalsNames;

  typename MatrixSetType::Pointer                 m_PointVectors;
  typename MatrixSetType::Pointer                 m_CellVectors;

  typename DataNameSetType::Pointer               m_PointVectorsNames;
  typename DataNameSetType::Pointer               m_CellVectorsNames;

  typename MatrixSetType::Pointer                 m_PointTextureCoordinates;
  typename MatrixSetType::Pointer                 m_CellTextureCoordinates;

  typename DataNameSetType::Pointer               m_PointTextureCoordinatesNames;
  typename DataNameSetType::Pointer               m_CellTextureCoordinatesNames;

  typename MatrixSetType::Pointer                 m_PointTensors;
  typename MatrixSetType::Pointer                 m_CellTensors;

  typename DataNameSetType::Pointer               m_PointTensorsNames;
  typename DataNameSetType::Pointer               m_CellTensorsNames;

  bool                                            m_BinaryData;
  std::ifstream                                   m_InputFile;

  typename MultiComponentScalarSetNamesType::Pointer m_MultiComponentScalarSetNames;
  typename MultiComponentScalarSetNamesType::Pointer m_CellMultiComponentScalarSetNames;

private:
  CaminoStreamlineFileReader( const Self& ); // purposely not implemented
  void operator=( const Self& ); // purposely not implemented

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkCaminoStreamlineFileReader.hxx"
#endif

#endif //_itkCaminoStreamlineFileReader_h
