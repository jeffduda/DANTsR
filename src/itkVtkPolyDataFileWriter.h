/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkVtkPolyDataFileWriter.h,v $
  Language:  C++
  Date:      $Date: 2009/03/04 23:10:58 $
  Version:   $Revision: 1.17 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkVtkPolyDataFileWriter_h
#define __itkVtkPolyDataFileWriter_h

#include "itkMesh.h"

#include "itkArray.h"
#include "itkImage.h"
#include "itkVectorContainer.h"

namespace itk {

/** \class VtkPolyDataFileWriter
 * \brief
 * Writes an itkMesh to a file in various txt file formats.
 *
 */
template <class TInputMesh>
class VtkPolyDataFileWriter : public Object
{
public:
  /** Standard "Self" typedef. */
  typedef VtkPolyDataFileWriter               Self;
  typedef Object                              Superclass;
  typedef SmartPointer<Self>                  Pointer;
  typedef SmartPointer<const Self>            ConstPointer;

  /** Method for creation through the object factory */
  itkNewMacro( Self );

  /** Write the Input mesh to the Output file.
   * Use either Update() or Write(). */
  void Update( void );
  void Write( void );

  /** Extract dimension from the output mesh. */
  itkStaticConstMacro( Dimension, unsigned int,
                       TInputMesh::PointType::Dimension );


  /** Run-time type information (and related methods). */
  itkTypeMacro( VtkPolyDataFileWriter, Object );

  /** Hold on to the type information specified by the template parameters. */
  typedef TInputMesh                             InputMeshType;
  typedef typename TInputMesh::Pointer           InputMeshPointer;
  typedef typename InputMeshType::MeshTraits     MeshTraits;
  typedef typename InputMeshType::Superclass     PointSetType;
  typedef typename InputMeshType::PointType      PointType;
  typedef typename MeshTraits::PixelType         PixelType;

  typedef typename Rcpp::NumericMatrix            MatrixType;
  typedef VectorContainer<long, MatrixType>       MatrixSetType;

  typedef Array<PixelType>                       MultiComponentScalarType;
  typedef Array<unsigned long>                   LineType;
  typedef VectorContainer<long,
    MultiComponentScalarType>                    MultiComponentScalarSetType;
  typedef VectorContainer<long,
    SmartPointer<MultiComponentScalarSetType> >   MultiComponentScalarMultiSetType;
  typedef VectorContainer<long, std::string >     MultiComponentScalarSetNamesType;

  typedef VectorContainer<long, std::string>      DataNameSetType;


  typedef typename InputMeshType::CellType      CellType;
  typedef typename CellType::CellAutoPointer    CellAutoPointer;

  typedef VectorContainer<long, LineType>        LineSetType;
  typedef Image<PixelType,
    itkGetStaticConstMacro( Dimension )>         VtkPolyDataImageType;
  typedef typename
    VtkPolyDataImageType::SizeType           ImageSizeType;
  typedef typename
    VtkPolyDataImageType::PointType          ImageOriginType;
  typedef typename
    VtkPolyDataImageType::SpacingType        ImageSpacingType;
  typedef typename
    VtkPolyDataImageType::DirectionType      ImageDirectionType;


  /** Set the Input */
  void SetInput( InputMeshType * input );

  void SetCellsAsPolygons( bool );
  void SetCellsAsLines( bool );

  /** Set/Get the name of the file where data are written. */
  itkSetStringMacro( FileName );
  itkGetStringMacro( FileName );

  itkSetMacro( WriteBinary, bool );
  itkGetMacro( WriteBinary, bool );

  /** Specify other attributes */
  itkSetMacro( Lines, typename LineSetType::Pointer );

  itkSetMacro( Polygons, typename LineSetType::Pointer );

  itkSetMacro( MultiComponentScalarSets,
    typename MultiComponentScalarMultiSetType::Pointer );

  itkSetMacro( MultiComponentScalarSetNames,
    typename MultiComponentScalarSetNamesType::Pointer );

  /** Specify image attributes if output is an image. */
  itkSetMacro( ImageSize, ImageSizeType );
  itkGetConstMacro( ImageSize, ImageSizeType );

  itkSetMacro( ImageOrigin, ImageOriginType );
  itkGetConstMacro( ImageOrigin, ImageOriginType );

  itkSetMacro( ImageSpacing, ImageSpacingType );
  itkGetConstMacro( ImageSpacing, ImageSpacingType );

  itkSetMacro( ImageDirection, ImageDirectionType );
  itkGetConstMacro( ImageDirection, ImageDirectionType );

  itkGetObjectMacro( PointNormals, MatrixSetType);
  itkGetObjectMacro( PointNormalsNames, DataNameSetType);

  MultiComponentScalarMultiSetType* GetMultiComponentScalarSets()
    { return this->m_MultiComponentScalarSets.GetPointer(); }

protected:
  VtkPolyDataFileWriter();
  virtual ~VtkPolyDataFileWriter();

  virtual void GenerateData();

  std::string                         m_FileName;
  InputMeshPointer                    m_Input;

  typename MultiComponentScalarMultiSetType::Pointer   m_MultiComponentScalarSets;
  typename LineSetType::Pointer                        m_Lines;
  typename LineSetType::Pointer                        m_Polygons;

  typename MatrixSetType::Pointer                 m_PointNormals;
  typename MatrixSetType::Pointer                 m_CellNormals;

  typename DataNameSetType::Pointer               m_PointNormalsNames;
  typename DataNameSetType::Pointer               m_CellNormalsNames;

  /**
   * If output is an image type, the attributes must be specified.
   */
  ImageSizeType                       m_ImageSize;
  ImageSpacingType                    m_ImageSpacing;
  ImageOriginType                     m_ImageOrigin;
  ImageDirectionType                  m_ImageDirection;

  typename MultiComponentScalarSetNamesType::Pointer m_MultiComponentScalarSetNames;

  void PrintSelf(std::ostream& os, Indent indent) const override;

  template <typename ValueType>
  void WriteVTKBinaryMatrix( MatrixType matrix ) {
    std::cout << "WriteVTKBinaryMatrix()" << std::endl;
    //ValueType * vec = ReadVTKBinaryData<ValueType>(normMatrix.rows()*normMatrix.cols());
    //for ( unsigned long i=0; i<normMatrix.rows()*normMatrix.cols(); i++ ) {
    //  normMatrix[i] = vec[i];
    //  }
    //delete [] vec;

/*
    float p;
    unsigned long totalSize = matrix.nrows()*matrix.ncols()
    float * data = new float [ totalSize ];

    unsigned long idx = 0;
    for ( unsigned long i=0; i<matrix.nrows(); i++)  {
      for (unsigned long j=0; j<matrix.ncols(); j++) {
        data[idx] = static_cast<float>( matrix(i,j) );
        ++idx;
      }
    }

    ByteSwapper<float>::SwapRangeFromSystemToBigEndian(data,totalSize);
    this->m_OutputFile.write( reinterpret_cast<char *>( data ), totalSize * sizeof(p) );
    delete [] data;
    */
  }

  void WriteVTKASCIIMatrix( MatrixType matrix ) {
    std::cout << "WriteVTKASCIIMatrix()" << std::endl;
    //double value=0;
    this->m_OutputFile << std::setprecision(10) << std::fixed;

    for (unsigned long i=0; i<matrix.rows(); i++) {
      for (unsigned long j=0; j<matrix.cols(); j++ ) {
        this->m_OutputFile << matrix(i,j) << " ";
      }
      this->m_OutputFile << std::endl;
    }
    //  this->m_InputFile >> value;
    //  matrix[i] = value;
    //}
    //std::string line;
    //std::getline( this->m_InputFile, line );
  }

  std::ofstream                                   m_OutputFile;

private:
  VtkPolyDataFileWriter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  bool m_WriteBinary;

  bool m_CellsAsPolygons;
  bool m_CellsAsLines;

  void WriteVTKFile();
  void WritePoints();
  void WritePointData();
  void WriteCellData();
  void WriteLines();
  void WritePolygons();
  void WriteCellsAs(std::string);
  void WriteVTKDataMatrix( MatrixType, std::string);


};

} // namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkVtkPolyDataFileWriter.hxx"
#endif

#endif
