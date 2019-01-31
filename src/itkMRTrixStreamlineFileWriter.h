#ifndef __itkMRTrixStreamlineFileWriter_h
#define __itkMRTrixStreamlineFileWriter_h

#include "itkMesh.h"
#include "itkArray.h"
#include "itkImage.h"
#include "itkVectorContainer.h"
#include "MRTrixHeader.h"

namespace itk {

/** \class MRTrixStreamlineFileWriter
 * \brief
 * Writes an itkMesh to a file in various txt file formats.
 *
 */
template <class TInputMesh, class TInputImage>
class  MRTrixStreamlineFileWriter : public Object
{
public:
  /** Standard "Self" typedef. */
  typedef MRTrixStreamlineFileWriter        Self;
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
  itkTypeMacro( MRTrixStreamlineFileWriter, Object );

  /** Hold on to the type information specified by the template parameters. */
  typedef TInputMesh                             InputMeshType;
  typedef typename TInputMesh::Pointer           InputMeshPointer;
  typedef typename InputMeshType::MeshTraits     MeshTraits;
  typedef typename InputMeshType::Superclass     PointSetType;
  typedef typename InputMeshType::PointType      PointType;
  typedef typename MeshTraits::PixelType         PixelType;
  typedef Array<PixelType>                       MultiComponentScalarType;
  typedef Array<unsigned long>                   LineType;
  typedef VectorContainer<long,
    MultiComponentScalarType>                    MultiComponentScalarSetType;
  typedef VectorContainer<long,
    SmartPointer<MultiComponentScalarSetType> >   MultiComponentScalarMultiSetType;
  typedef VectorContainer<long, std::string >     MultiComponentScalarSetNamesType;

  using ImageType = TInputImage;
  using ImagePointerType = typename ImageType::Pointer;


  typedef typename InputMeshType::CellType      CellType;
  typedef typename CellType::CellAutoPointer    CellAutoPointer;

  typedef VectorContainer<long, LineType>        LineSetType;

  typedef typename
    ImageType::SizeType           ImageSizeType;
  typedef typename
    ImageType::PointType          ImageOriginType;
  typedef typename
    ImageType::SpacingType        ImageSpacingType;
  typedef typename
    ImageType::DirectionType      ImageDirectionType;


  /** Set the Input */
  void SetInput( InputMeshType * input );

  void SetReferenceImage( ImageType * image );

  /** Set/Get the name of the file where data are written. */
  itkSetStringMacro( FileName );
  itkGetStringMacro( FileName );

  /** Specify other attributes */
  itkSetMacro( MultiComponentScalarSets,
    typename MultiComponentScalarMultiSetType::Pointer );

  itkSetMacro( MultiComponentScalarSetNames,
    typename MultiComponentScalarSetNamesType::Pointer );

  itkSetMacro( ReferenceImage, ImagePointerType );
  //itkGetConstObjectMacro( ReferenceImage, ImageType );

  /** Specify image attributes if output is an image. */
  itkSetMacro( ImageSize, ImageSizeType );
  itkGetConstMacro( ImageSize, ImageSizeType );

  itkSetMacro( ImageOrigin, ImageOriginType );
  itkGetConstMacro( ImageOrigin, ImageOriginType );

  itkSetMacro( ImageSpacing, ImageSpacingType );
  itkGetConstMacro( ImageSpacing, ImageSpacingType );

  itkSetMacro( ImageDirection, ImageDirectionType );
  itkGetConstMacro( ImageDirection, ImageDirectionType );

protected:
  MRTrixStreamlineFileWriter();
  virtual ~MRTrixStreamlineFileWriter();

  virtual void GenerateData();

  std::string                         m_FileName;
  InputMeshPointer                    m_Input;
  ImagePointerType                    m_ReferenceImage;

  typename MultiComponentScalarMultiSetType::Pointer   m_MultiComponentScalarSets;

  short int m_NScalars;
  short int m_NProperties;

    /**
   * If output is an image type, the attributes must be specified.
   */
  ImageSizeType                       m_ImageSize;
  ImageSpacingType                    m_ImageSpacing;
  ImageOriginType                     m_ImageOrigin;
  ImageDirectionType                  m_ImageDirection;

  typename MultiComponentScalarSetNamesType::Pointer m_MultiComponentScalarSetNames;

  void PrintSelf(std::ostream& os, Indent indent) const override;

private:
  MRTrixStreamlineFileWriter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  std::ofstream m_OutputFile;

  void WriteTckFile();
  void WriteTckHeader();
  void WriteTckTract(unsigned int i);

};

} // namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMRTrixStreamlineFileWriter.hxx"
#endif

#endif
