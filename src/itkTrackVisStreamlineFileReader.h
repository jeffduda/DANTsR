#ifndef __itkTrackVisStreamlineFileReader_h
#define __itkTrackVisStreamlineFileReader_h

#include "itkMesh.h"
#include "itkMeshSource.h"
#include "itkArray.h"
#include "itkImage.h"
#include "itkVectorContainer.h"
#include "itkPolyLineCell.h"
#include "trackVisHeader.h"

namespace itk {

/** \class TrackVisStreamlineFileReader
 * \brief
 * Writes an itkMesh to a file in various txt file formats.
 *
 */
template <class TOutputMesh>
class  TrackVisStreamlineFileReader : public MeshSource<TOutputMesh>
{
public:
  /** Standard "Self" typedef. */
  typedef TrackVisStreamlineFileReader        Self;
  typedef Object                              Superclass;
  typedef SmartPointer<Self>                  Pointer;
  typedef SmartPointer<const Self>            ConstPointer;

  /** Method for creation through the object factory */
  itkNewMacro( Self );

  /** Write the Input mesh to the Output file.
   * Use either Update() or Write(). */
  void Update( void ) ITK_OVERRIDE;
  void Read( void );

  /** Extract dimension from the output mesh. */
  itkStaticConstMacro( Dimension, unsigned int,
                       TOutputMesh::PointType::Dimension );


  /** Run-time type information (and related methods). */
  itkTypeMacro( TrackVisStreamlineFileReader, Object );

  /** Hold on to the type information specified by the template parameters. */
  typedef TOutputMesh                             OutputMeshType;
  typedef typename TOutputMesh::Pointer           OutputMeshPointer;
  typedef typename OutputMeshType::MeshTraits     MeshTraits;
  typedef typename OutputMeshType::Superclass     PointSetType;
  typedef typename OutputMeshType::PointType      PointType;
  typedef typename MeshTraits::PixelType         PixelType;
  typedef Array<PixelType>                       MultiComponentScalarType;
  typedef Array<unsigned long>                   LineType;
  typedef VectorContainer<long,
    MultiComponentScalarType>                    MultiComponentScalarSetType;
  typedef VectorContainer<long,
    SmartPointer<MultiComponentScalarSetType> >   MultiComponentScalarMultiSetType;
  typedef VectorContainer<long, std::string >     MultiComponentScalarSetNamesType;

  using ImageType = typename itk::Image<PixelType,3>;
  using ImagePointerType = typename ImageType::Pointer;

  typedef typename OutputMeshType::CellType       CellType;
  typedef typename CellType::CellAutoPointer      CellAutoPointer;
  typedef typename itk::PolyLineCell< CellType >  PolyLineCellType;

  typedef typename
    ImageType::SizeType           ImageSizeType;
  typedef typename
    ImageType::PointType          ImageOriginType;
  typedef typename
    ImageType::SpacingType        ImageSpacingType;
  typedef typename
    ImageType::DirectionType      ImageDirectionType;

  itkGetMacro( ReferenceImage, ImagePointerType );

  /** Set/Get the name of the file where data are written. */
  itkSetStringMacro( FileName );
  itkGetStringMacro( FileName );

  itkSetMacro( MultiComponentScalarSets,
    typename MultiComponentScalarMultiSetType::Pointer );

  itkSetMacro( MultiComponentScalarSetNames,
    typename MultiComponentScalarSetNamesType::Pointer );

  //itkSetMacro( ReferenceImage, ImagePointerType );
  //itkGetConstObjectMacro( ReferenceImage, ImageType );

  /** Specify image attributes if output is an image. */
  //itkSetMacro( ImageSize, ImageSizeType );
  //itkGetConstMacro( ImageSize, ImageSizeType );

  //itkSetMacro( ImageOrigin, ImageOriginType );
  //itkGetConstMacro( ImageOrigin, ImageOriginType );

  //itkSetMacro( ImageSpacing, ImageSpacingType );
  //itkGetConstMacro( ImageSpacing, ImageSpacingType );

  //itkSetMacro( ImageDirection, ImageDirectionType );
  //itkGetConstMacro( ImageDirection, ImageDirectionType );



  ImagePointerType m_ReferenceImage;

protected:
  TrackVisStreamlineFileReader();
  virtual ~TrackVisStreamlineFileReader();

  virtual void GenerateData() ITK_OVERRIDE;

  std::string                         m_FileName;

  typename MultiComponentScalarMultiSetType::Pointer   m_MultiComponentScalarSets;

  short int m_NScalars;
  short int m_NProperties;

  /**
   * If output is an image type, the attributes must be specified.
   */
  //ImageSizeType                       m_ImageSize;
  //ImageSpacingType                    m_ImageSpacing;
  //ImageOriginType                     m_ImageOrigin;
  //ImageDirectionType                  m_ImageDirection;

  typename MultiComponentScalarSetNamesType::Pointer m_MultiComponentScalarSetNames;

  void PrintSelf(std::ostream& os, Indent indent) const override;

private:
  TrackVisStreamlineFileReader(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  std::ifstream m_InputFile;

  void ReadTrkFile();

  TRACKVIS_HEADER_V2 ReadTrkHeader();

  void ReadTrkTract();

  void PrintTrkHeader( TRACKVIS_HEADER_V2 hdr );

};

} // namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTrackVisStreamlineFileReader.hxx"
#endif

#endif
