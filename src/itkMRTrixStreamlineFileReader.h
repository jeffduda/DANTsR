#ifndef __itkMRTrixStreamlineFileReader_h
#define __itkMRTrixStreamlineFileReader_h

#include "itkMesh.h"
#include "itkMeshSource.h"
#include "itkArray.h"
#include "itkImage.h"
#include "itkVectorContainer.h"
#include "itkPolyLineCell.h"

namespace itk {

/** \class MRTrixStreamlineFileReader
 * \brief
 * Writes an itkMesh to a file in various txt file formats.
 *
 */
template <class TOutputMesh>
class  MRTrixStreamlineFileReader : public MeshSource<TOutputMesh>
{
public:
  /** Standard "Self" typedef. */
  typedef MRTrixStreamlineFileReader          Self;
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
  itkTypeMacro( MRTrixStreamlineFileReader, Object );

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

  typedef typename std::map< std::string, std::string> HeaderType;

  typedef typename OutputMeshType::CellType       CellType;
  typedef typename CellType::CellAutoPointer      CellAutoPointer;
  typedef typename itk::PolyLineCell< CellType >  PolyLineCellType;

  /** Set/Get the name of the file where data are written. */
  itkSetStringMacro( FileName );
  itkGetStringMacro( FileName );

protected:
  MRTrixStreamlineFileReader();
  virtual ~MRTrixStreamlineFileReader();

  virtual void GenerateData() ITK_OVERRIDE;

  std::string                         m_FileName;

  HeaderType m_Header;

  std::ifstream m_InputFile;

  unsigned long m_SwapFlag;

  void PrintSelf(std::ostream& os, Indent indent) const override;

  void ReadTckFile();

  void ReadTckHeader();

  void ReadTckTracts(long dataSize, unsigned long count);

  void PrintTckHeader();

  std::string GetHeaderValue( std::string key, bool required=false );

private:
  MRTrixStreamlineFileReader(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

} // namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMRTrixStreamlineFileReader.hxx"
#endif

#endif
