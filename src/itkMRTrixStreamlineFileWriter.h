#ifndef __itkMRTrixStreamlineFileWriter_h
#define __itkMRTrixStreamlineFileWriter_h

#include "itkMesh.h"
#include "itkArray.h"
#include "itkImage.h"
#include "itkVectorContainer.h"

namespace itk {

/** \class MRTrixStreamlineFileWriter
 * \brief
 * Writes an itkMesh to a file in various txt file formats.
 *
 */
template <class TInputMesh>
class  MRTrixStreamlineFileWriter : public Object
{
public:
  /** Standard "Self" typedef. */
  typedef MRTrixStreamlineFileWriter          Self;
  typedef Object                              Superclass;
  typedef SmartPointer<Self>                  Pointer;
  typedef SmartPointer<const Self>            ConstPointer;

  /** Method for creation through the object factory */
  itkNewMacro( Self );

  /** Write the Input mesh to the Output file.
   * Use either Update() or Write(). */
  void Update( void );
  void Write( void );

  // FIXME - make enum
  void SetWriteFloat32()   { m_SwapFlag = 0; }
  void SetWriteFloat32BE() { m_SwapFlag = 1; }
  void SetWriteFloat32LE() { m_SwapFlag = 2; }

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
  typedef typename InputMeshType::CellType       CellType;
  typedef typename CellType::CellAutoPointer     CellAutoPointer;


  /** Set/Get the name of the file where data are written. */
  itkSetStringMacro( FileName );
  itkGetStringMacro( FileName );

  void SetInput( InputMeshType * input );

protected:
  MRTrixStreamlineFileWriter();
  virtual ~MRTrixStreamlineFileWriter();

  virtual void GenerateData();

  std::string                         m_FileName;
  InputMeshPointer                    m_Input;
  unsigned int                        m_SwapFlag;

  std::ofstream m_OutputFile;

  void WriteTckFile();
  void WriteTckHeader();
  void WriteTckTract(unsigned int i);

  void PrintSelf(std::ostream& os, Indent indent) const override;

private:
  MRTrixStreamlineFileWriter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented



};

} // namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMRTrixStreamlineFileWriter.hxx"
#endif

#endif
