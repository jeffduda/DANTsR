#ifndef __itkCaminoStreamlineFileWriter_h
#define __itkCaminoStreamlineFileWriter_h

#include "itkMesh.h"
#include "itkArray.h"
#include "itkImage.h"
#include "itkVectorContainer.h"

namespace itk {

/** \class CaminoStreamlineFileWriter
 * \brief
 * Writes an itkMesh to a file in various txt file formats.
 *
 */
template <class TInputMesh>
class  CaminoStreamlineFileWriter : public Object
{
public:
  /** Standard "Self" typedef. */
  typedef CaminoStreamlineFileWriter        Self;
  typedef Object                              Superclass;
  typedef SmartPointer<Self>                  Pointer;
  typedef SmartPointer<const Self>            ConstPointer;

  /** Method for creation through the object factory */
  itkNewMacro( Self );

  /** Write the Input mesh to the Output file.
   * Use either Update() or Write(). */
  void Update( void );
  void Write( void );

  void SetSeed( long i, unsigned long seed );

  /** Extract dimension from the output mesh. */
  itkStaticConstMacro( Dimension, unsigned int,
                       TInputMesh::PointType::Dimension );


  /** Run-time type information (and related methods). */
  itkTypeMacro( CaminoStreamlineFileWriter, Object );

  /** Hold on to the type information specified by the template parameters. */
  typedef TInputMesh                             InputMeshType;
  typedef typename TInputMesh::Pointer           InputMeshPointer;
  typedef typename InputMeshType::MeshTraits     MeshTraits;
  typedef typename InputMeshType::Superclass     PointSetType;
  typedef typename InputMeshType::PointType      PointType;
  typedef typename MeshTraits::PixelType         PixelType;
  typedef typename InputMeshType::CellType       CellType;
  typedef typename CellType::CellAutoPointer     CellAutoPointer;
  typedef VectorContainer<long, unsigned long>   SeedSetType;

  /** Set the Input */
  void SetInput( InputMeshType * input );

  /** Set/Get the name of the file where data are written. */
  itkSetStringMacro( FileName );
  itkGetStringMacro( FileName );

protected:
  CaminoStreamlineFileWriter();
  virtual ~CaminoStreamlineFileWriter();

  virtual void GenerateData();

  std::string                         m_FileName;
  InputMeshPointer                    m_Input;
  typename SeedSetType::Pointer       m_Seeds;

  void PrintSelf(std::ostream& os, Indent indent) const override;

private:
  CaminoStreamlineFileWriter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  std::ofstream m_OutputFile;

  bool m_WriteFloat;

  void WriteCaminoFile();
  template <typename PrecisionType> void WriteCaminoTract( unsigned int i );

};

} // namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkCaminoStreamlineFileWriter.hxx"
#endif

#endif
