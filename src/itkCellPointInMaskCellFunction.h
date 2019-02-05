#ifndef itkCellPointInMaskCellFunction_h
#define itkCellPointInMaskCellFunction_h

#include "itkCellFunction.h"

namespace itk
{
/** \class CellPointInMaskCellFunction
 * \brief Returns sum of length of all cell edges
 * This CellFunction returns the sum of the lengths for all of the edges
 * in a cell. The input mesh is set via method SetInputMesh().
 *
 * \ingroup CellFunctions
 *
 * \ingroup ITKMeshFunction
 */
template< typename TInputMesh, typename TInputImage >
class ITK_TEMPLATE_EXPORT CellPointInMaskCellFunction:
  public CellFunction< TInputMesh, bool >
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(CellPointInMaskCellFunction);

  /** Standard class type aliases. */
  using Self = CellPointInMaskCellFunction;
  using Superclass = CellFunction< TInputMesh, bool >;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  /** Run-time type information (and related methods). */
  itkTypeMacro(CellPointInMaskCellFunction, ImageFunction);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** InputMeshType type alias support */
  using InputMeshType = typename Superclass::InputMeshType;

  using CellType = typename InputMeshType::CellType;

  using CellAutoPointer = typename CellType::CellAutoPointer;

  using CellPointIterator = typename CellType::PointIdConstIterator;

  //using EdgeAutoPointer = typename CellType::EdgeAutoPointer;

  /** Typedef to describe the type of pixel. */
  using PixelType = typename TInputMesh::PixelType;

  /** Dimension underlying input image. */
  //static constexpr unsigned int Dimension = Superclass::Dimension;

  /** Point type alias support */
  using PointType = typename Superclass::PointType;

  /** Index type alias support */
  using IndexType = typename Superclass::IndexType;

  using PointIndexType = typename TInputMesh::PointIdentifier;

  using InputImageType = TInputImage;

  using ImagePointerType = typename InputImageType::Pointer;

  using ImagePointType = typename InputImageType::PointType;

  using ImageIndexType = typename InputImageType::IndexType;

  /** BinaryThreshold the image at a point position
   *
   * Returns true if the image intensity at the specified point position
   * satisfies the threshold criteria.  The point is assumed to lie within
   * the image buffer.
   *
   * ImageFunction::IsInsideBuffer() can be used to check bounds before
   * calling the method. */

  bool Evaluate(const IndexType & index) const override;

  void SetMask( ImagePointerType mask ) { m_Mask = mask; }


protected:
  CellPointInMaskCellFunction();
  ~CellPointInMaskCellFunction() override = default;
  void PrintSelf(std::ostream & os, Indent indent) const override;

  ImagePointerType m_Mask;

private:

};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkCellPointInMaskCellFunction.hxx"
#endif

#endif
