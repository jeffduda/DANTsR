#ifndef itkCellImageValueSummaryCellFunction_h
#define itkCellImageValueSummaryCellFunction_h

#include "itkCellFunction.h"

namespace itk
{
/** \class CellImageValueSummaryCellFunction
 * \brief Returns sum of length of all cell edges
 * This CellFunction returns the sum of the lengths for all of the edges
 * in a cell. The input mesh is set via method SetInputMesh().
 *
 * \ingroup CellFunctions
 *
 * \ingroup ITKMeshFunction
 */
template< typename TInputMesh, typename TInputImage, typename TOutput=float >
class ITK_TEMPLATE_EXPORT CellImageValueSummaryCellFunction:
  public CellFunction< TInputMesh, TOutput >
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(CellImageValueSummaryCellFunction);

  /** Standard class type aliases. */
  using Self = CellImageValueSummaryCellFunction;
  using Superclass = CellFunction< TInputMesh, TOutput >;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  /** Run-time type information (and related methods). */
  itkTypeMacro(CellImageValueSummaryCellFunction, ImageFunction);

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

  TOutput Evaluate(const IndexType & index) const override;

  void SetImage( ImagePointerType image ) { m_Image = image; }

  enum SummaryMeasure { MEAN, MEDIAN, MIN, MAX };

  // FIXME - change this to an enum
  void SetMeasure( SummaryMeasure measure ) { m_Measure=measure; }




protected:
  CellImageValueSummaryCellFunction();
  ~CellImageValueSummaryCellFunction() override = default;
  void PrintSelf(std::ostream & os, Indent indent) const override;

  ImagePointerType m_Image;

  SummaryMeasure m_Measure;

private:

};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkCellImageValueSummaryCellFunction.hxx"
#endif

#endif
