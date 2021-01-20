#ifndef itkCellCountImageFilter_h
#define itkCellCountImageFilter_h

#include "itkImageSource.h"

#include "itkAddImageFilter.h"
#include "itkPolyLineCell.h"
#include "itkMapContainer.h"
#include "itkVectorContainer.h"
#include "itkAutomaticTopologyMeshSource.h"
#include "itkPointSet.h"
#include "itkContinuousIndex.h"

#include <vector>

namespace itk
{
class Point1D
{
public:
  double m_X;
  int    m_Sign;

  Point1D(){}
  Point1D(const double p, const int s)
  {
    m_X = p;
    m_Sign = s;
  }

  Point1D(const Point1D & point)
  {
    m_X = point.m_X;
    m_Sign = point.m_Sign;
  }

  double getX() const
  {
    return m_X;
  }

  int  getSign() const
  {
    return m_Sign;
  }
};

/** \class CellCountImageFilter
 *
 * \brief 3D Rasterization algorithm Courtesy of Dr David Gobbi of Atamai Inc.

 * \author Leila Baghdadi, MICe, Hospital for Sick Childern, Toronto, Canada,
 * \ingroup ITKMesh
 */
template< typename TInputMesh, typename TOutputImage  >
class ITK_TEMPLATE_EXPORT CellCountImageFilter:public ImageSource< TOutputImage >
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(CellCountImageFilter);

  /** Standard class type aliases. */
  using Self = CellCountImageFilter;
  using Superclass = ImageSource< TOutputImage >;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  //using IndexType = typename TOutputImage::IndexType;
  //using SizeType = typename TOutputImage::SizeType;
  using OutputImageType = TOutputImage;
  using OutputImagePointer = typename OutputImageType::Pointer;
  using ValueType = typename OutputImageType::ValueType;
  //using SpacingType = typename OutputImageType::SpacingType;
  //using DirectionType = typename OutputImageType::DirectionType;
  using ContinuousIndexType = itk::ContinuousIndex< ValueType, OutputImageType::ImageDimension>;

  using AddFilterType = itk::AddImageFilter<OutputImageType, OutputImageType, OutputImageType>;
  using AddFilterPointerType = typename AddFilterType::Pointer;

  using SubsetType = typename std::vector<double>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(CellCountImageFilter, ImageSource);

  /** Superclass type alias. */
  using OutputImageRegionType = typename Superclass::OutputImageRegionType;

  /** Some convenient type alias. */
  using InputMeshType = TInputMesh;
  using InputMeshPointer = typename InputMeshType::Pointer;
  using InputPointType = typename InputMeshType::PointType;
  using InputPixelType = typename InputMeshType::PixelType;
  using InputCellTraitsType = typename InputMeshType::MeshTraits::CellTraits;
  using CellType = typename InputMeshType::CellType;
  using CellsContainerPointer = typename InputMeshType::CellsContainerPointer;
  using CellsContainerIterator = typename InputMeshType::CellsContainerIterator;
  using CellAutoPointer = typename CellType::CellAutoPointer;

  //using InfoImageType = itk::ImageBase< OutputImageType::ImageDimension >;
  using InfoImageType = OutputImageType;
  using InfoImagePointer = typename InfoImageType::Pointer;
  using IndexType = typename InfoImageType::IndexType;
  using SizeType = typename InfoImageType::SizeType;
  using SpacingType = typename InfoImageType::SpacingType;
  using DirectionType = typename InfoImageType::DirectionType;
  using PointType = typename InfoImageType::PointType;

  using InputPointsContainer = typename InputMeshType::PointsContainer;
  using InputPointsContainerPointer = typename InputPointsContainer::Pointer;
  using InputPointsContainerIterator = typename InputPointsContainer::Iterator;

  using PointSetType = itk::PointSet< double, 3 >;
  using PointsContainer = typename PointSetType::PointsContainer;

  //using PointType = itk::Point< double, 3 >;
  using Point2DType = itk::Point< double, 2 >;

  using DoubleArrayType = itk::Array< double >;

  using Point1DVector = std::vector< Point1D >;
  using Point1DArray = std::vector< std::vector< Point1D > >;

  using Point2DVector = std::vector< Point2DType >;
  using Point2DArray = std::vector< std::vector< Point2DType > >;

  using PointVector = std::vector< PointType >;
  using PointArray = std::vector< std::vector< PointType > >;

  using StencilIndexVector = std::vector< int >;
  /** Spacing (size of a pixel) of the output image. The
   * spacing is the geometric distance between image samples.
   * It is stored internally as double, but may be set from
   * float. \sa GetSpacing() */
  itkSetMacro(Spacing, SpacingType);
  virtual void SetSpacing(const double spacing[3]);

  virtual void SetSpacing(const float spacing[3]);

  itkGetConstReferenceMacro(Spacing, SpacingType);

  /** The Direction is a matix of direction cosines
   *  that specify the direction between samples.
   * */
  itkSetMacro(Direction, DirectionType);
  itkGetConstMacro(Direction, DirectionType);

  /** Set/Get the value for pixels inside the spatial object.
  * By default, this filter will return an image
  * If this "inside" value is changed to a non-null value,
  * the output produced by this filter will be a mask with inside/outside values
  * specified by the user. */
  itkSetMacro(InsideValue, ValueType);
  itkGetConstMacro(InsideValue, ValueType);

  /** Set/Get the value for pixels outside the spatial object.
  * By default, this filter will return an image
  * If this "outside" value is changed to a non-null value,
  * the output produced by this filter will be a mask with inside/outside values
  * specified by the user. */
  itkSetMacro(OutsideValue, ValueType);
  itkGetConstMacro(OutsideValue, ValueType);

  /** The origin of the output image. The origin is the geometric
   * coordinates of the index (0,0,...,0).  It is stored internally
   * as double but may be set from float.
   * \sa GetOrigin() */
  itkSetMacro(Origin, PointType);
  virtual void SetOrigin(const double origin[3]);

  virtual void SetOrigin(const float origin[3]);

  itkGetConstReferenceMacro(Origin, PointType);

  itkSetMacro(Subset, SubsetType);

  /** Set/Get Index */
  itkSetMacro(Index, IndexType);
  itkGetConstMacro(Index, IndexType);

  /** Set/Get Size */
  itkSetMacro(Size, SizeType);
  itkGetConstMacro(Size, SizeType);

  itkSetMacro(Target, bool);
  itkGetConstMacro(Target, bool);

  /** Set the mesh input of this process object.  */
  using Superclass::SetInput;
  void SetInput(InputMeshType *input);

  void SetInfoImage(InfoImageType *InfoImage)
  {
    if ( InfoImage != m_InfoImage )
      {
      this->Modified();
      m_InfoImage = InfoImage;
      }
  }


  /** Get the mesh input of this process object.  */
  InputMeshType * GetInput();

  InputMeshType * GetInput(unsigned int idx);

  /* Set the tolerance for doing spatial searches of the polydata. */
  itkSetMacro(Tolerance, double);
  itkGetConstMacro(Tolerance, double);

protected:
  CellCountImageFilter();
  ~CellCountImageFilter() override;

  void GenerateOutputInformation() override {}  // do nothing
  void GenerateData() override;

  InfoImageType *m_InfoImage;

  bool m_Target;

  IndexType m_Index;

  SizeType m_Size;

  SpacingType m_Spacing;

  PointType m_Origin;        //start value

  double m_Tolerance;

  ValueType m_InsideValue;
  ValueType m_OutsideValue;

  DirectionType m_Direction;

  StencilIndexVector m_StencilIndex;

  SubsetType m_Subset;

  void PrintSelf(std::ostream & os, Indent indent) const override;

private:

  //void CountCellVoxelHits( unsigned long );

  //OutputImagePointer m_CellMask;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkCellCountImageFilter.hxx"
#endif

#endif
