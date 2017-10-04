/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkLabeledImageToPointSetFilter_h
#define itkLabeledImageToPointSetFilter_h

#include "itkImageToMeshFilter.h"
#include "itkImageRegionIteratorWithIndex.h"

namespace itk
{
/** \class LabeledImageToPointSetFilter
 * \brief Generate a PointSet containing center points of all
 * non-zero voxels.
 *
 * LabeledImageToPointSetFilter takes a binary image as input
 * and generates a PointSet as output. The point set contains
 * points that are non-zero in the binary mask.
 *
 * This filter is intended to be used for generating seeds for
 * deterministic DTI tractography
 *
 * The filter is templated over the input image type and the
 * output mesh type. The only restriction is that the dimension
 * of points in the mesh should be equal to the input image dimension.
 *
 * \sa ReinitializeImageFilter
 * \sa PointSetToImageRegistrationMethod
 *
 * \ingroup ImageFilters  MeshFilters
 * \ingroup ITKLevelSets
 */
template< typename TInputImage, typename TOutputMesh >
class LabeledImageToPointSetFilter:
  public ImageToMeshFilter< TInputImage, TOutputMesh >
{
public:
  /** Standard class typedefs. */
  typedef LabeledImageToPointSetFilter Self;

  typedef ImageToMeshFilter< TInputImage, TOutputMesh >  Superclass;

  typedef SmartPointer< Self >                 Pointer;
  typedef SmartPointer< const Self >           ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(LabeledImageToPointSetFilter, ImageToMeshFilter);

  /** Some typedefs associated with the input images. */
  typedef TInputImage                           InputImageType;
  typedef typename InputImageType::ConstPointer InputImageConstPointer;
  typedef typename InputImageType::RegionType   InputImageRegionType;
  typedef typename InputImageType::PixelType    InputImagePixelType;

  typedef ImageRegionConstIteratorWithIndex< InputImageType > InputImageIterator;

  /** Some typedefs associated with the output mesh. */
  typedef TOutputMesh                                 OutputMeshType;
  typedef typename OutputMeshType::PointType          PointType;
  typedef typename OutputMeshType::Pointer            OutputMeshPointer;
  typedef typename OutputMeshType::ConstPointer       OutputMeshConstPointer;
  typedef typename OutputMeshType::PointsContainer    PointsContainer;
  typedef typename OutputMeshType::PointIdentifier    PointIdentifier;
  typedef typename PointsContainer::Pointer           PointsContainerPointer;
  typedef typename PointsContainer::Iterator          PointsContainerIterator;

  /** Image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Float image type to be used by the ReinitializeLevelSet image filter */
  typedef itk::Image< float,
                      itkGetStaticConstMacro(ImageDimension) >   RealImageType;

  /** The dimension of the output mesh. */
  itkStaticConstMacro(PointDimension, unsigned int,
                      TOutputMesh::PointDimension);

  /** Some typedefs associated with the output mesh. */
  void GenerateData(void) ITK_OVERRIDE;

  /** Some typedefs associated with the output mesh. */
  void GenerateOutputInformation(void) ITK_OVERRIDE;

  /** accept the input image */
  using Superclass::SetInput;
  void SetInput(const InputImageType *inputImage);

  itkSetMacro(Label, InputImagePixelType);
  itkGetMacro(Label, InputImagePixelType);

protected:
  LabeledImageToPointSetFilter();
  ~LabeledImageToPointSetFilter();
  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

private:
  LabeledImageToPointSetFilter(const LabeledImageToPointSetFilter &) ITK_DELETE_FUNCTION;
  void operator=(const LabeledImageToPointSetFilter &) ITK_DELETE_FUNCTION;

  InputImagePixelType   m_Label;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLabeledImageToPointSetFilter.hxx"
#endif

#endif
