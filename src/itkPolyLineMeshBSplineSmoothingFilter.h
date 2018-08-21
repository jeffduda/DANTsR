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
#ifndef itkPolyLineMeshBSplineSmoothingFilter_h
#define itkPolyLineMeshBSplineSmoothingFilter_h

#include "itkMeshToMeshFilter.h"
#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkBSplineControlPointImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkMesh.h"
#include "itkImage.h"

namespace itk
{
/** \class PolyLineMeshBSplineSmoothingFilter
 * \brief
 *
 * PolyLineMeshBSplineSmoothingFilter applies a transform to all the points
 * of a mesh.
 *
 * The additional content of the mesh is passed untouched. Including the
 * connectivity and the additional information contained on cells and points.
 *
 * Meshes that have added information like normal vector on the points, will
 * have to take care of transforming this data by other means.
 *
 * \ingroup MeshFilters
 * \ingroup ITKMesh
 */
template< typename TInputMesh, typename TOutputMesh >
class ITK_TEMPLATE_EXPORT PolyLineMeshBSplineSmoothingFilter:
  public MeshToMeshFilter< TInputMesh, TOutputMesh >
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(PolyLineMeshBSplineSmoothingFilter);

  /** Standard class type aliases. */
  using Self = PolyLineMeshBSplineSmoothingFilter;
  using Superclass = MeshToMeshFilter< TInputMesh, TOutputMesh >;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  using InputMeshType = TInputMesh;
  using OutputMeshType = TOutputMesh;
  using InputMeshPointer = typename InputMeshType::Pointer;
  using OutputMeshPointer = typename OutputMeshType::Pointer;

  using InputPointType = typename InputMeshType::PointType;
  using InputPixelType = typename InputMeshType::PixelType;
  using InputCellTraitsType = typename InputMeshType::MeshTraits::CellTraits;
  using CellType = typename InputMeshType::CellType;
  using CellsContainerPointer = typename InputMeshType::CellsContainerPointer;
  using CellsContainerIterator = typename InputMeshType::CellsContainerIterator;
  using CellAutoPointer = typename CellType::CellAutoPointer;

  using OutputPointType = typename OutputMeshType::PointType;
  using OutputPixelType = typename OutputMeshType::PixelType;
  //using OutputPixelValueType = typename OutputPixelType::ValueType;

  /** Type for representing coordinates. */
  using CoordRepType = typename InputMeshType::CoordRepType;

  using ParametricPixelType = typename itk::Vector<CoordRepType,InputMeshType::PointDimension>;
  //using ParametricPixelType = InputPointType;
  using ParametricMeshType = itk::Mesh<ParametricPixelType, 1>;
  using ParametricMeshPointer = typename ParametricMeshType::Pointer;
  using ParametricPointType = typename ParametricMeshType::PointType;
  using ParametricImageType = itk::Image<ParametricPixelType, 1>;
  using BSplineFilterType = BSplineScatteredDataPointSetToImageFilter<ParametricMeshType, ParametricImageType>;
  using ControlPointFilterType = BSplineControlPointImageFilter< ParametricImageType >;
  using BSplineImageWriterType = ImageFileWriter< ParametricImageType >;
  using ParametricImagePointer = typename ParametricImageType::Pointer;

  using ParametricMeshPointContainer = typename ParametricMeshType::PointsContainer;
  using ParmetricMeshPointsContainerPointer = typename ParametricMeshPointContainer::Pointer;

  itkSetMacro(UsePointDataAsParameter, bool);

  itkGetMacro(UsePointDataAsParameter, bool);
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(PolyLineMeshBSplineSmoothingFilter, MeshToMeshFilter);

protected:
  PolyLineMeshBSplineSmoothingFilter();
  ~PolyLineMeshBSplineSmoothingFilter() override {}
  void PrintSelf(std::ostream & os, Indent indent) const override;

  bool m_UsePointDataAsParameter;

  /** Generate Requested Data */
  void GenerateData() override;

};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPolyLineMeshBSplineSmoothingFilter.hxx"
#endif

#endif
