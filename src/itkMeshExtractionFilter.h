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
/*=========================================================================
 *
 *  Portions of this file are subject to the VTK Toolkit Version 3 copyright.
 *
 *  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
 *
 *  For complete copyright, license and disclaimer of warranty information
 *  please refer to the NOTICE file at the top of the ITK source tree.
 *
 *=========================================================================*/
#ifndef itkMeshExtractionFilter_h
#define itkMeshExtractionFilter_h

#include "itkIntTypes.h"
#include "itkMeshToMeshFilter.h"
#include "itkLineCell.h"
#include "itkPolygonCell.h"

namespace itk
{
/** \class MeshExtractionFilter
 * \brief Extract portions of a mesh.
 *
 * MeshExtractionFilter will extract portions of a mesh
 * from a list of point ids or cell ids
 *
 * \ingroup MeshFilters
 * \ingroup ITKMesh
 */

template< typename TInputMesh, typename TOutputMesh >
class ITK_TEMPLATE_EXPORT MeshExtractionFilter:
  public MeshToMeshFilter< TInputMesh, TOutputMesh >
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(MeshExtractionFilter);

  /**
   * Standard class type aliases.
   */
  using Self = MeshExtractionFilter;

  /**
   * Standard "Superclass" type alias.
   */
  using Superclass = MeshToMeshFilter< TInputMesh, TOutputMesh >;

  /**
   * Smart pointer type alias support
   */
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  /**
   * Method for creation through the object factory.
   */
  itkNewMacro(Self);

  /**
   * Convenient type alias for this filter.
   */
  using InputMeshType = TInputMesh;
  using OutputMeshType = TOutputMesh;
  using InputMeshConstPointer = typename TInputMesh::ConstPointer;
  using OutputMeshPointer = typename TOutputMesh::Pointer;

  static constexpr unsigned int PointDimension = TInputMesh::PointDimension;

  using InputMeshPointType = typename TInputMesh::PointType;
  using InputMeshPointIdentifier = typename TInputMesh::PointIdentifier;
  using InputMeshPointsContainerConstPointer = typename TInputMesh::PointsContainerConstPointer;
  using InputMeshCellsContainer = typename TInputMesh::CellsContainer;
  using InputMeshCellsContainerPointer = typename TInputMesh::CellsContainerPointer;
  using InputMeshCellsContainerConstPointer = typename TInputMesh::CellsContainerConstPointer;
  using InputMeshCellDataContainer = typename TInputMesh::CellDataContainer;
  using InputMeshCellDataContainerPointer = typename TInputMesh::CellDataContainerPointer;
  using InputMeshCellDataContainerConstPointer = typename TInputMesh::CellDataContainerConstPointer;
  using PointsContainerConstIterator = typename InputMeshType::PointsContainer::ConstIterator;
  using CellsContainerConstIterator = typename InputMeshType::CellsContainer::ConstIterator;
  using CellDataContainerConstIterator = typename InputMeshType::CellDataContainer::ConstIterator;
  using InputMeshCellPointer = typename TInputMesh::CellAutoPointer;
  using InputMeshPointIdConstIterator = typename TInputMesh::CellTraits::PointIdConstIterator;
  using InputMeshCellLinksContainerConstPointer = typename TInputMesh::CellLinksContainerConstPointer;
  using InputMeshCellLinksContainer = typename TInputMesh::PointCellLinksContainer;
  using InputMeshCellIdentifier = typename TInputMesh::CellIdentifier;

  using InputMeshCellType = typename TInputMesh::CellType;
  using InputMeshCellAutoPointer = typename InputMeshCellType::CellAutoPointer;

  using OutputMeshCellType = typename TOutputMesh::CellType;
  using OutputMeshCellAutoPointer = typename OutputMeshCellType::CellAutoPointer;

  using LineType = LineCell<InputMeshCellType>;
  using PolyType = PolygonCell<InputMeshCellType>;

  /**
   * Different modes of operation. Use these to specify
   * how to extract the regions.
   */
  enum { PointRegions = 0,
         CellRegions = 1 };

  /**
   * Methods specify mode of operation for the filter. Note that
   * some modes require additional information. For example,
   * SetExtractionModeToClosestPointRegion() also requires that
   * a point be defined.
   */
  itkSetMacro(ExtractionMode, int);
  itkGetConstMacro(ExtractionMode, int);

  void SetExtractionModeToPointRegions(void)
  {
    this->SetExtractionMode(Self::PointRegions);
  }

  void SetExtractionModeToCellRegions(void)
  {
    this->SetExtractionMode(Self::CellRegions);
  }

  /**
   * Initialize list of point ids/cell ids used to seed regions.
   */
  void InitializePointList(void)
  {
    this->Modified();
    m_PointList.clear();
  }

  /**
   * Initialize list of point ids/cell ids used to seed regions.
   */
  void InitializeCellList(void)
  {
    this->Modified();
    m_CellList.clear();
  }

  /**
   * Add a seed id (point or cell id). Note: ids are 0-offset.
   */
  void AddPointId(IdentifierType id)
  {
    this->Modified();
    m_PointList.push_back(id);
  }

  /**
   * Add a seed id (point or cell id). Note: ids are 0-offset.
   */
  void AddCellId(IdentifierType id)
  {
    this->Modified();
    m_CellList.push_back(id);
  }

  /**
   * Delete a seed id (point or cell id). Note: ids are 0-offset.
   */
  //void DeleteId(IdentifierType id);

  /**
   * Initialize list of region ids to extract.
   */
  //void InitializeSpecifiedRegionList(void)
  //{
  //  this->Modified();
  //  m_RegionList.clear();
  //}

  /**
   * Add a region id to extract. Note: ids are 0-offset.
   */
  //void AddSpecifiedRegion(IdentifierType id)
  //{
  //  this->Modified();
  //  m_RegionList.push_back(id);
  //}

  /**
   * Delete a region id to extract. Note: ids are 0-offset.
   */
  //void DeleteSpecifiedRegion(IdentifierType id);



protected:

  MeshExtractionFilter();
  ~MeshExtractionFilter() override {}

  void PrintSelf(std::ostream & os, Indent indent) const override;

  void GenerateData() override;

private:

  int                            m_ExtractionMode;
  std::vector< IdentifierType >  m_IdList;
  std::vector< IdentifierType >  m_RegionList;
  std::vector< SizeValueType >   m_RegionSizes;

  std::vector< IdentifierType >  m_PointList;
  std::vector< IdentifierType >  m_CellList;

  std::vector< OffsetValueType > m_Visited;
  SizeValueType                  m_NumberOfCellsInRegion;
  IdentifierType                 m_RegionNumber;

}; // class declaration
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMeshExtractionFilter.hxx"
#endif

#endif
