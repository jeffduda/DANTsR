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
#ifndef itkCellFunction_h
#define itkCellFunction_h

#include "itkFunctionBase.h"
#include "itkMesh.h"

namespace itk
{
/** \class CellFunction
 * \brief Evaluates a function of an mesh at specified position.
 *
 * CellFunction is a baseclass for all objects that evaluate
 * a function of an mesh at a cell index.
 * This class is templated over the input mesh type, adn the type
 * of the function output.
 *
 * The input mesh is set via method SetInputMesh().
 * Method EvaluateAtCellIndex evaluates the function at an
 * cell index.
 *
 * \sa Index
 *
 * \ingroup CellFunctions
 * \ingroup ITKCellFunction
 */
template<
  typename TInputMesh,
  typename TOutput>
class ITK_TEMPLATE_EXPORT CellFunction:
    public FunctionBase< typename TInputMesh::CellIdentifier, TOutput >
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(CellFunction);

  /** Dimension underlying input image. */
  //static constexpr unsigned int Dimension = TInputMesh::Dimension;

  /** Standard class type aliases. */
  using Self = CellFunction;

  using Superclass = FunctionBase<typename TInputMesh::CellIdentifier,
    TOutput>;

  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  /** Run-time type information (and related methods). */
  itkTypeMacro(CellFunction, FunctionBase);

  /** InputMeshType type alias support */
  using InputMeshType = TInputMesh;

  /** InputPixel type alias support */
  using InputPixelType = typename InputMeshType::PixelType;

  /** InputMeshPointer type alias support */
  using InputMeshConstPointer = typename InputMeshType::ConstPointer;

  /** OutputType type alias support */
  using OutputType = TOutput;

  /** Index Type. */
  using IndexType = typename TInputMesh::CellIdentifier;


  /** Point Type. */
  using PointType = typename InputMeshType::PointType;

  /** Set the input image.
   * \warning this method caches BufferedRegion information.
   * If the BufferedRegion has changed, user must call
   * SetInputMesh again to update cached values. */
  virtual void SetInputMesh(const InputMeshType *ptr);

  /** Get the input image. */
  const InputMeshType * GetInputMesh() const
  { return m_Mesh.GetPointer(); }

  /** Evaluate the function at specified Point position.
   * Subclasses must provide this method. */
  TOutput Evaluate(const IndexType & point) const override = 0;

  /** Check if an index is inside the image buffer.
   * We take into account the fact that each voxel has its
   * center at the integer coordinate and extends half way
   * to the next integer coordinate.
   * \warning For efficiency, no validity checking of
   * the input image is done. */
  virtual bool IsInsideBuffer(const IndexType & index) const
  {
    bool isInside = m_Mesh->GetCells()->IndexExists(index);
    return isInside;
  }

  itkGetConstReferenceMacro(StartIndex, IndexType);
  itkGetConstReferenceMacro(EndIndex, IndexType);

protected:
  CellFunction();
  ~CellFunction() override = default;
  void PrintSelf(std::ostream & os, Indent indent) const override;

  /** Const pointer to the input image. */
  InputMeshConstPointer m_Mesh;

  /** Cache some values for testing if indices are inside buffered region. */
  IndexType m_StartIndex;
  IndexType m_EndIndex;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkCellFunction.hxx"
#endif

#endif
