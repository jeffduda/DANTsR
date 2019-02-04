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
#ifndef itkEdgeLengthCellFunction_h
#define itkEdgeLengthCellFunction_h

#include "itkCellFunction.h"

namespace itk
{
/** \class EdgeLengthCellFunction
 * \brief Returns sum of length of all cell edges
 * This CellFunction returns the sum of the lengths for all of the edges
 * in a cell. The input mesh is set via method SetInputMesh().
 *
 * \ingroup CellFunctions
 *
 * \ingroup ITKImageFunction
 */
template< typename TInputMesh, typename TOutput = float >
class ITK_TEMPLATE_EXPORT EdgeLengthCellFunction:
  public CellFunction< TInputMesh, TOutput >
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(EdgeLengthCellFunction);

  /** Standard class type aliases. */
  using Self = EdgeLengthCellFunction;
  using Superclass = CellFunction< TInputMesh, TOutput >;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  /** Run-time type information (and related methods). */
  itkTypeMacro(EdgeLengthCellFunction, ImageFunction);

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

  /** BinaryThreshold the image at a point position
   *
   * Returns true if the image intensity at the specified point position
   * satisfies the threshold criteria.  The point is assumed to lie within
   * the image buffer.
   *
   * ImageFunction::IsInsideBuffer() can be used to check bounds before
   * calling the method. */

  TOutput Evaluate(const IndexType & index) const override;

protected:
  EdgeLengthCellFunction();
  ~EdgeLengthCellFunction() override = default;
  void PrintSelf(std::ostream & os, Indent indent) const override;

private:

};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkEdgeLengthCellFunction.hxx"
#endif

#endif
