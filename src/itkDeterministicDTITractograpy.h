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
#ifndef itkDeterministicDTITractography_h
#define itkDeterministicDTITractography_h

#include "vnl/vnl_matrix_fixed.h"
#include "itkMesh.h"
#include "itkImageToMeshFilter.h"
#include "itkTriangleCell.h"
#include "itkLineCell.h"
#include "itkCovariantVector.h"
#include "itkDefaultStaticMeshTraits.h"
#include "itkImageRegionConstIterator.h"

namespace itk
{
/** \class DeterministicDTITractography
 *
 *
 * \par
 * This class tries to construct a 3D mesh based on a vector image.
 *
 * \par PARAMETERS
 *  *
 * \par REFERENCE
 *
 *
 * \par INPUT
 * The input should be a 3D vector image.
 *
 * \ingroup ITKMesh
 */
template< typename TInputImage, typename TOutputMesh, typename TMaskImage >
class DeterministicDTITractography:public ImageToMeshFilter< TInputImage, TOutputMesh >
{
public:
  /** Standard "Self" typedef. */
  typedef DeterministicDTITractography                  Self;
  typedef ImageToMeshFilter< TInputImage, TOutputMesh > Superclass;
  typedef SmartPointer< Self >                          Pointer;
  typedef SmartPointer< const Self >                    ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(DeterministicDTITractography, ImageToMeshFilter);

  /** Hold on to the type information specified by the template parameters. */
  typedef TOutputMesh                         OutputMeshType;
  typedef typename OutputMeshType::MeshTraits OMeshTraits;
  typedef typename OutputMeshType::PointType  OPointType;
  typedef typename OMeshTraits::PixelType     OPixelType;

  /** Some convenient typedefs. */
  typedef typename OutputMeshType::Pointer                OutputMeshPointer;
  typedef typename OutputMeshType::CellTraits             CellTraits;
  typedef typename OutputMeshType::PointsContainerPointer PointsContainerPointer;
  typedef typename OutputMeshType::PointsContainer        PointsContainer;
  typedef typename OutputMeshType::CellsContainerPointer  CellsContainerPointer;
  typedef typename OutputMeshType::CellsContainer         CellsContainer;
  typedef CovariantVector< double, 2 >                    doubleVector;
  typedef CovariantVector< int, 2 >                       intVector;

  /** Define the triangular cell types which forms the surface of the model
   * and will be used in FEM application. */
  typedef CellInterface< OPixelType, CellTraits > TCellInterface;
  typedef LineCell< TCellInterface >              LineCell;
  typedef typename LineCell::SelfAutoPointer      LineCellAutoPointer;

  typedef TriangleCell< TCellInterface >          TriCell;
  typedef typename TriCell::SelfAutoPointer       TriCellAutoPointer;

  /** Input Image Type Definition. */
  typedef TInputImage                           InputImageType;
  typedef typename InputImageType::Pointer      InputImagePointer;
  typedef typename InputImageType::ConstPointer InputImageConstPointer;
  typedef typename InputImageType::PixelType    InputPixelType;
  typedef typename InputImageType::SpacingType  SpacingType;
  typedef typename InputImageType::PointType    OriginType;
  typedef typename InputImageType::RegionType   RegionType;
  typedef typename InputImageType::SizeType     SizeType;

  /** Type definition for the classified image index type. */
  typedef typename InputImageType::IndexType           InputImageIndexType;

  typedef ImageRegionConstIterator< InputImageType > InputImageIterator;

  typedef TMaskImage                        MaskImageType;
  typedef typename MaskImageType::Pointer   MaskImagePointer;

  typedef itk::IdentifierType                   IdentifierType;
  typedef itk::SizeValueType                    SizeValueType;

  itkSetMacro(ObjectValue, InputPixelType);

  itkGetMacro(Directed, bool);
  itkSetMacro(Directed, bool)

  itkGetConstMacro(NumberOfNodes, SizeValueType);
  itkGetConstMacro(NumberOfCells, SizeValueType);

  /** accept the input image */
  using Superclass::SetInput;
  virtual void SetInput(const InputImageType *inputImage);

  void SetSeeds(const OutputMeshType *seeds);

  itkSetObjectMacro(ValidRegion, MaskImagePointer);




protected:
  DeterministicDTITractography();
  ~DeterministicDTITractography();
  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

  void GenerateData() ITK_OVERRIDE;

  bool ValidPoint( OPointType point );

  PointsContainer Track( OPointType seed, bool forward);

  bool       m_Directed;

  virtual void GenerateOutputInformation() ITK_OVERRIDE {}  // do nothing ITK_OVERRIDE

private:
  DeterministicDTITractography(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;

  typedef typename InputImageType::SizeType InputImageSizeType;

  MaskImagePointer m_ValidRegion;

  /** temporary variables used in CreateMesh to avoid thousands of
   *  calls to GetInput() and GetOutput()
   */
  OutputMeshType       *m_OutputMesh;
  const InputImageType *m_InputImage;
  const OutputMeshType *m_SeedMesh;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDeterministicDTITractography.hxx"
#endif

#endif
