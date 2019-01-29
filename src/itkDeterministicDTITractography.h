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
#include "itkPolyLineCell.h"
#include "itkCovariantVector.h"
#include "itkDefaultStaticMeshTraits.h"
#include "itkImageRegionConstIterator.h"
#include "itkPolyLineCell.h"

namespace itk
{
/** \class DeterministicDTITractography
 *
 *
 * \par
 * This class performs particle tracking for points in a vector field. The
 * intended inputs are points in white matter to be tracked using a vector
 * field that represents the primary direction of diffusion in white matter
 * as measured with diffusion tensor magnetic resonance imaging.
 *
 * \par
 *
 * \par
 *
 * \par PARAMETERS
 *
 * \par REFERENCE
 *
*
 * \par INPUT
 * The first input should be a 3D vector field
 *
 * \ingroup ITKMesh
 */
template< typename TInputImage, typename TOutputMesh >
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
  itkTypeMacro(BinaryMask3DMeshSource, ImageToMeshFilter);

  /** Hold on to the type information specified by the template parameters. */
  using OutputMeshType = TOutputMesh;
  using OMeshTraits = typename OutputMeshType::MeshTraits;
  using OPointType = typename OutputMeshType::PointType;
  using OPixelType = typename OMeshTraits::PixelType;

  typedef typename OutputMeshType::CellType       CellType;
  typedef typename CellType::CellAutoPointer      CellAutoPointer;
  typedef typename itk::PolyLineCell< CellType >  PolyLineCellType;

  typedef itk::VectorContainer< unsigned long, unsigned long > SeedContainer;
  typedef typename SeedContainer::Pointer SeedContainerPointer;


  //typedef TOutputMesh                         OutputMeshType;
  //typedef typename OutputMeshType::MeshTraits OMeshTraits;
  //typedef typename OutputMeshType::PointType  OPointType;
  //typedef typename OMeshTraits::PixelType     OPixelType;
  typedef typename OutputMeshType::ConstPointer ConstMeshPointer;

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
  typedef CellInterface< OPixelType, CellTraits >  TCellInterface;
  typedef PolyLineCell< TCellInterface >           StreamlineCell;
  typedef typename StreamlineCell::SelfAutoPointer StreamlineCellAutoPointer;

  /** Input Image Type Definition. */
  typedef TInputImage                           InputImageType;
  typedef typename InputImageType::Pointer      InputImagePointer;
  typedef typename InputImageType::ConstPointer InputImageConstPointer;
  typedef typename InputImageType::PixelType    InputPixelType;
  typedef typename InputImageType::SpacingType  SpacingType;
  typedef typename InputImageType::PointType    PointType;
  typedef typename InputImageType::RegionType   RegionType;
  typedef typename InputImageType::SizeType     SizeType;

  typedef typename itk::Vector< typename InputImageType::PointValueType, InputImageType::ImageDimension> VectorType;

  typedef typename InputImageType::InternalPixelType ValueType;
  typedef typename itk::Image<ValueType, InputImageType::ImageDimension> MaskType;

  using MaskPointerType = typename MaskType::Pointer;
  using ConstMaskPointer = typename MaskType::ConstPointer;

  //using MaskPointerType = typename MaskType::Pointer;
  //using ConstMaskPointer = typename MaskType::ConstPointer;

  /** Type definition for the classified image index type. */
  typedef typename InputImageType::IndexType           InputImageIndexType;

  typedef ImageRegionConstIterator< InputImageType > InputImageIterator;

  typedef itk::IdentifierType                   IdentifierType;
  typedef itk::SizeValueType                    SizeValueType;

  itkGetObjectMacro(SeedOffsets, SeedContainer);

  itkSetMacro(MinimumNumberOfPoints, unsigned long);
  itkGetMacro(MinimumNumberOfPoints, unsigned long);

  itkSetMacro(MaximumNumberOfPoints, unsigned long);
  itkGetMacro(MaximumNumberOfPoints, unsigned long);

  itkSetMacro(StepSize, double);
  itkGetMacro(StepSize, double);

  itkSetMacro(RadianThreshold, double);
  itkGetMacro(RadianThreshold, double);

  itkSetMacro(DegreeThreshold, double);
  itkGetMacro(DegreeThreshold, double);

  itkSetMacro(ObjectValue, InputPixelType);

  itkGetConstMacro(NumberOfNodes, SizeValueType);
  itkGetConstMacro(NumberOfCells, SizeValueType);

  /** accept the input image */
  using Superclass::SetInput;

  void SetField(const InputImageType * inputImage);
  InputImageConstPointer GetField( void );

  void SetSeeds( const OutputMeshType * seeds);
  ConstMeshPointer GetSeeds( void );

  void SetMask( MaskPointerType mask )
    { m_Mask = mask; }

  ConstMaskPointer GetMask()
    { return m_Mask; }


  PointsContainerPointer TrackFiber( PointType point, bool forward );

  void SetRegionOfInterest( const RegionType & iRegion )
    {
    if( iRegion != m_RegionOfInterest )
      {
      this->m_RegionOfInterest = iRegion;
      this->m_RegionOfInterestProvidedByUser = true;
      this->Modified();
      }
    }

  itkGetConstReferenceMacro(RegionOfInterest, RegionType);

protected:
  DeterministicDTITractography();
  ~DeterministicDTITractography();
  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

  void GenerateData() ITK_OVERRIDE;


  bool       m_RegionOfInterestProvidedByUser;
  RegionType m_RegionOfInterest;

  virtual void GenerateOutputInformation() ITK_OVERRIDE {}  // do nothing ITK_OVERRIDE

private:
  DeterministicDTITractography(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;

  typedef typename InputImageType::SizeType InputImageSizeType;

  SizeValueType m_NumberOfNodes;
  SizeValueType m_NumberOfCells;

  ConstMeshPointer m_SeedMesh;
  ConstMaskPointer m_Mask;

  int m_NodeLimit;
  int m_CellLimit;
  int m_ImageWidth;
  int m_ImageHeight;
  int m_ImageDepth;

  unsigned long m_MinimumNumberOfPoints;
  unsigned long m_MaximumNumberOfPoints;

  double m_StepSize;
  double m_RadianThreshold;
  double m_DegreeThreshold;

  SeedContainerPointer m_SeedOffsets;

  unsigned char  m_PointFound;
  InputPixelType m_ObjectValue;

  OutputMeshPointer m_OutputMesh;
  InputImageConstPointer m_InputImage;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDeterministicDTITractography.hxx"
#endif

#endif
