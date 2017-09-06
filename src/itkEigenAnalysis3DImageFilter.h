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
#ifndef itkEigenAnalysis3DImageFilter_h
#define itkEigenAnalysis3DImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkSymmetricSecondRankTensor.h"

namespace itk
{
/** \class EigenAnalysis3DImageFilter
 * \brief Computes pixel-wise the eigen values and eigen vectors
 *        of a 3D symmetrical matrix.
 *
 * The filter expects an image with 6 component pixels [A,B,C,D,E,F] as the
 * the components of the matrix
 *
 *                    | A  B  C |
 *                    | B  D  E |
 *                    | C  E  F |
 *
 * The eigen values are stored in three output images, and the eigen
 * vectors are  stored in images using 3D vector as pixel type.
 *
 * \ingroup ShouldBeThreaded IntensityImageFilters
 * \ingroup ITKEigen
 */

template< typename TInputImage, typename TEigenValueImage, typename TEigenVectorImage >
class EigenAnalysis3DImageFilter:
  public ImageToImageFilter< TInputImage, TEigenValueImage >
{
public:
  /** Standard class typedefs. */
  typedef EigenAnalysis3DImageFilter                          Self;
  typedef ImageToImageFilter< TInputImage, TEigenValueImage > Superclass;
  typedef SmartPointer< Self >                                Pointer;
  typedef SmartPointer< const Self >                          ConstPointer;

  /** Run-time type information (and related methods).   */
  itkTypeMacro(EigenAnalysis3DImageFilter, ImageToImageFilter);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  typedef TInputImage InputImageType;
  typedef typename InputImageType::PixelType InputImagePixelType;
  typedef typename InputImagePixelType::ValueType InputImageValueType;

  typedef typename itk::SymmetricSecondRankTensor<InputImageValueType,3> TensorType;

  /** Typedef for the vector type representing the eigen vectors */
  typedef typename TEigenVectorImage::PixelType EigenVectorType;
  typedef typename EigenVectorType::ValueType   VectorComponentType;

  /** Superclass typedefs. */
  typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

  /** Some convenient typedefs. */
  typedef TEigenValueImage                          EigenValueImageType;
  typedef typename EigenValueImageType::Pointer     EigenValueImagePointer;
  typedef typename EigenValueImageType::RegionType  EigenValueImageRegionType;
  typedef typename EigenValueImageType::PixelType   EigenValueImagePixelType;
  typedef TEigenVectorImage                         EigenVectorImageType;
  typedef typename EigenVectorImageType::Pointer    EigenVectorImagePointer;
  typedef typename EigenVectorImageType::RegionType EigenVectorImageRegionType;
  typedef typename EigenVectorImageType::PixelType  EigenVectorImagePixelType;

  /** Image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Connect the image containting the elements [0,0]
   * of the input 3D matrix */
  //void SetInput(TInputImage *image);

  /** Get the Output image eigenvalue */
  EigenValueImageType * GetEigenValue(unsigned int);

  /** Get the Output image with the eigen vector */
  EigenVectorImageType * GetEigenVector(unsigned int);

  /** Get the Output image with the greatest eigenvalue */
  EigenValueImageType * GetMaxEigenValue();

  /** Get the Output image with the smallest eigenvalue */
  EigenValueImageType * GetMinEigenValue();

  /** Get the Output image with the eigen vector associated with
   * the greatest eigen value */
  EigenVectorImageType * GetMaxEigenVector();

  /**  Create the Output */
  typedef ProcessObject::DataObjectPointerArraySizeType DataObjectPointerArraySizeType;
  using Superclass::MakeOutput;
  DataObject::Pointer MakeOutput(DataObjectPointerArraySizeType idx) ITK_OVERRIDE;

#ifdef ITK_USE_CONCEPT_CHECKING
  // Begin concept checking
  itkConceptMacro( VectorComponentHasNumericTraitsCheck,
                   ( Concept::HasNumericTraits< VectorComponentType > ) );
  // End concept checking
#endif

protected:
  EigenAnalysis3DImageFilter();
  virtual ~EigenAnalysis3DImageFilter() {}

  void GenerateData(void) ITK_OVERRIDE;

private:
  EigenAnalysis3DImageFilter(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkEigenAnalysis3DImageFilter.hxx"
#endif

#endif
