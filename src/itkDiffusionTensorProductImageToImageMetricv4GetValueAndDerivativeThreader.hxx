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
#ifndef itkDiffusionTensorProductImageToImageMetricv4GetValueAndDerivativeThreader_hxx
#define itkDiffusionTensorProductImageToImageMetricv4GetValueAndDerivativeThreader_hxx

#include "itkDiffusionTensorProductImageToImageMetricv4GetValueAndDerivativeThreader.h"
#include "itkDefaultConvertPixelTraits.h"

namespace itk
{

template< typename TDomainPartitioner, typename TImageToImageMetric, typename TDiffusionTensorProductMetric >
typename DiffusionTensorProductImageToImageMetricv4GetValueAndDerivativeThreader< TDomainPartitioner, TImageToImageMetric, TDiffusionTensorProductMetric >::MeasureType
DiffusionTensorProductImageToImageMetricv4GetValueAndDerivativeThreader< TDomainPartitioner, TImageToImageMetric, TDiffusionTensorProductMetric >
::DiffusionTensorProduct(
    const FixedImagePixelType & dt1,
    const MovingImagePixelType & dt2 ) const
{
    MeasureType value = 0;
    // for each pair: (sqrt(x)(a,b,c) . sqrt(y)(d,e,f))^2
    //a*a d*d x y + 2 a b d e x y + 2 a c d f x y + b*b e*e x y + 2 b c e f x y + c*c f*f x y
    for (unsigned int i=0; i<3; i++) {
      for (unsigned int j=i; j<3; j++ ) {

        unsigned int iOffset = 3 + i*3;
        unsigned int jOffset = 3 + j*3;

        MeasureType x = dt1[i];
        MeasureType y = dt2[j];
        MeasureType a = dt1[iOffset];
        MeasureType b = dt1[iOffset+1];
        MeasureType c = dt1[iOffset+2];
        MeasureType d = dt2[jOffset];
        MeasureType e = dt2[jOffset+1];
        MeasureType f = dt2[jOffset+2];

        MeasureType pairValue = a*a*d*d*x*y
          + 2.0*a*b*d*e*x*y
          + 2.0*a*c*d*f*x*y
          + b*b*e*e*x*y
          + 2.0*b*c*e*f*x*y
          + c*c*f*f*x*y;

        value += pairValue;

        // double for off-diagonal elements
        if ( i != j ) {
          value += pairValue;
        }
      }
    }
    return value;
}

template< typename TDomainPartitioner, typename TImageToImageMetric, typename TDiffusionTensorProductMetric >
bool
DiffusionTensorProductImageToImageMetricv4GetValueAndDerivativeThreader< TDomainPartitioner, TImageToImageMetric, TDiffusionTensorProductMetric >
::ProcessPoint( const VirtualIndexType &,
                const VirtualPointType &           virtualPoint,
                const FixedImagePointType &,
                const FixedImagePixelType &        fixedImageValue,
                const FixedImageGradientType &,
                const MovingImagePointType &       ,
                const MovingImagePixelType &       movingImageValue,
                const MovingImageGradientType &    movingImageGradient,
                MeasureType &                      metricValueReturn,
                DerivativeType &                   localDerivativeReturn,
                const ThreadIdType                 threadId) const
{
  /** Only the voxelwise contribution given the point pairs. */
  FixedImagePixelType diff = fixedImageValue - movingImageValue;
  const unsigned int nComponents = NumericTraits<FixedImagePixelType>::GetLength( diff );
  metricValueReturn = NumericTraits<MeasureType>::ZeroValue();

  metricValueReturn = this->DiffusionTensorProduct( fixedImageValue, movingImageValue );
  //for ( unsigned int nc = 0; nc < nComponents; nc++ )
  //  {
  //  MeasureType diffC = DefaultConvertPixelTraits<FixedImagePixelType>::GetNthComponent(nc, diff);
  //  metricValueReturn += diffC*diffC;
  //  }

  if( ! this->GetComputeDerivative() )
    {
    return true;
    }

  /* Use a pre-allocated jacobian object for efficiency */
  typedef typename TImageToImageMetric::JacobianType & JacobianReferenceType;
  JacobianReferenceType jacobian = this->m_GetValueAndDerivativePerThreadVariables[threadId].MovingTransformJacobian;
  JacobianReferenceType jacobianPositional = this->m_GetValueAndDerivativePerThreadVariables[threadId].MovingTransformJacobianPositional;

  /** For dense transforms, this returns identity */
  this->m_Associate->GetMovingTransform()->
    ComputeJacobianWithRespectToParametersCachedTemporaries(virtualPoint,
                                                            jacobian,
                                                            jacobianPositional);

  for ( unsigned int par = 0; par < this->GetCachedNumberOfLocalParameters(); par++ )
    {
    localDerivativeReturn[par] = NumericTraits<DerivativeValueType>::ZeroValue();
    for ( unsigned int nc = 0; nc < nComponents; nc++ )
    {
      //MeasureType diffValue = DefaultConvertPixelTraits<FixedImagePixelType>::GetNthComponent(nc,diff);
      MeasureType dValue = NumericTraits<DerivativeValueType>::ZeroValue();
      for ( SizeValueType dim = 0; dim < ImageToImageMetricv4Type::MovingImageDimension; dim++ )
        {
        unsigned int idx = nc/4;      // which eigenvalue/vector pair
        unsigned int comp = nc-idx;
        unsigned int offset = idx*3;

        MeasureType x = fixedImageValue[idx];
        MeasureType y = movingImageValue[idx];
        MeasureType a = fixedImageValue[idx+1];
        MeasureType b = fixedImageValue[idx+2];
        MeasureType c = fixedImageValue[idx+3];
        MeasureType d = movingImageValue[idx+1];
        MeasureType e = movingImageValue[idx+2];
        MeasureType f = movingImageValue[idx+3];

        MeasureType dValue = NumericTraits<DerivativeValueType>::ZeroValue();

        switch(comp) {
          case 1:  //ddy
            dValue = a*a*d*d*x + 2.0*a*b*d*e*x + 2.0*a*c*d*f*x + b*b*e*e*x + 2.0*b*c*e*f*x + c*c*f*f*x;
            break;
          case 2: //ddd
            dValue = 2.0*a*a*d*x*y + 2.0*a*b*e*x*y + 2.0*a*c*f*x*y;
            break;
          case 3: //dde
            dValue = 2.0*a*b*d*x*y + 2.0*b*b*e*x*y + 2.0*b*c*f*x*y;
            break;
          case 4: //ddf
            dValue = 2.0*a*c*d*x*y + 2.0*b*c*e*x*y + 2.0*c*c*f*x*y;
            break;
          }

        localDerivativeReturn[par] += dValue * jacobian(dim, par) *
          DefaultConvertPixelTraits<MovingImageGradientType>::GetNthComponent(
            ImageToImageMetricv4Type::FixedImageDimension * nc + dim, movingImageGradient );
        }
      }
    }
  return true;
}

} // end namespace itk

#endif
