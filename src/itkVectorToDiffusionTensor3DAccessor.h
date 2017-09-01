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

#ifndef itkVectorToDiffusionTensor3DAccessor_h
#define itkVectorToDiffusionTensor3DAccessor_h
#include "itkDiffusionTensor3D.h"
#include "itkVariableLengthVector.h"

namespace itk
{
 namespace Accessor
 {
 template< typename T >
 class VectorToDiffusionTensor3DAccessor
 {
 public:
   typedef   VectorToDiffusionTensor3DAccessor Self;

   typedef   DiffusionTensor3D< T > ExternalType;

   typedef  VariableLengthVector< T > InternalType;

   inline void Set(InternalType & output, const ExternalType & input) const
   {
          output.SetData( input.GetDataPointer() );
   }

   inline ExternalType Get(const InternalType & input) const
   {
     ExternalType dti( input.GetDataPointer() );
     return dti;
   }

 private:
 };
 }  // end namespace Accessor
 }  // end namespace itk

#endif
