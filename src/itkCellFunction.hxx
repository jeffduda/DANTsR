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
#ifndef itkCellFunction_hxx
#define itkCellFunction_hxx

#include "itkCellFunction.h"

namespace itk
{
/**
 * Constructor
 */
template< typename TInputMesh, typename TOutput >
CellFunction< TInputMesh, TOutput >
::CellFunction()
{
  m_Mesh = nullptr;
  m_StartIndex = 0;
  m_EndIndex = 0;
}

/**
 * Standard "PrintSelf" method
 */
template< typename TInputMesh, typename TOutput >
void
CellFunction< TInputMesh, TOutput >
::PrintSelf(
  std::ostream & os,
  Indent indent) const
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "InputMesh: " << m_Mesh.GetPointer() << std::endl;
  os << indent << "StartIndex: " << m_StartIndex << std::endl;
  os << indent << "EndIndex: " << m_EndIndex << std::endl;
}

/**
 * Initialize by setting the input image
 */
template< typename TInputMesh, typename TOutput >
void
CellFunction< TInputMesh, TOutput >
::SetInputMesh(
  const InputMeshType *ptr)
{
  // set the input image
  m_Mesh = ptr;

  if ( ptr )
    {
    m_StartIndex = 0;
    m_EndIndex = m_Mesh->GetNumberOfCells()-1;
    }
}
} // end namespace itk

#endif
