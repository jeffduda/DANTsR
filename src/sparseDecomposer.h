

#ifndef __sparseDecomposer_h
#define __sparseDecomposer_h

#include <RcppDANTsR.h>
#include "itkImage.h"

namespace itk {

class sparseDecomposerOne
{
public:

  using Self = sparseDecomposerOne;
  using Pointer = SmartPointer<Self>
  using ConstPointer = SmartPointer<const Self>

  /** Method for creation through the object factory. */
  itkNewMacro( Self );



}

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "sparseDecomposer.hxx"
#endif
