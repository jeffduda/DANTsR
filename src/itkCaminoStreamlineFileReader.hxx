/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkCaminoStreamlineFileReader.hxx,v $
  Language:  C++
  Date:      $Date: 2009/03/13 19:48:16 $
  Version:   $Revision: 1.23 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkCaminoStreamlineFileReader_hxx
#define __itkCaminoStreamlineFileReader_hxx

#include "itkCaminoStreamlineFileReader.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLabelContourImageFilter.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkByteSwapper.h"

#include <fstream>
#include <stdio.h>
#include <string>

namespace itk {

//
// Constructor
//
template<class TOutputMesh>
CaminoStreamlineFileReader<TOutputMesh>
::CaminoStreamlineFileReader()
{
  //
  // Create the output
  //
  typename TOutputMesh::Pointer output = TOutputMesh::New();
  this->ProcessObject::SetNumberOfRequiredOutputs( 1 );
  this->ProcessObject::SetNthOutput( 0, output.GetPointer() );
}


template<class TOutputMesh>
void
CaminoStreamlineFileReader<TOutputMesh>
::GenerateData()
{
  if( this->m_FileName == "" )
    {
    itkExceptionMacro( "No input FileName" );
    return;
    }

  this->m_InputFile.open( this->m_FileName.c_str(), std::ifstream::in );
  if (!this->m_InputFile )
    {
    std::cerr << "Unable to open file: " << this->m_FileName << std::endl;
    return;
    }

  unsigned long nPoints=0;
  unsigned long nTracts=0;
  float n;
  while (!this->m_InputFile.eof()) {
    float nTractPoints[1];
    this->m_InputFile.read( reinterpret_cast< char * >( nTractPoints ), sizeof(n) );
    float skip[(int)nTractPoints[0]+1];
    this->m_InputFile.read( reinterpret_cast< char * >( skip ), (3*nTractPoints[0]+1)*sizeof(n) );
    nPoints += nTractPoints[0];
    nTracts++;
  }

  Rcpp::Rcout << "Read " << nTracts << " streamlines with " << nPoints << " points" << std::endl;

}


template<class TOutputMesh>
void
CaminoStreamlineFileReader<TOutputMesh>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "FileName: "
     << this->m_FileName << std::endl;
}

} //end of namespace itk


#endif
