#ifndef __RCPPANTSR_MESH_HPP
#define __RCPPANTSR_MESH_HPP

#include "itkMacro.h"
#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkVector.h"
#include "itkMesh.h"
#include "itkImageToImageMetricv4.h"
#include <Rcpp.h>

namespace Rcpp {

/*
template <> inline
SEXP wrap( const itk::ImageRegionIteratorWithIndex< itk::MeshType<PIXELTYPE,DIMENSION> > & iterator )
{
  typedef itk::MeshType<PIXELTYPE,DIMENSION>          MeshType;
  typedef MeshType::Pointer                           ImagePointerType;
  typedef itk::ImageRegionIteratorWithIndex<MeshType> IteratorType;

  IteratorType* rawPointer = new IteratorType( iterator );
  Rcpp::XPtr< IteratorType > xptr( rawPointer, true ) ;

  Rcpp::S4 itkImageIterator( std::string( "antsImageIterator" ) );
  itkImageIterator.slot( "pixeltype" ) = "PIXELTYPE";
  itkImageIterator.slot( "dimension" ) = DIMENSION;
  itkImageIterator.slot( "components" ) = 1;
  itkImageIterator.slot( "pointer" ) = xptr;

  return(wrap(itkImageIterator));
}
*/


template <> inline
SEXP wrap( const itk::Mesh< double, 2 >::Pointer & mesh )
{
  typedef itk::Mesh<double,2>             MeshType;
  typedef MeshType::Pointer               MeshPointerType;
  typedef Rcpp::XPtr<MeshPointerType>     XPtrType;

  MeshPointerType* rawPointer = new MeshPointerType( mesh );
  XPtrType xptr( rawPointer , true );

  Rcpp::S4 rMesh( std::string( "antsrMesh" ) );
  rMesh.slot( "precision" ) = "double";
  rMesh.slot( "dimension" ) = 2;
  rMesh.slot( "pointer" ) = xptr;

  return(wrap(rMesh));
}

template <> inline
SEXP wrap( const itk::Mesh< double, 3 >::Pointer & mesh )
{
  typedef itk::Mesh<double,3>             MeshType;
  typedef MeshType::Pointer               MeshPointerType;
  typedef Rcpp::XPtr<MeshPointerType>   XPtrType;

  MeshPointerType* rawPointer = new MeshPointerType( mesh );
  XPtrType xptr( rawPointer , true );

  Rcpp::S4 rMesh( std::string( "antsrMesh" ) );
  rMesh.slot( "precision" ) = "double";
  rMesh.slot( "dimension" ) = 3;
  rMesh.slot( "pointer" ) = xptr;

  return(wrap(rMesh));
}

template <> inline
SEXP wrap( const itk::Mesh< double, 4 >::Pointer & mesh )
{
  typedef itk::Mesh<double,4>             MeshType;
  typedef MeshType::Pointer               MeshPointerType;
  typedef Rcpp::XPtr<MeshPointerType>   XPtrType;

  MeshPointerType* rawPointer = new MeshPointerType( mesh );
  XPtrType xptr( rawPointer , true );

  Rcpp::S4 rMesh( std::string( "antsrMesh" ) );
  rMesh.slot( "precision" ) = "double";
  rMesh.slot( "dimension" ) = 4;
  rMesh.slot( "pointer" ) = xptr;

  return(wrap(rMesh));
}

template <> inline
SEXP wrap( const itk::Mesh< float, 2 >::Pointer & mesh )
{
  typedef itk::Mesh<float,2>             MeshType;
  typedef MeshType::Pointer               MeshPointerType;
  typedef Rcpp::XPtr<MeshPointerType>   XPtrType;

  MeshPointerType* rawPointer = new MeshPointerType( mesh );
  XPtrType xptr( rawPointer , true );

  Rcpp::S4 rMesh( std::string( "antsrMesh" ) );
  rMesh.slot( "precision" ) = "float";
  rMesh.slot( "dimension" ) = 2;
  rMesh.slot( "pointer" ) = xptr;

  return(wrap(rMesh));
}

template <> inline
SEXP wrap( const itk::Mesh< float, 3 >::Pointer & mesh )
{
  typedef itk::Mesh<float,3>             MeshType;
  typedef MeshType::Pointer               MeshPointerType;
  typedef Rcpp::XPtr<MeshPointerType>   XPtrType;

  MeshPointerType* rawPointer = new MeshPointerType( mesh );
  XPtrType xptr( rawPointer , true );

  Rcpp::S4 rMesh( std::string( "antsrMesh" ) );
  rMesh.slot( "precision" ) = "float";
  rMesh.slot( "dimension" ) = 3;
  rMesh.slot( "pointer" ) = xptr;

  return(wrap(rMesh));
}

template <> inline
SEXP wrap( const itk::Mesh< float, 4 >::Pointer & mesh )
{
  typedef itk::Mesh<float,4>             MeshType;
  typedef MeshType::Pointer               MeshPointerType;
  typedef Rcpp::XPtr<MeshPointerType>   XPtrType;

  MeshPointerType* rawPointer = new MeshPointerType( mesh );
  XPtrType xptr( rawPointer , true );

  Rcpp::S4 rMesh( std::string( "antsrMesh" ) );
  rMesh.slot( "precision" ) = "float";
  rMesh.slot( "dimension" ) = 4;
  rMesh.slot( "pointer" ) = xptr;

  return(wrap(rMesh));
}



template <> inline
itk::Mesh< double,2 >::Pointer as( SEXP mesh )
{
  const unsigned int Dim = 2;
  typedef itk::Mesh<double,Dim>  MeshType;
  typedef MeshType::Pointer      MeshPointerType;

  Rcpp::S4 rMesh( mesh );

  if (!rMesh.is( "antsrMesh") ||
      (Rcpp::as<std::string>(rMesh.slot("precision")) != "double") ||
      (Rcpp::as<int>(rMesh.slot("dimension")) != Dim) )
    {
    Rcpp::stop( "Invalid S4 object type");
    }

  XPtr<MeshPointerType> xptr( static_cast<SEXP>( rMesh.slot("pointer") ));
  return *xptr;
}

template <> inline
itk::Mesh< double,3 >::Pointer as( SEXP mesh )
{
  const unsigned int Dim = 3;
  typedef itk::Mesh<double,Dim>  MeshType;
  typedef MeshType::Pointer      MeshPointerType;

  Rcpp::S4 rMesh( mesh );

  if (!rMesh.is( "antsrMesh") ||
      (Rcpp::as<std::string>(rMesh.slot("precision")) != "double") ||
      (Rcpp::as<int>(rMesh.slot("dimension")) != Dim) )
    {
    Rcpp::stop( "Invalid S4 object type");
    }

  XPtr<MeshPointerType> xptr( static_cast<SEXP>( rMesh.slot("pointer") ));
  return *xptr;
}

template <> inline
itk::Mesh< double,4 >::Pointer as( SEXP mesh )
{
  const unsigned int Dim = 4;
  typedef itk::Mesh<double,Dim>  MeshType;
  typedef MeshType::Pointer      MeshPointerType;

  Rcpp::S4 rMesh( mesh );

  if (!rMesh.is( "antsrMesh") ||
      (Rcpp::as<std::string>(rMesh.slot("precision")) != "double") ||
      (Rcpp::as<int>(rMesh.slot("dimension")) != Dim) )
    {
    Rcpp::stop( "Invalid S4 object type");
    }

  XPtr<MeshPointerType> xptr( static_cast<SEXP>( rMesh.slot("pointer") ));
  return *xptr;
}

template <> inline
itk::Mesh< float,2 >::Pointer as( SEXP mesh )
{
  const unsigned int Dim = 2;
  typedef itk::Mesh<float,Dim>   MeshType;
  typedef MeshType::Pointer      MeshPointerType;

  Rcpp::S4 rMesh( mesh );

  if (!rMesh.is( "antsrMesh") ||
      (Rcpp::as<std::string>(rMesh.slot("precision")) != "float") ||
      (Rcpp::as<int>(rMesh.slot("dimension")) != Dim) )
    {
    Rcpp::stop( "Invalid S4 object type");
    }

  XPtr<MeshPointerType> xptr( static_cast<SEXP>( rMesh.slot("pointer") ));
  return *xptr;
}

template <> inline
itk::Mesh< float,3 >::Pointer as( SEXP mesh )
{
  const unsigned int Dim = 3;
  typedef itk::Mesh<float,Dim>   MeshType;
  typedef MeshType::Pointer      MeshPointerType;

  Rcpp::S4 rMesh( mesh );

  if (!rMesh.is( "antsrMesh") ||
      (Rcpp::as<std::string>(rMesh.slot("precision")) != "float") ||
      (Rcpp::as<int>(rMesh.slot("dimension")) != Dim) )
    {
    Rcpp::stop( "Invalid S4 object type");
    }

  XPtr<MeshPointerType> xptr( static_cast<SEXP>( rMesh.slot("pointer") ));
  return *xptr;
}

template <> inline
itk::Mesh< float,4 >::Pointer as( SEXP mesh )
{
  const unsigned int Dim = 4;
  typedef itk::Mesh<float,Dim>   MeshType;
  typedef MeshType::Pointer      MeshPointerType;

  Rcpp::S4 rMesh( mesh );

  if (!rMesh.is( "antsrMesh") ||
      (Rcpp::as<std::string>(rMesh.slot("precision")) != "float") ||
      (Rcpp::as<int>(rMesh.slot("dimension")) != Dim) )
    {
    Rcpp::stop( "Invalid S4 object type");
    }

  XPtr<MeshPointerType> xptr( static_cast<SEXP>( rMesh.slot("pointer") ));
  return *xptr;
}


}

#endif
