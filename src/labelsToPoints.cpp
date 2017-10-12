
#include <algorithm>
#include <vector>
#include <string>
#include <RcppDANTsR.h>

#include "itkImage.h"
#include "itkLabeledImageToPointSetFilter.h"
#include "itkMesh.h"

template< class ImageType, class MeshType >
SEXP labelsToPoints( SEXP r_img, SEXP r_label )
{

  typedef typename ImageType::PixelType PixelType;
  typedef typename ImageType::Pointer ImagePointerType;
  typedef typename MeshType::Pointer MeshPointerType;

  ImagePointerType image = Rcpp::as<ImagePointerType>(r_img);

  PixelType label = Rcpp::as< PixelType >( r_label );

  typedef itk::LabeledImageToPointSetFilter<ImageType, MeshType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(image);
  filter->SetLabel(label);
  filter->Update();

  typename MeshType::Pointer outMesh = filter->GetOutput();

  return( Rcpp::wrap<MeshPointerType>(outMesh) );

  /*
  Rcpp::NumericMatrix pts(filter->GetOutput()->GetNumberOfPoints(), ImageType::ImageDimension );
  for ( unsigned int i=0; i<filter->GetOutput()->GetNumberOfPoints(); i++)
  {
    for ( unsigned int j=0; j<ImageType::ImageDimension; j++)
    {
      pts(i,j) = filter->GetOutput()->GetPoint(i)[j];
    }
  }

  return( Rcpp::wrap(pts) );
  */
}



RcppExport SEXP labelsToPoints( SEXP r_img, SEXP r_label ) {
try
  {
  if( r_img == NULL || r_label == NULL )
    {
      Rcpp::Rcout << "Unspecified Arguments" << std::endl;
      return Rcpp::wrap( NA_REAL );
    }

  Rcpp::S4 antsimage( r_img ) ;
  std::string pixeltype = Rcpp::as< std::string >( antsimage.slot( "pixeltype" ) );
  unsigned int dimension = Rcpp::as< int >( antsimage.slot( "dimension" ) );
  unsigned int components = Rcpp::as< int >( antsimage.slot( "components" ) );
  //bool isVector = Rcpp::as<bool>( antsimage.slot("isVector") );

  if ( components != 1 ) {
    Rcpp::Rcout << "Input must have 1 components" << std::endl;
    return Rcpp::wrap(NA_REAL);
  }

  if ( ( dimension > 4) || (dimension < 2) ) {
    Rcpp::Rcout << "Invalid Dimension: must be 2,3, or 4" << std::endl;
    return Rcpp::wrap(NA_REAL);
  }

  if( pixeltype == "double" )
    {
    typedef double PixelType;
    if ( dimension == 2 )
      {
      const unsigned int Dimension = 2;
      typedef itk::Image< PixelType, Dimension >  ImageType;
      typedef itk::Mesh< PixelType, Dimension >   MeshType;
      return  labelsToPoints< ImageType, MeshType >( r_img, r_label );
      }
    else if ( dimension == 3 )
      {
      const unsigned int Dimension = 3;
      typedef itk::Image< PixelType, Dimension >  ImageType;
      typedef itk::Mesh< PixelType, Dimension >   MeshType;
      return  labelsToPoints< ImageType, MeshType >( r_img, r_label );
      }
    else if ( dimension == 4 )
      {
      const unsigned int Dimension = 4;
      typedef itk::Image< PixelType, Dimension >  ImageType;
      typedef itk::Mesh< PixelType, Dimension >   MeshType;
      return  labelsToPoints< ImageType, MeshType >( r_img, r_label );
      }
    }
  else if( pixeltype == "float" )
    {
    typedef float PixelType;
    if ( dimension == 2 )
      {
      const unsigned int Dimension = 2;
      typedef itk::Image< PixelType, Dimension >  ImageType;
      typedef itk::Mesh< PixelType, Dimension >   MeshType;
      return  labelsToPoints< ImageType, MeshType >( r_img, r_label );
      }
    else if ( dimension == 3 )
      {
      const unsigned int Dimension = 3;
      typedef itk::Image< PixelType, Dimension >  ImageType;
      typedef itk::Mesh< PixelType, Dimension >   MeshType;
      return  labelsToPoints< ImageType, MeshType >( r_img, r_label );
      }
    else if ( dimension == 4 )
      {
      const unsigned int Dimension = 4;
      typedef itk::Image< PixelType, Dimension >  ImageType;
      typedef itk::Mesh< PixelType, Dimension >   MeshType;
      return  labelsToPoints< ImageType, MeshType >( r_img, r_label );
      }
    }
  else if( pixeltype == "unsigned int" )
    {
    typedef unsigned int PixelType;
    if ( dimension == 2 )
      {
      const unsigned int Dimension = 2;
      typedef itk::Image< PixelType, Dimension >  ImageType;
      typedef itk::Mesh< float, Dimension >   MeshType;
      return  labelsToPoints< ImageType, MeshType >( r_img, r_label );
      }
    else if ( dimension == 3 )
      {
      const unsigned int Dimension = 3;
      typedef itk::Image< PixelType, Dimension >  ImageType;
      typedef itk::Mesh< float, Dimension >   MeshType;
      return  labelsToPoints< ImageType, MeshType >( r_img, r_label );
      }
    else if ( dimension == 4 )
      {
      const unsigned int Dimension = 4;
      typedef itk::Image< PixelType, Dimension >  ImageType;
      typedef itk::Mesh< float, Dimension >   MeshType;
      return  labelsToPoints< ImageType, MeshType >( r_img, r_label );
      }
    }
  else if( pixeltype == "unsigned char" )
    {
    typedef unsigned char PixelType;
    if ( dimension == 2 )
      {
      const unsigned int Dimension = 2;
      typedef itk::Image< PixelType, Dimension >  ImageType;
      typedef itk::Mesh< float, Dimension >   MeshType;
      return  labelsToPoints< ImageType, MeshType >( r_img, r_label );
      }
    else if ( dimension == 3 )
      {
      const unsigned int Dimension = 3;
      typedef itk::Image< PixelType, Dimension >  ImageType;
      typedef itk::Mesh< float, Dimension >   MeshType;
      return  labelsToPoints< ImageType, MeshType >( r_img, r_label );
      }
    else if ( dimension == 4 )
      {
      const unsigned int Dimension = 4;
      typedef itk::Image< PixelType, Dimension >  ImageType;
      typedef itk::Mesh< float, Dimension >   MeshType;
      return  labelsToPoints< ImageType, MeshType >( r_img, r_label );
      }
    }
  else
    {
    Rcpp::Rcout << "Unsupported PixelType" << std::endl ;
    return Rcpp::wrap( NA_REAL ) ;
    }

  }
catch( const std::exception& exc )
  {
    Rcpp::Rcout<< exc.what() << std::endl;
    return Rcpp::wrap( NA_REAL );
  }

// Never reached
return Rcpp::wrap( NA_REAL );

}
