
#include <algorithm>
#include <vector>
#include <string>
#include <RcppANTsR.h>

#include "itkImage.h"
#include "itkLabeledImageToPointSetFilter.h"
#include "itkMesh.h"
#include "itkPointSet.h"

template< class ImageType, class MeshType >
SEXP deterministicTracking( SEXP r_dfield, SEXP r_seeds )
{

  typedef typename ImageType::PixelType PixelType;
  typedef typename ImageType::Pointer ImagePointerType;

  ImagePointerType dfield = Rcpp::as<ImagePointerType>(r_dfield);

  //PixelType label = Rcpp::as< PixelType >( r_label );
  // convert point matrix to itk::Mesh

  Rcpp::NumericMatrix pts = Rcpp::as<Rcpp::NumericMatrix>( r_seeds );
  typename MeshType::Pointer mesh = MeshType::New();
  mesh->GetPoints()->Reserve( pts.nrow() );

  typename MeshType::PointType pt;

  for ( itk::IdentifierType i=0; i<pts.nrow(); i++)
  {
    for ( itk::IdentifierType j=0; j<ImageType::ImageDimension; j++)
    {
      pt[j] = mesh->GetPoints()->GetElement(i)[j];
    }
    mesh->GetPoints()->SetElement(i, pt);
  }



  return( Rcpp::wrap(1) );

}



RcppExport SEXP deterministicTracking( SEXP r_dfield, SEXP r_seeds ) {
try
  {
  if( r_dfield == NULL || r_seeds == NULL )
    {
      Rcpp::Rcout << "Unspecified Arguments" << std::endl;
      return Rcpp::wrap( NA_REAL );
    }

  Rcpp::S4 antsimage( r_dfield ) ;
  std::string pixeltype = Rcpp::as< std::string >( antsimage.slot( "pixeltype" ) );
  unsigned int dimension = Rcpp::as< int >( antsimage.slot( "dimension" ) );
  unsigned int components = Rcpp::as< int >( antsimage.slot( "components" ) );
  bool isVector = Rcpp::as<bool>( antsimage.slot("isVector") );

  if ( components != dimension ) {
    Rcpp::Rcout << "Input must have number of components equal to image dimension" << std::endl;
    return Rcpp::wrap(NA_REAL);
  }

  if ( ( dimension > 3) || (dimension < 2) ) {
    Rcpp::Rcout << "Invalid Dimension: must be 2,3,or 4" << std::endl;
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
      return  deterministicTracking< ImageType, MeshType >( r_dfield, r_seeds );
      }
    else if ( dimension == 3 )
      {
      const unsigned int Dimension = 3;
      typedef itk::Image< PixelType, Dimension >  ImageType;
      typedef itk::Mesh< PixelType, Dimension >   MeshType;
      return  deterministicTracking< ImageType, MeshType >( r_dfield, r_seeds );
      }
    else if ( dimension == 4 )
      {
      const unsigned int Dimension = 4;
      typedef itk::Image< PixelType, Dimension >  ImageType;
      typedef itk::Mesh< PixelType, Dimension >   MeshType;
      return  deterministicTracking< ImageType, MeshType >( r_dfield, r_seeds );
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
      return  deterministicTracking< ImageType, MeshType >( r_dfield, r_seeds );
      }
    else if ( dimension == 3 )
      {
      const unsigned int Dimension = 3;
      typedef itk::Image< PixelType, Dimension >  ImageType;
      typedef itk::Mesh< PixelType, Dimension >   MeshType;
      return  deterministicTracking< ImageType, MeshType >( r_dfield, r_seeds );
      }
    else if ( dimension == 4 )
      {
      const unsigned int Dimension = 4;
      typedef itk::Image< PixelType, Dimension >  ImageType;
      typedef itk::Mesh< PixelType, Dimension >   MeshType;
      return  deterministicTracking< ImageType, MeshType >( r_dfield, r_seeds );
      }
    }
  else
    {
    Rcpp::Rcout << "Unsupported PixelType - must be float or double" << std::endl ;
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
