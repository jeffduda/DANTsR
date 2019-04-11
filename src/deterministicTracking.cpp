
#include <algorithm>
#include <vector>
#include <string>
#include <RcppDANTsR.h>

#include "itkDeterministicDTITractography.h"
#include "itkImage.h"
#include "itkLabeledImageToPointSetFilter.h"
#include "itkMesh.h"
#include "itkPointSet.h"
#include "RcppITKObserver.h"

template< class VectorImageType, class MeshType >
SEXP deterministicTracking( SEXP r_dfield, SEXP r_seeds, SEXP r_mask )
{
  //using PixelType = typename VectorImageType::PixelType;
  using ImagePointerType = typename VectorImageType::Pointer;
  using MeshPointerType = typename MeshType::Pointer;
  using TrackerType = itk::DeterministicDTITractography< VectorImageType, MeshType >;
  using TrackerPointerType = typename TrackerType::Pointer;

  using MaskPointerType = typename TrackerType::MaskPointerType;
  using SeedPointerType = typename TrackerType::SeedContainerPointer;

  RcppITKObserver::Pointer myCommand = RcppITKObserver::New();

  ImagePointerType dfield = Rcpp::as<ImagePointerType>(r_dfield);
  MeshPointerType seedMesh = Rcpp::as<MeshPointerType>(r_seeds);
  MaskPointerType mask = Rcpp::as<MaskPointerType>(r_mask);
  TrackerPointerType tracker = TrackerType::New();
  tracker->SetInput( dfield );
  tracker->SetSeeds( seedMesh );
  tracker->SetMask( mask );
  tracker->SetMinimumNumberOfPoints(2);
  tracker->SetMaximumNumberOfPoints(2000);
  tracker->AddObserver(itk::ProgressEvent(), myCommand);
  tracker->Update();
  MeshPointerType outMesh = tracker->GetOutput();
  //Rcpp::Rcout << "Returned from tracker" << std::endl;

  SeedPointerType seedOffsets = tracker->GetSeedOffsets();

  Rcpp::NumericVector seedVector( seedOffsets->Size() );
  for (unsigned int i=0; i<seedOffsets->Size(); i++) {
    seedVector[i] = seedOffsets->GetElement(i);
  }


  //return( Rcpp::wrap(outMesh) );



  Rcpp::List list = Rcpp::List::create(Rcpp::Named("Mesh")=Rcpp::wrap(outMesh),
                                       Rcpp::Named("Seeds")=seedVector);
  return(Rcpp::wrap(list));


}


RcppExport SEXP deterministicTracking( SEXP r_dfield, SEXP r_seeds, SEXP r_mask ) {
try
  {
  if( r_dfield == nullptr || r_seeds == nullptr )
    {
      Rcpp::Rcout << "Unspecified Arguments" << std::endl;
      return Rcpp::wrap( NA_REAL );
    }

  Rcpp::S4 antsimage( r_dfield ) ;
  std::string pixeltype = Rcpp::as< std::string >( antsimage.slot( "pixeltype" ) );
  unsigned int dimension = Rcpp::as< int >( antsimage.slot( "dimension" ) );
  unsigned int components = Rcpp::as< int >( antsimage.slot( "components" ) );
  //bool isVector = Rcpp::as<bool>( antsimage.slot("isVector") );

  Rcpp::S4 antsrmesh( r_seeds );
  unsigned int meshDim = Rcpp::as< int >( antsrmesh.slot("dimension") );

  if ( meshDim != dimension ) {
    Rcpp::Rcout << "Image and mesh must have same dimensions" << std::endl;
    return( Rcpp::wrap(NA_REAL) );
  }

  if ( r_mask != NULL ) {
    Rcpp::S4 scalarimage( r_mask );
    std::string pixeltype2 = Rcpp::as< std::string >( scalarimage.slot( "pixeltype" ) );
    unsigned int dimension2 = Rcpp::as< int >( scalarimage.slot( "dimension" ) );
    //unsigned int components2 = Rcpp::as< int >( scalarimage.slot( "components" ) );
    bool isVector2 = Rcpp::as<bool>( scalarimage.slot("isVector") );

    if ( isVector2 || (pixeltype2 != pixeltype) || (dimension2 != dimension) ) {
      Rcpp::Rcout << "Mask must have same dimension and pixeltype as vector field" << std::endl;
      return Rcpp::wrap(NA_REAL);
    }
  }


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
      typedef itk::VectorImage< PixelType, Dimension >  VectorImageType;
      typedef itk::Mesh< PixelType, Dimension >   MeshType;
      //typedef itk::Image< PixelType, Dimension >  ImageType;
      return  deterministicTracking< VectorImageType, MeshType >( r_dfield, r_seeds, r_mask );
      }
    else if ( dimension == 3 )
      {
      const unsigned int Dimension = 3;
      typedef itk::VectorImage< PixelType, Dimension >  VectorImageType;
      typedef itk::Mesh< PixelType, Dimension >   MeshType;
      //typedef itk::Image< PixelType, Dimension >  ImageType;
      return  deterministicTracking< VectorImageType, MeshType >( r_dfield, r_seeds, r_mask );
      }
    else if ( dimension == 4 )
      {
      const unsigned int Dimension = 4;
      typedef itk::VectorImage< PixelType, Dimension >  VectorImageType;
      typedef itk::Mesh< PixelType, Dimension >   MeshType;
      //typedef itk::Image< PixelType, Dimension >  ImageType;
      return  deterministicTracking< VectorImageType, MeshType >( r_dfield, r_seeds, r_mask );
      }
    }
  else if( pixeltype == "float" )
    {
    typedef float PixelType;
    if ( dimension == 2 )
      {
      const unsigned int Dimension = 2;
      typedef itk::VectorImage< PixelType, Dimension >  VectorImageType;
      typedef itk::Mesh< PixelType, Dimension >   MeshType;
      //typedef itk::Image< PixelType, Dimension >  ImageType;
      return  deterministicTracking< VectorImageType, MeshType>( r_dfield, r_seeds, r_mask );
    }
    else if ( dimension == 3 )
      {
      const unsigned int Dimension = 3;
      typedef itk::VectorImage< PixelType, Dimension >  VectorImageType;
      typedef itk::Mesh< PixelType, Dimension >   MeshType;
      //typedef itk::Image< PixelType, Dimension >  ImageType;
      return  deterministicTracking< VectorImageType, MeshType >( r_dfield, r_seeds, r_mask );
      }
    else if ( dimension == 4 )
      {
      const unsigned int Dimension = 4;
      typedef itk::VectorImage< PixelType, Dimension >  VectorImageType;
      typedef itk::Mesh< PixelType, Dimension >   MeshType;
      //typedef itk::Image< PixelType, Dimension >  ImageType;
      return  deterministicTracking< VectorImageType, MeshType >( r_dfield, r_seeds, r_mask );
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
