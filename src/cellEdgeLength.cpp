#include <string>
#include "cellEdgeLength.h"

template< class MeshType >
SEXP cellEdgeLength( SEXP r_Mesh, SEXP r_index )
{
  using MeshPointerType = typename MeshType::Pointer;
  using FunctionType = itk::EdgeLengthCellFunction<MeshType, double>;
  using FunctionPointerType = typename FunctionType::Pointer;

  MeshPointerType mesh = Rcpp::as<MeshPointerType>(r_Mesh);
  if ( ! mesh.IsNotNull() )
    {
    Rcpp::stop("mesh not yet allocated");
    }

  //typename MeshType::CellIdentifier index = Rcpp::as<typename MeshType::CellIdentifier>(r_index)-1;
  Rcpp::NumericVector indices( r_index );
  Rcpp::NumericVector lengths( indices.size(), 0.0 );

  FunctionPointerType function = FunctionType::New();
  function->SetInputMesh( mesh );
  for (unsigned long i=0; i<indices.size(); i++ ) {
    lengths[i] = function->Evaluate( indices[i]-1 );
  }

  return( Rcpp::wrap(lengths) );

  //return Rcpp::wrap( function->Evaluate(index) );

}

RcppExport SEXP cellEdgeLength( SEXP r_Mesh, SEXP r_index )
{
try
{
  if( r_Mesh == nullptr )
    {
    Rcpp::Rcout << "Unspecified Argument" << std::endl ;
    return Rcpp::wrap( 1 ) ;
    }

  Rcpp::S4 Mesh( r_Mesh );
  unsigned int dimension = Rcpp::as< int >( Mesh.slot( "dimension" ) );
  std::string precision = Rcpp::as< std::string >( Mesh.slot( "precision" ) );

  if ( ( dimension < 2) || ( dimension > 3)  )
  {
    Rcpp::stop("Invalid mesh dimension");
  }

  if ( precision=="double") {
    using PrecisionType = double;
    if (dimension == 2)
      {
      using MeshType = itk::Mesh<PrecisionType,2>;
      return cellEdgeLength<MeshType>( r_Mesh, r_index );
      }
    else if (dimension == 3)
      {
        using MeshType = itk::Mesh<PrecisionType,3>;
        return cellEdgeLength<MeshType>( r_Mesh, r_index );
      }
    }
  else if ( precision=="float") {
    using PrecisionType = float;
    if (dimension == 2)
      {
      using MeshType = itk::Mesh<PrecisionType,2>;
      return cellEdgeLength<MeshType>( r_Mesh, r_index );
      }
    else if (dimension == 3)
      {
      using MeshType = itk::Mesh<PrecisionType,3>;
      return cellEdgeLength<MeshType>( r_Mesh, r_index );
      }
    }
  else
  {
    Rcpp::stop("Invalid precision");
  }

  // Never reached
  return(Rcpp::wrap(NA_REAL));

  }
catch( itk::ExceptionObject & err )
  {
  Rcpp::Rcout << "ITK ExceptionObject caught !" << std::endl;
  Rcpp::Rcout << err << std::endl;
  Rcpp::stop("ITK exception caught");
  }
catch( const std::exception& exc )
  {
  forward_exception_to_r( exc ) ;
  }
catch(...)
  {
	Rcpp::stop("c++ exception (unknown reason)");
  }
return Rcpp::wrap(NA_REAL); //not reached
}
