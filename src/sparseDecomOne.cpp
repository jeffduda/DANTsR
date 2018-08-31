/*=========================================================================


  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
  https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt
  for details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "itkMinimumMaximumImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkExtractImageFilter.h"
#include <vnl/vnl_random.h>
#include <vnl/vnl_trace.h>
#include <vnl/algo/vnl_ldl_cholesky.h>
#include <vnl/algo/vnl_qr.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/algo/vnl_generalized_eigensystem.h>
#include "antsSCCANObject.h"
#include <time.h>
#include "itkCSVNumericObjectFileWriter.h"
#include "itkCSVArray2DDataObject.h"
#include "itkCSVArray2DFileReader.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkSurfaceImageCurvature.h"
#include "itkImageFileWriter.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "ReadWriteData.h"
#include "sparseDecomOne.h"

/* sparseDecomOne
 * matrix: a matrix where rows are samples and columns are features
 * n_vecs: the number of components to return
 * sparseness: the fraction of each component with non-zero values
 */

RcppExport SEXP sparseDecomOne( SEXP r_matrix, SEXP r_nvecs, SEXP r_sparseness ) {

  unsigned int n_vecs = Rcpp::as<unsigned int>( r_nvecs );
  double sparseness = Rcpp::as<double>(r_sparseness);

  Rcpp::NumericMatrix X( r_matrix );
  std::vector<double> xdat =  Rcpp::as< std::vector<double> >( X );
  const double* _xdata = &xdat[0];
  MatrixType matrix( _xdata , X.cols(), X.rows()  );
  matrix = matrix.transpose();

  return( sparseDecomOneFunction( matrix, n_vecs, sparseness) );
}

SEXP sparseDecomOneFunction( MatrixType matrix,  unsigned int n_vecs, double sparseness )
{

  bool prior = false;

  RealType reconerr = 0;
  unsigned long clusterSize = 50;
  bool useL1 = false;
  bool keepPositive = false;


  // vector of correlation values
  VectorType canonicalCorrelations( n_vecs, 0.0);

  // vector holding sparseness level for each component
  VectorType sparsenessparams( n_vecs, sparseness );

  // If using image based prior regions
  // set sparseness by sparseness of prior regions
  /*
  if ( ( this->m_OriginalMatrixPriorROI.rows() > 0  ) &&
       ( this->m_OriginalMatrixPriorROI.cols() > 0  ) )
    {
    prior = true;
    Rcpp::Rcout << "Initialization: image-driven" << std::endl;
    this->m_MatrixPriorROI = this->m_OriginalMatrixPriorROI;
    n_vecs = this->m_OriginalMatrixPriorROI.rows();
    for( unsigned int x = 0; x < this->m_OriginalMatrixPriorROI.rows(); x++ )
      {
      VectorType   priorrow = this->m_OriginalMatrixPriorROI.get_row( x );
      RealType fnz = 0;
      for( unsigned int y = 0; y < this->m_OriginalMatrixPriorROI.cols(); y++ )
      	{
      	if(vnl_math_abs( priorrow( y ) ) > 1.e-8 )
          {
      	  fnz += 1;
      	  }
      	}
      if ( vnl_math_abs( this->m_FractionNonZeroP ) < 1.e-11 )
        sparsenessparams( x ) = 1.0 * (RealType) fnz / (RealType) this->m_OriginalMatrixPriorROI.cols();
      priorrow = priorrow/priorrow.two_norm();
      this->m_MatrixPriorROI.set_row( x , priorrow );
      }
    this->m_FractionNonZeroP = sparsenessparams.mean();
    if ( ! this->m_Silent )  std::cout <<"sparseness: " << sparsenessparams << std::endl;
    this->m_VariatesP = this->m_MatrixPriorROI.transpose();
    }
  */


  MatrixType matrixB( matrix.rows(), n_vecs );
  matrixB.fill( 0 );

  // Make this prerec - call scale(matrix) first in R
  //MatrixType matrixP = NormalizeMatrix<MatrixType>( matrix );

  VectorType icept( matrix.rows(), 0 );

  MatrixType variatesP;

  // Auto-select sparseness parameters
  /*
  if ( sparseness < 1.e-11 )
    { // estimate sparseness from PCA
      if ( ! this->m_Silent )  std::cout << " data-driven initialization from PCA " << std::endl;
      VectorType maxvals( this->m_MatrixP.cols() , 0 );
      MatrixType        cov = this->m_MatrixP * this->m_MatrixP.transpose();
      vnl_svd<RealType> eig( cov, 1.e-6 );
      this->m_VariatesP.set_size( this->m_MatrixP.cols(), n_vecs );
      this->m_VariatesP.fill( 0 );
      for( unsigned int i = 0; i < n_vecs; i++ )
	      {
    	  if( i < this->m_MatrixP.rows() )
    	    {
    	    VectorType u = eig.U().get_column( i );
    	    VectorType up = u * this->m_MatrixP;
    	    up = up / up.two_norm();
    	    this->m_VariatesP.set_column( i, up );
    	    }
	      }
      for( unsigned int i = 0; i < maxvals.size(); i++ )
	      {
        RealType maxval = 0;
    	  for( unsigned int j = 0; j < n_vecs; j++ )
    	    {
    	    RealType myval = vnl_math_abs( this->m_VariatesP( i, j ) );
    	    if (  myval > maxval )
    	      {
    	      maxvals( i ) = j;
    	      maxval = myval;
    	      }
    	    }
	      }
  sparsenessparams.fill( 0 );
  for( unsigned int i = 0; i < maxvals.size(); i++ )
  	{
  	unsigned int index = (unsigned int) ( maxvals( i ) + 0.5 );
  	sparsenessparams( index ) = sparsenessparams( index ) + 1;
  	}
  for( unsigned int i = 0; i < n_vecs; i++ )
  	{
  	sparsenessparams( i ) = 1.01 * ( sparsenessparams( i ) / this->m_MatrixP.cols() );
  	this->m_FractionNonZeroP = sparsenessparams( i );
  	VectorType vec = this->m_VariatesP.get_column( i );
  	this->SparsifyP( vec );
  	this->m_VariatesP.set_column( i, vec );
	   if ( ! this->m_Silent )  std::cout <<" Estimated-Sparseness " << i << " is " << sparsenessparams( i ) << std::endl;
	  }
    if ( ! this->m_Silent )  std::cout <<" Estimated-Sparseness Sum "  << sparsenessparams.sum() << std::endl;
    } // done estimate sparseness
  */



  // If no priors have been passed
  ImagePointerType mask = ITK_NULLPTR;

  if ( variatesP.rows() < 1  ) {
    Rcpp::Rcout << "Initialization: data-driven" << std::endl;

    variatesP.set_size( matrix.cols(), n_vecs );
    variatesP.fill( 0 );

    for( unsigned int i = 0; i < n_vecs; i++ )
      {
      VectorType initvec = InitializeVector( matrix, i + 1 ); // increment seed so each vec isn't the same?
      for( unsigned int j = 0; j < i; j++ )
        {
        initvec = OrthogonalizeVector( initvec, variatesP.get_column( j ) );
        }
      initvec = SpatiallySmoothVector( initvec, mask );

      // FIXME - why no re-center here?
      initvec = initvec / initvec.two_norm();
      variatesP.set_column( i, initvec );
      }
  }

  /** edit from here


  // now initialize B
  reconerr = this->SparseReconB( matrixB, icept  );
  this->SparseArnoldiSVD_Other( matrixB );
  if ( ! this->m_Silent )  std::cout << "begin : %var " << reconerr << std::endl;
  RealType matpfrobnorm = this->m_MatrixP.frobenius_norm();
  if ( false )
  for( unsigned int overit = 0; overit < this->m_MaximumNumberOfIterations; overit++ )
    {
    MatrixType vgrad = matrixB.transpose() * this->m_MatrixP -
      ( matrixB.transpose() * matrixB ) * this->m_VariatesP.transpose();
    this->m_VariatesP = this->m_VariatesP + vgrad.transpose();
    for(  unsigned int a = 0; a < n_vecs; a++ )
      {
      VectorType evec = this->m_VariatesP.get_column( a );
      this->SparsifyP( evec );
      this->m_VariatesP.set_column( a , evec );
      }
    reconerr = this->SparseReconB( matrixB, icept  );
    this->SparseArnoldiSVD_Other( matrixB );
    this->m_MatrixU = matrixB;
    if ( ! this->m_Silent )  std::cout << overit << ": %var " << reconerr << std::endl;
    }
  if ( true )
  for( unsigned int overit = 0; overit < this->m_MaximumNumberOfIterations; overit++ )
    {
       //a power iteration  method --- depends on the following
       //given any nonzero $z \in \mathbb{R}^n$, the Rayleigh quotient
       //$x^T X x / x^T x $ minimizes the function $\| \lambda x - X x \|^2 $
       //wrt $\lambda$.
       //so, if we find the vector x ( by sparse power iteration ) then we have a vector
       //that is a close approximation to the first eigenvector of X. If X is a residual
       //matrix then x is a good approximation of the $n^th$ eigenvector.

    VectorType   zero( this->m_MatrixP.cols(), 0 );
    VectorType   zerob( this->m_MatrixP.rows(), 0 );
    unsigned int a = 0;
    while(  a < n_vecs )
      {
      this->m_FractionNonZeroP = sparsenessparams( a );
      VectorType bvec = matrixB.get_column( a );
      matrixB.set_column( a, zerob );
      MatrixType tempMatrix = this->m_VariatesP;
      tempMatrix.set_column( a, zero );
      MatrixType partialmatrix = matrixB * tempMatrix.transpose();
      for(  unsigned int interc = 0; interc < this->m_MatrixP.rows(); interc++ )
	      {
        partialmatrix.set_row( interc,
          partialmatrix.get_row( interc ) + icept( interc ) );
 	      }
      partialmatrix = this->m_MatrixP - partialmatrix;
      this->m_CanonicalCorrelations[a] = 1 - ( partialmatrix.frobenius_norm() ) / matpfrobnorm;
      VectorType evec = this->m_VariatesP.get_column( a );
      this->m_VariatesP.set_column( a, zero );
      this->m_CanonicalCorrelations[a] = this->IHTPowerIteration(  partialmatrix,  evec, 5, a );    // 0 => a
      this->m_VariatesP.set_column( a, evec );
      matrixB.set_column( a, bvec );
      // update B matrix by linear regression
      reconerr = this->SparseReconB( matrixB, icept  );
      this->SparseArnoldiSVD_Other( matrixB );
      a++;
      } // while
    // update B matrix by linear regression
    reconerr = this->SparseReconB( matrixB, icept  );
    this->SparseArnoldiSVD_Other( matrixB );
    this->m_MatrixU = matrixB;
    if ( ! this->m_Silent )  std::cout << overit << ": %var " << reconerr << std::endl;
    if (  ( ! prior ) && ( overit == ( this->m_MaximumNumberOfIterations - 1 ) ) ) this->SortResults( n_vecs );
    }
  this->m_VariatesQ = matrixB;
  for( unsigned int i = 0; i < n_vecs; i++ )
    {
    VectorType v = this->m_VariatesP.get_column( i );
    if( v.min_value() < 0 )
      {
      this->m_VariatesP.set_column( i, v * ( -1 ) );
      }
    }


  return 1.0 / reconerr;
  */

  return( Rcpp::wrap(NA_REAL));



  /** a regression-based method
  for(  unsigned int a = 0; a < this->m_MatrixP.cols(); a++ )
    {
    VectorType x_i = this->m_MatrixP.get_column( a );
    VectorType lmsolv = this->m_VariatesP.get_row( a ); // good initialization should increase convergence speed
    (void) this->ConjGrad(  matrixB, lmsolv, x_i, 0, 10000 );    // A x = b
    VectorType x_recon = ( matrixB * lmsolv + this->m_Intercept );
    onenorm += x_i.one_norm() / this->m_MatrixP.rows();
    reconerr += ( x_i - x_recon ).one_norm() / this->m_MatrixP.rows();
    this->m_VariatesP.set_row( a , lmsolv );
    }
  */

  /* a gradient method
  E & = \langle X - U V^T , X - U V^T \rangle \\
  \partial E/\partial V & = \langle -U , X - U V^T \rangle \\
  & = -U^T X + U^T U V^T \rangle
  if ( overit == 0 ) this->m_VariatesP.fill( 0 );
  MatrixType vgrad = matrixB.transpose() * this->m_MatrixP -
    ( matrixB.transpose() * matrixB ) * this->m_VariatesP.transpose();
  //  this->m_VariatesP = this->m_VariatesP + vgrad.transpose() * 0.01;
  this->m_VariatesP = this->m_VariatesP + vgrad.transpose() ;
  for(  unsigned int a = 0; a < this->m_MatrixP.rows(); a++ )
    {
    VectorType evec = this->m_VariatesP.get_column( a );
    this->SparsifyP( evec );
    this->m_VariatesP.set_column( a, evec );
    }
  */


}

VectorType
InitializeVector( MatrixType p, unsigned long seed )
{
  VectorType w_p( p.columns() );

  w_p.fill(0);
  for( unsigned int its = 0; its < 1; its++ )
    {
    vnl_random randgen( seed ); /* FIXME - why this -> use constant seed to prevent weirdness */
    for( unsigned long i = 0; i < p.columns(); i++ )
      {
      if( seed > 0 )
        {
        w_p(i) = randgen.normal();
        // w_p(i)=randgen.drand32();
        }
      else
        {
        w_p(i) = 1.0;
        }
      }
    }
  w_p = w_p / p.columns();
  return w_p;
}

VectorType OrthogonalizeVector(VectorType Mvec, VectorType V, MatrixType* projecterM,  MatrixType* projecterV )
{
  if( ( !projecterM ) &&  ( !projecterV ) )
    {
    double ipv   = inner_product(V, V);
    if( ipv == 0 )
      {
      return Mvec;
      }
    double     ratio = inner_product(Mvec, V) / ipv;
    VectorType ortho = Mvec - V * ratio;
    return ortho;
    }
  else if( ( !projecterM ) && ( projecterV ) )
    {
    double     ratio = inner_product(Mvec, *projecterV * V) / inner_product(*projecterV * V, *projecterV * V);
    VectorType ortho = Mvec - V * ratio;
    return ortho;
    }
  else if( ( !projecterV ) && ( projecterM ) )
    {
    double     ratio = inner_product(*projecterM * Mvec, V) / inner_product(V, V);
    VectorType ortho = (*projecterM * Mvec) - V * ratio;
    return ortho;
    }
  else
    {
    double ratio = inner_product(*projecterM * Mvec, *projecterV * V) / inner_product(*projecterV * V,
                                                                                      *projecterV * V);
    VectorType ortho = Mvec - V * ratio;
    for( unsigned int i = 0; i < Mvec.size(); i++ )
      {
      if( Mvec(i) == 0 )
        {
        ortho(i) = 0;
        }
      }
    return ortho;
    }
}


VectorType
SpatiallySmoothVector( VectorType vec,  ImagePointerType mask, double smoothing, bool surface )
{
  if( mask.IsNull() || vnl_math_abs( smoothing ) < 1.e-9 )
    {
    return vec;
    }
  return vec;
}
/*

  RealType vecnorm = vec.two_norm();
  ImagePointer image = this->ConvertVariateToSpatialImage( vec, mask, false );
  RealType     spacingsize = 0;
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    RealType sp = mask->GetSpacing()[d];
    spacingsize += sp * sp;
    }
  spacingsize = sqrt( spacingsize );
  if ( this->m_Smoother  < 0.0  )
    {
    typedef itk::GradientAnisotropicDiffusionImageFilter< TInputImage,
      TInputImage > FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput( image );
    filter->SetNumberOfIterations( vnl_math_abs( this->m_Smoother ) );
    TRealType mytimestep = spacingsize / std::pow( 2.0 , static_cast<double>(ImageDimension+1) );
    TRealType reftimestep = 0.5 / std::pow( 2.0 , static_cast<double>(ImageDimension+1) );
    if ( mytimestep > reftimestep ) mytimestep = reftimestep;
    filter->SetTimeStep( mytimestep );
    filter->SetConductanceParameter( 0.25 ); // might need to change this
    filter->Update();
    VectorType gradvec = this->ConvertImageToVariate( filter->GetOutput(),  mask );
    return gradvec;
    }
  if ( ( surface )  && ( ImageDimension == 3 ) && ( false ) )
    {
  unsigned int sigma = ( unsigned int ) this->m_Smoother;
  if (  this->m_Smoother > 0.0001 )
    if ( sigma < 1 ) sigma = 1;
  typedef itk::SurfaceImageCurvature<TInputImage> ParamType;
  typename ParamType::Pointer Parameterizer = ParamType::New();
  Parameterizer->SetInputImage(mask);
  Parameterizer->SetFunctionImage(image);
  Parameterizer->SetNeighborhoodRadius( 1 );
  Parameterizer->SetSigma( 1.0 );
  Parameterizer->SetUseGeodesicNeighborhood(false);
  Parameterizer->SetUseLabel(false);
  Parameterizer->SetThreshold(0.5);
  Parameterizer->IntegrateFunctionOverSurface(true);
  for( unsigned int i = 0; i < sigma; i++ )
    {
    Parameterizer->IntegrateFunctionOverSurface(true);
    }
  VectorType svec = this->ConvertImageToVariate( Parameterizer->GetFunctionImage(),  mask );
  return svec;
    }
  else
    {
    typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType> dgf;
    typename dgf::Pointer filter = dgf::New();
    filter->SetUseImageSpacingOn();
    filter->SetVariance(  this->m_Smoother * spacingsize );
    filter->SetMaximumError( .01f );
    filter->SetInput( image );
    filter->Update();
    typename ImageType::Pointer imgout = filter->GetOutput();
    //    WriteImage<ImageType>( imgout , "ccaout_s.nii.gz" );
    VectorType gradvec = this->ConvertImageToVariate( imgout,  mask );
    gradvec = gradvec * vecnorm / gradvec.two_norm();
    //    if ( ! this->m_Silent )  std::cout << ImageDimension << " gvec " << gradvec[1] << " " << gradvec[100] << " " << gradvec[1000] << std::endl;
    return gradvec;
    }
}


/*
template <class TInputImage, class TRealType>
TRealType antsSCCANObject<TInputImage, TRealType>
::SparseReconB(   typename antsSCCANObject<TInputImage,
                                           TRealType>::MatrixType& matrixB,
                  typename antsSCCANObject<TInputImage,
                                           TRealType>::VectorType& icept )
{
  MatrixType temp( this->m_MatrixP );
  icept.fill( 0 );
  RealType meancorr = 0;
  matrixB.set_size( this->m_MatrixP.rows(), this->m_VariatesP.cols() );
  matrixB.fill( 0 );
  for(  unsigned int a = 0; a < this->m_MatrixP.rows(); a++ )
    {
    VectorType x_i = this->m_MatrixP.get_row( a );  // 1 by p
    VectorType lmsolv = matrixB.get_row( a );       // 1 by k ....   variatesp is  p by k ...
    // solve   ( p by k ) ( k by 1 ) = ( p by 1 )
    // which is    \| V U_j - X_j \|^2  where j is the j^th subject
    (void) this->ConjGrad(  this->m_VariatesP, lmsolv, x_i, 0, 100 ); // A x = b
    icept( a ) = this->m_Intercept;
    VectorType x_recon = ( this->m_VariatesP * lmsolv + icept( a ) );
    matrixB.set_row( a, lmsolv );
    RealType localcorr = this->PearsonCorr( x_recon, x_i  );
    temp.set_row( a, x_recon );
    meancorr += localcorr;
    }
  if( false )
    {
    std::vector<std::string> ColumnHeaders;
    for( unsigned int nv = 0; nv < this->m_MatrixP.cols(); nv++ )
      {
      std::stringstream ss;
      ss << (nv + 1);
      std::string stringout = ss.str();
      std::string colname = std::string("X") + stringout;
      ColumnHeaders.push_back( colname );
      }
    MatrixType recon = ( matrixB * this->m_VariatesP.transpose() );
    typedef itk::CSVNumericObjectFileWriter<double, 1, 1> WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( "Zrecon.csv" );
    writer->SetInput( &recon );
    writer->Write();
    }
  this->m_MatrixU = matrixB;
  RealType matpfrobnorm = this->m_MatrixP.frobenius_norm();
  RealType rr = ( temp - this->m_MatrixP ).frobenius_norm();
  return ( 1.0 - rr * rr / ( matpfrobnorm * matpfrobnorm ) );
}
*/
