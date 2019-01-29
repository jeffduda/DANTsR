#include <RcppDANTsR.h>
#include <vnl/vnl_random.h>
#include <vnl/vnl_trace.h>
#include <vnl/algo/vnl_ldl_cholesky.h>
#include <vnl/algo/vnl_qr.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/algo/vnl_generalized_eigensystem.h>
#include "itkImage.h"

/* sparseDecomOne
 * matrix: a matrix where rows are samples and columns are features
 * n_vecs: the number of components to return
 * sparseness: the fraction of each component with non-zero values
 */

using RealType = double;
using MatrixType = vnl_matrix<RealType>;
using VectorType = vnl_vector<RealType>;
using VariateType = MatrixType;
using DiagonalMatrixType = vnl_diag_matrix<RealType>;
using ImageType = itk::Image<double, 3>;
using ImagePointerType = typename ImageType::Pointer;


RcppExport SEXP sparseDecomOne( SEXP r_matrix, SEXP r_nvecs, SEXP r_sparseness );

SEXP sparseDecomOneFunction( MatrixType matrix,  unsigned int n_vecs, double sparseness );

VectorType InitializeVector( MatrixType p, unsigned long seed );

VectorType OrthogonalizeVector(VectorType Mvec, VectorType V, MatrixType* projecterM=ITK_NULLPTR, MatrixType* projecterV=ITK_NULLPTR );

VectorType SpatiallySmoothVector( VectorType, ImagePointerType, double smoothing, bool surface=true );
