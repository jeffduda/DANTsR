#include <RcppANTsR.h>
#include "itkImage.h"
#include "itkAffineTransform.h"
#include "itkTranslationTransform.h"
#include "itkEuler2DTransform.h"
#include "itkEuler3DTransform.h"
#include "itkCenteredAffineTransform.h"
/*
 * itkImageRegistrationMethod_UpdateFilterRun
 * sets up the registration, runs it, and returns results
 */
template<class RegistrationType>
SEXP antsrRegistrationUpdateFilterRun( SEXP r_registrationObject );

/*
 * itkImageRegistrationMethod_UpdateFilter
 * handles the templating over different registration methods
 */
template< class ImageType >
SEXP antsrRegistrationUpdateFilter( SEXP r_registrationObject );

/*
 * itkImageRegistrationMethod_UpdateFilter
 * handles the interfacing with R
 */
RcppExport SEXP antsrRegistrationUpdateFilter( SEXP r_registrationObject );


/**
 * Transform traits to generalize the rigid transform
 */
template <class TComputeType, unsigned int ImageDimension>
class RigidTransformTraits
{
  // Don't worry about the fact that the default option is the
  // affine Transform, that one will not actually be instantiated.
public:
typedef itk::AffineTransform<TComputeType, ImageDimension> TransformType;
};

template <>
class RigidTransformTraits<double, 2>
{
public:
typedef itk::Euler2DTransform<double> TransformType;
};

template <>
class RigidTransformTraits<float, 2>
{
public:
typedef itk::Euler2DTransform<float> TransformType;
};

template <>
class RigidTransformTraits<double, 3>
{
public:
  // typedef itk::VersorRigid3DTransform<double>    TransformType;
  // typedef itk::QuaternionRigidTransform<double>  TransformType;
typedef itk::Euler3DTransform<double> TransformType;
};

template <>
class RigidTransformTraits<float, 3>
{
public:
  // typedef itk::VersorRigid3DTransform<float>    TransformType;
  // typedef itk::QuaternionRigidTransform<float>  TransformType;
typedef itk::Euler3DTransform<float> TransformType;
};
