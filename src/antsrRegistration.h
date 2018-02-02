#include <RcppDANTsR.h>
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
SEXP antsrRegistrationMethodRun( SEXP r_registrationObject, SEXP r_fixed, SEXP r_moving );

template< class RegistrationType >
SEXP antsrRegistrationStage(SEXP r_registrationObject, SEXP r_fixed, SEXP r_moving );

template<class RegistrationType>
SEXP antsrRegistrationWithHelper(SEXP r_registrationObject, SEXP r_fixed, SEXP r_moving );

/*
 * itkImageRegistrationMethod_UpdateFilter
 * handles the templating over different registration methods
 */
template< class ImageType >
SEXP antsrRegistrationMethod( SEXP r_registrationObject, SEXP r_fixed, SEXP r_moving );

/*
 * itkImageRegistrationMethod_UpdateFilter
 * handles the interfacing with R
 */
RcppExport SEXP antsrRegistrationRun( SEXP r_registrationObject, SEXP r_fixed, SEXP r_moving );


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
