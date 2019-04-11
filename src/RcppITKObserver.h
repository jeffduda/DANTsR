#include "itkCommand.h"
#include <RcppDANTsR.h>

class RcppITKObserver : public itk::Command
{
  public:
    itkNewMacro( RcppITKObserver );

    void Execute(itk::Object * caller, const itk::EventObject & event) override
      {
      Execute( (const itk::Object *)caller, event);
      }

    void Execute(const itk::Object * caller, const itk::EventObject & event) override
      {
      if( ! itk::ProgressEvent().CheckEvent( &event ) )
        {
        return;
        }
      const auto * processObject =
        dynamic_cast< const itk::ProcessObject * >( caller );
      if( ! processObject )
        {
        return;
        }
      Rcpp::Rcout << "Progress: " << processObject->GetProgress() << std::endl;
      }
};
