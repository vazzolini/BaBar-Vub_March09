// Module: EvtItgSimpsonIntegrator.cc
//
// Description:
//      Abstraction of a generic function for use in integration methods elsewhere
//      in this package. 
//      (Stolen from  EvtGenModels package)
// Modified:
//      Alexei Volk  March 6, 2008
//

#ifndef XSLITGSIMPSONINTEGRATOR_HH
#define XSLITGSIMPSONINTEGRATOR_HH

//-------------
// C Headers --
//-------------
extern "C" {
}

#include "XSLItgAbsIntegrator.hh"

class XSLItgSimpsonIntegrator: public XSLItgAbsIntegrator {

public:
  
  XSLItgSimpsonIntegrator(const XSLItgAbsFunction &, double precision=1.0e-5, int maxLoop=20);

  virtual ~XSLItgSimpsonIntegrator( );
  
protected:
  
  virtual double evaluateIt(double , double) const;
  
private:
  
  double _precision;
  double _maxLoop;

  XSLItgSimpsonIntegrator();
  XSLItgSimpsonIntegrator( const XSLItgSimpsonIntegrator& );                //// Copy Constructor
  XSLItgSimpsonIntegrator& operator= ( const XSLItgSimpsonIntegrator& );    // Assignment op
  
};



#endif // XSLITGSIMPSONINTEGRATOR_HH
