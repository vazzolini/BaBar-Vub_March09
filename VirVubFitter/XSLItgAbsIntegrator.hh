// Module: XSLItgIntegrator.hh
//
// Description:
//      Simpson integrator (Stolen from  EvtGenModels package)
// Modified:
//      Alexei Volk  March 6, 2008
//

#ifndef XSLITGABSINTEGRATOR_HH
#define XSLITGABSINTEGRATOR_HH


#include "XSLItgAbsFunction.hh"

class XSLItgAbsIntegrator {

public:
  
  XSLItgAbsIntegrator(const XSLItgAbsFunction &);
  
  virtual ~XSLItgAbsIntegrator( );

  double evaluate(double lower, double upper) const;
 
  double normalisation() const;

protected:

   double trapezoid(double lower, double higher, int n, 
		   double &result) const;
  
  virtual double evaluateIt(double lower, double higher) const=0;
  
  double myFunction(double x) const {return _myFunction(x);}
 
private:
  
  const XSLItgAbsFunction &_myFunction;

  void boundsCheck(double &, double &) const;

  // Note: if your class needs a copy constructor or an assignment operator, 
  //  make one of the following public and implement it.
  XSLItgAbsIntegrator();
  XSLItgAbsIntegrator( const XSLItgAbsIntegrator& );                // Copy Constructor
  XSLItgAbsIntegrator& operator= ( const XSLItgAbsIntegrator& );    // Assignment op
 
};

#endif // XSLITGABSINTEGRATOR_HH
