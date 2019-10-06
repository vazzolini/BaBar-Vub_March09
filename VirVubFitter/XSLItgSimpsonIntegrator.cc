// Module: EvtItgSimpsonIntegrator.cc
//
// Description:
//      Abstraction of a generic function for use in integration methods elsewhere
//      in this package. 
//      (Stolen from  EvtGenModels package)
// Modified:
//      Alexei Volk  March 6, 2008
//

#include "XSLItgSimpsonIntegrator.hh"

//-------------
// C Headers --
//-------------
extern "C" {
}

//---------------
// C++ Headers --
//---------------

#include <math.h>
#include <iostream>


//-------------------------------
// Collaborating Class Headers --
//-------------------------------

#include "XSLItgAbsFunction.hh"
using std::endl;
using std::cout;


XSLItgSimpsonIntegrator::XSLItgSimpsonIntegrator(const XSLItgAbsFunction &theFunction, double precision, int maxLoop):
  XSLItgAbsIntegrator(theFunction),
  _precision(precision),
  _maxLoop(maxLoop)
{}


//--------------
// Destructor --
//--------------

XSLItgSimpsonIntegrator::~XSLItgSimpsonIntegrator()
{}

double
XSLItgSimpsonIntegrator::evaluateIt(double lower, double higher) const{
  
  // report(INFO,"XSLGen")<<"in evaluate"<<endl;
  int j;
  double result(0.0);
  double s, st, ost(0.0);
  for (j=1;j<4;j++) {
    st = trapezoid(lower, higher, j, result);
    s = (4.0 * st - ost)/3.0;
    ost=st;
  }

  double olds(s);
  st = trapezoid(lower, higher, j, result);
  s = (4.0 * st - ost)/3.0;

  if (fabs(s - olds) < _precision*fabs(olds) || (s==0.0 && olds==0.0))     return s;
  
  ost=st;

  for (j=5;j<_maxLoop;j++){

    st = trapezoid(lower, higher, j, result);
    s = (4.0 * st - ost)/3.0;
    
    if (fabs(s - olds) < _precision*fabs(olds) || (s==0.0 && olds==0.0))    return s;
    olds=s;
    ost=st;
  }
  
  cout << "Severe error in XSLItgSimpsonIntegrator.  Failed to converge after loop with 2**"
		 << _maxLoop << " calls to the integrand in." << endl;
  
  return 0.0;
    
}
