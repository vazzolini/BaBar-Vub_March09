// Module: XSLItgIntegrator.cc
//
// Description:
//      Simpson integrator (Stolen from  EvtGenModels package)
// Modified:
//      Alexei Volk  March 6, 2008
//
 

#include "XSLItgAbsIntegrator.hh"


//-------------
// C Headers --
//-------------
extern "C" {
}

#include <math.h>
#include <iostream>

#include "XSLItgAbsFunction.hh"
using std::endl;
using std::cout;

XSLItgAbsIntegrator::XSLItgAbsIntegrator(const XSLItgAbsFunction &theFunction):
  _myFunction(theFunction)
{}

XSLItgAbsIntegrator::~XSLItgAbsIntegrator()
{}
  
double 
XSLItgAbsIntegrator::normalisation() const {
  return evaluateIt(_myFunction.lowerRange(), _myFunction.upperRange());
}

double
XSLItgAbsIntegrator::evaluate(double lower, double upper) const{

  double newLower(lower), newUpper(upper);

  boundsCheck(newLower, newUpper);

  return evaluateIt(newLower, newUpper);
}

double 
XSLItgAbsIntegrator::trapezoid(double lower, double higher, int n, double &result) const {

  if (n==1) return 0.5*(higher-lower)*(_myFunction(lower) + _myFunction(higher));
  
  int it, j;
  
  for (it=1, j=1;j<n-1;j++) it <<=1;
  
  double itDouble(it);
  
  double sum(0.0);

  double deltaX((higher - lower)/itDouble);
  
  double x(lower + 0.5* deltaX);
    
  for (j=1;j<=it;j++){
    sum+=_myFunction(x);
    x+=deltaX;
  }
  
  result = 0.5*(result+(higher - lower)*sum/itDouble);

  return result;
}

void 
XSLItgAbsIntegrator::boundsCheck(double &lower, double &upper) const{

  if (lower < _myFunction.lowerRange() ) {
    cout << "Warning in XSLItgAbsIntegrator::evaluate.  Lower bound " << lower << " of integral " 
		    << " is less than lower bound " << _myFunction.lowerRange() 
		    << " of function.  No contribution from this range will be counted." << endl;
    lower = _myFunction.lowerRange();
  }

  if (upper > _myFunction.upperRange() ) {
    cout << "Warning in XSLItgAbsIntegrator::evaluate.  Upper bound " << upper << " of integral "
		    << " is greater than upper bound " << _myFunction.upperRange() 
		    << " of function.  No contribution from this range will be counted." << endl;  
    upper = _myFunction.upperRange();
  }

}
