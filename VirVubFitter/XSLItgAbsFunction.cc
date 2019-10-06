// Module: XslItgAbsFunction.cc
//
// Description:
//      Abstraction of a generic function for use in integration methods elsewhere
//      in this package. (Stolen from EvtGenModels package)
// Modified:
//      Alexei Volk  March 6, 2008
//


#include "XSLItgAbsFunction.hh"

//-------------
// C Headers --
//-------------
extern "C" {
}
#include "assert.h"
#include <iostream>

using std::endl;
using std::cout;

XSLItgAbsFunction::XSLItgAbsFunction(double lowerRange, double upperRange):
  _upperRange(upperRange),
  _lowerRange(lowerRange){}

XSLItgAbsFunction::~XSLItgAbsFunction( )
{}


double
XSLItgAbsFunction::value( double x) const{
  if (x >= _lowerRange && x <= _upperRange) return myFunction(x);
   cout << "Error in XSLItgAbsFunction::value.  Given co-ordinate " << x
                << " is outside of allowed range [" << _lowerRange << ", "
                << _upperRange << "].  Returning 0.0" << endl;
  return 0.0;  // Never get here
}
   
double 
XSLItgAbsFunction::operator()(double x) const{
  return myFunction(x);
}
