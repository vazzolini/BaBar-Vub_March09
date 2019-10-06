// Module: XslItgPtrFunction.cc
//
// Description:
//      Class describing a function with one vector of coefficients. 
//      (Stolen from EvtGenModels package).
//Modified: 
//      Alexei Volk       March 6, 2008
 
#include "XSLItgPtrFunction.hh"

//-------------
// C Headers --
//-------------
extern "C" {
}

//----------------
// Constructors --
//----------------
XSLItgPtrFunction::XSLItgPtrFunction( double (*theFunction)(double, const std::vector<double> &), double lowerRange, double upperRange, const std::vector<double> &coeffs1):
  XSLItgAbsFunction(lowerRange, upperRange),
  _myFunction(theFunction),
  _coeffs1(coeffs1)
{}


//--------------
// Destructor --
//--------------

XSLItgPtrFunction::~XSLItgPtrFunction( )
{}


double
XSLItgPtrFunction::myFunction(double x) const{
  return _myFunction(x, _coeffs1);
}

void
XSLItgPtrFunction::setCoeff(int vect, int which, double value)
{
  if (vect == 1) _coeffs1[which] = value;
}

double
XSLItgPtrFunction::getCoeff(int vect, int which)
{
  if (vect == 1) return _coeffs1[which];
  else {return 0;}
}
