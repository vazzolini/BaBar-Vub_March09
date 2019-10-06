// Module: XslItgPtrFunction.hh
//
// Description:
//      Class describing a function with one vector of coefficients. 
//      (Stolen from EvtGenModels package).
//Modified: 
//      Alexei Volk       March 6, 2008


#ifndef XSLITGPTRFUNCTION_HH
#define XSLITGPTRFUNCTION_HH

#include <vector>
#include "XSLItgAbsFunction.hh"

class XSLItgPtrFunction: public XSLItgAbsFunction {

public:

  XSLItgPtrFunction( double (*theFunction)(double, const std::vector<double> &),
		     double lowerRange, double upperRange, const std::vector<double> &coeffs1);
 
  virtual ~XSLItgPtrFunction( );

  virtual void setCoeff(int, int, double);
  virtual double getCoeff(int, int);

protected:
  
  virtual double myFunction(double x) const;
 
private:
 
  // Data members
  double (*_myFunction)(double x, const std::vector<double> & coeffs1);

  // Note: if your class needs a copy constructor or an assignment operator, 
  //  make one of the following public and implement it.
  XSLItgPtrFunction( const XSLItgPtrFunction& );                //// Copy Constructor
  XSLItgPtrFunction& operator= ( const XSLItgPtrFunction& );    // Assignment op
  std::vector<double> _coeffs1;

};

#endif // XSLITGPTRFUNCTION_HH
