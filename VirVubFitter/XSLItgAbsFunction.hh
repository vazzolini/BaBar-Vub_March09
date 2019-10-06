// Module: XslItgAbsFunction.hh
//
// Description:
//      Abstraction of a generic function for use in integration methods elsewhere
//      in this package. (Stolen from EvtGenModels package)
// Modified:
//      Alexei Volk  March 6, 2008
//

#ifndef XSLITGABSFUNCTION_HH
#define XSLITGABSFUNCTION_HH

//-------------
// C Headers --
//-------------
extern "C" {
}

class XSLItgAbsFunction {

public:

  // Constructors

  XSLItgAbsFunction(double lowerRange, double upperRange);

  // Destructor
  virtual ~XSLItgAbsFunction( );

  virtual double value( double x) const;

  virtual double operator()(double x) const;
  
  // Selectors (const)
  
  inline double upperRange() const {return _upperRange;}
  inline double lowerRange() const {return _lowerRange;}
  inline void   getRange(double &lower,double &upper) const { lower = _lowerRange; upper = _upperRange; } 
  virtual void setCoeff(int, int, double)=0;
  virtual double getCoeff(int, int)=0;
  
protected:
  
  virtual double myFunction(double x) const=0;
  void    setRange(double x1,double x2) { _lowerRange=x1; _upperRange=x2; };

private:
  
  double _upperRange;
  double _lowerRange;
 
  XSLItgAbsFunction( const XSLItgAbsFunction& );                // Copy Constructor
  XSLItgAbsFunction& operator= ( const XSLItgAbsFunction& );    // Assignment op

};

#endif // XSLITGABSFUNCTION_HH
