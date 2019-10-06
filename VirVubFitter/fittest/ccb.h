/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef CCB
#define CCB

#include "RooAbsPdf.hh"
#include "RooRealProxy.hh"
#include "RooAbsReal.hh"
#include "RooCBShape.hh" 
class ccb : public RooAbsPdf {
public:
  ccb(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _mean,
	      RooAbsReal& _sigma,
	      RooAbsReal& _alpha,
	      RooAbsReal& _n,
	      RooAbsReal& _endpoint);
  ccb(const ccb& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new ccb(*this,newname); }
  inline virtual ~ccb() { }

protected:

  RooRealProxy x ;
  RooRealProxy mean ;
  RooRealProxy sigma ;
  RooRealProxy alpha ;
  RooRealProxy n ;
  RooRealProxy endpoint ;
  Double_t evaluate() const ;

private:

  ClassDef(ccb,0) // Your description goes here...
};
 
#endif
