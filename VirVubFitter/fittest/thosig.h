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

#ifndef THOSIG
#define THOSIG

#include "RooAbsPdf.hh"
#include "RooRealProxy.hh"
#include "RooAbsReal.hh"
 
class thosig : public RooAbsPdf {
public:
  thosig(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _r,
	      RooAbsReal& _sigma_r1,
	      RooAbsReal& _x_c,
	      RooAbsReal& _sigma_r2,
	      RooAbsReal& _sigma_L,
	      RooAbsReal& _n,
	      RooAbsReal& _alpha);
  thosig(const thosig& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new thosig(*this,newname); }
  inline virtual ~thosig() { }

protected:

  RooRealProxy x ;
  RooRealProxy r ;
  RooRealProxy sigma_r1 ;
  RooRealProxy x_c ;
  RooRealProxy sigma_r2 ;
  RooRealProxy sigma_L ;
  RooRealProxy n ;
  RooRealProxy alpha ;
  
  Double_t evaluate() const ;

private:

  ClassDef(thosig,0) // Your description goes here...
};
 
#endif
