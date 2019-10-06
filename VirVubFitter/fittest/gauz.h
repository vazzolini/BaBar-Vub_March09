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

#ifndef GAUZ
#define GAUZ

#include "RooAbsPdf.hh"
#include "RooRealProxy.hh"
#include "RooAbsReal.hh"
 
class gauz : public RooAbsPdf {
public:
  gauz(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _x_max,
	      RooAbsReal& _sigma1,
	      RooAbsReal& _sigma2,
              RooAbsReal& _sigma3);
  gauz(const gauz& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new gauz(*this,newname); }
  inline virtual ~gauz() { }

protected:

  RooRealProxy x ;
  RooRealProxy x_max ;
  RooRealProxy sigma1 ;
  RooRealProxy sigma2 ;
  RooRealProxy sigma3 ;
  
  Double_t evaluate() const ;

private:

  ClassDef(gauz,0) // Your description goes here...
};
 
#endif
