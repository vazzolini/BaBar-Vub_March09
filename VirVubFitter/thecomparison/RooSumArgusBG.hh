/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: VirVubFitter                                                     *
 *    File: $Id: RooSumArgusBG.hh,v 1.2 2006/05/25 21:53:26 menges Exp $
 *                                                                           *
 * Authors:                                                                  *
 *   Wolfgang Menges, Queen Mary, University of London,                      *
 *                                            menges@slac.stanford.edu       *
 *                                                                           *
 * Copyright (c) 2006                                                        *
 *                                                                           *
 *****************************************************************************/
#ifndef ROO_SUM_ARGUS_BG
#define ROO_SUM_ARGUS_BG

#include <vector>

#include "RooFitCore/RooAbsPdf.hh"
#include "RooFitCore/RooRealProxy.hh"

class RooRealVar;
class RooAbsReal;

class RooSumArgusBG : public RooAbsPdf {
public:
  RooSumArgusBG(const char *name, const char *title, 
		RooAbsReal& _m, RooAbsReal& _m0, RooAbsReal& _c,
		std::vector<double> weight = std::vector<double>(), std::vector<double> endpoint = std::vector<double>());
  RooSumArgusBG(const char *name, const char *title, 
		RooAbsReal& _m, RooAbsReal& _m0, RooAbsReal& _c, RooAbsReal& _p);
  RooSumArgusBG(const RooSumArgusBG& other,const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooSumArgusBG(*this,newname); }
  inline virtual ~RooSumArgusBG() { }

protected:
  RooRealProxy m ;
  RooRealProxy m0 ;
  RooRealProxy c ;
  RooRealProxy p ;

  std::vector<double> _weight;
  std::vector<double> _endpoint;

  Double_t evaluate(const Double_t x, const Double_t endpoint) const;
  Double_t evaluate() const ;
//   void initGenerator();

private:
  ClassDef(RooSumArgusBG,0) // Sum of Argus background shape PDF with fixed endpoints
};

#endif
