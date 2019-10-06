/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: VirVubFitter                                                     *
 *    File: $Id: RooSumArgusBG.cc,v 1.1 2006/05/15 08:42:13 menges Exp $
 *                                                                           *
 * Authors:                                                                  *
 *   Wolfgang Menges, Queen Mary, University of London,                      *
 *                                            menges@slac.stanford.edu       *
 *                                                                           *
 * Copyright (c) 2006                                                        *
 *                                                                           *
 *****************************************************************************/

// -- CLASS DESCRIPTION [PDF] --
/*
  this class implements a sum of argus functions, where each has a 
  weight and a fixed endpoint
*/

#include "VirVubFitter/RooSumArgusBG.hh"

#include "iostream"
#include <math.h>

#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooRealConstant.hh"
#include "RooFitCore/RooMath.hh"
#include "TMath.h"

ClassImp(RooSumArgusBG)

RooSumArgusBG::RooSumArgusBG(const char *name, const char *title,
				 RooAbsReal& _m, RooAbsReal& _m0, RooAbsReal& _c,
				 std::vector<double> weight, std::vector<double> endpoint) :
  RooAbsPdf(name, title), 
  m("m","Mass",this,_m),
  m0("m0","Resonance mass",this,_m0),
  c("c","Slope parameter",this,_c),
  p("p","Power",this,(RooRealVar&)RooRealConstant::value(0.5))
{
  if (weight.empty() && endpoint.empty()) {
    _weight.push_back(1.);
    _endpoint.push_back(0.);
  } else if (weight.size() == endpoint.size()) {
    _weight = weight;
    _endpoint = endpoint;
  } else {
    std::cout << "Warning: weight vector does not match endpoint vector, reverting to single Argus!" << std::endl;
    _weight.push_back(1.);
    _endpoint.push_back(0.);
  }
}

RooSumArgusBG::RooSumArgusBG(const char *name, const char *title,
		       RooAbsReal& _m, RooAbsReal& _m0, RooAbsReal& _c, RooAbsReal& _p) :
  RooAbsPdf(name, title), 
  m("m","Mass",this,_m),
  m0("m0","Resonance mass",this,_m0),
  c("c","Slope parameter",this,_c),
  p("p","Power",this,_p)
{
}

RooSumArgusBG::RooSumArgusBG(const RooSumArgusBG& other, const char* name) :
  RooAbsPdf(other,name), 
  m("m",this,other.m), 
  m0("m0",this,other.m0), 
  c("c",this,other.c),
  p("p",this,other.p),
  _weight(other._weight),
  _endpoint(other._endpoint)
{
}

Double_t RooSumArgusBG::evaluate(const Double_t x, const Double_t endpoint) const {
  Double_t t = x/endpoint;
  if(t >= 1) return 0;

  Double_t u= 1 - t*t;
  //cout << "c = " << c << " result = " << m*TMath::Power(u,p)*exp(c*u) << endl ; 
  return m*TMath::Power(u,p)*exp(c*u);
}

Double_t RooSumArgusBG::evaluate() const {
  Double_t sum = 0.;
  for (size_t i = 0; i < _weight.size(); ++i) {
    sum += _weight[i]*evaluate(m, m0+_endpoint[i]);
  }  

  return sum;
}
