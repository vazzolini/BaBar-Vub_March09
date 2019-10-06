#ifndef mesData_h
#define mesData_h

#include <iostream>

#include "TString.h"

class mesData {

public :
  mesData();
  mesData(TString name, double sig, double sigE, double pu, double puE, double bg, double bgE);
  mesData(double sig);
  mesData(const mesData&);
  mesData& operator = ( const mesData& );
  virtual ~mesData();
  TString theName();
  double theSig();
  double theErrSig();
  double theBg();
  double theErrBg();
  double thePu();
  double theErrPu();

  double getSigma()  { return fSigma;}
  double getSigmaE() { return fSigmaE;}
  double getMean()   { return fMean;}
  double getMeanE()  { return fMeanE;}

  void setSigma(double x)  { fSigma = x;}
  void setSigmaE(double x) { fSigmaE = x;}
  void setMean(double x)   { fMean = x;}
  void setMeanE(double x)  { fMeanE = x;}

private :
  TString _name;
  double _sig;
  double _sigE;
  double _bg;
  double _bgE;
  double _pu;
  double _puE;

  double fSigma, fSigmaE;
  double fMean, fMeanE;

  // ----------------------------------------------------------------------
   ClassDef(mesData,1) 

};


#endif
