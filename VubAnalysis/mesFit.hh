#ifndef MESFIT_H
#define MESFIT_H

#include "TH1.h"
#include "TF1.h"
#include "TLatex.h"

#include "util.hh"

class mesFit: public TObject{

public :
  mesFit();
  mesFit(const char *s, const char *f = "gauss", int verbose = 1, int preFit = 0);
  mesFit(TH1D *h, const char *f = "gauss", int verbose = 1, int preFit = 0);
  mesFit(TH1D *h, const mesFit &fix, const char *f = "gauss", int verbose = 1);
  ~mesFit(); 

  void   hist(TH1D *h); 
  void   fitGauss(int printLevel = 0);
  void   fitCB(int printLevel = 0, int preFit = 0);
  void   print();
  void   dump();
  void   init();

  double getSig()    { return fSg;}
  double getSigE()   { return fSgE;}
  double getBg()     { return fBg;}
  double getBgE()    { return fBgE;}
  double getPur()    { return fPur;}
  double getPurE()   { return fPurE;}
  double getSb()     { return fSb;}
  double getSbE()    { return fSbE;}
  TH1D*  getHist()   { return fHist;}

  double getR()      { return fR;}
  double getRErr()   { return fRE;}

  double getSigma()  { return fSigma;}
  double getSigmaE() { return fSigmaE;}
  double getMean()   { return fMean;}
  double getMeanE()  { return fMeanE;}

  void   setSigma(double x)  { fSigma = x;}
  void   setSigmaE(double x) { fSigmaE = x;}
  void   setMean(double x)   { fMean = x;}
  void   setMeanE(double x)  { fMeanE = x;}

  void   setPrintLevel(int x = 15) { fPrintLevel = x;}
  void   setDontChange(Bool_t x = kTRUE) { fDontChange = x;}
  void   setTextSize(double x = 0.08);
  void   setMesMin(double x = 5.20) {fMesMin = x;}

  void   setTitleXY(const char *x, const char *y) {setTitleX(x); setTitleY(y);}
  void   setTitleX(const char *s) { fTitleX = TString(s);}
  void   setTitleY(const char *s) { fTitleY = TString(s);}
  void   setNdivX(int d = -405) {fNdivX = d;}
  void   setNdivY(int d = 306) {fNdivY = d;}
  void   setTitleSize(double x = 0.07) {fTitleSize = x;}
  void   setLabelSize(double x = 0.06) {fLabelSize = x;}
  void   setPadShrinkX(double x = 1.0) {fPadShrinkX = x;}
  void   setPadShrinkY(double x = 1.0) {fPadShrinkY = x;}
  void   setOffsetX(double x = 1.0) {fOffsetX = x;}
  void   setOffsetY(double x = 1.0) {fOffsetY = x;}

  void   setMarkerSize(double x = 1.0) {fMarkerSize = x;}
  void   setMarkerStyle(int x = 1) {fMarkerStyle = x;}
  void   setLineWidth(int x = 1) {fLineW = x;}

  void   setVerbose(int x = 1) {fVerbose = x;}

  void   fixSignal(double peak, double sigma, double xover = 0.9, double tail = 3.0); 
  void   fixSignal(const mesFit &); 
  void   limitSignal(double peak, double peakE, double sigma, double sigmaE, double xover, double xoverE, double tail, double tailE);
  void   limitSignal(const mesFit &); 

  double getPar(int ipar)  const   {return fPar[ipar];}
  double getParErr(int ipar) const {return fParErr[ipar];}

  void   setPar(int ipar, double value) {fPar[ipar] = value; fChangePar[ipar] = 1; fFixPar[ipar] = 0; fLimitPar[ipar] = 0;}
  void   fixPar(int ipar, double fix) {fPar[ipar]=fix; fParLo[ipar]=fix; fParHi[ipar]=fix; fChangePar[ipar]=1; fFixPar[ipar]=1;}
  void   limitPar(int ipar, double lo, double hi) {fParLo[ipar] = lo; fParHi[ipar] = hi; fLimitPar[ipar] = 1; }
  void   releasePar(int ipar) {fParLo[ipar] = 0.; fParHi[ipar] = 0.; fChangePar[ipar] = 0; fFixPar[ipar] = 0; fLimitPar[ipar] = 0;}
  void   releaseAllPar(); 

private :

  TString uniqueName(const char *);

  TH1D   *htemp, *fHist;
  TF1    *ftemp, *f0, *f1, *f2, *f3, *f4;

  double fPar[10], fParErr[10], fParLo[10], fParHi[10];
  double fParArgus[2], fParArgusE[2];
  int    fFixPar[10], fLimitPar[10], fChangePar[10];

  double fSg, fSgE; // signal     in signal region: fLo   .. fHi
  double fBg, fBgE; // background in signal region: fLo   .. fHi
  double fSb, fSbE; // integral of sideband:        fSbLo .. fSbHi
  double fPur, fPurE;
  double fR, fRE;

  double fSigma, fSigmaE;
  double fMean, fMeanE;
  double fArgus, fArgusE; 

  Bool_t fDontChange;
  int    fVerbose;
  double fLo, fHi;      // signal integration limits
  double fSbLo, fSbHi;  // side-band integration limits

  int    fPrintLevel;
  double fMesMinEvents, fMesMin;

  double fTxtSize, fTxtX, fTxtS;
  double fLabelSize, fTitleSize, fOffsetX, fOffsetY, fPadShrinkX, fPadShrinkY;
  TString fTitleX, fTitleY; 
  double fMarkerSize; 
  int    fMarkerStyle;
  int    fLineW, fNdivX, fNdivY;

  TLatex *tl; 
  char line[200];

  ClassDef(mesFit,1) // testing mesFit

};


#endif
