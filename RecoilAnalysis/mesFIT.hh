#ifndef MESFIT_H
#define MESFIT_H

#include "TH1.h"
#include "TF1.h"
#include "TLatex.h"

class mesFIT {

public :
  mesFIT(int histNameModifier=0);
  mesFIT(TH1D *h, const char *f = "cb", int setting = 0, int histNameModifier=0, int algorithm = 999, double en_beam=10.5795, double arg_uplim=5.26);
  void   hist(TH1D *h); 
  void   fitGauss(int printLevel = 0, int setting = 0);
  void   fitCB(int printLevel = 0, int setting = 0, int newAlgorithm=999);
  void   print();
  void   init(int pass=0);

  double getSig()    { return fSg;}
  double getSigE()   { return fSgE;}
  double getBg()     { return fBg;}
  double getBgE()    { return fBgE;}
  double getPur()    { return fPur;}
  double getPurE()   { return fPurE;}
  double getSb()     { return fSb;}
  double getSbE()    { return fSbE;}

  double getR()      { return fR;}
  double getRErr()   { return fRE;}

  double getSigma()  { return fSigma;}
  double getSigmaE() { return fSigmaE;}
  double getMean()   { return fMean;}
  double getMeanE()  { return fMeanE;}

  double getPar(int ipar) {return f3->GetParameter(ipar);}
  double getChi2() {return fchi2;}

  TH1D *GetResidual(int mode);

  //  void setFitMode(int mode) {fFitMode=mode;}

  void   setSigma(double x)  { fSigma = x;}
  void   setSigmaE(double x) { fSigmaE = x;}
  void   setMean(double x)   { fMean = x;}
  void   setMeanE(double x)  { fMeanE = x;}

  void   setPrintLevel(int x = 15) { fPrintLevel = x;}
  void   setDontChange(Bool_t x = kTRUE) { fDontChange = x;}
  void   setTextSize(double x = 0.08);
  void   setMesMin(double x = 5.2) {fMesMin = x;}
  void   setPars(double mesmin = 5.2, double messglow=5.27, double messghi=5.29, double mesbglow=5.21, double mesbghi=5.26, double ebeam=10.5795) {fMesMin = mesmin; fLo=messglow; fHi=messghi; fSbLo=mesbglow; fSbHi=mesbghi; fEbeam=ebeam;}

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

  void   setPar(int ipar, double value) {fPar[ipar] = value; fChangePar[ipar] = 1; fFixPar[ipar] = 0; fLimitPar[ipar] = 0;}
  void   fixPar(int ipar, double fix) {fPar[ipar] = fix; fParLo[ipar] = fix; fParHi[ipar] = fix; fFixPar[ipar] = 1; fChangePar[ipar] = 1;}
  void   limitPar(int ipar, double lo, double hi) {fParLo[ipar] = lo; fParHi[ipar] = hi; fLimitPar[ipar] = 1; }
  void   releasePar(int ipar) {fParLo[ipar] = 0.; fParHi[ipar] = 0.; fChangePar[ipar] = 0; fFixPar[ipar] = 0; fLimitPar[ipar] = 0;}

  void setTitles(TH1 *h, const char *sx, const char *sy, float size = 0.05, 
                 float xoff = 1.1, float yoff = 1.1, float lsize = 0.05, int font = 62);
  void setHist(TH1 *h, int color = kBlack, int symbol = 20, double size = 1., double width = 2.);
  void shrinkPad(double b = 0.1, double l = 0.1, double r = 0.1, double t = 0.1);

private :
  TH1D   *htemp, *fHist, *fitHist;
  TF1    *ftemp, *f0, *f1, *f2, *f3, *f4;

  int _histNameModifier;
  int _newAlgorithm;

  double fEbeam;

  int fFitMode;

  double fPar[10], fParLo[10], fParHi[10];
  double fPar_AR[2], fPar_ARLo[2], fPar_ARHi[2], fPar_ARE[2];
  double fchi2;
  int    fFixPar[10], fLimitPar[10], fChangePar[10];

  double fR, fRE;
  double fSg, fSgE; // signal     in signal region: fLo   .. fHi
  double fBg, fBgE; // background in signal region: fLo   .. fHi
  double fSb, fSbE; // integral of sideband:        fSbLo .. fSbHi
  double fPur, fPurE;

  double fSigma, fSigmaE;
  double fMean, fMeanE;

  double farg1, farg2;
  double farg1e, farg2e;

  Bool_t fDontChange;
  int    fVerbose;
  double fLo, fHi;      // signal integration limits
  double fSbLo, fSbHi;  // side-band integration limits

  int    fPrintLevel;
  double fMesMinEvents, fMesMin;

  double fTxtSize, fTxtX, fTxtS;
  double fLabelSize, fTitleSize, fOffsetX, fOffsetY, fPadShrinkX, fPadShrinkY;
  double fMarkerSize; 
  int    fMarkerStyle;
  int    fLineW;

  TLatex *tl; 
  char line[200];

//#ifndef __CINT__ 
  ClassDef(mesFIT,1) // testing mesFIT
//#endif
};


#endif
