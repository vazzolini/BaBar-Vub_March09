#ifndef VALHIST_H
#define VALHIST_H

#include "TObject.h"
#include <TH1.h>

#define MESMAX 40

class valHist {

public :
  valHist(); 

  // -- one mes histogram per bin 
  valHist(int base, const char *title, int nbins, double lo, double hi);

  // -- one mes histogram per group of bins
  valHist(int base, const char *title, int nbins, double lo, double hi, int nmes, int *bins);
  valHist(int base, TH1D *zb, int nmes = 1, int *bins = 0);

  void bookHist(int base, const char *title, int nbins, double lo, double hi, int nmes, int *bins);

  void setup(double *var, double *mes, Bool_t *signalBox, Bool_t *sideband, 
	     Bool_t *nocuts, Bool_t *allcuts, Bool_t *aocuts); 

  void fillHist(double w8 = 1.);
  void fillHist(int ntracks, float *tracks, int *goodTracks, double w8 = 1.);

  void print();


private:
  double *fpVar, *fpMes; 
  Bool_t *fpSignalBox, *fpSideband; 
  Bool_t *fpNoCuts, *fpSigCuts, *fpAoCuts; 
  
  int fBase; 
  int fNmes;
  int fVarBins[MESMAX];

  TH1D *a00,    *a10,    *a20;
  TH1D *a00sg,  *a10sg,  *a20sg; 
  TH1D *a00bg,  *a10bg,  *a20bg; 
  TH1D *a00mes[MESMAX], *a10mes[MESMAX], *a20mes[MESMAX]; 

};

#endif
