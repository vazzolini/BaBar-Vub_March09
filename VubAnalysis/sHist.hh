#ifndef SHIST_H
#define SHIST_H

#include <TObject.h>
#include <TH1.h>
#include <TLatex.h>
#include <TLine.h>

#define BINMAX 40

class sHist {

public :
  sHist(); 

  // -- one mes histogram per bin 
  sHist(int base, const char *tit, int nbins, double lo, double hi);
  sHist(int base);

  // -- to create and fill in a job:
  void  bookHist(int base, const char *tit, int nbins, double lo, double hi, int *bins);
  // -- to read in from a file:
  void bookHist(int base); 

  void  setup(double *var, double *mes);
  void  fillHist(double w8); 
  TH1D* signalHist(const char *s = "cb", int prefit = 0);

  void  print();

private:

  int fBase, fNmes;
  int fVarBins[BINMAX];

  double *fpVar, *fpMes;

  TH1D *fHist;
  double fMax; 

  TH1D *fhSum;
  TH1D *fhHist[BINMAX];

};

#endif
