#ifndef b2uClass_h
#define b2uClass_h

#include <vector>

#include "VirVubFitter/VirClass.hh" 

class b2uClass: public VirClass {

public :

  RooRealVar *VlepmxYes, *VmxYes, *VlepYaSe2, *VlepYaSe3;
  RooDataSet *datamcvubNC, *datamcvcbNC;
  RooDataSet *datamcvubtrue, *datamcvubfalse;
  RooDataSet *datamcvub_mx, *datamcvub_r;
  RooDataSet *datavubin_mx, *datavubin_r;
  double vubmcmx, errvubmcmx, vubmcr, errvubmcr, epsumx, errepsumx, epssel, errepssel, epsselbr, errepsselbr;
  double vubmcNC, errvubmcNC, vcbmcNC, errvcbmcNC;
  double N1, N2, N3, N1err, N2err, N3err;
  int fusemxfit;
  float fweights[10];

  ~b2uClass();
  b2uClass();
  b2uClass(TTree *tree, TString filename, int Sys, int q2F, int comb, int un, int unfB, double hiunfB, int me, int mu, int usemxfit, bool iscm2,bool varfit, const std::vector<float>& wFermivec, int newbin);

  void Loop(int isdata, int icat, int isMC, int nres, int comb);
  void b2uresultDumping();
  void b2ueff(int su);
  void b2upstar(int su);
  void b2umxBinning1d(double width);
  void b2uGetWeights(TString b2uwFile);
  int  b2unewrHistmx(double cmx);

};

#endif
