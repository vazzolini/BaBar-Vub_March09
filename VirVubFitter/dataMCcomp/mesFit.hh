#ifndef mesFit_hh
#define mesFit_hh

#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooDataSet.hh"
#include "RooFitCore/RooPlot.hh"
#include "TVector2.h"
#include <vector>

//class TVector2;
class RooAbsPdf;
//class RooRealVar;

enum ModelType { GaussFit = 0, ArgusAndCB = 1, ArgusAndThosig = 2};
enum parFitStdEnum { iMean=0, iSigma=1 , iAlpha=2, iN=3, iArgus=4 };
enum parFitThoEnum { iCutOff=5, iEndpoint, iThoSigR, iSigma_r1, iThoSigXc, iSigma_r2, iSigma_l, iThoSigN, iThoSigAlpha };

class mesFit : public TObject
{
public:

  mesFit(); //def ctor
  mesFit(ModelType,RooRealVar*,TFile*);
  ~mesFit(); //dtor
  TVector2* fitModel(RooDataSet*, int,vector<Double_t>&, const TString&,bool,bool);
  RooAbsPdf* createArgus(RooRealVar&);
  RooAbsPdf* createCB(RooRealVar&);
  RooAbsPdf* createThorsten(RooRealVar& mes);
  void setMesParamFile(const TString&);
  void readMesParamFile();

private:
  //  const RooDataSet *dataset;
  ModelType mesFitModel;
  RooRealVar *mes;
  TString mesparsetting;
  RooPlot *xframe;
  TFile *outfile;

  std::vector<double> mesNsl, mesdatacuts, mesvubcuts, mesvcbcuts, mesothcuts, mesNslMC, mesvcbMC, mesvubMC;
  std::vector<double> mesvubMCchop, mesvubMCall, mesvcbMCall, mesvubMClepteff, mesvubMCalleff, mesvuboutMC, mespstarMC, mespstarcuts;
    
  ClassDef(mesFit,1); //mES Fit Class
  
};

RooRealVar* getPointer(RooAbsPdf* pdf, const char* name);
Double_t getVal(RooAbsPdf* pdf, const char* name);
Double_t getError(RooAbsPdf* pdf, const char* name);

#endif
