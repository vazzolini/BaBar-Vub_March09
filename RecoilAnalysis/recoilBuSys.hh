#ifndef recoilBuSys_h
#define recoilBuSys_h

#include "../RecoilAnalysis/mesData.hh"

#include <TVector3.h>
#include <TLorentzVector.h>
#include <TChain.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <TObjString.h>
#include <TTree.h>
#include <TLatex.h>
#include <TRandom.h>

#include <TROOT.h>

class recoilBuSys {

private: 
  
  float _bweights[7], xv[7], xe[7], xt[7];


  void  readFile(const char *InputFileName);
  float randomized(float mean,float sig, TRandom& rndm, int seed);

public :
 // seed is the random number generator for the systematics. the default, 0, is NO smearing
  recoilBuSys(const char *InputFileName, int seed=0); // constructor for B->Dlnu decays (numbers hard coded)

  ~recoilBuSys(){}

  float weight(int bdecay);// for B decays
 
};


#endif
