#ifndef recoilDSys_h
#define recoilDSys_h

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
#include <Rtypes.h>

#define MAXNDMODES 100
#define MAXNIMODES 4

class recoilDSys {

private: 
  
  int _nPiMode[MAXNDMODES];
  int _nKMode[MAXNDMODES];
  int _nK0Mode[MAXNDMODES];
  int _nPi0Mode[MAXNDMODES];
  float _weightMode[MAXNDMODES];
  int _nmodes;
  float _b0weights[6];
  float _bchweights[6];

  float _b0WeightsAllSemilep;
  float _bchWeightsAllSemilep;

  void  readFile(const char *InputFileName, int seed=0, int iMode=0);
  float randomized(float mean,float sig, TRandom& rndm, int seed);
  float Brandomized(float mean, float lsig, float hsig, TRandom& rndm, int seed);

public :
 // seed is the random number generator for the systematics. the default, 0, is NO smearing
  recoilDSys(const char *InputFileName, int seed=0, int imode=0); // constructor for D decays (list from file)
  recoilDSys(int seed=0, int rel=14, int sysbdec=-99); // constructor for B->Dlnu decays (numbers hard coded for rel 14 (effectively SP4, SP5/SP6) and rel 18 (effectively SP8)

  virtual ~recoilDSys();

  void recoilDSys2(int num=0); // Fake constructor for B->Dlnu decays (read weights from file)
  void recoilDSys3(int num=0); // Fake constructor for D decays (read weights from file)

  float weight(int npi,int nk, int nk0, int npi0, int nlep, int imOde, int & mode);// for D decays  
  float weight(int bdecay);// for B decays

  float getAllSemilepWeight(bool); //to get the first weight.

  void reweightForSys(int); //Recomputes weight for B->Dlv systematic studies
  
  ClassDef(recoilDSys,1) 
};


#endif
