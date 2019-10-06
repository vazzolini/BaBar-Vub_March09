#ifndef VUBANALYSISCODE
#define VUBANALYSISCODE

#include "RecoilAnalysis/baseClass.hh"
#include "VubAnalysis/valHist.hh"

// ----------------------------------------------------------------------
class VubAnalysisCode: public baseClass {

protected:

  char    fDump4ParamFile[3][80];
  Bool_t  fisSkim;
  TTree  *fKsTree;
  Bool_t  fisDuplicate, fcontKs; 


  // ----------------------------------------------------------------------
  // -- recoil

  Int_t fEffCat, fNKshort, fBestKlind;
  Bool_t vubDepleted; 
  double totweight, totweightfRecoilTrkMult, totweightfRecoilNutMult;
  Int_t fRecoilNutMult80_160, fRecoilNutMult160_320, fRecoilNutMultfromB, 
        fRecoilNutMultfromB80_160, fRecoilNutMultfromB160_320, fRecoilPi0Mult;
  Double_t weightTrk[100], weightNeu[100];
  Double_t flMommin;
  Double_t fp4XminPi, fdeltaM, fwp4XminPi, fwdeltaM;
  Double_t fphotp4Xmin, fphotdeltaM;
  Int_t ifromBGam[80];
  Double_t tmpfMM2[15], tmpfNuT[15], fMNeupart[15];
  Double_t tmpfMxhadfit[15], tmpxmassResF[15], tmpfFxhadfit[15], tmpfTxhadfit[15], tmpfExhadfit[15];
  Int_t fcountNeu[15];

  // ----------------------------------------------------------------------
  // -- kinFit'ed quantities

  Double_t fBmassfit, fExhadfit, fMM2fit, fProbChi2, fChi2; 
  Int_t goodKshort[100], kshortLockGam[100], pi0LockGam[100];
  Int_t chargeGoodPi[60];
  Int_t nGoodPi;
  Double_t momentumGoodPi[60], thetaGoodPi[60], phiGoodPi[60], deltaMGoodPi[60];
  Int_t nGoodGam;
  Double_t momentumGoodGam[80], thetaGoodGam[80], phiGoodGam[80], deltaMGoodGam[80];
  Int_t convLockTrk[100], goodConv[100];

  Bool_t fOneLepton, 
    fGoodPRMM2, fGoodMM2, fGoodChargeCorr, fGoodNoHole[11],
    fGoodChargeCons, fGoodEvent, fGoodEventPh[15], fGoodEventPhNMNC[15];
  Bool_t fLowMx; 


  Bool_t fGoodAccLepton, fGoodAccElectron, fGoodAccMuon; 
  Bool_t fGoodLepton, fGoodElectron, fGoodMuon; 
  Bool_t fGoodRS,  fGoodWS,  fGoodRF,  fGoodWF; 
  Bool_t fGoodERS, fGoodEWS, fGoodERF, fGoodEWF; 
  Bool_t fGoodMRS, fGoodMWS, fGoodMRF, fGoodMWF; 

  Bool_t faoLepton, faoElectron, faoMuon; 
  Bool_t faoRS, faoERS, faoMRS; 
  Bool_t faoWS, faoEWS, faoMWS; 
  Bool_t faoRF, faoERF, faoMRF; 
  Bool_t faoWF, faoEWF, faoMWF; 

  Bool_t faoPRMM2,   faoPRMM2E,   faoPRMM2M; 
  Bool_t faoPRMM2RS, faoPRMM2ERS, faoPRMM2MRS; 
  Bool_t faoPRMM2WS, faoPRMM2EWS, faoPRMM2MWS; 
  Bool_t faoPRMM2RF, faoPRMM2ERF, faoPRMM2MRF; 
  Bool_t faoPRMM2WF, faoPRMM2EWF, faoPRMM2MWF; 

  Bool_t faoMM2,   faoMM2E,   faoMM2M; 
  Bool_t faoMM2RS, faoMM2ERS, faoMM2MRS; 
  Bool_t faoMM2WS, faoMM2EWS, faoMM2MWS; 
  Bool_t faoMM2RF, faoMM2ERF, faoMM2MRF; 
  Bool_t faoMM2WF, faoMM2EWF, faoMM2MWF; 

  Bool_t faoChargeCons, faoChargeConsE, faoChargeConsM;
  Bool_t faoChargeConsRS, faoChargeConsERS, faoChargeConsMRS;
  Bool_t faoChargeConsWS, faoChargeConsEWS, faoChargeConsMWS;
  Bool_t faoChargeConsRF, faoChargeConsERF, faoChargeConsMRF;
  Bool_t faoChargeConsWF, faoChargeConsEWF, faoChargeConsMWF;

  Int_t goodNr[100], goodWe[100], goodWk[100];
  float tlenTrack[100], ndchTrack[100], nsvtTrack[100], c2nTrack[100], dcaTrack[100], dcazTrack[100]; 

  // ----------------------------------------------------------------------
  // -- mcTruth

  Double_t fwCiuc, fxCiuc, fcsiCiuc, ftLepFit, ftLepG;
  Double_t fGwCiuc, fGxCiuc, fGcsiCiuc;
  Double_t fEwPwfit,fEwPw,fEwPwG;
  Int_t  fDnu;
  Int_t fWithKKbar; 
  Double_t fKplus;
  enum MCCategory { MCklong=1, MCkshortpippim=2, MCkshortpi0pi0=4,
                    MCkplus=8, MCkineff=16,  MCkmiss=32, MCcascade=64,
                    dummy1=128, dummy2=256, 
                    dummy3=512, dummy4=1024 };
  // ----------------------------------------------------------------------
  // -- Options

  Int_t fOptScaleKlongs,fOptKlongs,fOptMakeEventList;

  // ----------------------------------------------------------------------
  // -- Event

  double fLumiW8; 

  // ----------------------------------------------------------------------
  // -- MX study

#define MAXHIST 50

  valHist   *d1050[MAXHIST]; 
  valHist   *d1550[MAXHIST]; 
  valHist   *d1650[MAXHIST]; 
  valHist   *d1750[MAXHIST];

  valHist   *d2050[MAXHIST]; 
  valHist   *d2150[MAXHIST];   
  valHist   *d2250[MAXHIST];   
  valHist   *d2350[MAXHIST];   

  valHist   *d2550[MAXHIST]; 
  valHist   *d2650[MAXHIST]; 
  valHist   *d2750[MAXHIST];
  valHist   *d2850[MAXHIST]; 

  valHist   *d3050[MAXHIST]; 
  valHist   *d3150[MAXHIST]; 
  valHist   *d3550[MAXHIST]; 
  valHist   *d3650[MAXHIST];
  valHist   *d3750[MAXHIST];

  valHist   *d4050[MAXHIST]; 
  valHist   *d4150[MAXHIST]; 
  valHist   *d4250[MAXHIST];
  valHist   *d4350[MAXHIST];
  valHist   *d4450[MAXHIST];
  valHist   *d4550[MAXHIST];
  valHist   *d4650[MAXHIST];

  valHist   *d5050[MAXHIST]; // tracks
  valHist   *d5150[MAXHIST]; // photons
  valHist   *d5250[MAXHIST]; // charged kaons
  valHist   *d5350[MAXHIST]; // neutral kaons
  valHist   *d5450[MAXHIST]; // pions

  valHist   *d6050[MAXHIST]; 
  valHist   *d6150[MAXHIST];

  valHist   *d7000[MAXHIST], *d7050[MAXHIST];
  valHist   *d7100[MAXHIST], *d7150[MAXHIST];
  valHist   *d7200[MAXHIST], *d7250[MAXHIST];

  valHist   *d8000[MAXHIST], *d8050[MAXHIST];
  valHist   *d8100[MAXHIST], *d8150[MAXHIST];
  valHist   *d8200[MAXHIST]; 
  
public:

  VubAnalysisCode(TTree *tree=0,int isMC =0, int newFormat = 2);
  virtual void  bookHist(int dump);
  virtual void  mcTruth();
  void  doKlongScaling(int i) {fOptScaleKlongs = i;}
  void  runKlongs(int i) {fOptKlongs = i;}
  void  makeParam(char ifile[80]);
  void  makeEventList(int i) {fOptMakeEventList = i;}
  void  Skim(Double_t pCut, Int_t maxEvent, Int_t startEvent, Int_t isVerbose, const char *ITSfile);
  void  dumpEventList(const char *filename ="output.root");
  virtual void recoil();
  virtual void  Loop(Int_t maxEvent = 0, Int_t startEvent = 0, Int_t isVerbose = 0, Int_t lun = 0);
  void  mxCategory();
  void  dump4Param(int type, float P[28]);
  void  setLumiWeight(double w) {fLumiW8 = w;}
  void  maskKshorts(int i);
  void  maskPi0(int i);
  void  maskConversions(int modes=0 );
  void  fastBookHist(const char *name);
  void  fastFillHist(const char *name);
  void  fillMesHist(const char *dir, const char *le);
  void  fillRecoilHist(const char *);
  Double_t  kPlus(); 
  virtual void   Init(TTree *tree, int isMC);
  virtual void   initRest();
  virtual void   initVariables();
  virtual void   selectPhotons();
  Int_t  isKs2Pi0Dau(Int_t iGam);
  virtual void   selectTracks();

private:

  int ed;

};

#endif
