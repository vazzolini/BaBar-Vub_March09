#ifndef b2uNtp_h
#define b2uNtp_h

#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TVector3.h>

#include "RecoilAnalysis/PIDTable.hh"

class b2uNtp {

 public :

  void   doBreco(); 
  void   cleanupTracks(); 
  void   selectTracks(); 
  void   doSplitOffStudy();
  void   selectPhotons(); 
  void   selectKlongs(); 
  void   selectKshorts(); 
  void   selectKshortsZZ(); 

  void   recoil(); 
  void   fillPidMaps(); 
  void   findLeadingLepton(); 

  void   doBremRecovery();
  void   doPartialB0();
  void   doPartialB0bis();
  void   doPartialBp();
  void   doDeltaM();
  void   smearNeut(); 
  double clusterReCorrection(const double oldCorrectedEnergy, 
			     const double clusterPositionTheta) const ;
  double clusterCorrection(const double rawEnergy,
			   const double clusterPositionTheta, 
			   const bool   newCorr) const;

  //For the breco quality  
  void goodBreco(int lbrecoI);
  int getBflavor(int Bid);
  int getBflavor(int BRECOflavor, int charge);
  int getSeedDbyimc(int inimc);
  int getSeedDbymode(int inmodeB);

  void mcTruth(); 
  int  GetChg(int lundid); //For truth multiplicities
  void GetCat(int p, int VIndex, double dist); //For truth multiplicities
  void printLorentz(const TLorentzVector &);

  // -- temporary functions for tests
  void   doPi0(); 
  void   doKshorts(); // K0S ->pi+pi-
  void   doKszz();    // K0S ->pi0pi0
  void   doKlongs(); 
  void   doTracks(); 

  // -- Utility functions
  void   mk4Vector(TLorentzVector &p4, const Float_t p, const Float_t t, const Float_t f, const Float_t m); 
  void   mk3Vector(TVector3 &p3, const Float_t p, const Float_t t, const Float_t f); 

  double unwrapped(double phi1, double phi2); 
  int    isRecoed(int imc); 
  Bool_t isAncestor(int ancestor, int candidate); 
  Bool_t isStable(int id); 
  void   dumpGeneratorBlock(int b1 = -1, int b2 = -1);
  void   getMcDaughters(int imc, int &d1, int &d2, int &d3, int &d4, int &d5); 
  void   setupLeadingLepton(double lmass); 

  double kPlus(); 
  Bool_t isTruLepton(int i);
  Bool_t isTruEl(int i);
  Bool_t isTruMu(int i);
  Bool_t isTruTau(int i);

  Bool_t isRecEl(int i);    
  Bool_t isRecMu(int i);    
  Bool_t isRecLepton(int i);
  Bool_t isRecKaon(int i);

  Bool_t isPidKillEl(int i);
  Bool_t isPidKillMu(int i);
  Bool_t isPidKillKaon(int i);

  double lhKlong00(int i); 
  double lhKlong01(int i); 
  double lhKlong02(int i); 



  void   getPidTables(); 
  void   readCuts(TString filename);
  void   setup(const char *histFilename, int dumpLevel); 
  void   closeHistFile(); 


  // -- Event stuff
  int     fIsMC, fVerbose;  // ...
  int     fEvent;           // the event number in the chain
  int     fDump;            // control amount of output into root file
  TString fRunRange;        // Run 1, Run 2, Run 3, ...
  int     fFileChanged;     // a flag to allow a faster comparison if the underlying file has changed

  // -- Options
  int    fOptPidKillingEl, fOptPidKillingMu, fOptPidKillingKa, fOptSmearNeut; 


  // -- BRECO
  int    fChB;          // 1 for charged B
  int    fSeedMode;     // 0 - dc, 1 - dstar, 2 - d0, 3 - dstar0
  int    fIndexBestB;   // index to best BRECO candidate
  int    fBrecoOverlap; // mask to check for overlap with BRECO
  double fPurity, fIntPurity; 
  double fMes, fPcmsBreco; 
  int    fBrecoFlavor, fBrecoCharge, fBmode; 

  TLorentzVector f4Breco, f4BrecoNC, f4BrecoGen; 
  TLorentzVector f4Upsilon, f4Brecoil;
  TVector3 f3CmsBoost, f3UpsBoost;

  // -- Recoil
  int    fLeptonCharge, fLeptonIndex; 
  int    fNlepton, fElectron, fMuon; 
  int    fNEl, fNMu, flepmap, fBadReco; 
  double fPlab, fTlab, fTlabDR, fFlab, fFlabDR, fPcms, fTcms, fFcms, fEcms;
  double fpcmsTrklo;

  int    fB2uDepleted; 
  int    fNKs, fNKz, fNKl, fNKp;  // counters for KS->pi+pi-, KS->pi0pi0, KL, K+
  int    fQtot, fRecoilCharge, fRecoilTrkMult, fRecoilGamMult, brecoI;
  double fPxhad, fTxhad, fFxhad, fExhad, fMxhad, fMxtrk, fMxnut, fPcmsTrkLo; 
  double fPNu, fTNu, fCosTNu, fFNu, fMM2, fEmiss, fQ2, fMM2NC;
  double fMM1pr, fMM2pr, fMM3pr, fDeltaM, fCosBY; 
  double fOA1, fOA2, fOA3;

  TLorentzVector f4Lepton, f4Recoil, f4Xhad, f4Xnut, f4Xtrk, f4Neutrino, f4NeutrinoNC;
  TLorentzVector f4LeptonCMS;

  double fProbChi2, fChi2, fMxhadfit, fMM2fit, fQ2fit; 
  

  // -- Mask arrays
  int    fBrecoTrk[100], fBrecoPhoton[100]; 
  int    fRecEl[100], fRecMu[100], fRecKa[100]; 
  int    fGoodElectronTrk[100], fGoodMuonTrk[100], fGoodPionTrk[100], fGoodKaonTrk[100]; 
  int    fCleanGoodTrack[100]; 
  int    fGoodPhoton[100]; 
  int    fGoodKs[100]; 
  int    fGoodKl[100]; 


  int    fVub, fVcb, fOther;


  double fPvGen, fTvGen, fFvGen;
  double fPvcmsGen, fTvcmsGen, fFvcmsGen; 


  // split off quantities
  Double_t SOdth,SOdph,SOegam,SOptrk,SOthtrk,SOphitrk,SOthgam,SOphigam,SOecaltrk;
  Int_t SOcharge;
  Int_t splitOffGam[100];

  // ----------------------------------------------------------------------
  // -- mcTruth (TO BE REPLACED!!)
  TLorentzVector p4XhadGen, p4XtrkGen, p4XnutGen, p4XinvGen, p4MissGen, p4RecoilGen, p4LeptonGen, p4LeptonGenwPh, p4XhadGenwoPh; 
  Double_t fQ2Gen, fKplus;
  Double_t fMxhadGen, fPxhadGen, fTxhadGen, fFxhadGen, fExhadGen;
  Double_t fMxtrkGen, fPxtrkGen, fTxtrkGen, fFxtrkGen, fExtrkGen;
  Double_t fMxnutGen, fPxnutGen, fTxnutGen, fFxnutGen, fExnutGen;
  Double_t fwCiuc, fxCiuc, fcsiCiuc, ftLepFit, ftLepG;
  Double_t fGwCiuc, fGxCiuc, fGcsiCiuc;
  Double_t fEwPwfit,fEwPw,fEwPwG;
  Int_t fBbchgen, fQb;
  Double_t fPcmsGen, fTcmsGen, fFcmsGen, fEcmsGen;
  Int_t    fB1Index, fB2Index, fBVxb, fBVcb, fBVub, fBVxbTyp, fMxCategory;
  Int_t  fDpi, fDpiz, fDk,  fDks, fDCfDs, fD0CfDs, fDlep, fDgam, fBVSTyp;
  Int_t  fDkl, fDkspiopio, fDkspipi, fDnu;
  Double_t fctvGen,fctlGen,fchiGen;
  Double_t fMxhadGenwoPh, fPxhadGenwoPh, fTxhadGenwoPh, fFxhadGenwoPh, fExhadGenwoPh; //truth without bremsstr. photons

  Int_t fNLeptonMC, genLeptonCharge;
  Int_t fWithKKbar;
  Int_t fBrecoQual; //Breco quality
  Int_t fchgDau, fneuDau; //truth multiplicities
  // ----------------------------------------------------------------------


  // -- ROOT output
  TFile  *fHistFile;
  TTree  *fTree;


  // -- CUTS
  double PURITY, INTPURITY, IPURDC, IPURDSTAR, IPURD0, IPURDSTAR0;
  int    IDEL, IDMU, IDKA;
  double ELMOMLO, ELTHETALO, ELTHETAHI, ELBREMCONE; 
  double MUMOMLO, MUTHETALO, MUTHETAHI;
  double PIMOMLO, PITHETALO, PITHETAHI; 
  double KAMOMLO, KATHETALO, KATHETAHI;
  double GAMMAELO, GAMMATHETALO, GAMMATHETAHI;
  double SIGMANEUT;
  double TRACKKILL, NEUTKILL;

  int    TRACKSELECTION, PHOTONSELECTION, KLSELECTION, KSSELECTION, KSZZSELECTION, DOTRACKKILLING, DONEUTKILLING; 

  // -- Constants
  int    BZLUND, BPLUND; 
  double ELMASS, MUMASS, PIPMASS, PIZMASS, KAPMASS, KAZMASS, BPMASS, BZMASS, BMASS; 
  double BQMASS, A0;
  double DR; // = 180/pi

  // -- PidTables
  // five particles with two charges
  //  PIDTable *fPTel[10], *fPTmu[10], *fPTpi[10], *fPTka[10],  *fPTpr[10]; 
  // 3 DCH voltages ranges for SP3 and SP4
  //  TRKTable *fTT[5]; 
  PIDTable *fPTel[10], *fPTmu[10], *fPTka[10];
  double SHIFTELPIDTABLES, SHIFTMUPIDTABLES, SHIFTKAPIDTABLES, 
    SHIFTELMISTABLES, SHIFTMUMISTABLES, SHIFTKAMISTABLES; 
  char PIDTABLES[1000]; // contains the name of the file with all the PidTable



  // ----------------------------------------------------------------------
  // -- The rest is provided by ROOT. No need to edit. 
  // ----------------------------------------------------------------------

  b2uNtp(TTree *tree, int isMC);
  ~b2uNtp();
  Int_t  Cut(Int_t entry);
  Int_t  GetEntry(Int_t entry);
  Int_t  LoadTree(Int_t entry);
  void   Init(TTree *tree);
  void   Loop(int maxEvent = 0, int startEvent = 0); 
  Bool_t Notify();
  void   Show(Int_t entry = -1);


  // -- Tree definition 
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  //Declaration of leaves types
  Int_t           event;
  Int_t           runNumber;
  Int_t           platform;
  Int_t           partition;
  Int_t           upperID;
  Int_t           lowerID;
  Float_t         primVtxX;
  Float_t         primVtxY;
  Float_t         primVtxZ;
  Float_t         primVtxCovXX;
  Float_t         primVtxCovYY;
  Float_t         primVtxCovZZ;
  Float_t         primVtxCovXY;
  Float_t         primVtxCovYZ;
  Float_t         primVtxCovXZ;
  Float_t         primVtxChi2;
  Int_t           primVtxNdof;
  UChar_t         BCountFilter;
  UChar_t         DchTrig;
  UChar_t         EmcTrig;
  Float_t         R2All;
  Int_t           nGTLFid;
  Int_t           nChgFid;
  Float_t         eTotFid;
  Float_t         PrimVtxdr;
  Float_t         PrimVtxdz;
  Float_t         VtxProb;
  Float_t         beamSX;
  Float_t         beamSY;
  Float_t         beamSZ;
  Float_t         beamSCovXX;
  Float_t         beamSCovYY;
  Float_t         beamSCovZZ;
  Float_t         beamSCovXZ;
  Float_t         pxUps;
  Float_t         pyUps;
  Float_t         pzUps;
  Float_t         eUps;
  Int_t           nTrkTot;
  Float_t         W2;
  Float_t         FoxWol2;
  Float_t         FoxWol2Neu;
  Float_t         thrust;
  Float_t         thrustNeu;
  Int_t           nMc;
  Float_t         pMc[200];   //[nMc]
  Float_t         massMc[200];   //[nMc]
  Float_t         thetaMc[200];   //[nMc]
  Float_t         phiMc[200];   //[nMc]
  Int_t           idMc[200];   //[nMc]
  Int_t           mothMc[200];   //[nMc]
  UInt_t          nDauMc[200];   //[nMc]
  Float_t         xMc[200];   //[nMc]
  Float_t         yMc[200];   //[nMc]
  Float_t         zMc[200];   //[nMc]
  Int_t           nB0;
  Float_t         massB0[100];   //[nB0]
  Float_t         pB0[100];   //[nB0]
  Float_t         thB0[100];   //[nB0]
  Float_t         phiB0[100];   //[nB0]
  Float_t         errMassB0[100];   //[nB0]
  Float_t         m0B0[100];   //[nB0]
  Float_t         xB0[100];   //[nB0]
  Float_t         yB0[100];   //[nB0]
  Float_t         zB0[100];   //[nB0]
  Float_t         s2xB0[100];   //[nB0]
  Float_t         s2yB0[100];   //[nB0]
  Float_t         s2zB0[100];   //[nB0]
  Float_t         chi2B0[100];   //[nB0]
  Int_t           dofB0[100];   //[nB0]
  Int_t           stB0[100];   //[nB0]
  Int_t           ndauB0[100];   //[nB0]
  Int_t           MCB0[100];   //[nB0]
  Float_t         mseB0[100];   //[nB0]
  Float_t         mHatB0[100];   //[nB0]
  Float_t         deltaeB0[100];   //[nB0]
  Float_t         ThruB0[100];   //[nB0]
  Float_t         thThruB0[100];   //[nB0]
  Float_t         phiThruB0[100];   //[nB0]
  Float_t         cosTBB0[100];   //[nB0]
  Int_t           d1B0Index[100];   //[nB0]
  Int_t           d1B0Lund[100];   //[nB0]
  Int_t           d2B0Index[100];   //[nB0]
  Int_t           d2B0Lund[100];   //[nB0]
  Int_t           d3B0Index[100];   //[nB0]
  Int_t           d3B0Lund[100];   //[nB0]
  Int_t           d4B0Index[100];   //[nB0]
  Int_t           d4B0Lund[100];   //[nB0]
  Int_t           d5B0Index[100];   //[nB0]
  Int_t           d5B0Lund[100];   //[nB0]
  Int_t           d6B0Index[100];   //[nB0]
  Int_t           d6B0Lund[100];   //[nB0]
  Int_t           d7B0Index[100];   //[nB0]
  Int_t           d7B0Lund[100];   //[nB0]
  Int_t           modeB0[100];   //[nB0]
  Float_t         purB0[100];   //[nB0]
  Float_t         intpurB0[100];   //[nB0]
  Float_t         VtxXLepB0[100];   //[nB0]
  Float_t         VtxYLepB0[100];   //[nB0]
  Float_t         VtxZLepB0[100];   //[nB0]
  Float_t         VtxCovXXLepB0[100];   //[nB0]
  Float_t         VtxCovYYLepB0[100];   //[nB0]
  Float_t         VtxCovXYLepB0[100];   //[nB0]
  Float_t         VtxCovZZLepB0[100];   //[nB0]
  Float_t         VtxCovXZLepB0[100];   //[nB0]
  Float_t         VtxCovYZLepB0[100];   //[nB0]
  Float_t         VtxChiSqLepB0[100];   //[nB0]
  Float_t         VtxNDofLepB0[100];   //[nB0]
  Int_t           VtxStatLepB0[100];   //[nB0]
  Int_t           VtxNUsedLepB0[100];   //[nB0]
  Float_t         DocaLepB0[100];   //[nB0]
  Float_t         DocaErrLepB0[100];   //[nB0]
  Float_t         VtxXXB0[100];   //[nB0]
  Float_t         VtxYXB0[100];   //[nB0]
  Float_t         VtxZXB0[100];   //[nB0]
  Float_t         VtxCovXXXB0[100];   //[nB0]
  Float_t         VtxCovYYXB0[100];   //[nB0]
  Float_t         VtxCovXYXB0[100];   //[nB0]
  Float_t         VtxCovZZXB0[100];   //[nB0]
  Float_t         VtxCovXZXB0[100];   //[nB0]
  Float_t         VtxCovYZXB0[100];   //[nB0]
  Float_t         VtxChiSqXB0[100];   //[nB0]
  Float_t         VtxNDofXB0[100];   //[nB0]
  Int_t           VtxStatXB0[100];   //[nB0]
  Int_t           VtxNUsedXB0[100];   //[nB0]
  Float_t         VtxPXB0[100];   //[nB0]
  Float_t         VtxPhiXB0[100];   //[nB0]
  Float_t         VtxThetaXB0[100];   //[nB0]
  Float_t         ThrustXB0[100];   //[nB0]
  Float_t         ThrustXPhiB0[100];   //[nB0]
  Float_t         ThrustXThetaB0[100];   //[nB0]
  Float_t         MassPB0[100];   //[nB0]
  Float_t         MassPhiB0[100];   //[nB0]
  Float_t         MassThetaB0[100];   //[nB0]
  Float_t         Cov00B0[100];   //[nB0]
  Float_t         Cov10B0[100];   //[nB0]
  Float_t         Cov11B0[100];   //[nB0]
  Float_t         Cov20B0[100];   //[nB0]
  Float_t         Cov21B0[100];   //[nB0]
  Float_t         Cov22B0[100];   //[nB0]
  Float_t         Cov30B0[100];   //[nB0]
  Float_t         Cov31B0[100];   //[nB0]
  Float_t         Cov32B0[100];   //[nB0]
  Float_t         Cov33B0[100];   //[nB0]
  Int_t           nChB;
  Float_t         massChB[100];   //[nChB]
  Float_t         pChB[100];   //[nChB]
  Float_t         thChB[100];   //[nChB]
  Float_t         phiChB[100];   //[nChB]
  Float_t         errMassChB[100];   //[nChB]
  Float_t         m0ChB[100];   //[nChB]
  Float_t         xChB[100];   //[nChB]
  Float_t         yChB[100];   //[nChB]
  Float_t         zChB[100];   //[nChB]
  Float_t         s2xChB[100];   //[nChB]
  Float_t         s2yChB[100];   //[nChB]
  Float_t         s2zChB[100];   //[nChB]
  Float_t         chi2ChB[100];   //[nChB]
  Int_t           dofChB[100];   //[nChB]
  Int_t           stChB[100];   //[nChB]
  Int_t           ndauChB[100];   //[nChB]
  Int_t           MCChB[100];   //[nChB]
  Float_t         mseChB[100];   //[nChB]
  Float_t         mHatChB[100];   //[nChB]
  Float_t         deltaeChB[100];   //[nChB]
  Float_t         ThruChB[100];   //[nChB]
  Float_t         thThruChB[100];   //[nChB]
  Float_t         phiThruChB[100];   //[nChB]
  Float_t         cosTBChB[100];   //[nChB]
  Int_t           d1ChBIndex[100];   //[nChB]
  Int_t           d1ChBLund[100];   //[nChB]
  Int_t           d2ChBIndex[100];   //[nChB]
  Int_t           d2ChBLund[100];   //[nChB]
  Int_t           d3ChBIndex[100];   //[nChB]
  Int_t           d3ChBLund[100];   //[nChB]
  Int_t           d4ChBIndex[100];   //[nChB]
  Int_t           d4ChBLund[100];   //[nChB]
  Int_t           d5ChBIndex[100];   //[nChB]
  Int_t           d5ChBLund[100];   //[nChB]
  Int_t           d6ChBIndex[100];   //[nChB]
  Int_t           d6ChBLund[100];   //[nChB]
  Int_t           d7ChBIndex[100];   //[nChB]
  Int_t           d7ChBLund[100];   //[nChB]
  Int_t           modeChB[100];   //[nChB]
  Float_t         purChB[100];   //[nChB]
  Float_t         intpurChB[100];   //[nChB]
  Float_t         VtxXLepChB[100];   //[nChB]
  Float_t         VtxYLepChB[100];   //[nChB]
  Float_t         VtxZLepChB[100];   //[nChB]
  Float_t         VtxCovXXLepChB[100];   //[nChB]
  Float_t         VtxCovYYLepChB[100];   //[nChB]
  Float_t         VtxCovXYLepChB[100];   //[nChB]
  Float_t         VtxCovZZLepChB[100];   //[nChB]
  Float_t         VtxCovXZLepChB[100];   //[nChB]
  Float_t         VtxCovYZLepChB[100];   //[nChB]
  Float_t         VtxChiSqLepChB[100];   //[nChB]
  Float_t         VtxNDofLepChB[100];   //[nChB]
  Int_t           VtxStatLepChB[100];   //[nChB]
  Int_t           VtxNUsedLepChB[100];   //[nChB]
  Float_t         DocaLepChB[100];   //[nChB]
  Float_t         DocaErrLepChB[100];   //[nChB]
  Float_t         VtxXXChB[100];   //[nChB]
  Float_t         VtxYXChB[100];   //[nChB]
  Float_t         VtxZXChB[100];   //[nChB]
  Float_t         VtxCovXXXChB[100];   //[nChB]
  Float_t         VtxCovYYXChB[100];   //[nChB]
  Float_t         VtxCovXYXChB[100];   //[nChB]
  Float_t         VtxCovZZXChB[100];   //[nChB]
  Float_t         VtxCovXZXChB[100];   //[nChB]
  Float_t         VtxCovYZXChB[100];   //[nChB]
  Float_t         VtxChiSqXChB[100];   //[nChB]
  Float_t         VtxNDofXChB[100];   //[nChB]
  Int_t           VtxStatXChB[100];   //[nChB]
  Int_t           VtxNUsedXChB[100];   //[nChB]
  Float_t         VtxPXChB[100];   //[nChB]
  Float_t         VtxPhiXChB[100];   //[nChB]
  Float_t         VtxThetaXChB[100];   //[nChB]
  Float_t         ThrustXChB[100];   //[nChB]
  Float_t         ThrustXPhiChB[100];   //[nChB]
  Float_t         ThrustXThetaChB[100];   //[nChB]
  Float_t         MassPChB[100];   //[nChB]
  Float_t         MassPhiChB[100];   //[nChB]
  Float_t         MassThetaChB[100];   //[nChB]
  Float_t         Cov00ChB[100];   //[nChB]
  Float_t         Cov10ChB[100];   //[nChB]
  Float_t         Cov11ChB[100];   //[nChB]
  Float_t         Cov20ChB[100];   //[nChB]
  Float_t         Cov21ChB[100];   //[nChB]
  Float_t         Cov22ChB[100];   //[nChB]
  Float_t         Cov30ChB[100];   //[nChB]
  Float_t         Cov31ChB[100];   //[nChB]
  Float_t         Cov32ChB[100];   //[nChB]
  Float_t         Cov33ChB[100];   //[nChB]
  Int_t           nDstar;
  Float_t         massDstar[100];   //[nDstar]
  Float_t         pDstar[100];   //[nDstar]
  Float_t         thDstar[100];   //[nDstar]
  Float_t         phiDstar[100];   //[nDstar]
  Float_t         errMassDstar[100];   //[nDstar]
  Float_t         m0Dstar[100];   //[nDstar]
  Float_t         xDstar[100];   //[nDstar]
  Float_t         yDstar[100];   //[nDstar]
  Float_t         zDstar[100];   //[nDstar]
  Float_t         s2xDstar[100];   //[nDstar]
  Float_t         s2yDstar[100];   //[nDstar]
  Float_t         s2zDstar[100];   //[nDstar]
  Float_t         chi2Dstar[100];   //[nDstar]
  Int_t           dofDstar[100];   //[nDstar]
  Int_t           stDstar[100];   //[nDstar]
  Int_t           ndauDstar[100];   //[nDstar]
  Int_t           MCDstar[100];   //[nDstar]
  Int_t           d1DstarIndex[100];   //[nDstar]
  Int_t           d1DstarLund[100];   //[nDstar]
  Int_t           d2DstarIndex[100];   //[nDstar]
  Int_t           d2DstarLund[100];   //[nDstar]
  Int_t           nDstarBS;
  Float_t         massDstarBS[100];   //[nDstarBS]
  Float_t         chi2DstarBS[100];   //[nDstarBS]
  Int_t           dofDstarBS[100];   //[nDstarBS]
  Float_t         spixDstarBS[100];   //[nDstarBS]
  Float_t         spiyDstarBS[100];   //[nDstarBS]
  Float_t         spizDstarBS[100];   //[nDstarBS]
  Int_t           nDstar0;
  Float_t         massDstar0[100];   //[nDstar0]
  Float_t         pDstar0[100];   //[nDstar0]
  Float_t         thDstar0[100];   //[nDstar0]
  Float_t         phiDstar0[100];   //[nDstar0]
  Float_t         errMassDstar0[100];   //[nDstar0]
  Float_t         m0Dstar0[100];   //[nDstar0]
  Float_t         xDstar0[100];   //[nDstar0]
  Float_t         yDstar0[100];   //[nDstar0]
  Float_t         zDstar0[100];   //[nDstar0]
  Float_t         s2xDstar0[100];   //[nDstar0]
  Float_t         s2yDstar0[100];   //[nDstar0]
  Float_t         s2zDstar0[100];   //[nDstar0]
  Float_t         chi2Dstar0[100];   //[nDstar0]
  Int_t           dofDstar0[100];   //[nDstar0]
  Int_t           stDstar0[100];   //[nDstar0]
  Int_t           ndauDstar0[100];   //[nDstar0]
  Int_t           MCDstar0[100];   //[nDstar0]
  Int_t           d1Dstar0Index[100];   //[nDstar0]
  Int_t           d1Dstar0Lund[100];   //[nDstar0]
  Int_t           d2Dstar0Index[100];   //[nDstar0]
  Int_t           d2Dstar0Lund[100];   //[nDstar0]
  Int_t           nD0;
  Float_t         massD0[100];   //[nD0]
  Float_t         pD0[100];   //[nD0]
  Float_t         thD0[100];   //[nD0]
  Float_t         phiD0[100];   //[nD0]
  Float_t         errMassD0[100];   //[nD0]
  Float_t         m0D0[100];   //[nD0]
  Float_t         xD0[100];   //[nD0]
  Float_t         yD0[100];   //[nD0]
  Float_t         zD0[100];   //[nD0]
  Float_t         s2xD0[100];   //[nD0]
  Float_t         s2yD0[100];   //[nD0]
  Float_t         s2zD0[100];   //[nD0]
  Float_t         chi2D0[100];   //[nD0]
  Int_t           dofD0[100];   //[nD0]
  Int_t           stD0[100];   //[nD0]
  Int_t           ndauD0[100];   //[nD0]
  Int_t           MCD0[100];   //[nD0]
  Int_t           d1D0Index[100];   //[nD0]
  Int_t           d1D0Lund[100];   //[nD0]
  Int_t           d2D0Index[100];   //[nD0]
  Int_t           d2D0Lund[100];   //[nD0]
  Int_t           d3D0Index[100];   //[nD0]
  Int_t           d3D0Lund[100];   //[nD0]
  Int_t           d4D0Index[100];   //[nD0]
  Int_t           d4D0Lund[100];   //[nD0]
  Int_t           nChD;
  Float_t         massChD[100];   //[nChD]
  Float_t         pChD[100];   //[nChD]
  Float_t         thChD[100];   //[nChD]
  Float_t         phiChD[100];   //[nChD]
  Float_t         errMassChD[100];   //[nChD]
  Float_t         m0ChD[100];   //[nChD]
  Float_t         xChD[100];   //[nChD]
  Float_t         yChD[100];   //[nChD]
  Float_t         zChD[100];   //[nChD]
  Float_t         s2xChD[100];   //[nChD]
  Float_t         s2yChD[100];   //[nChD]
  Float_t         s2zChD[100];   //[nChD]
  Float_t         chi2ChD[100];   //[nChD]
  Int_t           dofChD[100];   //[nChD]
  Int_t           stChD[100];   //[nChD]
  Int_t           ndauChD[100];   //[nChD]
  Int_t           MCChD[100];   //[nChD]
  Int_t           d1ChDIndex[100];   //[nChD]
  Int_t           d1ChDLund[100];   //[nChD]
  Int_t           d2ChDIndex[100];   //[nChD]
  Int_t           d2ChDLund[100];   //[nChD]
  Int_t           d3ChDIndex[100];   //[nChD]
  Int_t           d3ChDLund[100];   //[nChD]
  Int_t           d4ChDIndex[100];   //[nChD]
  Int_t           d4ChDLund[100];   //[nChD]
  Int_t           nJpsi;
  Float_t         massJpsi[50];   //[nJpsi]
  Float_t         pJpsi[50];   //[nJpsi]
  Float_t         thJpsi[50];   //[nJpsi]
  Float_t         phiJpsi[50];   //[nJpsi]
  Float_t         errMassJpsi[50];   //[nJpsi]
  Float_t         m0Jpsi[50];   //[nJpsi]
  Float_t         xJpsi[50];   //[nJpsi]
  Float_t         yJpsi[50];   //[nJpsi]
  Float_t         zJpsi[50];   //[nJpsi]
  Float_t         s2xJpsi[50];   //[nJpsi]
  Float_t         s2yJpsi[50];   //[nJpsi]
  Float_t         s2zJpsi[50];   //[nJpsi]
  Float_t         chi2Jpsi[50];   //[nJpsi]
  Int_t           dofJpsi[50];   //[nJpsi]
  Int_t           stJpsi[50];   //[nJpsi]
  Int_t           ndauJpsi[50];   //[nJpsi]
  Int_t           MCJpsi[50];   //[nJpsi]
  Int_t           d1JpsiIndex[50];   //[nJpsi]
  Int_t           d1JpsiLund[50];   //[nJpsi]
  Int_t           d1JpsiGamIndex[50];   //[nJpsi]
  Int_t           d1JpsiGamBrIndex[50];   //[nJpsi]
  Int_t           d1JpsiGamNumBr[50];   //[nJpsi]
  Int_t           d2JpsiIndex[50];   //[nJpsi]
  Int_t           d2JpsiLund[50];   //[nJpsi]
  Int_t           d2JpsiGamIndex[50];   //[nJpsi]
  Int_t           d2JpsiGamBrIndex[50];   //[nJpsi]
  Int_t           d2JpsiGamNumBr[50];   //[nJpsi]
  Int_t           nKs;
  Float_t         massKs[50];   //[nKs]
  Float_t         pKs[50];   //[nKs]
  Float_t         thKs[50];   //[nKs]
  Float_t         phiKs[50];   //[nKs]
  Float_t         errMassKs[50];   //[nKs]
  Float_t         m0Ks[50];   //[nKs]
  Float_t         xKs[50];   //[nKs]
  Float_t         yKs[50];   //[nKs]
  Float_t         zKs[50];   //[nKs]
  Float_t         s2xKs[50];   //[nKs]
  Float_t         s2yKs[50];   //[nKs]
  Float_t         s2zKs[50];   //[nKs]
  Float_t         chi2Ks[50];   //[nKs]
  Int_t           dofKs[50];   //[nKs]
  Int_t           stKs[50];   //[nKs]
  Int_t           ndauKs[50];   //[nKs]
  Int_t           MCKs[50];   //[nKs]
  Int_t           d1KsIndex[50];   //[nKs]
  Int_t           d1KsLund[50];   //[nKs]
  Int_t           d2KsIndex[50];   //[nKs]
  Int_t           d2KsLund[50];   //[nKs]
  Int_t           nPi0;
  Float_t         massPi0[150];   //[nPi0]
  Float_t         pPi0[150];   //[nPi0]
  Float_t         thPi0[150];   //[nPi0]
  Float_t         phiPi0[150];   //[nPi0]
  Float_t         errMassPi0[150];   //[nPi0]
  Float_t         m0Pi0[150];   //[nPi0]
  Float_t         xPi0[150];   //[nPi0]
  Float_t         yPi0[150];   //[nPi0]
  Float_t         zPi0[150];   //[nPi0]
  Float_t         s2xPi0[150];   //[nPi0]
  Float_t         s2yPi0[150];   //[nPi0]
  Float_t         s2zPi0[150];   //[nPi0]
  Float_t         chi2Pi0[150];   //[nPi0]
  Int_t           dofPi0[150];   //[nPi0]
  Int_t           stPi0[150];   //[nPi0]
  Int_t           ndauPi0[150];   //[nPi0]
  Int_t           MCPi0[150];   //[nPi0]
  Int_t           d1Pi0Index[150];   //[nPi0]
  Int_t           d1Pi0Lund[150];   //[nPi0]
  Int_t           d2Pi0Index[150];   //[nPi0]
  Int_t           d2Pi0Lund[150];   //[nPi0]
  Int_t           nGConv;
  Float_t         massGConv[10];   //[nGConv]
  Float_t         pGConv[10];   //[nGConv]
  Float_t         thGConv[10];   //[nGConv]
  Float_t         phiGConv[10];   //[nGConv]
  Float_t         errMassGConv[10];   //[nGConv]
  Float_t         m0GConv[10];   //[nGConv]
  Float_t         xGConv[10];   //[nGConv]
  Float_t         yGConv[10];   //[nGConv]
  Float_t         zGConv[10];   //[nGConv]
  Float_t         s2xGConv[10];   //[nGConv]
  Float_t         s2yGConv[10];   //[nGConv]
  Float_t         s2zGConv[10];   //[nGConv]
  Float_t         chi2GConv[10];   //[nGConv]
  Int_t           dofGConv[10];   //[nGConv]
  Int_t           stGConv[10];   //[nGConv]
  Int_t           ndauGConv[10];   //[nGConv]
  Int_t           MCGConv[10];   //[nGConv]
  Int_t           d1GConvIndex[10];   //[nGConv]
  Int_t           d1GConvLund[10];   //[nGConv]
  Int_t           d2GConvIndex[10];   //[nGConv]
  Int_t           d2GConvLund[10];   //[nGConv]
  Int_t           nDalitz;
  Float_t         massDalitz[10];   //[nDalitz]
  Float_t         pDalitz[10];   //[nDalitz]
  Float_t         thDalitz[10];   //[nDalitz]
  Float_t         phiDalitz[10];   //[nDalitz]
  Float_t         errMassDalitz[10];   //[nDalitz]
  Float_t         m0Dalitz[10];   //[nDalitz]
  Float_t         xDalitz[10];   //[nDalitz]
  Float_t         yDalitz[10];   //[nDalitz]
  Float_t         zDalitz[10];   //[nDalitz]
  Float_t         s2xDalitz[10];   //[nDalitz]
  Float_t         s2yDalitz[10];   //[nDalitz]
  Float_t         s2zDalitz[10];   //[nDalitz]
  Float_t         chi2Dalitz[10];   //[nDalitz]
  Int_t           dofDalitz[10];   //[nDalitz]
  Int_t           stDalitz[10];   //[nDalitz]
  Int_t           ndauDalitz[10];   //[nDalitz]
  Int_t           MCDalitz[10];   //[nDalitz]
  Int_t           d1DalitzIndex[10];   //[nDalitz]
  Int_t           d1DalitzLund[10];   //[nDalitz]
  Int_t           d2DalitzIndex[10];   //[nDalitz]
  Int_t           d2DalitzLund[10];   //[nDalitz]
  Int_t           nTrk;
  Int_t           IfrLayTrk[60];   //[nTrk]
  Int_t           IfrNsTrk[60];   //[nTrk]
  UChar_t         IfrInnerTrk[60];   //[nTrk]
  UChar_t         IfrBarrelTrk[60];   //[nTrk]
  UChar_t         IfrFWDTrk[60];   //[nTrk]
  UChar_t         IfrBWDTrk[60];   //[nTrk]
  Float_t         IfrMeasIntLenTrk[60];   //[nTrk]
  Int_t           IfrFirstHitTrk[60];   //[nTrk]
  Int_t           IfrLastHitTrk[60];   //[nTrk]
  Float_t         lMomTrk[60];   //[nTrk]
  Float_t         ZMom42Trk[60];   //[nTrk]
  Float_t         ecalTrk[60];   //[nTrk]
  Float_t         ecalXTrk[60];   //[nTrk]
  Float_t         ecalYTrk[60];   //[nTrk]
  Float_t         ecalZTrk[60];   //[nTrk]
  Int_t           nCryTrk[60];   //[nTrk]
  Int_t           nBumpTrk[60];   //[nTrk]
  Float_t         ZMom20Trk[60];   //[nTrk]
  Float_t         secMomTrk[60];   //[nTrk]
  Float_t         s1s9Trk[60];   //[nTrk]
  Float_t         s9s25Trk[60];   //[nTrk]
  Float_t         erawTrk[60];   //[nTrk]
  Float_t         phiClusterTrk[60];   //[nTrk]
  Float_t         thetaClusterTrk[60];   //[nTrk]
  Float_t         covEETrk[60];   //[nTrk]
  Float_t         covTTTrk[60];   //[nTrk]
  Float_t         covPPTrk[60];   //[nTrk]
  Float_t         covRRTrk[60];   //[nTrk]
  Float_t         phicMatTrk[60];   //[nTrk]
  Float_t         trkcMatTrk[60];   //[nTrk]
  Int_t           nPidTrk[60];   //[nTrk]
  Int_t           emcStatusTrk[60];   //[nTrk]
  Float_t         phiAtEMCTrk[60];   //[nTrk]
  Float_t         thetaAtEMCTrk[60];   //[nTrk]
  Int_t           isvtTrk[60];   //[nTrk]
  Int_t           nsvtTrk[60];   //[nTrk]
  Int_t           fhitTrk[60];   //[nTrk]
  Int_t           ndchTrk[60];   //[nTrk]
  Int_t           lhitTrk[60];   //[nTrk]
  Float_t         tLenTrk[60];   //[nTrk]
  Int_t           ntdofTrk[60];   //[nTrk]
  Float_t         tproTrk[60];   //[nTrk]
  Float_t         tChi2Trk[60];   //[nTrk]
  Int_t           cPidTrk[60];   //[nTrk]
  Float_t         sfRangeTrk[60];   //[nTrk]
  UChar_t         trkFitStatusTrk[60];   //[nTrk]
  Int_t           chargeTrk[60];   //[nTrk]
  Float_t         momentumTrk[60];   //[nTrk]
  Float_t         ppcov00[60];   //[nTrk]
  Float_t         ppcov10[60];   //[nTrk]
  Float_t         ppcov11[60];   //[nTrk]
  Float_t         ppcov20[60];   //[nTrk]
  Float_t         ppcov21[60];   //[nTrk]
  Float_t         ppcov22[60];   //[nTrk]
  Float_t         xPocaTrk[60];   //[nTrk]
  Float_t         yPocaTrk[60];   //[nTrk]
  Float_t         zPocaTrk[60];   //[nTrk]
  Float_t         thetaTrk[60];   //[nTrk]
  Float_t         phiTrk[60];   //[nTrk]
  Int_t           muonIdTrk[60];   //[nTrk]
  Int_t           elecIdTrk[60];   //[nTrk]
  Int_t           kaonIdTrk[60];   //[nTrk]
  Int_t           pionIdTrk[60];   //[nTrk]
  Int_t           protonIdTrk[60];   //[nTrk]
  Int_t           idTrk[60];   //[nTrk]
  Int_t           IndexTrk[60];   //[nTrk]
  Int_t           IndexNtTrk[60];   //[nTrk]
  Int_t           B0RecTrk[60];   //[nTrk]
  Int_t           chBRecTrk[60];   //[nTrk]
  Int_t           nGam;
  Int_t           IfrLayGam[80];   //[nGam]
  Int_t           IfrNsGam[80];   //[nGam]
  UChar_t         IfrInnerGam[80];   //[nGam]
  UChar_t         IfrBarrelGam[80];   //[nGam]
  UChar_t         IfrFWDGam[80];   //[nGam]
  UChar_t         IfrBWDGam[80];   //[nGam]
  Float_t         IfrMeasIntLenGam[80];   //[nGam]
  Int_t           IfrFirstHitGam[80];   //[nGam]
  Int_t           IfrLastHitGam[80];   //[nGam]
  Float_t         IfrExpIntLenGam[80];   //[nGam]
  Float_t         IfrIntLenBeforeIronGam[80];   //[nGam]
  Float_t         IfrTrkMatchGam[80];   //[nGam]
  Float_t         IfrEmcMatchGam[80];   //[nGam]
  Int_t           IfrLastBarrelGam[80];   //[nGam]
  Float_t         IfrCLFitChi2Gam[80];   //[nGam]
  Int_t           IfrStrips0[80];   //[nGam]
  Int_t           IfrStrips1[80];   //[nGam]
  Int_t           IfrStrips2[80];   //[nGam]
  Int_t           IfrStrips3[80];   //[nGam]
  Int_t           IfrStrips4[80];   //[nGam]
  Int_t           IfrStrips5[80];   //[nGam]
  Int_t           IfrStrips6[80];   //[nGam]
  Int_t           IfrStrips7[80];   //[nGam]
  Int_t           IfrStrips8[80];   //[nGam]
  Int_t           IfrStrips9[80];   //[nGam]
  Int_t           IfrStrips10[80];   //[nGam]
  Int_t           IfrStrips11[80];   //[nGam]
  Int_t           IfrStrips12[80];   //[nGam]
  Int_t           IfrStrips13[80];   //[nGam]
  Int_t           IfrStrips14[80];   //[nGam]
  Int_t           IfrStrips15[80];   //[nGam]
  Int_t           IfrStrips16[80];   //[nGam]
  Int_t           IfrStrips17[80];   //[nGam]
  Int_t           IfrStrips18[80];   //[nGam]
  Int_t           IfrStrips19[80];   //[nGam]
  Float_t         lMomGam[80];   //[nGam]
  Float_t         ZMom42Gam[80];   //[nGam]
  Float_t         ecalGam[80];   //[nGam]
  Float_t         ecalXGam[80];   //[nGam]
  Float_t         ecalYGam[80];   //[nGam]
  Float_t         ecalZGam[80];   //[nGam]
  Int_t           nCryGam[80];   //[nGam]
  Int_t           nBumpGam[80];   //[nGam]
  Float_t         ZMom20Gam[80];   //[nGam]
  Float_t         secMomGam[80];   //[nGam]
  Float_t         s1s9Gam[80];   //[nGam]
  Float_t         s9s25Gam[80];   //[nGam]
  Float_t         erawGam[80];   //[nGam]
  Float_t         phiClusterGam[80];   //[nGam]
  Float_t         thetaClusterGam[80];   //[nGam]
  Float_t         covEEGam[80];   //[nGam]
  Float_t         covTTGam[80];   //[nGam]
  Float_t         covPPGam[80];   //[nGam]
  Float_t         covRRGam[80];   //[nGam]
  Int_t           emcStatusGam[80];   //[nGam]
  Float_t         thetaGam[80];   //[nGam]
  Float_t         phiGam[80];   //[nGam]
  Float_t         energyGam[80];   //[nGam]
  Int_t           idGam[80];   //[nGam]
  Int_t           IndexGam[80];   //[nGam]
  Int_t           IndexNtGam[80];   //[nGam]
  Int_t           B0RecGam[80];   //[nGam]
  Int_t           chBRecGam[80];   //[nGam]

  //List of branches
  TBranch        *b_event;   //!
  TBranch        *b_runNumber;   //!
  TBranch        *b_platform;   //!
  TBranch        *b_partition;   //!
  TBranch        *b_upperID;   //!
  TBranch        *b_lowerID;   //!
  TBranch        *b_primVtxX;   //!
  TBranch        *b_primVtxY;   //!
  TBranch        *b_primVtxZ;   //!
  TBranch        *b_primVtxCovXX;   //!
  TBranch        *b_primVtxCovYY;   //!
  TBranch        *b_primVtxCovZZ;   //!
  TBranch        *b_primVtxCovXY;   //!
  TBranch        *b_primVtxCovYZ;   //!
  TBranch        *b_primVtxCovXZ;   //!
  TBranch        *b_primVtxChi2;   //!
  TBranch        *b_primVtxNdof;   //!
  TBranch        *b_BCountFilter;   //!
  TBranch        *b_DchTrig;   //!
  TBranch        *b_EmcTrig;   //!
  TBranch        *b_R2All;   //!
  TBranch        *b_nGTLFid;   //!
  TBranch        *b_nChgFid;   //!
  TBranch        *b_eTotFid;   //!
  TBranch        *b_PrimVtxdr;   //!
  TBranch        *b_PrimVtxdz;   //!
  TBranch        *b_VtxProb;   //!
  TBranch        *b_beamSX;   //!
  TBranch        *b_beamSY;   //!
  TBranch        *b_beamSZ;   //!
  TBranch        *b_beamSCovXX;   //!
  TBranch        *b_beamSCovYY;   //!
  TBranch        *b_beamSCovZZ;   //!
  TBranch        *b_beamSCovXZ;   //!
  TBranch        *b_pxUps;   //!
  TBranch        *b_pyUps;   //!
  TBranch        *b_pzUps;   //!
  TBranch        *b_eUps;   //!
  TBranch        *b_nTrkTot;   //!
  TBranch        *b_W2;   //!
  TBranch        *b_FoxWol2;   //!
  TBranch        *b_FoxWol2Neu;   //!
  TBranch        *b_thrust;   //!
  TBranch        *b_thrustNeu;   //!
  TBranch        *b_nMc;   //!
  TBranch        *b_pMc;   //!
  TBranch        *b_massMc;   //!
  TBranch        *b_thetaMc;   //!
  TBranch        *b_phiMc;   //!
  TBranch        *b_idMc;   //!
  TBranch        *b_mothMc;   //!
  TBranch        *b_nDauMc;   //!
  TBranch        *b_xMc;   //!
  TBranch        *b_yMc;   //!
  TBranch        *b_zMc;   //!
  TBranch        *b_nB0;   //!
  TBranch        *b_massB0;   //!
  TBranch        *b_pB0;   //!
  TBranch        *b_thB0;   //!
  TBranch        *b_phiB0;   //!
  TBranch        *b_errMassB0;   //!
  TBranch        *b_m0B0;   //!
  TBranch        *b_xB0;   //!
  TBranch        *b_yB0;   //!
  TBranch        *b_zB0;   //!
  TBranch        *b_s2xB0;   //!
  TBranch        *b_s2yB0;   //!
  TBranch        *b_s2zB0;   //!
  TBranch        *b_chi2B0;   //!
  TBranch        *b_dofB0;   //!
  TBranch        *b_stB0;   //!
  TBranch        *b_ndauB0;   //!
  TBranch        *b_MCB0;   //!
  TBranch        *b_mseB0;   //!
  TBranch        *b_mHatB0;   //!
  TBranch        *b_deltaeB0;   //!
  TBranch        *b_ThruB0;   //!
  TBranch        *b_thThruB0;   //!
  TBranch        *b_phiThruB0;   //!
  TBranch        *b_cosTBB0;   //!
  TBranch        *b_d1B0Index;   //!
  TBranch        *b_d1B0Lund;   //!
  TBranch        *b_d2B0Index;   //!
  TBranch        *b_d2B0Lund;   //!
  TBranch        *b_d3B0Index;   //!
  TBranch        *b_d3B0Lund;   //!
  TBranch        *b_d4B0Index;   //!
  TBranch        *b_d4B0Lund;   //!
  TBranch        *b_d5B0Index;   //!
  TBranch        *b_d5B0Lund;   //!
  TBranch        *b_d6B0Index;   //!
  TBranch        *b_d6B0Lund;   //!
  TBranch        *b_d7B0Index;   //!
  TBranch        *b_d7B0Lund;   //!
  TBranch        *b_modeB0;   //!
  TBranch        *b_purB0;   //!
  TBranch        *b_intpurB0;   //!
  TBranch        *b_VtxXLepB0;   //!
  TBranch        *b_VtxYLepB0;   //!
  TBranch        *b_VtxZLepB0;   //!
  TBranch        *b_VtxCovXXLepB0;   //!
  TBranch        *b_VtxCovYYLepB0;   //!
  TBranch        *b_VtxCovXYLepB0;   //!
  TBranch        *b_VtxCovZZLepB0;   //!
  TBranch        *b_VtxCovXZLepB0;   //!
  TBranch        *b_VtxCovYZLepB0;   //!
  TBranch        *b_VtxChiSqLepB0;   //!
  TBranch        *b_VtxNDofLepB0;   //!
  TBranch        *b_VtxStatLepB0;   //!
  TBranch        *b_VtxNUsedLepB0;   //!
  TBranch        *b_DocaLepB0;   //!
  TBranch        *b_DocaErrLepB0;   //!
  TBranch        *b_VtxXXB0;   //!
  TBranch        *b_VtxYXB0;   //!
  TBranch        *b_VtxZXB0;   //!
  TBranch        *b_VtxCovXXXB0;   //!
  TBranch        *b_VtxCovYYXB0;   //!
  TBranch        *b_VtxCovXYXB0;   //!
  TBranch        *b_VtxCovZZXB0;   //!
  TBranch        *b_VtxCovXZXB0;   //!
  TBranch        *b_VtxCovYZXB0;   //!
  TBranch        *b_VtxChiSqXB0;   //!
  TBranch        *b_VtxNDofXB0;   //!
  TBranch        *b_VtxStatXB0;   //!
  TBranch        *b_VtxNUsedXB0;   //!
  TBranch        *b_VtxPXB0;   //!
  TBranch        *b_VtxPhiXB0;   //!
  TBranch        *b_VtxThetaXB0;   //!
  TBranch        *b_ThrustXB0;   //!
  TBranch        *b_ThrustXPhiB0;   //!
  TBranch        *b_ThrustXThetaB0;   //!
  TBranch        *b_MassPB0;   //!
  TBranch        *b_MassPhiB0;   //!
  TBranch        *b_MassThetaB0;   //!
  TBranch        *b_Cov00B0;   //!
  TBranch        *b_Cov10B0;   //!
  TBranch        *b_Cov11B0;   //!
  TBranch        *b_Cov20B0;   //!
  TBranch        *b_Cov21B0;   //!
  TBranch        *b_Cov22B0;   //!
  TBranch        *b_Cov30B0;   //!
  TBranch        *b_Cov31B0;   //!
  TBranch        *b_Cov32B0;   //!
  TBranch        *b_Cov33B0;   //!
  TBranch        *b_nChB;   //!
  TBranch        *b_massChB;   //!
  TBranch        *b_pChB;   //!
  TBranch        *b_thChB;   //!
  TBranch        *b_phiChB;   //!
  TBranch        *b_errMassChB;   //!
  TBranch        *b_m0ChB;   //!
  TBranch        *b_xChB;   //!
  TBranch        *b_yChB;   //!
  TBranch        *b_zChB;   //!
  TBranch        *b_s2xChB;   //!
  TBranch        *b_s2yChB;   //!
  TBranch        *b_s2zChB;   //!
  TBranch        *b_chi2ChB;   //!
  TBranch        *b_dofChB;   //!
  TBranch        *b_stChB;   //!
  TBranch        *b_ndauChB;   //!
  TBranch        *b_MCChB;   //!
  TBranch        *b_mseChB;   //!
  TBranch        *b_mHatChB;   //!
  TBranch        *b_deltaeChB;   //!
  TBranch        *b_ThruChB;   //!
  TBranch        *b_thThruChB;   //!
  TBranch        *b_phiThruChB;   //!
  TBranch        *b_cosTBChB;   //!
  TBranch        *b_d1ChBIndex;   //!
  TBranch        *b_d1ChBLund;   //!
  TBranch        *b_d2ChBIndex;   //!
  TBranch        *b_d2ChBLund;   //!
  TBranch        *b_d3ChBIndex;   //!
  TBranch        *b_d3ChBLund;   //!
  TBranch        *b_d4ChBIndex;   //!
  TBranch        *b_d4ChBLund;   //!
  TBranch        *b_d5ChBIndex;   //!
  TBranch        *b_d5ChBLund;   //!
  TBranch        *b_d6ChBIndex;   //!
  TBranch        *b_d6ChBLund;   //!
  TBranch        *b_d7ChBIndex;   //!
  TBranch        *b_d7ChBLund;   //!
  TBranch        *b_modeChB;   //!
  TBranch        *b_purChB;   //!
  TBranch        *b_intpurChB;   //!
  TBranch        *b_VtxXLepChB;   //!
  TBranch        *b_VtxYLepChB;   //!
  TBranch        *b_VtxZLepChB;   //!
  TBranch        *b_VtxCovXXLepChB;   //!
  TBranch        *b_VtxCovYYLepChB;   //!
  TBranch        *b_VtxCovXYLepChB;   //!
  TBranch        *b_VtxCovZZLepChB;   //!
  TBranch        *b_VtxCovXZLepChB;   //!
  TBranch        *b_VtxCovYZLepChB;   //!
  TBranch        *b_VtxChiSqLepChB;   //!
  TBranch        *b_VtxNDofLepChB;   //!
  TBranch        *b_VtxStatLepChB;   //!
  TBranch        *b_VtxNUsedLepChB;   //!
  TBranch        *b_DocaLepChB;   //!
  TBranch        *b_DocaErrLepChB;   //!
  TBranch        *b_VtxXXChB;   //!
  TBranch        *b_VtxYXChB;   //!
  TBranch        *b_VtxZXChB;   //!
  TBranch        *b_VtxCovXXXChB;   //!
  TBranch        *b_VtxCovYYXChB;   //!
  TBranch        *b_VtxCovXYXChB;   //!
  TBranch        *b_VtxCovZZXChB;   //!
  TBranch        *b_VtxCovXZXChB;   //!
  TBranch        *b_VtxCovYZXChB;   //!
  TBranch        *b_VtxChiSqXChB;   //!
  TBranch        *b_VtxNDofXChB;   //!
  TBranch        *b_VtxStatXChB;   //!
  TBranch        *b_VtxNUsedXChB;   //!
  TBranch        *b_VtxPXChB;   //!
  TBranch        *b_VtxPhiXChB;   //!
  TBranch        *b_VtxThetaXChB;   //!
  TBranch        *b_ThrustXChB;   //!
  TBranch        *b_ThrustXPhiChB;   //!
  TBranch        *b_ThrustXThetaChB;   //!
  TBranch        *b_MassPChB;   //!
  TBranch        *b_MassPhiChB;   //!
  TBranch        *b_MassThetaChB;   //!
  TBranch        *b_Cov00ChB;   //!
  TBranch        *b_Cov10ChB;   //!
  TBranch        *b_Cov11ChB;   //!
  TBranch        *b_Cov20ChB;   //!
  TBranch        *b_Cov21ChB;   //!
  TBranch        *b_Cov22ChB;   //!
  TBranch        *b_Cov30ChB;   //!
  TBranch        *b_Cov31ChB;   //!
  TBranch        *b_Cov32ChB;   //!
  TBranch        *b_Cov33ChB;   //!
  TBranch        *b_nDstar;   //!
  TBranch        *b_massDstar;   //!
  TBranch        *b_pDstar;   //!
  TBranch        *b_thDstar;   //!
  TBranch        *b_phiDstar;   //!
  TBranch        *b_errMassDstar;   //!
  TBranch        *b_m0Dstar;   //!
  TBranch        *b_xDstar;   //!
  TBranch        *b_yDstar;   //!
  TBranch        *b_zDstar;   //!
  TBranch        *b_s2xDstar;   //!
  TBranch        *b_s2yDstar;   //!
  TBranch        *b_s2zDstar;   //!
  TBranch        *b_chi2Dstar;   //!
  TBranch        *b_dofDstar;   //!
  TBranch        *b_stDstar;   //!
  TBranch        *b_ndauDstar;   //!
  TBranch        *b_MCDstar;   //!
  TBranch        *b_d1DstarIndex;   //!
  TBranch        *b_d1DstarLund;   //!
  TBranch        *b_d2DstarIndex;   //!
  TBranch        *b_d2DstarLund;   //!
  TBranch        *b_nDstarBS;   //!
  TBranch        *b_massDstarBS;   //!
  TBranch        *b_chi2DstarBS;   //!
  TBranch        *b_dofDstarBS;   //!
  TBranch        *b_spixDstarBS;   //!
  TBranch        *b_spiyDstarBS;   //!
  TBranch        *b_spizDstarBS;   //!
  TBranch        *b_nDstar0;   //!
  TBranch        *b_massDstar0;   //!
  TBranch        *b_pDstar0;   //!
  TBranch        *b_thDstar0;   //!
  TBranch        *b_phiDstar0;   //!
  TBranch        *b_errMassDstar0;   //!
  TBranch        *b_m0Dstar0;   //!
  TBranch        *b_xDstar0;   //!
  TBranch        *b_yDstar0;   //!
  TBranch        *b_zDstar0;   //!
  TBranch        *b_s2xDstar0;   //!
  TBranch        *b_s2yDstar0;   //!
  TBranch        *b_s2zDstar0;   //!
  TBranch        *b_chi2Dstar0;   //!
  TBranch        *b_dofDstar0;   //!
  TBranch        *b_stDstar0;   //!
  TBranch        *b_ndauDstar0;   //!
  TBranch        *b_MCDstar0;   //!
  TBranch        *b_d1Dstar0Index;   //!
  TBranch        *b_d1Dstar0Lund;   //!
  TBranch        *b_d2Dstar0Index;   //!
  TBranch        *b_d2Dstar0Lund;   //!
  TBranch        *b_nD0;   //!
  TBranch        *b_massD0;   //!
  TBranch        *b_pD0;   //!
  TBranch        *b_thD0;   //!
  TBranch        *b_phiD0;   //!
  TBranch        *b_errMassD0;   //!
  TBranch        *b_m0D0;   //!
  TBranch        *b_xD0;   //!
  TBranch        *b_yD0;   //!
  TBranch        *b_zD0;   //!
  TBranch        *b_s2xD0;   //!
  TBranch        *b_s2yD0;   //!
  TBranch        *b_s2zD0;   //!
  TBranch        *b_chi2D0;   //!
  TBranch        *b_dofD0;   //!
  TBranch        *b_stD0;   //!
  TBranch        *b_ndauD0;   //!
  TBranch        *b_MCD0;   //!
  TBranch        *b_d1D0Index;   //!
  TBranch        *b_d1D0Lund;   //!
  TBranch        *b_d2D0Index;   //!
  TBranch        *b_d2D0Lund;   //!
  TBranch        *b_d3D0Index;   //!
  TBranch        *b_d3D0Lund;   //!
  TBranch        *b_d4D0Index;   //!
  TBranch        *b_d4D0Lund;   //!
  TBranch        *b_nChD;   //!
  TBranch        *b_massChD;   //!
  TBranch        *b_pChD;   //!
  TBranch        *b_thChD;   //!
  TBranch        *b_phiChD;   //!
  TBranch        *b_errMassChD;   //!
  TBranch        *b_m0ChD;   //!
  TBranch        *b_xChD;   //!
  TBranch        *b_yChD;   //!
  TBranch        *b_zChD;   //!
  TBranch        *b_s2xChD;   //!
  TBranch        *b_s2yChD;   //!
  TBranch        *b_s2zChD;   //!
  TBranch        *b_chi2ChD;   //!
  TBranch        *b_dofChD;   //!
  TBranch        *b_stChD;   //!
  TBranch        *b_ndauChD;   //!
  TBranch        *b_MCChD;   //!
  TBranch        *b_d1ChDIndex;   //!
  TBranch        *b_d1ChDLund;   //!
  TBranch        *b_d2ChDIndex;   //!
  TBranch        *b_d2ChDLund;   //!
  TBranch        *b_d3ChDIndex;   //!
  TBranch        *b_d3ChDLund;   //!
  TBranch        *b_d4ChDIndex;   //!
  TBranch        *b_d4ChDLund;   //!
  TBranch        *b_nJpsi;   //!
  TBranch        *b_massJpsi;   //!
  TBranch        *b_pJpsi;   //!
  TBranch        *b_thJpsi;   //!
  TBranch        *b_phiJpsi;   //!
  TBranch        *b_errMassJpsi;   //!
  TBranch        *b_m0Jpsi;   //!
  TBranch        *b_xJpsi;   //!
  TBranch        *b_yJpsi;   //!
  TBranch        *b_zJpsi;   //!
  TBranch        *b_s2xJpsi;   //!
  TBranch        *b_s2yJpsi;   //!
  TBranch        *b_s2zJpsi;   //!
  TBranch        *b_chi2Jpsi;   //!
  TBranch        *b_dofJpsi;   //!
  TBranch        *b_stJpsi;   //!
  TBranch        *b_ndauJpsi;   //!
  TBranch        *b_MCJpsi;   //!
  TBranch        *b_d1JpsiIndex;   //!
  TBranch        *b_d1JpsiLund;   //!
  TBranch        *b_d1JpsiGamIndex;   //!
  TBranch        *b_d1JpsiGamBrIndex;   //!
  TBranch        *b_d1JpsiGamNumBr;   //!
  TBranch        *b_d2JpsiIndex;   //!
  TBranch        *b_d2JpsiLund;   //!
  TBranch        *b_d2JpsiGamIndex;   //!
  TBranch        *b_d2JpsiGamBrIndex;   //!
  TBranch        *b_d2JpsiGamNumBr;   //!
  TBranch        *b_nKs;   //!
  TBranch        *b_massKs;   //!
  TBranch        *b_pKs;   //!
  TBranch        *b_thKs;   //!
  TBranch        *b_phiKs;   //!
  TBranch        *b_errMassKs;   //!
  TBranch        *b_m0Ks;   //!
  TBranch        *b_xKs;   //!
  TBranch        *b_yKs;   //!
  TBranch        *b_zKs;   //!
  TBranch        *b_s2xKs;   //!
  TBranch        *b_s2yKs;   //!
  TBranch        *b_s2zKs;   //!
  TBranch        *b_chi2Ks;   //!
  TBranch        *b_dofKs;   //!
  TBranch        *b_stKs;   //!
  TBranch        *b_ndauKs;   //!
  TBranch        *b_MCKs;   //!
  TBranch        *b_d1KsIndex;   //!
  TBranch        *b_d1KsLund;   //!
  TBranch        *b_d2KsIndex;   //!
  TBranch        *b_d2KsLund;   //!
  TBranch        *b_nPi0;   //!
  TBranch        *b_massPi0;   //!
  TBranch        *b_pPi0;   //!
  TBranch        *b_thPi0;   //!
  TBranch        *b_phiPi0;   //!
  TBranch        *b_errMassPi0;   //!
  TBranch        *b_m0Pi0;   //!
  TBranch        *b_xPi0;   //!
  TBranch        *b_yPi0;   //!
  TBranch        *b_zPi0;   //!
  TBranch        *b_s2xPi0;   //!
  TBranch        *b_s2yPi0;   //!
  TBranch        *b_s2zPi0;   //!
  TBranch        *b_chi2Pi0;   //!
  TBranch        *b_dofPi0;   //!
  TBranch        *b_stPi0;   //!
  TBranch        *b_ndauPi0;   //!
  TBranch        *b_MCPi0;   //!
  TBranch        *b_d1Pi0Index;   //!
  TBranch        *b_d1Pi0Lund;   //!
  TBranch        *b_d2Pi0Index;   //!
  TBranch        *b_d2Pi0Lund;   //!
  TBranch        *b_nGConv;   //!
  TBranch        *b_massGConv;   //!
  TBranch        *b_pGConv;   //!
  TBranch        *b_thGConv;   //!
  TBranch        *b_phiGConv;   //!
  TBranch        *b_errMassGConv;   //!
  TBranch        *b_m0GConv;   //!
  TBranch        *b_xGConv;   //!
  TBranch        *b_yGConv;   //!
  TBranch        *b_zGConv;   //!
  TBranch        *b_s2xGConv;   //!
  TBranch        *b_s2yGConv;   //!
  TBranch        *b_s2zGConv;   //!
  TBranch        *b_chi2GConv;   //!
  TBranch        *b_dofGConv;   //!
  TBranch        *b_stGConv;   //!
  TBranch        *b_ndauGConv;   //!
  TBranch        *b_MCGConv;   //!
  TBranch        *b_d1GConvIndex;   //!
  TBranch        *b_d1GConvLund;   //!
  TBranch        *b_d2GConvIndex;   //!
  TBranch        *b_d2GConvLund;   //!
  TBranch        *b_nDalitz;   //!
  TBranch        *b_massDalitz;   //!
  TBranch        *b_pDalitz;   //!
  TBranch        *b_thDalitz;   //!
  TBranch        *b_phiDalitz;   //!
  TBranch        *b_errMassDalitz;   //!
  TBranch        *b_m0Dalitz;   //!
  TBranch        *b_xDalitz;   //!
  TBranch        *b_yDalitz;   //!
  TBranch        *b_zDalitz;   //!
  TBranch        *b_s2xDalitz;   //!
  TBranch        *b_s2yDalitz;   //!
  TBranch        *b_s2zDalitz;   //!
  TBranch        *b_chi2Dalitz;   //!
  TBranch        *b_dofDalitz;   //!
  TBranch        *b_stDalitz;   //!
  TBranch        *b_ndauDalitz;   //!
  TBranch        *b_MCDalitz;   //!
  TBranch        *b_d1DalitzIndex;   //!
  TBranch        *b_d1DalitzLund;   //!
  TBranch        *b_d2DalitzIndex;   //!
  TBranch        *b_d2DalitzLund;   //!
  TBranch        *b_nTrk;   //!
  TBranch        *b_IfrLayTrk;   //!
  TBranch        *b_IfrNsTrk;   //!
  TBranch        *b_IfrInnerTrk;   //!
  TBranch        *b_IfrBarrelTrk;   //!
  TBranch        *b_IfrFWDTrk;   //!
  TBranch        *b_IfrBWDTrk;   //!
  TBranch        *b_IfrMeasIntLenTrk;   //!
  TBranch        *b_IfrFirstHitTrk;   //!
  TBranch        *b_IfrLastHitTrk;   //!
  TBranch        *b_lMomTrk;   //!
  TBranch        *b_ZMom42Trk;   //!
  TBranch        *b_ecalTrk;   //!
  TBranch        *b_ecalXTrk;   //!
  TBranch        *b_ecalYTrk;   //!
  TBranch        *b_ecalZTrk;   //!
  TBranch        *b_nCryTrk;   //!
  TBranch        *b_nBumpTrk;   //!
  TBranch        *b_ZMom20Trk;   //!
  TBranch        *b_secMomTrk;   //!
  TBranch        *b_s1s9Trk;   //!
  TBranch        *b_s9s25Trk;   //!
  TBranch        *b_erawTrk;   //!
  TBranch        *b_phiClusterTrk;   //!
  TBranch        *b_thetaClusterTrk;   //!
  TBranch        *b_covEETrk;   //!
  TBranch        *b_covTTTrk;   //!
  TBranch        *b_covPPTrk;   //!
  TBranch        *b_covRRTrk;   //!
  TBranch        *b_phicMatTrk;   //!
  TBranch        *b_trkcMatTrk;   //!
  TBranch        *b_nPidTrk;   //!
  TBranch        *b_emcStatusTrk;   //!
  TBranch        *b_phiAtEMCTrk;   //!
  TBranch        *b_thetaAtEMCTrk;   //!
  TBranch        *b_isvtTrk;   //!
  TBranch        *b_nsvtTrk;   //!
  TBranch        *b_fhitTrk;   //!
  TBranch        *b_ndchTrk;   //!
  TBranch        *b_lhitTrk;   //!
  TBranch        *b_tLenTrk;   //!
  TBranch        *b_ntdofTrk;   //!
  TBranch        *b_tproTrk;   //!
  TBranch        *b_tChi2Trk;   //!
  TBranch        *b_cPidTrk;   //!
  TBranch        *b_sfRangeTrk;   //!
  TBranch        *b_trkFitStatusTrk;   //!
  TBranch        *b_chargeTrk;   //!
  TBranch        *b_momentumTrk;   //!
  TBranch        *b_ppcov00;   //!
  TBranch        *b_ppcov10;   //!
  TBranch        *b_ppcov11;   //!
  TBranch        *b_ppcov20;   //!
  TBranch        *b_ppcov21;   //!
  TBranch        *b_ppcov22;   //!
  TBranch        *b_xPocaTrk;   //!
  TBranch        *b_yPocaTrk;   //!
  TBranch        *b_zPocaTrk;   //!
  TBranch        *b_thetaTrk;   //!
  TBranch        *b_phiTrk;   //!
  TBranch        *b_muonIdTrk;   //!
  TBranch        *b_elecIdTrk;   //!
  TBranch        *b_kaonIdTrk;   //!
  TBranch        *b_pionIdTrk;   //!
  TBranch        *b_protonIdTrk;   //!
  TBranch        *b_idTrk;   //!
  TBranch        *b_IndexTrk;   //!
  TBranch        *b_IndexNtTrk;   //!
  TBranch        *b_B0RecTrk;   //!
  TBranch        *b_chBRecTrk;   //!
  TBranch        *b_nGam;   //!
  TBranch        *b_IfrLayGam;   //!
  TBranch        *b_IfrNsGam;   //!
  TBranch        *b_IfrInnerGam;   //!
  TBranch        *b_IfrBarrelGam;   //!
  TBranch        *b_IfrFWDGam;   //!
  TBranch        *b_IfrBWDGam;   //!
  TBranch        *b_IfrMeasIntLenGam;   //!
  TBranch        *b_IfrFirstHitGam;   //!
  TBranch        *b_IfrLastHitGam;   //!
  TBranch        *b_IfrExpIntLenGam;   //!
  TBranch        *b_IfrIntLenBeforeIronGam;   //!
  TBranch        *b_IfrTrkMatchGam;   //!
  TBranch        *b_IfrEmcMatchGam;   //!
  TBranch        *b_IfrLastBarrelGam;   //!
  TBranch        *b_IfrCLFitChi2Gam;   //!
  TBranch        *b_IfrStrips0;   //!
  TBranch        *b_IfrStrips1;   //!
  TBranch        *b_IfrStrips2;   //!
  TBranch        *b_IfrStrips3;   //!
  TBranch        *b_IfrStrips4;   //!
  TBranch        *b_IfrStrips5;   //!
  TBranch        *b_IfrStrips6;   //!
  TBranch        *b_IfrStrips7;   //!
  TBranch        *b_IfrStrips8;   //!
  TBranch        *b_IfrStrips9;   //!
  TBranch        *b_IfrStrips10;   //!
  TBranch        *b_IfrStrips11;   //!
  TBranch        *b_IfrStrips12;   //!
  TBranch        *b_IfrStrips13;   //!
  TBranch        *b_IfrStrips14;   //!
  TBranch        *b_IfrStrips15;   //!
  TBranch        *b_IfrStrips16;   //!
  TBranch        *b_IfrStrips17;   //!
  TBranch        *b_IfrStrips18;   //!
  TBranch        *b_IfrStrips19;   //!
  TBranch        *b_lMomGam;   //!
  TBranch        *b_ZMom42Gam;   //!
  TBranch        *b_ecalGam;   //!
  TBranch        *b_ecalXGam;   //!
  TBranch        *b_ecalYGam;   //!
  TBranch        *b_ecalZGam;   //!
  TBranch        *b_nCryGam;   //!
  TBranch        *b_nBumpGam;   //!
  TBranch        *b_ZMom20Gam;   //!
  TBranch        *b_secMomGam;   //!
  TBranch        *b_s1s9Gam;   //!
  TBranch        *b_s9s25Gam;   //!
  TBranch        *b_erawGam;   //!
  TBranch        *b_phiClusterGam;   //!
  TBranch        *b_thetaClusterGam;   //!
  TBranch        *b_covEEGam;   //!
  TBranch        *b_covTTGam;   //!
  TBranch        *b_covPPGam;   //!
  TBranch        *b_covRRGam;   //!
  TBranch        *b_emcStatusGam;   //!
  TBranch        *b_thetaGam;   //!
  TBranch        *b_phiGam;   //!
  TBranch        *b_energyGam;   //!
  TBranch        *b_idGam;   //!
  TBranch        *b_IndexGam;   //!
  TBranch        *b_IndexNtGam;   //!
  TBranch        *b_B0RecGam;   //!
  TBranch        *b_chBRecGam;   //!

};

// ----------------------------------------------------------------------
inline void b2uNtp::mk4Vector(TLorentzVector &p4, const float p, const float t, const float f, const float m) {
  p4.SetXYZM(p*TMath::Sin(t)*TMath::Cos(f), p*TMath::Sin(t)*TMath::Sin(f), p*TMath::Cos(t), m);
}

// ----------------------------------------------------------------------
inline void b2uNtp::mk3Vector(TVector3 &p3, const float p, const float t, const float f) {
  p3.SetXYZ(p*TMath::Sin(t)*TMath::Cos(f), p*TMath::Sin(t)*TMath::Sin(f), p*TMath::Cos(t));
}


// ----------------------------------------------------------------------
inline double b2uNtp::unwrapped(double phi1, double phi2) {
  double dphi = phi1 - phi2;
  while (dphi >  3.1415) dphi -= 6.2830;
  while (dphi < -3.1415) dphi += 6.2830;
  return TMath::Abs(dphi); 
}


// ----------------------------------------------------------------------
// -- The meaning of the bits in the bitmap is documented in  
// http://www.slac.stanford.edu/BFROOT/www/Physics/Analysis/AWG/Charmonium/CharmUser/bitmap.html
// xxxxIdTrk & 1  -> noCal/MIP
// xxxxIdTrk & 2  -> veryLoose
// xxxxIdTrk & 4  -> loose
// xxxxIdTrk & 8  -> tight
// xxxxIdTrk & 16 -> veryTight
// xxxxIdTrk & 32 -> likelihood 
inline Bool_t b2uNtp::isRecEl(int i)     {if (fRecEl[i] == 1) {return kTRUE;} else {return kFALSE;} }
inline Bool_t b2uNtp::isRecMu(int i)     {if (fRecMu[i] == 1) {return kTRUE;} else {return kFALSE;} }
inline Bool_t b2uNtp::isRecLepton(int i) {return (isRecEl(i) || isRecMu(i));}
inline Bool_t b2uNtp::isRecKaon(int i)   {if (fRecKa[i] == 1) {return kTRUE;} else {return kFALSE;} }

// -- PID for generator block
inline Bool_t b2uNtp::isTruEl(Int_t i)  {return (TMath::Abs(idMc[i]) == 11);}
inline Bool_t b2uNtp::isTruMu(Int_t i)  {return (TMath::Abs(idMc[i]) == 13);}
inline Bool_t b2uNtp::isTruTau(Int_t i) {return (TMath::Abs(idMc[i]) == 15);}
inline Bool_t b2uNtp::isTruLepton(Int_t i) {return (isTruEl(i) || isTruMu(i));}

#endif

