#ifndef baseClass_h
#define baseClass_h

#include <iostream>

#include <TROOT.h>
#include <TVector3.h>
#include <TEventList.h>
#include <TLorentzVector.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TList.h>
#include <TMap.h>
#include <TH1.h>

#include "PIDTable.hh"
#include "TRKTable.hh"
#include "recoilDSys.hh"

class baseClass {

protected:

  TEventList *fToBeCopied;
  TFile  *fHistFile;
  TTree  *fTree;
  TTree  *fGTree;
  TTree  *fSOTree;// split offs tree (tmp)

  long fReturnLog[100];
  TString fReturnString[100];
  Int_t   fSkimcount;
  Int_t   fDump, fLkt, fIVal, fDupcount, fNewFormat;
  
  Bool_t fIsMC;

  Bool_t  fIsMakeParam;
  TString fCutFile; 

  int fChB;

  // ----------------------------------------------------------------------
  // -- Binning and constants

  Double_t DR; 
  int PCBIN, PLBIN, FBIN, B0LUND, CHBLUND, XBIN;
  double PCMAX, PLMAX, XMAX;
  Double_t ELMASS, MUMASS, PIPMASS, PIZMASS, KAPMASS, KAZMASS, BPMASS, BZMASS, BMASS, BQMASS, A0;
  int findPro, findUps;




  // ----------------------------------------------------------------------
  // -- CUTS
  // -- Note: additional cuts may require modifications in initRest(), dumpCuts() and readCuts()

  Double_t MESSIGNALLO, MESSIGNALHI, MESSIDEBANDLO, MESSIDEBANDHI; 
  Double_t DESIGNALLO,  DESIGNALHI, MESSIGNALBANDLO, MESSIGNALBANDHI;
  Double_t PURITY, INTPURITY, IPURDSTAR, IPURDC, IPURDSTAR0, IPURD0; 
  Double_t PCMSLO, PLABLO, TLABLO, TLABHI, ELMOMLO, MUMOMLO, KAMOMLO; 
  Int_t NLEPTON, IDEL, IDMU, IDKA; 
  Int_t DOBDECWEIGHT, DODDECWEIGHT; 

  double SMEARELPIDTABLES, SMEARMUPIDTABLES, SMEARKAPIDTABLES, 
    SMEARELMISTABLES, SMEARMUMISTABLES, SMEARKAMISTABLES, 
    SMEARTRKPX, SMEARTRKPY, 
    SIGMANEUT, SHIFTNEUT, 
    TRACKSELECTION, PHOTONSELECTION;
  int DOTRACKKILLING;
  double SLOWPIONKILL, TRACKKILL; 
  char PIDTABLES[1000];
  char TRKTABLES[1000];
  char WEIGHTNEU[1000];
  char WEIGHTCHG[1000];

  Double_t PRMM2, MM2LO, MM2HI; 
  Bool_t   REQCHARGECORR;
  Double_t REQTOTALCHARGE; 
  Double_t KSPIPLO, KSPIPHI, KSPIPRLO, KSPIZLO, KSPIZHI; 
  Double_t PTLO, GAMMAELO, GAMMAEHI; 



  // ----------------------------------------------------------------------
  // -- Options

  Int_t fOptBlind, fOptGammas, fOptCategories, fOptKlongs, fOptMakeEventList, fOptFilterK0sPi0Pi0, 
    fOptPidKilling, fOptPidKillingKaon, fOptPidKillingEl, fOptPidKillingMu, fOptSmearTracks, fOptSmearNeut;
  Int_t fVerbose;



  // ----------------------------------------------------------------------
  // -- Event

  int fEvent;

  TString fRunRange;   // Run 1, Run 2, Run 3, ...
  int fFileChanged; // a flag to allow a faster comparison if the underlying file has changed

  TRKTable *fTT[5]; 
  PIDTable *fPTel[10], *fPTmu[10], *fPTka[10];
  double fEvtW8;   // filled with B/D BF reweighting, ...
  int bestB0;
  int indexbestB;
  Int_t brecoOverlap;
  TLorentzVector p4Upsilon; 
  TVector3 cmsBoost, upsBoost;

  recoilDSys *Dvar;
  recoilDSys *Bsem;


  // ----------------------------------------------------------------------
  // -- mcTruth

  Int_t fRunnumber, fLower, fUpper; 

  TLorentzVector p4XhadGen, p4MissGen, p4RecoilGen, p4LeptonGen; 
  Double_t fQ2Gen;
  Double_t fMxhadGen, fPxhadGen, fTxhadGen, fFxhadGen, fExhadGen;

  Int_t fBbchgen;
  Double_t fPcmsGen, fTcmsGen, fFcmsGen, fEcmsGen;
  Int_t    fB1Index, fB2Index, fBVxb, fBVcb, fBVub, fBVxbTyp, fMxCategory;

  Int_t  fDpi, fDpiz, fDk, fDkmiss,  fDks, fDCfDs, fD0CfDs, fDlep, fDgam, fBVSTyp;
  Int_t  fDkl, fDkspiopio, fDkspipi;


  // ----------------------------------------------------------------------
  // -- breco
  TLorentzVector p4Brecoil, p4Breco, p4BrecoNC, p4BrecoGen, p4BrecoGenSmear; 
  Double_t fPcmsBreco;
  Int_t lundB, fSeedMode, fBmode; 
  Bool_t mesSideband, mesSignalband, signalBox, tsdump; 
  Int_t fNLeptonMC, genLeptonCharge;
  Int_t fBrecoCharge, fBrecoFlavor;

  Int_t    fVub, fVcb, fOther, fBrecoMc;
  Int_t fNlep; // -- Number of leptons in the recoil (all)
  Double_t fPurity, fIntPurity; 
  Double_t fMes, fDeltaE, fBmass; 





  // ----------------------------------------------------------------------
  // -- recoil
  Double_t fQtot, fQrecoil, fNtracks, fNneutrals, frealNLepton, frealNKp, frealNKshort;
  Int_t fRecoilCharge, fLeptonCharge , fRecoilTrkMult, fRecoilNutMult;

  Int_t fNLepton, fNEl, fNMu, fNKp; 
  Bool_t fIsPrompt, fIsCascade;
  Bool_t fElectron, fMuon; 

  Double_t fGammaMax,fENeu,fEPiz;
  Double_t fQ2, fQ2NC, fQ2Fit, fQ2Res;
  Double_t fPNu, fTNu, fCosTNu, fFNu, fEmiss, fMM2,fMinKMom,fMaxKMom;
  Double_t fPNuNC, fTNuNC, fFNuNC, fEmissNC, fMM2NC;
  Double_t fPlab, fTlab, fTlabDR, fFlab, fFlabDR, fPcms, fTcms, fFcms, fEcms, fPups;
  Double_t fPxhad, fTxhad, fFxhad, fExhad, fMxhad, fMxhadRes; 
  Double_t fPxhadchg, fTxhadchg, fFxhadchg, fExhadchg, fMxhadchg;
  Double_t fPrecoil, fTrecoil, fFrecoil, fErecoil, fMrecoil;
  Double_t fPAllev, fTAllev, fFAllev, fEAllev, fMAllev;

  Int_t fcountChg;

  Double_t fMCharpart, fTNeupart[15], fTCharpart, fENeupart[15], fECharpart,  fPNeupart[15], fPCharpart;
  Double_t tmpfMxhad[15], tmpxmassRes[15], tmpfFxhad[15], tmpfTxhad[15], tmpfExhad[15];
  Double_t theweighttrk[20], theweightneu[20];
  Double_t brecosig[6000], brecobkg[6000], brecointpur[6000];  

  // split off quantities (tmp)
  Double_t SOdth,SOdph,SOegam,SOptrk,SOthtrk,SOphitrk,SOthgam,SOphigam,SOecaltrk;
  Int_t SOcharge;

  // ----------------------------------------------------------------------
  // -- kinFit'ed quantities

  Int_t fParametrization; 
  Double_t fMxhadfit;
  Double_t fMxhadfitRes; 

  TLorentzVector p4ChargPart;
  TLorentzVector p4Xhad, p4Miss, p4Recoil, p4LeptonLab, p4LeptonCms, p4LeptonUps;

  Int_t recEl[100], recMu[100], recKa[100];
  Int_t splitOffGam[100], goodTrack[100], goodChargedKaon[100], goodHadron[100], 
        goodPion[100], goodPhoton[100], loopTrack[100], kshortLockTrk[100];

  // ----------------------------------------------------------------------
  // -- MX study


   TTree          *fChain;   //pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //current Tree number in a TChain

   Int_t           platform;
   Int_t           partition;

   Int_t           event;
   Int_t           nnpi0;
   Int_t           nnks;
   Int_t           nntrk;
   Int_t           nnpar;
   Int_t           IndexTrk[60];
   Int_t           idTrk[60];
   Float_t         momentumTrk[60];
   Float_t         thetaTrk[60];
   Float_t         phiTrk[60];
   Int_t           nTrk;
   Int_t           runNumber;
   Int_t           upperID;
   Int_t           lowerID;
   Int_t           nGam;
   Int_t           IndexGam[80];
   Int_t           d1ChBIndex[100];
   Int_t           d1ChBLund[100];
   Int_t           d2ChBIndex[100];
   Int_t           d2ChBLund[100];
   Int_t           d3ChBIndex[100];
   Int_t           d3ChBLund[100];
   Int_t           d4ChBIndex[100];
   Int_t           d4ChBLund[100];
   Int_t           d5ChBIndex[100];
   Int_t           d5ChBLund[100];
   Int_t           d6ChBIndex[100];
   Int_t           d6ChBLund[100];
   Int_t           d7ChBIndex[100];
   Int_t           d7ChBLund[100];
   Int_t           modeChB[100];
   Float_t         purChB[100];
   Float_t         intpurChB[100];
   Int_t           d1B0Index[100];
   Int_t           d1B0Lund[100];
   Int_t           d2B0Index[100];
   Int_t           d2B0Lund[100];
   Int_t           d3B0Index[100];
   Int_t           d3B0Lund[100];
   Int_t           d4B0Index[100];
   Int_t           d4B0Lund[100];
   Int_t           d5B0Index[100];
   Int_t           d5B0Lund[100];
   Int_t           d6B0Index[100];
   Int_t           d6B0Lund[100];
   Int_t           d7B0Index[100];
   Int_t           d7B0Lund[100];
   Float_t         ppcov00[60];
   Float_t         ppcov10[60];
   Float_t         ppcov11[60];
   Float_t         ppcov20[60];
   Float_t         ppcov21[60];
   Float_t         ppcov22[60];
   Int_t           nMc;
   Float_t         pMc[200];
   Float_t         massMc[200];
   Float_t         thetaMc[200];
   Float_t         phiMc[200];
   Int_t           idMc[200];
   Int_t           mothMc[200];
   UInt_t          nDauMc[200];
   Float_t         energyGam[80];
   Float_t         thetaGam[80];
   Float_t         xPocaTrk[60];
   Float_t         yPocaTrk[60];
   Float_t         zPocaTrk[60];
   Float_t         beamSX;
   Float_t         beamSY;
   Float_t         beamSZ;
   Int_t           isvtTrk[60];
   Int_t           nsvtTrk[60];
   Int_t           fhitTrk[60];
   Int_t           ndchTrk[60];
   Int_t           lhitTrk[60];
   Float_t         tLenTrk[60];
   Int_t           ntdofTrk[60];
   Float_t         tproTrk[60];
   Float_t         tChi2Trk[60];
   Int_t           cPidTrk[60];
   Int_t           nCryGam[80];
   Float_t         lMomGam[80];
   Float_t         phiGam[80];
   Float_t         s1s9Gam[80];
   Float_t         s9s25Gam[80];
   Int_t           muonIdTrk[60];
   Int_t           elecIdTrk[60];
   Int_t           kaonIdTrk[60];
   Int_t           pionIdTrk[60];
   Int_t           nB0;
   Float_t         massB0[100];
   Float_t         pB0[100];
   Float_t         thB0[100];
   Float_t         phiB0[100];
   Float_t         errMassB0[100];
   Float_t         mseB0[100];
   Float_t         mseChB[100];
   Float_t         ecalGam[80];
   Float_t         ecalXGam[80];
   Float_t         ecalYGam[80];
   Float_t         ecalZGam[80];
   Int_t           modeB0[100];
   Int_t           nChB;
   Float_t         secMomGam[80];
   Float_t         ZMom20Gam[80];
   Float_t         ZMom42Gam[80];
   Int_t           IfrStrips0[80];
   Int_t           IfrStrips1[80];
   Int_t           IfrStrips2[80];
   Int_t           IfrStrips3[80];
   Int_t           IfrStrips4[80];
   Int_t           IfrStrips5[80];
   Int_t           IfrStrips6[80];
   Int_t           IfrStrips7[80];
   Int_t           IfrStrips8[80];
   Int_t           IfrStrips9[80];
   Int_t           IfrStrips10[80];
   Int_t           IfrStrips11[80];
   Int_t           IfrStrips12[80];
   Int_t           IfrStrips13[80];
   Int_t           IfrStrips14[80];
   Int_t           IfrStrips15[80];
   Int_t           IfrStrips16[80];
   Int_t           IfrStrips17[80];
   Int_t           IfrStrips18[80];
   Int_t           IfrStrips19[80];
   Int_t           IfrFirstHitGam[80];
   Int_t           IfrLastHitGam[80];
   Int_t           IfrLayGam[80];
   Int_t           IfrNsGam[80];
   Int_t           chargeTrk[60];
   Float_t         primVtxX;
   Float_t         primVtxY;
   Float_t         primVtxZ;
   Int_t           B0RecTrk[60];
   Int_t           chBRecTrk[60];
   Int_t           nPi0;
   Int_t           d1Pi0Index[150];
   Int_t           d1Pi0Lund[150];
   Int_t           d2Pi0Index[150];
   Int_t           d2Pi0Lund[150];
   Float_t         pxUps;
   Float_t         pyUps;
   Float_t         pzUps;
   Float_t         eUps;
   Float_t         VtxXXB0[100];
   Float_t         VtxYXB0[100];
   Float_t         VtxZXB0[100];
   Float_t         VtxXLepB0[100];
   Float_t         VtxYLepB0[100];
   Float_t         VtxZLepB0[100];
   Float_t         VtxCovXXLepB0[100];
   Float_t         VtxCovYYLepB0[100];
   Float_t         VtxCovXYLepB0[100];
   Float_t         VtxCovZZLepB0[100];
   Float_t         VtxCovXZLepB0[100];
   Float_t         VtxCovYZLepB0[100];
   Float_t         VtxCovXXXB0[100];
   Float_t         VtxCovYYXB0[100];
   Float_t         VtxCovXYXB0[100];
   Float_t         VtxCovZZXB0[100];
   Float_t         VtxCovXZXB0[100];
   Float_t         VtxCovYZXB0[100];
   Int_t           B0RecGam[80];
   Int_t           chBRecGam[80];
   Float_t         VtxXXChB[100];
   Float_t         VtxYXChB[100];
   Float_t         VtxZXChB[100];
   Float_t         VtxXLepChB[100];
   Float_t         VtxYLepChB[100];
   Float_t         VtxZLepChB[100];
   Float_t         VtxCovXXLepChB[100];
   Float_t         VtxCovYYLepChB[100];
   Float_t         VtxCovXYLepChB[100];
   Float_t         VtxCovZZLepChB[100];
   Float_t         VtxCovXZLepChB[100];
   Float_t         VtxCovYZLepChB[100];
   Float_t         VtxChiSqLepChB[100];
   Float_t         VtxNDofLepChB[100];
   Int_t           VtxStatLepChB[100];
   Int_t           VtxNUsedLepChB[100];
   Int_t           d1KsIndex[50];
   Int_t           d1KsLund[50];
   Int_t           d2KsIndex[50];
   Int_t           d2KsLund[50];
   Float_t         massKs[50];
   Float_t         pKs[50];
   Float_t         thKs[50];
   Float_t         phiKs[50];
   Float_t         errMassKs[50];
   Float_t         m0Ks[50];
   Float_t         xKs[50];
   Float_t         yKs[50];
   Float_t         zKs[50];
   Float_t         VtxCovXXXChB[100];
   Float_t         VtxCovYYXChB[100];
   Float_t         VtxCovXYXChB[100];
   Float_t         VtxCovZZXChB[100];
   Float_t         VtxCovXZXChB[100];
   Float_t         VtxCovYZXChB[100];
   Float_t         VtxChiSqXChB[100];
   Float_t         phiAtEMCTrk[60];
   Float_t         thetaAtEMCTrk[60];
   Float_t         ecalTrk[60];
   Float_t         xMc[200];
   Float_t         yMc[200];
   Float_t         zMc[200];
   Float_t         MassPB0[100];
   Float_t         MassPhiB0[100];
   Float_t         MassThetaB0[100];
   Float_t         massChB[100];
   Float_t         pChB[100];
   Float_t         thChB[100];
   Float_t         phiChB[100];
   Float_t         beamSCovXX;
   Float_t         beamSCovYY;
   Float_t         beamSCovZZ;
   Float_t         beamSCovXZ;
   Int_t           nTrkTot;
   Float_t         W2;
   Float_t         FoxWol2;
   Float_t         FoxWol2Neu;
   Float_t         thrust;
   Float_t         thrustNeu;
   Float_t         m0B0[100];
   Float_t         xB0[100];
   Float_t         yB0[100];
   Float_t         zB0[100];
   Float_t         s2xB0[100];
   Float_t         s2yB0[100];
   Float_t         s2zB0[100];
   Float_t         chi2B0[100];
   Int_t           dofB0[100];
   Int_t           stB0[100];
   Int_t           ndauB0[100];
   Int_t           MCB0[100];
   Float_t         mHatB0[100];
   Float_t         deltaeB0[100];
   Float_t         ThruB0[100];
   Float_t         thThruB0[100];
   Float_t         phiThruB0[100];
   Float_t         cosTBB0[100];
   Float_t         purB0[100];
   Float_t         intpurB0[100];
   Float_t         VtxChiSqLepB0[100];
   Float_t         VtxNDofLepB0[100];
   Int_t           VtxStatLepB0[100];
   Int_t           VtxNUsedLepB0[100];
   Float_t         DocaLepB0[100];
   Float_t         DocaErrLepB0[100];
   Float_t         VtxChiSqXB0[100];
   Float_t         VtxNDofXB0[100];
   Int_t           VtxStatXB0[100];
   Int_t           VtxNUsedXB0[100];
   Float_t         VtxPXB0[100];
   Float_t         VtxPhiXB0[100];
   Float_t         VtxThetaXB0[100];
   Float_t         ThrustXB0[100];
   Float_t         ThrustXPhiB0[100];
   Float_t         ThrustXThetaB0[100];
   Float_t         Cov00B0[100];
   Float_t         Cov10B0[100];
   Float_t         Cov11B0[100];
   Float_t         Cov20B0[100];
   Float_t         Cov21B0[100];
   Float_t         Cov22B0[100];
   Float_t         Cov30B0[100];
   Float_t         Cov31B0[100];
   Float_t         Cov32B0[100];
   Float_t         Cov33B0[100];
   Float_t         errMassChB[100];
   Float_t         m0ChB[100];
   Float_t         xChB[100];
   Float_t         yChB[100];
   Float_t         zChB[100];
   Float_t         s2xChB[100];
   Float_t         s2yChB[100];
   Float_t         s2zChB[100];
   Float_t         chi2ChB[100];
   Int_t           dofChB[100];
   Int_t           stChB[100];
   Int_t           ndauChB[100];
   Int_t           MCChB[100];
   Float_t         mHatChB[100];
   Float_t         deltaeChB[100];
   Float_t         ThruChB[100];
   Float_t         thThruChB[100];
   Float_t         phiThruChB[100];
   Float_t         cosTBChB[100];
   Float_t         MassPChB[100];
   Float_t         MassPhiChB[100];
   Float_t         MassThetaChB[100];
   Float_t         DocaLepChB[100];
   Float_t         DocaErrLepChB[100];
   Float_t         VtxNDofXChB[100];
   Int_t           VtxStatXChB[100];
   Int_t           VtxNUsedXChB[100];
   Float_t         VtxPXChB[100];
   Float_t         VtxPhiXChB[100];
   Float_t         VtxThetaXChB[100];
   Float_t         ThrustXChB[100];
   Float_t         ThrustXPhiChB[100];
   Float_t         ThrustXThetaChB[100];
   Float_t         Cov00ChB[100];
   Float_t         Cov10ChB[100];
   Float_t         Cov11ChB[100];
   Float_t         Cov20ChB[100];
   Float_t         Cov21ChB[100];
   Float_t         Cov22ChB[100];
   Float_t         Cov30ChB[100];
   Float_t         Cov31ChB[100];
   Float_t         Cov32ChB[100];
   Float_t         Cov33ChB[100];
   Int_t           nDstar;
   Float_t         massDstar[100];
   Float_t         pDstar[100];
   Float_t         thDstar[100];
   Float_t         phiDstar[100];
   Float_t         errMassDstar[100];
   Float_t         m0Dstar[100];
   Float_t         xDstar[100];
   Float_t         yDstar[100];
   Float_t         zDstar[100];
   Float_t         s2xDstar[100];
   Float_t         s2yDstar[100];
   Float_t         s2zDstar[100];
   Float_t         chi2Dstar[100];
   Int_t           dofDstar[100];
   Int_t           stDstar[100];
   Int_t           ndauDstar[100];
   Int_t           MCDstar[100];
   Int_t           d1DstarIndex[100];
   Int_t           d1DstarLund[100];
   Int_t           d2DstarIndex[100];
   Int_t           d2DstarLund[100];
   Int_t           nDstarBS;
   Float_t         massDstarBS[100];
   Float_t         chi2DstarBS[100];
   Int_t           dofDstarBS[100];
   Int_t           nDstar0;
   Float_t         massDstar0[100];
   Float_t         pDstar0[100];
   Float_t         thDstar0[100];
   Float_t         phiDstar0[100];
   Float_t         errMassDstar0[100];
   Float_t         m0Dstar0[100];
   Float_t         xDstar0[100];
   Float_t         yDstar0[100];
   Float_t         zDstar0[100];
   Float_t         s2xDstar0[100];
   Float_t         s2yDstar0[100];
   Float_t         s2zDstar0[100];
   Float_t         chi2Dstar0[100];
   Int_t           dofDstar0[100];
   Int_t           stDstar0[100];
   Int_t           ndauDstar0[100];
   Int_t           MCDstar0[100];
   Int_t           d1Dstar0Index[100];
   Int_t           d1Dstar0Lund[100];
   Int_t           d2Dstar0Index[100];
   Int_t           d2Dstar0Lund[100];
   Int_t           nD0;
   Float_t         massD0[100];
   Float_t         pD0[100];
   Float_t         thD0[100];
   Float_t         phiD0[100];
   Float_t         errMassD0[100];
   Float_t         m0D0[100];
   Float_t         xD0[100];
   Float_t         yD0[100];
   Float_t         zD0[100];
   Float_t         s2xD0[100];
   Float_t         s2yD0[100];
   Float_t         s2zD0[100];
   Float_t         chi2D0[100];
   Int_t           dofD0[100];
   Int_t           stD0[100];
   Int_t           ndauD0[100];
   Int_t           MCD0[100];
   Int_t           d1D0Index[100];
   Int_t           d1D0Lund[100];
   Int_t           d2D0Index[100];
   Int_t           d2D0Lund[100];
   Int_t           d3D0Index[100];
   Int_t           d3D0Lund[100];
   Int_t           d4D0Index[100];
   Int_t           d4D0Lund[100];
   Int_t           nChD;
   Float_t         massChD[100];
   Float_t         pChD[100];
   Float_t         thChD[100];
   Float_t         phiChD[100];
   Float_t         errMassChD[100];
   Float_t         m0ChD[100];
   Float_t         xChD[100];
   Float_t         yChD[100];
   Float_t         zChD[100];
   Float_t         s2xChD[100];
   Float_t         s2yChD[100];
   Float_t         s2zChD[100];
   Float_t         chi2ChD[100];
   Int_t           dofChD[100];
   Int_t           stChD[100];
   Int_t           ndauChD[100];
   Int_t           MCChD[100];
   Int_t           d1ChDIndex[100];
   Int_t           d1ChDLund[100];
   Int_t           d2ChDIndex[100];
   Int_t           d2ChDLund[100];
   Int_t           d3ChDIndex[100];
   Int_t           d3ChDLund[100];
   Int_t           d4ChDIndex[100];
   Int_t           d4ChDLund[100];
   Float_t         s2xKs[50];
   Float_t         s2yKs[50];
   Float_t         s2zKs[50];
   Float_t         chi2Ks[50];
   Int_t           dofKs[50];
   Int_t           stKs[50];
   Int_t           ndauKs[50];
   Float_t         massPi0[150];
   Float_t         pPi0[150];
   Float_t         thPi0[150];
   Float_t         phiPi0[150];
   Float_t         errMassPi0[150];
   Float_t         m0Pi0[150];
   Float_t         xPi0[150];
   Float_t         yPi0[150];
   Float_t         zPi0[150];
   Float_t         s2xPi0[150];
   Float_t         s2yPi0[150];
   Float_t         s2zPi0[150];
   Float_t         chi2Pi0[150];
   Int_t           dofPi0[150];
   Int_t           stPi0[150];
   Int_t           ndauPi0[150];
   Int_t           MCPi0[150];
   Int_t           nGConv;
   Float_t         massGConv[10];
   Float_t         pGConv[10];
   Float_t         thGConv[10];
   Float_t         phiGConv[10];
   Float_t         errMassGConv[10];
   Float_t         m0GConv[10];
   Float_t         xGConv[10];
   Float_t         yGConv[10];
   Float_t         zGConv[10];
   Float_t         s2xGConv[10];
   Float_t         s2yGConv[10];
   Float_t         s2zGConv[10];
   Float_t         chi2GConv[10];
   Int_t           dofGConv[10];
   Int_t           stGConv[10];
   Int_t           ndauGConv[10];
   Int_t           MCGConv[10];
   Int_t           d1GConvIndex[10];
   Int_t           d1GConvLund[10];
   Int_t           d2GConvIndex[10];
   Int_t           d2GConvLund[10];
   Int_t           nDalitz;
   Float_t         massDalitz[10];
   Float_t         pDalitz[10];
   Float_t         thDalitz[10];
   Float_t         phiDalitz[10];
   Float_t         errMassDalitz[10];
   Float_t         m0Dalitz[10];
   Float_t         xDalitz[10];
   Float_t         yDalitz[10];
   Float_t         zDalitz[10];
   Float_t         s2xDalitz[10];
   Float_t         s2yDalitz[10];
   Float_t         s2zDalitz[10];
   Float_t         chi2Dalitz[10];
   Int_t           dofDalitz[10];
   Int_t           stDalitz[10];
   Int_t           ndauDalitz[10];
   Int_t           MCDalitz[10];
   Int_t           d1DalitzIndex[10];
   Int_t           d1DalitzLund[10];
   Int_t           d2DalitzIndex[10];
   Int_t           d2DalitzLund[10];
   Int_t           nJpsi;
   Float_t         massJpsi[10];
   Float_t         pJpsi[10];
   Float_t         thJpsi[10];
   Float_t         phiJpsi[10];
   Float_t         errMassJpsi[10];
   Float_t         m0Jpsi[10];
   Float_t         xJpsi[10];
   Float_t         yJpsi[10];
   Float_t         zJpsi[10];
   Float_t         s2xJpsi[10];
   Float_t         s2yJpsi[10];
   Float_t         s2zJpsi[10];
   Float_t         chi2Jpsi[10];
   Int_t           dofJpsi[10];
   Int_t           stJpsi[10];
   Int_t           ndauJpsi[10];
   Int_t           MCJpsi[10];
   Int_t           d1JpsiIndex[10];
   Int_t           d1JpsiLund[10];
   Int_t           d2JpsiIndex[10];
   Int_t           d2JpsiLund[10];
   Int_t           IfrLayTrk[60];
   Int_t           IfrNsTrk[60];
   UChar_t         IfrInnerTrk[60];
   UChar_t         IfrBarrelTrk[60];
   UChar_t         IfrFWDTrk[60];
   UChar_t         IfrBWDTrk[60];
   Float_t         IfrMeasIntLenTrk[60];
   Int_t           IfrFirstHitTrk[60];
   Int_t           IfrLastHitTrk[60];
   Float_t         lMomTrk[60];
   Float_t         ZMom42Trk[60];
   Float_t         ecalXTrk[60];
   Float_t         ecalYTrk[60];
   Float_t         ecalZTrk[60];
   Int_t           nCryTrk[60];
   Int_t           nBumpTrk[60];
   Float_t         ZMom20Trk[60];
   Float_t         secMomTrk[60];
   Float_t         s1s9Trk[60];
   Float_t         s9s25Trk[60];
   Float_t         erawTrk[60];
   Float_t         phiClusterTrk[60];
   Float_t         thetaClusterTrk[60];
   Float_t         covEETrk[60];
   Float_t         covTTTrk[60];
   Float_t         covPPTrk[60];
   Float_t         covRRTrk[60];
   Float_t         phicMatTrk[60];
   Float_t         trkcMatTrk[60];
   Int_t           nPidTrk[60];
   Int_t           emcStatusTrk[60];
   Int_t           MCKs[50];
   Float_t         sfRangeTrk[60];
   UChar_t         trkFitStatusTrk[60];
   Int_t           IndexNtTrk[60];
   UChar_t         IfrInnerGam[80];
   UChar_t         IfrBarrelGam[80];
   UChar_t         IfrFWDGam[80];
   UChar_t         IfrBWDGam[80];
   Float_t         IfrMeasIntLenGam[80];
   Float_t         IfrExpIntLenGam[80];
   Float_t         IfrTrkMatchGam[80];
   Float_t         IfrEmcMatchGam[80];
   Int_t           IfrLastBarrelGam[80];
   Float_t         IfrCLFitChi2Gam[80];
   Float_t         IfrIntLenBeforeIronGam[80];
   Int_t           nBumpGam[80];
   Float_t         erawGam[80];
   Float_t         phiClusterGam[80];
   Float_t         thetaClusterGam[80];
   Float_t         covEEGam[80];
   Float_t         covTTGam[80];
   Float_t         covPPGam[80];
   Float_t         covRRGam[80];
   Int_t           emcStatusGam[80];
   Int_t           idGam[80];
   Int_t           IndexNtGam[80];
   Int_t           nKs;



//List of branches
   TBranch        *b_event;
   TBranch        *b_runNumber;
   TBranch        *b_platform;
   TBranch        *b_partition;
   TBranch        *b_upperID;
   TBranch        *b_lowerID;
   TBranch        *b_beamSX;
   TBranch        *b_beamSY;
   TBranch        *b_beamSZ;
   TBranch        *b_primVtxX;
   TBranch        *b_primVtxY;
   TBranch        *b_primVtxZ;
   TBranch        *b_beamSCovXX;
   TBranch        *b_beamSCovYY;
   TBranch        *b_beamSCovZZ;
   TBranch        *b_beamSCovXZ;
   TBranch        *b_pxUps;
   TBranch        *b_pyUps;
   TBranch        *b_pzUps;
   TBranch        *b_eUps;
   TBranch        *b_nTrkTot;
   TBranch        *b_W2;
   TBranch        *b_FoxWol2;
   TBranch        *b_FoxWol2Neu;
   TBranch        *b_thrust;
   TBranch        *b_thrustNeu;
   TBranch        *b_nMc;
   TBranch        *b_pMc;
   TBranch        *b_massMc;
   TBranch        *b_thetaMc;
   TBranch        *b_phiMc;
   TBranch        *b_idMc;
   TBranch        *b_mothMc;
   TBranch        *b_nDauMc;
   TBranch        *b_xMc;
   TBranch        *b_yMc;
   TBranch        *b_zMc;
   TBranch        *b_nB0;
   TBranch        *b_massB0;
   TBranch        *b_pB0;
   TBranch        *b_thB0;
   TBranch        *b_phiB0;
   TBranch        *b_errMassB0;
   TBranch        *b_m0B0;
   TBranch        *b_xB0;
   TBranch        *b_yB0;
   TBranch        *b_zB0;
   TBranch        *b_s2xB0;
   TBranch        *b_s2yB0;
   TBranch        *b_s2zB0;
   TBranch        *b_chi2B0;
   TBranch        *b_dofB0;
   TBranch        *b_stB0;
   TBranch        *b_ndauB0;
   TBranch        *b_MCB0;
   TBranch        *b_mseB0;
   TBranch        *b_mHatB0;
   TBranch        *b_deltaeB0;
   TBranch        *b_ThruB0;
   TBranch        *b_thThruB0;
   TBranch        *b_phiThruB0;
   TBranch        *b_cosTBB0;
   TBranch        *b_d1B0Index;
   TBranch        *b_d1B0Lund;
   TBranch        *b_d2B0Index;
   TBranch        *b_d2B0Lund;
   TBranch        *b_d3B0Index;
   TBranch        *b_d3B0Lund;
   TBranch        *b_d4B0Index;
   TBranch        *b_d4B0Lund;
   TBranch        *b_d5B0Index;
   TBranch        *b_d5B0Lund;
   TBranch        *b_d6B0Index;
   TBranch        *b_d6B0Lund;
   TBranch        *b_d7B0Index;
   TBranch        *b_d7B0Lund;
   TBranch        *b_modeB0;
   TBranch        *b_purB0;
   TBranch        *b_intpurB0;
   TBranch        *b_VtxXLepB0;
   TBranch        *b_VtxYLepB0;
   TBranch        *b_VtxZLepB0;
   TBranch        *b_VtxCovXXLepB0;
   TBranch        *b_VtxCovYYLepB0;
   TBranch        *b_VtxCovXYLepB0;
   TBranch        *b_VtxCovZZLepB0;
   TBranch        *b_VtxCovXZLepB0;
   TBranch        *b_VtxCovYZLepB0;
   TBranch        *b_VtxChiSqLepB0;
   TBranch        *b_VtxNDofLepB0;
   TBranch        *b_VtxStatLepB0;
   TBranch        *b_VtxNUsedLepB0;
   TBranch        *b_DocaLepB0;
   TBranch        *b_DocaErrLepB0;
   TBranch        *b_VtxXXB0;
   TBranch        *b_VtxYXB0;
   TBranch        *b_VtxZXB0;
   TBranch        *b_VtxCovXXXB0;
   TBranch        *b_VtxCovYYXB0;
   TBranch        *b_VtxCovXYXB0;
   TBranch        *b_VtxCovZZXB0;
   TBranch        *b_VtxCovXZXB0;
   TBranch        *b_VtxCovYZXB0;
   TBranch        *b_VtxChiSqXB0;
   TBranch        *b_VtxNDofXB0;
   TBranch        *b_VtxStatXB0;
   TBranch        *b_VtxNUsedXB0;
   TBranch        *b_VtxPXB0;
   TBranch        *b_VtxPhiXB0;
   TBranch        *b_VtxThetaXB0;
   TBranch        *b_ThrustXB0;
   TBranch        *b_ThrustXPhiB0;
   TBranch        *b_ThrustXThetaB0;
   TBranch        *b_MassPB0;
   TBranch        *b_MassPhiB0;
   TBranch        *b_MassThetaB0;
   TBranch        *b_Cov00B0;
   TBranch        *b_Cov10B0;
   TBranch        *b_Cov11B0;
   TBranch        *b_Cov20B0;
   TBranch        *b_Cov21B0;
   TBranch        *b_Cov22B0;
   TBranch        *b_Cov30B0;
   TBranch        *b_Cov31B0;
   TBranch        *b_Cov32B0;
   TBranch        *b_Cov33B0;
   TBranch        *b_nChB;
   TBranch        *b_massChB;
   TBranch        *b_pChB;
   TBranch        *b_thChB;
   TBranch        *b_phiChB;
   TBranch        *b_errMassChB;
   TBranch        *b_m0ChB;
   TBranch        *b_xChB;
   TBranch        *b_yChB;
   TBranch        *b_zChB;
   TBranch        *b_s2xChB;
   TBranch        *b_s2yChB;
   TBranch        *b_s2zChB;
   TBranch        *b_chi2ChB;
   TBranch        *b_dofChB;
   TBranch        *b_stChB;
   TBranch        *b_ndauChB;
   TBranch        *b_MCChB;
   TBranch        *b_mseChB;
   TBranch        *b_mHatChB;
   TBranch        *b_deltaeChB;
   TBranch        *b_ThruChB;
   TBranch        *b_thThruChB;
   TBranch        *b_phiThruChB;
   TBranch        *b_cosTBChB;
   TBranch        *b_d1ChBIndex;
   TBranch        *b_d1ChBLund;
   TBranch        *b_d2ChBIndex;
   TBranch        *b_d2ChBLund;
   TBranch        *b_d3ChBIndex;
   TBranch        *b_d3ChBLund;
   TBranch        *b_d4ChBIndex;
   TBranch        *b_d4ChBLund;
   TBranch        *b_d5ChBIndex;
   TBranch        *b_d5ChBLund;
   TBranch        *b_d6ChBIndex;
   TBranch        *b_d6ChBLund;
   TBranch        *b_d7ChBIndex;
   TBranch        *b_d7ChBLund;
   TBranch        *b_modeChB;
   TBranch        *b_purChB;
   TBranch        *b_intpurChB;
   TBranch        *b_VtxXLepChB;
   TBranch        *b_VtxYLepChB;
   TBranch        *b_VtxZLepChB;
   TBranch        *b_VtxCovXXLepChB;
   TBranch        *b_VtxCovYYLepChB;
   TBranch        *b_VtxCovXYLepChB;
   TBranch        *b_VtxCovZZLepChB;
   TBranch        *b_VtxCovXZLepChB;
   TBranch        *b_VtxCovYZLepChB;
   TBranch        *b_VtxChiSqLepChB;
   TBranch        *b_VtxNDofLepChB;
   TBranch        *b_VtxStatLepChB;
   TBranch        *b_VtxNUsedLepChB;
   TBranch        *b_DocaLepChB;
   TBranch        *b_DocaErrLepChB;
   TBranch        *b_VtxXXChB;
   TBranch        *b_VtxYXChB;
   TBranch        *b_VtxZXChB;
   TBranch        *b_VtxCovXXXChB;
   TBranch        *b_VtxCovYYXChB;
   TBranch        *b_VtxCovXYXChB;
   TBranch        *b_VtxCovZZXChB;
   TBranch        *b_VtxCovXZXChB;
   TBranch        *b_VtxCovYZXChB;
   TBranch        *b_VtxChiSqXChB;
   TBranch        *b_VtxNDofXChB;
   TBranch        *b_VtxStatXChB;
   TBranch        *b_VtxNUsedXChB;
   TBranch        *b_VtxPXChB;
   TBranch        *b_VtxPhiXChB;
   TBranch        *b_VtxThetaXChB;
   TBranch        *b_ThrustXChB;
   TBranch        *b_ThrustXPhiChB;
   TBranch        *b_ThrustXThetaChB;
   TBranch        *b_MassPChB;
   TBranch        *b_MassPhiChB;
   TBranch        *b_MassThetaChB;
   TBranch        *b_Cov00ChB;
   TBranch        *b_Cov10ChB;
   TBranch        *b_Cov11ChB;
   TBranch        *b_Cov20ChB;
   TBranch        *b_Cov21ChB;
   TBranch        *b_Cov22ChB;
   TBranch        *b_Cov30ChB;
   TBranch        *b_Cov31ChB;
   TBranch        *b_Cov32ChB;
   TBranch        *b_Cov33ChB;
   TBranch        *b_nDstar;
   TBranch        *b_massDstar;
   TBranch        *b_pDstar;
   TBranch        *b_thDstar;
   TBranch        *b_phiDstar;
   TBranch        *b_errMassDstar;
   TBranch        *b_m0Dstar;
   TBranch        *b_xDstar;
   TBranch        *b_yDstar;
   TBranch        *b_zDstar;
   TBranch        *b_s2xDstar;
   TBranch        *b_s2yDstar;
   TBranch        *b_s2zDstar;
   TBranch        *b_chi2Dstar;
   TBranch        *b_dofDstar;
   TBranch        *b_stDstar;
   TBranch        *b_ndauDstar;
   TBranch        *b_MCDstar;
   TBranch        *b_d1DstarIndex;
   TBranch        *b_d1DstarLund;
   TBranch        *b_d2DstarIndex;
   TBranch        *b_d2DstarLund;
   TBranch        *b_nDstarBS;
   TBranch        *b_massDstarBS;
   TBranch        *b_chi2DstarBS;
   TBranch        *b_dofDstarBS;
   TBranch        *b_spixDstarBS;
   TBranch        *b_spiyDstarBS;
   TBranch        *b_spizDstarBS;
   TBranch        *b_nDstar0;
   TBranch        *b_massDstar0;
   TBranch        *b_pDstar0;
   TBranch        *b_thDstar0;
   TBranch        *b_phiDstar0;
   TBranch        *b_errMassDstar0;
   TBranch        *b_m0Dstar0;
   TBranch        *b_xDstar0;
   TBranch        *b_yDstar0;
   TBranch        *b_zDstar0;
   TBranch        *b_s2xDstar0;
   TBranch        *b_s2yDstar0;
   TBranch        *b_s2zDstar0;
   TBranch        *b_chi2Dstar0;
   TBranch        *b_dofDstar0;
   TBranch        *b_stDstar0;
   TBranch        *b_ndauDstar0;
   TBranch        *b_MCDstar0;
   TBranch        *b_d1Dstar0Index;
   TBranch        *b_d1Dstar0Lund;
   TBranch        *b_d2Dstar0Index;
   TBranch        *b_d2Dstar0Lund;
   TBranch        *b_nD0;
   TBranch        *b_massD0;
   TBranch        *b_pD0;
   TBranch        *b_thD0;
   TBranch        *b_phiD0;
   TBranch        *b_errMassD0;
   TBranch        *b_m0D0;
   TBranch        *b_xD0;
   TBranch        *b_yD0;
   TBranch        *b_zD0;
   TBranch        *b_s2xD0;
   TBranch        *b_s2yD0;
   TBranch        *b_s2zD0;
   TBranch        *b_chi2D0;
   TBranch        *b_dofD0;
   TBranch        *b_stD0;
   TBranch        *b_ndauD0;
   TBranch        *b_MCD0;
   TBranch        *b_d1D0Index;
   TBranch        *b_d1D0Lund;
   TBranch        *b_d2D0Index;
   TBranch        *b_d2D0Lund;
   TBranch        *b_d3D0Index;
   TBranch        *b_d3D0Lund;
   TBranch        *b_d4D0Index;
   TBranch        *b_d4D0Lund;
   TBranch        *b_nChD;
   TBranch        *b_massChD;
   TBranch        *b_pChD;
   TBranch        *b_thChD;
   TBranch        *b_phiChD;
   TBranch        *b_errMassChD;
   TBranch        *b_m0ChD;
   TBranch        *b_xChD;
   TBranch        *b_yChD;
   TBranch        *b_zChD;
   TBranch        *b_s2xChD;
   TBranch        *b_s2yChD;
   TBranch        *b_s2zChD;
   TBranch        *b_chi2ChD;
   TBranch        *b_dofChD;
   TBranch        *b_stChD;
   TBranch        *b_ndauChD;
   TBranch        *b_MCChD;
   TBranch        *b_d1ChDIndex;
   TBranch        *b_d1ChDLund;
   TBranch        *b_d2ChDIndex;
   TBranch        *b_d2ChDLund;
   TBranch        *b_d3ChDIndex;
   TBranch        *b_d3ChDLund;
   TBranch        *b_d4ChDIndex;
   TBranch        *b_d4ChDLund;
   TBranch        *b_nKs;
   TBranch        *b_massKs;
   TBranch        *b_pKs;
   TBranch        *b_thKs;
   TBranch        *b_phiKs;
   TBranch        *b_errMassKs;
   TBranch        *b_m0Ks;
   TBranch        *b_xKs;
   TBranch        *b_yKs;
   TBranch        *b_zKs;
   TBranch        *b_s2xKs;
   TBranch        *b_s2yKs;
   TBranch        *b_s2zKs;
   TBranch        *b_chi2Ks;
   TBranch        *b_dofKs;
   TBranch        *b_stKs;
   TBranch        *b_ndauKs;
   TBranch        *b_MCKs;
   TBranch        *b_d1KsIndex;
   TBranch        *b_d1KsLund;
   TBranch        *b_d2KsIndex;
   TBranch        *b_d2KsLund;
   TBranch        *b_nPi0;
   TBranch        *b_massPi0;
   TBranch        *b_pPi0;
   TBranch        *b_thPi0;
   TBranch        *b_phiPi0;
   TBranch        *b_errMassPi0;
   TBranch        *b_m0Pi0;
   TBranch        *b_xPi0;
   TBranch        *b_yPi0;
   TBranch        *b_zPi0;
   TBranch        *b_s2xPi0;
   TBranch        *b_s2yPi0;
   TBranch        *b_s2zPi0;
   TBranch        *b_chi2Pi0;
   TBranch        *b_dofPi0;
   TBranch        *b_stPi0;
   TBranch        *b_ndauPi0;
   TBranch        *b_MCPi0;
   TBranch        *b_d1Pi0Index;
   TBranch        *b_d1Pi0Lund;
   TBranch        *b_d2Pi0Index;
   TBranch        *b_d2Pi0Lund;
   TBranch        *b_nGConv;
   TBranch        *b_massGConv;
   TBranch        *b_pGConv;
   TBranch        *b_thGConv;
   TBranch        *b_phiGConv;
   TBranch        *b_errMassGConv;
   TBranch        *b_m0GConv;
   TBranch        *b_xGConv;
   TBranch        *b_yGConv;
   TBranch        *b_zGConv;
   TBranch        *b_s2xGConv;
   TBranch        *b_s2yGConv;
   TBranch        *b_s2zGConv;
   TBranch        *b_chi2GConv;
   TBranch        *b_dofGConv;
   TBranch        *b_stGConv;
   TBranch        *b_ndauGConv;
   TBranch        *b_MCGConv;
   TBranch        *b_d1GConvIndex;
   TBranch        *b_d1GConvLund;
   TBranch        *b_d2GConvIndex;
   TBranch        *b_d2GConvLund;
   TBranch        *b_nDalitz;
   TBranch        *b_massDalitz;
   TBranch        *b_pDalitz;
   TBranch        *b_thDalitz;
   TBranch        *b_phiDalitz;
   TBranch        *b_errMassDalitz;
   TBranch        *b_m0Dalitz;
   TBranch        *b_xDalitz;
   TBranch        *b_yDalitz;
   TBranch        *b_zDalitz;
   TBranch        *b_s2xDalitz;
   TBranch        *b_s2yDalitz;
   TBranch        *b_s2zDalitz;
   TBranch        *b_chi2Dalitz;
   TBranch        *b_dofDalitz;
   TBranch        *b_stDalitz;
   TBranch        *b_ndauDalitz;
   TBranch        *b_MCDalitz;
   TBranch        *b_d1DalitzIndex;
   TBranch        *b_d1DalitzLund;
   TBranch        *b_d2DalitzIndex;
   TBranch        *b_d2DalitzLund;
   TBranch        *b_nJpsi;
   TBranch        *b_massJpsi;
   TBranch        *b_pJpsi;
   TBranch        *b_thJpsi;
   TBranch        *b_phiJpsi;
   TBranch        *b_errMassJpsi;
   TBranch        *b_m0Jpsi;
   TBranch        *b_xJpsi;
   TBranch        *b_yJpsi;
   TBranch        *b_zJpsi;
   TBranch        *b_s2xJpsi;
   TBranch        *b_s2yJpsi;
   TBranch        *b_s2zJpsi;
   TBranch        *b_chi2Jpsi;
   TBranch        *b_dofJpsi;
   TBranch        *b_stJpsi;
   TBranch        *b_ndauJpsi;
   TBranch        *b_MCJpsi;
   TBranch        *b_d1JpsiIndex;
   TBranch        *b_d1JpsiLund;
   TBranch        *b_d2JpsiIndex;
   TBranch        *b_d2JpsiLund;
   TBranch        *b_nTrk;
   TBranch        *b_IfrLayTrk;
   TBranch        *b_IfrNsTrk;
   TBranch        *b_IfrInnerTrk;
   TBranch        *b_IfrBarrelTrk;
   TBranch        *b_IfrFWDTrk;
   TBranch        *b_IfrBWDTrk;
   TBranch        *b_IfrMeasIntLenTrk;
   TBranch        *b_IfrFirstHitTrk;
   TBranch        *b_IfrLastHitTrk;
   TBranch        *b_lMomTrk;
   TBranch        *b_ZMom42Trk;
   TBranch        *b_ecalTrk;
   TBranch        *b_ecalXTrk;
   TBranch        *b_ecalYTrk;
   TBranch        *b_ecalZTrk;
   TBranch        *b_nCryTrk;
   TBranch        *b_nBumpTrk;
   TBranch        *b_ZMom20Trk;
   TBranch        *b_secMomTrk;
   TBranch        *b_s1s9Trk;
   TBranch        *b_s9s25Trk;
   TBranch        *b_erawTrk;
   TBranch        *b_phiClusterTrk;
   TBranch        *b_thetaClusterTrk;
   TBranch        *b_covEETrk;
   TBranch        *b_covTTTrk;
   TBranch        *b_covPPTrk;
   TBranch        *b_covRRTrk;
   TBranch        *b_phicMatTrk;
   TBranch        *b_trkcMatTrk;
   TBranch        *b_nPidTrk;
   TBranch        *b_emcStatusTrk;
   TBranch        *b_phiAtEMCTrk;
   TBranch        *b_thetaAtEMCTrk;
   TBranch        *b_isvtTrk;
   TBranch        *b_nsvtTrk;
   TBranch        *b_fhitTrk;
   TBranch        *b_ndchTrk;
   TBranch        *b_lhitTrk;
   TBranch        *b_tLenTrk;
   TBranch        *b_ntdofTrk;
   TBranch        *b_tproTrk;
   TBranch        *b_tChi2Trk;
   TBranch        *b_cPidTrk;
   TBranch        *b_sfRangeTrk;
   TBranch        *b_trkFitStatusTrk;
   TBranch        *b_chargeTrk;
   TBranch        *b_momentumTrk;
   TBranch        *b_ppcov00;
   TBranch        *b_ppcov10;
   TBranch        *b_ppcov11;
   TBranch        *b_ppcov20;
   TBranch        *b_ppcov21;
   TBranch        *b_ppcov22;
   TBranch        *b_xPocaTrk;
   TBranch        *b_yPocaTrk;
   TBranch        *b_zPocaTrk;
   TBranch        *b_thetaTrk;
   TBranch        *b_phiTrk;
   TBranch        *b_muonIdTrk;
   TBranch        *b_elecIdTrk;
   TBranch        *b_kaonIdTrk;
   TBranch        *b_pionIdTrk;
   TBranch        *b_idTrk;
   TBranch        *b_IndexTrk;
   TBranch        *b_IndexNtTrk;
   TBranch        *b_B0RecTrk;
   TBranch        *b_chBRecTrk;
   TBranch        *b_nGam;
   TBranch        *b_IfrLayGam;
   TBranch        *b_IfrNsGam;
   TBranch        *b_IfrInnerGam;
   TBranch        *b_IfrBarrelGam;
   TBranch        *b_IfrFWDGam;
   TBranch        *b_IfrBWDGam;
   TBranch        *b_IfrMeasIntLenGam;
   TBranch        *b_IfrFirstHitGam;
   TBranch        *b_IfrLastHitGam;
   TBranch        *b_IfrExpIntLenGam;
   TBranch        *b_IfrIntLenBeforeIronGam;
   TBranch        *b_IfrTrkMatchGam;
   TBranch        *b_IfrEmcMatchGam;
   TBranch        *b_IfrLastBarrelGam;
   TBranch        *b_IfrCLFitChi2Gam;
   TBranch        *b_IfrStrips0;
   TBranch        *b_IfrStrips1;
   TBranch        *b_IfrStrips2;
   TBranch        *b_IfrStrips3;
   TBranch        *b_IfrStrips4;
   TBranch        *b_IfrStrips5;
   TBranch        *b_IfrStrips6;
   TBranch        *b_IfrStrips7;
   TBranch        *b_IfrStrips8;
   TBranch        *b_IfrStrips9;
   TBranch        *b_IfrStrips10;
   TBranch        *b_IfrStrips11;
   TBranch        *b_IfrStrips12;
   TBranch        *b_IfrStrips13;
   TBranch        *b_IfrStrips14;
   TBranch        *b_IfrStrips15;
   TBranch        *b_IfrStrips16;
   TBranch        *b_IfrStrips17;
   TBranch        *b_IfrStrips18;
   TBranch        *b_IfrStrips19;
   TBranch        *b_lMomGam;
   TBranch        *b_ZMom42Gam;
   TBranch        *b_ecalGam;
   TBranch        *b_ecalXGam;
   TBranch        *b_ecalYGam;
   TBranch        *b_ecalZGam;
   TBranch        *b_nCryGam;
   TBranch        *b_nBumpGam;
   TBranch        *b_ZMom20Gam;
   TBranch        *b_secMomGam;
   TBranch        *b_s1s9Gam;
   TBranch        *b_s9s25Gam;
   TBranch        *b_erawGam;
   TBranch        *b_phiClusterGam;
   TBranch        *b_thetaClusterGam;
   TBranch        *b_covEEGam;
   TBranch        *b_covTTGam;
   TBranch        *b_covPPGam;
   TBranch        *b_covRRGam;
   TBranch        *b_emcStatusGam;
   TBranch        *b_thetaGam;
   TBranch        *b_phiGam;
   TBranch        *b_energyGam;
   TBranch        *b_idGam;
   TBranch        *b_IndexGam;
   TBranch        *b_IndexNtGam;
   TBranch        *b_B0RecGam;
   TBranch        *b_chBRecGam;


public :

  baseClass(TTree *tree=0,int isMC =0, int newFormat = 2);
  virtual ~baseClass();
  baseClass(TString, TString, int isMC=0, int newFormat = 2);

  virtual Int_t  Cut(Int_t entry);
  virtual void   Loop(Int_t maxEvent = 0, Int_t startEvent = 0, Int_t isVerbose = 0, Int_t lun = 0);
  virtual Int_t  GetEntry(Int_t entry);
  virtual Int_t  LoadTree(Int_t entry);
  virtual void   Init(TTree *tree, int isMC);
  virtual Bool_t Notify();
  virtual void   Show(Int_t entry = -1);

  TFile*    openHistFile(TString name);
  virtual void      closeHistFile();
  virtual void      bookHist(int dump = 0);
  virtual void      switchOffReading(const char *s);
  virtual void      redoPidKilling(int i) {fOptPidKilling = i;}
  virtual void      redoPidKillingKaon(int i) {fOptPidKillingKaon = i;}
  virtual void      redoPidKillingEl(int i) {fOptPidKillingEl = i;}
  virtual void      redoPidKillingMu(int i) {fOptPidKillingMu = i;}
  virtual void      runGammas(int i) {fOptGammas = i;}
  virtual void      doTrackSmearing(int i) {fOptSmearTracks = i;}
  virtual void      doNeutSmearing(int i) {fOptSmearNeut = i;}
  virtual void      filterK0s(int i) {fOptFilterK0sPi0Pi0 = i;}
  virtual void      runCategories(int i) {fOptCategories = i;}
  virtual void      runBlind(int i) {fOptBlind = i;}

  virtual void      initRest();
  virtual void      initVariables();
  virtual void      readCuts(TString filename, int dump = 1);
  virtual void      getPidTables();
  virtual void      getTrkTables();
  virtual void      dumpCuts();
  void      readweightchg();
  void      readweightneu();
  Double_t  getweightchg(double mom);   
  Double_t  getweightneu(double energy);   
  virtual Bool_t    isAncestor(int ancestor, int candidate); 
  virtual Int_t     isRecoed(int imc); 
  virtual void      readintpur();
  virtual void      timestamp(const char * lun);
  virtual void      fillPidMaps();
  virtual int       skipBadBreco();  
  Int_t compChgBreco(int nBs);
  virtual void      smearTracks(); 
  virtual void      smearNeut(); 
  virtual void      selectTracks(); 
  virtual void      selectPhotons(); 
  virtual void      breco();
  virtual void      findbestB();
  virtual void      calculateEvtW8(); 
  virtual double    getBsysweight(int decType,int thevub); 
  virtual double    getDsysweight(int decDpi,int decDk,int decDks,int decDpiz,int decDlep,int decImode,int thevub); 
  double clusterReCorrection(const double oldCorrectedEnergy, 
			     const double clusterPositionTheta) const ;
  double clusterCorrection(const double rawEnergy,
			   const double clusterPositionTheta, 
			   const bool   newCorr) const;
  virtual double KLlikeEMC(int iCand)const;
  virtual double KLlikeIFR(int iCand)const;
  virtual void      cleanGoodTracks(int what = 7); 
  virtual void      recoil();
  virtual void      doSplitOffStudy(bool book=false);
  virtual void      mcTruth();
  virtual void      dumpOneB(int b1 = -1); 
  virtual void      dumpGeneratorBlock(int b1 = -1, int b2 = -1); // prints overlap information if indices to the two B are passed

  virtual Double_t  fermi(double kp, double m, double a);
  virtual Double_t  fw8(double kp, double m, double a);

  virtual void mk4Vector(TLorentzVector &p4, const Float_t p, const Float_t t, const Float_t f, const Float_t m); 
  virtual void mk3Vector(TVector3 &p3, const Float_t p, const Float_t t, const Float_t f); 

  virtual Bool_t    isPrompt(Int_t iTrk);
  virtual Bool_t    isCascade(Int_t iTrk);
  virtual void      printLorentz(const TLorentzVector &);
  Int_t     bestKsIndex(Int_t isVerbose);

  virtual Double_t  getTrkW8(int i);

  virtual Bool_t    isTruLepton(int i);
  virtual Bool_t    isTruEl(int i);
  virtual Bool_t    isTruMu(int i);
  virtual Bool_t    isTruTau(int i);

  virtual Bool_t    isRecEl(int i);    
  virtual Bool_t    isRecMu(int i);    
  virtual Bool_t    isRecLepton(int i);
  virtual Bool_t    isRecKaon(int i);

  virtual Bool_t    isPidKillEl(Int_t i);
  virtual Bool_t    isPidKillMu(Int_t i);
  virtual Bool_t    isPidKillKaon(Int_t i);

  virtual Bool_t    isMcEl(int i);    
  virtual Bool_t    isMcMu(int i);    
  virtual Bool_t    isMcLepton(int i);
  virtual Bool_t    isMcKaon(int i);

  virtual Bool_t    isTrkOfBreco(int i);
  virtual Bool_t    isGamOfBreco(int i);
};

inline Bool_t baseClass::isTrkOfBreco(int i) {
  Int_t Brectrktmp;
  Brectrktmp = B0RecTrk[i];
  if(fChB) Brectrktmp = chBRecTrk[i];
  if ((Brectrktmp&brecoOverlap))  {return kTRUE;} else {return kFALSE;};
}

inline Bool_t baseClass::isGamOfBreco(int i) {
  Int_t Brecgamtmp;
  Brecgamtmp = B0RecGam[i];
  if(fChB) Brecgamtmp = chBRecGam[i];
  if (Brecgamtmp&brecoOverlap)  {return kTRUE;} else {return kFALSE;}
}

inline Bool_t baseClass::isPrompt(int i) {
  if (!fIsMC) return kFALSE; 
  if (IndexTrk[i]-1 > nMc) return kFALSE;
  int lundMom = idMc[mothMc[IndexTrk[i]-1]-1];
  return ((TMath::Abs(lundMom) == 511) || (TMath::Abs(lundMom) == 521));
} 

inline Bool_t baseClass::isCascade(int i) {
  if (!fIsMC) return kFALSE; 
  if (IndexTrk[i]-1 > nMc) return kFALSE;
  int lundMom = idMc[mothMc[IndexTrk[i]-1]-1];
  int did = ((lundMom - (lundMom%100))/100)%10;
  return (TMath::Abs(did) == 4);
} 

// ----------------------------------------------------------------------
inline void baseClass::mk4Vector(TLorentzVector &p4, const Float_t p, const Float_t t, const Float_t f, const Float_t m) {
  p4.SetXYZM(p*TMath::Sin(t)*TMath::Cos(f), p*TMath::Sin(t)*TMath::Sin(f), p*TMath::Cos(t), m);
}

inline void baseClass::mk3Vector(TVector3 &p3, const Float_t p, const Float_t t, const Float_t f) {
  p3.SetXYZ(p*TMath::Sin(t)*TMath::Cos(f), p*TMath::Sin(t)*TMath::Sin(f), p*TMath::Cos(t));
}



// -- PID options
// --------------
// Note that so far nor PID killing has been run on MC ...
// http://www.slac.stanford.edu/BFROOT/www/Physics/Analysis/AWG/Charmonium/CharmUser/bitmap.html
// xxxxIdTrk & 1  -> noCal/MIP
// xxxxIdTrk & 2  -> veryLoose
// xxxxIdTrk & 4  -> loose
// xxxxIdTrk & 8  -> tight
// xxxxIdTrk & 16 -> veryTight
// xxxxIdTrk & 32 -> likelihood (to be added, especially for electrons)
 
// -- PID with reco selectors 
inline Bool_t baseClass::isRecEl(Int_t i)     {if (recEl[i] == 1) {return kTRUE;} else {return kFALSE;} }
inline Bool_t baseClass::isRecMu(Int_t i)     {if (recMu[i] == 1) {return kTRUE;} else {return kFALSE;} }
inline Bool_t baseClass::isRecLepton(Int_t i) {return (isRecEl(i) || isRecMu(i));}
inline Bool_t baseClass::isRecKaon(Int_t i)   {if (recKa[i] == 1) {return kTRUE;} else {return kFALSE;} }

// -- PID with MC truth  matched information
inline Bool_t baseClass::isMcEl(Int_t i)     {return (TMath::Abs(idTrk[i]) == 11);}
inline Bool_t baseClass::isMcMu(Int_t i)     {return (TMath::Abs(idTrk[i]) == 13);}
inline Bool_t baseClass::isMcLepton(Int_t i) {return (isMcEl(i) || isMcMu(i));}
inline Bool_t baseClass::isMcKaon(Int_t i) {return (TMath::Abs(idTrk[i]) == 321);}

// -- PID for generator block
inline Bool_t baseClass::isTruEl(Int_t i)  {return (TMath::Abs(idMc[i]) == 11);}
inline Bool_t baseClass::isTruMu(Int_t i)  {return (TMath::Abs(idMc[i]) == 13);}
inline Bool_t baseClass::isTruTau(Int_t i) {return (TMath::Abs(idMc[i]) == 15);}
inline Bool_t baseClass::isTruLepton(Int_t i) {return (isTruEl(i) || isTruMu(i));}

// -- Shape function
inline Double_t baseClass::fermi(double kp, double m, double a) {
  // FIXME: This probably messes up the normalization
  if (kp<(BMASS-m)) {
    return TMath::Power((1 - (kp/(BMASS-m))), a) * TMath::Exp((1+a) * (kp/(BMASS-m))); 
  } 
  return 1.;
}
inline Double_t baseClass::fw8(double kp, double mb, double a) {
  return fermi(kp, mb, a) / fermi(kp, BQMASS, A0);
}

// ----------------------------------------------------------------------
// -- HV Ranges according to 
// http://www.slac.stanford.edu/BFROOT/www/Physics/TrackEfficTaskForce/Recipe/GoodTracksLoose-Run2.html
// ----------------------------------------------------
// HV   DATA            SP3              SP4
// 1900 9931   14461    300000 399999    600000 679999
// 1960 14462  17106    400000 499999    680000 749999
// 1930 18190  ...                       770000 1009999
// ----------------------------------------------------
inline Double_t baseClass::getTrkW8(int i) {
  double w8(1.);
  double pt = momentumTrk[i]*TMath::Sin(thetaTrk[i]);
  if ((fRunnumber >= 300000) && (fRunnumber <= 399999)) {
    // -- SP3 at 1900 V
    w8 = fTT[0]->w8R(pt, nTrk, thetaTrk[i], phiTrk[i]);
  } else if ((fRunnumber >= 400000) && (fRunnumber <= 499999)) {
    // -- SP3 at 1960 V
    w8 = fTT[1]->w8R(pt, nTrk, thetaTrk[i], phiTrk[i]);
  } else if ((fRunnumber >= 600000) && (fRunnumber <= 684999)) {
    // -- SP4 at 1900 V
    w8 = fTT[2]->w8R(pt, nTrk, thetaTrk[i], phiTrk[i]);
  } else if ((fRunnumber >= 680000) && (fRunnumber <= 749999)) {
    // -- SP4 at 1960 V
    w8 = fTT[3]->w8R(pt, nTrk, thetaTrk[i], phiTrk[i]);
  } else if ((fRunnumber >= 770000) && (fRunnumber <= 1009999)) {
    // -- SP4 at 1930 V
    w8 = fTT[4]->w8R(pt, nTrk, thetaTrk[i], phiTrk[i]);
  }
  //  cout << fRunnumber << "  " << pt << "  " << nTrk << "  " << thetaTrk[i] << "  " << phiTrk[i] << " -> w8 = " << w8 << endl;
  return w8;
}


#endif
