#include "recoilBase.hh"

#include <fstream>
#include <iostream>
#include <iomanip>

#include "TH2.h"
#include "TH1.h"
#include <TMap.h>
#include <TObjString.h>
#include "TLorentzVector.h"

using namespace std;

//-------------------------------------------------
//FORTRAN STUFF
//-------------------------------------------------
#if 0
extern "C" 
{
  // The fortran generator:
  int BILTYP;
  float BCVAL[4];
  float BP_REC[16];
  float BP_FIT[16];
  float BCHI2T, BPROBCHI2;
  int   BIERR, BISMEAR, BISV;
  int abcfit_interface_vub_(int *BISMEAR, int *BILTYP, float *BCVAL, float *BP_REC, float *BP_FIT, float *BCHI2T,float *BPROBCHI2,int *BIERR, int *BISV);
}
#endif

// ----------------------------------------------------------------------
void  recoilBase::printLorentz(const TLorentzVector &p) {
  char line[200]; 
  sprintf(line, "p = (%7.5f, %7.5f, %7.5f, %7.5f), m = %7.5f",  p.X(), p.Y(), p.Z(), p.E(), p.Mag()); 
  cout << line;
}


// ----------------------------------------------------------------------
recoilBase::recoilBase(TTree *tree, int isMC, int mcasdata, int doInit, int newFormat) {
  fNewFormat = newFormat; 
  cout <<"recoilBase constructor tree isMC"<<isMC<<endl;
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/u/ec/ursl/d/output/breco-110201/mc-b0cocktail-20.root");
    if (!f) {
      f = new TFile("/u/ec/ursl/d/output/breco-110201/mc-b0cocktail-20.root");
    }
    tree = (TTree*)gDirectory->Get("h9");
    
  }
  if (doInit) {
    Init(tree,isMC,mcasdata);
    fToBeCopied = new TEventList("toBeCopied", "Events to be copied", 1000);
    //ADDED CB
    initRest();
  }
}

// ----------------------------------------------------------------------
recoilBase::~recoilBase() {
  cout <<"recoilBase destructor"<<endl;
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}


// ----------------------------------------------------------------------
recoilBase::recoilBase(TString filename, TString treename,int isMC, int mcasdata, int newFormat) {
  fNewFormat = newFormat; 
  cout<<"recoilBase constructor filename treename isMC"<<isMC<<endl;
  TFile *f = new TFile(filename);
  TTree *tree = (TTree*)f->Get(treename);
  if (!tree) { 
    cout << "Did not find " << treename << " in file " << filename << endl;
    f->ls();
  } else {
    Init(tree,isMC,mcasdata);
  }
  initRest();
}

// ---------------------------------------------------------------------- 

  void recoilBase::initRest() {
  DR = 57.29577951;
  PCBIN = 30;
  PCMAX = 3.0;
  PLBIN = 40;
  PLMAX = 4.0;
  FBIN = 36;
  B0LUND = 511;
  CHBLUND = 521;

  XMAX = 5.;
  XBIN = 50;

  ELMASS  = 0.000511;
  MUMASS  = 0.10567;

  PIPMASS = 0.13957; 
  PIZMASS = 0.13498;
  KAPMASS = 0.49368;
  KAZMASS = 0.49767;
  BPMASS  = 5.2790;
  BZMASS  = 5.2794;

  fOptBlind = fOptGammas = fOptCategories = fOptKlongs = fOptMakeEventList = 0;

  // -- CUTS
  DESIGNALLO = -0.1; 
  DESIGNALHI =  0.1; 

  MESSIGNALLO   =  5.27; 
  MESSIGNALHI   =  5.29; 

  MESSIDEBANDLO = 5.20; 
  MESSIDEBANDHI = 5.26; 

  MESSIGNALBANDLO = 5.27; 
  MESSIGNALBANDHI = 5.29; 

  PURITY     = 0.;
  INTPURITY  = 0.;
  IPURDSTAR  = 0.;
  IPURDC     = 0.;
  IPURDSTAR0 = 0.;
  IPURD0     = 0.;

  TLABLO = 20.6;  // BAD 90 
  TLABHI = 135.9;

  PCMSLO = 1.0;
  PLABLO = 0.5;
  NLEPTON = 1;
  ELMOMLO = 1.0;
  MUMOMLO = 1.0;
  //  IDEL = IDMU = IDKA = 16; 
  IDEL = IDMU = 16; 
  IDKA = 8; 

  MM2LO = -9999.;
  MM2HI = 0.5;

  REQCHARGECORR  = kTRUE;
  REQTOTALCHARGE = kTRUE;

  KSPIPLO = KSPIZLO = 0.2;
  KSPIPHI = KSPIZHI = 0.8;

  GAMMAELO = 0.080;
  GAMMAEHI = 100.0;
  
  fVerbose = 0;
  fHistFile = 0;
  fDump = 0;

  Int_t i(0);
  for (i = 0; i < 100; ++i) { 
    kshortLockTrk[i] = 0; 
    kshortLockGam[i] = 0; 
    goodKshort[i] = 0;
  }

  READALL=1;
}

// ----------------------------------------------------------------------
  void recoilBase::readCuts(TString filename, int dump) {
  char  buffer[200];
  fCutFile = filename;
  sprintf(buffer, "%s", filename.Data());
  ifstream is(buffer);
  char CutName[100];
  float CutValue;
  int ok(0);
  while (is.getline(buffer, 200, '\n')) {
    ok = 0;
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);
    // -- breco 
    if (!strcmp(CutName, "mesSignalLo")) {MESSIGNALLO   = CutValue; ok = 1;}
    if (!strcmp(CutName, "mesSignalHi")) {MESSIGNALHI   = CutValue; ok = 1;}
    if (!strcmp(CutName, "mesSidebandLo")) {MESSIDEBANDLO = CutValue; ok = 1;}
    if (!strcmp(CutName, "mesSidebandHi")) {MESSIDEBANDHI = CutValue; ok = 1;}
    if (!strcmp(CutName, "mesSignalbandLo")) {MESSIGNALBANDLO = CutValue; ok = 1;}
    if (!strcmp(CutName, "mesSignalbandHi")) {MESSIGNALBANDHI = CutValue; ok = 1;}
    if (!strcmp(CutName, "deSignalLo")) {DESIGNALLO = CutValue; ok = 1;}
    if (!strcmp(CutName, "deSignalHi")) {DESIGNALHI = CutValue; ok = 1;}
    if (!strcmp(CutName, "purity")) {PURITY = CutValue; ok = 1;}
    if (!strcmp(CutName, "intPurity")) {INTPURITY = CutValue; ok = 1;}
    if (!strcmp(CutName, "ipurDstar"))  {IPURDSTAR = CutValue; ok = 1;}
    if (!strcmp(CutName, "ipurDc"))     {IPURDC = CutValue; ok = 1;}
    if (!strcmp(CutName, "ipurDstar0")) {IPURDSTAR0 = CutValue; ok = 1;}
    if (!strcmp(CutName, "ipurD0"))     {IPURD0 = CutValue; ok = 1;}

    // -- Lepton 
    if (!strcmp(CutName, "pcmsLo")) {PCMSLO = CutValue; ok = 1;}
    if (!strcmp(CutName, "plabLo")) {PLABLO = CutValue; ok = 1;}
    if (!strcmp(CutName, "tlabLo")) {TLABLO = CutValue; ok = 1;}
    if (!strcmp(CutName, "tlabHi")) {TLABHI = CutValue; ok = 1;}
    if (!strcmp(CutName, "nLepton")) {NLEPTON = int(CutValue); ok = 1;}
    if (!strcmp(CutName, "elmomLo")) {ELMOMLO = CutValue; ok = 1;}
    if (!strcmp(CutName, "mumomLo")) {MUMOMLO = CutValue; ok = 1;}
    if (!strcmp(CutName, "idEl")) {IDEL = int(CutValue); ok = 1;}
    if (!strcmp(CutName, "idMu")) {IDMU = int(CutValue); ok = 1;}
    if (!strcmp(CutName, "idKa")) {IDKA = int(CutValue); ok = 1;}

    // -- recoil
    if (!strcmp(CutName, "mm2Lo")) {MM2LO = CutValue; ok = 1;}
    if (!strcmp(CutName, "mm2Hi")) {MM2HI = CutValue; ok = 1;}
    if (!strcmp(CutName, "reqChargeCoor")) {REQCHARGECORR = (CutValue > 0.5) ?  kTRUE: kFALSE; ok = 1;}
    if (!strcmp(CutName, "reqTotalCharge")) {REQTOTALCHARGE = (CutValue > 0.5) ?  kTRUE: kFALSE; ok = 1;}
    if (!strcmp(CutName, "kspipLo")) {KSPIPLO = CutValue; ok = 1;}
    if (!strcmp(CutName, "kspipHi")) {KSPIPHI = CutValue; ok = 1;}
    if (!strcmp(CutName, "kspizLo")) {KSPIZLO = CutValue; ok = 1;}
    if (!strcmp(CutName, "kspizHi")) {KSPIZHI = CutValue; ok = 1;}
    if (!strcmp(CutName, "gammaELo")) {GAMMAELO = CutValue; ok = 1;}
    if (!strcmp(CutName, "gammaEHi")) {GAMMAEHI = CutValue; ok = 1;}

    if (ok == 0)  cout << "==> recoilBase::readCuts() Error: Don't know about variable " << CutName << endl;
  }

  if (dump == 1) dumpCuts();

}


// ----------------------------------------------------------------------
Bool_t recoilBase::isAncestor(int ancestor, int cand) {
  // ancestor is in C++ standard, e.g. counting from 0
  // cand is in Fortran standard, i.e. counting from 1
  int mom(cand); 
  while (mom > 0) {
    //    cout<<TMath::Abs(idMc[mom])<<endl;
    mom = mothMc[mom]-1; 
    //    cout << "   now looking at " << mom << " which is a " << idMc[mom] << endl;
    if (mom == ancestor) return true;
  }
  return false;
}

// ----------------------------------------------------------------------
Int_t recoilBase::isRecoed(int imc) {
  int result(-1), i(0); 
  for (i = 0; i < nTrk; ++i) {
    if (imc == (IndexTrk[i]-1)) {
      result = i;
      break;
    }
  }
  if (result > -1) return result;
  for (i = 0; i < nGam; ++i) {
    if (imc == (IndexGam[i]-1)) {
      result = i;
      break;
    }
  }
  return result;
}


// ----------------------------------------------------------------------
void recoilBase::dumpCuts() {
  cout << "====================================" << endl;
  cout << "Cut file " << fCutFile << endl; 
  cout << "------------------------------------" << endl;
  cout << "mesSignal:      " << MESSIGNALLO << " ... " << MESSIGNALHI << endl;
  cout << "mesSideband:    " << MESSIDEBANDLO << " ... " << MESSIDEBANDHI << endl;
  cout << "mesSignalband:  " << MESSIGNALBANDLO << " ... " << MESSIGNALBANDHI << endl;
  cout << "deSignal:       " << DESIGNALLO << " ... " << DESIGNALHI << endl;
  cout << "purity:         " << PURITY << endl;
  cout << "intPurity:      " << INTPURITY << endl;
  cout << "ipurDstar:      " << IPURDSTAR << endl;
  cout << "ipurDc:         " << IPURDC << endl;
  cout << "ipurDstar0:     " << IPURDSTAR0 << endl;
  cout << "ipurD0:         " << IPURD0 << endl;

  cout << "pcmsLo:         " << PCMSLO << endl;
  cout << "plabLo:         " << PLABLO << endl;
  cout << "tlabLo:         " << TLABLO << endl;
  cout << "tlabHi:         " << TLABHI << endl;
  cout << "nLepton:        " << NLEPTON << endl;
  cout << "elmomLo:        " << ELMOMLO << endl;
  cout << "mumomLo:        " << MUMOMLO << endl;
  cout << "idEl:           " << IDEL << endl;
  cout << "idMu:           " << IDMU << endl;
  cout << "idKa:           " << IDKA << endl;

  cout << "mm2Lo:          " << MM2LO << endl;
  cout << "mm2Hi:          " << MM2HI << endl;
  cout << "reqChargeCoor:  " << (REQCHARGECORR ?  "true": "false") << endl;
  cout << "reqTotalCharge: " << (REQTOTALCHARGE ? "true": "false") << endl;
  cout << "kspip:          " << KSPIPLO << " ... " << KSPIPHI << endl;
  cout << "kspiz:          " << KSPIZLO << " ... " << KSPIZHI << endl;
  cout << "GammaE:         " << GAMMAELO << " ... " << GAMMAEHI << endl;
  cout << "====================================" << endl;
}


// ----------------------------------------------------------------------
void recoilBase::breco(int chbcand, int isVer) {

  static Bool_t first(kTRUE); 
  char name[100], title[100];
  TH1D *h;
  TH2D *h2;

  if (first) {
    first = kFALSE;
    
    fHistFile->cd();
    fHistFile->mkdir("breco", "breco");
    fHistFile->cd("breco");

    sprintf(name, "h100");  sprintf(title, "Breco modes");  h = new TH1D(name, title, 20000, 0., 20000.); 
    sprintf(name, "f100");  sprintf(title, "Breco modes with purity <=0");  h = new TH1D(name, title, 20000, 0., 20000.); 
    sprintf(name, "h101");  sprintf(title, "purity B0");  h = new TH1D(name, title, 500, 0., 1.); 
    sprintf(name, "h102");  sprintf(title, "intpurity B0");  h = new TH1D(name, title, 500, 0., 100.); 
    sprintf(name, "h111");  sprintf(title, "purity B+");  h = new TH1D(name, title, 500, 0., 1.); 
    sprintf(name, "h112");  sprintf(title, "intpurity B+");  h = new TH1D(name, title, 500, 0., 100.); 
    
    sprintf(name, "h200");  sprintf(title, "mes");  h = new TH1D(name, title, 40, 5.2, 5.3); 
    sprintf(name, "h201");  sprintf(title, "mes, new");  h = new TH1D(name, title, 80, 5.1, 5.3); 
    sprintf(name, "h202");  sprintf(title, "mes, old");  h = new TH1D(name, title, 80, 5.1, 5.3); 
    sprintf(name, "h300");  sprintf(title, "mes dE all");  h2 = new TH2D(name, title, 40, 5.2, 5.3, 40, -0.2, 0.2); 

    sprintf(name, "h400");  sprintf(title, "p B, new");  h = new TH1D(name, title, 50, 0., .5); 
    sprintf(name, "h401");  sprintf(title, "p B, old");  h = new TH1D(name, title, 50, 0., .5); 
    
    sprintf(name, "mesallDupli");  sprintf(title, "mes All Dupli");  h = new TH1D(name, title, 40, 5.2, 5.3); 
    sprintf(name, "deallDupli");  sprintf(title, "delta E All Dupli");  h = new TH1D(name, title, 40, -0.1, 0.1); 
    sprintf(name, "mesdeDupli");  sprintf(title, "mes dE all Dupli");  h2 = new TH2D(name, title, 40, 5.2, 5.3, 40, -0.2, 0.2); 

    int i; 
    for (i = 0; i < 600; ++i) {
      sprintf(name, "h%d", 11000+i);  sprintf(title, "mes, mode = %d", 11000+i);  h = new TH1D(name, title, 40, 5.2, 5.3); 
      sprintf(name, "h%d", 12000+i);  sprintf(title, "mes, mode = %d", 12000+i);  h = new TH1D(name, title, 40, 5.2, 5.3); 
      sprintf(name, "h%d", 13000+i);  sprintf(title, "mes, mode = %d", 13000+i);  h = new TH1D(name, title, 40, 5.2, 5.3); 
      sprintf(name, "h%d", 14000+i);  sprintf(title, "mes, mode = %d", 14000+i);  h = new TH1D(name, title, 40, 5.2, 5.3); 

      sprintf(name, "ae%d", 11000+i);  sprintf(title, "all events, mode = %d", 11000+i);  h = new TH1D(name, title, 40, 5.2, 5.3); 
      sprintf(name, "ae%d", 12000+i);  sprintf(title, "all events, mode = %d", 12000+i);  h = new TH1D(name, title, 40, 5.2, 5.3); 
      sprintf(name, "ae%d", 13000+i);  sprintf(title, "all events, mode = %d", 13000+i);  h = new TH1D(name, title, 40, 5.2, 5.3); 
      sprintf(name, "ae%d", 14000+i);  sprintf(title, "all events, mode = %d", 14000+i);  h = new TH1D(name, title, 40, 5.2, 5.3); 

      sprintf(name, "ac%d", 11000+i);  sprintf(title, "all cuts, mode = %d", 11000+i);  h = new TH1D(name, title, 40, 5.2, 5.3); 
      sprintf(name, "ac%d", 12000+i);  sprintf(title, "all cuts, mode = %d", 12000+i);  h = new TH1D(name, title, 40, 5.2, 5.3); 
      sprintf(name, "ac%d", 13000+i);  sprintf(title, "all cuts, mode = %d", 13000+i);  h = new TH1D(name, title, 40, 5.2, 5.3); 
      sprintf(name, "ac%d", 14000+i);  sprintf(title, "all cuts, mode = %d", 14000+i);  h = new TH1D(name, title, 40, 5.2, 5.3); 

    }

  }

  fHistFile->cd("breco");
  
  // -- Reco quantities
  TLorentzVector tmpp4Breco, tmpp4BrecoNC;
  int tmpdauB;
  int tmpmodeB;
  if(chbcand) {
    fBrecoMc = MCChB[indexbestB];
    fMes = mseChB[indexbestB];
    fDeltaE = deltaeChB[indexbestB];
    tmpdauB = d1ChBLund[indexbestB];
    tmpmodeB = modeChB[indexbestB];
    mk4Vector(tmpp4Breco, MassPChB[indexbestB], MassThetaChB[indexbestB], MassPhiChB[indexbestB], BPMASS); 
    mk4Vector(tmpp4BrecoNC, pChB[indexbestB], thChB[indexbestB], phiChB[indexbestB], massChB[indexbestB]); 
  } else {
    fBrecoMc = MCB0[indexbestB];
    fMes = mseB0[indexbestB];
    fDeltaE = deltaeB0[indexbestB];
    tmpdauB = d1B0Lund[indexbestB];
    tmpmodeB = modeB0[indexbestB];
    mk4Vector(tmpp4Breco, MassPB0[indexbestB], MassThetaB0[indexbestB], MassPhiB0[indexbestB], BZMASS); 
    mk4Vector(tmpp4BrecoNC, pB0[indexbestB], thB0[indexbestB], phiB0[indexbestB], massB0[indexbestB]); 
  }
  double pB_old=tmpp4BrecoNC.Vect().Mag();
  
  // this is used throughout
  p4Upsilon = TLorentzVector(pxUps, pyUps, pzUps, eUps); 
  upsBoost = TVector3(pxUps, pyUps, pzUps);
  upsBoost.SetMag(upsBoost.Mag()/eUps);
  
  tmpp4BrecoNC.Boost(-upsBoost); 
  double pB=tmpp4BrecoNC.Vect().Mag();
  // mes = sqrt(E*^2 - p*^2) 
  double E=p4Upsilon.Mag();
  double myMes=sqrt(E*E/4 - pB*pB);
  double oldMes=fMes;
  // cout << "Mes (E=" << E << ") : (calculated - ntuple = difference) - (" << myMes << " - " << fMes << " = " << myMes-fMes << endl;
  fMes=myMes;
  
  fBmode = tmpmodeB;


  // BetaCoreTools/BtaExclusiveDecayList.hh
  // the following is an arbitrary definition
  // fSeedMode = 0 dc
  // fSeedMode = 1 dstar
  // fSeedMode = 2 d0
  // fSeedMode = 3 dstar0
  if ((11000 <= tmpmodeB) &&  (tmpmodeB < 12000)) {
    fSeedMode = 3;
  } else if  ((12000 <= tmpmodeB) &&  (tmpmodeB < 13000)) {
    fSeedMode = 0;
  } else if  ((13000 <= tmpmodeB) &&  (tmpmodeB < 14000)) {
    fSeedMode = 1;
  } else if  ((14000 <= tmpmodeB) &&  (tmpmodeB < 15000)) {
    fSeedMode = 2;
  } else {
    fSeedMode = -1;
  }

  
  // PUTTING IN DANIELE's FIX for D0->Kspipi
  // docs: http://www.slac.stanford.edu/cgi-bin/lwgate/VUB-RECOIL/archives/vub-recoil.200207/Author/article-93.html
  // NOTE: fBrecoCharge is the charge of the BRECO candidate
  // ----- fBrecoFlavor is the sign of a prompt lepton if the BRECO decayed semileptonically
  // tmpdauB needs -1 in front because B+(521) -> anti-D0(-421) nomenclature
  fBrecoCharge= 0;
  fBrecoFlavor = 1;    

  if (chbcand) {
    if(d2ChBLund[indexbestB]!=0&&d2ChBLund[indexbestB]!=111&&d2ChBLund[indexbestB]!=310){ fBrecoCharge=fBrecoCharge+(TMath::Abs(d2ChBLund[indexbestB])/d2ChBLund[indexbestB]);}    
    if(d3ChBLund[indexbestB]!=0&&d3ChBLund[indexbestB]!=111&&d3ChBLund[indexbestB]!=310){ fBrecoCharge=fBrecoCharge+(TMath::Abs(d3ChBLund[indexbestB])/d3ChBLund[indexbestB]);}            
    if(d4ChBLund[indexbestB]!=0&&d4ChBLund[indexbestB]!=111&&d4ChBLund[indexbestB]!=310){ fBrecoCharge=fBrecoCharge+(TMath::Abs(d4ChBLund[indexbestB])/d4ChBLund[indexbestB]);}            
    if(d5ChBLund[indexbestB]!=0&&d5ChBLund[indexbestB]!=111&&d5ChBLund[indexbestB]!=310){ fBrecoCharge=fBrecoCharge+(TMath::Abs(d5ChBLund[indexbestB])/d5ChBLund[indexbestB]);}            
    if(d6ChBLund[indexbestB]!=0&&d6ChBLund[indexbestB]!=111&&d6ChBLund[indexbestB]!=310){ fBrecoCharge=fBrecoCharge+(TMath::Abs(d6ChBLund[indexbestB])/d6ChBLund[indexbestB]);}            
    if(d7ChBLund[indexbestB]!=0&&d7ChBLund[indexbestB]!=111&&d7ChBLund[indexbestB]!=310){ fBrecoCharge=fBrecoCharge+(TMath::Abs(d7ChBLund[indexbestB])/d7ChBLund[indexbestB]);}    
    fBrecoFlavor=fBrecoCharge;        
  } else {
    fBrecoCharge=0;
    fBrecoFlavor=-1*(TMath::Abs(d1B0Lund[indexbestB])/d1B0Lund[indexbestB]);
  }
  // END IMPLEMENTING DANIELE'S FIX


  if ((MESSIDEBANDLO < fMes) && (fMes < MESSIDEBANDHI)) {
    mesSideband = kTRUE; 
  } else {
    mesSideband = kFALSE; 
  }

  if ((MESSIGNALBANDLO < fMes) && (fMes < MESSIGNALBANDHI)) {
    mesSignalband = kTRUE; 
  } else {
    mesSignalband = kFALSE; 
  }

  if ((MESSIGNALLO < fMes) && (fMes < MESSIGNALHI)
      && (DESIGNALLO < fDeltaE) && (fDeltaE < DESIGNALHI)
      ) {
    signalBox = kTRUE;
  } else {
    signalBox = kFALSE;
  }

  if ((MESSIDEBANDLO < fMes) && (fMes < MESSIGNALHI)
      && (DESIGNALLO < fDeltaE) && (fDeltaE < DESIGNALHI)
      ) {
    tsdump = kTRUE;
  } else {
    tsdump = kFALSE;
  }



  fHistFile->cd("breco");

  ((TH1D*)gDirectory->Get("h100"))->Fill(tmpmodeB);
  ((TH1D*)gDirectory->Get("h200"))->Fill(fMes);
  ((TH1D*)gDirectory->Get("h201"))->Fill(fMes);
  ((TH1D*)gDirectory->Get("h202"))->Fill(oldMes);
  ((TH2D*)gDirectory->Get("h300"))->Fill(fMes, fDeltaE);

  ((TH1D*)gDirectory->Get("h400"))->Fill(pB);
  ((TH1D*)gDirectory->Get("h401"))->Fill(pB_old);

  if(fisDuplicate) {
    ((TH1D*)gDirectory->Get("deallDupli"))->Fill(fDeltaE);
    ((TH1D*)gDirectory->Get("mesallDupli"))->Fill(fMes);
    ((TH1D*)gDirectory->Get("mesdeDupli"))->Fill(fMes, fDeltaE);
  }

  sprintf(name, "h%d", fBmode);  fillHist(name, name, template_MES, fMes);

  double tmpIpurB, tmppurB;
  if(chbcand) {
    tmpIpurB = intpurChB[indexbestB];
    tmppurB = purChB[indexbestB];
    ((TH1D*)gDirectory->Get("h111"))->Fill(tmppurB);
    ((TH1D*)gDirectory->Get("h112"))->Fill(tmpIpurB);
  } else {
    tmpIpurB = intpurB0[indexbestB];
    tmppurB = purB0[indexbestB];
    ((TH1D*)gDirectory->Get("h101"))->Fill(tmppurB);
    ((TH1D*)gDirectory->Get("h102"))->Fill(tmpIpurB);
  }
  
#if 0
    int SMEARBRECO = 1;
    if (SMEARBRECO) {
      int BILTYP; 
      float BCVAL[4] = {0., 0., 0., 1110.};
      float BP_REC[16] = {
	p4BrecoGen.Px(), p4BrecoGen.Py(), p4BrecoGen.Pz(), p4BrecoGen.E(), 
	0., 0., 0., 0., 
	0., 0., 0., 0., 
	0., 0., 0., 0.};
      float BP_FIT[16];
      float BCHI2T, BPROBCHI2; 
      
      int   BISMEAR(1);
      int   BIERR;
      int   BISV = isVer;      
      int i1 = abcfit_interface_vub_(&BISMEAR,&BILTYP,BCVAL,BP_REC,BP_FIT,&BCHI2T,&BPROBCHI2,&BIERR, &BISV);

      p4BrecoGenSmear.SetXYZM(BP_FIT[0], BP_FIT[1], BP_FIT[2], BP_FIT[3]); 

      if (signalBox) {
	fHistFile->cd("sgallevents");
      } else {
	fHistFile->cd("bgallevents");
      }

      ((TH1D*)gDirectory->Get("s100"))->Fill(p4BrecoGen.Vect().Mag(), 1.);
      ((TH1D*)gDirectory->Get("s101"))->Fill(p4BrecoNC.Vect().Mag(), 1.);
      ((TH1D*)gDirectory->Get("s102"))->Fill(p4BrecoGenSmear.Vect().Mag(), 1.);
      ((TH1D*)gDirectory->Get("s103"))->Fill(p4BrecoNC.Vect().Mag() - p4BrecoGenSmear.Vect().Mag(), 1.);
      ((TH1D*)gDirectory->Get("s104"))->Fill(p4BrecoNC.Vect().Mag() - p4BrecoGen.Vect().Mag(), 1.);
      ((TH1D*)gDirectory->Get("s105"))->Fill(p4BrecoGen.Vect().Mag() - p4BrecoGenSmear.Vect().Mag(), 1.);

    }
#endif

  
  fHistFile->cd();
}  
  
// ----------------------------------------------------------------------
void recoilBase::testConversion() {

  static Bool_t first(kTRUE); 
  char name[100], title[100];
  TH1D *h;
  TH2D *h2;

  if (first) {
    first = kFALSE;
    
    fHistFile->cd();
    fHistFile->mkdir("conversion", "conversion");
    fHistFile->cd("conversion");
    sprintf(name, "h100");  sprintf(title, "mass");  h = new TH1D(name, title, 100, 0., 0.1); 
    sprintf(name, "h101");  sprintf(title, "mass0");  h = new TH1D(name, title, 100, 0., 0.1); 
    sprintf(name, "h102");  sprintf(title, "daughter ID");  h = new TH1D(name, title, 1000, -500., 500.); 
    sprintf(name, "h103");  sprintf(title, "radius");  h = new TH1D(name, title, 200, 0., 20.); 
    sprintf(name, "h104");  sprintf(title, "rec electrons");  h = new TH1D(name, title, 10, 0., 10.); 
    //    sprintf(name, "h105");  sprintf(title, "chi2");  h = new TH1D(name, title, 100, 0., 1.); 
    sprintf(name, "h200");  sprintf(title, "p1 vs p2");  h2 = new TH2D(name, title, 40, 0., 2., 40, 0., 2.); 
    sprintf(name, "h201");  sprintf(title, "x vs y ");  h2 = new TH2D(name, title, 100, -10., 10., 100, -10., 10.); 
  }

  fHistFile->cd("conversion");
  Int_t i;
  Double_t mass, m0, radius;
  Double_t p1, p2, x, y; 
  Int_t id1, id2, ip1, ip2; 
  for (i = 0; i < nGConv; ++i) {
    mass = massGConv[i];
    m0 = m0GConv[i];
    id1 = d1GConvLund[i]; 
    ip1 = d1GConvIndex[i]; 
    id2 = d2GConvLund[i]; 
    ip2 = d2GConvIndex[i]; 
    p1 = momentumTrk[ip1];
    p2 = momentumTrk[ip2];
    x = xGConv[i];
    y = yGConv[i]; 
    radius = TMath::Sqrt(x*x + y*y); 
    
    ((TH1D*)gDirectory->Get("h100"))->Fill(mass);
    ((TH1D*)gDirectory->Get("h101"))->Fill(m0);
    ((TH1D*)gDirectory->Get("h102"))->Fill(id1);
    ((TH1D*)gDirectory->Get("h102"))->Fill(id2);
    ((TH1D*)gDirectory->Get("h103"))->Fill(radius);

    ((TH1D*)gDirectory->Get("h104"))->Fill(0.);
    if (isRecEl(ip1))  ((TH1D*)gDirectory->Get("h104"))->Fill(1.);
    if (isRecEl(ip2))  ((TH1D*)gDirectory->Get("h104"))->Fill(2.);
    if ((isRecEl(ip1)) && (isRecEl(ip2))) ((TH1D*)gDirectory->Get("h104"))->Fill(3.);

    if (isRecEl(ip1) || isRecEl(ip2)) {
      if (mass < 0.010) {
        ((TH2D*)gDirectory->Get("h200"))->Fill(p1, p2);
        ((TH2D*)gDirectory->Get("h201"))->Fill(x, y);
      }
    }
  }

}

// ----------------------------------------------------------------------
void recoilBase::maskPi0(int mode) {
  Int_t i(0);
  for (i = 0; i < nGam; ++i) { pi0LockGam[i] = 0; }
  
  for (i = 0; i < nPi0; ++i) {
    if ((TMath::Abs(d1Pi0Lund[i]) == 22) && (TMath::Abs(d2Pi0Lund[i]) == 22)) {
      int ga1 =  d1Pi0Index[i]-1;
      int ga2 =  d2Pi0Index[i]-1;
      pi0LockGam[ga1] = 1;
      pi0LockGam[ga2] = 1;
      //        cout << "pi0->gg mass = " << m0Pi0[i] << " from  "  << ga1 << "  " << ga2 << endl;

    }
    else {
      // this happens. need to investigate why!?x
      //        cout << "pi0 decay mass = " << m0Pi0[i] << "   " << d1Pi0Lund[i] << "," << d2Pi0Lund[i] << ") " << endl;
    }

  }
     
}

// ----------------------------------------------------------------------
void recoilBase::maskConversions( int mode ) 
{
  Int_t i(0),ip1,ip2;
  Double_t mass;
  
  for (i = 0; i < 100; ++i) {
    goodConv[i] = 0;
    convLockTrk[i] = -1;
  }
  
  for (i = 0; i < nGConv; ++i) {
    mass = massGConv[i];
    ip1 = d1GConvIndex[i] - 1; 
    ip2 = d2GConvIndex[i] - 1; 
    if ((isRecEl(ip1) || isRecEl(ip2)) && (mass<CONVMAXMEE)) {
      goodConv[i] = 1;
      convLockTrk[ip1] = i;
      convLockTrk[ip2] = i;
    }
  }
}
  
// ----------------------------------------------------------------------
void recoilBase::maskKshorts(int mode, int Ver) {
  static Bool_t first(kTRUE);
  if (first) {
    first = kFALSE;
    fHistFile->cd();
    fHistFile->mkdir("kshorts", "kshorts");
    fHistFile->cd("kshorts");
    TH1D *h;
    TH2D *h2;
    char name[100], title[100];
    sprintf(name, "mcTruth");  sprintf(title, "mass all Kshorts MC Truth");  h = new TH1D(name, title, 4, 0., 4.); 
    sprintf(name, "mcTotal");  sprintf(title, "mass all Kshorts MC Total");  h = new TH1D(name, title, 4, 0., 4.); 
    sprintf(name, "mcBreco");  sprintf(title, "mass all Kshorts MC BRECO");  h = new TH1D(name, title, 4, 0., 4.); 
    sprintf(name, "mcRecoil");  sprintf(title, "mass all Kshorts MC Recoil");  h = new TH1D(name, title, 4, 0., 4.); 
    // KS->pi+pi-
    sprintf(name, "e100");  sprintf(title, "Kshort mass with electrons");  h = new TH1D(name, title, 100, 0.46, 0.54); 
    sprintf(name, "m100");  sprintf(title, "Kshort mass with muons");  h = new TH1D(name, title, 100, 0.46, 0.54); 
    sprintf(name, "k100");  sprintf(title, "Kshort mass with kaons");  h = new TH1D(name, title, 100, 0.46, 0.54); 
    sprintf(name, "h100");  sprintf(title, "mass all Kshorts");  h = new TH1D(name, title, 100, 0.46, 0.54); 
    sprintf(name, "h101");  sprintf(title, "mass Kshorts with one BRECO overlap");  h = new TH1D(name, title, 100, 0.46, 0.54); 
    sprintf(name, "h102");  sprintf(title, "mass Kshorts with two BRECO overlaps");  h = new TH1D(name, title, 100, 0.46, 0.54); 
    sprintf(name, "h103");  sprintf(title, "mass Kshorts with no  BRECO overlaps");  h = new TH1D(name, title, 100, 0.46, 0.54); 
    sprintf(name, "h104");  sprintf(title, "mass of selected Kshorts");  h = new TH1D(name, title, 100, 0.46, 0.54); 
    sprintf(name, "h105");  sprintf(title, "mass of Kshorts (signalBox)");  h = new TH1D(name, title, 100, 0.46, 0.54); 
    sprintf(name, "h106");  sprintf(title, "mass of selected Kshorts (MCT matched)");  h = new TH1D(name, title, 100, 0.46, 0.54); 

    sprintf(name, "nkcharged");  sprintf(title, "KS -> pi+pi-");  h = new TH1D(name, title, 20, 0., 20.); 
    sprintf(name, "nkneutral");  sprintf(title, "KS -> pi0pi0");  h = new TH1D(name, title, 20, 0., 20.); 

    // KS->pi0pi0
    sprintf(name, "mcpi0egamma");  sprintf(title, "gamma energy for pi0");  h = new TH1D(name, title, 50, 0.0, 1.0); 
    sprintf(name, "mcpi0m100");  sprintf(title, "pi0 mass");  h = new TH1D(name, title, 70, 0.100, 0.170); 
    sprintf(name, "mcgaE_pi0M");  sprintf(title, " ");  h2 = new TH2D(name, title, 50, 0.0, 1.0, 70, 0.100, 0.170); 
    sprintf(name, "mcpi0p0");     sprintf(title, "pi0 momentum");  h = new TH1D(name, title, 50, 0.0, 1.0); 
    sprintf(name, "mcksp0");  sprintf(title, "Ks momentum");  h = new TH1D(name, title, 50, 0., 2.5); 

    sprintf(name, "secmom");  sprintf(title, "secmom for pi0");  h = new TH1D(name, title, 50, 0.0, 0.02); 
    sprintf(name, "lmom");  sprintf(title, "LAT for pi0");  h = new TH1D(name, title, 50, 0.0, 1.0); 
    sprintf(name, "ncry");  sprintf(title, "NCRY for pi0");  h = new TH1D(name, title, 20, 0.0, 20.0); 
    sprintf(name, "nbump");  sprintf(title, "NBUMP for pi0");  h = new TH1D(name, title, 20, 0.0, 20.0); 

    sprintf(name, "gaE_pi0M");  sprintf(title, " ");  h2 = new TH2D(name, title, 50, 0.0, 1.0, 70, 0.100, 0.170); 
    sprintf(name, "pi0egamma");  sprintf(title, "gamma energy for pi0");  h = new TH1D(name, title, 50, 0.0, 1.0); 
    sprintf(name, "pi0m100");  sprintf(title, "pi0 mass");  h = new TH1D(name, title, 70, 0.100, 0.170); 
    sprintf(name, "pi0p0");     sprintf(title, "pi0 momentum");  h = new TH1D(name, title, 50, 0.0, 1.0); 
    sprintf(name, "ksp0");    sprintf(title, "Ks momentum");  h = new TH1D(name, title, 50, 0., 2.5); 

    sprintf(name, "i100");  sprintf(title, "Kshort ->pi0pi0 mass");  h = new TH1D(name, title, 100, 0.46, 0.54); 
    sprintf(name, "mci100");  sprintf(title, "Kshort ->pi0pi0 mass");  h = new TH1D(name, title, 100, 0.46, 0.54); 

    sprintf(name, "i101");  sprintf(title, "Kshort ->pi0pi0 mass");  h = new TH1D(name, title, 100, 0.46, 0.54); 
    sprintf(name, "mci101");  sprintf(title, "Kshort ->pi0pi0 mass");  h = new TH1D(name, title, 100, 0.46, 0.54); 

    sprintf(name, "i102");  sprintf(title, "Kshort ->pi0pi0 mass");  h = new TH1D(name, title, 100, 0.46, 0.54); 
    sprintf(name, "mci102");  sprintf(title, "Kshort ->pi0pi0 mass");  h = new TH1D(name, title, 100, 0.46, 0.54); 

    sprintf(name, "i103");  sprintf(title, "Kshort ->pi0pi0 mass");  h = new TH1D(name, title, 100, 0.46, 0.54); 
    sprintf(name, "mci103");  sprintf(title, "Kshort ->pi0pi0 mass");  h = new TH1D(name, title, 100, 0.46, 0.54); 
  }
  fHistFile->cd("kshorts");

  Int_t i(0), j(0), goodKcharged(0), goodKneutral(0);
  Int_t partOfKs[100];
  Int_t Brectrktmp, bmctmp; 
  brecoOverlap=2;

  // -- MCT kshort counter
  if (fIsMC) {
    int Brectrktmp(0); 
    for (Int_t imc = 0; imc < nMc; ++imc) {
      if (TMath::Abs(idMc[imc]) == 310) {
	bmctmp = MCB0[indexbestB];  if(fBrecoCharge != 0) bmctmp = MCChB[indexbestB];
	if (bmctmp-1 >=0) {
	  if (isAncestor(bmctmp-1, imc)) {
	    ((TH1D*)gDirectory->Get("mcTruth"))->Fill(1.);
	  } else {
	    ((TH1D*)gDirectory->Get("mcTruth"))->Fill(2.);
	  }
	} else {
	    ((TH1D*)gDirectory->Get("mcTruth"))->Fill(0.);
	}
	int cntTot(0), cntBreco(0), cntRecoil(0);
	for (int itrk = 0 ; itrk < nTrk; ++itrk) {
	  if (mothMc[IndexTrk[itrk]-1]-1 == imc) {
	    ++cntTot;
	    Brectrktmp = B0RecTrk[itrk]; if(fBrecoCharge != 0) Brectrktmp = chBRecTrk[itrk];
	    if (Brectrktmp  == brecoOverlap)  ++cntBreco;
	    if (Brectrktmp  != brecoOverlap)  ++cntRecoil;
	  }
	  ((TH1D*)gDirectory->Get("mcTotal"))->Fill(cntTot);
	  ((TH1D*)gDirectory->Get("mcBreco"))->Fill(cntBreco);
	  ((TH1D*)gDirectory->Get("mcRecoil"))->Fill(cntRecoil);
	}
      }
    }
  }
  for (i = 0; i < 100; ++i) { 
    kshortLockTrk[i] = 0; 
    kshortLockGam[i] = 0; 
    partOfKs[i] = 0; 
    goodKshort[i] = 1;
    goodWe[i] = 0;
    goodWk[i] = 0;
    goodNr[i] = 0;
  }
  if (Ver) cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  for (i = 0; i < nKs; ++i) {
    if (Ver) cout << "Ks[" << i << "] -> (" << d1KsIndex[i]-1 << "," << d2KsIndex[i]-1 << ") " 
		  << "with daughters " << d1KsLund[i] << " " << d2KsLund[i]
		  << " and mass = " << massKs[i] << endl;
    partOfKs[d1KsIndex[i]-1] += 1; 
    partOfKs[d2KsIndex[i]-1] += 1; 
  }
  
  // =============
  // -- KS->pi+pi-
  // =============
  // -- Two passes over KS block: (1) Reject trivial cases: Ks with kaons and electron, overlaps with BRECO
  for (i = 0; i < nKs; ++i) {
    if (TMath::Abs(d1KsLund[i]) != 211) continue;
    if (TMath::Abs(d2KsLund[i]) != 211) continue;
    int pi1 = d1KsIndex[i]-1;
    int pi2 = d2KsIndex[i]-1; 
    double mass = massKs[i];
    double residual = TMath::Abs(KAZMASS - mass); 
    ((TH1D*)gDirectory->Get("h100"))->Fill(mass);
    if (isRecMu(pi1) || isRecMu(pi2)) ((TH1D*)gDirectory->Get("m100"))->Fill(mass);	
    // -- Veto KS candidates with electrons or kaons (note: KS candidates with muons peak at m(KS))
    if (isRecEl(pi1)) {
      if (Ver) cout << "Rejecting Ks[" << i <<"] due to identified electron " << pi1 << endl;
      goodKshort[i] = 0;
      goodWe[i] += 1;
      ((TH1D*)gDirectory->Get("e100"))->Fill(mass);	
      continue;
    }
    if (isRecEl(pi2)) {
      if (Ver) cout << "Rejecting Ks[" << i <<"] due to identified electron " << pi2 << endl;
      goodKshort[i] = 0;
      goodWe[i] += 1;
      ((TH1D*)gDirectory->Get("e100"))->Fill(mass);	
      continue;
    }
    if (isRecKaon(pi1)) {
      if (Ver) cout << "Rejecting Ks[" << i <<"] due to identified kaon " << pi1 << endl;
      goodKshort[i] = 0;
      goodWk[i] += 1;
      ((TH1D*)gDirectory->Get("k100"))->Fill(mass);	
      continue;
    }
    if (isRecKaon(pi2)) {
      if (Ver) cout << "Rejecting Ks[" << i <<"] due to identified kaon " << pi2 << endl;
      goodKshort[i] = 0;
      goodWk[i] += 1;
      ((TH1D*)gDirectory->Get("k100"))->Fill(mass);	
      continue;
    }
    // -- Check against overlap with BRECO candidate
    int overlap(0);
    Brectrktmp = B0RecTrk[pi1];
    if(fBrecoCharge != 0) Brectrktmp = chBRecTrk[pi1];
    if (Brectrktmp  == brecoOverlap) {
      if(Ver) cout << "Rejecting Ks[" << i << "] due to pi1 at " << pi1 << " overlaps with BRECO " << endl;
      overlap++;
    }
    Brectrktmp = B0RecTrk[pi2];
    if(fBrecoCharge != 0) Brectrktmp = chBRecTrk[pi2];
    if (Brectrktmp  == brecoOverlap) {
      if(Ver) cout << "Rejecting Ks[" << i << "] due to pi2 at " << pi2 << " overlaps with BRECO " << endl;
      overlap++;
    }
    if (overlap > 0) {
      goodKshort[i] = 0;
      goodNr[i] = overlap;
      if (overlap == 1) ((TH1D*)gDirectory->Get("h101"))->Fill(mass);
      if (overlap == 2) ((TH1D*)gDirectory->Get("h102"))->Fill(mass);
      continue;
    } else {
      ((TH1D*)gDirectory->Get("h103"))->Fill(mass);
    }
  } // first pass over KS block

  // -- Second Pass: Reject Ks whose daughters are part of a 'better' KS
  for (i = 0; i < nKs; ++i) {
    if (goodKshort[i] == 0) continue;
    if (TMath::Abs(d1KsLund[i]) != 211) continue;
    if (TMath::Abs(d2KsLund[i]) != 211) continue;
    int takeThis(1); 
    int pi1 = d1KsIndex[i]-1;
    int pi2 = d2KsIndex[i]-1; 
    double mass = massKs[i];
    double residual = TMath::Abs(KAZMASS - mass); 
    if (partOfKs[pi1] > 1) {
      for (j = 0; j < nKs; ++j) {
	if (i == j) continue;
	if (goodKshort[j] == 0) continue;
	// -- Compare residuals
	if ((d1KsIndex[j]-1 == pi1) || (d2KsIndex[j]-1 == pi1)) {
	  if (TMath::Abs(KAZMASS - massKs[j]) < residual) {
	    if (Ver) cout << "Rejecting KS[" << i << "] -> (" << pi1 << "," << pi2 << ")"
			  << " with mass = " << mass 
			  << ", since KS[" << j << "] -> (" << d1KsIndex[j]-1 << "," << d2KsIndex[j]-1 << ")"
			  << " has better mass: " << massKs[j] << endl;
	    takeThis = 0; 
	    goodKshort[i] = 0;
	    break;
	  }
	}
      }
    }
    if (partOfKs[pi2] > 1) {
      for (j = 0; j < nKs; ++j) {
	if (i == j) continue;
	// -- do not compare against KS overlapping with BRECO
	if (goodKshort[j] == 0) continue;
	// -- Compare residuals
	if ((d1KsIndex[j]-1 == pi2) || (d2KsIndex[j]-1 == pi2)) {
	  if (TMath::Abs(KAZMASS - massKs[j]) < residual) {
	    if (Ver) cout << "Rejecting KS[" << i << "] -> (" << pi1 << "," << pi2 << ")"
			  << " mass = " << mass 
			  << ", since KS[" << j << "] -> (" << d1KsIndex[j]-1 << "," << d2KsIndex[j]-1 << ")"
			  << " has better mass: " << massKs[j] << endl;
	    takeThis = 0; 
	    goodKshort[i] = 0;
	    break;
	  }
	}
      }
    }
    if (takeThis == 1) {
      ((TH1D*)gDirectory->Get("h104"))->Fill(mass);
      if (signalBox) ((TH1D*)gDirectory->Get("h105"))->Fill(mass);
      goodKshort[i] = 1;
      kshortLockTrk[pi1] = 1;
      kshortLockTrk[pi2] = 1;      
      if (fIsMC && (TMath::Abs(idMc[MCKs[i]-1])==310)) ((TH1D*)gDirectory->Get("h106"))->Fill(mass);
      if(Ver) cout << "Locking KS[" << i << "] -> ("  << pi1 << "," << pi2 << ") mass = " << mass << endl;
      ++goodKcharged;
    }
  } // second pass over KS block


  //  return;

  // =============
  // -- KS->pi0pi0
  // =============
  for (i = 0; i < nKs; ++i) {
    if (TMath::Abs(d1KsLund[i]) != 111) continue;
    if (TMath::Abs(d2KsLund[i]) != 111) continue;
    int pi1 = d1KsIndex[i]-1;
    int pi2 = d2KsIndex[i]-1; 
    double pi0m1 = m0Pi0[pi1];
    double pi0m2 = m0Pi0[pi2];  
    ((TH1D*)gDirectory->Get("pi0m100"))->Fill(pi0m1);
    ((TH1D*)gDirectory->Get("pi0m100"))->Fill(pi0m2);
    ((TH1D*)gDirectory->Get("pi0p0"))->Fill(pPi0[pi1]);
    ((TH1D*)gDirectory->Get("pi0p0"))->Fill(pPi0[pi2]);

    int g[4] = {d1Pi0Index[pi1]-1, d2Pi0Index[pi1]-1, d1Pi0Index[pi2]-1, g[3] = d2Pi0Index[pi2]-1};
    int goodGammas(1);
    for (int ig = 0; ig < 4; ++ig) {
      ((TH1D*)gDirectory->Get("pi0egamma"))->Fill(energyGam[g[ig]]);
      ((TH1D*)gDirectory->Get("secmom"))->Fill(secMomGam[g[ig]]);
      ((TH1D*)gDirectory->Get("lmom"))->Fill(lMomGam[g[ig]]);
      ((TH1D*)gDirectory->Get("ncry"))->Fill(nCryGam[g[ig]]);
      ((TH1D*)gDirectory->Get("nbump"))->Fill(nBumpGam[g[ig]]);
      if (energyGam[g[ig]] < 0.010) goodGammas = 0;
      if (nCryGam[g[ig]] < 3) goodGammas = 0;
      if (nBumpGam[g[ig]] < 1) goodGammas = 0;
    }    

    if (goodGammas == 0) continue;

    ((TH2D*)gDirectory->Get("gaE_pi0M"))->Fill(energyGam[g[0]], pi0m1);
    ((TH2D*)gDirectory->Get("gaE_pi0M"))->Fill(energyGam[g[1]], pi0m1);
    ((TH2D*)gDirectory->Get("gaE_pi0M"))->Fill(energyGam[g[2]], pi0m2);
    ((TH2D*)gDirectory->Get("gaE_pi0M"))->Fill(energyGam[g[3]], pi0m2);

    // -- Check for MC truth matched pi0
    int mcpi1 = MCPi0[pi1]-1;
    int mcpi2 = MCPi0[pi2]-1;
    if (mcpi1 > nMc) {
      cout << " corrupt index" << endl;
      continue;
    }
    if (mcpi2 > nMc) {
      cout << " corrupt index" << endl;
      continue;
    }

    if (mcpi1 > -1 && idMc[mcpi1] == 111) {
      ((TH1D*)gDirectory->Get("mcpi0egamma"))->Fill(energyGam[g[0]]);
      ((TH1D*)gDirectory->Get("mcpi0egamma"))->Fill(energyGam[g[1]]);
      ((TH1D*)gDirectory->Get("mcpi0m100"))->Fill(pi0m1);

      ((TH2D*)gDirectory->Get("mcgaE_pi0M"))->Fill(energyGam[g[0]], pi0m1);
      ((TH2D*)gDirectory->Get("mcgaE_pi0M"))->Fill(energyGam[g[1]], pi0m1);
    }

    if (mcpi2 > -1 && idMc[mcpi2] == 111) {
      ((TH1D*)gDirectory->Get("mcpi0egamma"))->Fill(energyGam[g[2]]);
      ((TH1D*)gDirectory->Get("mcpi0egamma"))->Fill(energyGam[g[3]]);
      ((TH1D*)gDirectory->Get("mcpi0m100"))->Fill(pi0m2);

      ((TH2D*)gDirectory->Get("mcgaE_pi0M"))->Fill(energyGam[g[2]], pi0m2);
      ((TH2D*)gDirectory->Get("mcgaE_pi0M"))->Fill(energyGam[g[3]], pi0m2);
    }

    double mass = massKs[i];
    double residual = TMath::Abs(KAZMASS - mass); 
    ((TH1D*)gDirectory->Get("i100"))->Fill(mass);
    ((TH1D*)gDirectory->Get("ksp0"))->Fill(pKs[i]);
    int mcTruthMatched(0);
    if (mcpi1 > -1 && idMc[mcpi1] == 111 && mcpi2 > -1 && idMc[mcpi2] == 111) mcTruthMatched = 1;
    if (mcTruthMatched) { 
      ((TH1D*)gDirectory->Get("mci100"))->Fill(mass);
      ((TH1D*)gDirectory->Get("mcpi0p0"))->Fill(pPi0[pi1]);
      ((TH1D*)gDirectory->Get("mcpi0p0"))->Fill(pPi0[pi2]);
      ((TH1D*)gDirectory->Get("mcksp0"))->Fill(pKs[i]);
    }

    if (energyGam[g[0]] > 0.1 && energyGam[g[1]] > 0.1 && energyGam[g[2]] > 0.1 && energyGam[g[3]] > 0.1) {
      ((TH1D*)gDirectory->Get("i101"))->Fill(mass);
      if (mcTruthMatched == 1) {
	((TH1D*)gDirectory->Get("mci101"))->Fill(mass);
      }

      if (pPi0[pi1] > 0.1 && pPi0[pi2] > 0.1) {
	((TH1D*)gDirectory->Get("i102"))->Fill(mass);
	if (mcTruthMatched == 1) {
	  ((TH1D*)gDirectory->Get("mci102"))->Fill(mass);
	}

	if (pi0m1 > 0.124 && pi0m1 < 0.144 && pi0m2 > 0.124 && pi0m2 < 0.144) {
	  ((TH1D*)gDirectory->Get("i103"))->Fill(mass);
	  if (mcTruthMatched == 1) {
	    ((TH1D*)gDirectory->Get("mci103"))->Fill(mass);
	  }
	}
      }
    }
    
    ++goodKneutral;
  }

  ((TH1D*)gDirectory->Get("nkcharged"))->Fill(goodKcharged);
  ((TH1D*)gDirectory->Get("nkneutral"))->Fill(goodKneutral);
  

}



// ----------------------------------------------------------------------

void recoilBase::mcGam() {

  for(int je=0;je<80;je++) {
    ifromBGam[je] = 0;
  }
  
  int idB[2];

  int s = 0;

  idB[0] = -100;
  idB[1] = -100;
  
  for (int f = 0; f < nMc; f++) {
    if ((TMath::Abs(idMc[f]) == 511) || (TMath::Abs(idMc[f]) == 521) ) {
      idB[s] = f;
      s++;
    }
  }

  for (int h = 0; h < nGam; h++) {    
    if(idGam[h] == 22) {
      if (isAncestor(idB[0] ,h) ||  isAncestor(idB[1],h)) { 
	ifromBGam[h] = 1; 
      }
    }
  }
}


void recoilBase::timestamp( int lun ) {

  //  This function prints out eventIDs and run numbers in the official
  //  BaBar format (e.g. eventID in hex as aa:bb:xxxxxx/yyyyyyyy:Z, 
  //  run number as an 8-digit decimal number following "run=" token
  //  Output is written to LUN <lun>

  int k;
  int temp;
  int iplat;
  Int_t ipart;
  char cplat[2];
  char cpart[8],tspart[8];
  int nplat,npart;
  
  const char * array[17]={"G","H","J","K","L","M","N","P","Q","R","S",
			  "T","U","V","W","X","Y"};

  if(lun == 0) return;

  iplat=platform;
  ipart = partition;

  //      iplat=127
  //      partition=8388607             ! 0x7fffff

  // implement wildcard behavior of platform
  if ( iplat == 0 ) {
    sprintf(cplat,"%s","\0");
    nplat=1;
  }  else {
    sprintf(cplat,"%X%s",iplat,"\0");
    nplat=2;
  }

  // implement wildcard behavior of partitionMask

  if ( partition == 0 ) {
    sprintf(cpart,"%s","*\0");
    sprintf(tspart,"%s","0\0");
    npart=1;
    //   Need to redo twos complement interpretation (MC has negative partition id)
    //  i.e. map -2 ->FFFFFFFE , -1 -> FFFFFFFF, etc.
  } else if (partition < 0) {
    //    temp = (Int_t) (pow(2,32) + partition);
    sprintf(cpart,"%X%s", temp,"\0");
    sprintf(tspart,"%X%s", temp,"\0");
    npart=8;
  }  else {
    sprintf(cpart,"%X", partition);
    sprintf(tspart,"%X", partition);
    npart=6;
  }

  // checksum
  k= (TMath::Abs(iplat+ipart+upperID+lowerID)%17) +1;

  // now print
  //  if(fIntPurity < 30) cout<<" Ts file must be set up "<<fIntPurity<<endl;
  if(fIntPurity > 80) lun = lun+2;
  if(fIntPurity > 50) lun = lun+2;
  if(fIntPurity > 10) lun = lun+2;

  char name[100];
  char name2[100];
  sprintf(name,"outp3/ts_file.%i",lun);
  sprintf(name2,"outp3/ts_file.%i",lun+1);
  //  ofstream outFile(name,ios::app|ios::uppercase);
  //  ofstream outFile2(name2,ios::app|ios::uppercase);
  ofstream outFile(name);
  ofstream outFile2(name2);
  outFile.setf(ios::uppercase);
  outFile2.setf(ios::uppercase);
  outFile<<hex<<setw(nplat)<<setfill('0')<<cplat<<":"<<setw(npart)
	 <<setfill('0')<<cpart<<":"<<setw(8)<<setfill('0')<<upperID
	 <<"/"<<setw(8)<<setfill('0')<<lowerID<<":"<<array[k-1]<<"  run="<<dec<<setw(8)
	 <<setfill('0')<<runNumber<<endl;
  outFile2<<hex<<setw(nplat)<<setfill('0')<<cplat<<":"<<setw(npart)
	 <<setfill('0')<<tspart<<":"<<setw(8)<<setfill('0')<<upperID
	 <<"/"<<setw(8)<<setfill('0')<<lowerID<<":"<<array[k-1]<<"  run="<<dec<<setw(8)
	 <<setfill('0')<<runNumber<<"  "<<dec<<fPurity<<" "
	  <<fBmode<<" "<<endl;
}

// ----------------------------------------------------------------------
Int_t recoilBase::bestKsIndex(Int_t isVerbose)
{
  Double_t bestdm(1000);
  Int_t ibest(-1);
  
  if (isVerbose)
    cout << "----------------" << endl;
  for (Int_t i=0; i<nKs ; i++) {
    if (goodKshort[i] == 1) {
      if (isVerbose>0) 
	cout << "Good KS Index: " << i << "  Mass : " << massKs[i] << endl;
      Double_t currdm = TMath::Abs(KAZMASS - massKs[i]);
      if (currdm<bestdm) {
        bestdm=currdm;
        ibest=i;
      }
    }
  }
  if (isVerbose>0) {
    cout << "N. KS = " << nKs << endl
	 << "Good KS: = " << fNKshort << endl
	 << "KS index = " << ibest << endl;
  }
  return ibest;
}

  
// ----------------------------------------------------------------------

void recoilBase::GamStudy( int mode, TLorentzVector &m4Xhad) {

  //  cout<<m4Xhad.Mag()<<endl;
  //  cout<<"New Gamma Study ************** New event"<<endl;
  Int_t Brecgamtmp ;
  Int_t brecoOv(2); // this is flag for overlap with BRECO candidate  
  Int_t bmctmp;
  fHistFile->cd();
  fHistFile->cd("Photon");
  int i;
  for(int je=0;je<15;je++) {
    fcountNeu[je] = 0;
    pGS[je] = m4Xhad;
  }

  int inp1, inp2;
  bool contKs =kFALSE;

  if (fIsMC){
    for (i = 0; i < nGam; ++i) {
      Brecgamtmp = B0RecGam[i];
      if(mode) Brecgamtmp = chBRecGam[i];
      if(Brecgamtmp == brecoOv) continue;
      if(ecalGam[i] < 0) continue;
      
      ((TH1D*)gDirectory->Get("EnGam"))->Fill(energyGam[i], 1.);
      ((TH1D*)gDirectory->Get("lMomGam"))->Fill( lMomGam[i], 1.);
      ((TH1D*)gDirectory->Get("ZMom42Gam"))->Fill( ZMom42Gam[i], 1.);
      ((TH2D*)gDirectory->Get("E vs th"))->Fill(TMath::Cos(thetaGam[i]),energyGam[i]);
      ((TH2D*)gDirectory->Get("E vs phi"))->Fill(energyGam[i],phiGam[i]);
      
      // return kFALSE; 
      //      cout<<"tr with initial"<<endl;
      int MomId = 1;

      if(idGam[i] == 22) {
	((TH1D*)gDirectory->Get("EnGamMP"))->Fill(energyGam[i], 1.);  

	if( isAncestor(fB1Index,i) || isAncestor(fB2Index,i) ) {
	  ((TH1D*)gDirectory->Get("EGam Mt"))->Fill(energyGam[i], 1.);
	  if(fVub == 1) ((TH1D*)gDirectory->Get("EGam Mt vub"))->Fill(energyGam[i], 1.);
	  if(fVcb == 1) ((TH1D*)gDirectory->Get("EGam Mt vcb"))->Fill(energyGam[i], 1.);
	} else {
	  ((TH1D*)gDirectory->Get("EGam Un"))->Fill(energyGam[i], 1.);
	  if(fVub == 1) ((TH1D*)gDirectory->Get("EGam Un vub"))->Fill(energyGam[i], 1.);
	  if(fVcb == 1) ((TH1D*)gDirectory->Get("EGam Un vcb"))->Fill(energyGam[i], 1.);
	}
	
	for (int ik = 0; ik < nMc; ik++) {
	  if(TMath::Abs(idMc[ik]) == 310) {
	    bool isfirPI = kTRUE;
	    int Npi = 0;
	    for (int ip = 0; ip < nMc; ip++) {
	      if((TMath::Abs(idMc[ip]) == 111) && (mothMc[ip]-1 == ik) && Npi <2) {
		if(isfirPI) {
		  inp1 = ip;
		  isfirPI = kFALSE;
		} else {
		  inp2 = ip;
		}
		Npi++;
	      } else if ((TMath::Abs(idMc[ip]) == 211) && (mothMc[ip]-1 == ik) && Npi <2) {
		if(isfirPI) {
		  inp1 = ip;
		  isfirPI = kFALSE;
		} else {
		  inp2 = ip;
		}
		Npi++;
	      }
	    }

	    //Event with Ks -> pi0pi0 ; Ks -> pi+pi- 

	    if(inp2 > 0 && inp1 > 0 && Npi == 2) contKs = kTRUE;

	  }
	}

	if(contKs) {
	  ((TH1D*)gDirectory->Get("EnGamKs"))->Fill(energyGam[i], 1.);
	} else {
	  ((TH1D*)gDirectory->Get("EnGamNoKs"))->Fill(energyGam[i], 1.);
	}	

	int indexmom;
	bool isFir = kTRUE;
	while(MomId != 0) { 
	  if(isFir) {
	    indexmom = mothMc[i]-1;
	    isFir = kFALSE;
	  }
	  if ((indexmom > 0) && (indexmom < nMc)) {
	    MomId = idMc[indexmom]; //Id of the mother
	  } else {
	    if(indexmom == 0) findUps++;
	    MomId = 0;
	    findPro++;
	  }
	  if( TMath::Abs(MomId) >= 1000000000 ) {
	    ((TH1D*)gDirectory->Get("ThetaGam Un"))->Fill(TMath::Cos(thetaGam[i]), 1.);
	    ((TH1D*)gDirectory->Get("PhiGam Un"))->Fill(phiGam[i], 1.);
	    MomId = 0;
	  }
	  if( (TMath::Abs(MomId) == 511) || (TMath::Abs(MomId) == 521) ) { //B daughter
	    ((TH1D*)gDirectory->Get("EnGamMB vub"))->Fill(energyGam[i], 1.);
	    if(energyGam[i] <= 0.080 ) {
	      ((TH1D*)gDirectory->Get("EnGam (MB) cut"))->Fill(energyGam[i], 1.);
	      ((TH1D*)gDirectory->Get("lMomGam (MB) cut"))->Fill( lMomGam[i], 1.);
	      ((TH1D*)gDirectory->Get("ZMom42Gam (MB) cut"))->Fill( ZMom42Gam[i], 1.);
	      ((TH1D*)gDirectory->Get("PhiGam (MB) cut"))->Fill(phiGam[i], 1.);
	    } else {
	      ((TH1D*)gDirectory->Get("EnGam (MB)"))->Fill(energyGam[i], 1.);
	      ((TH1D*)gDirectory->Get("lMomGam (MB)"))->Fill( lMomGam[i], 1.);
	      ((TH1D*)gDirectory->Get("ZMom42Gam (MB)"))->Fill( ZMom42Gam[i], 1.);
	      ((TH1D*)gDirectory->Get("PhiGam (MB)"))->Fill(phiGam[i], 1.);
	    }
	    ((TH2D*)gDirectory->Get("E vs theta (MB)"))->Fill(TMath::Cos(thetaGam[i]),energyGam[i]);
	    MomId = 0;
	  } else if ( (MomId == 0) ||(MomId == 70553))  { //I.M. daughter
	    if(MomId == 70553) cout<<"Upsi ev."<<endl;
	    ((TH1D*)gDirectory->Get("EnGamMA"))->Fill(energyGam[i], 1.);
	    if(energyGam[i] <= 0.080 ) {
	      ((TH1D*)gDirectory->Get("EnGam (MA) cut"))->Fill(energyGam[i], 1.);
	      ((TH1D*)gDirectory->Get("lMomGam (MA) cut"))->Fill( lMomGam[i], 1.);
	      ((TH1D*)gDirectory->Get("ZMom42Gam (MA) cut"))->Fill( ZMom42Gam[i], 1.);
	      ((TH1D*)gDirectory->Get("PhiGam (MA) cut"))->Fill(phiGam[i], 1.);
	    } else {
	      ((TH1D*)gDirectory->Get("EnGam (MA)"))->Fill(energyGam[i], 1.);
	      ((TH1D*)gDirectory->Get("lMomGam (MA)"))->Fill( lMomGam[i], 1.);
	      ((TH1D*)gDirectory->Get("ZMom42Gam (MA)"))->Fill( ZMom42Gam[i], 1.);
	      ((TH1D*)gDirectory->Get("PhiGam (MA)"))->Fill(phiGam[i], 1.);
	    }
	    ((TH2D*)gDirectory->Get("E vs theta (MA)"))->Fill(TMath::Cos(thetaGam[i]),energyGam[i]);
	    MomId = 0;
	  } else {
	    indexmom = mothMc[indexmom]-1;
	  }
	} //end while
      } else {
	((TH1D*)gDirectory->Get("EnGamNMP"))->Fill(energyGam[i], 1.);
	if(energyGam[i] <= 0.080 ) {
	  ((TH1D*)gDirectory->Get("EnGam (NMP) cut"))->Fill(energyGam[i], 1.);
	  ((TH1D*)gDirectory->Get("lMomGam (NMP) cut"))->Fill( lMomGam[i], 1.);
	  ((TH1D*)gDirectory->Get("ZMom42Gam (NMP) cut"))->Fill( ZMom42Gam[i], 1.);
	  ((TH1D*)gDirectory->Get("PhiGam (NMP) cut"))->Fill(phiGam[i], 1.);
	} else {
	  ((TH1D*)gDirectory->Get("EnGam (NMP)"))->Fill(energyGam[i], 1.);
	  ((TH1D*)gDirectory->Get("lMomGam (NMP)"))->Fill( lMomGam[i], 1.);
	  ((TH1D*)gDirectory->Get("ZMom42Gam (NMP)"))->Fill( ZMom42Gam[i], 1.);
	  ((TH1D*)gDirectory->Get("PhiGam (NMP)"))->Fill(phiGam[i], 1.);
	}
	((TH2D*)gDirectory->Get("E vs theta (NMP)"))->Fill(TMath::Cos(thetaGam[i]),energyGam[i]);
      }
      
      //#Crystal distribution
      int ind1pi0GS(0), ind2pi0GS(0);
      for(int jj=0;jj<nPi0;jj++) {
	ind1pi0GS =   d1Pi0Index[jj] - 1;
	ind2pi0GS =   d2Pi0Index[jj] - 1;
	if(i == ind1pi0GS || i == ind2pi0GS ) {
	  if(energyGam[i] <= 0.080 ) {
	    ((TH1D*)gDirectory->Get("CrySoftpi"))->Fill(nCryGam[i], 1.);
	  } else {
	    ((TH1D*)gDirectory->Get("CryHardpi"))->Fill(nCryGam[i], 1.);
	  }
	}
      }
      if(energyGam[i] <= 0.080 ) {
	((TH1D*)gDirectory->Get("CrySoft"))->Fill(nCryGam[i], 1.);
      } else {
	((TH1D*)gDirectory->Get("CryHard"))->Fill(nCryGam[i], 1.);
      }
      
      //Global distributions
      
      if(energyGam[i] <= 0.080 ) {
	((TH2D*)gDirectory->Get("E vs phi cut"))->Fill(energyGam[i],phiGam[i]);
	((TH1D*)gDirectory->Get("ThetaGam"))->Fill(TMath::Cos(thetaGam[i]), 1.);
	((TH1D*)gDirectory->Get("PhiGam"))->Fill(phiGam[i], 1.);
      } else {
	((TH1D*)gDirectory->Get("ThetaGamcut"))->Fill(TMath::Cos(thetaGam[i]), 1.);
      }
      
      if(energyGam[i] <= 0.010) continue;
      
    }
    
    //-----------------------------------------------------------
    TLorentzVector p4tph(0., 0., 0., 0.);

    mk4Vector(p4tph, energyGam[i], thetaGam[i], phiGam[i], 0.);

    //    cout<< m4Xhad.Mag() <<" m4 in gamstudy"<< tmp4Xhad.Mag() <<endl;

    //No cut vector
    //Starting the cut study
    for(int jk=0;jk<10;jk++) {
      //      if(energyGam[i] <= 0.008) cout<<"En of cand = "<<energyGam[i]<<endl;
      if(energyGam[i] >= (jk * 0.010 + 0.020))  {
	//	cout<<i<<" th Ph "<<energyGam[i]<<" <- En ; the -> "<<thetaGam[i]<<endl;
	pGS[jk] += p4tph;
        fcountNeu[jk]++;          
      }
    }
    int tmpind1pi0(0), tmpind2pi0(0);
    for(int jl=0;jl<nPi0;jl++) {
      tmpind1pi0 =   d1Pi0Index[jl] - 1;
      tmpind2pi0 =   d2Pi0Index[jl] - 1;
      if(i == tmpind1pi0 || i == tmpind2pi0 ) {
	pGS[10]    += p4tph;
        fcountNeu[10]++;          
      }
    }
    if(energyGam[i] >= 0.08 && lMomGam[i] >= 0.01) {
      pGS[11]    += p4tph;
      fcountNeu[11]++;          
    }

    if(energyGam[i] >= 0.080 && fcountNeu[12]<2) { 
      pGS[12]    += p4tph;
      fcountNeu[12]++;          
    }

    if( (energyGam[i] >= 0.080) || ((TMath::Cos(thetaGam[i]) <= 0.8)  && (TMath::Cos(thetaGam[i]) >= -0.7)) ) {
      pGS[13]  += p4tph;
      fcountNeu[13]++;          
    }
    if((energyGam[i] >= 0.050 && TMath::Cos(thetaGam[i]) > 0.8 ) ||
       (energyGam[i] >= 0.030 && TMath::Cos(thetaGam[i]) < -0.7) || 
       ((TMath::Cos(thetaGam[i])>=-0.7) && (TMath::Cos(thetaGam[i])<= 0.8))) {
      pGS[14]    += p4tph;
      fcountNeu[14]++;          
    }
  }

  for(int jf=0; jf<15; jf++) {
    tmpfMxhad[jf] = pGS[jf].Mag();
    tmpfTxhad[jf] = pGS[jf].Theta(); 
    tmpfFxhad[jf] = pGS[jf].Phi(); 
    tmpfExhad[jf] = pGS[jf].E();
  }

}

// ----------------------------------------------------------------------
// void recoilBase::recoil(int chbcand, int isVer, Bool_t fisDuplicate) 
// VALERY: REMOVED, should be in the subclass!
// ----------------------------------------------------------------------

// ----------------------------------------------------------------------
void  recoilBase::mxCategory() {

//    char name[100], title[100];
//    TH1D *h;
//    TH2D *h2;
 
//    if (ini == 1) {
//      fHistFile->cd();
//      fHistFile->mkdir(dir, dir);
//      fHistFile->cd(dir);

//      sprintf(name, "u100");  sprintf(title, "categories, Vub enh. ");  h = new TH1D(name, title, 20, 0., 20.); 
//      sprintf(name, "c100");  sprintf(title, "categories, Vub depl.");  h = new TH1D(name, title, 20, 0., 20.); 

//      sprintf(name, "u1000");  sprintf(title, "Mx - MxGen, Vub enh. ");  h2 = new TH2D(name, title, 20, 0., 20., 100, -2., 2.); 
//      sprintf(name, "c1000");  sprintf(title, "Mx - MxGen, Vub depl.");  h2 = new TH2D(name, title, 20, 0., 20., 100, -2., 2.); 

//      sprintf(name, "u1001");  sprintf(title, "Mxfit - MxGen, Vub enh. ");  h2 = new TH2D(name, title, 20, 0., 20., 100, -2., 2.); 
//      sprintf(name, "c1001");  sprintf(title, "Mxfit - MxGen, Vub depl.");  h2 = new TH2D(name, title, 20, 0., 20., 100, -2., 2.); 

//      sprintf(name, "u101");  sprintf(title, "KL momentum, Vub enh. ");  h = new TH1D(name, title, 50, 0., 5.); 
//      sprintf(name, "c101");  sprintf(title, "KL momentum, Vub depl. ");  h = new TH1D(name, title, 50, 0., 5.); 

//      sprintf(name, "u102");  sprintf(title, "n momentum, Vub enh. ");  h = new TH1D(name, title, 50, 0., 5.); 
//      sprintf(name, "c102");  sprintf(title, "n momentum, Vub depl.");  h = new TH1D(name, title, 50, 0., 5.); 

//      sprintf(name, "u103");  sprintf(title, "lost reco momentum, Vub enh. ");  h = new TH1D(name, title, 50, 0., 5.); 
//      sprintf(name, "c103");  sprintf(title, "lost reco momentum, Vub depl.");  h = new TH1D(name, title, 50, 0., 5.); 

//      sprintf(name, "u104");  sprintf(title, "casc momentum, Vub enh. ");  h = new TH1D(name, title, 50, 0., 5.); 
//      sprintf(name, "c104");  sprintf(title, "casc momentum, Vub depl.");  h = new TH1D(name, title, 50, 0., 5.); 

//      sprintf(name, "u105");  sprintf(title, "K+ momentum, Vub enh. ");  h = new TH1D(name, title, 50, 0., 5.); 
//      sprintf(name, "c105");  sprintf(title, "K+ momentum, Vub depl.");  h = new TH1D(name, title, 50, 0., 5.); 



//      sprintf(name, "u201");  sprintf(title, "KL mass, Vub enh. ");  h = new TH1D(name, title, 50, 0., 5.); 
//      sprintf(name, "c201");  sprintf(title, "KL mass, Vub depl. ");  h = new TH1D(name, title, 50, 0., 5.); 

//      sprintf(name, "u202");  sprintf(title, "n mass, Vub enh. ");  h = new TH1D(name, title, 50, 0., 5.); 
//      sprintf(name, "c202");  sprintf(title, "n mass, Vub depl.");  h = new TH1D(name, title, 50, 0., 5.); 

//      sprintf(name, "u203");  sprintf(title, "lost reco mass, Vub enh. ");  h = new TH1D(name, title, 50, 0., 5.); 
//      sprintf(name, "c203");  sprintf(title, "lost reco mass, Vub depl.");  h = new TH1D(name, title, 50, 0., 5.); 

//      sprintf(name, "u204");  sprintf(title, "casc mass, Vub enh. ");  h = new TH1D(name, title, 50, 0., 5.); 
//      sprintf(name, "c204");  sprintf(title, "casc mass, Vub depl.");  h = new TH1D(name, title, 50, 0., 5.); 

//      sprintf(name, "u205");  sprintf(title, "K+ mass, Vub enh. ");  h = new TH1D(name, title, 50, 0., 5.); 
//      sprintf(name, "c205");  sprintf(title, "K+ mass, Vub depl.");  h = new TH1D(name, title, 50, 0., 5.); 

//    }

  int overlap(-1), recoed(-1), aid(0);
  Bool_t cascade(kFALSE), klong(kFALSE), kshort(kFALSE), neutron(kFALSE), lostk(kFALSE), kplus(kFALSE);
  TLorentzVector p(0., 0., 0., 0.), plostmax(0., 0., 0., 0.)
    , pkl(0., 0., 0., 0.), pks(0., 0., 0., 0.), pk(0., 0., 0., 0.), pn(0., 0., 0., 0.), pl(0., 0., 0., 0.)
    , pmiss(0., 0., 0., 0.);

  Int_t bla(0);
  for (int i = 0; i < nMc; ++i) {
    mk4Vector(p, pMc[i], thetaMc[i], phiMc[i], massMc[i]);
    // -- recoil B is #2 ????TRUE on generic as well ????
    if (isAncestor(fB1Index, i)) {
      overlap = 1;
    } else if (isAncestor(fB2Index, i)) {
      overlap = 2;
    } else {
      overlap = 0;
    }
    if (overlap != 2) continue;
    
    // -- Cascade
    if ((overlap == 2) && (isTruEl(i)||isTruMu(i)||isTruTau(i))) {
      aid = TMath::Abs(idMc[mothMc[i]-1]);
      int did = ((aid - (aid%100))/100)%10;
      if (TMath::Abs(did) == 4) {
	pl += p; 
	cascade = kTRUE;
      } 
    }
    // -- Charged Kaons not identified? 
    bla = isRecoed(i);
    if (bla > -1) {
      recoed = 1;
      if ((p.Vect().Mag() > 0.3) && (TMath::Abs(idMc[i]) == 321)) {
	Bool_t ka = isRecKaon(bla);
	if (!ka) {
	  lostk = kTRUE;
	  pk = p;
	}
      }
    } else {
      pmiss += p;
      if (p.Vect().Mag() > plostmax.Vect().Mag()) {
	plostmax = p; 
      }
      recoed = 0; 
    }
    // -- charged Kaons
    if ((overlap == 2) && (TMath::Abs(idMc[i]) == 321)) {
      kplus = kTRUE;
      pkl += p;
    }
    // -- KSHORT
    if ((overlap == 2) && (TMath::Abs(idMc[i]) == 310)) {
      kshort = kTRUE;
      pks += p;
    }
    // -- KLONG
    if ((overlap == 2) && (TMath::Abs(idMc[i]) == 130)) {
      klong = kTRUE;
      pkl += p;
    }
    // -- neutrons from HQ decays
    if ((overlap == 2) && (TMath::Abs(idMc[i]) == 2112)) {
      aid = TMath::Abs(idMc[mothMc[i]-1]);
      int did = TMath::Abs(((aid - (aid%100))/100)%10);
      if ((did == 4) || (did == 5)) {
	neutron = kTRUE;
	pn += p;
      }
    }
  }
  

  fMxCategory = 0;
  if (cascade)  fMxCategory += MCcascade;  // 1
  if (klong)    fMxCategory += MCklong;    // 2 
  if (kshort)   fMxCategory += MCkshort;   // 4
  if (neutron)  fMxCategory += MCneutron;  // 8
  if (lostk)    fMxCategory += MClostk;    // 16
  if (kplus)    fMxCategory += MCkplus;    // 32
  if (klong && cascade) fMxCategory+= MCcascadeklong; //64
  if (klong && !cascade) fMxCategory+= MCnotcascadeklong; //128
  if (kshort && cascade) fMxCategory+= MCcascadekshort; //256
  if (kshort && !cascade) fMxCategory+= MCnotcascadekshort; //512
  if (!klong && !kshort && cascade) fMxCategory+= MCcascadenotk; //1024
  

//    fHistFile->cd(dir);
//    if (vubDepleted) {sprintf(name, "c"); } else { sprintf(name, "u"); }
//    double lostPartMom = plostmax.Vect().Mag();
//    double xmassRes = fMxhad - fMxhadGen;
//    double xmassfitRes = fMxhadfit - fMxhadGen;

//    sprintf(title, "%s100", name); ((TH1D*)gDirectory->Get(title))->Fill(0.);
//    sprintf(title, "%s1000",name); ((TH2D*)gDirectory->Get(title))->Fill(0., xmassRes);
//    sprintf(title, "%s1001",name); ((TH2D*)gDirectory->Get(title))->Fill(0., xmassfitRes);
//    // -- KLONG
//    if (klong) {
//      sprintf(title, "%s100", name); ((TH1D*)gDirectory->Get(title))->Fill(1.);
//      sprintf(title, "%s1000",name); ((TH2D*)gDirectory->Get(title))->Fill(1., xmassRes);
//      sprintf(title, "%s1001",name); ((TH2D*)gDirectory->Get(title))->Fill(1., xmassfitRes);
//      sprintf(title, "%s101", name); ((TH1D*)gDirectory->Get(title))->Fill(pkl.Vect().Mag());
//      sprintf(title, "%s201", name); ((TH1D*)gDirectory->Get(title))->Fill(pkl.Mag());
//    }
//    // -- neutron
//    if (neutron) {
//      sprintf(title, "%s100", name); ((TH1D*)gDirectory->Get(title))->Fill(2.);
//      sprintf(title, "%s1000",name); ((TH2D*)gDirectory->Get(title))->Fill(2., xmassRes);
//      sprintf(title, "%s1001",name); ((TH2D*)gDirectory->Get(title))->Fill(2., xmassfitRes);
//      sprintf(title, "%s102", name); ((TH1D*)gDirectory->Get(title))->Fill(pn.Vect().Mag());
//      sprintf(title, "%s202", name); ((TH1D*)gDirectory->Get(title))->Fill(pn.Mag());
//    }
//    // -- lost particles
//    if (lostPartMom > 0.5) {
//      sprintf(title, "%s100", name); ((TH1D*)gDirectory->Get(title))->Fill(3.);
//      sprintf(title, "%s1000",name); ((TH2D*)gDirectory->Get(title))->Fill(3., xmassRes);
//      sprintf(title, "%s1001",name); ((TH2D*)gDirectory->Get(title))->Fill(3., xmassfitRes);
//      sprintf(title, "%s103", name); ((TH1D*)gDirectory->Get(title))->Fill(pmiss.Vect().Mag());
//      sprintf(title, "%s203", name); ((TH1D*)gDirectory->Get(title))->Fill(pmiss.Mag());
//    }
//    // -- cascade decays
//    if (cascade) {
//      sprintf(title, "%s100", name); ((TH1D*)gDirectory->Get(title))->Fill(4.);
//      sprintf(title, "%s1000",name); ((TH2D*)gDirectory->Get(title))->Fill(4., xmassRes);
//      sprintf(title, "%s1001",name); ((TH2D*)gDirectory->Get(title))->Fill(4., xmassfitRes);
//      sprintf(title, "%s104", name); ((TH1D*)gDirectory->Get(title))->Fill(pl.Vect().Mag());
//      sprintf(title, "%s204", name); ((TH1D*)gDirectory->Get(title))->Fill(pl.Mag());
//    }
//    // -- un-identified charged kaons
//    if (lostk) {
//      sprintf(title, "%s100", name); ((TH1D*)gDirectory->Get(title))->Fill(5.);
//      sprintf(title, "%s1000",name); ((TH2D*)gDirectory->Get(title))->Fill(5., xmassRes);
//      sprintf(title, "%s1001",name); ((TH2D*)gDirectory->Get(title))->Fill(5., xmassfitRes);
//      sprintf(title, "%s105", name); ((TH1D*)gDirectory->Get(title))->Fill(pk.Vect().Mag());
//      sprintf(title, "%s205", name); ((TH1D*)gDirectory->Get(title))->Fill(pk.Mag());
//    }
                  
}

// ----------------------------------------------------------------------
void recoilBase::fillMesHist(const char *dir) {

  fHistFile->cd(dir);
  Int_t l; 
  char le[10];
  for (l = 0; l < 3; ++l) {
    if (l == 0) {
      if (!fElectron) continue;
      sprintf(le, "e");
    }
    if (l == 1) {
      if (!fMuon) continue;
      sprintf(le, "m");
    }
    if (l == 2) {
      sprintf(le, "a");
    }    

    char name[200];
    sprintf(name, "mes%s%d", le, 1);  if (fGoodAccLepton) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
    sprintf(name, "mes%s%d", le, 2);  if (fGoodChargeCorr) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
    sprintf(name, "mes%s%d", le, 3);  if (fGoodLepton) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
    sprintf(name, "mes%s%d", le, 4);  if (fOneLepton) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
    sprintf(name, "mes%s%d", le, 5);  if (fGoodChargeCons) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
    sprintf(name, "mes%s%d", le, 6);  if (fGoodMM2) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
    sprintf(name, "mes%s%d", le, 7);  if (fGoodEvent) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
    
    sprintf(name, "mes%s%d", le, 11);      
    if (fGoodAccLepton) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
    sprintf(name, "mes%s%d", le, 12);      
    if (fGoodAccLepton&&fGoodChargeCorr) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
    sprintf(name, "mes%s%d", le, 13);      
    if (fGoodAccLepton&&fGoodChargeCorr&&fGoodLepton) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
    sprintf(name, "mes%s%d", le, 14);      
    if (fGoodAccLepton&&fGoodChargeCorr&&fGoodLepton&&fOneLepton) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
    sprintf(name, "mes%s%d", le, 15);      
    if (fGoodAccLepton&&fGoodChargeCorr&&fGoodLepton&&fOneLepton&&fGoodChargeCons) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
    sprintf(name, "mes%s%d", le, 16);      
    if (fGoodAccLepton&&fGoodChargeCorr&&fGoodLepton&&fOneLepton&&fGoodChargeCons&&fGoodMM2) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
    sprintf(name, "mes%s%d", le, 17);      
    if (fGoodEvent) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
  }
  fHistFile->cd();
  return;
}

// ----------------------------------------------------------------------
void recoilBase::fillRecoilHist(const char *dir, int chbcand) {
  
  fHistFile->cd(dir);

  Int_t i, l; 
  double pnuRes = fPNu - p4MissGen.Vect().Mag();
  double tnuRes = TMath::Cos(fTNu) - TMath::Cos(p4MissGen.Theta());
  double q2Res  = fQ2 - fQ2Gen;
  
  double xmassRes    = fMxhad - fMxhadGen;
  double xmassfitRes = fMxhadfit - fMxhadGen;
  if (fOptGammas) {
    for(int jh =0; jh<15; jh++) {
      tmpxmassResF[jh] = tmpfMxhadfit[jh]- fMxhadGen; 
      tmpxmassRes[jh] = tmpfMxhad[jh]- fMxhadGen; 
    }
  }

  // -- Kshorts
  for (Int_t iks = 0; iks < nKs; ++iks) {
    double mass = massKs[iks];
    if (goodKshort[iks] == 0) continue;   // skip kshorts whose daughters make up better combinations
    ((TH1D*)gDirectory->Get("ks100"))->Fill(mass);
    if (fGoodLepton)     ((TH1D*)gDirectory->Get("ks101"))->Fill(mass);
    if (fGoodMM2)        ((TH1D*)gDirectory->Get("ks102"))->Fill(mass);
    if (fGoodChargeCorr) ((TH1D*)gDirectory->Get("ks103"))->Fill(mass);
    if (fGoodChargeCons) ((TH1D*)gDirectory->Get("ks104"))->Fill(mass);
    if (fGoodEvent)      ((TH1D*)gDirectory->Get("ks105"))->Fill(mass);
    
    ((TH1D*)gDirectory->Get("ks200"))->Fill(pKs[iks]);
    ((TH1D*)gDirectory->Get("ks201"))->Fill(thKs[iks]*DR);
  }

  // -- charged kaons
  for (i = 0; i < nTrk; ++i) {
    if (goodChargedKaon[i] == 0)  continue;
    if (fGoodEvent) {
      ((TH1D*)gDirectory->Get("kp200"))->Fill(momentumTrk[i]);
      ((TH1D*)gDirectory->Get("kp201"))->Fill(thetaTrk[i]*DR);
    }
  }

  char le[10];
  for (l = 0; l < 3; ++l) {
    if (l == 0) {
      if (!fElectron) continue;
      sprintf(le, "e");
    }
    if (l == 1) {
      if (!fMuon) continue;
      sprintf(le, "m");
    }
    if (l == 2) {
      sprintf(le, "a");
    }    

    for (i = 0; i < 7; ++i) {
      if ((i == 1) && (fVcb != 1)) continue;
      if ((i == 2) && (fVub != 1)) continue;
      if ((i == 3) && (fBVxbTyp != 1)) continue;
      if ((i == 4) && (fBVxbTyp != 2)) continue;
      if ((i == 5) && (fBVxbTyp != 3)) continue;
      if ((i == 6) && (fOther < 1)) continue;

      char name[100];

      //Photon study
      if (fOptGammas) {
	
	gDirectory->cd("CutPlots");
	
	for(int jp =0; jp<15; jp++) {
	  if (fGoodEventPhNMNC[jp]) {
	    sprintf(name, "bin_nc_Mnu%d", jp); ((TH1D*)gDirectory->Get(name))->Fill(tmpfMM2[jp],1.);
	    sprintf(name, "nc_Mnu%d", jp); ((TH1D*)gDirectory->Get(name))->Fill(tmpfMM2[jp],1.);
	    sprintf(name, "nc_ResMX%d", jp); ((TH1D*)gDirectory->Get(name))->Fill(tmpxmassRes[jp],1.);
	    sprintf(name, "bin_nc_ResMX%d", jp); ((TH1D*)gDirectory->Get(name))->Fill(tmpxmassRes[jp],1.);
	    sprintf(name, "nc_ResMXF%d", jp); ((TH1D*)gDirectory->Get(name))->Fill(tmpxmassResF[jp],1.);
	    sprintf(name, "bin_nc_ResMXF%d", jp); ((TH1D*)gDirectory->Get(name))->Fill(tmpxmassResF[jp],1.);
	  }
	}
	
	for(int jp =0; jp<15; jp++) {
	  if (fGoodEventPh[jp]) {
	    sprintf(name, "Mnu%d", jp); ((TH1D*)gDirectory->Get(name))->Fill(tmpfMM2[jp],1.);
	    sprintf(name, "bin_Mnu%d", jp); ((TH1D*)gDirectory->Get(name))->Fill(tmpfMM2[jp],1.);
	    sprintf(name, "ResMX%d", jp); ((TH1D*)gDirectory->Get(name))->Fill(tmpxmassRes[jp],1.);
	    sprintf(name, "bin_ResMX%d", jp); ((TH1D*)gDirectory->Get(name))->Fill(tmpxmassRes[jp],1.);
	    sprintf(name, "ResMXF%d", jp); ((TH1D*)gDirectory->Get(name))->Fill(tmpxmassResF[jp],1.);
	    sprintf(name, "bin_ResMXF%d", jp); ((TH1D*)gDirectory->Get(name))->Fill(tmpxmassResF[jp],1.);
	    
	    for(int io =0; io<10; io++) {
	      if (fGoodNoHole[(io+1)]) {
		sprintf(name,"MNH%d%d", jp,io); ((TH1D*)gDirectory->Get(name))->Fill(tmpfMM2[jp],1.);
		sprintf(name,"bin_MNH%d%d", jp,io); ((TH1D*)gDirectory->Get(name))->Fill(tmpfMM2[jp],1.);   
		sprintf(name,"ResNH%d%d", jp,io); ((TH1D*)gDirectory->Get(name))->Fill(tmpxmassRes[jp],1.);
		sprintf(name,"bin_ResNH%d%d", jp,io); ((TH1D*)gDirectory->Get(name))->Fill(tmpxmassRes[jp],1.);
		sprintf(name,"MxNH%d%d", jp,io); ((TH1D*)gDirectory->Get(name))->Fill(tmpfMxhad[jp],1.);
		sprintf(name,"bin_MxNH%d%d", jp,io); ((TH1D*)gDirectory->Get(name))->Fill(tmpfMxhad[jp],1.);
		sprintf(name,"ResNHF%d%d", jp,io); ((TH1D*)gDirectory->Get(name))->Fill(tmpxmassResF[jp],1.);
		sprintf(name,"bin_ResNHF%d%d", jp,io); ((TH1D*)gDirectory->Get(name))->Fill(tmpxmassResF[jp],1.);
		sprintf(name,"MxNHF%d%d", jp,io); ((TH1D*)gDirectory->Get(name))->Fill(tmpfMxhadfit[jp],1.);
		sprintf(name,"bin_MxNHF%d%d", jp,io); ((TH1D*)gDirectory->Get(name))->Fill(tmpfMxhadfit[jp],1.);
	      }
	      
	    }     
	  }
	} 
      }

      fHistFile->cd(dir);    
      gDirectory->cd();
      

      ((TH1D*)gDirectory->Get("nc_mNuSqNC"))->Fill(fMM2NC);
      
      if (fGoodAccLepton) {
	sprintf(name, "%s%d", le, 1000+i);  ((TH1D*)gDirectory->Get(name))->Fill(fPcms, 1.);
	sprintf(name, "%s%d", le, 1100+i);   if (fIsPrompt) ((TH1D*)gDirectory->Get(name))->Fill(fPcms, 1.);
	sprintf(name, "%s%d", le, 1200+i);   if (fIsCascade) ((TH1D*)gDirectory->Get(name))->Fill(fPcms, 1.);
	sprintf(name, "%s%d", le, 1300+i);   ((TH1D*)gDirectory->Get(name))->Fill(TMath::Cos(fTcms), 1.);
	sprintf(name, "%s%d", le, 1400+i);   ((TH1D*)gDirectory->Get(name))->Fill(fPups, 1.);

	// -- all other cuts
	if (faoLepton) {
	  sprintf(name, "%s%d", le, 1020+i);  ((TH1D*)gDirectory->Get(name))->Fill(fPcms, 1.);
	  sprintf(name, "%s%d", le, 1120+i);   if (fIsPrompt) ((TH1D*)gDirectory->Get(name))->Fill(fPcms, 1.);
	  sprintf(name, "%s%d", le, 1220+i);   if (fIsCascade) ((TH1D*)gDirectory->Get(name))->Fill(fPcms, 1.);
	  sprintf(name, "%s%d", le, 1320+i);   ((TH1D*)gDirectory->Get(name))->Fill(TMath::Cos(fTcms), 1.);
	  sprintf(name, "%s%d", le, 1420+i);   ((TH1D*)gDirectory->Get(name))->Fill(fPups, 1.);
	}
	if (faoMM2) {
	  sprintf(name, "%s%d", le, 3020+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMM2);
	  if (fEffCat > 0) sprintf(name, "%s%d", le, 3020 + fEffCat*100 + i);   ((TH1D*)gDirectory->Get(name))->Fill(fMM2);
	  //  	  sprintf(name, "%s%d", le, 3120+i);   ((TH1D*)gDirectory->Get(name))->Fill(TMath::Cos(fTNu));
	  //  	  sprintf(name, "%s%d", le, 3220+i);   ((TH1D*)gDirectory->Get(name))->Fill(fPNu);
	  sprintf(name, "%s%d", le, 3820+i);   ((TH1D*)gDirectory->Get(name))->Fill(fQ2);
	  sprintf(name, "%s%d", le, 9120+i);   ((TH1D*)gDirectory->Get(name))->Fill(pnuRes);
	  sprintf(name, "%s%d", le, 9220+i);   ((TH1D*)gDirectory->Get(name))->Fill(tnuRes);
	  sprintf(name, "%s%d", le, 9320+i);   ((TH1D*)gDirectory->Get(name))->Fill(q2Res);
	}
	if (faoChargeCons) {
	  sprintf(name, "%s%d", le, 4020+i);   ((TH1D*)gDirectory->Get(name))->Fill(fRecoilTrkMult);
	  sprintf(name, "%s%d", le, 4120+i);   ((TH1D*)gDirectory->Get(name))->Fill(fRecoilNutMult);
	  sprintf(name, "%s%d", le, 4220+i);   ((TH1D*)gDirectory->Get(name))->Fill(fRecoilCharge);
	  sprintf(name, "%s%d", le, 4320+i);   ((TH1D*)gDirectory->Get(name))->Fill(fRecoilCharge+fBrecoCharge);
	}
	if (fGoodEvent) { // this is just for consistency
	  sprintf(name, "%s%d", le, 2020+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhad);
	  sprintf(name, "%s%d", le, 2120+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhad);
	  sprintf(name, "%s%d", le, 2220+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhadGen);
	  sprintf(name, "%s%d", le, 2320+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhad*fMxhad);
	  sprintf(name, "%s%d", le, 2420+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhadfit);
	  sprintf(name, "%s%d", le, 2520+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhadfit);
	  sprintf(name, "%s%d", le, 3920+i);   ((TH1D*)gDirectory->Get(name))->Fill(fGammaMax);
	  sprintf(name, "%s%d", le, 9020+i);   ((TH1D*)gDirectory->Get(name))->Fill(xmassRes);
	  sprintf(name, "%s%d", le, 9420+i);   ((TH1D*)gDirectory->Get(name))->Fill(xmassfitRes);
	}
      }

      if (fGoodLepton) {
	sprintf(name, "%s%d", le, 2000+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhad);
	sprintf(name, "%s%d", le, 2100+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhad);
	sprintf(name, "%s%d", le, 2200+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhadGen);
	sprintf(name, "%s%d", le, 2300+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhad*fMxhad);
	sprintf(name, "%s%d", le, 2400+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhadfit);
	sprintf(name, "%s%d", le, 2500+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhadfit);
	sprintf(name, "%s%d", le, 3000+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMM2);
	if (fEffCat > 0) sprintf(name, "%s%d", le, 3000 + fEffCat*100 + i);   ((TH1D*)gDirectory->Get(name))->Fill(fMM2);
	sprintf(name, "%s%d", le, 3800+i);   ((TH1D*)gDirectory->Get(name))->Fill(fQ2);
        sprintf(name, "%s%d", le, 3900+i);   ((TH1D*)gDirectory->Get(name))->Fill(fGammaMax);
	sprintf(name, "%s%d", le, 4000+i);   ((TH1D*)gDirectory->Get(name))->Fill(fRecoilTrkMult);
	sprintf(name, "%s%d", le, 4100+i);   ((TH1D*)gDirectory->Get(name))->Fill(fRecoilNutMult);
	sprintf(name, "%s%d", le, 4200+i);   ((TH1D*)gDirectory->Get(name))->Fill(fRecoilCharge);
	sprintf(name, "%s%d", le, 4300+i);   ((TH1D*)gDirectory->Get(name))->Fill(fRecoilCharge+fBrecoCharge);
	
	sprintf(name, "%s%d", le, 9000+i);   ((TH1D*)gDirectory->Get(name))->Fill(xmassRes);
	sprintf(name, "%s%d", le, 9100+i);   ((TH1D*)gDirectory->Get(name))->Fill(pnuRes);
	sprintf(name, "%s%d", le, 9200+i);   ((TH1D*)gDirectory->Get(name))->Fill(tnuRes);
	sprintf(name, "%s%d", le, 9300+i);   ((TH1D*)gDirectory->Get(name))->Fill(q2Res);
	sprintf(name, "%s%d", le, 9400+i);   ((TH1D*)gDirectory->Get(name))->Fill(xmassfitRes);
          
	if (fOptCategories) {
	  Int_t bit=1;
	  for (Int_t k=1; k<=numOfCategories; k++) {
	    Bool_t goodEv;
	    if (k==3 || k==9 || k==10)
	      goodEv = fGoodEventKS;
	    else
	      goodEv = fGoodEvent;
	    
	    if (goodEv) {
	      if ( fMxCategory & bit ) {
		sprintf(name, "%s%d", le, 5000+i+k*10);
		((TH1D*)gDirectory->Get(name))->Fill(fCatMxhad[k]);
		sprintf(name, "%s%d", le, 7000+i+k*10);
		((TH1D*)gDirectory->Get(name))->Fill(fCatMxhadFit[k]);
		sprintf(name, "%s%d", le, 5000+i);
		((TH1D*)gDirectory->Get(name))->Fill(k);
		if (fCatMxhad[k]<1.5) {
		  sprintf(name, "%s%d", le, 6000+i);
		  ((TH1D*)gDirectory->Get(name))->Fill(k);
		}
		if (fCatMxhadFit[k]<1.5) {
		  sprintf(name, "%s%d", le, 8000+i);
		  ((TH1D*)gDirectory->Get(name))->Fill(k);
		}
	      }
	      bit = bit << 1 ;
	    }
	    
	    if (fMxhad > 0) {
	      sprintf(name, "%s%d", le, 5600+i);
	      ((TH1D*)gDirectory->Get(name))->Fill(fMxCategory);
	      sprintf(name, "%s%d", le, 7600+i);
	      ((TH1D*)gDirectory->Get(name))->Fill(fMxCategory);
	    }
	    
	    if (fNKshort>0) {
	      sprintf(name, "%s%d", le, 5500+i);
	      ((TH1D*)gDirectory->Get(name))->Fill(fCatMxhad[3]);
	      sprintf(name, "%s%d", le, 7500+i);
	      ((TH1D*)gDirectory->Get(name))->Fill(fCatMxhadFit[3]);
	      
	      if (0) {
		cout << "step i: " << i << endl;
		cout << "fCatMxhad[3] = " << fCatMxhad[3] << endl;
		cout << "fMxhad = " << fMxhad << endl;
		cout << "fnkshort" << fNKshort << endl;
	      }
	    }
	  }
	}
      }

      // -- after cuts
      if (fGoodEvent) {
	((TH1D*)gDirectory->Get("mNuSqNC"))->Fill(fMM2NC);
        sprintf(name, "%s%d", le, 1010+i);   ((TH1D*)gDirectory->Get(name))->Fill(fPcms);
	sprintf(name, "%s%d", le, 1110+i);   if (fIsPrompt) ((TH1D*)gDirectory->Get(name))->Fill(fPcms);
	sprintf(name, "%s%d", le, 1210+i);   if (fIsCascade) ((TH1D*)gDirectory->Get(name))->Fill(fPcms);
	sprintf(name, "%s%d", le, 1310+i);   ((TH1D*)gDirectory->Get(name))->Fill(TMath::Cos(fTcms));
	sprintf(name, "%s%d", le, 1410+i);   ((TH1D*)gDirectory->Get(name))->Fill(fPups);

        sprintf(name, "%s%d", le, 2010+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhad);
	sprintf(name, "%s%d", le, 2110+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhad);
        sprintf(name, "%s%d", le, 2210+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhadGen);
        sprintf(name, "%s%d", le, 2310+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhad*fMxhad);
	sprintf(name, "%s%d", le, 2410+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhadfit);
	sprintf(name, "%s%d", le, 2510+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhadfit);

        sprintf(name, "%s%d", le, 3010+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMM2);
        if (fEffCat > 0) sprintf(name, "%s%d", le, 3010 + fEffCat*100 + i);   ((TH1D*)gDirectory->Get(name))->Fill(fMM2);
 	//          sprintf(name, "%s%d", le, 3110+i);   ((TH1D*)gDirectory->Get(name))->Fill(TMath::Cos(fTNu));
	//          sprintf(name, "%s%d", le, 3210+i);   ((TH1D*)gDirectory->Get(name))->Fill(fPNu);
        sprintf(name, "%s%d", le, 3810+i);   ((TH1D*)gDirectory->Get(name))->Fill(fQ2);
        sprintf(name, "%s%d", le, 3910+i);   ((TH1D*)gDirectory->Get(name))->Fill(fGammaMax);

        sprintf(name, "%s%d", le, 4010+i);   ((TH1D*)gDirectory->Get(name))->Fill(fRecoilTrkMult);
        sprintf(name, "%s%d", le, 4110+i);   ((TH1D*)gDirectory->Get(name))->Fill(fRecoilNutMult);
	sprintf(name, "%s%d", le, 4210+i);   ((TH1D*)gDirectory->Get(name))->Fill(fRecoilCharge);
	sprintf(name, "%s%d", le, 4310+i);   ((TH1D*)gDirectory->Get(name))->Fill(fRecoilCharge+fBrecoCharge);

        sprintf(name, "%s%d", le, 9010+i);   ((TH1D*)gDirectory->Get(name))->Fill(xmassRes);
        sprintf(name, "%s%d", le, 9110+i);   ((TH1D*)gDirectory->Get(name))->Fill(pnuRes);
        sprintf(name, "%s%d", le, 9210+i);   ((TH1D*)gDirectory->Get(name))->Fill(tnuRes);
        sprintf(name, "%s%d", le, 9310+i);   ((TH1D*)gDirectory->Get(name))->Fill(q2Res);
        sprintf(name, "%s%d", le, 9410+i);   ((TH1D*)gDirectory->Get(name))->Fill(xmassfitRes);
      }
    }
  }

  fHistFile->cd();
}

// ----------------------------------------------------------------------
Int_t recoilBase::compChgBreco( int nBs , int chbcand ) {
  int chg(0);
  int tmpdBLUND[7] = { d1B0Lund[nBs],
		       d2B0Lund[nBs],
		       d3B0Lund[nBs],
		       d4B0Lund[nBs],
		       d5B0Lund[nBs],
		       d6B0Lund[nBs],
		       d7B0Lund[nBs] };
  if(chbcand) {
    tmpdBLUND[0] =  d1ChBLund[nBs];
    tmpdBLUND[1] =  d2ChBLund[nBs];
    tmpdBLUND[2] =  d3ChBLund[nBs];
    tmpdBLUND[3] =  d4ChBLund[nBs];
    tmpdBLUND[4] =  d5ChBLund[nBs];
    tmpdBLUND[5] =  d6ChBLund[nBs];
    tmpdBLUND[6] =  d7ChBLund[nBs] ;
  }

  for(int i = 0; i<7; i++) {

    if(tmpdBLUND[i] == 211  || tmpdBLUND[i] == 321   || tmpdBLUND[i] == 411  || tmpdBLUND[i] == 413)     chg++;
    if(tmpdBLUND[i] == -211 || tmpdBLUND[i] == -321  || tmpdBLUND[i] == -411 || tmpdBLUND[i] == -413)    chg--;
    // if(TMath::Abs(tmpdBLUND[i]) != 211 && TMath::Abs(tmpdBLUND[i]) != 321 && TMath::Abs(tmpdBLUND[i]) != 411 && TMath::Abs(tmpdBLUND[i]) != 413 ) cout << tmpdBLUND[i] << endl;
  }
  return chg;
}

// ----------------------------------------------------------------------
void recoilBase::mcTruth( int chbcand ) {

  if (!fIsMC) return;
  
  Bool_t isSemilep(kFALSE);
  Int_t ipB[] = {-1, -1};
  Int_t typB[] = {7, 7};
  Int_t cnt(-1);

  Int_t Imc(-1), Jmc(-1), Kmc(-1); 
  fB1Index = fB2Index = fBVxb = -99;

  // -- Find indices of B's
  for (Int_t imc = 0; imc < nMc; ++imc) {
    if ((TMath::Abs(idMc[imc]) == 511) || (TMath::Abs(idMc[imc]) == 521)) {
      cnt++;
      ipB[cnt] = imc;
      if (fB1Index < 0) {
	fB1Index = imc; 
      } else {
	fB2Index = imc; 
      }
    }
    if (cnt == 1) break;
  }

  // -- Determine event type
  genLeptonCharge = 0;
  Int_t ib;
  Bool_t isCharm(kFALSE);
  for (ib = 0; ib < 2; ++ib) { 
    for (Int_t imc = 0; imc < nMc; ++imc) {
      if ((mothMc[imc]-1) == ipB[ib]) {

	// -- D*0
	if (TMath::Abs(idMc[imc]) == 423) {
	  isSemilep = kFALSE;
	  for (Int_t s = 0; s < nMc; ++s) {
	    if (((mothMc[s]) == (mothMc[imc])) && isTruLepton(s)) {
	      isSemilep = kTRUE;
	      break;
	    }
	  }
	  if (!isSemilep) {
	    //	    cout << "Found type = 5 had B decay with s = " << s << " -> " idMc[s] << endl;
	    typB[ib] = 5;
	  } else {
	    //	    cout << "Found type = 5 sl B decay with s = " << s << " -> " idMc[s] << endl;
	  }
	} // D*0

	// -- D0
	if (TMath::Abs(idMc[imc]) == 421) {
	  isSemilep = kFALSE;
	  for (Int_t s = 0; s < nMc; ++s) {
	    if ((mothMc[s] == mothMc[imc]) && isTruLepton(s)) {
	      isSemilep = kTRUE;
	      break;
	    }
	  }
	  if (!isSemilep) {
	    //	    cout << "Found type = 6 had B decay with s = " << s << " -> " idMc[s] << endl;
	    typB[ib] = 6;
	  } else {
	    //	    cout << "Found type = 6 sl B decay with s = " << s << " -> " idMc[s] << endl;
	  }
	} // D0

	// -- D*+
	if (TMath::Abs(idMc[imc]) == 413) {
	  isSemilep = kFALSE;
	  for (Int_t s = 0; s < nMc; ++s) {
  	    if ((mothMc[s] == mothMc[imc]) && isTruLepton(s)) {
  	      isSemilep = kTRUE;
	      break;
  	    }
	  }
	  if (!isSemilep) {
	    //	    cout << "Found type = 3 had B decay with s = " << s << " -> " idMc[s] << endl;
	    typB[ib] = 3;
	  } else {
	    //	    cout << "Found type = 3 sl B decay with s = " << s << " -> " idMc[s] << endl;
	  }
	} // D*+

	// -- D+
	if (TMath::Abs(idMc[imc]) == 411) {
	  isSemilep = kFALSE;
	  for (Int_t s = 0; s < nMc; ++s) {
	    if ((mothMc[s] == mothMc[imc]) && isTruLepton(s)) {
	      isSemilep = kTRUE;
	      break;
	    }
	  }
	  if (!isSemilep) {
	    //	    cout << "Found type = 4 had B decay with s = " << s << " -> " idMc[s] << endl;
	    typB[ib] = 4;
	  } else {
	    //	    cout << "Found type = 4 sl B decay with s = " << s << " -> " idMc[s] << endl;
	  }
	} // D+


	// -- sl Decay: Vcb or Vub? 
	if (isTruLepton(imc)) {
	  Imc = imc;
	  Int_t idmc(0);
	  genLeptonCharge = -1*idMc[imc]/TMath::Abs(idMc[imc]);
	  typB[ib] = 1; 
	  fBVxb = ipB[ib];
	  for (Int_t s = 0; s < nMc; ++s) {
	    idmc = TMath::Abs(idMc[s]);
	    if (mothMc[s] == mothMc[imc]) {
	      if (Jmc < 0)  Jmc = s; 
	      else if (Kmc < 0)  Kmc = s; 
	    }
	    if ((mothMc[s] == mothMc[imc]) 
		&& (((idmc >= 400) && (idmc < 500))
		    || ((idmc >= 10400) && (idmc < 10500))
		    || ((idmc >= 20400) && (idmc < 20500))
		    || ((idmc >= 30400) && (idmc < 30500))
		    )) {
	      typB[ib] = 2; 
	      isCharm = kTRUE; 
	      break;
	    }
	  }
	}
      }
    }
  }

  // -- Determine special modes of sl B decay: 
  //    1 D+ e nu (gamma)
  //    2 Dstar e nu (gamma)
  //    3 D{+, *, **} pi e nu (gamma)
  //    4 VUB

  // THE FOLLOWING IS WRONG BECAUSE IT DOES NOT ACCOUNT FOR POSSIBLE BREMSSTRAHLUNG. 
  // YOU NEED TO REDUCE THE NUMBER OF B_DAUGHTERS BY ONE IF ONE OF THE DAUGHTERS IS A PHOTON.
  fBVxbTyp = -1; 
  for (Int_t ss = 0; ss < nMc; ++ss) {
    if (mothMc[ss]-1 == fBVxb) {
      int idmc = TMath::Abs(idMc[ss]);
      if (((idmc == 411) || (idmc == 421)) &&(nDauMc[fBVxb] == 3)) {
	fBVxbTyp = 1; 
	break;
      }
      else if (((idmc == 413) || (idmc == 423)) &&  (nDauMc[fBVxb] == 3)) {
	fBVxbTyp = 2; 
	break;
      }
      else if (isCharm) {
	  fBVxbTyp = 3; 
      }	
      if (idmc == 41) {
	fBVxbTyp = 4; 
      }	
    }
  }

//    if (fBVxbTyp > 0) {
//      cout << "---------------------------------------------------------------------- " << fBVxbTyp << " " << fBVxb << endl;
//      dumpGeneratorBlock();
//    }

  //    cout << "The two B are at " << fB1Index << " and " << fB2Index << " and sl decay of " << fBVxb << endl;
  //    if (fBVxb > -1) {
  //      dumpGeneratorBlock();
  //    }
  
  Int_t nfullyrecoDstar0(0), nfullyrecoDstar(0), nfullyrecoDc(0), nfullyrecoD0(0);
  fVub = fVcb = fOther = 0; 
  for (Int_t j = 0; j < 2; ++j) {
    if (typB[j] == 5) {
      nfullyrecoDstar0++;
    } else if (typB[j] == 6) {
      nfullyrecoD0++;
    } else if (typB[j] == 3) {
      nfullyrecoDstar++;
    } else if (typB[j] == 4) {
      nfullyrecoDc++;
    } else if (typB[j] == 2) {
      //      cout << "idMc[Imc] = " << idMc[Imc] << " mom = " << idMc[mothMc[Imc]-1] 
      //           << " sisters " << idMc[Jmc] << "  " <<idMc[Kmc] << endl;
      fVcb++;
    } else if (typB[j] == 1) {
      //      cout << "vub idMc[Imc] = " << idMc[Imc]<<" mom = "<<idMc[mothMc[Imc]-1]
      //           <<" sisters " << idMc[Jmc] << " " <<idMc[Kmc] << endl;
      //        if (idMc[Jmc] == 41) {
      //  	dumpGeneratorBlock(); 
      //        }
      fVub++;
    } else {
      fOther++;
    }
  }

  fHistFile->cd("mcTruth");
  ((TH1D*)gDirectory->Get("h100"))->Fill(fVub, 1.);
  ((TH1D*)gDirectory->Get("h101"))->Fill(fVcb, 1.);
  ((TH1D*)gDirectory->Get("h102"))->Fill(nfullyrecoDstar, 1.);
  ((TH1D*)gDirectory->Get("h103"))->Fill(nfullyrecoDc, 1.);
  ((TH1D*)gDirectory->Get("h104"))->Fill(nfullyrecoDstar0, 1.);
  ((TH1D*)gDirectory->Get("h105"))->Fill(nfullyrecoD0, 1.);
  ((TH1D*)gDirectory->Get("h106"))->Fill(fOther, 1.);
  ((TH1D*)gDirectory->Get("h107"))->Fill(fBVxbTyp, 1.);
  ((TH2D*)gDirectory->Get("h9009"))->Fill(typB[0], typB[1]);
  

  // -- BRECO, recoiling B, Y(4S)
  int seedMode(1);
  TLorentzVector p4Null(0., 0., 0., 0.);
  TLorentzVector pMcBreco(p4Null), pMcRecoil(p4Null), pMcUpsilon(p4Null);
  TLorentzVector p(p4Null); 
  Int_t nElectron(0), nMuon(0), nTau(0); 
  Int_t icountmc(0);

  fNLeptonMC = 0;
  p4LeptonGen = p4Null;
  p4MissGen.SetXYZM(0., 0., 0., 0.);
  p4XhadGen.SetXYZM(0., 0., 0., 0.);
  p4RecoilGen.SetXYZM(0., 0., 0., 0.);

  for (int imc = 0; imc < nMc; ++imc) {
    if (idMc[imc] == 70553) {
      mk4Vector(pMcUpsilon, pMc[imc], thetaMc[imc], phiMc[imc], massMc[imc]);
      break;
    }
  }

  for (ib = 0; ib < 2; ++ib) {
    // -- Breco
    if (typB[ib] >= seedMode+2) {
      int ibk = ipB[ib];
      mk4Vector(pMcBreco, pMc[ibk], thetaMc[ibk], phiMc[ibk], massMc[ibk]);
    }
    // -- recoil
    if ((typB[ib] == 1) || (typB[ib] == 2)) {
      int ibk = ipB[ib];
      mk4Vector(pMcRecoil, pMc[ibk], thetaMc[ibk], phiMc[ibk], massMc[ibk]);
    }
  }

  // -- Recoil
  for (ib = 0; ib < 2; ++ib) {
    if (typB[ib] > 2) continue;
    
    TVector3 boostVector(-pMcRecoil.Vect());  
    boostVector.SetMag(boostVector.Mag()/pMcRecoil.E());

    //    cout<<    ipB[ib] << " <- B index ; nMC particles ->  "<< nMc << endl;
    for (int imc = 0; imc < nMc; ++imc) {
      if (mothMc[imc]-1 == ipB[ib]) {
	mk4Vector(p, pMc[imc], thetaMc[imc], phiMc[imc], massMc[imc]);
	p.Boost(boostVector); 
	
	if (isTruLepton(imc)) {
	  fNLeptonMC++;
	  p4LeptonGen = p;
	}
	if (isTruEl(imc)) nElectron++;
	if (isTruMu(imc)) nMuon++;
	if (isTruTau(imc)) nTau++;
	
	int idmc = TMath::Abs(idMc[imc]);
	if ((idmc < 11) || (idmc > 16)) {   // if not a lepton nor a neutrino
	  ((TH1D*)gDirectory->Get("h700"))->Fill(idMc[imc], 1.);
	  p4XhadGen += p;
	  icountmc++;
	}
	if ((idmc != 12) && (idmc != 14) && (idmc != 16)) p4RecoilGen += p;
 	if ((idmc == 12) || (idmc == 14) || (idmc == 16)) p4MissGen += p;
      } 
    }
  }
  
  fQ2Gen   = p4MissGen*p4LeptonGen;

  fPcmsGen = p4LeptonGen.Vect().Mag();
  fTcmsGen = p4LeptonGen.Theta();
  fFcmsGen = p4LeptonGen.Phi();
  fEcmsGen = p4LeptonGen.E();

  fPxhadGen = p4XhadGen.Vect().Mag();
  fTxhadGen = p4XhadGen.Theta();
  fFxhadGen = p4XhadGen.Phi();
  fExhadGen = p4XhadGen.E();
  fMxhadGen = p4XhadGen.Mag();

  if (fVub + fVcb > 0) {
    ((TH1D*)gDirectory->Get("h77000"))->Fill(fNLeptonMC, 1.);
    if (fNLeptonMC == 1) {
      ((TH1D*)gDirectory->Get("h121000"))->Fill(p4LeptonGen.Vect().Mag(), 1.);
      ((TH1D*)gDirectory->Get("h123000"))->Fill(fMxhadGen, 1.);

      if (fVcb == 1) {
	((TH1D*)gDirectory->Get("h121005"))->Fill(p4LeptonGen.Vect().Mag(), 1.);
	((TH1D*)gDirectory->Get("h123005"))->Fill(fMxhadGen, 1.);
	((TH1D*)gDirectory->Get("h124005"))->Fill(icountmc, 1.);
      }      
      
      if (fVub == 1) {
	((TH1D*)gDirectory->Get("h121006"))->Fill(p4LeptonGen.Vect().Mag(), 1.);
	((TH1D*)gDirectory->Get("h123006"))->Fill(fMxhadGen, 1.);
	((TH1D*)gDirectory->Get("h124006"))->Fill(icountmc, 1.);
      }      
      
    }
  }

}

// ----------------------------------------------------------------------
void recoilBase::dumpGeneratorBlock(int b1, int b2) {
  char line[200];
  if (b1 < 0) {
    for (int i = 0; i < nMc; ++i) {
      sprintf(line, "%3d %+6d mom(%3d) ndau(%3d) p=%5.3f, t=%5.3f f=%+5.3f v=(%+7.3f,%+7.3f,%+7.3f)", 
	      i, idMc[i], mothMc[i]-1, nDauMc[i],
	      pMc[i], thetaMc[i], phiMc[i],
	      xMc[i], yMc[i], zMc[i]);
      cout << line << endl;
    }
  } else {
    int overlap(-1);
    char recoed[2];
    for (int i = 0; i < nMc; ++i) {
      if (isAncestor(b1, i)) {
	overlap = 1;
      } else if (isAncestor(b2, i)) {
	overlap = 2;
      } else {
	overlap = 0;
      }
      if (isRecoed(i) > -1) {
	sprintf(recoed, "r");
      } else {
	sprintf(recoed, " ");
      }
      sprintf(line, "%3d %+6d mom(%3d) flag(%1d%s) p=%5.3f, t=%5.3f f=%+5.3f v=(%+7.3f,%+7.3f,%+7.3f)", 
	      i, idMc[i], mothMc[i]-1, overlap, recoed,
	      pMc[i], thetaMc[i], phiMc[i],
	      xMc[i], yMc[i], zMc[i]);
      cout << line << endl;
    }
  }
}

// ----------------------------------------------------------------------
void recoilBase::Loop(Int_t maxEvent, Int_t startEvent, Int_t isVerbose, Int_t lun) {
  int step(1000);
  int ischB(0);
  findPro = findUps = 0;
  // 

  double tmpMassPB, tmpMassThetaB, tmpMassPhiB ;
  double tmpPB, tmpThetaB, tmpPhiB ;
  double tmpPgen, tmpThetagen, tmpPhigen, tmpMassgen ;
  double tmpMB, tmpBevM;
  //


  if (fChain == 0) return;
  Int_t nentries = Int_t(fChain->GetEntries());
  if (maxEvent == 0) maxEvent = nentries;
  if (nentries < 1) {
    cout << "Found no entries in " << fChain->GetName() << endl;
  } else {
    cout << "Found " << nentries << " entries in tree " << fChain->GetName() << endl;
  }

  if (startEvent > 0) {
    cout << "Will start at event " << startEvent << endl;
    if (startEvent+maxEvent >  nentries) {
      cout << "Requested " << maxEvent << " events, but will run only to end of chain"  << endl;
      maxEvent = nentries - startEvent; 
    }
  }

  Int_t nbytes = 0, nb = 0;

  for (Int_t jentry = startEvent; jentry < startEvent+maxEvent; jentry++) {
    if (isVerbose) cout << "->  new event " << endl;
    fEvent = jentry;
    // in case of a TChain, ientry is the entry number in the current file
    tsdump = kFALSE;
    Int_t ientry = LoadTree(jentry); 
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%step == 0) cout << "Event " << jentry << endl;
    if (ientry == 0) cout << "File " << fChain->GetCurrentFile()->GetName() << endl;


     // --jump the events (reco Bch) with wrong B-D flavor correlation
    if(nB0 == 0 && !(modeChB[indexbestB]-14000>299 && modeChB[indexbestB]-14000<399) && !(modeChB[indexbestB]-15000>299 && modeChB[indexbestB]-15000<399) && !(modeChB[indexbestB]-11000>299 && modeChB[indexbestB]-11000<399)){
      int fBrecoChargetmp= 0;
      if(d2ChBLund[indexbestB]!=0&&d2ChBLund[indexbestB]!=111&&d2ChBLund[indexbestB]!=310) fBrecoChargetmp+=TMath::Abs(d2ChBLund[indexbestB])/d2ChBLund[indexbestB];
      if(d3ChBLund[indexbestB]!=0&&d3ChBLund[indexbestB]!=111&&d3ChBLund[indexbestB]!=310) fBrecoChargetmp+=TMath::Abs(d3ChBLund[indexbestB])/d3ChBLund[indexbestB];
      if(d4ChBLund[indexbestB]!=0&&d4ChBLund[indexbestB]!=111&&d4ChBLund[indexbestB]!=310) fBrecoChargetmp+=TMath::Abs(d4ChBLund[indexbestB])/d4ChBLund[indexbestB];
      if(d5ChBLund[indexbestB]!=0&&d5ChBLund[indexbestB]!=111&&d5ChBLund[indexbestB]!=310) fBrecoChargetmp+=TMath::Abs(d5ChBLund[indexbestB])/d5ChBLund[indexbestB];
      if(d6ChBLund[indexbestB]!=0&&d6ChBLund[indexbestB]!=111&&d6ChBLund[indexbestB]!=310) fBrecoChargetmp+=TMath::Abs(d6ChBLund[indexbestB])/d6ChBLund[indexbestB];
      if(d7ChBLund[indexbestB]!=0&&d7ChBLund[indexbestB]!=111&&d7ChBLund[indexbestB]!=310) fBrecoChargetmp+=TMath::Abs(d7ChBLund[indexbestB])/d7ChBLund[indexbestB];
      cout << "jump the events (reco Bch) with wrong B-D flavor correlation: " << fBrecoChargetmp << endl;
      if(-1*(TMath::Abs(d1ChBLund[indexbestB])/d1ChBLund[indexbestB]) != fBrecoChargetmp) {
	continue;
      }
    }
    
 
    if (fOptMakeEventList) {
      if (jentry%2 == 0) {
	//      cout << "--> Copying " << jentry << " into fToBeCopied " << bla << endl;
	fToBeCopied->Enter(jentry);
      }
    }

    if (nPi0 == 100) { 
      cout << "Event with pi0 maximum -- skipping due to Linux problem" << endl;
      continue;
    }

    if (nB0 == 0) {
      ischB = 1;
    }

    if ((nB0 + nChB) > 1) {
      if(massB0[indexbestB] > 0 ) {
 	ischB = 0;
      } else if (massChB[indexbestB] > 0 ) {
 	ischB = 1;
      } else {
 	cout <<massChB[indexbestB]<< " something fishy " << massB0[indexbestB] <<endl;
 	continue;
      }
    }

    // -- MonteCarlo Truth
    int brecoI(-99); 
    mcTruth(ischB);
    if (fIsMC) {
      if (fBVxb == fB1Index) {
	brecoI = fB2Index;
      } else {
	brecoI = fB1Index;
      }
      tmpPgen = pMc[brecoI];
      tmpThetagen = thetaMc[brecoI];
      tmpPhigen = phiMc[brecoI];
      tmpMassgen = massMc[brecoI];
      mk4Vector(p4BrecoGen, tmpPgen, tmpThetagen, tmpPhigen, tmpMassgen);
    }

    // -- Reco quantities
    if(ischB == 0) {
      tmpMassPB =MassPB0[indexbestB];
      tmpMassThetaB =MassThetaB0[indexbestB];
      tmpMassPhiB =MassPhiB0[indexbestB];
      tmpPB = pB0[indexbestB];
      tmpBevM = massB0[indexbestB];
      tmpThetaB =thB0[indexbestB];
      tmpPhiB =phiB0[indexbestB];
      tmpMB = BZMASS;
    } else {
      tmpMassPB =MassPChB[indexbestB];
      tmpMassThetaB =MassThetaChB[indexbestB];
      tmpMassPhiB =MassPhiChB[indexbestB];
      tmpPB = pChB[indexbestB];
      tmpBevM = massChB[indexbestB];
      tmpThetaB =thChB[indexbestB];
      tmpPhiB =phiChB[indexbestB];
      tmpMB = BPMASS;
    }
    mk4Vector(p4Breco, tmpMassPB , tmpMassThetaB, tmpMassPhiB, tmpMB); 
    mk4Vector(p4BrecoNC, tmpPB , tmpThetaB, tmpPhiB, tmpBevM); 

    // -- Some event quantities
    fRunnumber = runNumber;
    fLower = lowerID; 
    fUpper = upperID; 

    p4Brecoil = p4Upsilon - p4Breco; 
    cmsBoost.SetXYZ(p4Brecoil.Px(), p4Brecoil.Py(), p4Brecoil.Pz());
    cmsBoost.SetMag(cmsBoost.Mag()/p4Brecoil.E());

    breco(ischB, isVerbose);
    if (fOptCategories) mxCategory();
    maskKshorts(1, isVerbose);
    //     maskPi0(1);
    maskConversions();

    recoil(ischB, isVerbose, fisDuplicate);
    //    testConversion();
    
    if (fOptKlongs) KLStudy();

    if (fDump > 0) {
      if  (fDump & 1) {
	if (fOptGammas) fGTree->Fill();
      }
      if ((fDump & 2) && (fVub == 1)) {
	if (isVerbose) cout << " fill fTree 2" << endl;
	fTree->Fill();
      }
      if ((fDump & 4) && ((fPcms > 0.) || (fVub == 1) || (fVcb == 1))) {
	if (isVerbose) cout << " fill fTree 4" << endl;
	fTree->Fill();
      }
      if (fDump & 8) {
	if (isVerbose) cout << " fill fTree 8 " << endl;
	fTree->Fill();
      }

    }
    //    cout<<"new event"<<endl;
    
    if(tsdump) timestamp(lun);

    if( isVerbose ) {   cout << " jentry = " << jentry << "  mes = " << fMes << " pcms = " << fPcms << " mass = " << p4LeptonLab.Mag()
			     << " xhadmass = " << fMxhad << " mm2 = " << fMM2NC << " xhadfitmass = " << fMxhadfit
			     << " vub = " << fVub << "  vcb = " << fVcb << endl;
    }
  }
  
  if (fOptGammas)  cout<<findPro<<" "<<findUps<<endl;

}

// ----------------------------------------------------------------------
void recoilBase::LoopKill(const char *Ifile, Int_t maxEvent, Int_t startEvent, Int_t isVerbose, Int_t lun) {

  int step(1000);
  fDupcount = 0;

  if (fChain == 0) return;
  Int_t nentries = Int_t(fChain->GetEntries());
  if (maxEvent == 0) maxEvent = nentries;
  if (nentries < 1) {
    cout << "Found no entries in " << fChain->GetName() << endl;
  } else {
    cout << "Found " << nentries << " entries in tree " << fChain->GetName() << endl;
  }

  if (startEvent > 0) {
    cout << "Will start at event " << startEvent << endl;
    if (startEvent+maxEvent >  nentries) {
      cout << "Requested " << maxEvent << " events, but will run only to end of chain"  << endl;
      maxEvent = nentries - startEvent; 
    }
  }

  Int_t nbytes = 0, nb = 0;

  read_killTab(Ifile);

  for (Int_t jentry = startEvent; jentry < startEvent+maxEvent; jentry++) {
    if (isVerbose) cout << "->  new event " << endl;
    fEvent = jentry;
    // in case of a TChain, ientry is the entry number in the current file
    tsdump = kFALSE;
    Int_t ientry = LoadTree(jentry); 
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%step == 0) cout << "Event " << jentry << endl;
    if (ientry == 0) cout << "File " << fChain->GetCurrentFile()->GetName() << endl;

    fisDuplicate = kFALSE;
    //    cout<<upperID<< " Up ; Low "<<lowerID<<endl;

    if(Ifile != " ")  fisDuplicate = kill_dupli(upperID,lowerID);

    if (fOptMakeEventList) {
      if (!fisDuplicate) {
	//	cout << "--> Copying " << jentry << endl;
	fToBeCopied->Enter(jentry);
      } else {
	cout<<"Duplicate not copied"<<endl;
      }
    }

    fIVal = 0;

    const char *Arstr[64] = { 
      "1", "2" ,"3", "4",       
      "11", "12" ,"13", "14",   
      "21", "22" ,"23", "24",   
      "31", "32" ,"33", "34",
      "41", "42" ,"43", "44",
      "51", "52" ,"53", "54",
      "61", "62" ,"63", "64",
      "71", "72" ,"73", "74",
      "81", "82" ,"83", "84",
      "91", "92" ,"93", "94",
      "101", "102" ,"103", "104",
      "111", "112" ,"113", "114",
      "121", "122" ,"123", "124",
      "131", "132" ,"133", "134",
      "141", "142" ,"143", "144",
      "151", "152" ,"153", "154" };
    
    if(fisDuplicate) {
      cout<<" Duplicated event!!!!   Recoil part will not be computed"<<endl;
      for (int iVal= 0; iVal < 64; iVal++) {
	//	cout<<Arstr[iVal]<<" "<<fValMap<<endl;
	if (strcmp(fValMap,Arstr[iVal]) == 0) {
	  fIVal = iVal + 10;
	}
      }
      fDupcount++;
    }


  }
  
  if(Ifile != " ") cout<<fDupcount<<" Num. of duplicates"<<endl;

}
  

// ----------------------------------------------------------------------
TFile* recoilBase::openHistFile(TString name) {
  fHistFile = new TFile(name.Data(), "RECREATE");
  fHistFile->cd();
  cout << "Opened " << fHistFile->GetName() << endl;
  return fHistFile;
}


// ----------------------------------------------------------------------
void recoilBase::closeHistFile() {
  cout << "Writing " << fHistFile->GetName() << endl;
  fHistFile->cd();
  fHistFile->Write();
  fHistFile->Close();
  delete fHistFile;

}

// ----------------------------------------------------------------------
void recoilBase::dumpEventList(const char *filename) {
  TFile *f = new TFile(filename, "RECREATE");
  fToBeCopied->Print();
  fChain->SetEventList(fToBeCopied);
  TTree *small = fChain->CopyTree("");
  small->Write();
  // small->Print();
}


// ----------------------------------------------------------------------
void recoilBase::bookHist(int dump) {
  TH1 *h;
  cout << " entering recoilBase bookHist " << endl;
  char name[100], title[100];
  char title1[100], variable[100];
  if (!fHistFile) {
    cout << "Call recoilBase::openHistFile(...) before booking histograms" << endl;
  } else {
    fHistFile->cd();
  }    

  if (dump > 0) {
    cout << "Booking events tree" << endl;
    fDump = dump;
    fTree = new TTree("events", "events"); 
    fTree->Branch("run",    &fRunnumber, "run/I");
    fTree->Branch("lower",  &fLower, "lower/I");
    fTree->Branch("upper",  &fUpper, "upper/I");

    // -- breco
    fTree->Branch("bmass",      &fBmass, "bmass/D");
    fTree->Branch("bmassfit",   &fBmassfit, "bmassfit/D");
    fTree->Branch("sbox",       &signalBox, "sbox/b");
    fTree->Branch("mes",        &fMes, "mes/D");
    fTree->Branch("de",         &fDeltaE, "de/D");
    fTree->Branch("pur",        &fPurity, "pur/D");
    fTree->Branch("intpur",     &fIntPurity, "intpur/D");
    fTree->Branch("brecoflav",  &fBrecoFlavor, "brecoflav/I");
    fTree->Branch("brecocharge",&fBrecoCharge , "brecocharge/I");
    fTree->Branch("brecomc",    &fBrecoMc,  "brecomc/I");

    // -- generator block
    fTree->Branch("mxhadgen", &fMxhadGen, "mxhadgen/D");
    fTree->Branch("pcmsgen",  &fPcmsGen, "pcmsgen/D");
    fTree->Branch("tcmsgen",  &fTcmsGen, "tcmsgen/D");
    fTree->Branch("fcmsgen",  &fFcmsGen, "fcmsgen/D");
    fTree->Branch("ecmsgen",  &fEcmsGen, "ecmsgen/D");
    fTree->Branch("pxhadgen", &fPxhadGen, "pxhadgen/D");
    fTree->Branch("txhadgen", &fTxhadGen, "txhadgen/D");
    fTree->Branch("fxhadgen", &fFxhadGen, "fxhadgen/D");
    fTree->Branch("exhadgen", &fExhadGen, "exhadgen/D");
    fTree->Branch("GoodEvent", &fGoodEvent, "GoodEvent/b");
    fTree->Branch("isDupli", &fisDuplicate, "isDupli/b");
    fTree->Branch("ValMap", &fIVal, "ValMap/I");

    fTree->Branch("vub",    &fVub, "vub/I");
    fTree->Branch("vcb",    &fVcb, "vcb/I");
    fTree->Branch("vxbtyp", &fBVxbTyp, "vxbtyp/I");
    fTree->Branch("other",  &fOther, "other/I");

    // -- recoil
    fTree->Branch("xcharge", &fRecoilCharge, "xcharge/I"); // note that this is the entire recoil, not xhad!
    fTree->Branch("pxhad",   &fPxhad, "pxhad/D");
    fTree->Branch("txhad",   &fTxhad, "txhad/D");
    fTree->Branch("fxhad",   &fFxhad, "fxhad/D");
    fTree->Branch("exhad",   &fExhad, "exhad/D");
    fTree->Branch("mxhad",   &fMxhad, "mxhad/D");
    fTree->Branch("gmax",    &fGammaMax,"gmax/D");
    fTree->Branch("mxhadfit",&fMxhadfit, "mxhadfit/D");

    // -- lepton
    fTree->Branch("lcharge", &fLeptonCharge ,"lcharge/I");
    fTree->Branch("plab",    &fPlab, "plab/D");
    fTree->Branch("tlab",    &fTlab, "tlab/D");
    fTree->Branch("flab",    &fFlab, "flab/D");
    fTree->Branch("pcms",    &fPcms, "pcms/D");
    fTree->Branch("tcms",    &fTcms, "tcms/D");
    fTree->Branch("fcms",    &fFcms, "fcms/D");
    fTree->Branch("ecms",    &fEcms, "ecms/D");

    fTree->Branch("nle",  &fNLepton, "nle/I");
    fTree->Branch("nel",  &fNEl, "nel/I");
    fTree->Branch("nmu",  &fNMu, "nmu/I");
    fTree->Branch("nchg",  &fRecoilTrkMult, "nchg/I");
    fTree->Branch("nneu",  &fRecoilNutMult, "nneu/I");
    fTree->Branch("nneu80_160",  &fRecoilNutMult80_160, "nneu80_160/I");
    fTree->Branch("nneu160_320",  &fRecoilNutMult160_320, "nneu160_320/I");
    fTree->Branch("nneufromB",  &fRecoilNutMultfromB, "nneufromB/I");
    fTree->Branch("nneufromB80_160",  &fRecoilNutMultfromB80_160, "nneufromB80_160/I");
    fTree->Branch("nneufromB160_320",  &fRecoilNutMultfromB160_320, "nneufromB160_320/I");
    fTree->Branch("nkp",  &fNKp, "nkp/I");
    fTree->Branch("nks",  &fNKshort, "nks/I");

    // -- neutrino
    fTree->Branch("pnu",    &fPNu, "pnu/D");
    fTree->Branch("tnu",    &fTNu, "tnu/D");
    fTree->Branch("fnu",    &fFNu, "fnu/D");
    fTree->Branch("mm2",    &fMM2, "mm2/D");
    fTree->Branch("mm2nc",  &fMM2NC, "mm2nc/D");
    fTree->Branch("mm2fit", &fMM2fit, "mm2fit/D");

    if (fOptGammas) {  
      fGTree = new TTree("photons", "photons"); 
      for(int ju =0; ju<15; ju++) {
	sprintf(title1, "G_NuM%d", ju);    sprintf(variable, "G_NuM%d%s", ju, "/D");
	fGTree->Branch(title1,    &tmpfMM2[ju], variable);
	sprintf(title1, "G_ResMx%d", ju);    sprintf(variable, "G_ResMx%d%s", ju, "/D");
	fGTree->Branch(title1,    &tmpxmassRes[ju], variable);
	sprintf(title1, "G_ResMxF%d", ju);    sprintf(variable, "G_ResMxF%d%s", ju, "/D");
	fGTree->Branch(title1,    &tmpxmassResF[ju], variable);
	sprintf(title1, "G_MxHad%d", ju);    sprintf(variable, "G_MxHad%d%s", ju, "/D");
	fGTree->Branch(title1,    &tmpfMxhad[ju], variable);
	sprintf(title1, "G_ThetaxHad%d", ju);    sprintf(variable, "G_ThetaxHad%d%s", ju, "/D");
	fGTree->Branch(title1,    &tmpfTxhad[ju], variable);
	sprintf(title1, "G_PhixHad%d", ju);    sprintf(variable, "G_PhixHad%d%s", ju, "/D");
	fGTree->Branch(title1,    &tmpfFxhad[ju], variable);
	sprintf(title1, "G_EnxHad%d", ju);    sprintf(variable, "G_EnxHad%d%s", ju, "/D");
	fGTree->Branch(title1,    &tmpfExhad[ju], variable);
	sprintf(title1, "G_MxHadF%d", ju);    sprintf(variable, "G_MxHadF%d%s", ju, "/D");
	fGTree->Branch(title1,    &tmpfMxhadfit[ju], variable);
	sprintf(title1, "G_ThetaxHadF%d", ju);    sprintf(variable, "G_ThetaxHadF%d%s", ju, "/D");
	fGTree->Branch(title1,    &tmpfTxhadfit[ju], variable);
	sprintf(title1, "G_PhixHadF%d", ju);    sprintf(variable, "G_PhixHadF%d%s", ju, "/D");
	fGTree->Branch(title1,    &tmpfFxhadfit[ju], variable);
	sprintf(title1, "G_EnxHadF%d", ju);    sprintf(variable, "G_EnxHadF%d%s", ju, "/D");
	fGTree->Branch(title1,    &tmpfExhadfit[ju], variable);
	sprintf(title1, "G_CountNeu%d", ju);    sprintf(variable, "G_CountNeu%d%s", ju, "/I");
	fGTree->Branch(title1,    &fcountNeu[ju], variable);
	sprintf(title1, "G_MassNeu%d", ju);    sprintf(variable, "G_MassNeu%d%s", ju, "/D");
	fGTree->Branch(title1,    &fMNeupart[ju], variable);
	sprintf(title1, "G_EnNeu%d", ju);    sprintf(variable, "G_EnNeu%d%s", ju, "/D");
	fGTree->Branch(title1,    &fENeupart[ju], variable);
	sprintf(title1, "G_ThNeu%d", ju);    sprintf(variable, "G_ThNeu%d%s", ju, "/D");
	fGTree->Branch(title1,    &fTNeupart[ju], variable);
	sprintf(title1, "G_PhiNeu%d", ju);    sprintf(variable, "G_PhiNeu%d%s", ju, "/D");
	fGTree->Branch(title1,    &fPNeupart[ju], variable);
	sprintf(title1, "GoodEvent%d", ju);    sprintf(variable, "GoodEvent%d%s", ju, "/b");
	fGTree->Branch(title1,    &fGoodEventPh[ju], variable);
	
      }
      
      for(int juu =0; juu<10; juu++) {
	sprintf(title1, "GoodNoHol%d", juu);    sprintf(variable, "GoodNoHol%d%s", juu, "/b");
	fGTree->Branch(title1, &fGoodNoHole[juu], variable);
      }
      
      fGTree->Branch("Probchi2fit", &fProbChi2, "Probchi2fit/D");
      fGTree->Branch("Mm2fit", &fMM2fit, "Mm2fit/D");
      fGTree->Branch("Mxhadfit",&fMxhadfit, "Mxhadfit/D");
      fGTree->Branch("Bmassfit",   &fBmassfit, "Bmassfit/D");
      fGTree->Branch("mes",        &fMes, "mes/D");
      fGTree->Branch("MxHadGen", &fMxhadGen, "MxHadGen/D");
      fGTree->Branch("NumofChgpart", &fcountChg, "NumofChgpart/I");
      fGTree->Branch("GoodEvent", &fGoodEvent, "GoodEvent/b");
      fGTree->Branch("isDupli", &fisDuplicate, "isDupli/b");
      fGTree->Branch("ValMap", &fIVal, "ValMap/I");
      fGTree->Branch("ChargeCorr", &fGoodChargeCorr, "ChargeCorr/b");
      fGTree->Branch("ChargeCons", &fGoodChargeCons, "ChargeCons/b");
      fGTree->Branch("GoodLep", &fGoodLepton, "GoodLep/b");
      fGTree->Branch("Masschpa", &fMCharpart, "Masschpa/D");
      fGTree->Branch("Enchpa", &fECharpart, "Enchpa/D");
      fGTree->Branch("Thchpa", &fTCharpart, "Thchpa/D");
      fGTree->Branch("Phichpa", &fPCharpart, "Phichpa/D");
      fGTree->Branch("sbox",    &signalBox,   "sbox/b");
      fGTree->Branch("sideb",   &mesSideband, "sideb/b");
      //    fGTree->Branch("ChgMom", fChargPart, "ChgMom[4]/D");

    }

  }

  fHistFile->mkdir("mcTruth", "mcTruth");
  fHistFile->cd("mcTruth");
  sprintf(name, "h100");  sprintf(title, "nVub");  h = new TH1D(name, title, 3, 0., 3.); 
  sprintf(name, "h101");  sprintf(title, "nVcb");  h = new TH1D(name, title, 3, 0., 3.); 
  sprintf(name, "h102");  sprintf(title, "n fully reco Dstar");  h = new TH1D(name, title, 3, 0., 3.); 
  sprintf(name, "h103");  sprintf(title, "n fully reco Dc");  h = new TH1D(name, title, 3, 0., 3.); 
  sprintf(name, "h104");  sprintf(title, "n fully reco Dstar0");  h = new TH1D(name, title, 3, 0., 3.); 
  sprintf(name, "h105");  sprintf(title, "n fully reco D0");  h = new TH1D(name, title, 3, 0., 3.); 
  sprintf(name, "h106");  sprintf(title, "nOther");  h = new TH1D(name, title, 3, 0., 3.); 
  sprintf(name, "h107");  sprintf(title, "fVxbTyp");  h = new TH1D(name, title, 11, -1., 10.); 
  sprintf(name, "h700");  sprintf(title, "id mc");  h = new TH1D(name, title, 100, -600., 600.); 
  sprintf(name, "h9009");  sprintf(title, "type event MC");  h = new TH2D(name, title, 7, 1., 8., 7, 1., 8.); 
  
  sprintf(name, "h77000");  sprintf(title, "number of leptons");  h = new TH1D(name, title, 4, 0., 4.); 
  
  sprintf(name, "h121000");  sprintf(title, "generator pcms e/mu");  h = new TH1D(name, title, PCBIN, 0., PCMAX); 
  sprintf(name, "h121005");  sprintf(title, "generator pcms e/mu Vcb"); h = new TH1D(name, title, PCBIN, 0., PCMAX); 
  sprintf(name, "h121006");  sprintf(title, "generator pcms e/mu Vub"); h = new TH1D(name, title, PCBIN, 0., PCMAX); 

  sprintf(name, "h123000");  sprintf(title, "generator Mx one lepton");  h = new TH1D(name, title, XBIN, 0., XMAX); 
  sprintf(name, "h123005");  sprintf(title, "generator Mx one lepton Vcb"); h = new TH1D(name, title, XBIN, 0., XMAX); 
  sprintf(name, "h123006");  sprintf(title, "generator Mx one lepton Vub"); h = new TH1D(name, title, XBIN, 0., XMAX); 

  sprintf(name, "h124000");  sprintf(title, "generator N hadrons with one lepton"); h = new TH1D(name, title, 20, 0., 20.); 
  sprintf(name, "h124005");  sprintf(title, "generator N hadrons with one lepton Vcb"); h = new TH1D(name, title, 20, 0., 20.); 
  sprintf(name, "h124006");  sprintf(title, "generator N hadrons with one lepton Vub"); h = new TH1D(name, title, 20, 0., 20.); 
  
  if (fOptGammas) {
    fHistFile->cd();
    cout<<"making the directory"<<endl;
    fHistFile->mkdir("Photon", "Photon");     fHistFile->cd("Photon");
    
    //Photon block
    sprintf(name, "ThetaGam");  sprintf(title, "Cos(Theta) (E<80 MeV)");  h = new TH1D(name, title, 50, -1., 1.0); h->Sumw2();
    sprintf(name, "ThetaGam Un");  sprintf(title, "Cos(Theta) (E<80 MeV)");  h = new TH1D(name, title, 50, -1., 1.0); h->Sumw2();
    sprintf(name, "EGam Un");  sprintf(title, "E(gam) MeV");  h = new TH1D(name, title, 100, 0., 1.0); h->Sumw2();
    sprintf(name, "EGam Mt");  sprintf(title, "E(gam) MeV");  h = new TH1D(name, title, 100, 0., 1.0); h->Sumw2();
    sprintf(name, "EGam Un vub");  sprintf(title, "E(gam) MeV");  h = new TH1D(name, title, 100, 0., 1.0); h->Sumw2();
    sprintf(name, "EGam Mt vub");  sprintf(title, "E(gam) MeV");  h = new TH1D(name, title, 100, 0., 1.0); h->Sumw2();
    sprintf(name, "EGam Un vcb");  sprintf(title, "E(gam) MeV");  h = new TH1D(name, title, 100, 0., 1.0); h->Sumw2();
    sprintf(name, "EGam Mt vcb");  sprintf(title, "E(gam) MeV");  h = new TH1D(name, title, 100, 0., 1.0); h->Sumw2();
    sprintf(name, "ThetaGamcut");  sprintf(title, "Cos(Theta) (E>80 MeV)");  h = new TH1D(name, title, 50, -1., 1.0); h->Sumw2();
    sprintf(name, "PhiGam");  sprintf(title, "Phi of the recoil gamma");  h = new TH1D(name, title, 100, -3.15, 3.15); h->Sumw2();
    sprintf(name, "PhiGam Un");  sprintf(title, "Phi of the recoil gamma");  h = new TH1D(name, title, 100, -3.15, 3.15); h->Sumw2();
    sprintf(name, "CrySoft");  sprintf(title, "# of crystal all gammas");  h = new TH1D(name, title, 30, -5., 25.); h->Sumw2(); h->SetMarkerStyle(3); h->SetMarkerSize(2);
    sprintf(name, "CryHard");  sprintf(title, "# of crystal all gammas");  h = new TH1D(name, title, 30, -5., 25.); h->Sumw2();h->SetMarkerStyle(3); h->SetMarkerSize(2);
    sprintf(name, "CrySoftpi");  sprintf(title, "# of crystal for pi0's gammas");  h = new TH1D(name, title, 30, -5., 25.); h->Sumw2();h->SetMarkerStyle(3); h->SetMarkerSize(2);
    sprintf(name, "CryHardpi");  sprintf(title, "# of crystal for pi0's gammas");  h = new TH1D(name, title, 30, -5., 25.); h->Sumw2();h->SetMarkerStyle(3); h->SetMarkerSize(2);
    sprintf(name, "EnGam");  sprintf(title, "Energy");  h = new TH1D(name, title, 100, 0., 1.); h->Sumw2();
    sprintf(name, "EnGamKs");  sprintf(title, "Energy");  h = new TH1D(name, title, 100, 0., 1.); h->Sumw2();
    sprintf(name, "EnGamNoKs");  sprintf(title, "Energy");  h = new TH1D(name, title, 100, 0., 1.); h->Sumw2();
    sprintf(name, "lMomGam");  sprintf(title, "l mom ");  h = new TH1D(name, title, 100, 0., 0.4); h->Sumw2();
    sprintf(name, "ZMom42Gam");  sprintf(title, "Z42 mom ");  h = new TH1D(name, title, 80, 0., 0.2); h->Sumw2();
    sprintf(name, "Theta vs E");  sprintf(title, "Theta vs E");  h = new TH2D(name, title, 100, -1., 1., 100, 0., 3.);
    sprintf(name, "E vs phi");  sprintf(title, "E vs phi");  h = new TH2D(name, title, 80, 0., 2., 80, -3.15, 3.15);
    sprintf(name, "E vs phi cut");  sprintf(title, "E vs phi cut");  h = new TH2D(name, title, 80, 0., 0.1, 80, -3.15, 3.15);
    sprintf(name, "E vs th");  sprintf(title, "E vs theta");  h = new TH2D(name, title, 80, -1. , 1. , 80, 0.0, 0.08);

    sprintf(name, "EnGam (MA)");  sprintf(title, "Energy");  h = new TH1D(name, title, 80, 0.08, 1.); h->Sumw2();
    sprintf(name, "EnGamNMP");  sprintf(title, "Energy");  h = new TH1D(name, title, 100, 0., 1.); h->Sumw2();
    sprintf(name, "EnGamMP");  sprintf(title, "Energy");  h = new TH1D(name, title, 100, 0., 1.); h->Sumw2();
    sprintf(name, "EnGamMA");  sprintf(title, "Energy");  h = new TH1D(name, title, 100, 0., 1.); h->Sumw2();
    sprintf(name, "EnGamMB vub");  sprintf(title, "Energy");  h = new TH1D(name, title, 100, 0., 1.); h->Sumw2();
    sprintf(name, "EnGamMB vcb");  sprintf(title, "Energy");  h = new TH1D(name, title, 100, 0., 1.); h->Sumw2();
    sprintf(name, "EnGamMB vub cut");  sprintf(title, "Energy");  h = new TH1D(name, title, 100, 0.07, 1.); h->Sumw2();
    sprintf(name, "EnGamMB vcb cut");  sprintf(title, "Energy");  h = new TH1D(name, title, 100, 0.07, 1.); h->Sumw2();
    sprintf(name, "lMomGam (MA)");  sprintf(title, "l mom ");  h = new TH1D(name, title, 100, 0., 0.4); h->Sumw2();
    sprintf(name, "ZMom42Gam (MA)");  sprintf(title, "Z42 mom ");  h = new TH1D(name, title, 100, 0., 0.2); h->Sumw2();
    sprintf(name, "PhiGam (MA)");  sprintf(title, "Phi of the recoil gamma");  h = new TH1D(name, title, 100, -3.15, 3.15); h->Sumw2();
    
    sprintf(name, "EnGam (MA) cut");  sprintf(title, "Energy");  h = new TH1D(name, title, 100, 0., 0.08); h->Sumw2();
    sprintf(name, "lMomGam (MA) cut");  sprintf(title, "l mom ");  h = new TH1D(name, title, 100, 0., 0.4); h->Sumw2();
    sprintf(name, "ZMom42Gam (MA) cut");  sprintf(title, "Z42 mom ");  h = new TH1D(name, title, 100, 0., 0.2); h->Sumw2();
    sprintf(name, "PhiGam (MA) cut");  sprintf(title, "Phi of the recoil gamma");  h = new TH1D(name, title, 100, -3.15, 3.15); h->Sumw2();
    sprintf(name, "E vs theta (MA)");  sprintf(title, "Theta vs E");  h = new TH2D(name, title, 100, -1., 1., 100, 0., 2.);
    
    sprintf(name, "EnGam (MB)");  sprintf(title, "Energy");  h = new TH1D(name, title, 100, 0.08, 1.); h->Sumw2();
    sprintf(name, "lMomGam (MB)");  sprintf(title, "l mom ");  h = new TH1D(name, title, 100, 0., 0.4); h->Sumw2();
    sprintf(name, "ZMom42Gam (MB)");  sprintf(title, "Z42 mom ");  h = new TH1D(name, title, 100, 0., 0.2); h->Sumw2();
    sprintf(name, "PhiGam (MB)");  sprintf(title, "Phi of the recoil gamma");  h = new TH1D(name, title, 100, -3.15, 3.15); h->Sumw2();

    sprintf(name, "EnGam (MB) cut");  sprintf(title, "Energy");  h = new TH1D(name, title, 100, 0., 0.08); h->Sumw2();
    sprintf(name, "lMomGam (MB) cut");  sprintf(title, "l mom ");  h = new TH1D(name, title, 100, 0., 0.4); h->Sumw2();
    sprintf(name, "ZMom42Gam (MB) cut");  sprintf(title, "Z42 mom ");  h = new TH1D(name, title, 100, 0., 0.2); h->Sumw2();
    sprintf(name, "PhiGam (MB) cut");  sprintf(title, "Phi of the recoil gamma");  h = new TH1D(name, title, 100, -3.15, 3.15); h->Sumw2();
    sprintf(name, "E vs theta (MB)");  sprintf(title, "Theta vs E");  h = new TH2D(name, title, 100, -1., 1., 100, 0., 2.);

    sprintf(name, "EnGam (NMP)");  sprintf(title, "Energy");  h = new TH1D(name, title, 100, 0.08, 1.); h->Sumw2();
    sprintf(name, "lMomGam (NMP)");  sprintf(title, "l mom ");  h = new TH1D(name, title, 100, 0., 0.4); h->Sumw2();
    sprintf(name, "ZMom42Gam (NMP)");  sprintf(title, "Z42 mom ");  h = new TH1D(name, title, 100, 0., 0.2); h->Sumw2();
    sprintf(name, "PhiGam (NMP)");  sprintf(title, "Phi of the recoil gamma");  h = new TH1D(name, title, 100, -3.15, 3.15); h->Sumw2();
    
    sprintf(name, "EnGam (NMP) cut");  sprintf(title, "Energy");  h = new TH1D(name, title, 100, 0., 0.08); h->Sumw2();
    sprintf(name, "lMomGam (NMP) cut");  sprintf(title, "l mom ");  h = new TH1D(name, title, 100, 0., 0.4); h->Sumw2();
    sprintf(name, "ZMom42Gam (NMP) cut");  sprintf(title, "Z42 mom ");  h = new TH1D(name, title, 100, 0., 0.2); h->Sumw2();
    sprintf(name, "PhiGam (NMP) cut");  sprintf(title, "Phi of the recoil gamma");  h = new TH1D(name, title, 100, -3.15, 3.15); h->Sumw2();
    sprintf(name, "E vs theta (NMP)");  sprintf(title, "Theta vs E");  h = new TH2D(name, title, 100, -1., 1., 100, 0., 2.);
    
    fHistFile->cd();
        
    sprintf(name, "EREcoil");  sprintf(title, "Energy recoil");  h = new TH1D(name, title, 100, 2., 8.0); h->Sumw2();
    sprintf(name, "PREcoil");  sprintf(title, "Momentum recoil");  h = new TH1D(name, title, 100, 0., 5.0); h->Sumw2();
    sprintf(name, "FREcoil");  sprintf(title, "Phi recoil");  h = new TH1D(name, title, 100, -3.15, 3.15); h->Sumw2();
    sprintf(name, "TREcoil");  sprintf(title, "Theta recoil");  h = new TH1D(name, title, 100, 0., 2.); h->Sumw2();
    sprintf(name, "MREcoil");  sprintf(title, "Mass recoil");  h = new TH1D(name, title, 100, 3., 7.5); h->Sumw2();
    
    sprintf(name, "EAllev");  sprintf(title, "Energy all event");  h = new TH1D(name, title, 100, 8., 16.0); h->Sumw2();
    sprintf(name, "PAllev");  sprintf(title, "Momentum all event");  h = new TH1D(name, title, 100, 2., 10.0); h->Sumw2();
    sprintf(name, "FAllev");  sprintf(title, "Phi all event");  h = new TH1D(name, title, 100, -3.15, 3.15); h->Sumw2();
    sprintf(name, "TAllev");  sprintf(title, "Theta all event");  h = new TH1D(name, title, 100, 0., 2.); h->Sumw2();
    sprintf(name, "MAllev");  sprintf(title, "Mass all event");  h = new TH1D(name, title, 100, 7.5, 12.5); h->Sumw2();
    
    sprintf(name, "EREcoilsig");  sprintf(title, "Energy recoil");  h = new TH1D(name, title, 100, 2., 8.0); h->Sumw2();
    sprintf(name, "PREcoilsig");  sprintf(title, "Momentum recoil");  h = new TH1D(name, title, 100, 0., 5.0); h->Sumw2();
    sprintf(name, "FREcoilsig");  sprintf(title, "Phi recoil");  h = new TH1D(name, title, 100, -3.15, 3.15); h->Sumw2();
    sprintf(name, "TREcoilsig");  sprintf(title, "Theta recoil");  h = new TH1D(name, title, 100, 0., 2.); h->Sumw2();
    sprintf(name, "MREcoilsig");  sprintf(title, "Mass recoil");  h = new TH1D(name, title, 100, 3., 7.5); h->Sumw2();
    
    sprintf(name, "EAllevsig");  sprintf(title, "Energy all event");  h = new TH1D(name, title, 100, 8., 16.0); h->Sumw2();
    sprintf(name, "PAllevsig");  sprintf(title, "Momentum all event");  h = new TH1D(name, title, 100, 2., 10.0); h->Sumw2();
    sprintf(name, "FAllevsig");  sprintf(title, "Phi all event");  h = new TH1D(name, title, 100, -3.15, 3.15); h->Sumw2();
    sprintf(name, "TAllevsig");  sprintf(title, "Theta all event");  h = new TH1D(name, title, 100, 0., 2.); h->Sumw2();
    sprintf(name, "MAllevsig");  sprintf(title, "Mass all event");  h = new TH1D(name, title, 100, 7.5, 12.5); h->Sumw2();
  }
  
    
  fHistFile->cd();
  
  sprintf(name, "mesallevents");  sprintf(title, "mes All");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "deallevents");  sprintf(title, "delta E All");  h = new TH1D(name, title, 40, -0.1, 0.1); 
  sprintf(name, "mesalleventsA");  sprintf(title, "mes Vub enhanced");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesalleventsB");  sprintf(title, "mes Vub enhanced");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesalleventsC");  sprintf(title, "mes Vub enhanced");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  
  sprintf(name, "mesallcuts");  sprintf(title, "mes All Cuts");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "deallcuts");  sprintf(title, "delta E All Cuts");  h = new TH1D(name, title, 40, -0.1, 0.1); 
  sprintf(name, "mesallcutsA");  sprintf(title, "mes Vub enhanced");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesallcutsB");  sprintf(title, "mes Vub enhanced");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesallcutsC");  sprintf(title, "mes Vub enhanced");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  
  sprintf(name, "mesallel");  sprintf(title, "mes All el");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "deallel");  sprintf(title, "delta E All el");  h = new TH1D(name, title, 40, -0.1, 0.1); 
  sprintf(name, "mesallelA");  sprintf(title, "mes Vub enhanced");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesallelB");  sprintf(title, "mes Vub enhanced");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesallelC");  sprintf(title, "mes Vub enhanced");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  
  sprintf(name, "mesallmu");  sprintf(title, "mes All mu");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "deallmu");  sprintf(title, "delta E All mu");  h = new TH1D(name, title, 40, -0.1, 0.1); 
  sprintf(name, "mesallmuA");  sprintf(title, "mes Vub enhanced");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesallmuB");  sprintf(title, "mes Vub enhanced");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesallmuC");  sprintf(title, "mes Vub enhanced");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  
  sprintf(name, "mesrecoil");  sprintf(title, "mes Recoil");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "derecoil");  sprintf(title, "delta E Recoil");  h = new TH1D(name, title, 40, -0.1, 0.1); 
  sprintf(name, "mesrecoilA");  sprintf(title, "mes Vub enhanced");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesrecoilB");  sprintf(title, "mes Vub enhanced");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesrecoilC");  sprintf(title, "mes Vub enhanced");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  
  sprintf(name, "mesall");  sprintf(title, "mes All");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "deall");  sprintf(title, "delta E All");  h = new TH1D(name, title, 40, -0.1, 0.1); 
  sprintf(name, "mesallA");  sprintf(title, "mes Vub enhanced");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesallB");  sprintf(title, "mes Vub enhanced");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesallC");  sprintf(title, "mes Vub enhanced");  h = new TH1D(name, title, 40, 5.2, 5.3); 
    
  sprintf(name, "mesvcb");  sprintf(title, "mes Vub depleted");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "devcb");  sprintf(title, "delta Vub depleted");  h = new TH1D(name, title, 40, -0.1, 0.1); 
  sprintf(name, "mesvcbA");  sprintf(title, "mes Vub enhanced");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesvcbB");  sprintf(title, "mes Vub enhanced");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesvcbC");  sprintf(title, "mes Vub enhanced");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  
  sprintf(name, "mesvub");  sprintf(title, "mes Vub enhanced");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "devub");  sprintf(title, "delta Vub enhanced");  h = new TH1D(name, title, 40, -0.1, 0.1); 
  sprintf(name, "mesvubA");  sprintf(title, "mes Vub enhanced");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesvubB");  sprintf(title, "mes Vub enhanced");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesvubC");  sprintf(title, "mes Vub enhanced");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  
  
  sprintf(name,"p100");  sprintf(title, "prob chi2");  h = new TH1D(name, title, 100, 0., 1.); h->Sumw2();
  sprintf(name,"p101");  sprintf(title, "chi2");       h = new TH1D(name, title, 100, 0., 100.); h->Sumw2();
  
  int no(0), le(0), i(0);
  char lc[10];
  for (i = 0; i < 9; ++i) {
    if (i==0) {fHistFile->mkdir("recoil","recoil");  fHistFile->cd("recoil"); gDirectory->mkdir("CutPlots","CutPlots"); }
    if (i==1) {fHistFile->mkdir("bgall","Sideband all");  fHistFile->cd("bgall"); gDirectory->mkdir("CutPlots","CutPlots"); }
    if (i==2) {fHistFile->mkdir("sgall","Signalbox all");  fHistFile->cd("sgall"); gDirectory->mkdir("CutPlots","CutPlots"); }
    if (i==3) {fHistFile->mkdir("bgvcb","Sideband Vub depleted");fHistFile->cd("bgvcb"); gDirectory->mkdir("CutPlots","CutPlots"); }
    if (i==4) {fHistFile->mkdir("sgvcb","Signalbox Vub depleted");fHistFile->cd("sgvcb"); gDirectory->mkdir("CutPlots","CutPlots"); }
    if (i==5) {fHistFile->mkdir("bgvub","Sideband Vub enhanced");fHistFile->cd("bgvub"); gDirectory->mkdir("CutPlots","CutPlots");}
    if (i==6) {fHistFile->mkdir("sgvub","Signalbox Vub enhanced");fHistFile->cd("sgvub"); gDirectory->mkdir("CutPlots","CutPlots");}
    if (i==7) {fHistFile->mkdir("bgallevents","Sideband");fHistFile->cd("bgallevents");gDirectory->mkdir("CutPlots","CutPlots"); }
    if (i==8) {fHistFile->mkdir("sgallevents","Signalbox");fHistFile->cd("sgallevents"); gDirectory->mkdir("CutPlots","CutPlots");}
    
    
    sprintf(name,"ks100");  sprintf(title, "mass kshorts");  h = new TH1D(name, title, 100, 0.46, 0.54); h->Sumw2();
    sprintf(name,"ks101");  sprintf(title, "mass kshorts goodLepton");  h = new TH1D(name, title, 100, 0.46, 0.54); h->Sumw2();
    sprintf(name,"ks102");  sprintf(title, "mass kshorts goodMM2");  h = new TH1D(name, title, 100, 0.46, 0.54); h->Sumw2();
    sprintf(name,"ks103");  sprintf(title, "mass kshorts chargeCorr");  h = new TH1D(name, title, 100, 0.46, 0.54); h->Sumw2();
    sprintf(name,"ks104");  sprintf(title, "mass kshorts chargeCons");  h = new TH1D(name, title, 100, 0.46, 0.54); h->Sumw2();
    sprintf(name,"ks105");  sprintf(title, "mass kshorts goodEvent");  h = new TH1D(name, title, 100, 0.46, 0.54); h->Sumw2();
  
    sprintf(name,"ks200");  sprintf(title, "plab kshorts");  h = new TH1D(name, title, 50, 0., 5.); h->Sumw2();
    sprintf(name,"ks201");  sprintf(title, "tlab kshorts");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();

    sprintf(name,"kp200");  sprintf(title, "plab k");  h = new TH1D(name, title, 50, 0., 5.); h->Sumw2();
    sprintf(name,"kp201");  sprintf(title, "tlab k");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
    
    sprintf(name, "s100");  sprintf(title, "gen Breco momentum");  h = new TH1D(name, title, 100, 2.5, 3.5); 
    sprintf(name, "s101");  sprintf(title, "reco Breco momentum");  h = new TH1D(name, title, 100, 2.5, 3.5); 
    sprintf(name, "s102");  sprintf(title, "smeared Breco momentum");  h = new TH1D(name, title, 100, 2.5, 3.5); 
    sprintf(name, "s103");  sprintf(title, "reco-smeared Breco momentum");  h = new TH1D(name, title, 100, -0.5, 0.5); 
    sprintf(name, "s104");  sprintf(title, "reco-generated Breco momentum");  h = new TH1D(name, title, 100, -0.5, 0.5); 
    sprintf(name, "s105");  sprintf(title, "generated-smeared Breco momentum");  h = new TH1D(name, title, 100, -0.5, 0.5); 
      
    
    
    if (i == 0) {
      sprintf(name,"elSelBits");  sprintf(title, "electron selector bits"); h = new TH1D(name, title, 33, -1., 32.);
      sprintf(name,"muSelBits");  sprintf(title, "muon     selector bits"); h = new TH1D(name, title, 33, -1., 32.);
      sprintf(name,"kaSelBits");  sprintf(title, "kaon     selector bits"); h = new TH1D(name, title, 33, -1., 32.);
    }      
    
    for (le = 0; le < 3; ++le) {
      if (le == 0) sprintf(lc, "a");
      if (le == 1) sprintf(lc, "e");
      if (le == 2) sprintf(lc, "m");
      
      // -- mes histograms for cut studies
      no = 1; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      no = 2; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      no = 3; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      no = 4; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      no = 5; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      no = 6; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      no = 7; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      
      no =11; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      no =12; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      no =13; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      no =14; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      no =15; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      no =16; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      no =17; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      
      // -- Lepton
      no = 1000;
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "pcms all e/mu");      h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "pcms all e/mu Vcb");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "pcms all e/mu Vub");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "pcms all e/mu D    ");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "pcms all e/mu D*   ");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "pcms all e/mu D(*)x");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "pcms all e/mu other");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      
      no = 1010;		     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "pcms all e/mu");      h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "pcms all e/mu Vcb");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "pcms all e/mu Vub");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "pcms all e/mu D");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "pcms all e/mu D*");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "pcms all e/mu D(*)x");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "pcms all e/mu other");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      
      no = 1020;		     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "pcms all e/mu");      h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "pcms all e/mu Vcb");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "pcms all e/mu Vub");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "pcms all e/mu D");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "pcms all e/mu D*");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "pcms all e/mu D(*)x");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "pcms all e/mu other");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      
      no = 1100;         
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "pcms prompt e/mu");      h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "pcms prompt e/mu Vcb");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "pcms prompt e/mu Vub");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "pcms prompt e/mu D    ");h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "pcms prompt e/mu D*   ");h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "pcms prompt e/mu D(*)x");h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "pcms prompt e/mu other");h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();

      no = 1110;         
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "pcms prompt e/mu");      h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "pcms prompt e/mu Vcb");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "pcms prompt e/mu Vub");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "pcms prompt e/mu D    ");h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "pcms prompt e/mu D*   ");h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "pcms prompt e/mu D(*)x");h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "pcms prompt e/mu other");h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();

      no = 1120;         
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "pcms prompt e/mu");      h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "pcms prompt e/mu Vcb");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "pcms prompt e/mu Vub");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "pcms prompt e/mu D    ");h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "pcms prompt e/mu D*   ");h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "pcms prompt e/mu D(*)x");h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "pcms prompt e/mu other");h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();

      no = 1200;           
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "pcms cascad e/mu");      h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "pcms cascad e/mu Vcb");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "pcms cascad e/mu Vub");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "pcms cascad e/mu D    ");h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "pcms cascad e/mu D*   ");h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "pcms cascad e/mu D(*)x");h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "pcms cascad e/mu other");h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();

      no = 1210;           
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "pcms cascad e/mu");      h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "pcms cascad e/mu Vcb");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "pcms cascad e/mu Vub");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "pcms cascad e/mu D    ");h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "pcms cascad e/mu D*   ");h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "pcms cascad e/mu D(*)x");h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "pcms cascad e/mu other");h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();

      no = 1220;           
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "pcms cascad e/mu");      h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "pcms cascad e/mu Vcb");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "pcms cascad e/mu Vub");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "pcms cascad e/mu D    ");h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "pcms cascad e/mu D*   ");h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "pcms cascad e/mu D(*)x");h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "pcms cascad e/mu other");h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();

      no = 1300;           
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "cos#theta e/mu");      h = new TH1D(name, title, 50, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "cos#theta e/mu Vcb");  h = new TH1D(name, title, 50, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "cos#theta e/mu Vub");  h = new TH1D(name, title, 50, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "cos#theta e/mu D    ");  h = new TH1D(name, title, 50, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "cos#theta e/mu D*   ");  h = new TH1D(name, title, 50, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "cos#theta e/mu D(*)x");  h = new TH1D(name, title, 50, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "cos#theta e/mu other");  h = new TH1D(name, title, 50, -1., 1.); h->Sumw2();

      no = 1310;           
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "cos#theta e/mu");      h = new TH1D(name, title, 50, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "cos#theta e/mu Vcb");  h = new TH1D(name, title, 50, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "cos#theta e/mu Vub");  h = new TH1D(name, title, 50, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "cos#theta e/mu D    ");  h = new TH1D(name, title, 50, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "cos#theta e/mu D*   ");  h = new TH1D(name, title, 50, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "cos#theta e/mu D(*)x");  h = new TH1D(name, title, 50, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "cos#theta e/mu other");  h = new TH1D(name, title, 50, -1., 1.); h->Sumw2();

      no = 1320;           
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "cos#theta e/mu");      h = new TH1D(name, title, 50, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "cos#theta e/mu Vcb");  h = new TH1D(name, title, 50, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "cos#theta e/mu Vub");  h = new TH1D(name, title, 50, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "cos#theta e/mu D    ");  h = new TH1D(name, title, 50, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "cos#theta e/mu D*   ");  h = new TH1D(name, title, 50, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "cos#theta e/mu D(*)x");  h = new TH1D(name, title, 50, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "cos#theta e/mu other");  h = new TH1D(name, title, 50, -1., 1.); h->Sumw2();
				     
      no = 1400;		     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "pups all e/mu");      h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "pups all e/mu Vcb");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "pups all e/mu Vub");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "pups all e/mu D    ");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "pups all e/mu D*   ");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "pups all e/mu D(*)x");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "pups all e/mu other");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
				     
      no = 1410;		     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "pups all e/mu");      h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "pups all e/mu Vcb");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "pups all e/mu Vub");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "pups all e/mu D");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "pups all e/mu D*");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "pups all e/mu D(*)x");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "pups all e/mu other");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();

      no = 1420;		     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "pups all e/mu");      h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "pups all e/mu Vcb");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "pups all e/mu Vub");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "pups all e/mu D");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "pups all e/mu D*");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "pups all e/mu D(*)x");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "pups all e/mu other");  h = new TH1D(name, title, PCBIN, 0., PCMAX); h->Sumw2();
				     
				     
      // -- Recoil		     
      no = 2000;		     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "xhad mass");     h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "xhad mass Vcb"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "xhad mass Vub"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "xhad mass D    "); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "xhad mass D*   "); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "xhad mass D(*)x"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "xhad mass other"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
				     
      no = 2010;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "xhad mass");     h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "xhad mass Vcb"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "xhad mass Vub"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "xhad mass D    "); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "xhad mass D*   "); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "xhad mass D(*)x"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "xhad mass other"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();

      no = 2020;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "xhad mass");     h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "xhad mass Vcb"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "xhad mass Vub"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "xhad mass D    "); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "xhad mass D*   "); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "xhad mass D(*)x"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "xhad mass other"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
				     
      no = 2100;		     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "xhad mass");     h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "xhad mass Vcb"); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "xhad mass Vub"); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "xhad mass D    "); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "xhad mass D*   "); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "xhad mass D(*)x"); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "xhad mass other"); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      				     
      no = 2110;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "xhad mass");     h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "xhad mass Vcb"); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "xhad mass Vub"); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "xhad mass D    "); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "xhad mass D*   "); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "xhad mass D(*)x"); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "xhad mass other"); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();

      no = 2120;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "xhad mass");     h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "xhad mass Vcb"); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "xhad mass Vub"); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "xhad mass D    "); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "xhad mass D*   "); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "xhad mass D(*)x"); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "xhad mass other"); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
				     
      no = 2200;		     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "xhadgen mass");     h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "xhadgen mass Vcb"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "xhadgen mass Vub"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "xhadgen mass D    "); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "xhadgen mass D*   "); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "xhadgen mass D(*)x"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "xhadgen mass other"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      				     
      no = 2210;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "xhadgen mass");     h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "xhadgen mass Vcb"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "xhadgen mass Vub"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "xhadgen mass D    "); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "xhadgen mass D*   "); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "xhadgen mass D(*)x"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "xhadgen mass other"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();

      no = 2220;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "xhadgen mass");     h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "xhadgen mass Vcb"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "xhadgen mass Vub"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "xhadgen mass D    "); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "xhadgen mass D*   "); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "xhadgen mass D(*)x"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "xhadgen mass other"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
				     
      no = 2300;		     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "xhad mass**2");     h = new TH1D(name, title, XBIN, -4., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "xhad mass**2 Vcb"); h = new TH1D(name, title, XBIN, -4., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "xhad mass**2 Vub"); h = new TH1D(name, title, XBIN, -4., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "xhad mass**2 D    "); h = new TH1D(name, title, XBIN, -4., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "xhad mass**2 D*   "); h = new TH1D(name, title, XBIN, -4., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "xhad mass**2 D(*)x"); h = new TH1D(name, title, XBIN, -4., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "xhad mass**2 other"); h = new TH1D(name, title, XBIN, -4., 10.); h->Sumw2();
      				     
      no = 2310;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "xhad mass");     h = new TH1D(name, title, XBIN, -4., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "xhad mass Vcb"); h = new TH1D(name, title, XBIN, -4., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "xhad mass Vub"); h = new TH1D(name, title, XBIN, -4., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "xhad mass D    "); h = new TH1D(name, title, XBIN, -4., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "xhad mass D*   "); h = new TH1D(name, title, XBIN, -4., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "xhad mass D(*)x"); h = new TH1D(name, title, XBIN, -4., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "xhad mass other"); h = new TH1D(name, title, XBIN, -4., 10.); h->Sumw2();

      no = 2320;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "xhad mass");     h = new TH1D(name, title, XBIN, -4., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "xhad mass Vcb"); h = new TH1D(name, title, XBIN, -4., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "xhad mass Vub"); h = new TH1D(name, title, XBIN, -4., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "xhad mass D    "); h = new TH1D(name, title, XBIN, -4., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "xhad mass D*   "); h = new TH1D(name, title, XBIN, -4., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "xhad mass D(*)x"); h = new TH1D(name, title, XBIN, -4., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "xhad mass other"); h = new TH1D(name, title, XBIN, -4., 10.); h->Sumw2();
				     
				     
      no = 2400;		     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "xhad mass fit");     h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "xhad mass fit Vcb"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "xhad mass fit Vub"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "xhad mass fit D    "); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "xhad mass fit D*   "); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "xhad mass fit D(*)x"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "xhad mass fit other"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
				     
      no = 2410;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "xhad mass fit");     h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "xhad mass fit Vcb"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "xhad mass fit Vub"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "xhad mass fit D    "); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "xhad mass fit D*   "); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "xhad mass fit D(*)x"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "xhad mass fit other"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();

      no = 2420;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "xhad mass fit");     h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "xhad mass fit Vcb"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "xhad mass fit Vub"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "xhad mass fit D    "); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "xhad mass fit D*   "); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "xhad mass fit D(*)x"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "xhad mass fit other"); h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
				     
      no = 2500;		     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "xhad mass fit");     h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "xhad mass fit Vcb"); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "xhad mass fit Vub"); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "xhad mass fit D    "); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "xhad mass fit D*   "); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "xhad mass fit D(*)x"); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "xhad mass fit other"); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      				     
      no = 2510;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "xhad mass fit");     h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "xhad mass fit Vcb"); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "xhad mass fit Vub"); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "xhad mass fit D    "); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "xhad mass fit D*   "); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "xhad mass fit D(*)x"); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "xhad mass fit other"); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();

      no = 2520;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "xhad mass fit");     h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "xhad mass fit Vcb"); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "xhad mass fit Vub"); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "xhad mass fit D    "); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "xhad mass fit D*   "); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "xhad mass fit D(*)x"); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "xhad mass fit other"); h = new TH1D(name, title, XBIN/10, 0., XMAX); h->Sumw2();

      no = 3000;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "mm2");     h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "mm2 Vcb"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "mm2 Vub"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "mm2 D    "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "mm2 D*   "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "mm2 D(*)x"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "mm2 other"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
				     
      no = 3010;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "mm2");     h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "mm2 Vcb"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "mm2 Vub"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "mm2 D    "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "mm2 D*   "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "mm2 D(*)x"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "mm2 other"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();

      no = 3020;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "mm2");     h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "mm2 Vcb"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "mm2 Vub"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "mm2 D    "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "mm2 D*   "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "mm2 D(*)x"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "mm2 other"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();

      no = 3100;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "mm2");     h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "mm2 Vcb"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "mm2 Vub"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "mm2 D    "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "mm2 D*   "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "mm2 D(*)x"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "mm2 other"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
				     
      no = 3110;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "mm2");     h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "mm2 Vcb"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "mm2 Vub"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "mm2 D    "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "mm2 D*   "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "mm2 D(*)x"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "mm2 other"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();

      no = 3120;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "mm2");     h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "mm2 Vcb"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "mm2 Vub"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "mm2 D    "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "mm2 D*   "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "mm2 D(*)x"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "mm2 other"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();

      no = 3200;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "mm2");     h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "mm2 Vcb"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "mm2 Vub"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "mm2 D    "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "mm2 D*   "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "mm2 D(*)x"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "mm2 other"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
				     
      no = 3210;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "mm2");     h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "mm2 Vcb"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "mm2 Vub"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "mm2 D    "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "mm2 D*   "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "mm2 D(*)x"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "mm2 other"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();

      no = 3220;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "mm2");     h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "mm2 Vcb"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "mm2 Vub"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "mm2 D    "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "mm2 D*   "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "mm2 D(*)x"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "mm2 other"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();

      no = 3300;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "mm2");     h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "mm2 Vcb"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "mm2 Vub"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "mm2 D    "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "mm2 D*   "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "mm2 D(*)x"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "mm2 other"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
				     
      no = 3310;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "mm2");     h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "mm2 Vcb"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "mm2 Vub"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "mm2 D    "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "mm2 D*   "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "mm2 D(*)x"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "mm2 other"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();

      no = 3320;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "mm2");     h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "mm2 Vcb"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "mm2 Vub"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "mm2 D    "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "mm2 D*   "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "mm2 D(*)x"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "mm2 other"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();

      no = 3400;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "mm2");     h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "mm2 Vcb"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "mm2 Vub"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "mm2 D    "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "mm2 D*   "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "mm2 D(*)x"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "mm2 other"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
				     
      no = 3410;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "mm2");     h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "mm2 Vcb"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "mm2 Vub"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "mm2 D    "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "mm2 D*   "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "mm2 D(*)x"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "mm2 other"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();

      no = 3420;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "mm2");     h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "mm2 Vcb"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "mm2 Vub"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "mm2 D    "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "mm2 D*   "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "mm2 D(*)x"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "mm2 other"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();

      no = 3500;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "mm2");     h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "mm2 Vcb"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "mm2 Vub"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "mm2 D    "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "mm2 D*   "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "mm2 D(*)x"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "mm2 other"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
				     
      no = 3510;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "mm2");     h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "mm2 Vcb"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "mm2 Vub"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "mm2 D    "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "mm2 D*   "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "mm2 D(*)x"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "mm2 other"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();

      no = 3520;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "mm2");     h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "mm2 Vcb"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "mm2 Vub"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "mm2 D    "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "mm2 D*   "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "mm2 D(*)x"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "mm2 other"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();

      no = 3600;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "mm2");     h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "mm2 Vcb"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "mm2 Vub"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "mm2 D    "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "mm2 D*   "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "mm2 D(*)x"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "mm2 other"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
				     
      no = 3610;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "mm2");     h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "mm2 Vcb"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "mm2 Vub"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "mm2 D    "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "mm2 D*   "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "mm2 D(*)x"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "mm2 other"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();

      no = 3620;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "mm2");     h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "mm2 Vcb"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "mm2 Vub"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "mm2 D    "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "mm2 D*   "); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "mm2 D(*)x"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "mm2 other"); h = new TH1D(name, title, XBIN, -XMAX, XMAX); h->Sumw2();

      				     
//        no = 3100;           	     
//        sprintf(name,"%s%d",lc,no+0); sprintf(title, "cos#theta neutrino");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+1); sprintf(title, "cos#theta neutrino");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+2); sprintf(title, "cos#theta neutrino");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+3); sprintf(title, "cos#theta neutrino");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+4); sprintf(title, "cos#theta neutrino");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+5); sprintf(title, "cos#theta neutrino");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+6); sprintf(title, "cos#theta neutrino");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
				     
//        no = 3110;           	     
//        sprintf(name,"%s%d",lc,no+0); sprintf(title, "cos#theta neutrino");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+1); sprintf(title, "cos#theta neutrino");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+2); sprintf(title, "cos#theta neutrino");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+3); sprintf(title, "cos#theta neutrino");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+4); sprintf(title, "cos#theta neutrino");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+5); sprintf(title, "cos#theta neutrino");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+6); sprintf(title, "cos#theta neutrino");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();

//        no = 3120;           	     
//        sprintf(name,"%s%d",lc,no+0); sprintf(title, "cos#theta neutrino");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+1); sprintf(title, "cos#theta neutrino");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+2); sprintf(title, "cos#theta neutrino");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+3); sprintf(title, "cos#theta neutrino");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+4); sprintf(title, "cos#theta neutrino");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+5); sprintf(title, "cos#theta neutrino");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+6); sprintf(title, "cos#theta neutrino");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
      				     

//        no = 3200;           	     
//        sprintf(name,"%s%d",lc,no+0); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+1); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+2); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+3); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+4); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+5); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+6); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
				     
//        no = 3210;           	     
//        sprintf(name,"%s%d",lc,no+0); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+1); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+2); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+3); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+4); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+5); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+6); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();

//        no = 3220;           	     
//        sprintf(name,"%s%d",lc,no+0); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+1); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+2); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+3); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+4); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+5); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
//        sprintf(name,"%s%d",lc,no+6); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();

				     
      no = 3800;           	     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "Q2");  h = new TH1D(name, title, 30, 0., 15.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "Q2");  h = new TH1D(name, title, 30, 0., 15.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "Q2");  h = new TH1D(name, title, 30, 0., 15.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "Q2");  h = new TH1D(name, title, 30, 0., 15.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "Q2");  h = new TH1D(name, title, 30, 0., 15.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "Q2");  h = new TH1D(name, title, 30, 0., 15.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "Q2");  h = new TH1D(name, title, 30, 0., 15.); h->Sumw2();
				     
      no = 3810;           	     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "Q2");  h = new TH1D(name, title, 30, 0., 15.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "Q2");  h = new TH1D(name, title, 30, 0., 15.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "Q2");  h = new TH1D(name, title, 30, 0., 15.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "Q2");  h = new TH1D(name, title, 30, 0., 15.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "Q2");  h = new TH1D(name, title, 30, 0., 15.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "Q2");  h = new TH1D(name, title, 30, 0., 15.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "Q2");  h = new TH1D(name, title, 30, 0., 15.); h->Sumw2();

      no = 3820;           	     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "Q2");  h = new TH1D(name, title, 30, 0., 15.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "Q2");  h = new TH1D(name, title, 30, 0., 15.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "Q2");  h = new TH1D(name, title, 30, 0., 15.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "Q2");  h = new TH1D(name, title, 30, 0., 15.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "Q2");  h = new TH1D(name, title, 30, 0., 15.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "Q2");  h = new TH1D(name, title, 30, 0., 15.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "Q2");  h = new TH1D(name, title, 30, 0., 15.); h->Sumw2();


      no = 3900;           	     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "emax");  h = new TH1D(name, title, 30, 0., 3.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "emax");  h = new TH1D(name, title, 30, 0., 3.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "emax");  h = new TH1D(name, title, 30, 0., 3.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "emax");  h = new TH1D(name, title, 30, 0., 3.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "emax");  h = new TH1D(name, title, 30, 0., 3.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "emax");  h = new TH1D(name, title, 30, 0., 3.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "emax");  h = new TH1D(name, title, 30, 0., 3.); h->Sumw2();

      no = 3910;           	     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "emax");  h = new TH1D(name, title, 30, 0., 3.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "emax");  h = new TH1D(name, title, 30, 0., 3.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "emax");  h = new TH1D(name, title, 30, 0., 3.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "emax");  h = new TH1D(name, title, 30, 0., 3.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "emax");  h = new TH1D(name, title, 30, 0., 3.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "emax");  h = new TH1D(name, title, 30, 0., 3.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "emax");  h = new TH1D(name, title, 30, 0., 3.); h->Sumw2();

      no = 3920;           	     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "emax");  h = new TH1D(name, title, 30, 0., 3.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "emax");  h = new TH1D(name, title, 30, 0., 3.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "emax");  h = new TH1D(name, title, 30, 0., 3.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "emax");  h = new TH1D(name, title, 30, 0., 3.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "emax");  h = new TH1D(name, title, 30, 0., 3.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "emax");  h = new TH1D(name, title, 30, 0., 3.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "emax");  h = new TH1D(name, title, 30, 0., 3.); h->Sumw2();
      				     
      no = 4000;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "nTrk");      h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "nTrk Vcb");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "nTrk Vub");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "nTrk D    ");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "nTrk D*   ");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "nTrk D(*)x");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "nTrk other");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();

      no = 4010;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "nTrk");      h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "nTrk Vcb");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "nTrk Vub");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "nTrk D    ");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "nTrk D*   ");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "nTrk D(*)x");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "nTrk other");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
				     
      no = 4020;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "nTrk");      h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "nTrk Vcb");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "nTrk Vub");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "nTrk D    ");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "nTrk D*   ");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "nTrk D(*)x");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "nTrk other");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      				     
      no = 4100;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "nNut");      h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "nNut Vcb");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "nNut Vub");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "nNut D    ");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "nNut D*   ");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "nNut D(*)x");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "nNut other");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
				     
      no = 4110;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "nNut");      h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "nNut Vcb");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "nNut Vub");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "nNut D    ");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "nNut D*   ");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "nNut D(*)x");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "nNut other");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();

      no = 4120;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "nNut");      h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "nNut Vcb");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "nNut Vub");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "nNut D    ");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "nNut D*   ");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "nNut D(*)x");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "nNut other");  h = new TH1D(name, title, 20, 0., 20.); h->Sumw2();

      				     
      no = 4200;           	     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "chargeRecoil");      h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "chargeRecoil Vcb");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2(); 
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "chargeRecoil Vub");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "chargeRecoil D    ");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "chargeRecoil D*   ");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "chargeRecoil D(*)x");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "chargeRecoil other");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();

      no = 4210;           	     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "chargeRecoil");      h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "chargeRecoil Vcb");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2(); 
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "chargeRecoil Vub");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "chargeRecoil D    ");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "chargeRecoil D*   ");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "chargeRecoil D(*)x");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "chargeRecoil other");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();

      no = 4220;           	     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "chargeRecoil");      h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "chargeRecoil Vcb");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2(); 
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "chargeRecoil Vub");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "chargeRecoil D    ");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "chargeRecoil D*   ");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "chargeRecoil D(*)x");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "chargeRecoil other");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();

      no = 4300;           	     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "charge total");      h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "charge total Vcb");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2(); 
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "charge total Vub");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "charge total D    ");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "charge total D*   ");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "charge total D(*)x");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "charge total other");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();
				     			    
      no = 4310;           	     			    
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "charge total");      h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "charge total Vcb");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2(); 
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "charge total Vub");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "charge total D    ");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "charge total D*   ");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "charge total D(*)x");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "charge total other");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();
							    
      no = 4320;           	     			    
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "charge total");      h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "charge total Vcb");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2(); 
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "charge total Vub");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "charge total D    ");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "charge total D*   ");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "charge total D(*)x");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "charge total other");  h = new TH1D(name, title, 10, -5., 5.); h->Sumw2();

      if (fOptCategories) {

	static const char* mxcattitle[11] = { "Cascade ","KL ","KS ","neutron ","misid K", "charged K", "KL & cascade", "KL no cascade ", "KS  & cascade", "Ks no cascade", "cascade not K"};
	static const char* sample[7] = {" ", " Vcb", " Vub", " D", " D*", "D(*)x", "other"};
	  
	no = 5000;

	for (Int_t k=0; k<7; k++) {
	  for (Int_t j=0; j<11; j++) {
	      
	    sprintf(name,"%s%d",lc,no+k+j*10+10); sprintf(title, "%s xhad mass %s",mxcattitle[j], sample[k]);
	      
	    h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
	  }
	  sprintf(name,"%s%d",lc,no+k+500); sprintf(title, "Best KS subtracted xhad mass %s", sample[k]);
	  h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
	  sprintf(name,"%s%d",lc,no+k+600); sprintf(title, "MC Category bitmap%s", sample[k]);
	  h = new TH1D(name, title, 65, -0.5, 64.5); h->Sumw2();
	  sprintf(name,"%s%d",lc,no+k); sprintf(title, "MC Category %s", sample[k]);
	  h = new TH1D(name, title, 12, -0.5, 11.5); h->Sumw2();
	}
	  
	no=6000;
	  
	for (Int_t k=0; k<7; k++) {
	  sprintf(name,"%s%d",lc,no+k); sprintf(title, "MC Category %s", sample[k]);
	  h = new TH1D(name, title, 12, -0.5, 11.5); h->Sumw2();
	}
	  
	no=7000;
	for (Int_t k=0; k<7; k++) {
	  for (Int_t j=0; j<11; j++) {
	      
	    sprintf(name,"%s%d",lc,no+k+j*10+10); sprintf(title, "%s xhad fit mass %s",mxcattitle[j], sample[k]);
	      
	    h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
	  }
	  sprintf(name,"%s%d",lc,no+k+500); sprintf(title, "Best KS subtracted xhad mass fit %s", sample[k]);
	  h = new TH1D(name, title, XBIN, 0., XMAX); h->Sumw2();
	  sprintf(name,"%s%d",lc,no+k+600); sprintf(title, "MC Category bitmap%s", sample[k]);
	  h = new TH1D(name, title, 65, -0.5, 64.5); h->Sumw2();
	  sprintf(name,"%s%d",lc,no+k); sprintf(title, "MC Category %s", sample[k]);
	  h = new TH1D(name, title, 12, -0.5, 11.5); h->Sumw2();
	}
	  
	no=8000;
	  
	for (Int_t k=0; k<7; k++) {
	  sprintf(name,"%s%d",lc,no+k); sprintf(title, "MC Category %s", sample[k]);
	  h = new TH1D(name, title, 12, -0.5, 11.5); h->Sumw2();
	}
      }
	  
      no = 9000;
          

      sprintf(name,"%s%d",lc,no+0); sprintf(title, "res xhad mass");      h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "res xhad mass Vcb");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2(); 
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "res xhad mass Vub");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "res xhad mass D    ");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "res xhad mass D*   ");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "res xhad mass D(*)x");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "res xhad mass other");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
				     
      no = 9010;           	     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "res xhad mass");      h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "res xhad mass Vcb");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2(); 
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "res xhad mass Vub");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "res xhad mass D    ");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "res xhad mass D*   ");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "res xhad mass D(*)x");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "res xhad mass other");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();

      no = 9020;           	     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "res xhad mass");      h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "res xhad mass Vcb");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2(); 
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "res xhad mass Vub");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "res xhad mass D    ");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "res xhad mass D*   ");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "res xhad mass D(*)x");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "res xhad mass other");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();

      no = 9100;           	     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "res pmiss");      h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "res pmiss Vcb");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2(); 
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "res pmiss Vub");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "res pmiss D    ");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "res pmiss D*   ");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "res pmiss D(*)x");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "res pmiss other");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();

      no = 9110;           	     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "res pmiss");      h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "res pmiss Vcb");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2(); 
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "res pmiss Vub");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "res pmiss D    ");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "res pmiss D*   ");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "res pmiss D(*)x");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "res pmiss other");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();

      no = 9120;           	     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "res pmiss");      h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "res pmiss Vcb");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2(); 
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "res pmiss Vub");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "res pmiss D    ");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "res pmiss D*   ");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "res pmiss D(*)x");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "res pmiss other");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();

      no = 9200;           	     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "res tmiss");      h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "res tmiss Vcb");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2(); 
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "res tmiss Vub");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "res tmiss D    ");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "res tmiss D*   ");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "res tmiss D(*)x");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "res tmiss other");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();

      no = 9210;           	     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "res tmiss");      h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "res tmiss Vcb");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2(); 
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "res tmiss Vub");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "res tmiss D    ");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "res tmiss D*   ");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "res tmiss D(*)x");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "res tmiss other");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();

      no = 9220;           	     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "res tmiss");      h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "res tmiss Vcb");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2(); 
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "res tmiss Vub");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "res tmiss D    ");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "res tmiss D*   ");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "res tmiss D(*)x");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "res tmiss other");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();

      no = 9300;           	     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "res Q2");      h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "res Q2 Vcb");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2(); 
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "res Q2 Vub");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "res Q2 D    ");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "res Q2 D*   ");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "res Q2 D(*)x");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "res Q2 other");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();

      no = 9310;           	     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "res Q2");      h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "res Q2 Vcb");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2(); 
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "res Q2 Vub");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "res Q2 D    ");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "res Q2 D*   ");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "res Q2 D(*)x");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "res Q2 other");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();

      no = 9320;           	     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "res Q2");      h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "res Q2 Vcb");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2(); 
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "res Q2 Vub");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "res Q2 D    ");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "res Q2 D*   ");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "res Q2 D(*)x");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "res Q2 other");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();


      no = 9400;
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "res xhad fitted mass");      h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "res xhad fitted mass Vcb");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2(); 
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "res xhad fitted mass Vub");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "res xhad fitted mass D    ");h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "res xhad fitted mass D*   ");h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "res xhad fitted mass D(*)x");h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "res xhad fitted mass other");h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
				     
      no = 9410;           	     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "res xhad fitted mass");      h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "res xhad fitted mass Vcb");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2(); 
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "res xhad fitted mass Vub");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "res xhad fitted mass D    ");h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "res xhad fitted mass D*   ");h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "res xhad fitted mass D(*)x");h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "res xhad fitted mass other");h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();

      no = 9420;           	     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "res xhad fitted mass");      h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "res xhad fitted mass Vcb");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2(); 
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "res xhad fitted mass Vub");  h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "res xhad fitted mass D    ");h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "res xhad fitted mass D*   ");h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "res xhad fitted mass D(*)x");h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "res xhad fitted mass other");h = new TH1D(name, title, 50, -5., 5.); h->Sumw2();
				     

    }

    sprintf(name, "NuTheta");  sprintf(title, "Cos(Theta) angle of neutrino momentum");  h = new TH1D(name, title, 50, -1., 1.0); h->Sumw2();
    sprintf(name, "NuPhi");  sprintf(title, "Phi angle of neutrino momentum");  h = new TH1D(name, title, 50, -3.15, 3.15); h->Sumw2();
    sprintf(name, "mNuSqNC");  sprintf(title, "Missing mass squared");  h = new TH1D(name, title, 100, -3., 8.0); h->Sumw2();
    sprintf(name, "nc_mNuSqNC");  sprintf(title, "Missing mass squared (no cut)");  h = new TH1D(name, title, 100, -3., 8.0); h->Sumw2();
    sprintf(name, "mes all");  sprintf(title, "mes distribution");  h = new TH1D(name, title, 50, 5.2, 5.3); 
    
    if (fOptGammas) {
      gDirectory->cd("CutPlots");
      for(int jv =0; jv<15; jv++) {
	sprintf(name, "Mnu%d", jv);  sprintf(title, "Missing mass squared");  h = new TH1D(name, title, 4, -1.8, 3.0); h->Sumw2();
	sprintf(name, "bin_Mnu%d", jv);  sprintf(title, "Missing mass squared");  h = new TH1D(name, title, 100, -4., 8.0); h->Sumw2();
	sprintf(name, "nc_Mnu%d", jv);  sprintf(title, "Missing mass squared");  h = new TH1D(name, title, 4, -1.8, 3.0); h->Sumw2();
	sprintf(name, "bin_nc_Mnu%d", jv);  sprintf(title, "Missing mass squared");  h = new TH1D(name, title, 100, -4., 8.0); h->Sumw2();
	//Reso block
	sprintf(name, "ResMX%d",jv);  sprintf(title, "Reso Mx");  h = new TH1D(name, title, 100, -4.0, 4.0); h->Sumw2();
	sprintf(name, "nc_ResMX%d",jv);  sprintf(title, "nc_Reso Mx");  h = new TH1D(name, title, 100, -4.0, 4.0); h->Sumw2();
	sprintf(name, "bin_ResMX%d",jv);  sprintf(title, "Reso Mx");  h = new TH1D(name, title, 8, -4.0, 4.0); h->Sumw2();
	sprintf(name, "bin_nc_ResMX%d",jv);  sprintf(title, "nc_Reso Mx");  h = new TH1D(name, title, 8, -4.0, 4.0); h->Sumw2();
	sprintf(name, "ResMXF%d",jv);  sprintf(title, "Reso Mx");  h = new TH1D(name, title, 100, -4.0, 4.0); h->Sumw2();
	sprintf(name, "nc_ResMXF%d",jv);  sprintf(title, "nc_Reso Mx");  h = new TH1D(name, title, 100, -4.0, 4.0); h->Sumw2();
	sprintf(name, "bin_ResMXF%d",jv);  sprintf(title, "Reso Mx");  h = new TH1D(name, title, 8, -4.0, 4.0); h->Sumw2();
	sprintf(name, "bin_nc_ResMXF%d",jv);  sprintf(title, "nc_Reso Mx");  h = new TH1D(name, title, 8, -4.0, 4.0); h->Sumw2();
	// cout<<"Problems with sub directories"<<endl;
	
	for(int iv =0; iv<10; iv++) {
	  sprintf(name, "MNH%d%d", jv,iv);  sprintf(title, "Missing mass squared");  h = new TH1D(name, title, 4, -1.8, 3.0); h->Sumw2();
	  sprintf(name, "bin_MNH%d%d", jv,iv);  sprintf(title, "Missing mass squared");  h = new TH1D(name, title, 100, -4., 8.0); h->Sumw2();
	  //Reso block
	  sprintf(name, "ResNH%d%d",jv,iv);  sprintf(title, "Reso Mx");  h = new TH1D(name, title, 100, -2.0, 2.0); h->Sumw2();
	  sprintf(name, "bin_ResNH%d%d",jv,iv);  sprintf(title, "Reso Mx");  h = new TH1D(name, title, 16, -2.0, 2.0); h->Sumw2();
	  sprintf(name, "MxNH%d%d",jv,iv);  sprintf(title, "Mx");  h = new TH1D(name, title, 100, 0.0, 4.5); h->Sumw2();
	  sprintf(name, "bin_MxNH%d%d",jv,iv);  sprintf(title, "Mx");  h = new TH1D(name, title, 3, 0.0, 4.5); h->Sumw2();
	  //Reso block Fitted
	  sprintf(name, "ResNHF%d%d",jv,iv);  sprintf(title, "Reso Mx");  h = new TH1D(name, title, 100, -2.0, 2.0); h->Sumw2();
	  sprintf(name, "bin_ResNHF%d%d",jv,iv);  sprintf(title, "Reso Mx");  h = new TH1D(name, title, 16, -2.0, 2.0); h->Sumw2();
	  sprintf(name, "MxNHF%d%d",jv,iv);  sprintf(title, "Mx");  h = new TH1D(name, title, 100, 0.0, 4.5); h->Sumw2();
	    sprintf(name, "bin_MxNHF%d%d",jv,iv);  sprintf(title, "Mx");  h = new TH1D(name, title, 3, 0.0, 4.5); h->Sumw2();
	}
	
	//	fHistFile->cd();
      }
    }
    fHistFile->cd();
    gDirectory->cd();
  }
  bookKLHisto();
}



// ----------------------------------------------------------------------
void recoilBase::switchOffReading(const char *s) {
  fChain->SetBranchStatus(s,0);
}


// ----------------------------------------------------------------------
Int_t recoilBase::GetEntry(Int_t entry) {
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

// ----------------------------------------------------------------------
Int_t recoilBase::LoadTree(Int_t entry) {
// Set the environment to read one entry
   if (!fChain) return -5;
   Int_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->IsA() != TChain::Class()) return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

// ----------------------------------------------------------------------
void recoilBase::Init(TTree *tree, int isMC, int mcasdata) {
//   Set branch addresses
  
  if (isMC) { 
    fIsMC = kTRUE; 
  } else {
    fIsMC = kFALSE; 
  }
   if (tree == 0) return;
   fChain    = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event",&event);

   READALL=1;

   if(READALL==1) {
     fChain->SetBranchAddress("runNumber",&runNumber);
     fChain->SetBranchAddress("platform",&platform);
     fChain->SetBranchAddress("partition",&partition);
     fChain->SetBranchAddress("upperID",&upperID);
     fChain->SetBranchAddress("lowerID",&lowerID);
   }

   fChain->SetBranchAddress("primVtxX",&primVtxX);
   fChain->SetBranchAddress("primVtxY",&primVtxY);
   fChain->SetBranchAddress("primVtxZ",&primVtxZ);
   fChain->SetBranchAddress("beamSX",&beamSX);
   fChain->SetBranchAddress("beamSY",&beamSY);
   fChain->SetBranchAddress("beamSZ",&beamSZ);

   if(READALL==1) {   
     fChain->SetBranchAddress("beamSCovXX",&beamSCovXX);
     fChain->SetBranchAddress("beamSCovYY",&beamSCovYY);
     fChain->SetBranchAddress("beamSCovZZ",&beamSCovZZ);
     fChain->SetBranchAddress("beamSCovXZ",&beamSCovXZ);
   }

   fChain->SetBranchAddress("pxUps",&pxUps);
   fChain->SetBranchAddress("pyUps",&pyUps);
   fChain->SetBranchAddress("pzUps",&pzUps);
   fChain->SetBranchAddress("eUps",&eUps);
   fChain->SetBranchAddress("nTrkTot",&nTrkTot);

   if(READALL==1) {   
     fChain->SetBranchAddress("W2",&W2);
     fChain->SetBranchAddress("FoxWol2",&FoxWol2);
     fChain->SetBranchAddress("FoxWol2Neu",&FoxWol2Neu);
     fChain->SetBranchAddress("thrust",&thrust);
     fChain->SetBranchAddress("thrustNeu",&thrustNeu);
   }

   //MC block
   if(isMC || mcasdata) { 
     fChain->SetBranchAddress("nMc",&nMc);
     fChain->SetBranchAddress("pMc",pMc);
     fChain->SetBranchAddress("massMc",massMc);
     fChain->SetBranchAddress("thetaMc",thetaMc);
     fChain->SetBranchAddress("phiMc",phiMc);
     fChain->SetBranchAddress("idMc",idMc);
     fChain->SetBranchAddress("mothMc",mothMc);
     fChain->SetBranchAddress("nDauMc",nDauMc);
     fChain->SetBranchAddress("xMc",xMc);
     fChain->SetBranchAddress("yMc",yMc);
     fChain->SetBranchAddress("zMc",zMc);
   }
   fChain->SetBranchAddress("nB0",&nB0);
   fChain->SetBranchAddress("massB0",massB0);
   fChain->SetBranchAddress("pB0",pB0);
   fChain->SetBranchAddress("thB0",thB0);
   fChain->SetBranchAddress("phiB0",phiB0);

   if(READALL==1) {   
     fChain->SetBranchAddress("errMassB0",errMassB0);
     fChain->SetBranchAddress("m0B0",m0B0);
     fChain->SetBranchAddress("xB0",xB0);
     fChain->SetBranchAddress("yB0",yB0);
     fChain->SetBranchAddress("zB0",zB0);
     fChain->SetBranchAddress("s2xB0",s2xB0);
     fChain->SetBranchAddress("s2yB0",s2yB0);
     fChain->SetBranchAddress("s2zB0",s2zB0);
     fChain->SetBranchAddress("chi2B0",chi2B0);
     fChain->SetBranchAddress("dofB0",dofB0);
     fChain->SetBranchAddress("stB0",stB0);
     fChain->SetBranchAddress("ndauB0",ndauB0);
     fChain->SetBranchAddress("mHatB0",mHatB0);
     fChain->SetBranchAddress("ThruB0",ThruB0);
     fChain->SetBranchAddress("thThruB0",thThruB0);
     fChain->SetBranchAddress("phiThruB0",phiThruB0);
     fChain->SetBranchAddress("cosTBB0",cosTBB0);
   }
   if(isMC || mcasdata)   fChain->SetBranchAddress("MCB0",MCB0);
   fChain->SetBranchAddress("mseB0",mseB0);
   fChain->SetBranchAddress("deltaeB0",deltaeB0);
   
   fChain->SetBranchAddress("d1B0Index",d1B0Index);
   fChain->SetBranchAddress("d1B0Lund",d1B0Lund);
   fChain->SetBranchAddress("d2B0Index",d2B0Index);
   fChain->SetBranchAddress("d2B0Lund",d2B0Lund);
   fChain->SetBranchAddress("d3B0Index",d3B0Index);
   fChain->SetBranchAddress("d3B0Lund",d3B0Lund);
   fChain->SetBranchAddress("d4B0Index",d4B0Index);
   fChain->SetBranchAddress("d4B0Lund",d4B0Lund);
   fChain->SetBranchAddress("d5B0Index",d5B0Index);
   fChain->SetBranchAddress("d5B0Lund",d5B0Lund);
   fChain->SetBranchAddress("d6B0Index",d6B0Index);
   fChain->SetBranchAddress("d6B0Lund",d6B0Lund);
   fChain->SetBranchAddress("d7B0Index",d7B0Index);
   fChain->SetBranchAddress("d7B0Lund",d7B0Lund);
   fChain->SetBranchAddress("modeB0",modeB0);
   fChain->SetBranchAddress("purB0",purB0);
   fChain->SetBranchAddress("intpurB0",intpurB0);

   if(READALL==1) {   
     fChain->SetBranchAddress("VtxXLepB0",VtxXLepB0);
     fChain->SetBranchAddress("VtxYLepB0",VtxYLepB0);
     fChain->SetBranchAddress("VtxZLepB0",VtxZLepB0);
     fChain->SetBranchAddress("VtxCovXXLepB0",VtxCovXXLepB0);
     fChain->SetBranchAddress("VtxCovYYLepB0",VtxCovYYLepB0);
     fChain->SetBranchAddress("VtxCovXYLepB0",VtxCovXYLepB0);
     fChain->SetBranchAddress("VtxCovZZLepB0",VtxCovZZLepB0);
     fChain->SetBranchAddress("VtxCovXZLepB0",VtxCovXZLepB0);
     fChain->SetBranchAddress("VtxCovYZLepB0",VtxCovYZLepB0);
     fChain->SetBranchAddress("VtxChiSqLepB0",VtxChiSqLepB0);
     fChain->SetBranchAddress("VtxNDofLepB0",VtxNDofLepB0);
     fChain->SetBranchAddress("VtxStatLepB0",VtxStatLepB0);
     fChain->SetBranchAddress("VtxNUsedLepB0",VtxNUsedLepB0);
     fChain->SetBranchAddress("DocaLepB0",DocaLepB0);
     fChain->SetBranchAddress("DocaErrLepB0",DocaErrLepB0);
     fChain->SetBranchAddress("VtxXXB0",VtxXXB0);
     fChain->SetBranchAddress("VtxYXB0",VtxYXB0);
     fChain->SetBranchAddress("VtxZXB0",VtxZXB0);
     fChain->SetBranchAddress("VtxCovXXXB0",VtxCovXXXB0);
     fChain->SetBranchAddress("VtxCovYYXB0",VtxCovYYXB0);
     fChain->SetBranchAddress("VtxCovXYXB0",VtxCovXYXB0);
     fChain->SetBranchAddress("VtxCovZZXB0",VtxCovZZXB0);
     fChain->SetBranchAddress("VtxCovXZXB0",VtxCovXZXB0);
     fChain->SetBranchAddress("VtxCovYZXB0",VtxCovYZXB0);
     fChain->SetBranchAddress("VtxChiSqXB0",VtxChiSqXB0);
     fChain->SetBranchAddress("VtxNDofXB0",VtxNDofXB0);
     fChain->SetBranchAddress("VtxStatXB0",VtxStatXB0);
     fChain->SetBranchAddress("VtxNUsedXB0",VtxNUsedXB0);
     fChain->SetBranchAddress("VtxPXB0",VtxPXB0);
     fChain->SetBranchAddress("VtxPhiXB0",VtxPhiXB0);
     fChain->SetBranchAddress("VtxThetaXB0",VtxThetaXB0);
     fChain->SetBranchAddress("ThrustXB0",ThrustXB0);
     fChain->SetBranchAddress("ThrustXPhiB0",ThrustXPhiB0);
     fChain->SetBranchAddress("ThrustXThetaB0",ThrustXThetaB0);
   }

   fChain->SetBranchAddress("MassPB0",MassPB0);
   fChain->SetBranchAddress("MassPhiB0",MassPhiB0);
   fChain->SetBranchAddress("MassThetaB0",MassThetaB0);

   if(READALL==1) {   
     fChain->SetBranchAddress("Cov00B0",Cov00B0);
     fChain->SetBranchAddress("Cov10B0",Cov10B0);
     fChain->SetBranchAddress("Cov11B0",Cov11B0);
     fChain->SetBranchAddress("Cov20B0",Cov20B0);
     fChain->SetBranchAddress("Cov21B0",Cov21B0);
     fChain->SetBranchAddress("Cov22B0",Cov22B0);
     fChain->SetBranchAddress("Cov30B0",Cov30B0);
     fChain->SetBranchAddress("Cov31B0",Cov31B0);
     fChain->SetBranchAddress("Cov32B0",Cov32B0);
     fChain->SetBranchAddress("Cov33B0",Cov33B0);
   }

   fChain->SetBranchAddress("nChB",&nChB);
   fChain->SetBranchAddress("massChB",massChB);
   fChain->SetBranchAddress("pChB",pChB);
   fChain->SetBranchAddress("thChB",thChB);
   fChain->SetBranchAddress("phiChB",phiChB);

   if(READALL==1) {   
     fChain->SetBranchAddress("errMassChB",errMassChB);
     fChain->SetBranchAddress("m0ChB",m0ChB);
     fChain->SetBranchAddress("xChB",xChB);
     fChain->SetBranchAddress("yChB",yChB);
     fChain->SetBranchAddress("zChB",zChB);
     fChain->SetBranchAddress("s2xChB",s2xChB);
     fChain->SetBranchAddress("s2yChB",s2yChB);
     fChain->SetBranchAddress("s2zChB",s2zChB);
     fChain->SetBranchAddress("chi2ChB",chi2ChB);
     fChain->SetBranchAddress("dofChB",dofChB);
     fChain->SetBranchAddress("stChB",stChB);
     fChain->SetBranchAddress("ndauChB",ndauChB);
     fChain->SetBranchAddress("mHatChB",mHatChB);
     fChain->SetBranchAddress("ThruChB",ThruChB);
     fChain->SetBranchAddress("thThruChB",thThruChB);
     fChain->SetBranchAddress("phiThruChB",phiThruChB);
     fChain->SetBranchAddress("cosTBChB",cosTBChB);
   }
   if(isMC || mcasdata)   fChain->SetBranchAddress("MCChB",MCChB);
   fChain->SetBranchAddress("mseChB",mseChB);
   fChain->SetBranchAddress("deltaeChB",deltaeChB);

   fChain->SetBranchAddress("d1ChBIndex",d1ChBIndex);
   fChain->SetBranchAddress("d1ChBLund",d1ChBLund);
   fChain->SetBranchAddress("d2ChBIndex",d2ChBIndex);
   fChain->SetBranchAddress("d2ChBLund",d2ChBLund);
   fChain->SetBranchAddress("d3ChBIndex",d3ChBIndex);
   fChain->SetBranchAddress("d3ChBLund",d3ChBLund);
   fChain->SetBranchAddress("d4ChBIndex",d4ChBIndex);
   fChain->SetBranchAddress("d4ChBLund",d4ChBLund);
   fChain->SetBranchAddress("d5ChBIndex",d5ChBIndex);
   fChain->SetBranchAddress("d5ChBLund",d5ChBLund);
   fChain->SetBranchAddress("d6ChBIndex",d6ChBIndex);
   fChain->SetBranchAddress("d6ChBLund",d6ChBLund);
   fChain->SetBranchAddress("d7ChBIndex",d7ChBIndex);
   fChain->SetBranchAddress("d7ChBLund",d7ChBLund);
   fChain->SetBranchAddress("modeChB",modeChB);
   fChain->SetBranchAddress("purChB",purChB);
   fChain->SetBranchAddress("intpurChB",intpurChB);

   if(READALL==1) {   
     fChain->SetBranchAddress("VtxXLepChB",VtxXLepChB);
     fChain->SetBranchAddress("VtxYLepChB",VtxYLepChB);
     fChain->SetBranchAddress("VtxZLepChB",VtxZLepChB);
     fChain->SetBranchAddress("VtxCovXXLepChB",VtxCovXXLepChB);
     fChain->SetBranchAddress("VtxCovYYLepChB",VtxCovYYLepChB);
     fChain->SetBranchAddress("VtxCovXYLepChB",VtxCovXYLepChB);
     fChain->SetBranchAddress("VtxCovZZLepChB",VtxCovZZLepChB);
     fChain->SetBranchAddress("VtxCovXZLepChB",VtxCovXZLepChB);
     fChain->SetBranchAddress("VtxCovYZLepChB",VtxCovYZLepChB);
     fChain->SetBranchAddress("VtxChiSqLepChB",VtxChiSqLepChB);
     fChain->SetBranchAddress("VtxNDofLepChB",VtxNDofLepChB);
     fChain->SetBranchAddress("VtxStatLepChB",VtxStatLepChB);
     fChain->SetBranchAddress("VtxNUsedLepChB",VtxNUsedLepChB);
     fChain->SetBranchAddress("DocaLepChB",DocaLepChB);
     fChain->SetBranchAddress("DocaErrLepChB",DocaErrLepChB);
     fChain->SetBranchAddress("VtxXXChB",VtxXXChB);
     fChain->SetBranchAddress("VtxYXChB",VtxYXChB);
     fChain->SetBranchAddress("VtxZXChB",VtxZXChB);
     fChain->SetBranchAddress("VtxCovXXXChB",VtxCovXXXChB);
     fChain->SetBranchAddress("VtxCovYYXChB",VtxCovYYXChB);
     fChain->SetBranchAddress("VtxCovXYXChB",VtxCovXYXChB);
     fChain->SetBranchAddress("VtxCovZZXChB",VtxCovZZXChB);
     fChain->SetBranchAddress("VtxCovXZXChB",VtxCovXZXChB);
     fChain->SetBranchAddress("VtxCovYZXChB",VtxCovYZXChB);
     fChain->SetBranchAddress("VtxChiSqXChB",VtxChiSqXChB);
     fChain->SetBranchAddress("VtxNDofXChB",VtxNDofXChB);
     fChain->SetBranchAddress("VtxStatXChB",VtxStatXChB);
     fChain->SetBranchAddress("VtxNUsedXChB",VtxNUsedXChB);
     fChain->SetBranchAddress("VtxPXChB",VtxPXChB);
     fChain->SetBranchAddress("VtxPhiXChB",VtxPhiXChB);
     fChain->SetBranchAddress("VtxThetaXChB",VtxThetaXChB);
     fChain->SetBranchAddress("ThrustXChB",ThrustXChB);
     fChain->SetBranchAddress("ThrustXPhiChB",ThrustXPhiChB);
     fChain->SetBranchAddress("ThrustXThetaChB",ThrustXThetaChB);
   }

   fChain->SetBranchAddress("MassPChB",MassPChB);
   fChain->SetBranchAddress("MassPhiChB",MassPhiChB);
   fChain->SetBranchAddress("MassThetaChB",MassThetaChB);

   if(READALL==1) {   
     fChain->SetBranchAddress("Cov00ChB",Cov00ChB);
     fChain->SetBranchAddress("Cov10ChB",Cov10ChB);
     fChain->SetBranchAddress("Cov11ChB",Cov11ChB);
     fChain->SetBranchAddress("Cov20ChB",Cov20ChB);
     fChain->SetBranchAddress("Cov21ChB",Cov21ChB);
     fChain->SetBranchAddress("Cov22ChB",Cov22ChB);
     fChain->SetBranchAddress("Cov30ChB",Cov30ChB);
     fChain->SetBranchAddress("Cov31ChB",Cov31ChB);
     fChain->SetBranchAddress("Cov32ChB",Cov32ChB);
     fChain->SetBranchAddress("Cov33ChB",Cov33ChB);
     fChain->SetBranchAddress("nDstar",&nDstar);
     fChain->SetBranchAddress("massDstar",massDstar);
     fChain->SetBranchAddress("pDstar",pDstar);
     fChain->SetBranchAddress("thDstar",thDstar);
     fChain->SetBranchAddress("phiDstar",phiDstar);
     fChain->SetBranchAddress("errMassDstar",errMassDstar);
     fChain->SetBranchAddress("m0Dstar",m0Dstar);
     fChain->SetBranchAddress("xDstar",xDstar);
     fChain->SetBranchAddress("yDstar",yDstar);
     fChain->SetBranchAddress("zDstar",zDstar);
     fChain->SetBranchAddress("s2xDstar",s2xDstar);
     fChain->SetBranchAddress("s2yDstar",s2yDstar);
     fChain->SetBranchAddress("s2zDstar",s2zDstar);
     fChain->SetBranchAddress("chi2Dstar",chi2Dstar);
     fChain->SetBranchAddress("dofDstar",dofDstar);
     fChain->SetBranchAddress("stDstar",stDstar);
     fChain->SetBranchAddress("ndauDstar",ndauDstar);
     if(isMC || mcasdata)   fChain->SetBranchAddress("MCDstar",MCDstar);
     fChain->SetBranchAddress("d1DstarIndex",d1DstarIndex);
     fChain->SetBranchAddress("d1DstarLund",d1DstarLund);
     fChain->SetBranchAddress("d2DstarIndex",d2DstarIndex);
     fChain->SetBranchAddress("d2DstarLund",d2DstarLund);
     fChain->SetBranchAddress("nDstarBS",&nDstarBS);
     fChain->SetBranchAddress("massDstarBS",massDstarBS);
     fChain->SetBranchAddress("chi2DstarBS",chi2DstarBS);
     fChain->SetBranchAddress("dofDstarBS",dofDstarBS);
     fChain->SetBranchAddress("spixDstarBS",spixDstarBS);
     fChain->SetBranchAddress("spiyDstarBS",spiyDstarBS);
     fChain->SetBranchAddress("spizDstarBS",spizDstarBS);
     fChain->SetBranchAddress("nDstar0",&nDstar0);
     fChain->SetBranchAddress("massDstar0",massDstar0);
     fChain->SetBranchAddress("pDstar0",pDstar0);
     fChain->SetBranchAddress("thDstar0",thDstar0);
     fChain->SetBranchAddress("phiDstar0",phiDstar0);
     fChain->SetBranchAddress("errMassDstar0",errMassDstar0);
     fChain->SetBranchAddress("m0Dstar0",m0Dstar0);
     fChain->SetBranchAddress("xDstar0",xDstar0);
     fChain->SetBranchAddress("yDstar0",yDstar0);
     fChain->SetBranchAddress("zDstar0",zDstar0);
     fChain->SetBranchAddress("s2xDstar0",s2xDstar0);
     fChain->SetBranchAddress("s2yDstar0",s2yDstar0);
     fChain->SetBranchAddress("s2zDstar0",s2zDstar0);
     fChain->SetBranchAddress("chi2Dstar0",chi2Dstar0);
     fChain->SetBranchAddress("dofDstar0",dofDstar0);
     fChain->SetBranchAddress("stDstar0",stDstar0);
     fChain->SetBranchAddress("ndauDstar0",ndauDstar0);
     if(isMC || mcasdata)   fChain->SetBranchAddress("MCDstar0",MCDstar0);
     fChain->SetBranchAddress("d1Dstar0Index",d1Dstar0Index);
     fChain->SetBranchAddress("d1Dstar0Lund",d1Dstar0Lund);
     fChain->SetBranchAddress("d2Dstar0Index",d2Dstar0Index);
     fChain->SetBranchAddress("d2Dstar0Lund",d2Dstar0Lund);
     fChain->SetBranchAddress("nD0",&nD0);
     fChain->SetBranchAddress("massD0",massD0);
     fChain->SetBranchAddress("pD0",pD0);
     fChain->SetBranchAddress("thD0",thD0);
     fChain->SetBranchAddress("phiD0",phiD0);
     fChain->SetBranchAddress("errMassD0",errMassD0);
     fChain->SetBranchAddress("m0D0",m0D0);
     fChain->SetBranchAddress("xD0",xD0);
     fChain->SetBranchAddress("yD0",yD0);
     fChain->SetBranchAddress("zD0",zD0);
     fChain->SetBranchAddress("s2xD0",s2xD0);
     fChain->SetBranchAddress("s2yD0",s2yD0);
     fChain->SetBranchAddress("s2zD0",s2zD0);
     fChain->SetBranchAddress("chi2D0",chi2D0);
     fChain->SetBranchAddress("dofD0",dofD0);
     fChain->SetBranchAddress("stD0",stD0);
     fChain->SetBranchAddress("ndauD0",ndauD0);
     if(isMC || mcasdata)   fChain->SetBranchAddress("MCD0",MCD0);
     fChain->SetBranchAddress("d1D0Index",d1D0Index);
     fChain->SetBranchAddress("d1D0Lund",d1D0Lund);
     fChain->SetBranchAddress("d2D0Index",d2D0Index);
     fChain->SetBranchAddress("d2D0Lund",d2D0Lund);
     fChain->SetBranchAddress("d3D0Index",d3D0Index);
     fChain->SetBranchAddress("d3D0Lund",d3D0Lund);
     fChain->SetBranchAddress("d4D0Index",d4D0Index);
     fChain->SetBranchAddress("d4D0Lund",d4D0Lund);
     fChain->SetBranchAddress("nChD",&nChD);
     fChain->SetBranchAddress("massChD",massChD);
     fChain->SetBranchAddress("pChD",pChD);
     fChain->SetBranchAddress("thChD",thChD);
     fChain->SetBranchAddress("phiChD",phiChD);
     fChain->SetBranchAddress("errMassChD",errMassChD);
     fChain->SetBranchAddress("m0ChD",m0ChD);
     fChain->SetBranchAddress("xChD",xChD);
     fChain->SetBranchAddress("yChD",yChD);
     fChain->SetBranchAddress("zChD",zChD);
     fChain->SetBranchAddress("s2xChD",s2xChD);
     fChain->SetBranchAddress("s2yChD",s2yChD);
     fChain->SetBranchAddress("s2zChD",s2zChD);
     fChain->SetBranchAddress("chi2ChD",chi2ChD);
     fChain->SetBranchAddress("dofChD",dofChD);
     fChain->SetBranchAddress("stChD",stChD);
     fChain->SetBranchAddress("ndauChD",ndauChD);
     if(isMC || mcasdata)   fChain->SetBranchAddress("MCChD",MCChD);
     fChain->SetBranchAddress("d1ChDIndex",d1ChDIndex);
     fChain->SetBranchAddress("d1ChDLund",d1ChDLund);
     fChain->SetBranchAddress("d2ChDIndex",d2ChDIndex);
     fChain->SetBranchAddress("d2ChDLund",d2ChDLund);
     fChain->SetBranchAddress("d3ChDIndex",d3ChDIndex);
     fChain->SetBranchAddress("d3ChDLund",d3ChDLund);
     fChain->SetBranchAddress("d4ChDIndex",d4ChDIndex);
     fChain->SetBranchAddress("d4ChDLund",d4ChDLund);
     fChain->SetBranchAddress("nKs",&nKs);
     fChain->SetBranchAddress("massKs",massKs);
     fChain->SetBranchAddress("pKs",pKs);
     fChain->SetBranchAddress("thKs",thKs);
     fChain->SetBranchAddress("phiKs",phiKs);
     fChain->SetBranchAddress("errMassKs",errMassKs);
     fChain->SetBranchAddress("m0Ks",m0Ks);
     fChain->SetBranchAddress("xKs",xKs);
     fChain->SetBranchAddress("yKs",yKs);
     fChain->SetBranchAddress("zKs",zKs);
     fChain->SetBranchAddress("s2xKs",s2xKs);
     fChain->SetBranchAddress("s2yKs",s2yKs);
     fChain->SetBranchAddress("s2zKs",s2zKs);
     fChain->SetBranchAddress("chi2Ks",chi2Ks);
     fChain->SetBranchAddress("dofKs",dofKs);
     fChain->SetBranchAddress("stKs",stKs);
     fChain->SetBranchAddress("ndauKs",ndauKs);
     if(isMC || mcasdata)   fChain->SetBranchAddress("MCKs",MCKs);
     fChain->SetBranchAddress("d1KsIndex",d1KsIndex);
     fChain->SetBranchAddress("d1KsLund",d1KsLund);
     fChain->SetBranchAddress("d2KsIndex",d2KsIndex);
     fChain->SetBranchAddress("d2KsLund",d2KsLund);
     fChain->SetBranchAddress("nPi0",&nPi0);
     fChain->SetBranchAddress("massPi0",massPi0);
     fChain->SetBranchAddress("pPi0",pPi0);
     fChain->SetBranchAddress("thPi0",thPi0);
     fChain->SetBranchAddress("phiPi0",phiPi0);
     fChain->SetBranchAddress("errMassPi0",errMassPi0);
     fChain->SetBranchAddress("m0Pi0",m0Pi0);
     fChain->SetBranchAddress("xPi0",xPi0);
     fChain->SetBranchAddress("yPi0",yPi0);
     fChain->SetBranchAddress("zPi0",zPi0);
     fChain->SetBranchAddress("s2xPi0",s2xPi0);
     fChain->SetBranchAddress("s2yPi0",s2yPi0);
     fChain->SetBranchAddress("s2zPi0",s2zPi0);
     fChain->SetBranchAddress("chi2Pi0",chi2Pi0);
     fChain->SetBranchAddress("dofPi0",dofPi0);
     fChain->SetBranchAddress("stPi0",stPi0);
     fChain->SetBranchAddress("ndauPi0",ndauPi0);
     if(isMC || mcasdata)   fChain->SetBranchAddress("MCPi0",MCPi0);
     fChain->SetBranchAddress("d1Pi0Index",d1Pi0Index);
     fChain->SetBranchAddress("d1Pi0Lund",d1Pi0Lund);
     fChain->SetBranchAddress("d2Pi0Index",d2Pi0Index);
     fChain->SetBranchAddress("d2Pi0Lund",d2Pi0Lund);
   }

   fChain->SetBranchAddress("nGConv",&nGConv);
   fChain->SetBranchAddress("massGConv",massGConv);
   fChain->SetBranchAddress("pGConv",pGConv);
   fChain->SetBranchAddress("thGConv",thGConv);
   fChain->SetBranchAddress("phiGConv",phiGConv);
   fChain->SetBranchAddress("errMassGConv",errMassGConv);
   fChain->SetBranchAddress("m0GConv",m0GConv);
   fChain->SetBranchAddress("xGConv",xGConv);
   fChain->SetBranchAddress("yGConv",yGConv);
   fChain->SetBranchAddress("zGConv",zGConv);
   fChain->SetBranchAddress("s2xGConv",s2xGConv);
   fChain->SetBranchAddress("s2yGConv",s2yGConv);
   fChain->SetBranchAddress("s2zGConv",s2zGConv);
   fChain->SetBranchAddress("chi2GConv",chi2GConv);
   fChain->SetBranchAddress("dofGConv",dofGConv);
   fChain->SetBranchAddress("stGConv",stGConv);
   fChain->SetBranchAddress("ndauGConv",ndauGConv);
   if(isMC || mcasdata)   fChain->SetBranchAddress("MCGConv",MCGConv);
   fChain->SetBranchAddress("d1GConvIndex",d1GConvIndex);
   fChain->SetBranchAddress("d1GConvLund",d1GConvLund);
   fChain->SetBranchAddress("d2GConvIndex",d2GConvIndex);
   fChain->SetBranchAddress("d2GConvLund",d2GConvLund);

   if (fNewFormat) {
     fChain->SetBranchAddress("nDalitz",&nDalitz);
     fChain->SetBranchAddress("massDalitz",massDalitz);
     fChain->SetBranchAddress("pDalitz",pDalitz);
     fChain->SetBranchAddress("thDalitz",thDalitz);
     fChain->SetBranchAddress("phiDalitz",phiDalitz);
     fChain->SetBranchAddress("errMassDalitz",errMassDalitz);
     fChain->SetBranchAddress("m0Dalitz",m0Dalitz);
     fChain->SetBranchAddress("xDalitz",xDalitz);
     fChain->SetBranchAddress("yDalitz",yDalitz);
     fChain->SetBranchAddress("zDalitz",zDalitz);
     fChain->SetBranchAddress("s2xDalitz",s2xDalitz);
     fChain->SetBranchAddress("s2yDalitz",s2yDalitz);
     fChain->SetBranchAddress("s2zDalitz",s2zDalitz);
     fChain->SetBranchAddress("chi2Dalitz",chi2Dalitz);
     fChain->SetBranchAddress("dofDalitz",dofDalitz);
     fChain->SetBranchAddress("stDalitz",stDalitz);
     fChain->SetBranchAddress("ndauDalitz",ndauDalitz);
     if(isMC || mcasdata)   fChain->SetBranchAddress("MCDalitz",MCDalitz);
     fChain->SetBranchAddress("d1DalitzIndex",d1DalitzIndex);
     fChain->SetBranchAddress("d1DalitzLund",d1DalitzLund);
     fChain->SetBranchAddress("d2DalitzIndex",d2DalitzIndex);
     fChain->SetBranchAddress("d2DalitzLund",d2DalitzLund);
     
     fChain->SetBranchAddress("nJpsi",&nJpsi);
     fChain->SetBranchAddress("massJpsi",massJpsi);
     fChain->SetBranchAddress("pJpsi",pJpsi);
     fChain->SetBranchAddress("thJpsi",thJpsi);
     fChain->SetBranchAddress("phiJpsi",phiJpsi);
     fChain->SetBranchAddress("errMassJpsi",errMassJpsi);
     fChain->SetBranchAddress("m0Jpsi",m0Jpsi);
     fChain->SetBranchAddress("xJpsi",xJpsi);
     fChain->SetBranchAddress("yJpsi",yJpsi);
     fChain->SetBranchAddress("zJpsi",zJpsi);
     fChain->SetBranchAddress("s2xJpsi",s2xJpsi);
     fChain->SetBranchAddress("s2yJpsi",s2yJpsi);
     fChain->SetBranchAddress("s2zJpsi",s2zJpsi);
     fChain->SetBranchAddress("chi2Jpsi",chi2Jpsi);
     fChain->SetBranchAddress("dofJpsi",dofJpsi);
     fChain->SetBranchAddress("stJpsi",stJpsi);
     fChain->SetBranchAddress("ndauJpsi",ndauJpsi);
     if(isMC || mcasdata)   fChain->SetBranchAddress("MCJpsi",MCJpsi);
     fChain->SetBranchAddress("d1JpsiIndex",d1JpsiIndex);
     fChain->SetBranchAddress("d1JpsiLund",d1JpsiLund);
     fChain->SetBranchAddress("d2JpsiIndex",d2JpsiIndex);
     fChain->SetBranchAddress("d2JpsiLund",d2JpsiLund);   
   }

   fChain->SetBranchAddress("nTrk",&nTrk);

   if(READALL==1) {   
     fChain->SetBranchAddress("IfrLayTrk",IfrLayTrk);
     fChain->SetBranchAddress("IfrNsTrk",IfrNsTrk);
     fChain->SetBranchAddress("IfrInnerTrk",IfrInnerTrk);
     fChain->SetBranchAddress("IfrBarrelTrk",IfrBarrelTrk);
     fChain->SetBranchAddress("IfrFWDTrk",IfrFWDTrk);
     fChain->SetBranchAddress("IfrBWDTrk",IfrBWDTrk);
     fChain->SetBranchAddress("IfrMeasIntLenTrk",IfrMeasIntLenTrk);
     fChain->SetBranchAddress("IfrFirstHitTrk",IfrFirstHitTrk);
     fChain->SetBranchAddress("IfrLastHitTrk",IfrLastHitTrk);
     fChain->SetBranchAddress("lMomTrk",lMomTrk);
     fChain->SetBranchAddress("ZMom42Trk",ZMom42Trk);
     fChain->SetBranchAddress("ecalTrk",ecalTrk);
     fChain->SetBranchAddress("ecalXTrk",ecalXTrk);
     fChain->SetBranchAddress("ecalYTrk",ecalYTrk);
     fChain->SetBranchAddress("ecalZTrk",ecalZTrk);
     fChain->SetBranchAddress("nCryTrk",nCryTrk);
     fChain->SetBranchAddress("nBumpTrk",nBumpTrk);
     fChain->SetBranchAddress("ZMom20Trk",ZMom20Trk);
     fChain->SetBranchAddress("secMomTrk",secMomTrk);
     fChain->SetBranchAddress("s1s9Trk",s1s9Trk);
     fChain->SetBranchAddress("s9s25Trk",s9s25Trk);
     fChain->SetBranchAddress("erawTrk",erawTrk);
     fChain->SetBranchAddress("phiClusterTrk",phiClusterTrk);
     fChain->SetBranchAddress("thetaClusterTrk",thetaClusterTrk);
     fChain->SetBranchAddress("covEETrk",covEETrk);
     fChain->SetBranchAddress("covTTTrk",covTTTrk);
     fChain->SetBranchAddress("covPPTrk",covPPTrk);
     fChain->SetBranchAddress("covRRTrk",covRRTrk);
     fChain->SetBranchAddress("phicMatTrk",phicMatTrk);
     fChain->SetBranchAddress("trkcMatTrk",trkcMatTrk);
     fChain->SetBranchAddress("nPidTrk",nPidTrk);
     fChain->SetBranchAddress("emcStatusTrk",emcStatusTrk);
     fChain->SetBranchAddress("phiAtEMCTrk",phiAtEMCTrk);
     fChain->SetBranchAddress("thetaAtEMCTrk",thetaAtEMCTrk);
     fChain->SetBranchAddress("isvtTrk",isvtTrk);
     fChain->SetBranchAddress("nsvtTrk",nsvtTrk);
     fChain->SetBranchAddress("fhitTrk",fhitTrk);
     fChain->SetBranchAddress("ndchTrk",ndchTrk);
     fChain->SetBranchAddress("lhitTrk",lhitTrk);
     fChain->SetBranchAddress("tLenTrk",tLenTrk);
     fChain->SetBranchAddress("ntdofTrk",ntdofTrk);
     fChain->SetBranchAddress("tproTrk",tproTrk);
     fChain->SetBranchAddress("tChi2Trk",tChi2Trk);
     fChain->SetBranchAddress("cPidTrk",cPidTrk);
     fChain->SetBranchAddress("sfRangeTrk",sfRangeTrk);
     fChain->SetBranchAddress("trkFitStatusTrk",trkFitStatusTrk);
   }
   
   fChain->SetBranchAddress("chargeTrk",chargeTrk);
   fChain->SetBranchAddress("momentumTrk",momentumTrk);
   
   if(READALL==1) {   
     fChain->SetBranchAddress("ppcov00",ppcov00);
     fChain->SetBranchAddress("ppcov10",ppcov10);
     fChain->SetBranchAddress("ppcov11",ppcov11);
     fChain->SetBranchAddress("ppcov20",ppcov20);
     fChain->SetBranchAddress("ppcov21",ppcov21);
     fChain->SetBranchAddress("ppcov22",ppcov22);
   }
   
   fChain->SetBranchAddress("xPocaTrk",xPocaTrk);
   fChain->SetBranchAddress("yPocaTrk",yPocaTrk);
   fChain->SetBranchAddress("zPocaTrk",zPocaTrk);
   fChain->SetBranchAddress("thetaTrk",thetaTrk);
   fChain->SetBranchAddress("phiTrk",phiTrk);
   fChain->SetBranchAddress("muonIdTrk",muonIdTrk);
   fChain->SetBranchAddress("elecIdTrk",elecIdTrk);
   fChain->SetBranchAddress("kaonIdTrk",kaonIdTrk);
   fChain->SetBranchAddress("pionIdTrk",pionIdTrk);
   if (isMC || mcasdata) {
     fChain->SetBranchAddress("idTrk",idTrk);
     fChain->SetBranchAddress("IndexTrk",IndexTrk);
     fChain->SetBranchAddress("IndexNtTrk",IndexNtTrk);
   }
   fChain->SetBranchAddress("B0RecTrk",B0RecTrk);
   fChain->SetBranchAddress("chBRecTrk",chBRecTrk);

   if(READALL==1) {   
     fChain->SetBranchAddress("nGam",&nGam);
     fChain->SetBranchAddress("IfrLayGam",IfrLayGam);
     fChain->SetBranchAddress("IfrNsGam",IfrNsGam);
     fChain->SetBranchAddress("IfrInnerGam",IfrInnerGam);
     fChain->SetBranchAddress("IfrBarrelGam",IfrBarrelGam);
     fChain->SetBranchAddress("IfrFWDGam",IfrFWDGam);
     fChain->SetBranchAddress("IfrBWDGam",IfrBWDGam);
     fChain->SetBranchAddress("IfrMeasIntLenGam",IfrMeasIntLenGam);
     fChain->SetBranchAddress("IfrFirstHitGam",IfrFirstHitGam);
     fChain->SetBranchAddress("IfrLastHitGam",IfrLastHitGam);
     fChain->SetBranchAddress("IfrExpIntLenGam",IfrExpIntLenGam);
     fChain->SetBranchAddress("IfrIntLenBeforeIronGam",IfrIntLenBeforeIronGam);
     fChain->SetBranchAddress("IfrTrkMatchGam",IfrTrkMatchGam);
     fChain->SetBranchAddress("IfrEmcMatchGam",IfrEmcMatchGam);
     fChain->SetBranchAddress("IfrLastBarrelGam",IfrLastBarrelGam);
     fChain->SetBranchAddress("IfrCLFitChi2Gam",IfrCLFitChi2Gam);
     fChain->SetBranchAddress("IfrStrips0",IfrStrips0);
     fChain->SetBranchAddress("IfrStrips1",IfrStrips1);
     fChain->SetBranchAddress("IfrStrips2",IfrStrips2);
     fChain->SetBranchAddress("IfrStrips3",IfrStrips3);
     fChain->SetBranchAddress("IfrStrips4",IfrStrips4);
     fChain->SetBranchAddress("IfrStrips5",IfrStrips5);
     fChain->SetBranchAddress("IfrStrips6",IfrStrips6);
     fChain->SetBranchAddress("IfrStrips7",IfrStrips7);
     fChain->SetBranchAddress("IfrStrips8",IfrStrips8);
     fChain->SetBranchAddress("IfrStrips9",IfrStrips9);
     fChain->SetBranchAddress("IfrStrips10",IfrStrips10);
     fChain->SetBranchAddress("IfrStrips11",IfrStrips11);
     fChain->SetBranchAddress("IfrStrips12",IfrStrips12);
     fChain->SetBranchAddress("IfrStrips13",IfrStrips13);
     fChain->SetBranchAddress("IfrStrips14",IfrStrips14);
     fChain->SetBranchAddress("IfrStrips15",IfrStrips15);
     fChain->SetBranchAddress("IfrStrips16",IfrStrips16);
     fChain->SetBranchAddress("IfrStrips17",IfrStrips17);
     fChain->SetBranchAddress("IfrStrips18",IfrStrips18);
     fChain->SetBranchAddress("IfrStrips19",IfrStrips19);
     fChain->SetBranchAddress("lMomGam",lMomGam);
     fChain->SetBranchAddress("ZMom42Gam",ZMom42Gam);
     fChain->SetBranchAddress("ecalGam",ecalGam);
     fChain->SetBranchAddress("ecalXGam",ecalXGam);
     fChain->SetBranchAddress("ecalYGam",ecalYGam);
     fChain->SetBranchAddress("ecalZGam",ecalZGam);
     fChain->SetBranchAddress("nCryGam",nCryGam);
     fChain->SetBranchAddress("nBumpGam",nBumpGam);
     fChain->SetBranchAddress("ZMom20Gam",ZMom20Gam);
     fChain->SetBranchAddress("secMomGam",secMomGam);
     fChain->SetBranchAddress("s1s9Gam",s1s9Gam);
     fChain->SetBranchAddress("s9s25Gam",s9s25Gam);
     fChain->SetBranchAddress("erawGam",erawGam);
     fChain->SetBranchAddress("phiClusterGam",phiClusterGam);
     fChain->SetBranchAddress("thetaClusterGam",thetaClusterGam);
     fChain->SetBranchAddress("covEEGam",covEEGam);
     fChain->SetBranchAddress("covTTGam",covTTGam);
     fChain->SetBranchAddress("covPPGam",covPPGam);
     fChain->SetBranchAddress("covRRGam",covRRGam);
     fChain->SetBranchAddress("emcStatusGam",emcStatusGam);
     fChain->SetBranchAddress("thetaGam",thetaGam);
     fChain->SetBranchAddress("phiGam",phiGam);
     fChain->SetBranchAddress("energyGam",energyGam);
     if(isMC || mcasdata){
       fChain->SetBranchAddress("idGam",idGam);
       fChain->SetBranchAddress("IndexGam",IndexGam);
       fChain->SetBranchAddress("IndexNtGam",IndexNtGam);
     }
     fChain->SetBranchAddress("B0RecGam",B0RecGam);
     fChain->SetBranchAddress("chBRecGam",chBRecGam);
   }

   Notify();
}

// ----------------------------------------------------------------------
Bool_t recoilBase::Notify() {
//   called when loading a new file
//   get branch pointers
   b_event = fChain->GetBranch("event");
   b_runNumber = fChain->GetBranch("runNumber");
   b_platform = fChain->GetBranch("platform");
   b_partition = fChain->GetBranch("partition");
   b_upperID = fChain->GetBranch("upperID");
   b_lowerID = fChain->GetBranch("lowerID");
   b_beamSX = fChain->GetBranch("beamSX");
   b_beamSY = fChain->GetBranch("beamSY");
   b_beamSZ = fChain->GetBranch("beamSZ");
   b_beamSCovXX = fChain->GetBranch("beamSCovXX");
   b_beamSCovYY = fChain->GetBranch("beamSCovYY");
   b_beamSCovZZ = fChain->GetBranch("beamSCovZZ");
   b_beamSCovXZ = fChain->GetBranch("beamSCovXZ");
   b_pxUps = fChain->GetBranch("pxUps");
   b_pyUps = fChain->GetBranch("pyUps");
   b_pzUps = fChain->GetBranch("pzUps");
   b_eUps = fChain->GetBranch("eUps");
   b_nTrkTot = fChain->GetBranch("nTrkTot");
   b_W2 = fChain->GetBranch("W2");
   b_FoxWol2 = fChain->GetBranch("FoxWol2");
   b_FoxWol2Neu = fChain->GetBranch("FoxWol2Neu");
   b_thrust = fChain->GetBranch("thrust");
   b_thrustNeu = fChain->GetBranch("thrustNeu");
   b_nMc = fChain->GetBranch("nMc");
   b_pMc = fChain->GetBranch("pMc");
   b_massMc = fChain->GetBranch("massMc");
   b_thetaMc = fChain->GetBranch("thetaMc");
   b_phiMc = fChain->GetBranch("phiMc");
   b_idMc = fChain->GetBranch("idMc");
   b_mothMc = fChain->GetBranch("mothMc");
   b_nDauMc = fChain->GetBranch("nDauMc");
   b_xMc = fChain->GetBranch("xMc");
   b_yMc = fChain->GetBranch("yMc");
   b_zMc = fChain->GetBranch("zMc");
   b_nB0 = fChain->GetBranch("nB0");
   b_massB0 = fChain->GetBranch("massB0");
   b_pB0 = fChain->GetBranch("pB0");
   b_thB0 = fChain->GetBranch("thB0");
   b_phiB0 = fChain->GetBranch("phiB0");
   b_errMassB0 = fChain->GetBranch("errMassB0");
   b_m0B0 = fChain->GetBranch("m0B0");
   b_xB0 = fChain->GetBranch("xB0");
   b_yB0 = fChain->GetBranch("yB0");
   b_zB0 = fChain->GetBranch("zB0");
   b_s2xB0 = fChain->GetBranch("s2xB0");
   b_s2yB0 = fChain->GetBranch("s2yB0");
   b_s2zB0 = fChain->GetBranch("s2zB0");
   b_chi2B0 = fChain->GetBranch("chi2B0");
   b_dofB0 = fChain->GetBranch("dofB0");
   b_stB0 = fChain->GetBranch("stB0");
   b_ndauB0 = fChain->GetBranch("ndauB0");
   b_MCB0 = fChain->GetBranch("MCB0");
   b_mseB0 = fChain->GetBranch("mseB0");
   b_mHatB0 = fChain->GetBranch("mHatB0");
   b_deltaeB0 = fChain->GetBranch("deltaeB0");
   b_ThruB0 = fChain->GetBranch("ThruB0");
   b_thThruB0 = fChain->GetBranch("thThruB0");
   b_phiThruB0 = fChain->GetBranch("phiThruB0");
   b_cosTBB0 = fChain->GetBranch("cosTBB0");
   b_d1B0Index = fChain->GetBranch("d1B0Index");
   b_d1B0Lund = fChain->GetBranch("d1B0Lund");
   b_d2B0Index = fChain->GetBranch("d2B0Index");
   b_d2B0Lund = fChain->GetBranch("d2B0Lund");
   b_d3B0Index = fChain->GetBranch("d3B0Index");
   b_d3B0Lund = fChain->GetBranch("d3B0Lund");
   b_d4B0Index = fChain->GetBranch("d4B0Index");
   b_d4B0Lund = fChain->GetBranch("d4B0Lund");
   b_d5B0Index = fChain->GetBranch("d5B0Index");
   b_d5B0Lund = fChain->GetBranch("d5B0Lund");
   b_d6B0Index = fChain->GetBranch("d6B0Index");
   b_d6B0Lund = fChain->GetBranch("d6B0Lund");
   b_d7B0Index = fChain->GetBranch("d7B0Index");
   b_d7B0Lund = fChain->GetBranch("d7B0Lund");
   b_modeB0 = fChain->GetBranch("modeB0");
   b_purB0 = fChain->GetBranch("purB0");
   b_intpurB0 = fChain->GetBranch("intpurB0");
   b_VtxXLepB0 = fChain->GetBranch("VtxXLepB0");
   b_VtxYLepB0 = fChain->GetBranch("VtxYLepB0");
   b_VtxZLepB0 = fChain->GetBranch("VtxZLepB0");
   b_VtxCovXXLepB0 = fChain->GetBranch("VtxCovXXLepB0");
   b_VtxCovYYLepB0 = fChain->GetBranch("VtxCovYYLepB0");
   b_VtxCovXYLepB0 = fChain->GetBranch("VtxCovXYLepB0");
   b_VtxCovZZLepB0 = fChain->GetBranch("VtxCovZZLepB0");
   b_VtxCovXZLepB0 = fChain->GetBranch("VtxCovXZLepB0");
   b_VtxCovYZLepB0 = fChain->GetBranch("VtxCovYZLepB0");
   b_VtxChiSqLepB0 = fChain->GetBranch("VtxChiSqLepB0");
   b_VtxNDofLepB0 = fChain->GetBranch("VtxNDofLepB0");
   b_VtxStatLepB0 = fChain->GetBranch("VtxStatLepB0");
   b_VtxNUsedLepB0 = fChain->GetBranch("VtxNUsedLepB0");
   b_DocaLepB0 = fChain->GetBranch("DocaLepB0");
   b_DocaErrLepB0 = fChain->GetBranch("DocaErrLepB0");
   b_VtxXXB0 = fChain->GetBranch("VtxXXB0");
   b_VtxYXB0 = fChain->GetBranch("VtxYXB0");
   b_VtxZXB0 = fChain->GetBranch("VtxZXB0");
   b_VtxCovXXXB0 = fChain->GetBranch("VtxCovXXXB0");
   b_VtxCovYYXB0 = fChain->GetBranch("VtxCovYYXB0");
   b_VtxCovXYXB0 = fChain->GetBranch("VtxCovXYXB0");
   b_VtxCovZZXB0 = fChain->GetBranch("VtxCovZZXB0");
   b_VtxCovXZXB0 = fChain->GetBranch("VtxCovXZXB0");
   b_VtxCovYZXB0 = fChain->GetBranch("VtxCovYZXB0");
   b_VtxChiSqXB0 = fChain->GetBranch("VtxChiSqXB0");
   b_VtxNDofXB0 = fChain->GetBranch("VtxNDofXB0");
   b_VtxStatXB0 = fChain->GetBranch("VtxStatXB0");
   b_VtxNUsedXB0 = fChain->GetBranch("VtxNUsedXB0");
   b_VtxPXB0 = fChain->GetBranch("VtxPXB0");
   b_VtxPhiXB0 = fChain->GetBranch("VtxPhiXB0");
   b_VtxThetaXB0 = fChain->GetBranch("VtxThetaXB0");
   b_ThrustXB0 = fChain->GetBranch("ThrustXB0");
   b_ThrustXPhiB0 = fChain->GetBranch("ThrustXPhiB0");
   b_ThrustXThetaB0 = fChain->GetBranch("ThrustXThetaB0");
   b_MassPB0 = fChain->GetBranch("MassPB0");
   b_MassPhiB0 = fChain->GetBranch("MassPhiB0");
   b_MassThetaB0 = fChain->GetBranch("MassThetaB0");
   b_Cov00B0 = fChain->GetBranch("Cov00B0");
   b_Cov10B0 = fChain->GetBranch("Cov10B0");
   b_Cov11B0 = fChain->GetBranch("Cov11B0");
   b_Cov20B0 = fChain->GetBranch("Cov20B0");
   b_Cov21B0 = fChain->GetBranch("Cov21B0");
   b_Cov22B0 = fChain->GetBranch("Cov22B0");
   b_Cov30B0 = fChain->GetBranch("Cov30B0");
   b_Cov31B0 = fChain->GetBranch("Cov31B0");
   b_Cov32B0 = fChain->GetBranch("Cov32B0");
   b_Cov33B0 = fChain->GetBranch("Cov33B0");
   b_nChB = fChain->GetBranch("nChB");
   b_massChB = fChain->GetBranch("massChB");
   b_pChB = fChain->GetBranch("pChB");
   b_thChB = fChain->GetBranch("thChB");
   b_phiChB = fChain->GetBranch("phiChB");
   b_errMassChB = fChain->GetBranch("errMassChB");
   b_m0ChB = fChain->GetBranch("m0ChB");
   b_xChB = fChain->GetBranch("xChB");
   b_yChB = fChain->GetBranch("yChB");
   b_zChB = fChain->GetBranch("zChB");
   b_s2xChB = fChain->GetBranch("s2xChB");
   b_s2yChB = fChain->GetBranch("s2yChB");
   b_s2zChB = fChain->GetBranch("s2zChB");
   b_chi2ChB = fChain->GetBranch("chi2ChB");
   b_dofChB = fChain->GetBranch("dofChB");
   b_stChB = fChain->GetBranch("stChB");
   b_ndauChB = fChain->GetBranch("ndauChB");
   b_MCChB = fChain->GetBranch("MCChB");
   b_mseChB = fChain->GetBranch("mseChB");
   b_mHatChB = fChain->GetBranch("mHatChB");
   b_deltaeChB = fChain->GetBranch("deltaeChB");
   b_ThruChB = fChain->GetBranch("ThruChB");
   b_thThruChB = fChain->GetBranch("thThruChB");
   b_phiThruChB = fChain->GetBranch("phiThruChB");
   b_cosTBChB = fChain->GetBranch("cosTBChB");
   b_d1ChBIndex = fChain->GetBranch("d1ChBIndex");
   b_d1ChBLund = fChain->GetBranch("d1ChBLund");
   b_d2ChBIndex = fChain->GetBranch("d2ChBIndex");
   b_d2ChBLund = fChain->GetBranch("d2ChBLund");
   b_d3ChBIndex = fChain->GetBranch("d3ChBIndex");
   b_d3ChBLund = fChain->GetBranch("d3ChBLund");
   b_d4ChBIndex = fChain->GetBranch("d4ChBIndex");
   b_d4ChBLund = fChain->GetBranch("d4ChBLund");
   b_d5ChBIndex = fChain->GetBranch("d5ChBIndex");
   b_d5ChBLund = fChain->GetBranch("d5ChBLund");
   b_d6ChBIndex = fChain->GetBranch("d6ChBIndex");
   b_d6ChBLund = fChain->GetBranch("d6ChBLund");
   b_d7ChBIndex = fChain->GetBranch("d7ChBIndex");
   b_d7ChBLund = fChain->GetBranch("d7ChBLund");
   b_modeChB = fChain->GetBranch("modeChB");
   b_purChB = fChain->GetBranch("purChB");
   b_intpurChB = fChain->GetBranch("intpurChB");
   b_VtxXLepChB = fChain->GetBranch("VtxXLepChB");
   b_VtxYLepChB = fChain->GetBranch("VtxYLepChB");
   b_VtxZLepChB = fChain->GetBranch("VtxZLepChB");
   b_VtxCovXXLepChB = fChain->GetBranch("VtxCovXXLepChB");
   b_VtxCovYYLepChB = fChain->GetBranch("VtxCovYYLepChB");
   b_VtxCovXYLepChB = fChain->GetBranch("VtxCovXYLepChB");
   b_VtxCovZZLepChB = fChain->GetBranch("VtxCovZZLepChB");
   b_VtxCovXZLepChB = fChain->GetBranch("VtxCovXZLepChB");
   b_VtxCovYZLepChB = fChain->GetBranch("VtxCovYZLepChB");
   b_VtxChiSqLepChB = fChain->GetBranch("VtxChiSqLepChB");
   b_VtxNDofLepChB = fChain->GetBranch("VtxNDofLepChB");
   b_VtxStatLepChB = fChain->GetBranch("VtxStatLepChB");
   b_VtxNUsedLepChB = fChain->GetBranch("VtxNUsedLepChB");
   b_DocaLepChB = fChain->GetBranch("DocaLepChB");
   b_DocaErrLepChB = fChain->GetBranch("DocaErrLepChB");
   b_VtxXXChB = fChain->GetBranch("VtxXXChB");
   b_VtxYXChB = fChain->GetBranch("VtxYXChB");
   b_VtxZXChB = fChain->GetBranch("VtxZXChB");
   b_VtxCovXXXChB = fChain->GetBranch("VtxCovXXXChB");
   b_VtxCovYYXChB = fChain->GetBranch("VtxCovYYXChB");
   b_VtxCovXYXChB = fChain->GetBranch("VtxCovXYXChB");
   b_VtxCovZZXChB = fChain->GetBranch("VtxCovZZXChB");
   b_VtxCovXZXChB = fChain->GetBranch("VtxCovXZXChB");
   b_VtxCovYZXChB = fChain->GetBranch("VtxCovYZXChB");
   b_VtxChiSqXChB = fChain->GetBranch("VtxChiSqXChB");
   b_VtxNDofXChB = fChain->GetBranch("VtxNDofXChB");
   b_VtxStatXChB = fChain->GetBranch("VtxStatXChB");
   b_VtxNUsedXChB = fChain->GetBranch("VtxNUsedXChB");
   b_VtxPXChB = fChain->GetBranch("VtxPXChB");
   b_VtxPhiXChB = fChain->GetBranch("VtxPhiXChB");
   b_VtxThetaXChB = fChain->GetBranch("VtxThetaXChB");
   b_ThrustXChB = fChain->GetBranch("ThrustXChB");
   b_ThrustXPhiChB = fChain->GetBranch("ThrustXPhiChB");
   b_ThrustXThetaChB = fChain->GetBranch("ThrustXThetaChB");
   b_MassPChB = fChain->GetBranch("MassPChB");
   b_MassPhiChB = fChain->GetBranch("MassPhiChB");
   b_MassThetaChB = fChain->GetBranch("MassThetaChB");
   b_Cov00ChB = fChain->GetBranch("Cov00ChB");
   b_Cov10ChB = fChain->GetBranch("Cov10ChB");
   b_Cov11ChB = fChain->GetBranch("Cov11ChB");
   b_Cov20ChB = fChain->GetBranch("Cov20ChB");
   b_Cov21ChB = fChain->GetBranch("Cov21ChB");
   b_Cov22ChB = fChain->GetBranch("Cov22ChB");
   b_Cov30ChB = fChain->GetBranch("Cov30ChB");
   b_Cov31ChB = fChain->GetBranch("Cov31ChB");
   b_Cov32ChB = fChain->GetBranch("Cov32ChB");
   b_Cov33ChB = fChain->GetBranch("Cov33ChB");
   b_nDstar = fChain->GetBranch("nDstar");
   b_massDstar = fChain->GetBranch("massDstar");
   b_pDstar = fChain->GetBranch("pDstar");
   b_thDstar = fChain->GetBranch("thDstar");
   b_phiDstar = fChain->GetBranch("phiDstar");
   b_errMassDstar = fChain->GetBranch("errMassDstar");
   b_m0Dstar = fChain->GetBranch("m0Dstar");
   b_xDstar = fChain->GetBranch("xDstar");
   b_yDstar = fChain->GetBranch("yDstar");
   b_zDstar = fChain->GetBranch("zDstar");
   b_s2xDstar = fChain->GetBranch("s2xDstar");
   b_s2yDstar = fChain->GetBranch("s2yDstar");
   b_s2zDstar = fChain->GetBranch("s2zDstar");
   b_chi2Dstar = fChain->GetBranch("chi2Dstar");
   b_dofDstar = fChain->GetBranch("dofDstar");
   b_stDstar = fChain->GetBranch("stDstar");
   b_ndauDstar = fChain->GetBranch("ndauDstar");
   b_MCDstar = fChain->GetBranch("MCDstar");
   b_d1DstarIndex = fChain->GetBranch("d1DstarIndex");
   b_d1DstarLund = fChain->GetBranch("d1DstarLund");
   b_d2DstarIndex = fChain->GetBranch("d2DstarIndex");
   b_d2DstarLund = fChain->GetBranch("d2DstarLund");
   b_nDstarBS = fChain->GetBranch("nDstarBS");
   b_massDstarBS = fChain->GetBranch("massDstarBS");
   b_chi2DstarBS = fChain->GetBranch("chi2DstarBS");
   b_dofDstarBS = fChain->GetBranch("dofDstarBS");
   b_spixDstarBS = fChain->GetBranch("spixDstarBS");
   b_spiyDstarBS = fChain->GetBranch("spiyDstarBS");
   b_spizDstarBS = fChain->GetBranch("spizDstarBS");
   b_nDstar0 = fChain->GetBranch("nDstar0");
   b_massDstar0 = fChain->GetBranch("massDstar0");
   b_pDstar0 = fChain->GetBranch("pDstar0");
   b_thDstar0 = fChain->GetBranch("thDstar0");
   b_phiDstar0 = fChain->GetBranch("phiDstar0");
   b_errMassDstar0 = fChain->GetBranch("errMassDstar0");
   b_m0Dstar0 = fChain->GetBranch("m0Dstar0");
   b_xDstar0 = fChain->GetBranch("xDstar0");
   b_yDstar0 = fChain->GetBranch("yDstar0");
   b_zDstar0 = fChain->GetBranch("zDstar0");
   b_s2xDstar0 = fChain->GetBranch("s2xDstar0");
   b_s2yDstar0 = fChain->GetBranch("s2yDstar0");
   b_s2zDstar0 = fChain->GetBranch("s2zDstar0");
   b_chi2Dstar0 = fChain->GetBranch("chi2Dstar0");
   b_dofDstar0 = fChain->GetBranch("dofDstar0");
   b_stDstar0 = fChain->GetBranch("stDstar0");
   b_ndauDstar0 = fChain->GetBranch("ndauDstar0");
   b_MCDstar0 = fChain->GetBranch("MCDstar0");
   b_d1Dstar0Index = fChain->GetBranch("d1Dstar0Index");
   b_d1Dstar0Lund = fChain->GetBranch("d1Dstar0Lund");
   b_d2Dstar0Index = fChain->GetBranch("d2Dstar0Index");
   b_d2Dstar0Lund = fChain->GetBranch("d2Dstar0Lund");
   b_nD0 = fChain->GetBranch("nD0");
   b_massD0 = fChain->GetBranch("massD0");
   b_pD0 = fChain->GetBranch("pD0");
   b_thD0 = fChain->GetBranch("thD0");
   b_phiD0 = fChain->GetBranch("phiD0");
   b_errMassD0 = fChain->GetBranch("errMassD0");
   b_m0D0 = fChain->GetBranch("m0D0");
   b_xD0 = fChain->GetBranch("xD0");
   b_yD0 = fChain->GetBranch("yD0");
   b_zD0 = fChain->GetBranch("zD0");
   b_s2xD0 = fChain->GetBranch("s2xD0");
   b_s2yD0 = fChain->GetBranch("s2yD0");
   b_s2zD0 = fChain->GetBranch("s2zD0");
   b_chi2D0 = fChain->GetBranch("chi2D0");
   b_dofD0 = fChain->GetBranch("dofD0");
   b_stD0 = fChain->GetBranch("stD0");
   b_ndauD0 = fChain->GetBranch("ndauD0");
   b_MCD0 = fChain->GetBranch("MCD0");
   b_d1D0Index = fChain->GetBranch("d1D0Index");
   b_d1D0Lund = fChain->GetBranch("d1D0Lund");
   b_d2D0Index = fChain->GetBranch("d2D0Index");
   b_d2D0Lund = fChain->GetBranch("d2D0Lund");
   b_d3D0Index = fChain->GetBranch("d3D0Index");
   b_d3D0Lund = fChain->GetBranch("d3D0Lund");
   b_d4D0Index = fChain->GetBranch("d4D0Index");
   b_d4D0Lund = fChain->GetBranch("d4D0Lund");
   b_nChD = fChain->GetBranch("nChD");
   b_massChD = fChain->GetBranch("massChD");
   b_pChD = fChain->GetBranch("pChD");
   b_thChD = fChain->GetBranch("thChD");
   b_phiChD = fChain->GetBranch("phiChD");
   b_errMassChD = fChain->GetBranch("errMassChD");
   b_m0ChD = fChain->GetBranch("m0ChD");
   b_xChD = fChain->GetBranch("xChD");
   b_yChD = fChain->GetBranch("yChD");
   b_zChD = fChain->GetBranch("zChD");
   b_s2xChD = fChain->GetBranch("s2xChD");
   b_s2yChD = fChain->GetBranch("s2yChD");
   b_s2zChD = fChain->GetBranch("s2zChD");
   b_chi2ChD = fChain->GetBranch("chi2ChD");
   b_dofChD = fChain->GetBranch("dofChD");
   b_stChD = fChain->GetBranch("stChD");
   b_ndauChD = fChain->GetBranch("ndauChD");
   b_MCChD = fChain->GetBranch("MCChD");
   b_d1ChDIndex = fChain->GetBranch("d1ChDIndex");
   b_d1ChDLund = fChain->GetBranch("d1ChDLund");
   b_d2ChDIndex = fChain->GetBranch("d2ChDIndex");
   b_d2ChDLund = fChain->GetBranch("d2ChDLund");
   b_d3ChDIndex = fChain->GetBranch("d3ChDIndex");
   b_d3ChDLund = fChain->GetBranch("d3ChDLund");
   b_d4ChDIndex = fChain->GetBranch("d4ChDIndex");
   b_d4ChDLund = fChain->GetBranch("d4ChDLund");
   b_nKs = fChain->GetBranch("nKs");
   b_massKs = fChain->GetBranch("massKs");
   b_pKs = fChain->GetBranch("pKs");
   b_thKs = fChain->GetBranch("thKs");
   b_phiKs = fChain->GetBranch("phiKs");
   b_errMassKs = fChain->GetBranch("errMassKs");
   b_m0Ks = fChain->GetBranch("m0Ks");
   b_xKs = fChain->GetBranch("xKs");
   b_yKs = fChain->GetBranch("yKs");
   b_zKs = fChain->GetBranch("zKs");
   b_s2xKs = fChain->GetBranch("s2xKs");
   b_s2yKs = fChain->GetBranch("s2yKs");
   b_s2zKs = fChain->GetBranch("s2zKs");
   b_chi2Ks = fChain->GetBranch("chi2Ks");
   b_dofKs = fChain->GetBranch("dofKs");
   b_stKs = fChain->GetBranch("stKs");
   b_ndauKs = fChain->GetBranch("ndauKs");
   b_MCKs = fChain->GetBranch("MCKs");
   b_d1KsIndex = fChain->GetBranch("d1KsIndex");
   b_d1KsLund = fChain->GetBranch("d1KsLund");
   b_d2KsIndex = fChain->GetBranch("d2KsIndex");
   b_d2KsLund = fChain->GetBranch("d2KsLund");
   b_nPi0 = fChain->GetBranch("nPi0");
   b_massPi0 = fChain->GetBranch("massPi0");
   b_pPi0 = fChain->GetBranch("pPi0");
   b_thPi0 = fChain->GetBranch("thPi0");
   b_phiPi0 = fChain->GetBranch("phiPi0");
   b_errMassPi0 = fChain->GetBranch("errMassPi0");
   b_m0Pi0 = fChain->GetBranch("m0Pi0");
   b_xPi0 = fChain->GetBranch("xPi0");
   b_yPi0 = fChain->GetBranch("yPi0");
   b_zPi0 = fChain->GetBranch("zPi0");
   b_s2xPi0 = fChain->GetBranch("s2xPi0");
   b_s2yPi0 = fChain->GetBranch("s2yPi0");
   b_s2zPi0 = fChain->GetBranch("s2zPi0");
   b_chi2Pi0 = fChain->GetBranch("chi2Pi0");
   b_dofPi0 = fChain->GetBranch("dofPi0");
   b_stPi0 = fChain->GetBranch("stPi0");
   b_ndauPi0 = fChain->GetBranch("ndauPi0");
   b_MCPi0 = fChain->GetBranch("MCPi0");
   b_d1Pi0Index = fChain->GetBranch("d1Pi0Index");
   b_d1Pi0Lund = fChain->GetBranch("d1Pi0Lund");
   b_d2Pi0Index = fChain->GetBranch("d2Pi0Index");
   b_d2Pi0Lund = fChain->GetBranch("d2Pi0Lund");
   b_nGConv = fChain->GetBranch("nGConv");
   b_massGConv = fChain->GetBranch("massGConv");
   b_pGConv = fChain->GetBranch("pGConv");
   b_thGConv = fChain->GetBranch("thGConv");
   b_phiGConv = fChain->GetBranch("phiGConv");
   b_errMassGConv = fChain->GetBranch("errMassGConv");
   b_m0GConv = fChain->GetBranch("m0GConv");
   b_xGConv = fChain->GetBranch("xGConv");
   b_yGConv = fChain->GetBranch("yGConv");
   b_zGConv = fChain->GetBranch("zGConv");
   b_s2xGConv = fChain->GetBranch("s2xGConv");
   b_s2yGConv = fChain->GetBranch("s2yGConv");
   b_s2zGConv = fChain->GetBranch("s2zGConv");
   b_chi2GConv = fChain->GetBranch("chi2GConv");
   b_dofGConv = fChain->GetBranch("dofGConv");
   b_stGConv = fChain->GetBranch("stGConv");
   b_ndauGConv = fChain->GetBranch("ndauGConv");
   b_MCGConv = fChain->GetBranch("MCGConv");
   b_d1GConvIndex = fChain->GetBranch("d1GConvIndex");
   b_d1GConvLund = fChain->GetBranch("d1GConvLund");
   b_d2GConvIndex = fChain->GetBranch("d2GConvIndex");
   b_d2GConvLund = fChain->GetBranch("d2GConvLund");
   if (fNewFormat) {
     b_nDalitz = fChain->GetBranch("nDalitz");
     b_massDalitz = fChain->GetBranch("massDalitz");
     b_pDalitz = fChain->GetBranch("pDalitz");
     b_thDalitz = fChain->GetBranch("thDalitz");
     b_phiDalitz = fChain->GetBranch("phiDalitz");
     b_errMassDalitz = fChain->GetBranch("errMassDalitz");
     b_m0Dalitz = fChain->GetBranch("m0Dalitz");
     b_xDalitz = fChain->GetBranch("xDalitz");
     b_yDalitz = fChain->GetBranch("yDalitz");
     b_zDalitz = fChain->GetBranch("zDalitz");
     b_s2xDalitz = fChain->GetBranch("s2xDalitz");
     b_s2yDalitz = fChain->GetBranch("s2yDalitz");
     b_s2zDalitz = fChain->GetBranch("s2zDalitz");
     b_chi2Dalitz = fChain->GetBranch("chi2Dalitz");
     b_dofDalitz = fChain->GetBranch("dofDalitz");
     b_stDalitz = fChain->GetBranch("stDalitz");
     b_ndauDalitz = fChain->GetBranch("ndauDalitz");
     b_MCDalitz = fChain->GetBranch("MCDalitz");
     b_d1DalitzIndex = fChain->GetBranch("d1DalitzIndex");
     b_d1DalitzLund = fChain->GetBranch("d1DalitzLund");
     b_d2DalitzIndex = fChain->GetBranch("d2DalitzIndex");
     b_d2DalitzLund = fChain->GetBranch("d2DalitzLund");
     b_nJpsi = fChain->GetBranch("nJpsi");
     b_massJpsi = fChain->GetBranch("massJpsi");
     b_pJpsi = fChain->GetBranch("pJpsi");
     b_thJpsi = fChain->GetBranch("thJpsi");
     b_phiJpsi = fChain->GetBranch("phiJpsi");
     b_errMassJpsi = fChain->GetBranch("errMassJpsi");
     b_m0Jpsi = fChain->GetBranch("m0Jpsi");
     b_xJpsi = fChain->GetBranch("xJpsi");
     b_yJpsi = fChain->GetBranch("yJpsi");
     b_zJpsi = fChain->GetBranch("zJpsi");
     b_s2xJpsi = fChain->GetBranch("s2xJpsi");
     b_s2yJpsi = fChain->GetBranch("s2yJpsi");
     b_s2zJpsi = fChain->GetBranch("s2zJpsi");
     b_chi2Jpsi = fChain->GetBranch("chi2Jpsi");
     b_dofJpsi = fChain->GetBranch("dofJpsi");
     b_stJpsi = fChain->GetBranch("stJpsi");
     b_ndauJpsi = fChain->GetBranch("ndauJpsi");
     b_MCJpsi = fChain->GetBranch("MCJpsi");
     b_d1JpsiIndex = fChain->GetBranch("d1JpsiIndex");
     b_d1JpsiLund = fChain->GetBranch("d1JpsiLund");
     b_d2JpsiIndex = fChain->GetBranch("d2JpsiIndex");
     b_d2JpsiLund = fChain->GetBranch("d2JpsiLund");
   }
   b_nTrk = fChain->GetBranch("nTrk");
   b_IfrLayTrk = fChain->GetBranch("IfrLayTrk");
   b_IfrNsTrk = fChain->GetBranch("IfrNsTrk");
   b_IfrInnerTrk = fChain->GetBranch("IfrInnerTrk");
   b_IfrBarrelTrk = fChain->GetBranch("IfrBarrelTrk");
   b_IfrFWDTrk = fChain->GetBranch("IfrFWDTrk");
   b_IfrBWDTrk = fChain->GetBranch("IfrBWDTrk");
   b_IfrMeasIntLenTrk = fChain->GetBranch("IfrMeasIntLenTrk");
   b_IfrFirstHitTrk = fChain->GetBranch("IfrFirstHitTrk");
   b_IfrLastHitTrk = fChain->GetBranch("IfrLastHitTrk");
   b_lMomTrk = fChain->GetBranch("lMomTrk");
   b_ZMom42Trk = fChain->GetBranch("ZMom42Trk");
   b_ecalTrk = fChain->GetBranch("ecalTrk");
   b_ecalXTrk = fChain->GetBranch("ecalXTrk");
   b_ecalYTrk = fChain->GetBranch("ecalYTrk");
   b_ecalZTrk = fChain->GetBranch("ecalZTrk");
   b_nCryTrk = fChain->GetBranch("nCryTrk");
   b_nBumpTrk = fChain->GetBranch("nBumpTrk");
   b_ZMom20Trk = fChain->GetBranch("ZMom20Trk");
   b_secMomTrk = fChain->GetBranch("secMomTrk");
   b_s1s9Trk = fChain->GetBranch("s1s9Trk");
   b_s9s25Trk = fChain->GetBranch("s9s25Trk");
   b_erawTrk = fChain->GetBranch("erawTrk");
   b_phiClusterTrk = fChain->GetBranch("phiClusterTrk");
   b_thetaClusterTrk = fChain->GetBranch("thetaClusterTrk");
   b_covEETrk = fChain->GetBranch("covEETrk");
   b_covTTTrk = fChain->GetBranch("covTTTrk");
   b_covPPTrk = fChain->GetBranch("covPPTrk");
   b_covRRTrk = fChain->GetBranch("covRRTrk");
   b_phicMatTrk = fChain->GetBranch("phicMatTrk");
   b_trkcMatTrk = fChain->GetBranch("trkcMatTrk");
   b_nPidTrk = fChain->GetBranch("nPidTrk");
   b_emcStatusTrk = fChain->GetBranch("emcStatusTrk");
   b_phiAtEMCTrk = fChain->GetBranch("phiAtEMCTrk");
   b_thetaAtEMCTrk = fChain->GetBranch("thetaAtEMCTrk");
   b_isvtTrk = fChain->GetBranch("isvtTrk");
   b_nsvtTrk = fChain->GetBranch("nsvtTrk");
   b_fhitTrk = fChain->GetBranch("fhitTrk");
   b_ndchTrk = fChain->GetBranch("ndchTrk");
   b_lhitTrk = fChain->GetBranch("lhitTrk");
   b_tLenTrk = fChain->GetBranch("tLenTrk");
   b_ntdofTrk = fChain->GetBranch("ntdofTrk");
   b_tproTrk = fChain->GetBranch("tproTrk");
   b_tChi2Trk = fChain->GetBranch("tChi2Trk");
   b_cPidTrk = fChain->GetBranch("cPidTrk");
   b_sfRangeTrk = fChain->GetBranch("sfRangeTrk");
   b_trkFitStatusTrk = fChain->GetBranch("trkFitStatusTrk");
   b_chargeTrk = fChain->GetBranch("chargeTrk");
   b_momentumTrk = fChain->GetBranch("momentumTrk");
   b_ppcov00 = fChain->GetBranch("ppcov00");
   b_ppcov10 = fChain->GetBranch("ppcov10");
   b_ppcov11 = fChain->GetBranch("ppcov11");
   b_ppcov20 = fChain->GetBranch("ppcov20");
   b_ppcov21 = fChain->GetBranch("ppcov21");
   b_ppcov22 = fChain->GetBranch("ppcov22");
   b_xPocaTrk = fChain->GetBranch("xPocaTrk");
   b_yPocaTrk = fChain->GetBranch("yPocaTrk");
   b_zPocaTrk = fChain->GetBranch("zPocaTrk");
   b_thetaTrk = fChain->GetBranch("thetaTrk");
   b_phiTrk = fChain->GetBranch("phiTrk");
   b_muonIdTrk = fChain->GetBranch("muonIdTrk");
   b_elecIdTrk = fChain->GetBranch("elecIdTrk");
   b_kaonIdTrk = fChain->GetBranch("kaonIdTrk");
   b_pionIdTrk = fChain->GetBranch("pionIdTrk");
   b_idTrk = fChain->GetBranch("idTrk");
   b_IndexTrk = fChain->GetBranch("IndexTrk");
   b_IndexNtTrk = fChain->GetBranch("IndexNtTrk");
   b_B0RecTrk = fChain->GetBranch("B0RecTrk");
   b_chBRecTrk = fChain->GetBranch("chBRecTrk");
   b_nGam = fChain->GetBranch("nGam");
   b_IfrLayGam = fChain->GetBranch("IfrLayGam");
   b_IfrNsGam = fChain->GetBranch("IfrNsGam");
   b_IfrInnerGam = fChain->GetBranch("IfrInnerGam");
   b_IfrBarrelGam = fChain->GetBranch("IfrBarrelGam");
   b_IfrFWDGam = fChain->GetBranch("IfrFWDGam");
   b_IfrBWDGam = fChain->GetBranch("IfrBWDGam");
   b_IfrMeasIntLenGam = fChain->GetBranch("IfrMeasIntLenGam");
   b_IfrFirstHitGam = fChain->GetBranch("IfrFirstHitGam");
   b_IfrLastHitGam = fChain->GetBranch("IfrLastHitGam");
   b_IfrExpIntLenGam = fChain->GetBranch("IfrExpIntLenGam");
   b_IfrIntLenBeforeIronGam = fChain->GetBranch("IfrIntLenBeforeIronGam");
   b_IfrTrkMatchGam = fChain->GetBranch("IfrTrkMatchGam");
   b_IfrEmcMatchGam = fChain->GetBranch("IfrEmcMatchGam");
   b_IfrLastBarrelGam = fChain->GetBranch("IfrLastBarrelGam");
   b_IfrCLFitChi2Gam = fChain->GetBranch("IfrCLFitChi2Gam");
   b_IfrStrips0 = fChain->GetBranch("IfrStrips0");
   b_IfrStrips1 = fChain->GetBranch("IfrStrips1");
   b_IfrStrips2 = fChain->GetBranch("IfrStrips2");
   b_IfrStrips3 = fChain->GetBranch("IfrStrips3");
   b_IfrStrips4 = fChain->GetBranch("IfrStrips4");
   b_IfrStrips5 = fChain->GetBranch("IfrStrips5");
   b_IfrStrips6 = fChain->GetBranch("IfrStrips6");
   b_IfrStrips7 = fChain->GetBranch("IfrStrips7");
   b_IfrStrips8 = fChain->GetBranch("IfrStrips8");
   b_IfrStrips9 = fChain->GetBranch("IfrStrips9");
   b_IfrStrips10 = fChain->GetBranch("IfrStrips10");
   b_IfrStrips11 = fChain->GetBranch("IfrStrips11");
   b_IfrStrips12 = fChain->GetBranch("IfrStrips12");
   b_IfrStrips13 = fChain->GetBranch("IfrStrips13");
   b_IfrStrips14 = fChain->GetBranch("IfrStrips14");
   b_IfrStrips15 = fChain->GetBranch("IfrStrips15");
   b_IfrStrips16 = fChain->GetBranch("IfrStrips16");
   b_IfrStrips17 = fChain->GetBranch("IfrStrips17");
   b_IfrStrips18 = fChain->GetBranch("IfrStrips18");
   b_IfrStrips19 = fChain->GetBranch("IfrStrips19");
   b_lMomGam = fChain->GetBranch("lMomGam");
   b_ZMom42Gam = fChain->GetBranch("ZMom42Gam");
   b_ecalGam = fChain->GetBranch("ecalGam");
   b_ecalXGam = fChain->GetBranch("ecalXGam");
   b_ecalYGam = fChain->GetBranch("ecalYGam");
   b_ecalZGam = fChain->GetBranch("ecalZGam");
   b_nCryGam = fChain->GetBranch("nCryGam");
   b_nBumpGam = fChain->GetBranch("nBumpGam");
   b_ZMom20Gam = fChain->GetBranch("ZMom20Gam");
   b_secMomGam = fChain->GetBranch("secMomGam");
   b_s1s9Gam = fChain->GetBranch("s1s9Gam");
   b_s9s25Gam = fChain->GetBranch("s9s25Gam");
   b_erawGam = fChain->GetBranch("erawGam");
   b_phiClusterGam = fChain->GetBranch("phiClusterGam");
   b_thetaClusterGam = fChain->GetBranch("thetaClusterGam");
   b_covEEGam = fChain->GetBranch("covEEGam");
   b_covTTGam = fChain->GetBranch("covTTGam");
   b_covPPGam = fChain->GetBranch("covPPGam");
   b_covRRGam = fChain->GetBranch("covRRGam");
   b_emcStatusGam = fChain->GetBranch("emcStatusGam");
   b_thetaGam = fChain->GetBranch("thetaGam");
   b_phiGam = fChain->GetBranch("phiGam");
   b_energyGam = fChain->GetBranch("energyGam");
   b_idGam = fChain->GetBranch("idGam");
   b_IndexGam = fChain->GetBranch("IndexGam");
   b_IndexNtGam = fChain->GetBranch("IndexNtGam");
   b_B0RecGam = fChain->GetBranch("B0RecGam");
   b_chBRecGam = fChain->GetBranch("chBRecGam");
   return kTRUE;
}

// ----------------------------------------------------------------------

void recoilBase::Show(Int_t entry) {
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

// ----------------------------------------------------------------------
Int_t recoilBase::Cut(Int_t entry) {
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void recoilBase::KLStudy()
{
  static Bool_t firstTime(kTRUE);
  
  if (firstTime) {
    char name[100], title[100];
    TH1 *h;
    fHistFile->mkdir("klstudyks", "klstudyks");
    fHistFile->cd("klstudyks");
    // book the histograms

    sprintf(name, "goodks");  sprintf(title, "goodks");  h = new TH1D(name, title, 4, -.5, 3.5);
    sprintf(name, "ksubmx");  sprintf(title, "ksubmx");  h = new TH1D(name, title, 100, 0., 5.);
    sprintf(name, "ksubdemx");  sprintf(title, "ksubdemx");  h = new TH1D(name, title, 100, -2., 2.);
    sprintf(name, "ksp");  sprintf(title, "ksp");  h = new TH1D(name, title, 100, 0., 5.);
    sprintf(name, "ksth");  sprintf(title, "ksth");  h = new TH1D(name, title, 100, 0., 3.1416);
    sprintf(name, "ksph");  sprintf(title, "ksph");  h = new TH1D(name, title, 100, 0., 6.2832);
    sprintf(name, "kspdm");  sprintf(title, "kspdm");  h = new TH2D(name, title, 100, 0., 5., 100, -2., 2.);
    sprintf(name, "ksthdm");  sprintf(title, "ksthdm");  h = new TH2D(name, title, 100, 0., 3.1416, 100, -2., 2.);
    sprintf(name, "ksphdm");  sprintf(title, "ksphdm");  h = new TH2D(name, title, 100, 0., 6.2832, 100, -2., 2.);
    fHistFile->cd();
  }
  firstTime = kFALSE;
  

  // Take the good  neutral kaons
  

  Double_t submassx,bestdm(-10.0);
  TLorentzVector p4X,p4tmp(0,0,0,0),subp4X(0,0,0,0),p4ks(0,0,0,0);
  int nGoodKS(0);

  mk4Vector( p4X, fPxhad, fTxhad, fFxhad, fMxhad);
  for (Int_t i=0; i<nKs; ++i) {
    if ( goodKshort[i] == 1 ) {
      nGoodKS++;
      mk4Vector( p4tmp, pKs[i], thKs[i], phiKs[i], KAZMASS);
      
      if (bestdm>TMath::Abs(massKs[i]-KAZMASS)) {
        bestdm = TMath::Abs(massKs[i]-KAZMASS);
        p4ks=p4tmp;
        subp4X = p4X - p4tmp;
        submassx = subp4X.M(); 
      }
    }
  }
  fHistFile->cd("klstudyks");
  TH1D * h = (TH1D*)gDirectory->Get("goodks");
  h->Fill(nGoodKS);

  
  if (nGoodKS) {
    TH1D * h2 = (TH1D*)gDirectory->Get("ksubmx");
    TH1D * h3 = (TH1D*)gDirectory->Get("ksubdemx");
    TH1D * h4 = (TH1D*)gDirectory->Get("ksp");
    TH1D * h5 = (TH1D*)gDirectory->Get("ksth");
    TH1D * h6 = (TH1D*)gDirectory->Get("ksph");
    TH2D * h7 = (TH2D*)gDirectory->Get("kspdm");
    TH2D * h8 = (TH2D*)gDirectory->Get("ksthdm");
    TH2D * h9 = (TH2D*)gDirectory->Get("ksphdm");
    h2->Fill(submassx);     
    h3->Fill(submassx-fMxhad); 
    h4->Fill(p4ks.P()); 
    h5->Fill(p4ks.Theta()); 
    h6->Fill(p4ks.Phi()); 
    h7->Fill(bestdm,p4ks.P()); 
    h8->Fill(bestdm,p4ks.Theta()); 
    h9->Fill(bestdm,p4ks.Phi()); 
  }
  
  fHistFile->cd();
  
  if (!fIsMC) return;

 // since this is run AFTER mctruth 
 // fB1Index and fB2Index point to the 2 B
 // fBVxb points to semilept B  
 // fBVxbTyp (D+,Dstar, D pi, Vub)
 // anyway fVub and fVcb contains the number of b->c and b-> semil. decays
 // fOther otherwise
 
 if (fBVxb<0) return;
 Int_t klcand,indexKL(-1),nKL(0);
 Bool_t isKL = kFALSE;
 Double_t pKL(-1000.);
 
 for (klcand=0; klcand<nMc; klcand++) {
   if (idMc[klcand] == 130) {
     if (isAncestor(fBVxb, klcand)) {
       isKL = kTRUE;
       nKL++;
       if (pKL < pMc[klcand]) 
         pKL = pMc[indexKL = klcand];
     }
   }
 }

// Retrive and fill histograms 
 fHistFile->cd("klstudy");
 
 TH1D * h1 = (TH1D*)gDirectory->Get("klklp");
 TH1D * h2 = (TH1D*)gDirectory->Get("klmx");
 TH1D * h3 = (TH1D*)gDirectory->Get("klnmxlo");
 TH2D * h4 = (TH2D*)gDirectory->Get("klmm2mx");
 TH1D * h5 = (TH1D*)gDirectory->Get("klmm2");
 TH1D * h6 = (TH1D*)gDirectory->Get("klmm2lo");
 TH1D * h7 = (TH1D*)gDirectory->Get("klklplo");

 TH1D * h12 = (TH1D*)gDirectory->Get("noklmx");
 TH1D * h13 = (TH1D*)gDirectory->Get("noklnmxlo");
 TH2D * h14 = (TH2D*)gDirectory->Get("noklmm2mx");
 TH1D * h15 = (TH1D*)gDirectory->Get("noklmm2");
 TH1D * h16 = (TH1D*)gDirectory->Get("noklmm2lo");

 TH1D * h22 = (TH1D*)gDirectory->Get("tagklmx");
 TH1D * h26 = (TH1D*)gDirectory->Get("tagklmm2lo");
 TH1D * h32 = (TH1D*)gDirectory->Get("tagnoklmx");
 TH1D * h36 = (TH1D*)gDirectory->Get("tagnoklmm2lo");

 Int_t index(0);
 Bool_t KLtag = klselection(index);
 

 if (isKL && fGoodEvent ) {
   h1->Fill(pKL);
   h2->Fill(fMxhad);
//   cout << 1 << endl;
   if (KLtag) h22->Fill(fMxhad);
//   cout << 2 << endl;
   h4->Fill(fMM2,fMxhad);
   h5->Fill(fMM2);
   if (fMxhad<1.5) {
     h3->Fill(nKL);
     h6->Fill(fMM2); 
//     cout << 3 << endl;
     if (KLtag) h26->Fill(fMM2);
//     cout << 4 << endl;
     h7->Fill(pKL);
   }
   
 }
 
 if(!isKL && fGoodEvent) {
   h12->Fill(fMxhad);
//   cout << 5 << endl;
   if (KLtag) h32->Fill(fMxhad);
//   cout << 6 << endl;
   h14->Fill(fMM2,fMxhad);
   h15->Fill(fMM2);
   if (fMxhad<1.5) { 
     h13->Fill(0.0);
     h16->Fill(fMM2);
//     cout << 7<< endl;
     if (KLtag) h36->Fill(fMM2);
//     cout << 8 << endl;
   }

 }
 
 fHistFile->cd();
}

void recoilBase::bookKLHisto()
{
  char name[100], title[100];
  TH1 *h;

  fHistFile->mkdir("klstudy", "klstudy");  fHistFile->cd("klstudy");

  sprintf (name, "klklp");  sprintf(title, "klklp"); 
  h = new TH1D(name, title, 70, 0.0, 3.5);

  sprintf (name, "klmx");  sprintf(title, "klmx"); 
  h = new TH1D(name, title, 100, 0.0, 5.0);

  sprintf (name, "klnmxlo");  sprintf(title, "klnmxlo"); 
  h = new TH1D(name, title, 4, -0.5, 3.5);

  sprintf (name, "klmm2mx");  sprintf(title, "klmm2mx"); 
  h = new TH2D(name, title, 200, -1.0, 1.0, 100, 0.0, 5.0);

  sprintf (name, "klmm2");  sprintf(title, "klmm2"); 
  h = new TH1D(name, title,  50, -1.0, 1.0);

  sprintf (name, "klmm2lo");  sprintf(title, "klmm2lo"); 
  h = new TH1D(name, title,  50, -1.0, 1.0);

  sprintf (name, "klklplo");  sprintf(title, "klklplo"); 
  h = new TH1D(name, title, 75, 0.0, 3.5);

  sprintf (name, "noklmx");  sprintf(title, "noklmx"); 
  h = new TH1D(name, title,  100, 0.0, 5.0);

  sprintf (name, "noklnmxlo");  sprintf(title, "noklnmxlo"); 
  h = new TH1D(name, title,  4, -0.5, 3.5);

  sprintf (name, "noklmm2mx");  sprintf(title, "noklmm2mx"); 
  h = new TH2D(name, title, 200, -1.0, 1.0, 100, 0.0, 5.0);

  sprintf (name, "noklmm2");  sprintf(title, "noklmm2"); 
  h = new TH1D(name, title, 200, -1.0, 1.0);
  
  sprintf (name, "noklmm2lo");  sprintf(title, "noklmm2lo"); 
  h = new TH1D(name, title, 200, -1.0, 1.0);
  
  sprintf (name, "tagklmx");  sprintf(title, "tagklmx"); 
  h = new TH1D(name, title, 100, 0.0, 5.0);

  sprintf (name, "tagklmm2lo");  sprintf(title, "tagklmm2lo"); 
  h = new TH1D(name, title,  50, -1.0, 1.0);

  sprintf (name, "tagnoklmx");  sprintf(title, "tagnoklmx"); 
  h = new TH1D(name, title,  100, 0.0, 5.0);

  sprintf (name, "tagnoklmm2lo");  sprintf(title, "tagnoklmm2lo"); 
  h = new TH1D(name, title, 200, -1.0, 1.0);
  
  fHistFile->cd();
  
  
}
      
Bool_t recoilBase::klselection( Int_t& klindex )
{
  Bool_t selection = kFALSE;
  
  for (Int_t i=0; i<nGam ; i++) {
    if (ecalGam[i]<= 0) {
      if ((IfrFWDGam[i])&&IfrLastHitGam[i]>14) continue;
      Int_t nlay = IfrLayGam[i];
      if ( nlay>2 ) 
        if ( minTrkClusterOpen(thetaGam[i], phiGam[i]) > .500 ) {
          selection = kTRUE;
          klindex = i;
          break ;
        }
    }
  }
  return selection;
}
  
Double_t recoilBase::minTrkClusterOpen(Double_t theta, Double_t phi)  
{
  Double_t min=-10.0;
  
  TVector3 v1(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
  for (Int_t i=0; i<nTrk; i++) {
    if ( momentumTrk[i]>.750) {
      Double_t theta2 = thetaAtEMCTrk[i];
      Double_t phi2 = phiAtEMCTrk[i];
      
      TVector3 v2(sin(theta2)*cos(phi2),sin(theta2)*sin(phi2),cos(theta2));
      Double_t cosangle = v1 * v2;
      if (min<cosangle) min=cosangle;
    }
  }
  return acos(min);
}

TList * recoilBase::createTrueKLList() 
{
  
  Int_t i ;
  TList * list = new TList;
  
  for (i=0; i<nMc; i++) {
    if (idMc[i] == 130) {
      TLorentzVector * p4tmp = new TLorentzVector;
      mk4Vector( *p4tmp, pMc[i], thetaMc[i], phiMc[i], massMc[i]);
      list->Add( p4tmp );
    }
  }
  
  return list;
}

Bool_t recoilBase::angularMatch(TList& list,  Double_t recTheta, Double_t recPhi) 
{
  Bool_t match(kFALSE);
  
  TIterator*  next = list.MakeIterator();
  TLorentzVector* p4;
  while ((p4 = ((TLorentzVector*) (*next)())) != 0  ) {
    if ( TMath::Abs(p4->Phi()-recPhi) < .200 ) 
      if ( TMath::Abs(p4->Theta()-recTheta) < .200 ) {
        match = kTRUE;
        break ;
      }
  }
  delete next;
  return match;
}
  
void recoilBase::P4XhadKS(TLorentzVector& p4Miss, TLorentzVector& pX,
                         Int_t isVerbose)
{
  
  if (fMxhad>0) {
//    TLorentzVector p4X(0,0,0,0);
//    mk4Vector( p4X, fPxhad, fTxhad, fFxhad, fMxhad);

    Int_t iks = bestKsIndex(0);
    if (iks>=0) {
      TLorentzVector ksp4;
      mk4Vector( ksp4, pKs[iks], thKs[iks], phiKs[iks], KAZMASS);
      pX = p4Xhad - ksp4;
      p4Miss = p4Upsilon - p4Breco - pX;
    }
  }
  if ( isVerbose )
    cout << "ks sub p4x = " << pX.P() << endl;
}

void recoilBase::MxhadchK(Double_t& mass, Double_t& massfit)
{
 
  Double_t M(0);
  
  if (fMxhad>0) {
    TLorentzVector p4X(0,0,0,0),p4(0,0,0,0);
    mk4Vector( p4X, fPxhad, fTxhad, fFxhad, fMxhad);
    for (Int_t i=0; i<nTrk; i++) {
      if (isRecKaon(i)) {
        Double_t p  = momentumTrk[i];
        Double_t ek = sqrt( p*p + KAPMASS*KAPMASS);
        Double_t epi= sqrt( p*p + PIPMASS*PIPMASS);
        Double_t EX = sqrt(fMxhad*fMxhad + fPxhad*fPxhad);
        Double_t E  = EX - ek + epi;
        M = sqrt(E*E - fPxhad*fPxhad);
        continue;
      }
    }
  }
  
  mass=M;
  massfit=M;
}

Double_t recoilBase::MxhadKS(Double_t& missMass,TLorentzVector& pMiss, 
                            TLorentzVector& pX, Int_t isVerbose)
{
  static Double_t counter1 = 0;
  static Double_t counter2 = 0;

  P4XhadKS(p4Miss,pX,isVerbose);
  Double_t M = pX.M();
  missMass = p4Miss.M();
  
  counter1++;
  if (M>0) 
    counter2++;
  if (isVerbose)
    cout << "KS p =" << pX.P() << "x = " << pX.X() << "y= " << pX.Y() 
         << "z= " << pX.Z() << "fMxhad = " << fMxhad 
         << "fnkshort: " << fNKshort << endl;
  
  if (isVerbose)
    cout << "counter1 = " << counter1 << "counter2= " << counter2 << endl;
  return M;

}

Bool_t recoilBase::kill_dupli(int Up, int Low)
{
  // char Bts[100];
  char iUpLo[100];
  bool tmpDupli = kFALSE;

  fValMap = " ";

  sprintf(iUpLo,"%s%X%s%X","00",Up,"/",Low);
  TMapIter *IterMap = new TMapIter(map);
  TObject *NextItem;
  IterMap->Reset();
  NextItem = IterMap->Next();
  while (NextItem != 0) {
    if( NextItem != 0)  {
      //      cout<<iUpLo<<" "<<((TObjString*)(map->GetValue(NextItem)))->GetString()<<endl;
      if(  strcmp(iUpLo,
		  ((TObjString*)(map->GetValue(NextItem)))->GetString())
	   == 0 ) {
	fValMap = ((TObjString*)(map2->GetValue(NextItem)))->GetString();
	tmpDupli = kTRUE;
	break;
      }
      NextItem = IterMap->Next();
    }
  }
  return tmpDupli;
}

void recoilBase::read_killTab(const char *lsfile) {

  if(lsfile == "") return; 

  cout<<" Reading this file for killing: "<<lsfile<<endl;

  int count(0);
  char Ifile[100], Ifile1[100], Ibuff[200], lsbuff[200];
  char ts[30];
  char val[100], tmpfl[100];

  fLkt = 0;
  map = new TMap(10,10);

  map2 = new TMap(10,10);
  TObjString *flag;

  ifstream lsev(lsfile);
  while (lsev.getline(lsbuff, 200, '\n')) {
    sscanf(lsbuff,"%s", Ifile1);
    count++;
    TString *val1 = new TString(Ifile1);
    if ( val1->Contains("a_cand") ) {
      if  ( (val1->Contains("dzdz")) ) {
	sprintf(tmpfl,"1");
      } else if  ( (val1->Contains("dzdc")) || (val1->Contains("dcdz")) ){
	sprintf(tmpfl,"2");
      } else if  ( (val1->Contains("dzds")) || (val1->Contains("dsdz")) ){
	sprintf(tmpfl,"3");
      } else if  ( (val1->Contains("dzdsz")) ||  (val1->Contains("dszdz")) ){
	sprintf(tmpfl,"4");
      }

    } else if ( val1->Contains("b_cand") ) { 
	
      if  ( (val1->Contains("dcdz")) || (val1->Contains("dzdc")) ){
	sprintf(tmpfl,"41");
      } else if  ( (val1->Contains("dcdc")) ){
	sprintf(tmpfl,"42");
      } else if  ( (val1->Contains("dcds")) || (val1->Contains("dsdc")) ){
	sprintf(tmpfl,"43");
      } else if  ( (val1->Contains("dcdsz")) || (val1->Contains("dszdc")) ){
	sprintf(tmpfl,"44");
      }

    } else if ( val1->Contains("c_cand") ) {

      if  ( (val1->Contains("dsdz")) || (val1->Contains("dzds")) ){
	sprintf(tmpfl,"81");
      } else if  ( (val1->Contains("dsdc")) ||  (val1->Contains("dcds")) ){
	sprintf(tmpfl,"82");
      } else if  ( (val1->Contains("dsds")) ){
	sprintf(tmpfl,"83");
      } else if  ( (val1->Contains("dsdsz")) || (val1->Contains("dszds")) ){
	sprintf(tmpfl,"84");
      }
      
    } else if ( val1->Contains("d_cand") ) {

      if  ( (val1->Contains("dszdz")) || (val1->Contains("dzdsz")) ){
	sprintf(tmpfl,"121");
      } else if  ( (val1->Contains("dszdc")) || (val1->Contains("dcdsz")) ){
	sprintf(tmpfl,"122");
      } else if  ( (val1->Contains("dszds")) || (val1->Contains("dsdsz")) ){
	sprintf(tmpfl,"123");
      } else if  ( (val1->Contains("dszdsz")) ) {
	sprintf(tmpfl,"124");
      }

    } else if ( val1->Contains("e_cand") ) {

      if  ( (val1->Contains("dzdz")) ) {
	sprintf(tmpfl,"11");
      } else if  ( (val1->Contains("dzdc")) || (val1->Contains("dcdz")) ){
	sprintf(tmpfl,"12");
      } else if  ( (val1->Contains("dzds")) || (val1->Contains("dsdz")) ){
	sprintf(tmpfl,"13");
      } else if  ( (val1->Contains("dzdsz")) ||  (val1->Contains("dszdz")) ){
	sprintf(tmpfl,"14");
      }

    } else if ( val1->Contains("f_cand") ) {

      if  ( (val1->Contains("dcdz")) || (val1->Contains("dzdc")) ){
	sprintf(tmpfl,"51");
      } else if  ( (val1->Contains("dcdc")) ){
	sprintf(tmpfl,"52");
      } else if  ( (val1->Contains("dcds")) || (val1->Contains("dsdc")) ){
	sprintf(tmpfl,"53");
      } else if  ( (val1->Contains("dcdsz")) || (val1->Contains("dszdc")) ){
	sprintf(tmpfl,"54");
      }

    } else if ( val1->Contains("g_cand") ) {

      if  ( (val1->Contains("dsdz")) || (val1->Contains("dzds")) ){
	sprintf(tmpfl,"91");
      } else if  ( (val1->Contains("dsdc")) ||  (val1->Contains("dcds")) ){
	sprintf(tmpfl,"92");
      } else if  ( (val1->Contains("dsds")) ){
	sprintf(tmpfl,"93");
      } else if  ( (val1->Contains("dsdsz")) || (val1->Contains("dszds")) ){
	sprintf(tmpfl,"94");
      }
      
    } else if ( val1->Contains("h_cand") ) {

      if  ( (val1->Contains("dszdz")) || (val1->Contains("dzdsz")) ){
	sprintf(tmpfl,"131");
      } else if  ( (val1->Contains("dszdc")) || (val1->Contains("dcdsz")) ){
	sprintf(tmpfl,"132");
      } else if  ( (val1->Contains("dszds")) || (val1->Contains("dsdsz")) ){
	sprintf(tmpfl,"133");
      } else if  ( (val1->Contains("dszdsz")) ) {
	sprintf(tmpfl,"134");
      }

    } else if ( val1->Contains("i_cand") ) {

      if  ( (val1->Contains("dzdz")) ) {
	sprintf(tmpfl,"21");
      } else if  ( (val1->Contains("dzdc")) || (val1->Contains("dcdz")) ){
	sprintf(tmpfl,"22");
      } else if  ( (val1->Contains("dzds")) || (val1->Contains("dsdz")) ){
	sprintf(tmpfl,"23");
      } else if  ( (val1->Contains("dzdsz")) ||  (val1->Contains("dszdz")) ){
	sprintf(tmpfl,"24");
      }

    } else if ( val1->Contains("l_cand") ) {

      if  ( (val1->Contains("dcdz")) || (val1->Contains("dzdc")) ){
	sprintf(tmpfl,"61");
      } else if  ( (val1->Contains("dcdc")) ){
	sprintf(tmpfl,"62");
      } else if  ( (val1->Contains("dcds")) || (val1->Contains("dsdc")) ){
	sprintf(tmpfl,"63");
      } else if  ( (val1->Contains("dcdsz")) || (val1->Contains("dszdc")) ){
	sprintf(tmpfl,"64");
      }

    } else if ( val1->Contains("m_cand") ) {

      if  ( (val1->Contains("dsdz")) || (val1->Contains("dzds")) ){
	sprintf(tmpfl,"101");
      } else if  ( (val1->Contains("dsdc")) ||  (val1->Contains("dcds")) ){
	sprintf(tmpfl,"102");
      } else if  ( (val1->Contains("dsds")) ){
	sprintf(tmpfl,"103");
      } else if  ( (val1->Contains("dsdsz")) || (val1->Contains("dszds")) ){
	sprintf(tmpfl,"104");
      }

    } else if ( val1->Contains("n_cand") ) {

      if  ( (val1->Contains("dszdz")) || (val1->Contains("dzdsz")) ){
	sprintf(tmpfl,"141");
      } else if  ( (val1->Contains("dszdc")) || (val1->Contains("dcdsz")) ){
	sprintf(tmpfl,"142");
      } else if  ( (val1->Contains("dszds")) || (val1->Contains("dsdsz")) ){
	sprintf(tmpfl,"143");
      } else if  ( (val1->Contains("dszdsz")) ) {
	sprintf(tmpfl,"144");
      }

    } else if ( val1->Contains("o_cand") ) {

      if  ( (val1->Contains("dzdz")) ) {
	sprintf(tmpfl,"31");
      } else if  ( (val1->Contains("dzdc")) || (val1->Contains("dcdz")) ){
	sprintf(tmpfl,"32");
      } else if  ( (val1->Contains("dzds")) || (val1->Contains("dsdz")) ){
	sprintf(tmpfl,"33");
      } else if  ( (val1->Contains("dzdsz")) ||  (val1->Contains("dszdz")) ){
	sprintf(tmpfl,"34");
      }

    } else if ( val1->Contains("p_cand") ) {

      if  ( (val1->Contains("dcdz")) || (val1->Contains("dzdc")) ){
	sprintf(tmpfl,"71");
      } else if  ( (val1->Contains("dcdc")) ){
	sprintf(tmpfl,"72");
      } else if  ( (val1->Contains("dcds")) || (val1->Contains("dsdc")) ){
	sprintf(tmpfl,"73");
      } else if  ( (val1->Contains("dcdsz")) || (val1->Contains("dszdc")) ){
	sprintf(tmpfl,"74");
      }

    } else if ( val1->Contains("q_cand") ) {

      if  ( (val1->Contains("dsdz")) || (val1->Contains("dzds")) ){
	sprintf(tmpfl,"111");
      } else if  ( (val1->Contains("dsdc")) ||  (val1->Contains("dcds")) ){
	sprintf(tmpfl,"112");
      } else if  ( (val1->Contains("dsds")) ){
	sprintf(tmpfl,"113");
      } else if  ( (val1->Contains("dsdsz")) || (val1->Contains("dszds")) ){
	sprintf(tmpfl,"114");
      }

    } else if ( val1->Contains("r_cand") ) {

      if  ( (val1->Contains("dszdz")) || (val1->Contains("dzdsz")) ){
	sprintf(tmpfl,"151");
      } else if  ( (val1->Contains("dszdc")) || (val1->Contains("dcdsz")) ){
	sprintf(tmpfl,"152");
      } else if  ( (val1->Contains("dszds")) || (val1->Contains("dsdsz")) ){
	sprintf(tmpfl,"153");
      } else if  ( (val1->Contains("dszdsz")) ) {
	sprintf(tmpfl,"154");
      }

    }  else {
      cout<<"Something really fishy"<<endl;
    }

    flag = new TObjString(tmpfl);

    //Tipical: m103dsdsz_0_d_cand
    sprintf(Ifile,"%s%s","kill_stu/test12345/",Ifile1);
    ifstream Iev(Ifile);
    TObject *num;
    num = (TObject *)10;
    TObjString *stval;
    while (Iev.getline(Ibuff, 200, '\n')) {
      sprintf(val,"%d",fLkt);
      stval = new TObjString(val);
      sscanf(Ibuff,"%s", ts);
      TObjString *st = new TObjString(ts);
      map->Add(stval,st);
      map2->Add(stval,flag);
      fLkt++;
    }
  }
  cout << count<< " num of files"<<endl;
}

void recoilBase::KinFit(TLorentzVector& p4XReco, TLorentzVector& p4Nu,
                       TLorentzVector& p4XFit, TLorentzVector& p4NuFit)
{
  int BILTYP; 
  if (fElectron) {
    BILTYP = 2; 
  } else if (fMuon) {
    BILTYP = 1; 
  } 
  
  float BCVAL[4] = {pxUps, pyUps, pzUps, eUps};
  float BP_REC[28] = {
    p4BrecoNC.Px(), p4BrecoNC.Py(), p4BrecoNC.Pz(), p4BrecoNC.E(), 
    p4LeptonLab.Px(), p4LeptonLab.Py(), p4LeptonLab.Pz(), p4LeptonLab.E(), 
    p4XReco.Px(), p4XReco.Py(), p4XReco.Pz(), p4XReco.E(), 
    p4Nu.Px(), p4Nu.Py(), p4Nu.Pz(), p4Nu.E()};
  float BP_FIT[28];
  float BCHI2T, BPROBCHI2; 
  int   BISMEAR(0);
  int   BIERR;
  int BISV = kFALSE;

#if 0
  abcfit_interface_vub_(&BISMEAR,&BILTYP,BCVAL,BP_REC,BP_FIT,&BCHI2T,&BPROBCHI2,&BIERR, &BISV);
#endif

  p4XFit.SetXYZT(BP_FIT[8], BP_FIT[9], BP_FIT[10], BP_FIT[11]);
  p4NuFit.SetXYZT(BP_FIT[12], BP_FIT[13], BP_FIT[14], BP_FIT[15]); 
}


// -------------- general histogram filling functions ---------------------
void recoilBase::fillHist(const char *hname, const char *htitle, 
	      TH1D *templateHist, double value, double weight) {
      TH1D *h=(TH1D*)gDirectory->Get(hname);
      if(h==NULL) {
	h=new TH1D(*templateHist);
	h->SetName(hname);
	h->SetTitle(htitle);
	h->Sumw2();
      }
      h->Fill(value, weight);
}
void recoilBase::fillHist(const char *hname, const char *htitle, 
			  int nbins, double min, double max, double value, double weight) {
      TH1D *h=(TH1D*)gDirectory->Get(hname);
      if(h==NULL) {
	h=new TH1D(hname, htitle, nbins, min, max);
	h->Sumw2();
      }
      h->Fill(value, weight);
}



// -- PID with reco selectors 
 Bool_t recoilBase::isRecEl(Int_t i)     {return (TMath::Abs(elecIdTrk[i]) & IDEL);}
 Bool_t recoilBase::isRecMu(Int_t i)     {return (TMath::Abs(muonIdTrk[i]) & IDMU);}
 Bool_t recoilBase::isRecLepton(Int_t i) {return (isRecEl(i) || isRecMu(i));}
 Bool_t recoilBase::isRecKaon(Int_t i)   {return (TMath::Abs(kaonIdTrk[i]) & IDKA);} 

// -- PID with MC truth  matched information
 Bool_t recoilBase::isMcEl(Int_t i)     {return (TMath::Abs(idTrk[i]) == 11);}
 Bool_t recoilBase::isMcMu(Int_t i)     {return (TMath::Abs(idTrk[i]) == 13);}
 Bool_t recoilBase::isMcLepton(Int_t i) {return (isMcEl(i) || isMcMu(i));}
 Bool_t recoilBase::isMcPi(Int_t i)     {return (TMath::Abs(idTrk[i]) == 211);}
 Bool_t recoilBase::isMcK(Int_t i)     {return (TMath::Abs(idTrk[i]) == 321);}

// -- PID for generator block
 Bool_t recoilBase::isTruPion(Int_t i)  {return (TMath::Abs(idMc[i]) == 211);}
 Bool_t recoilBase::isTruEl(Int_t i)  {return (TMath::Abs(idMc[i]) == 11);}
 Bool_t recoilBase::isTruNuEl(Int_t i)  {return (TMath::Abs(idMc[i]) == 12);}
 Bool_t recoilBase::isTruMu(Int_t i)  {return (TMath::Abs(idMc[i]) == 13);}
 Bool_t recoilBase::isTruNuMu(Int_t i)  {return (TMath::Abs(idMc[i]) == 14);}
 Bool_t recoilBase::isTruTau(Int_t i) {return (TMath::Abs(idMc[i]) == 15);}
 Bool_t recoilBase::isTruLepton(Int_t i) {return (isTruEl(i) || isTruMu(i));}

