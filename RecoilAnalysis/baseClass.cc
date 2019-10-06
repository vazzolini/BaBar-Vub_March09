#include "baseClass.hh"

#include <fstream>
#include <iomanip>

#include "TH2.h"
#include "TH1.h"
#include "TProfile.h"
#include <TMap.h>
#include <TExMap.h>
#include <TObjString.h>
#include "TLorentzVector.h"
#include "TRandom.h"

using namespace std;

#include "splitOff.icc"
#include "mcTruth.icc"
#include "bookHist.icc"
#include "init.icc"

Double_t p_energy_loss_corrected(Double_t pin, Double_t dip, Int_t itype);

// ----------------------------------------------------------------------
baseClass::baseClass(TTree *tree, int isMC, int newFormat) {
  fNewFormat = newFormat; 
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("csx-vubmix-new2002.root");
    if (!f) {
      f = new TFile("csx-vubmix-new2002.root");
    }
    tree = (TTree*)gDirectory->Get("h9");
    
  }
  
  Init(tree,isMC); 
  fToBeCopied = new TEventList("toBeCopied", "Events to be copied", 1000);
  //ADDED CB
  initRest();
}

// ----------------------------------------------------------------------
baseClass::~baseClass() {
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}



// ----------------------------------------------------------------------
void baseClass::readCuts(TString filename, int dump) {
  char  buffer[200];
  fCutFile = filename;
  cout << "Reading " << filename.Data() << " for cut settings" << endl;
  sprintf(buffer, "%s", filename.Data());
  ifstream is(buffer);
  char CutName[100], tablefile[1000];
  float CutValue;
  int ok(0);
  for(int i=0; i<20; i++){
    theweighttrk[i] = 1;
    theweightneu[i] = 1;
  }

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
    if (!strcmp(CutName, "kamomLo")) {KAMOMLO = CutValue; ok = 1;}
    if (!strcmp(CutName, "idEl")) {IDEL = int(CutValue); ok = 1;}
    if (!strcmp(CutName, "idMu")) {IDMU = int(CutValue); ok = 1;}
    if (!strcmp(CutName, "idKa")) {IDKA = int(CutValue); ok = 1;}

    // -- recoil
    if (!strcmp(CutName, "prmm2")) {PRMM2 = CutValue; ok = 1;}
    if (!strcmp(CutName, "mm2Lo")) {MM2LO = CutValue; ok = 1;}
    if (!strcmp(CutName, "mm2Hi")) {MM2HI = CutValue; ok = 1;}
    if (!strcmp(CutName, "reqChargeCoor")) {REQCHARGECORR = (CutValue > 0.5) ?  kTRUE: kFALSE; ok = 1;}
    if (!strcmp(CutName, "reqTotalCharge")) {REQTOTALCHARGE = CutValue; ok = 1;}
    if (!strcmp(CutName, "kspipLo")) {KSPIPLO = CutValue; ok = 1;}
    if (!strcmp(CutName, "kspipHi")) {KSPIPHI = CutValue; ok = 1;}
    if (!strcmp(CutName, "kspipRlo")) {KSPIPRLO = CutValue; ok = 1;}
    if (!strcmp(CutName, "kspizLo")) {KSPIZLO = CutValue; ok = 1;}
    if (!strcmp(CutName, "kspizHi")) {KSPIZHI = CutValue; ok = 1;}
    if (!strcmp(CutName, "ptLo"))    {PTLO = CutValue; ok = 1;}
    if (!strcmp(CutName, "gammaELo")) {GAMMAELO = CutValue; ok = 1;}
    if (!strcmp(CutName, "gammaEHi")) {GAMMAEHI = CutValue; ok = 1;}
    if (!strcmp(CutName, "smearTrkPx"))  {SMEARTRKPX = CutValue; ok = 1;}
    if (!strcmp(CutName, "smearTrkPy"))  {SMEARTRKPY = CutValue; ok = 1;}
    if (!strcmp(CutName, "sigmaNeut")) {SIGMANEUT  = CutValue; ok = 1;}
    if (!strcmp(CutName, "shiftNeut")) {SHIFTNEUT  = CutValue; ok = 1;}
    if (!strcmp(CutName, "trackSelection")) {TRACKSELECTION = CutValue; ok = 1;}
    if (!strcmp(CutName, "trackKilling"))   {DOTRACKKILLING = (int)CutValue; ok = 1;}
    if (!strcmp(CutName, "slowPionKill"))   {SLOWPIONKILL = CutValue; ok = 1;}
    if (!strcmp(CutName, "trackKill")   )   {TRACKKILL = CutValue; ok = 1;}
    if (!strcmp(CutName, "kFitPara")) {fParametrization = CutValue; ok = 1;}
    if (!strcmp(CutName, "photonSelection")) {PHOTONSELECTION = CutValue; ok = 1;}

    if (!strcmp(CutName, "DoBdecayWeight")) {DOBDECWEIGHT = CutValue; ok = 1;}
    if (!strcmp(CutName, "DoDdecayWeight")) {DODDECWEIGHT = CutValue; ok = 1;}

    if (!strcmp(CutName, "photonSelection")) {PHOTONSELECTION = CutValue; ok = 1;}

    // -- Pidmaps and trk efficiency tables
    if (!strcmp(CutName, "smearElPidTables")) {SMEARELPIDTABLES = CutValue; ok = 1;}
    if (!strcmp(CutName, "smearMuPidTables")) {SMEARMUPIDTABLES = CutValue; ok = 1;}
    if (!strcmp(CutName, "smearKaPidTables")) {SMEARKAPIDTABLES = CutValue; ok = 1;}
    if (!strcmp(CutName, "smearElMisTables")) {SMEARELMISTABLES = CutValue; ok = 1;}
    if (!strcmp(CutName, "smearMuMisTables")) {SMEARMUMISTABLES = CutValue; ok = 1;}
    if (!strcmp(CutName, "smearKaMisTables")) {SMEARKAMISTABLES = CutValue; ok = 1;}
    if (!strcmp(CutName, "PidTables")) {
      sscanf(buffer, "%s %s", CutName, tablefile);
      sprintf(PIDTABLES, "%s", tablefile); ok = 1;
      getPidTables();
    } 
    if (!strcmp(CutName, "TrkTables")) {
      sscanf(buffer, "%s %s", CutName, tablefile);
      sprintf(TRKTABLES, "%s", tablefile); ok = 1;
      getTrkTables();
    }


    if (ok == 0)  cout << "==> baseClass::readCuts() Error: Don't know about variable " << CutName << endl;
  }

  if (dump == 1) dumpCuts();

  readintpur();
}



// ----------------------------------------------------------------------
void baseClass::dumpCuts() {
  cout << "====================================" << endl;
  cout << "Cut file " << fCutFile << endl; 
  cout << "------------------------------------" << endl;
  cout << "mesSignal:         " << MESSIGNALLO << " ... " << MESSIGNALHI << endl;
  cout << "mesSideband:       " << MESSIDEBANDLO << " ... " << MESSIDEBANDHI << endl;
  cout << "mesSignalband:     " << MESSIGNALBANDLO << " ... " << MESSIGNALBANDHI << endl;
  cout << "deSignal:          " << DESIGNALLO << " ... " << DESIGNALHI << endl;
  cout << "purity:            " << PURITY << endl;
  cout << "intPurity:         " << INTPURITY << endl;
  cout << "ipurDstar:         " << IPURDSTAR << endl;
  cout << "ipurDc:            " << IPURDC << endl;
  cout << "ipurDstar0:        " << IPURDSTAR0 << endl;
  cout << "ipurD0:            " << IPURD0 << endl;

  cout << "pcmsLo:            " << PCMSLO << endl;
  cout << "plabLo:            " << PLABLO << endl;
  cout << "tlabLo:            " << TLABLO << endl;
  cout << "tlabHi:            " << TLABHI << endl;
  cout << "nLepton:           " << NLEPTON << endl;
  cout << "elmomLo:           " << ELMOMLO << endl;
  cout << "mumomLo:           " << MUMOMLO << endl;
  cout << "kamomLo:           " << KAMOMLO << endl;
  cout << "idEl:              " << IDEL << endl;
  cout << "idMu:              " << IDMU << endl;
  cout << "idKa:              " << IDKA << endl;

  cout << "killTracks:        " << (DOTRACKKILLING & 1 ? 1 : 0) << " with probability " << SLOWPIONKILL << endl;
  cout << "killSlowPions:     " << (DOTRACKKILLING & 2 ? 1 : 0) << " with probability " << TRACKKILL <<  endl;
  cout << "smearTrkPx:        " << SMEARTRKPX << endl;
  cout << "smearTrkPy:        " << SMEARTRKPY << endl;
  cout << "kfitPara:          " << fParametrization << endl;
  cout << "trackSelection:    " << TRACKSELECTION << endl;
  cout << "doTrackKilling:    " << DOTRACKKILLING << endl;
  cout << "photonSelection:   " << PHOTONSELECTION << endl;
  cout << "smearElPidTables:  " << SMEARELPIDTABLES << endl;
  cout << "smearMuPidTables:  " << SMEARMUPIDTABLES << endl;
  cout << "smearKaPidTables:  " << SMEARKAPIDTABLES << endl;
  cout << "smearElMisTables:  " << SMEARELMISTABLES << endl;
  cout << "smearMuMisTables:  " << SMEARMUMISTABLES << endl;
  cout << "smearKaMisTables:  " << SMEARKAMISTABLES << endl;


  cout << "prmm2:             " << PRMM2 << endl;
  cout << "mm2Lo:             " << MM2LO << endl;
  cout << "mm2Hi:             " << MM2HI << endl;
  cout << "reqChargeCoor:     " << (REQCHARGECORR ?  "true": "false") << endl;
  cout << "reqTotalCharge:    " << REQTOTALCHARGE << endl;
  cout << "kspip:             " << KSPIPLO << " ... " << KSPIPHI << endl;
  cout << "kspiz:             " << KSPIZLO << " ... " << KSPIZHI << endl;
  cout << "kspipRlo:          " << KSPIPRLO << endl;
  cout << "pTLo:              " << PTLO << endl;
  cout << "GammaE:            " << GAMMAELO << " ... " << GAMMAEHI << endl;

  cout << "DoBdecayWeight:    " << DOBDECWEIGHT << endl;
  cout << "DoDdecayWeight:    " << DODDECWEIGHT << endl;

  cout << "====================================" << endl;
}


// ----------------------------------------------------------------------
void baseClass::getTrkTables() {
  cout << "Reading TrkTables from " << TRKTABLES << endl;
  ifstream is(TRKTABLES);
  char tableName[1000], buffer[100], fname[100];
  int ok(0);
  while (is.getline(buffer, 200, '\n')) {
    if (buffer[0] == '#') {continue;}
    sscanf(buffer, "%s", tableName);
    cout << tableName << endl;
    sprintf(fname, "trk%d", ok);
    //    fTT[ok].flush();
    fTT[ok] = new TRKTable(tableName, fname);
    ++ok;
  }
}    



// ----------------------------------------------------------------------
void baseClass::getPidTables() {
  cout << "Reading PidTables from " << PIDTABLES << "  " << SMEARELPIDTABLES << "  " << SMEARELMISTABLES << endl;
  ifstream is(PIDTABLES);
  char tableName[1000], selector[100], buffer[1000], fname[200];
  int el(0), mu(0), ka(0), source(0), sink(0);
  while (is.getline(buffer, 200, '\n')) {
    if (buffer[0] == '#') {continue;}
    sscanf(buffer, "%s %d %s %d", tableName, &source, selector, &sink);
    if (TMath::Abs(sink) == 11) {
      cout << "Electron Table " << el << "  ->  " << tableName << endl;
      sprintf(fname, "%s%d_%d", selector, sink, source); 
      //      fPTel[el].flush();
      fPTel[el] = new PIDTable(tableName, fname);
      if (TMath::Abs(source) == TMath::Abs(sink)) {
	cout << "shifting el pidtables relative by " << SMEARELPIDTABLES << endl;
	if (TMath::Abs(SMEARELPIDTABLES) > 0.001) fPTel[el]->shiftRel(SMEARELPIDTABLES);
      } else {
	cout << "shifting el misdtables relative by " << SMEARELMISTABLES << endl;
	if (TMath::Abs(SMEARELMISTABLES) > 0.001) fPTel[el]->shiftRel(SMEARELMISTABLES);
      }
      el++;
    }
    if (TMath::Abs(sink) == 13) {
      cout << "Muon Table " << mu << "  ->  " << tableName << endl;
      sprintf(fname, "%s%d_%d", selector, sink, source); 
      //      fPTmu[mu].flush();
      fPTmu[mu] = new PIDTable(tableName, fname);
      if (TMath::Abs(source) == TMath::Abs(sink)) {
	cout << "shifting mu pidtables relative by " << SMEARMUPIDTABLES << endl;
	if (TMath::Abs(SMEARMUPIDTABLES) > 0.001) fPTmu[mu]->shiftRel(SMEARMUPIDTABLES);
      } else {
	cout << "shifting mu mistables relative by " << SMEARMUMISTABLES << endl;
	if (TMath::Abs(SMEARMUMISTABLES) > 0.001) fPTmu[mu]->shiftRel(SMEARMUMISTABLES);
      }
      mu++;
    }
    if (TMath::Abs(sink) == 321) {
      cout << "Kaon Table " << ka << "  ->  " << tableName << endl;
      sprintf(fname, "%s%d_%d", selector, sink, source); 
      //      fPTka[ka].flush();
      fPTka[ka] = new PIDTable(tableName, fname);
      if (TMath::Abs(source) == TMath::Abs(sink)) {
	cout << "shifting ka pidtables relative by " << SMEARKAPIDTABLES << endl;
	if (TMath::Abs(SMEARKAPIDTABLES) > 0.001) fPTka[ka]->shiftRel(SMEARKAPIDTABLES);
      } else {
	cout << "shifting ka mistables relative by " << SMEARKAMISTABLES << endl;
	if (TMath::Abs(SMEARKAMISTABLES) > 0.001) fPTka[ka]->shiftRel(SMEARKAMISTABLES);
      }
      ka++;
    }
  }
}

// ----------------------------------------------------------------------
void baseClass::readintpur(){    

   char buffer[200];
   float bmode, dmode, sig, bkg, pur, sb;
   ifstream is("tables/tablepurity.dat");
   int mode;
   while (is.getline(buffer, 200, '\n')) {
     if (buffer[0] == '#') {continue;}
     sscanf(buffer, "%f %f %f %f %f %f", &bmode, &dmode, &sig, &bkg, &pur, &sb);
     mode = (dmode+100) * 100 + bmode-10000;
     brecosig[mode] = sig;	
     brecobkg[mode] = bkg;
     brecointpur[mode] = pur; 
     
   }
	
}

// ----------------------------------------------------------------------
Bool_t baseClass::isAncestor(int ancestor, int cand) {
  int mom(cand); 
  while (mom > 0) {
    mom = mothMc[mom]-1; 
    if (mom == ancestor) return true;
  }
  return false;
}

// ----------------------------------------------------------------------
Int_t baseClass::isRecoed(int imc) {
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
void  baseClass::printLorentz(const TLorentzVector &p) {
  char line[200]; 
  sprintf(line, "p = (%7.5f, %7.5f, %7.5f, %7.5f), m = %7.5f",  p.X(), p.Y(), p.Z(), p.E(), p.Mag()); 
  cout << line;
}


// --jump the events (reco Bch) with wrong B-D flavor correlation)
int baseClass::skipBadBreco() {
  int result(0);
  if(bestB0 == 0 && !(modeChB[indexbestB]-14000>299 && modeChB[indexbestB]-14000<399) && !(modeChB[indexbestB]-15000>299 && modeChB[indexbestB]-15000<399) && !(modeChB[indexbestB]-11000>299 && modeChB[indexbestB]-11000<399)){
    int fBrecoChargetmp= 0;
    if(d2ChBLund[indexbestB]!=0&&d2ChBLund[indexbestB]!=111&&d2ChBLund[indexbestB]!=310){ fBrecoChargetmp=fBrecoChargetmp+(TMath::Abs(d2ChBLund[indexbestB])/d2ChBLund[indexbestB]);}
    if(d3ChBLund[indexbestB]!=0&&d3ChBLund[indexbestB]!=111&&d3ChBLund[indexbestB]!=310){ fBrecoChargetmp=fBrecoChargetmp+(TMath::Abs(d3ChBLund[indexbestB])/d3ChBLund[indexbestB]);}
    if(d4ChBLund[indexbestB]!=0&&d4ChBLund[indexbestB]!=111&&d4ChBLund[indexbestB]!=310){ fBrecoChargetmp=fBrecoChargetmp+(TMath::Abs(d4ChBLund[indexbestB])/d4ChBLund[indexbestB]);}
    if(d5ChBLund[indexbestB]!=0&&d5ChBLund[indexbestB]!=111&&d5ChBLund[indexbestB]!=310){ fBrecoChargetmp=fBrecoChargetmp+(TMath::Abs(d5ChBLund[indexbestB])/d5ChBLund[indexbestB]);}
    if(d6ChBLund[indexbestB]!=0&&d6ChBLund[indexbestB]!=111&&d6ChBLund[indexbestB]!=310){ fBrecoChargetmp=fBrecoChargetmp+(TMath::Abs(d6ChBLund[indexbestB])/d6ChBLund[indexbestB]);}
    if(d7ChBLund[indexbestB]!=0&&d7ChBLund[indexbestB]!=111&&d7ChBLund[indexbestB]!=310){ fBrecoChargetmp=fBrecoChargetmp+(TMath::Abs(d7ChBLund[indexbestB])/d7ChBLund[indexbestB]);}  
    if(-1*(TMath::Abs(d1ChBLund[indexbestB])/d1ChBLund[indexbestB]) != fBrecoChargetmp) {
      result = 1;
    }
  }
  return result; 
}


// ----------------------------------------------------------------------
void baseClass::smearTracks() { 
  static Bool_t first(kTRUE);
  if (0 == fOptSmearTracks) return;
  if (fRunnumber < 100000) return; // no smearing for data
  //  if (fRunnumber > 500000) return; // no smearing for SP4

  TH1D *h;
  TDirectory *old = gDirectory; 
  fHistFile->cd();

  if (first) { 
    first = kFALSE; 
    cout << "smearing Tracks with s = " << SMEARTRKPX << endl;
    h = new TH1D("r0", "prec - pgen", 100, -0.1, 0.1);
    h = new TH1D("r2", "prec - pgen 2", 100, -0.1, 0.1);

    h = new TH1D("t0", "trec - tgen", 100, -0.1, 0.1);
    h = new TH1D("t2", "trec - tgen 2", 100, -0.1, 0.1);

    h = new TH1D("f0", "frec - fgen", 100, -0.1, 0.1);
    h = new TH1D("f2", "frec - fgen 2", 100, -0.1, 0.1);

    h = new TH1D("s0", "ptrec - ptgen", 100, -0.1, 0.1);
    h = new TH1D("s2", "ptrec - ptgen 2", 100, -0.1, 0.1);

    h = new TH1D("a0", "1/ptrec - 1/ptgen", 100, -0.1, 0.1);
    h = new TH1D("a2", "1/ptrec - 1/ptgen 2", 100, -0.1, 0.1);
  }

  // smearing is probably happening in pT, 061002 moving to p smearing
  // non-diagonal covariances are ignored at the moment 
  // the additional smearing should increase the end resolution by factors SMEARTRKPX and SMEARTRKPY: 
  //     sigma(tooSmall)**2 + sigma(additional)**2 = s**2 * sigma(tooSmall)**2 
  // ->  sigma(additional)**2 =  s**2 * sigma(tooSmall)**2 - sigma(tooSmall)**2
  double p(0.), p2(0.), px(0.), py(0.), pz(0.), pt(0.), npt(0.), npx(0.), npy(0.), npz(0.), np(0.);
  double t(0.), f(0.), nnt(0.), nnf(0.);
  double shift(0.);
  if (fVerbose) cout << "== Start track list in smearTracks() ==" << endl;
  for (int itrk = 0 ; itrk < nTrk; ++itrk) {
    px  = momentumTrk[itrk]*TMath::Sin(thetaTrk[itrk])*TMath::Cos(phiTrk[itrk]);
    py  = momentumTrk[itrk]*TMath::Sin(thetaTrk[itrk])*TMath::Sin(phiTrk[itrk]);
    pz  = momentumTrk[itrk]*TMath::Cos(thetaTrk[itrk]);
    p   = momentumTrk[itrk];
    t   = thetaTrk[itrk];
    f   = phiTrk[itrk];
    p2  = p*p;
    pt  = TMath::Sqrt(px*px + py*py);
    npx = px;
    npy = py;
    np  = p;
    if (fVerbose) cout << itrk << " unsmeared: " << momentumTrk[itrk]
		       << " (" << thetaTrk[itrk] << ", " << phiTrk[itrk] << ") ";

    double xsadd  = SMEARTRKPX*SMEARTRKPX*ppcov00[itrk] - ppcov00[itrk];
    if (xsadd > 0.) xsadd = TMath::Sqrt(xsadd);
    double xsmear = gRandom->Gaus(shift, xsadd); 
    if (fVerbose) cout <<ppcov00[itrk] << " .. " <<SMEARTRKPX <<" 77 " <<xsmear<<endl;
    if (TMath::Abs(px) > 1.e-3) { npx =  px + xsmear; }

    double ysadd  = SMEARTRKPX*SMEARTRKPX*ppcov11[itrk] - ppcov11[itrk];
    if (ysadd > 0.) ysadd = TMath::Sqrt(ysadd);
    double ysmear = gRandom->Gaus(shift, ysadd); 
    if (TMath::Abs(py) > 1.e-3) { npy = py + ysmear; }

    double zsadd  = SMEARTRKPX*SMEARTRKPX*ppcov22[itrk] - ppcov22[itrk];
    if (zsadd > 0.) zsadd = TMath::Sqrt(zsadd);
    double zsmear = gRandom->Gaus(shift, zsadd); 
    if (TMath::Abs(pz) > 1.e-3) { npz = pz + zsmear; }
    
    npt =  TMath::Sqrt(npx*npx + npy*npy);
    np = TMath::Sqrt(npx*npx + npy*npy + pz*pz);
    double nnp  = TMath::Sqrt(npx*npx + npy*npy + npz*npz);
    double nnpt = TMath::Sqrt(npx*npx + npy*npy);

    nnt =  (npx == 0.0 && npy == 0.0 && npz == 0.0 ? 0.0 : TMath::ATan2(TMath::Sqrt(npx*npx+npy*npy),npz));
    nnf =  (npx == 0.0 && npy == 0.0 ? 0.0 : TMath::ATan2(npy,npx));
    thetaTrk[itrk] = nnt;
    phiTrk[itrk]   = nnf;
    momentumTrk[itrk] = nnp;

    if (fVerbose) cout << " -> smeared: " << momentumTrk[itrk] 
		       << " (" << thetaTrk[itrk] << ", " << phiTrk[itrk] << ") " << endl;

    int imc = IndexTrk[itrk] - 1;
    double mpt(0.);
    if (imc > 0 && imc < nMc) {
      mpt = pMc[imc]*TMath::Sin(thetaMc[imc]);
      if (mpt > 0.2) {
        ((TH1D*)fHistFile->Get("r0"))->Fill(p - pMc[imc]);
        ((TH1D*)fHistFile->Get("r2"))->Fill(nnp - pMc[imc]);

        ((TH1D*)fHistFile->Get("t0"))->Fill(t - thetaMc[imc]);
        ((TH1D*)fHistFile->Get("t2"))->Fill(nnt - thetaMc[imc]);

        ((TH1D*)fHistFile->Get("f0"))->Fill(f - phiMc[imc]);
        ((TH1D*)fHistFile->Get("f2"))->Fill(nnf - phiMc[imc]);

        ((TH1D*)fHistFile->Get("s0"))->Fill(pt - mpt);
        ((TH1D*)fHistFile->Get("s2"))->Fill(nnpt - mpt);
        ((TH1D*)fHistFile->Get("a0"))->Fill(1./pt - 1./mpt);
        ((TH1D*)fHistFile->Get("a2"))->Fill(1./nnpt - 1./mpt);
      }
    }    

  }
  if (fVerbose) cout << "== End track list in smearTracks() ==" << endl;
  old->cd();
}

// ----------------------------------------------------------------------
void baseClass::smearNeut() { 
  static Bool_t first(kTRUE);
  if (0 == fOptSmearNeut) return;
  if (fRunnumber < 100000) return; // no smearing for data
  if (fFileChanged) { 
    if (fRunRange.Contains("Run 1")) {
      fOptSmearNeut = 1; 
      cout << "Smearing/ Shifting Neutrals with RUN1 recipe "  << endl;      
    } else if (fRunRange.Contains("Run 2a")) {
      fOptSmearNeut = 2; 
      cout << "Smearing/ Shifting Neutrals with RUN2 recipe "  << endl;      
    } else if (fRunRange.Contains("Run 2b")) {
      fOptSmearNeut = 2; 
      cout << "Smearing/ Shifting Neutrals with RUN2 recipe "  << endl;      
    } else {
      cout << "Warning: No run range determined. Taking RUN2 as default!" << endl;
      fOptSmearNeut = 2; 
    }
  }

  for (int i = 0; i < nGam; ++i) {
    double tempene = clusterReCorrection(energyGam[i], thetaGam[i]);
    SIGMANEUT = 0;
    if(energyGam[i]< 0.1 ) {
	SIGMANEUT =  0.03;
    } else if(energyGam[i]< 0.3 ) {
	SIGMANEUT =  0.026;
    } else if(energyGam[i]< 0.6 ) {
        SIGMANEUT =  (fOptSmearNeut==1) ? 0.024: 0.016;
    } else if(energyGam[i]< 1. ) {     
        SIGMANEUT =  (fOptSmearNeut==1) ? 0.020: 0.016;
    }
    SIGMANEUT *= energyGam[i];

    tempene = gRandom->Gaus(tempene, SIGMANEUT); 
    TH1D *h;
    TDirectory *old = gDirectory;

    fHistFile->cd();
    if (first) {      
      first = kFALSE;
      h = new TH1D("energypre", " ", 100, 0., 4.);
      h = new TH1D("energypost", " ", 100,  0., 4.);
      h = new TH1D("energysmear", " ", 100,  -0.01, .01);
    }
    ((TH1D*)gDirectory->Get("energypre"))->Fill(energyGam[i]);
    ((TH1D*)gDirectory->Get("energysmear"))->Fill(energyGam[i]-tempene);

    energyGam[i] = tempene;
    ((TH1D*)gDirectory->Get("energypost"))->Fill(energyGam[i]);


  }  
  first = kFALSE;

}

// ----------------------------------------------------------------------
// level = 0  CT
//         1  CTACC
//         2  GTVL
//         3  GTVLACC
//         4  GTL
//         5  GTLACC
//         6  GTVLACC && DCH hits for pt > 0.2
//         7  GTVLACC && DCH hits for pt > 0.2 && looper removal

void baseClass::selectTracks() {
  static Bool_t first(kTRUE);
  static int killedTracks(0);
  static int killedSlowPions(0); 
  if (first == kTRUE) {
    first = kFALSE;
    cout << "-> Selecting tracks  at level " << TRACKSELECTION << endl;
    cout << "-> " << PTLO << " < pT " << endl; 
    if ((DOTRACKKILLING > 0) && (fRunnumber > 100000)) {
      if (DOTRACKKILLING &1) cout << "-> Killing tracks with flat probability " << TRACKKILL << endl;
      if (DOTRACKKILLING &2) cout << "-> Killing in addition slow pions with flat probability " << SLOWPIONKILL << endl;
    } else {
      cout << "-> No track killing " << endl;
    }
  }
  Bool_t ct(kTRUE), gtvl(kFALSE), gtl(kFALSE), acc(kFALSE);
  Double_t pt(0.), dca(0.), dcaz(0.), rand(0.);
  if (fVerbose) cout << "== Start track list in selectTracks() ==" << endl;
  for (int i = 0; i < nTrk; ++i) {
    goodTrack[i] = goodHadron[i] = goodChargedKaon[i] = goodPion[i] = 0; 
    gtvl = gtl = acc = 0;
    pt = momentumTrk[i]*sin(thetaTrk[i]);

    dcaz = zPocaTrk[i] - beamSZ;
    dca  = TMath::Sqrt((xPocaTrk[i]-beamSX)*(xPocaTrk[i]-beamSX) + (yPocaTrk[i]-beamSY)*(yPocaTrk[i]-beamSY));

    ct = kTRUE;
    acc = ((thetaTrk[i] > 0.410) && (thetaTrk[i] < 2.54));
    gtvl = ((pt > 0.0) 
	    && (momentumTrk[i] < 10.0)
	    && (tproTrk[i] >= 0.) 
	    //  	    && ((xPocaTrk[i]*xPocaTrk[i] + yPocaTrk[i]*yPocaTrk[i]) <= 1.5*1.5)
	    //  	    && (zPocaTrk[i] >= -10.0 && zPocaTrk[i] <= 10.0));
	    && (dca  <= 1.5)
	    && (TMath::Abs(dcaz) <= 10.0));
    gtl  = ((pt > 0.1) 
	    && (momentumTrk[i] <= 10.0)
	    && (ndchTrk[i] >= 12)
	    && (tproTrk[i] >= 0.) 
	    //  	    && ((xPocaTrk[i]*xPocaTrk[i] + yPocaTrk[i]*yPocaTrk[i]) <= 1.5*1.5)
	    //  	    && (zPocaTrk[i] >= -10.0 && zPocaTrk[i] <= 10.0));
	    && (dca  <= 1.5)
	    && (TMath::Abs(dcaz) <= 10.0));

    if (TRACKSELECTION == 0) { 
      if (ct) goodTrack[i] = 1; 
    } 
    if (TRACKSELECTION == 1) {
      if (ct && acc) goodTrack[i] = 1;
    }
    if (TRACKSELECTION == 2) {
      if (gtvl) goodTrack[i] = 1;
    }
    if (TRACKSELECTION == 3) {
      if (gtvl && acc) goodTrack[i] = 1;
    }
    if (TRACKSELECTION == 4) {
      if (gtl) goodTrack[i] = 1;
    }
    if (TRACKSELECTION == 5) {
      if (gtl && acc) goodTrack[i] = 1;
    }

    // if it is supposed to have dch hits it should have them (loose because of SVT only track resolution)...
    if (TRACKSELECTION == 6) {
      if (acc && gtvl&& (ndchTrk[i] > 0 || pt< 0.2 )) goodTrack[i] = 1;
    }
    
    if (TRACKSELECTION == 7) {
      if (acc && gtvl&& (ndchTrk[i] > 0 || pt< 0.2 )) goodTrack[i] = 1;
    }


    if (pt < PTLO)  goodTrack[i] = 0;

    // http://www.slac.stanford.edu/BFROOT/www/Physics/TrackEfficTaskForce/TrackingTaskForce-2001.html
    // states: 
    // For GoodTracksVeryLoose and ChargedTracks, if you are working with a very low multiplicity sample 
    // (less than 5 tracks per event) you can set the systematic to 0.5% per track.
    // Otherwise, the systematic must be raised to 1.3% per track since the multiplicity dependence 
    // of these efficiencies is not well known. 
    //
    // For the special case of slow pions (momenta below about 200 MeV) 
    // a systematic uncertainty of 1.6% should be used.


    // -- Track killing in MC: 
    //    DOTRACKKILLING&1: All tracks with flat probability 1.3%
    //    DOTRACKKILLING&2: In addition, slow pions can be killed with 1.6%
    if (fRunnumber > 100000) {
      if (DOTRACKKILLING & 1) {
	rand = gRandom->Rndm(); 
	if (rand < TRACKKILL) {
	  ++killedTracks; 
	  goodTrack[i] = 0; 
	  if (0 == killedTracks%1000) cout << "killed " << killedTracks << " tracks so far (in the entire event)" << endl;
	}      
      }
      if (DOTRACKKILLING & 2) {
        if ((pt < 0.200) && (rand < SLOWPIONKILL) && (goodTrack[i] == 1)) {
          goodTrack[i] = 0; 
	  ++killedSlowPions; 
	  if (0 == killedSlowPions%100) cout << "killed " << killedSlowPions << " slow pions so far (in the entire event)" << endl;
        }
      }
    }

    if (fVerbose) cout << i << "  " << momentumTrk[i] << "  " << goodTrack[i] << endl;
  }

  if (TRACKSELECTION == 7) {
    cleanGoodTracks(); 
  }


  // -- Initialize all particle flag arrays with track selection result
  for (int j = 0; j < nTrk; ++j) {
    goodHadron[j] = goodTrack[j]; 
    goodChargedKaon[j] = goodTrack[j]; 
    goodPion[j] = goodTrack[j]; 
  }

  if (fVerbose) cout << "== End track list in selectTracks() ==" << endl;
}



// ----------------------------------------------------------------------
// level = 0 E_gamma > 0.080 GeV
//         1 GPL
//         2 GPLACC
//         3 GPD
//         4 GPDACC
//         5 Cuts on: low/high E, low LAT
//         6 Cuts on: low/high E, low/high LAT, low s9s25
//         7 Cuts on: low/high E, low/high LAT, low s9s25, veto unmatched clusters
void baseClass::selectPhotons() {
  fHistFile->cd();
  static Bool_t first(kTRUE);
  if (first == kTRUE) {
    first = kFALSE;
    cout << "-> Selecting photons at level " << PHOTONSELECTION << endl;
  }
  Bool_t gpl(kTRUE), gpd(kFALSE), acc(kFALSE), superric(kFALSE), klsel(kFALSE), ric(kFALSE), ricTight(kFALSE), gg(kFALSE);
  for (int i = 0; i < nGam; ++i) {
    goodPhoton[i] = 0;
    gpl = gpd = acc = ric = kFALSE;
    acc = ((thetaGam[i] > 0.410) && (thetaGam[i] < 2.54));
    gpl = ((energyGam[i] >= 0.030)
           && (nCryGam[i] >= 1.)
           && (lMomGam[i] <= 0.8));
    gpd = ((energyGam[i] >= 0.100)
	   && (nCryGam[i] >= 1.)
	   && (lMomGam[i] <= 0.8));

    TLorentzVector p4Gam(0., 0., 0., 0.);
    mk4Vector(p4Gam,energyGam[i],thetaGam[i],phiGam[i],0);
    if(p4Brecoil.T()<=0)cout <<" photon selection screwed because recoil B energy is "<<p4Brecoil.T()<<endl;
    p4Gam.Boost(-p4Brecoil.BoostVector());

    superric = energyGam[i] >= 0.08 && p4Gam.T()<2.8 && lMomGam[i]>0.05 && lMomGam[i]<0.5 && s9s25Gam[i]>0.9;
    klsel = energyGam[i] >= 0.08 && p4Gam.T()<2.8 && KLlikeEMC(i)<0;
        
    ric = energyGam[i] >= 0.08 && energyGam[i]<4. &&  lMomGam[i]>0.05;
    ricTight = energyGam[i] >= 0.08&& energyGam[i]<4. && lMomGam[i]>0.05 && lMomGam[i]<0.5 && s9s25Gam[i]>0.9;

    gg = ((energyGam[i] >= 0.030)
	  && (nCryGam[i] >= 1.)
	  && (lMomGam[i] <= 0.8)
	  && (p4Gam.T() < 2.8)
	  && (splitOffGam[i] == 0)
	  );

    if (PHOTONSELECTION == 0) { 
      if (energyGam[i] > 0.080) goodPhoton[i] = 1; 
    } 
    if (PHOTONSELECTION == 1) {
      if (gg && acc) goodPhoton[i] = 1;
      continue;
    }
    if (PHOTONSELECTION == 2) {
      if (gpl && acc) goodPhoton[i] = 1;
    }
    if (PHOTONSELECTION == 3) {
      if (gpd) goodPhoton[i] = 1;
    }
    if (PHOTONSELECTION == 4) {
      if (gpd && acc) goodPhoton[i] = 1;
    }
    if (PHOTONSELECTION == 5) {
      if (ric && acc)
	 goodPhoton[i] = 1;
    }
    if (PHOTONSELECTION == 6) {
      if (ricTight && acc)
         goodPhoton[i] = 1;
    }
    if (PHOTONSELECTION == 7) {
      if (ricTight && acc &&  splitOffGam[i]==0 )
         goodPhoton[i] = 1;
    }
    if (PHOTONSELECTION == 8) {
      if (superric && acc &&  splitOffGam[i]==0 )
         goodPhoton[i] = 1;
    }
    if (PHOTONSELECTION == 9) {
      if (klsel && acc &&  splitOffGam[i]==0 )
         goodPhoton[i] = 1;
    }

    if (fVerbose) cout << i << "  " << energyGam[i] <<" "<<thetaGam[i]<<" "<<phiGam[i]<<p4Gam.T()<<" "<<lMomGam[i]<<" "<<s9s25Gam[i]<<" "<<splitOffGam[i]<< "  " << goodPhoton[i] << endl;
  }
}



// ----------------------------------------------------------------------
Bool_t baseClass::isPidKillEl(int i) {
  if ((idTrk[i]) == -211)  return fPTel[4]->idR(momentumTrk[i], thetaTrk[i], phiTrk[i]);
  if ((idTrk[i]) == 211)   return fPTel[5]->idR(momentumTrk[i], thetaTrk[i], phiTrk[i]);
  if ((idTrk[i]) == -321)  return fPTel[6]->idR(momentumTrk[i], thetaTrk[i], phiTrk[i]);
  if ((idTrk[i]) == 321)   return fPTel[7]->idR(momentumTrk[i], thetaTrk[i], phiTrk[i]);
  if ((idTrk[i]) == 11)    return fPTel[0]->idR(momentumTrk[i], thetaTrk[i], phiTrk[i]);
  if ((idTrk[i]) == -11)   return fPTel[1]->idR(momentumTrk[i], thetaTrk[i], phiTrk[i]);
  if ((idTrk[i]) == 13)    return fPTel[2]->idR(momentumTrk[i], thetaTrk[i], phiTrk[i]);
  if ((idTrk[i]) == -13)   return fPTel[3]->idR(momentumTrk[i], thetaTrk[i], phiTrk[i]);
  if ((idTrk[i]) == -2212) return fPTel[8]->idR(momentumTrk[i], thetaTrk[i], phiTrk[i]);
  if ((idTrk[i]) == 2212)  return fPTel[9]->idR(momentumTrk[i], thetaTrk[i], phiTrk[i]);
  return kFALSE;
}

// ----------------------------------------------------------------------
Bool_t baseClass::isPidKillMu(int i) {
  if ((idTrk[i]) == -211)  return fPTmu[4]->idR(momentumTrk[i], thetaTrk[i], phiTrk[i]);
  if ((idTrk[i]) == 211)   return fPTmu[5]->idR(momentumTrk[i], thetaTrk[i], phiTrk[i]);
  if ((idTrk[i]) == -321)  return fPTmu[6]->idR(momentumTrk[i], thetaTrk[i], phiTrk[i]);
  if ((idTrk[i]) == 321)   return fPTmu[7]->idR(momentumTrk[i], thetaTrk[i], phiTrk[i]);
  if ((idTrk[i]) == 11)    return kFALSE; // p.d. electrons don't fake muons
  if ((idTrk[i]) == -11)   return kFALSE; // p.d. electrons don't fake muons
  if ((idTrk[i]) == 13)    return fPTmu[2]->idR(momentumTrk[i], thetaTrk[i], phiTrk[i]);
  if ((idTrk[i]) == -13)   return fPTmu[3]->idR(momentumTrk[i], thetaTrk[i], phiTrk[i]);
  if ((idTrk[i]) == -2212) return fPTmu[8]->idR(momentumTrk[i], thetaTrk[i], phiTrk[i]);
  if ((idTrk[i]) == 2212)  return fPTmu[9]->idR(momentumTrk[i], thetaTrk[i], phiTrk[i]);
  return kFALSE;
}

// ----------------------------------------------------------------------
Bool_t baseClass::isPidKillKaon(int i) {
  if ((idTrk[i]) == -211)  return fPTka[4]->idR(momentumTrk[i], thetaTrk[i], phiTrk[i]);
  if ((idTrk[i]) == 211)   return fPTka[5]->idR(momentumTrk[i], thetaTrk[i], phiTrk[i]);
  if ((idTrk[i]) == -321)  return fPTka[6]->idR(momentumTrk[i], thetaTrk[i], phiTrk[i]);
  if ((idTrk[i]) == 321)   return fPTka[7]->idR(momentumTrk[i], thetaTrk[i], phiTrk[i]);
  if ((idTrk[i]) == 11)    return fPTka[0]->idR(momentumTrk[i], thetaTrk[i], phiTrk[i]);
  if ((idTrk[i]) == -11)   return fPTka[1]->idR(momentumTrk[i], thetaTrk[i], phiTrk[i]);
  if ((idTrk[i]) == 13)    return fPTka[2]->idR(momentumTrk[i], thetaTrk[i], phiTrk[i]);
  if ((idTrk[i]) == -13)   return fPTka[3]->idR(momentumTrk[i], thetaTrk[i], phiTrk[i]);
  if ((idTrk[i]) == -2212) return fPTka[8]->idR(momentumTrk[i], thetaTrk[i], phiTrk[i]);
  if ((idTrk[i]) == 2212)  return fPTka[9]->idR(momentumTrk[i], thetaTrk[i], phiTrk[i]);
  return kFALSE;	
}



// ----------------------------------------------------------------------
void baseClass::fillPidMaps() {
  int i(0);
  static bool first(kTRUE), firstK(kTRUE), firstE(kTRUE), firstM(kTRUE);

  // -- Protection against running PidKilling on Data
  if ((fRunnumber < 100000) && (fOptPidKilling || fOptPidKillingKaon || fOptPidKillingEl || fOptPidKillingMu)) {
    if (first) cout << "-> Resetting PidKilling options, since this is data" << endl;
    fOptPidKilling = kFALSE;
    fOptPidKillingKaon = kFALSE;
    fOptPidKillingEl = kFALSE;
    fOptPidKillingMu = kFALSE;
  }

  //    cout << "----------------------------------------------------------------------" << endl;
  Bool_t gtl;
  Double_t pt(0.), dca(0.), dcaz(0.);
  if (fOptPidKilling > 0) {
    if (first) { 
      first = kFALSE; 
      cout << "-> Running with general pidkilling" << endl;
    }
    for (i = 0; i < nTrk; ++i) {
      recEl[i] = 0; 
      recMu[i] = 0; 
      recKa[i] = 0; 
      if (isPidKillKaon(i)) recKa[i] = 1; 
      if (momentumTrk[i] < KAMOMLO) recKa[i] = 0;
      // -- Tighter track requirements for leptons
      pt = momentumTrk[i]*sin(thetaTrk[i]);
      dcaz = zPocaTrk[i] - beamSZ;
      dca  = TMath::Sqrt((xPocaTrk[i]-beamSX)*(xPocaTrk[i]-beamSX) + (yPocaTrk[i]-beamSY)*(yPocaTrk[i]-beamSY));
      gtl = ((pt > 0.1) 
	     && (momentumTrk[i] <= 10.0)
	     && (ndchTrk[i] >= 12)
	     && (tproTrk[i] >= 0.) 
	     && (dca  <= 1.5)
	     && (TMath::Abs(dcaz) <= 10.0));
      if (kFALSE == gtl) continue;
      if (isPidKillEl(i))   recEl[i] = 1; 
      if (isPidKillMu(i))   recMu[i] = 1; 
    } 
  } else {
    if (first) { 
      first = kFALSE; 
      cout << "-> Running without pidkilling" << endl;
    }
    for (i = 0; i < nTrk; ++i) {
      recKa[i] = 0; 
      recEl[i] = 0;
      recMu[i] = 0;
      // -- Kaons
      if (TMath::Abs(kaonIdTrk[i]) & IDKA) {
	recKa[i] = 1; 
	if (momentumTrk[i] < KAMOMLO) recKa[i] = 0; 
      } else {
	recKa[i] = 0; 
      }
      if(fOptPidKillingKaon > 0){
	if (firstK) { 
	  firstK = kFALSE; 
	  cout << "-> Running with Kaon pidkilling" << endl;	
	}
	if (isPidKillKaon(i)) {
	  recKa[i] = 1; 
	  if (momentumTrk[i] < KAMOMLO) recKa[i] = 0;
	} else {
	  recKa[i] = 0;
	}
      } 
      // -- Tighter track requirements for leptons
      pt = momentumTrk[i]*sin(thetaTrk[i]);
      dcaz = zPocaTrk[i] - beamSZ;
      dca  = TMath::Sqrt((xPocaTrk[i]-beamSX)*(xPocaTrk[i]-beamSX) + (yPocaTrk[i]-beamSY)*(yPocaTrk[i]-beamSY));
      gtl = ((pt > 0.1) 
	     && (momentumTrk[i] <= 10.0)
	     && (ndchTrk[i] >= 12)
	     && (tproTrk[i] >= 0.) 
	     && (dca  <= 1.5)
	     && (TMath::Abs(dcaz) <= 10.0));
      
      //      cout <<" track " << i << " ele " << elecIdTrk[i] << " " << IDEL << " mu " <<muonIdTrk[i]<<" " << IDMU <<" kaon " << kaonIdTrk[i]<<" " <<IDKA<<endl;
      if (kFALSE == gtl) continue;

      if (TMath::Abs(elecIdTrk[i]) & IDEL) {recEl[i] = 1; } 
      if ((TMath::Abs(muonIdTrk[i]) & IDMU) && ((TMath::Abs(kaonIdTrk[i]) & IDKA) == 0))  {recMu[i] = 1; } 

      if (fOptPidKillingEl > 0) {
	if (firstE) { 
	  firstE = kFALSE; 
	  cout << "-> Running with Electron pidkilling" << endl;	
	}
	if (isPidKillEl(i)) {recEl[i] = 1;} else {recEl[i] = 0;}
      } 

      if (fOptPidKillingMu > 0) {
	if (firstM) { 
	  firstM = kFALSE; 
	  cout << "-> Running with Muon pidkilling" << endl;	
	}
	if (isPidKillMu(i)) {recMu[i] = 1;} else {recMu[i] = 0;}

      } 

    }
  }
  
  if (fVerbose) {
    cout << "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << endl;
    for (i = 0; i < nTrk; ++i) {
      if (isRecEl(i))   cout << "Track " << i << " is an electron with p = " << momentumTrk[i] << endl;
      if (isRecMu(i))   cout << "Track " << i << " is a  muon with p = " << momentumTrk[i] << endl;
      if (isRecKaon(i)) cout << "Track " << i << " is a kaon with p = " << momentumTrk[i] << endl;
    }
    cout << "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << endl;
  }    
}

void baseClass::breco() {

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
    
    sprintf(name, "h200");  sprintf(title, "mes");  h = new TH1D(name, title, 40, 5.2, 5.3); 
    sprintf(name, "h300");  sprintf(title, "mes dE all");  h2 = new TH2D(name, title, 40, 5.2, 5.3, 40, -0.2, 0.2); 
    
    sprintf(name, "mesallDupli");  sprintf(title, "mes All Dupli");  h = new TH1D(name, title, 40, 5.2, 5.3); 
    sprintf(name, "deallDupli");  sprintf(title, "delta E All Dupli");  h = new TH1D(name, title, 40, -0.1, 0.1); 
    sprintf(name, "mesdeDupli");  sprintf(title, "mes dE all Dupli");  h2 = new TH2D(name, title, 40, 5.2, 5.3, 40, -0.2, 0.2); 

  }

  
  // #trx #pi0 #ks
  nnpi0=0; nntrk=0; nnks=0; nnpar=0;
  if(!fChB) {
    //B0
    if(d1B0Lund[indexbestB]==111) {nnpi0++;}
    else if(d1B0Lund[indexbestB]==310) {nnks++;}
    else if(d1B0Lund[indexbestB]!=0) {nntrk++;}
    if(d2B0Lund[indexbestB]==111) {nnpi0++;}
    else if(d2B0Lund[indexbestB]==310) {nnks++;}
    else if(d2B0Lund[indexbestB]!=0) {nntrk++;}
    if(d3B0Lund[indexbestB]==111) {nnpi0++;}
    else if(d3B0Lund[indexbestB]==310) {nnks++;}
    else if(d3B0Lund[indexbestB]!=0) {nntrk++;}
    if(d4B0Lund[indexbestB]==111) {nnpi0++;}
    else if(d4B0Lund[indexbestB]==310) {nnks++;}
    else if(d4B0Lund[indexbestB]!=0) {nntrk++;}
    if(d5B0Lund[indexbestB]==111) {nnpi0++;}
    else if(d5B0Lund[indexbestB]==310) {nnks++;}
    else if(d5B0Lund[indexbestB]!=0) {nntrk++;}
    if(d6B0Lund[indexbestB]==111) {nnpi0++;}
    else if(d6B0Lund[indexbestB]==310) {nnks++;}
    else if(d6B0Lund[indexbestB]!=0) {nntrk++;}
    if(d7B0Lund[indexbestB]==111) {nnpi0++;}
    else if(d7B0Lund[indexbestB]==310) {nnks++;}
    else if(d7B0Lund[indexbestB]!=0) {nntrk++;}
  } else {  
    //ChB
    if(d1ChBLund[indexbestB]==111) {nnpi0++;}
    else if(d1ChBLund[indexbestB]==310) {nnks++;}
    else if(d1ChBLund[indexbestB]!=0) {nntrk++;}
    if(d2ChBLund[indexbestB]==111) {nnpi0++;}
    else if(d2ChBLund[indexbestB]==310) {nnks++;}
    else if(d2ChBLund[indexbestB]!=0) {nntrk++;}
    if(d3ChBLund[indexbestB]==111) {nnpi0++;}
    else if(d3ChBLund[indexbestB]==310) {nnks++;}
    else if(d3ChBLund[indexbestB]!=0) {nntrk++;}
    if(d4ChBLund[indexbestB]==111) {nnpi0++;}
    else if(d4ChBLund[indexbestB]==310) {nnks++;}
    else if(d4ChBLund[indexbestB]!=0) {nntrk++;}
    if(d5ChBLund[indexbestB]==111) {nnpi0++;}
    else if(d5ChBLund[indexbestB]==310) {nnks++;}
    else if(d5ChBLund[indexbestB]!=0) {nntrk++;}
    if(d6ChBLund[indexbestB]==111) {nnpi0++;}
    else if(d6ChBLund[indexbestB]==310) {nnks++;}
    else if(d6ChBLund[indexbestB]!=0) {nntrk++;}
    if(d7ChBLund[indexbestB]==111) {nnpi0++;}
    else if(d7ChBLund[indexbestB]==310) {nnks++;}
    else if(d7ChBLund[indexbestB]!=0) {nntrk++;}
  }
  nnpar=nnpi0+nnks+nntrk;

  

  fHistFile->cd("breco");
	
  fBrecoMc = MCB0[indexbestB];
  fMes = mseB0[indexbestB];
  fDeltaE = deltaeB0[indexbestB];
  int tmpdauB = d1B0Lund[indexbestB];
  int tmpmodeB = modeB0[indexbestB];
  if(fChB) {
    fBrecoMc = MCChB[indexbestB];
    fMes = mseChB[indexbestB];
    fDeltaE = deltaeChB[indexbestB];
    tmpdauB = d1ChBLund[indexbestB];
    tmpmodeB = modeChB[indexbestB];
  }

  fBmode = tmpmodeB;

  // BetaCoreTools/BtaExclusiveDecayList.hh
  // the following is an arbitrary definition
  // fSeedMode = 0 dc
  // fSeedMode = 1 dstar
  // fSeedMode = 2 d0
  // fSeedMode = 3 dstar0
  if ((11000 <= tmpmodeB) &&  (tmpmodeB < 12000)) {
    fSeedMode = 2;
  } else if  ((12000 <= tmpmodeB) &&  (tmpmodeB < 13000)) {
    fSeedMode = 0;
  } else if  ((13000 <= tmpmodeB) &&  (tmpmodeB < 14000)) {
    fSeedMode = 1;
  } else if  ((14000 <= tmpmodeB) &&  (tmpmodeB < 16000)) {
    fSeedMode = 3;
  } else {
    fSeedMode = -1;
  }

  fBrecoCharge= 0;
  fBrecoFlavor = 1;  

  if (fChB) {
     if(d2ChBLund[indexbestB]!=0&&d2ChBLund[indexbestB]!=111&&d2ChBLund[indexbestB]!=310){ fBrecoCharge=fBrecoCharge+(TMath::Abs(d2ChBLund[indexbestB])/d2ChBLund[indexbestB]);}	
     if(d3ChBLund[indexbestB]!=0&&d3ChBLund[indexbestB]!=111&&d3ChBLund[indexbestB]!=310){ fBrecoCharge=fBrecoCharge+(TMath::Abs(d3ChBLund[indexbestB])/d3ChBLund[indexbestB]);}
     if(d4ChBLund[indexbestB]!=0&&d4ChBLund[indexbestB]!=111&&d4ChBLund[indexbestB]!=310){ fBrecoCharge=fBrecoCharge+(TMath::Abs(d4ChBLund[indexbestB])/d4ChBLund[indexbestB]);}
     if(d5ChBLund[indexbestB]!=0&&d5ChBLund[indexbestB]!=111&&d5ChBLund[indexbestB]!=310){ fBrecoCharge=fBrecoCharge+(TMath::Abs(d5ChBLund[indexbestB])/d5ChBLund[indexbestB]);}
     if(d6ChBLund[indexbestB]!=0&&d6ChBLund[indexbestB]!=111&&d6ChBLund[indexbestB]!=310){ fBrecoCharge=fBrecoCharge+(TMath::Abs(d6ChBLund[indexbestB])/d6ChBLund[indexbestB]);}
     if(d7ChBLund[indexbestB]!=0&&d7ChBLund[indexbestB]!=111&&d7ChBLund[indexbestB]!=310){ fBrecoCharge=fBrecoCharge+(TMath::Abs(d7ChBLund[indexbestB])/d7ChBLund[indexbestB]);}	
     fBrecoFlavor=fBrecoCharge;	
  }else{
     fBrecoCharge=0;
     fBrecoFlavor=-1*(TMath::Abs(d1B0Lund[indexbestB])/d1B0Lund[indexbestB]);
  }

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
  
  fHistFile->cd();
}

void baseClass::findbestB()
{
    double tmpintpur=-999;
    indexbestB=-999;
    bestB0=1;

    static int first(1); 
    TH1D *h; 
    if (first && (fHistFile != 0)) {
      char title[200];
      char name[200];
      first = 0;

      fHistFile->cd();
      sprintf(name, "meszero");  sprintf(title, "mes zero");  h = new TH1D(name, title, 40, 5.2, 5.3);
      sprintf(name, "mesbest");  sprintf(title, "mes best");  h = new TH1D(name, title, 40, 5.2, 5.3);

      fHistFile->mkdir("bestb", "bestb");
      fHistFile->cd("bestb");
      sprintf(name, "h101");  sprintf(title, "int purity B0");  h = new TH1D(name, title, 500, 0., 1.); 
      sprintf(name, "h111");  sprintf(title, "int purity B+");  h = new TH1D(name, title, 500, 0., 1.); 
    }
    
    for (int iB0=0; iB0<nB0; iB0++){
      int mode = modeB0[iB0]-10000;
      if(brecointpur[mode]>tmpintpur) {
	indexbestB = iB0;    
        tmpintpur = brecointpur[mode];
      }
    }
    
    for (int iChB=0; iChB<nChB; iChB++){
      int mode = modeChB[iChB]-10000;
      if(brecointpur[mode]>tmpintpur) {
	indexbestB = iChB;    
        tmpintpur = brecointpur[mode];
	bestB0=0;
      }
    }

    fChB = 0;
    if (bestB0 == 0) {
      fChB = 1;
    }

    brecoOverlap = 1; // this is flag for overlap with BRECO candidate
    if(indexbestB == 1) brecoOverlap = 2;

    if(fHistFile!=0){
      fHistFile->cd();
      ((TH1D*)gDirectory->Get("meszero"))->Fill(mseB0[0]);
      ((TH1D*)gDirectory->Get("meszero"))->Fill(mseChB[0]);
      if(bestB0) { 
        ((TH1D*)gDirectory->Get("mesbest"))->Fill(mseB0[indexbestB]);
      }else{
        ((TH1D*)gDirectory->Get("mesbest"))->Fill(mseChB[indexbestB]);
      }
    }

    if (fHistFile != 0) {
      fHistFile->cd("bestb");
      if(bestB0) {
	((TH1D*)gDirectory->Get("h101"))->Fill(tmpintpur);
      } else {
	((TH1D*)gDirectory->Get("h111"))->Fill(tmpintpur);
      }
      fHistFile->cd();
    }

}

// ----------------------------------------------------------------------
void baseClass::calculateEvtW8() {
  static int first(1);
  
  if (first) {
    first = 0; 

    TH1D *h;
    char name[100], title[100];
    fHistFile->cd();
    sprintf(name,"bw8");sprintf(title,"B decay weights"); h = new TH1D(name, title, 200, 0., 2.0); h->Sumw2();
    sprintf(name,"dw8");sprintf(title,"B decay weights"); h = new TH1D(name, title, 200, 0., 2.0); h->Sumw2();
    sprintf(name,"ew8");sprintf(title,"Evt weights");     h = new TH1D(name, title, 200, 0., 2.0); h->Sumw2();
  }

  int dImode(0); 
  double b = getBsysweight(fBVxbTyp,fVub); 
  ((TH1D*)gDirectory->Get("bw8"))->Fill(b);
  double d = getDsysweight(fDpi,fDk,fDks,fDpiz,fDlep,dImode,fVub);
  ((TH1D*)gDirectory->Get("dw8"))->Fill(d);
  fEvtW8 = b * d * fEvtW8; 
  ((TH1D*)gDirectory->Get("ew8"))->Fill(b);
}

// ----------------------------------------------------------------------
double baseClass::getBsysweight(int decType,int thevub) {
  double theweight;
  theweight = 1.;  
  if (fIsMC) {
    if(DOBDECWEIGHT) theweight = Bsem->weight(decType); 
    if(thevub) theweight = 1.;
  }
  return theweight;
}

// ----------------------------------------------------------------------
double baseClass::getDsysweight(int decDpi,int decDk,int decDks,int decDpiz,int decDlep,int decImode,int thevub) {
  double theweight;
  theweight = 1.;  
  int bla(1);
  if (fIsMC) {
    if(DODDECWEIGHT) theweight = Dvar->weight(decDpi,decDk,decDks,decDpiz,decDlep,decImode,bla);    
    if(thevub) theweight = 1.;
  }
  return theweight;
}



//------------------------------------------------------------------------------//
//
// Description: Function for the Correction of the Cluster Energy
//              => returns the corrected cluster energy
//
//   clusterCorrection(double rawEnergy, double clusterPostionTheta, bool newCorr=true)
//     => returns corrected energy with
//        - constants from March    2002 ("new") if newCorr=true (default)
//        - constants from November 1999 ("old") if newCorr=false
//
//   clusterReCorrection(oldCorrectedEnergy, clusterPostionTheta)
//     => convert old corrected energy
//             to new corrected energy
//
//              rawEnergy           - raw cluster energy (in GeV)
//              clusterPostionTheta - theta of the cluster positon (in rad)
//              oldCorrectedEnergy  - old corrected cluster energy (in GeV)
//
// Author : Enrico Maly  25 Apr 2002
//
//------------------------------------------------------------------------------//

double
baseClass::clusterCorrection(const double rawEnergy,
                  const double clusterPositionTheta, 
                  const bool   newCorr) const
{

// constants from April 2002
 double constants[19] = {
   +9.046e-01,
   +1.243e-02,
   +7.785e-03,
   -2.178e-03,
   +1.620e-02,
   -5.686e-03,
   +2.063e-04,
   +1.408e-01,
   -1.133e-01,
   +2.214e-02,
   +8.409e-03,
   -1.650e-02,
   +5.301e-03,
   +4.998e-02,
   -2.964e-02,
   +9.960e-01,
   -9.097e-02,
   +4.351e-02,
   -5.891e-03
 };
 
 // constants from November 1999
 double oldConstants[19] = {
   +1.024e-00,
   -1.093e-01,
   +4.528e-02,
   -5.959e-03,
   +3.955e-04,
   +4.375e-04,
   -2.855e-04,
   +1.643e-02,
   -1.881e-02,
   +4.838e-03,
   +1.583e-02,
   -1.680e-02,
   +5.074e-03,
   +1.989e-02,
   -1.907e-02,
   +1.133e-00,
   -2.420e-01,
   +9.396e-02,
   -1.142e-02,
 };

  const double lgE = log10(1000*rawEnergy);
  const double cosTh = cos(clusterPositionTheta);

  double *c;
  if (newCorr) c=constants;
  else c=oldConstants;

  double eCorr=rawEnergy;

  // barrel section
  if (cosTh<=0.892)
    {

      const double p0 = c[0]+c[1]*lgE+c[2]*lgE*lgE+c[3]*pow(lgE,3);
      const double p1 = c[4]+c[5]*lgE+c[6]*lgE*lgE;
      const double p2 = c[7]+c[8]*lgE+c[9]*lgE*lgE;
      const double p3 = c[10]+c[11]*lgE+c[12]*lgE*lgE;
      const double p4 = c[13]+c[14]*lgE;

      const double correctionForBarrel = p0+p1*cosTh+p2*cosTh*cosTh+p3*pow(cosTh,3)+p4*pow(cosTh,4);
      eCorr = rawEnergy/correctionForBarrel;

    }
  // endcap section
  else if (cosTh>0.892)
    {
      const double correctionForEndcap = c[15]+c[16]*lgE+c[17]*lgE*lgE+c[18]*pow(lgE,3);
      eCorr = rawEnergy/correctionForEndcap;
    }

  return eCorr;

}
double
baseClass::clusterReCorrection(const double oldCorrectedEnergy, 
                    const double clusterPositionTheta) const
{

  double eCorr=oldCorrectedEnergy;
  const double eRaw=eCorr*eCorr/clusterCorrection(eCorr,clusterPositionTheta,false);
  eCorr=clusterCorrection(eRaw,clusterPositionTheta,true);

  return eCorr;

}
	
double	
baseClass::KLlikeEMC(int iCand) const
{	
  //-------- Background PDFs for the data of RUN2 --------------------------
  double back_ecal[10] = { 0.5275, 0.21681, 0.10349, 0.05486, 0.03238, 0.02183, 0.01579, 0.01172, 0.00884, 0.00677};
  double back_lat[10] = { 0.09198, 0.26943, 0.28318, 0.13907, 0.06998, 0.04745, 0.0368, 0.03117, 0.02151, 0.00942};
  double back_crystals[10] = { 0.03202, 0.42276, 0.34644, 0.11175, 0.03933, 0.02358, 0.01408, 0.00618, 0.00261, 0.00126};
  double back_secmom[10] = { 0.43028, 0.38517, 0.07238, 0.03521, 0.0262, 0.01545, 0.01215, 0.00955, 0.00753, 0.00608};
  double back_zern20[10] = { 0.00745, 0.00894, 0.01011, 0.01139, 0.01581, 0.02165, 0.03045, 0.05589, 0.25287, 0.58544};
  double back_zern42[10] = { 0.66203, 0.19298, 0.05003, 0.03441, 0.02507, 0.01716, 0.00951, 0.00573, 0.00273, 0.34011E-03};
  double back_s1s9[8] = { 0., 0.93556E-05, 0.00472, 0.04841, 0.14427, 0.20776, 0.24495, 0.34988};
  double back_s9s25[6] = { 0., 0., 0.29512E-03, 0.01094, 0.14129, 0.84748};
  //----------------------------------------------------------------------------
  
  //----------- Signal PDFs for the data of RUN2 --------------------------
  double sig_ecal[10] = { 0.23928, 0.33592, 0.22081, 0.10395, 0.05966, 0.01193, 0.01193, 0.76422E-03, 0.76422E-03, 0.76422E-03};
  double sig_lat[10] = { 0.1364, 0.11018, 0.11018, 0.10434, 0.10614, 0.1407, 0.12379, 0.10709, 0.03059, 0.03059};
  double sig_crystals[10] = { 0.10184, 0.24081, 0.2948, 0.19563, 0.10597, 0.04163, 0.00644, 0.00644, 0.00644, 0.00644};
  double sig_secmom[10] = { 0.22105, 0.22105, 0.18421, 0.12933, 0.10712, 0.05343, 0.04217, 0.02627, 0.00769, 0.00769};
  double sig_zern20[10] = { 0.00746, 0.00746, 0.00746, 0.02248, 0.02248, 0.05591, 0.11991, 0.22323, 0.26681, 0.26681};
  double sig_zern42[10] = { 0.33813, 0.3211, 0.15336, 0.12222, 0.01716, 0.01716, 0.00772, 0.00772, 0.00772, 0.00772};
  double sig_s1s9[8] = { 0.01492, 0.11158, 0.19429, 0.23105, 0.18548, 0.13799, 0.06235, 0.06235};
  double sig_s9s25[6] = {0.01065, 0.01065, 0.05403, 0.12455, 0.22319, 0.57693};
  //--------------------------------------------------------------------------------
  
  //**** Inizialization of the ****
  //*****various contributions ****
  //**to the likelihood function***
  double like1 = 0.;
  double like2 = 0.;
  double like3 = 0.;
  double like4 = 0.;
  double like5 = 0.;
  double like6 = 0.;
  double like7 = 0.;
  double like8 = 0.;
  //********************
  
  double likeli = 0.;
  double step = 0.;
  int n_chan = 0;
  
  
  //___________ For each variable, access to the correct values of _________
  //______________ the signal and background PDFs and evaluate _____________
  //_________ the corresponding contribution to the likelihood function ____
  n_chan = 10;
  step = 0.2;
  for(int ja = 1; ja <= n_chan; ja++) {
    int jja = ja-1;
    if(ecalGam[iCand]>(0.2+(ja-1)*step) && ecalGam[iCand]<(0.2+ja*step)) {
      if(back_ecal[jja] != 0. && sig_ecal[jja] != 0.) {
	like1 = log(float(sig_ecal[jja])) - log(float(back_ecal[jja]));
      }
    }
  }
  
  step = 0.1;
  for(int jb = 1; jb <= n_chan; jb++){
    int jjb = jb-1;
    if(lMomGam[iCand]>=((jb-1)*step) && lMomGam[iCand]<(jb*step)) {
      if(back_lat[jjb] != 0. && sig_lat[jjb] != 0.) {
	like2 = log(float(sig_lat[jjb])) - log(float(back_lat[jjb]));
      }
    }
  }
  
  step = 5.;
  for(int jc = 1; jc <= n_chan; jc++){
    int jjc = jc-1;
    if(nCryGam[iCand]>=((jc-1)*step) && nCryGam[iCand]<(jc*step)) {
      if(back_crystals[jjc] != 0. && sig_crystals[jjc] != 0.) {
	like3 = log(float(sig_crystals[jjc])) - log(float(back_crystals[jjc]));
      }
    }
  }
  
  step = 0.001;
  for(int jd = 1; jd <= n_chan; jd++){
    int jjd = jd-1;
    if(secMomGam[iCand]>=((jd-1)*step) && secMomGam[iCand]<(jd*step)) {
      if(back_secmom[jjd] != 0. && sig_secmom[jjd] != 0.) {
	like4 = log(float(sig_secmom[jjd])) - log(float(back_secmom[jjd]));
      }
    }
  }
  
  step = 0.1;
  for(int je = 1; je <= n_chan; je++){
    int jje = je-1;
    if(ZMom20Gam[iCand]>=((je-1)*step) && ZMom20Gam[iCand]<(je*step)) {
      if(back_zern20[jje] != 0. && sig_zern20[jje] != 0.) {
	like5 = log(float(sig_zern20[jje])) - log(float(back_zern20[jje]));
      }
    }
  }
  
  step = 0.05;
  for(int jf = 1; jf <= n_chan; jf++){
    int jjf = jf-1;
    if(ZMom42Gam[iCand]>=((jf-1)*step) && ZMom42Gam[iCand]<(jf*step)) {
      if(back_zern42[jjf] != 0. && sig_zern42[jjf] != 0.) {
	like6 = log(float(sig_zern42[jjf])) - log(float(back_zern42[jjf]));
      }
    }
  }
  
  step = 0.1;
  n_chan = 8;
  for(int jg = 1; jg <= n_chan; jg++){
    int jjg = jg-1;
    if(s1s9Gam[iCand]>=(0.2+(jg-1)*step) && s1s9Gam[iCand]<(0.2+jg*step)) {
      if(back_s1s9[jjg] != 0. && sig_s1s9[jjg] != 0.) {
	like7 = log(float(sig_s1s9[jjg])) - log(float(back_s1s9[jjg]));
      }
    }
  }
  
  step = 0.1;
  n_chan = 6;
  for(int jh = 1; jh <= n_chan; jh++){
    int jjh = jh-1;
    if(s9s25Gam[iCand]>=(0.4+(jh-1)*step) && s9s25Gam[iCand]<(0.4+jh*step)) {
      if(back_s9s25[jjh] != 0. && sig_s9s25[jjh] != 0.) {
	like8 = log(float(sig_s9s25[jjh])) - log(float(back_s9s25[jjh]));
      }
    }
  }
  //____________________________________________________________________________
  
  likeli = like1 + like2 + like3 + like4 + like5 + like6 + like7 + like8;
  
  return likeli;
  
}

double
baseClass::KLlikeIFR(int iCand) const
{
  //-------- Background PDFs for the data of RUN2 --------------------------
    double back_first[13] = {0.60094, 0.10427, 0.07344, 0.03872, 0.03763, 0.02899, 0.02053, 0.01668, 0.01419, 0.01786, 0.01649, 0.01215, 0.01811};
    double back_last[16] = {0.10622, 0.11292, 0.11276, 0.09837, 0.10227, 0.08403, 0.05166, 0.04824, 0.0444, 0.03653, 0.02859, 0.02598, 0.02811, 0.01636, 0.02523, 0.07834};
    double back_layers[10] = {0.39475, 0.31462, 0.13947, 0.07263, 0.03439, 0.0188, 0.01093, 0.0069, 0.00449, 0.00301};
    double back_length[18] = {0.24839, 0.20055, 0.11651, 0.07728, 0.04762, 0.03555, 0.02201,0.02003, 0.01516, 0.01234, 0.00989, 0.00989, 0.00918, 0.00636, 0.00874, 0.00991, 0.06256, 0.06256};
    double back_strips[6] = {0.30616, 0.41523, 0.13778, 0.07169, 0.0426, 0.02655};
    double back_multip[6] = {0.1402, 0.50683, 0.23513, 0.0678, 0.03268, 0.01735};
    double back_sigma[5] = {0.49318, 0.31015, 0.12636, 0.05229, 0.01802};
    //----------------------------------------------------------------------------

  //----------- Signal PDFs for the data of RUN2 --------------------------
    double sig_first[13] = {0.16, 0.16, 0.16, 0.09481, 0.09481, 0.09481, 0.02844, 0.02844, 0.02844, 0.03756, 0.03756, 0.03756, 0.03756};
    double sig_last[16] = {0.05955, 0.05955, 0.05955, 0.05955, 0.09936, 0.09936, 0.09936, 0.04142, 0.04142, 0.04142, 0.07989, 0.07989, 0.07989, 0.03327, 0.03327, 0.03327};
    double sig_layers[10] = {0.18847, 0.21875, 0.15265, 0.15265, 0.15265, 0.02421, 0.02421, 0.02421, 0.02421, 0.02421};
    double sig_length[18] = {0.20804, 0.18532, 0.13347, 0.13347, 0.06371, 0.06371, 0.06371, 0.0195, 0.0195, 0.0195, 0.0195, 0.01011, 0.01011, 0.01011, 0.01011, 0.01011, 0.01011, 0.01011};
    double sig_strips[6] = {0.04651, 0.42515, 0.21216, 0.15645, 0.08475, 0.07498};
    double sig_multip[6] = {0.0903, 0.49604, 0.31731, 0.03212, 0.03212, 0.03212};
    double sig_sigma[5] = {0.20203, 0.33329, 0.2966, 0.08404, 0.08404};
    //--------------------------------------------------------------------------

  //**** Inizialization of the ****
  //*****various contributions ****
  //**to the likelihood function***
    double like1 = 0.;
    double like2 = 0.;
    double like3 = 0.;
    double like4 = 0.;
    double like5 = 0.;
    double like6 = 0.;
    double like7 = 0.;
    //******************************

    double likeli = 0.;
    double step = 0.;
    int n_chan = 0;


  //___________ For each variable, access to the correct values of _________
  //______________ the signal and background PDFs and evaluate _____________
  //_________ the corresponding contribution to the likelihood function ____
    n_chan = 13;
    for(int ja = 1; ja <= n_chan; ja++) {
       if(IfrFirstHitGam[iCand] == ja) {
       int jja = ja-1;
         if(back_first[jja] != 0. && sig_first[jja] != 0.) {
           like1 = log(float(sig_first[jja])) - log(float(back_first[jja]));
         }
       }
    }

    n_chan = 16;
    for(int jb = 1; jb <= n_chan; jb++){
       int jjb = jb-1;
       if(IfrLastHitGam[iCand] == jb) {
         if(back_last[jjb] != 0. && sig_last[jjb] != 0.) {
           like2 = log(float(sig_last[jjb])) - log(float(back_last[jjb]));
         }
       }
    }

    n_chan = 10;
    for(int jc = 1; jc <= n_chan; jc++){
       int jjc = jc-1;
       if(IfrLayGam[iCand] == jc) {
         if(back_layers[jjc] != 0. && sig_layers[jjc] != 0.) {
           like3 = log(float(sig_layers[jjc])) - log(float(back_layers[jjc]));
         }
       }
    }

    n_chan = 18;
    int length = IfrLastHitGam[iCand] - IfrFirstHitGam[iCand] + 1;
    for(int jd = 1; jd <= n_chan; jd++){
       int jjd = jd-1;  
       if(length == jd) {
         if(back_length[jjd] != 0. && sig_length[jjd] != 0.) {
           like4 = log(float(sig_length[jjd])) - log(float(back_length[jjd]));
         }
       }
    }

    n_chan = 6;
    step = 10.;
    for(int je = 1; je <= n_chan; je++){
       int jje = je-1;
       if(IfrNsGam[iCand] >= ((je -1)*step) && IfrNsGam[iCand] < (je*step)) {
         if(back_strips[jje] != 0. && sig_strips[jje] != 0.) {
           like5 = log(float(sig_strips[jje])) - log(float(back_strips[jje]));
         }
       }
    }
        
    n_chan = 6;
    step = 3.;
    double multip = IfrNsGam[iCand]/IfrLayGam[iCand];
    for(int jf = 1; jf <= n_chan; jf++){
       int jjf = jf-1;
       if(multip >= ((jf -1)*step) && multip < (jf*step)) {
         if(back_multip[jjf] != 0. && sig_multip[jjf] != 0.) {
           like6 = log(float(sig_multip[jjf])) - log(float(back_multip[jjf]));
         }
       }
    }

    step = 2.;
    n_chan = 5;
    //------- Sigma -------------------------
    double sqr_diff = 0.;
    int n_lay= 0;
    double mult_sum = 0.;
    // uglyyyyyyyyyyyyyyyyyy
    int nStrips[20];
    nStrips[0]=IfrStrips0[iCand];
    nStrips[1]=IfrStrips1[iCand];
    nStrips[2]=IfrStrips2[iCand];
    nStrips[3]=IfrStrips3[iCand];
    nStrips[4]=IfrStrips4[iCand];
    nStrips[5]=IfrStrips5[iCand];
    nStrips[6]=IfrStrips6[iCand];
    nStrips[7]=IfrStrips7[iCand];
    nStrips[8]=IfrStrips8[iCand];
    nStrips[9]=IfrStrips9[iCand];
    nStrips[10]=IfrStrips10[iCand];
    nStrips[11]=IfrStrips11[iCand];
    nStrips[12]=IfrStrips12[iCand];
    nStrips[13]=IfrStrips13[iCand];
    nStrips[14]=IfrStrips14[iCand];
    nStrips[15]=IfrStrips15[iCand];
    nStrips[16]=IfrStrips16[iCand];
    nStrips[17]=IfrStrips17[iCand];
    nStrips[18]=IfrStrips18[iCand];
    nStrips[19]=IfrStrips19[iCand];

    for(int iii= IfrFirstHitGam[iCand]; iii<=IfrLastHitGam[iCand]; iii++) {
       if(iii != 0) {

         if(nStrips[iii] != 0) {
           mult_sum = mult_sum + nStrips[iii];
           n_lay = n_lay +1;
           sqr_diff = sqr_diff + ((nStrips[iii] - multip)*(nStrips[iii] - multip));
         } // nStrips(iii) != 0
       } // iii != 0
    } // end loop firstHit()<iii<lastHit()
    //cout<<"    sqr_diff: "<<sqr_diff<< endl;
    double sigma_multip = sqrt(sqr_diff/(n_lay - 1)); 
    if(n_lay == 1) sigma_multip = sqrt(sqr_diff);
    //-----------------------------------------------
    for(int jg = 1; jg <= n_chan; jg++){
       int jjg = jg-1;
       if(sigma_multip >= ((jg-1)*step) && sigma_multip < (jg*step)) {
         if(back_sigma[jjg] != 0. && sig_sigma[jjg] != 0.) {
           like7 = log(float(sig_sigma[jjg])) - log(float(back_sigma[jjg]));
         }
       }
    }
    //___________________________________________________________________________

    likeli = like2 + like3 + like4 + like5 + like7;

    return likeli;
}




void baseClass::timestamp(const char * lun) {

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

  if(lun == "") return;

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

  char name[100];
  char name2[100];

  sprintf(name,"%s",lun);
  sprintf(name2,"%s%s",lun,"MI");


  ofstream outFile(name,ios::app);
  ofstream outFile2(name2,ios::app);
  outFile.setf(ios::uppercase);
  outFile2.setf(ios::uppercase);
  int pur(0);// tmp fix
  outFile<<hex<<setw(nplat)<<setfill('0')<<cplat<<":"<<setw(npart)
	 <<setfill('0')<<cpart<<":"<<setw(8)<<setfill('0')<<upperID
	 <<"/"<<setw(8)<<setfill('0')<<lowerID<<":"<<array[k-1]<<"  run="<<dec<<setw(8)
	 <<setfill('0')<<runNumber<<endl;
  outFile2<<hex<<setw(nplat)<<setfill('0')<<cplat<<":"<<setw(npart)
	 <<setfill('0')<<tspart<<":"<<setw(8)<<setfill('0')<<upperID
	 <<"/"<<setw(8)<<setfill('0')<<lowerID<<":"<<array[k-1]<<"  run="<<dec<<setw(8)
	 <<setfill('0')<<runNumber<<"  "<<dec<<pur<<" "
	 <<fBmode<<" "<<endl;
}


// ----------------------------------------------------------------------
void baseClass::cleanGoodTracks(int what) {

  static int doHist(0);

  int osLooper  = 1;
  int ssLooper  = 1; 
  int ssParallel= 1; 
  int print(0); 

  // -- ca. 5 sigma
  const double phiSSmatch = 0.220; 
  const double phiOSmatch = 0.190; 
  const double thetaSSmatch = 0.215;  // from central part fit
  const double thetaOSmatch = 0.300; 
  const double ptmatch = 0.1; 

  const double paraphimatch = 0.1; 
  const double parathetamatch = 0.1; 
  const double paraptmatch = 0.05; 

  static Bool_t first(kTRUE);
  int osmult[100], ssfmult[100], sssmult[100];
  for (int i= 0; i < 100; ++i) {
    osmult[i] = ssfmult[i] = sssmult[i]  = 0; 
  }


  TH1D *h1; 
  TH2D *h2; 
  if (doHist && (first == kTRUE)) {
    first = kFALSE;
    fHistFile->cd();
    fHistFile->mkdir("tracks");
    fHistFile->cd("tracks"); 

    h1 = new TH1D("g0",   "Nevent", 10, 0., 10.);
    h1 = new TH1D("g1",   "Ntrk", 30, 0., 30.);
    h1 = new TH1D("g10",  "opposite sign multiplicity", 6, 0., 6.);
    h1 = new TH1D("g11",  "slow same sign multiplicity", 6, 0., 6.);
    h1 = new TH1D("g12",  "fast same sign multiplicity", 6, 0., 6.);

    h1 = new TH1D("g100", "all tracks, pt < 0.18", 100, 0., 3.15);
    h1 = new TH1D("g101", "goodTracks, pt < 0.18", 100, 0., 3.15);
    h1 = new TH1D("g102", "goodTracks after removal, pt < 0.18", 100, 0., 3.15);
    h2 = new TH2D("g103", "goodTracks after removal, pt < 0.18", 100, 0., 3.15, 100, -3.15, 3.15);
    h2 = new TH2D("g104", "goodTracks after removal, pt < 0.18", 100, 0., 3.15, 100, 0., 0.2);

    h1 = new TH1D("g110", "all good tracks", 100, 0., 3.15);
    h1 = new TH1D("g111", "all good tracks", 100, 0., 1.0);
    h2 = new TH2D("g112", "all good tracks", 100, 0., 3.15, 100, 0., 1.0);

    h1 = new TH1D("g120", "osLooper", 100, 0., 3.15);
    h1 = new TH1D("g121", "osLooper", 100, 0., 1.0);
    h2 = new TH2D("g122", "osLooper", 100, 0., 3.15, 100, 0., 1.0);

    h1 = new TH1D("g130", "ssLooper slow", 100, 0., 3.15);
    h1 = new TH1D("g131", "ssLooper slow", 100, 0., 1.0);
    h2 = new TH2D("g132", "ssLooper slow", 100, 0., 3.15, 100, 0., 3.15);

    h1 = new TH1D("g140", "ss parallel", 100, 0., 3.15);
    h1 = new TH1D("g141", "ss Parallel", 100, 0., 1.0);
    h2 = new TH2D("g142", "ss Parallel", 100, 0., 3.15, 100, 0., 1.0);


    // --  os loopers
    h1 = new TH1D("h100", "os delta phi,    pt < 0.18", 500, -6.50, 6.50);     
    h1 = new TH1D("k100", "os delta phi,    pt < 0.18", 500, -6.50, 6.50);     
    h2 = new TH2D("H100", "os delta phi vs. pt", 50, 0., 0.2, 50, -6.50, 6.50);     

    h1 = new TH1D("h101", "os delta theta,  pt < 0.18", 500, -2.00, 2.00);     
    h1 = new TH1D("k101", "os delta theta,  pt < 0.18", 500, -2.00, 2.00);     
    h2 = new TH2D("H101", "os delta theta vs. pt", 50, 0., 0.2, 50, -2.0, 2.0);     

    h1 = new TH1D("h102", "os delta pt,    pt < 0.18", 500, -.50, .50);     
    h1 = new TH1D("k102", "os delta pt,    pt < 0.18", 500, -.50, .50);     
    h2 = new TH2D("H102", "os delta pt vs. pt", 50, 0., 0.2, 50, -0.5, 0.5);     

    h1 = new TH1D("h103", "os delta tan(dip),  pt < 0.18", 500, -2.00, 2.00);     
    h1 = new TH1D("k103", "os delta tan(dip),  pt < 0.18", 500, -2.00, 2.00);     
    h2 = new TH2D("H103", "os delta tan(dip) vs. pt", 50, 0., 0.2, 50, -2.0, 2.0);     

    // --  ss loopers
    h1 = new TH1D("h200", "ss delta phi,    pt < 0.18", 500, -6.50, 6.50);     
    h1 = new TH1D("k200", "ss delta phi,    pt < 0.18", 500, -6.50, 6.50);     
    h2 = new TH2D("H200", "ss delta phi vs. pt", 50, 0., 0.2, 50, -6.50, 6.50);     

    h1 = new TH1D("h201", "ss delta theta,  pt < 0.18", 500, -2.00, 2.00);     
    h1 = new TH1D("k201", "ss delta theta,  pt < 0.18", 500, -2.00, 2.00);     
    h2 = new TH2D("H201", "ss delta theta vs. pt", 50, 0., 0.2, 50, -2.0, 2.0);     

    h1 = new TH1D("h202", "ss delta pt,    pt < 0.18", 500, -.50, .50);     
    h1 = new TH1D("k202", "ss delta pt,    pt < 0.18", 500, -.50, .50);     
    h2 = new TH2D("H202", "ss delta pt vs. pt", 50, 0., 0.2, 50, -0.5, 0.5);     

    h1 = new TH1D("h203", "os delta tan(dip),  pt < 0.18", 500, -2.00, 2.00);     
    h1 = new TH1D("k203", "os delta tan(dip),  pt < 0.18", 500, -2.00, 2.00);     
    h2 = new TH2D("H203", "os delta tan(dip) vs. pt", 50, 0., 0.2, 50, -2.0, 2.0);     

    // --  ss parallel tracks
    h1 = new TH1D("h300", "ssPara delta phi", 500, -6.50, 6.50);     
    h1 = new TH1D("k300", "ssPara delta phi", 500, -6.50, 6.50);     
    h2 = new TH2D("H300", "ssPara delta phi vs. pt", 50, 0., 2.0, 50, -6.50, 6.50);     

    h1 = new TH1D("h301", "ssPara delta theta", 500, -2.00, 2.00);     
    h1 = new TH1D("k301", "ssPara delta theta", 500, -2.00, 2.00);     
    h2 = new TH2D("H301", "ssPara delta theta vs. pt", 100, 0., 2.0, 50, -2.0, 2.0);     

    h1 = new TH1D("h302", "ssPara delta pt", 500, -1.50, 1.50);     
    h1 = new TH1D("k302", "ssPara delta pt", 500, -1.50, 1.50);     
    h2 = new TH2D("H302", "ssPara delta pt vs. pt", 50, 0., 0.2, 50, -0.5, 0.5);     

    fHistFile->cd(); 
  }

  if (doHist) fHistFile->cd("tracks"); 

  double pt(-99.), ndch(-99.), nsvt(-99.), charge(-99.), dca(-99.), dcaz(-99.), dz(-99.);
  double phi(-99.), theta(-99.); 
  int kill(0); 

  Bool_t histSSPhi, histOSPhi, histSSTheta, histOSTheta, histPt; 

  double pi = TMath::Pi(); 

  for (int i = 0; i < nTrk; ++i) {
    pt     = momentumTrk[i]*sin(thetaTrk[i]);
    nsvt   = nsvtTrk[i]; 
    ndch   = ndchTrk[i]; 
    charge = chargeTrk[i]; 
    phi    = phiTrk[i]; 
    theta  = thetaTrk[i];
    
    dz   = zPocaTrk[i] - primVtxZ;
    dcaz = zPocaTrk[i] - beamSZ;
    dca  = TMath::Sqrt((xPocaTrk[i]-beamSX)*(xPocaTrk[i]-beamSX) + (yPocaTrk[i]-beamSY)*(yPocaTrk[i]-beamSY));

    if (goodTrack[i] == 0) continue; 

    for (int j = 0; j < nTrk; ++j) {
      if (j == i) continue; 
      if (goodTrack[j] == 0) continue; 
      kill = j; 

      histSSPhi   = (TMath::Abs(phi - phiTrk[j]) < phiSSmatch); 
      histOSPhi   = (TMath::Abs(TMath::Abs(phi - phiTrk[j]) - pi) < phiOSmatch); 
      histSSTheta = (TMath::Abs(theta - thetaTrk[j]) < thetaSSmatch); 
      histOSTheta = (TMath::Abs(theta - thetaTrk[j]) < thetaOSmatch); 
      histPt      = (TMath::Abs(pt - momentumTrk[j]*sin(thetaTrk[j])) < ptmatch); 

      // -- in central theta: opposite sign loopers at opposite phi
      if ((charge*chargeTrk[j] < 0.)
	  && (theta > 1.4) && (theta < 1.7) 
	  && (pt < 0.18)
	  ) {

	// -- killing decision
	if ((charge*chargeTrk[j] < 0.) && (pt < 0.18)
	    && (TMath::Abs(TMath::Abs(phi - phiTrk[j]) - pi) < phiOSmatch)
	    && (TMath::Abs(theta - thetaTrk[j]) < thetaOSmatch)
	    && (TMath::Abs(pt - momentumTrk[j]*sin(thetaTrk[j])) < 0.1)
	    ) {
	  
	  if (TMath::Abs(zPocaTrk[j] - primVtxZ) < TMath::Abs(dz)) {
	    kill = i; 
	    osmult[j] += 1;
	  } else {
	    kill = j; 
	    osmult[i] += 1;
	  }
	  print += 1; 
	  if (osLooper) {
	    goodTrack[kill] = 0; 
	    loopTrack[kill] = 1;
	    continue; 
	  }
	}
      }


      // -- in central theta:  same sign loopers at the same phi
      if ((charge*chargeTrk[j] > 0.)
	  &&(theta < 1.7) && (theta > 1.4)
	  && (pt < 0.18)  
	  ) {

	// -- killing decision
	if ((charge*chargeTrk[j] > 0.)  && (pt < 0.18)
	    && (TMath::Abs(phi - phiTrk[j]) < phiSSmatch)
	    && (TMath::Abs(theta - thetaTrk[j]) < thetaSSmatch)
	    && (TMath::Abs(pt - momentumTrk[j]*sin(thetaTrk[j])) < 0.1)
	    ) {
	  
	  if (TMath::Abs(zPocaTrk[j] - primVtxZ) < TMath::Abs(dz)) {
	    kill = i; 
	    sssmult[j] += 1;
	  } else {
	    kill = j; 
	    sssmult[i] += 1;
	  }
	  print += 2; 
	  if (ssLooper) {
	    goodTrack[kill] = 0; 
	    loopTrack[kill] = 2;
	    continue; 
	  }
	}
      }


      // -- parallel tracks
      if ((charge*chargeTrk[j] > 0.)
	  && (pt > 0.)	 && (pt < 0.35)	  
	  ) {

	if ((charge*chargeTrk[j] > 0.) && (pt > 0.) && (pt < 0.35)
	    && (TMath::Abs(theta - thetaTrk[j]) < parathetamatch)
	    && (TMath::Abs(phi - phiTrk[j]) < paraphimatch)
	    && (TMath::Abs(pt - momentumTrk[j]*sin(thetaTrk[j])) < paraptmatch)
	    ) {

	  // -- killing decision
	  if (ndch < ndchTrk[j]) {
	    kill = i; 
	    ssfmult[j] += 1;
	  } else {
	    kill = j; 
	    ssfmult[i] += 1;
	  }
	  print += 4; 
	  if (ssParallel) {
	    goodTrack[kill] = 0; 
	    loopTrack[kill] = 3;
	  }
	}
      }

    }

  }

  if (doHist) {
    for (int j = 0; j < nTrk; ++j) {
      if (momentumTrk[j]*sin(thetaTrk[j]) >  0.18) continue; 

    }
  }

  if (0) { // print == 7) {
    for (int j = 0; j < nTrk; ++j) {
      cout << "  " << Form("%2d", j)
	   << Form("  pt:%+4.3f",chargeTrk[j]*momentumTrk[j]*sin(thetaTrk[j]))
	   << Form("  t: %4.3f", thetaTrk[j]) 
	   << Form("  f: %+4.3f", phiTrk[j]) 
	   << Form("  d: %3.0f", double(ndchTrk[j])) 
	   << Form("  s: %3.0f", double(nsvtTrk[j]))
	   << Form("  z: %+8.3f", zPocaTrk[j] - primVtxZ)
	   << Form("  z: %+8.3f", zPocaTrk[j] - beamSZ)
	   << Form("  d: %8.4f", TMath::Sqrt((xPocaTrk[j]-beamSX)*(xPocaTrk[j]-beamSX) + (yPocaTrk[j]-beamSY)*(yPocaTrk[j]-beamSY)))
	   << Form("  l: %5.4f", tLenTrk[j])
	   << endl;
    }
    
    cout << "event: " << fEvent << " in run " << fRunnumber << Form(":  %x/%x", fUpper, fLower) << endl;
    cout << "----------------------------------------------------------------------" << endl;
    
}

}


// ----------------------------------------------------------------------
void baseClass::recoil() {

  //this is a sample recoil method ... kinda like vub's recoil

  if (fVerbose) cout << "Starting with recoil()" << endl;

  fHistFile->cd("recoil");

  // -- Find lepton with highest p*
  TLorentzVector p4l; 
  Double_t mass(0.), pmax(0.), plab(0.), pcms(0.), lmass(0.), mNuSq(0.), mNuSqNC(0.);
  Int_t lmax(-1),  nB(0);
  double tmppurB(0.), tmpIpurB(0.);
    
  //Assign the correct Lund to the B reco
  fLeptonCharge = 0;
  Int_t tmpblund;
  tmpblund = B0LUND;
  nB = nB0;
  int modeB = modeB0[indexbestB];
  if(bestB0==0) {
    modeB = modeChB[indexbestB]; 
    tmpblund = CHBLUND;
    nB = nChB;
  }  
  tmpIpurB = brecointpur[modeB-10000];
  tmppurB = brecosig[modeB-10000]/(brecosig[modeB-10000]+brecobkg[modeB-10000]);
  
  fPurity = tmppurB;
  fIntPurity = tmpIpurB;

  if (tmppurB < 0.1) {
    if (fReturnLog[1] == 0) fReturnString[1] = TString(Form("purity too low"));
  }
  if (fVerbose) cout << "Survived cut on breco purity" << endl;

  if ((fSeedMode == 0) && (tmppurB < IPURDC)) {
    if (fReturnLog[2] == 0) fReturnString[2] = TString(Form("Dc, ipur too low"));
    fReturnLog[2]++;
    return;
  }

  if ((fSeedMode == 1) && (tmppurB < IPURDSTAR)) {
    if (fReturnLog[3] == 0) fReturnString[3] = TString(Form("D*, ipur too low"));
    fReturnLog[3]++;
    return;
  }
  if ((fSeedMode == 2) && (tmppurB < IPURD0)) {
    if (fReturnLog[4] == 0) fReturnString[4] = TString(Form("D0, ipur too low"));
    fReturnLog[4]++;
    return;
  }
  if ((fSeedMode == 3) && (tmppurB < IPURDSTAR0)) {
    if (fReturnLog[5] == 0) fReturnString[5] = TString(Form("D*0, ipur too low"));
    fReturnLog[5]++;
    return;
  }

  if (fVerbose) cout << "Survived cuts on breco integrated purities" << endl;

  fHistFile->cd("recoil");
  

  // -- Find leading lepton
  Int_t i;
  fNLepton = fNKp = fNlep = 0;
  for (i = 0; i < nTrk; ++i) {

    if (isTrkOfBreco(i)) continue;


    if (momentumTrk[i] < PLABLO) continue;
    if (goodTrack[i] == 0) continue;



    if (isRecLepton(i)) {  // do PID
      mass = MUMASS; 
      if (isRecEl(i)) mass = ELMASS; // electrons override muons
      mk4Vector(p4l, momentumTrk[i], thetaTrk[i], phiTrk[i], mass);
      p4l.Boost(-cmsBoost);
      plab = momentumTrk[i];
      pcms = p4l.Vect().Mag(); 

      if (isRecEl(i) && (pcms > ELMOMLO)) {
	fNlep++;
      } else if (isRecMu(i) && (pcms > MUMOMLO)) {
	fNlep++;
      }

      if (pcms > pmax) {
	pmax = pcms;
	lmax = i;
	lmass = mass;
	fLeptonCharge = chargeTrk[i];
      }
    }
  }

  fNLepton = fNlep;
  frealNLepton = fNLepton; 
  // -- Set up leading lepton
  fElectron = fMuon = kFALSE;
  fNEl = fNMu = 0;
  p4LeptonLab.SetXYZM(0., 0., 0., 0.);
  p4LeptonCms.SetXYZM(0., 0., 0., 0.);
  fIsPrompt = fIsCascade = kFALSE; 
  fPlab = fTlab = fTlabDR = fFlabDR = fPcms = fTcms = fFcms = fEcms = fGammaMax = -99.;
  if (lmax > -1) {
    if (fVerbose) cout << "Found a leading lepton" << endl;
    mk4Vector(p4LeptonLab, momentumTrk[lmax], thetaTrk[lmax], phiTrk[lmax], lmass);
    mk4Vector(p4LeptonCms, momentumTrk[lmax], thetaTrk[lmax], phiTrk[lmax], lmass);
    mk4Vector(p4LeptonUps, momentumTrk[lmax], thetaTrk[lmax], phiTrk[lmax], lmass);
    p4LeptonCms.Boost(-cmsBoost);
    p4LeptonUps.Boost(-upsBoost);
    fPlab   = p4LeptonLab.Vect().Mag();
    fTlab   = p4LeptonLab.Theta();
    fTlabDR = p4LeptonLab.Theta()*DR;
    fFlab   = p4LeptonLab.Phi();
    fFlabDR = p4LeptonLab.Phi()*DR;
    fPcms   = p4LeptonCms.Vect().Mag();
    fTcms   = p4LeptonCms.Theta();
    fFcms   = p4LeptonCms.Phi();
    fEcms   = p4LeptonCms.E();
    fPups   = p4LeptonUps.Vect().Mag();
    if (isRecEl(lmax)) {
      fElectron = kTRUE;
      fNEl = 1;
    }
    if (isRecMu(lmax)) {
      fMuon = kTRUE;
      fNMu = 1;
    }
    // -- remove overlap: misid'ed electrons (PidKilling) are automatically also misid'ed muons
    if (fElectron && fMuon) {
      fNMu = 0; 
      fMuon = kFALSE;
    }

    if (isPrompt(lmax)) {
      fIsPrompt = kTRUE; 
      fIsCascade = kFALSE; 
    } else if (isCascade(lmax)) {
      fIsPrompt = kFALSE; 
      fIsCascade = kTRUE; 
    }
  }

  fHistFile->cd();
  ((TH1D*)gDirectory->Get("deallevents"))->Fill(fDeltaE);
  ((TH1D*)gDirectory->Get("mesallevents"))->Fill(fMes);
  if (0 != fBrecoCharge) {
    ((TH1D*)gDirectory->Get("mesalleventsBch"))->Fill(fMes);
  } else {
    ((TH1D*)gDirectory->Get("mesalleventsBnu"))->Fill(fMes);
  }
  if (fSeedMode == 0) ((TH1D*)gDirectory->Get("mesalleventsS0"))->Fill(fMes);
  else if (fSeedMode == 1) ((TH1D*)gDirectory->Get("mesalleventsS1"))->Fill(fMes);
  else if (fSeedMode == 2) ((TH1D*)gDirectory->Get("mesalleventsS2"))->Fill(fMes);
  else if (fSeedMode == 3) ((TH1D*)gDirectory->Get("mesalleventsS3"))->Fill(fMes);



  if ((fPlab > 0.5) && (fPcms > 0.5)
      && (TLABLO < fTlab*DR) && (fTlab*DR < TLABHI)) {
  } else {
    fPcms = -96.; // reset fPcms so that the event is not dumped into 'events'
    if (fReturnLog[13] == 0) fReturnString[13] = TString("No Lepton in event");
    fReturnLog[13]++;
    return;
  }


  // -- Recoil calculation
  // ---------------------

  p4Xhad.SetXYZM(0., 0., 0., 0.);
  p4Recoil.SetXYZM(0., 0., 0., 0.);
  TLorentzVector p4t(0., 0., 0., 0.), p4BRecoilNC(0., 0., 0., 0.);

  // -- calculate recoil and Xhad (= recoil - hardest lepton)
  // --------------------------------------------------------

  // -- Tracks
  for (i = 0; i < nTrk; ++i) {
      if (isTrkOfBreco(i)) {
        goodTrack[i] = goodHadron[i] = goodChargedKaon[i] = goodPion[i] = 0; 
        continue;
      }


    fRecoilCharge += chargeTrk[i];
    ++fRecoilTrkMult;

    if (isRecEl(i)) {  
      if (i==lmax) {
        mk4Vector(p4t, momentumTrk[i], thetaTrk[i], phiTrk[i], ELMASS);
      }else{
        mk4Vector(p4t, momentumTrk[i], thetaTrk[i], phiTrk[i], PIPMASS);
      }

    } else if (isRecKaon(i) && (momentumTrk[i] > KAMOMLO)) {  
      // -- correct kaon tracks for wrong tracking hypothesis
      double momentum = p_energy_loss_corrected(momentumTrk[i], 0.5*TMath::Pi() - thetaTrk[i], 1);
      mk4Vector(p4t, momentum, thetaTrk[i], phiTrk[i], KAPMASS);
      momentumTrk[i] = momentum;
    } else if (isRecMu(i)) {
      mk4Vector(p4t, momentumTrk[i], thetaTrk[i], phiTrk[i], MUMASS);
    } else {
      mk4Vector(p4t, momentumTrk[i], thetaTrk[i], phiTrk[i], PIPMASS);
    }

    if (fVerbose) cout << i 
		       << "  " << momentumTrk[i] 
		       << "  " << thetaTrk[i] 
		       << "  " << phiTrk[i] 
		       << "  " << chargeTrk[i] 
		       << "  " << goodTrack[i]
		       << "  " << goodHadron[i]
		       << "  " << goodChargedKaon[i]
		       << "  " << goodPion[i]
		       << "  " << kshortLockTrk[i] 
		       << endl;

    p4Recoil += p4t;

    // -- Skip tracks which are leading leptons for the rest
    if (i == lmax) continue;
    p4Xhad += p4t;
    fcountChg++;
    p4ChargPart = p4Xhad;

  }
 // define mxhadchg: mxhad reconstructed with charged tracks only
  // (without the neutral ones):

  fPxhadchg   = p4Xhad.Vect().Mag(); 
  fTxhadchg   = p4Xhad.Theta(); 
  fFxhadchg   = p4Xhad.Phi(); 
  fExhadchg   = p4Xhad.E();
  fMxhadchg   = p4Xhad.Mag();


  if (fVerbose) cout << "== End track list in recoil() ==" << endl;
  frealNKp = fNKp; 

  // -- Photons
  int nEffGam(0); 

  fENeu=fEPiz=0;
  for (i = 0; i < nGam; ++i) {

    if (isGamOfBreco(i)) {
      goodPhoton[i] = 0; 
    }


    if (ecalGam[i] < 0) goodPhoton[i] = 0; 
    if (goodPhoton[i] == 0) continue;
       
    // test test rescale photon energy to get agreement ('the function from E. Maly)
    //    energyGam[i]=clusterCorrection(energyGam[i], thetaGam[i],true);
    //    if(fRunnumber>40000)energyGam[i]=1.09*energyGam[i];

    mk4Vector(p4t, energyGam[i], thetaGam[i], phiGam[i], 0.);
    if (energyGam[i] < 0.16) ++nEffGam;
 
    ++fRecoilNutMult;

    fENeu+=p4t.T();

    p4Xhad += p4t;

    p4Recoil += p4t;
    p4t.Boost(-cmsBoost);
    if (p4t.Vect().Mag() > fGammaMax) fGammaMax = p4t.Vect().Mag();
  }

  fPxhad   = p4Xhad.Vect().Mag(); 
  fTxhad   = p4Xhad.Theta(); 
  fFxhad   = p4Xhad.Phi(); 
  fExhad   = p4Xhad.E();
  fMxhad   = p4Xhad.Mag();

  fPrecoil   = p4Recoil.Vect().Mag(); 
  fTrecoil   = p4Recoil.Theta(); 
  fFrecoil   = p4Recoil.Phi(); 
  fErecoil   = p4Recoil.E();
  fMrecoil   = p4Recoil.Mag();
  
  TLorentzVector p4Allev = p4Breco + p4Recoil;

  fPAllev   = p4Allev.Vect().Mag(); 
  fTAllev   = p4Allev.Theta(); 
  fFAllev   = p4Allev.Phi(); 
  fEAllev   = p4Allev.E();
  fMAllev   = p4Allev.Mag();

  fQtot = fRecoilCharge + fBrecoCharge;  
  fNtracks = double(fRecoilTrkMult);
  fNneutrals = double(fRecoilNutMult);


  fHistFile->cd();


  p4Brecoil   = p4Upsilon - p4Breco;
  p4BRecoilNC = p4Upsilon - p4BrecoNC;
  TLorentzVector p4Neutrino = p4Upsilon - p4Breco - p4Recoil;
  TLorentzVector pfourneu[15]; 

  fMCharpart   = p4ChargPart.Mag();
  fTCharpart   = p4ChargPart.Theta(); 
  fPCharpart   = p4ChargPart.Phi(); 
  fECharpart   = p4ChargPart.E();

  fPNu   = p4Neutrino.Vect().Mag(); 
  fTNu   = p4Neutrino.Theta(); 
  fCosTNu= p4Neutrino.CosTheta(); 
  fFNu   = p4Neutrino.Phi(); 
  fMM2   = p4Neutrino.Mag2();
  fEmiss = p4Neutrino.E();
  fQ2    = 2.*p4Neutrino*p4LeptonLab;

  TLorentzVector p4NeutrinoNC = p4BRecoilNC - p4Xhad - p4LeptonLab;

  fPNuNC   = p4NeutrinoNC.Vect().Mag(); 
  fTNuNC   = p4NeutrinoNC.Theta(); 
  fFNuNC   = p4NeutrinoNC.Phi(); 
  fMM2NC   = p4NeutrinoNC.Mag2();
  fEmissNC = p4NeutrinoNC.E();
  fQ2NC    = p4NeutrinoNC*p4LeptonLab;

  mNuSq = p4Neutrino * p4Neutrino;
  mNuSqNC = p4NeutrinoNC * p4NeutrinoNC;

  
  fHistFile->cd();

  fQ2Fit = 0.;


  // -- Fill histograms
  // ------------------
  fMxhadRes    = fMxhad - fMxhadGen; 
  fMxhadfitRes = fMxhadfit - fMxhadGen; 
  fQ2Res       = fQ2 - fQ2Gen; 

  fHistFile->cd(); 

  // -- Determine Event weight
  // -------------------------
  if (fIsMC) calculateEvtW8(); 

  fHistFile->cd();
  
}

// ----------------------------------------------------------------------
void baseClass::dumpOneB(int b1) {
  char line[200];

  for (int i = 0; i < nMc; ++i) {
    if ((i == b1) || isAncestor(b1, i)) {
      sprintf(line, "%3d %+6d mom(%3d) ndau(%3d) p=%5.3f, t=%5.3f f=%+5.3f m = %+9.6f v=(%+7.3f,%+7.3f,%+7.3f)", 
	      i, idMc[i], mothMc[i]-1, nDauMc[i],
	      pMc[i], thetaMc[i], phiMc[i], massMc[i], 
	      xMc[i], yMc[i], zMc[i]);
      cout << line << endl;
    }
  }
}

// ----------------------------------------------------------------------
// prints overlap information if indices to the two B are passed
void baseClass::dumpGeneratorBlock(int b1, int b2) {
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
void baseClass::Loop(Int_t maxEvent, Int_t startEvent, Int_t isVerbose, Int_t lun) {

  int step(1000);
  findPro = findUps = 0;
  fVerbose = isVerbose; 
  double tmpMassPB, tmpMassThetaB, tmpMassPhiB ;
  double tmpPB, tmpThetaB, tmpPhiB ;
  double tmpPgen, tmpThetagen, tmpPhigen, tmpMassgen ;
  double tmpMB, tmpBevM;

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
  const char *pChar; 

  int oldrunnumber(0); 
  TString rfile(fHistFile->GetName()); cout << rfile<< endl;
  rfile.ReplaceAll(".root", ".runs"); 
  ofstream fRUN(rfile); 

  for (Int_t jentry = startEvent; jentry < startEvent+maxEvent; jentry++) {


    if (fReturnLog[0] == 0) fReturnString[0] = TString("Loop event counter");
    fReturnLog[0]++;
    fEvent = jentry;
    // in case of a TChain, ientry is the entry number in the current file
    Int_t ientry = LoadTree(jentry); 
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (fVerbose) cout << "->  new event " <<upperID<< " : " <<lowerID<< endl;
    if (jentry%step == 0) cout << "Event " << jentry << endl;
    if (ientry == 0) {
      cout << "File " << fChain->GetCurrentFile()->GetName(); 
      fFileChanged = 1; // file has changed
      if (pChar = strstr(fChain->GetCurrentFile()->GetName(), "2000")) {
	fRunRange = TString("Run 1"); 
	cout << " Run 1: 2000" << endl;
      } else if (pChar = strstr(fChain->GetCurrentFile()->GetName(), "2001")) {
	fRunRange = TString("Run 2a"); 
	cout << " Run 2a: 2001" << endl;
      } else if (pChar = strstr(fChain->GetCurrentFile()->GetName(), "2002")) {
	fRunRange = TString("Run 2b"); 
	cout << " Run 2b: 2002" << endl;
      }
      else {
	fRunRange = TString("undefined"); 
	cout << " Runrange ?" << endl;
      }
    } else {
      fFileChanged = 0;  // staying in same file
    }

    // -- Initialize event
    initVariables();

    if (fRunnumber != oldrunnumber) {
      fRUN << fRunnumber << endl;
      oldrunnumber = fRunnumber; 
    }

    findbestB();
    
    if (skipBadBreco()) {
      if (fReturnLog[6] == 0) fReturnString[6] = TString("Bad B-D flavor correlation");
      fReturnLog[6]++;
      continue;
    }

    fillPidMaps();
    if (fOptSmearTracks > 0) {
      smearTracks();
    }
    if (fOptSmearNeut > 0) {
      smearNeut();
    }


    int brecoI(-99);

    if (fIsMC) {
      mcTruth();
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
    if(fChB == 0) {
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

    p4Upsilon = TLorentzVector(pxUps, pyUps, pzUps, eUps); 
    upsBoost = TVector3(pxUps, pyUps, pzUps);
    upsBoost.SetMag(upsBoost.Mag()/eUps);

    p4Brecoil = p4Upsilon - p4Breco; 
    cmsBoost = p4Brecoil.BoostVector();

    // track and neutrals selection must be after the computation of the recoil B, which is used
    selectTracks();

    // splitoff studies --- needed for photon selection
    doSplitOffStudy();
    selectPhotons();

    TLorentzVector p4t = p4Breco; 
    p4t.Boost(-upsBoost);
    fPcmsBreco = p4t.Vect().Mag(); 

    if (fVerbose) cout << "CALLING breco()" << endl;
    breco();

    if (fVerbose) cout << "CALLING recoil()" << endl;
    recoil();

    if (fDump > 0) {
      if  (fDump & 1) {
	if (fOptGammas) fGTree->Fill();
	if (isVerbose)	cout<<"filling Gammas"<<endl;
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


    if(fVerbose) cout << " jentry = " << jentry << "  mes = " << fMes << " pcms = " << fPcms << " mass = " << p4LeptonLab.Mag()
		      << " xhadmass = " << fMxhad << " mm2 = " << fMM2NC << " xhadfitmass = " << fMxhadfit
		      << " vub = " << fVub << "  vcb = " << fVcb << endl;
  }

  cout << "----------------------------------------------------------------------" << endl;
  if (fOptGammas)  cout<<findPro<<" "<<findUps<<endl;
}

// ----------------------------------------------------------------------
TFile* baseClass::openHistFile(TString name) {
  cout << name << endl;
  fHistFile = new TFile(name.Data(), "RECREATE");
  fHistFile->cd();
  cout << "Opened " << fHistFile->GetName() << endl;
  return fHistFile;
}


// ----------------------------------------------------------------------
void baseClass::closeHistFile() {
  cout << "Return Log: " << endl;
  for (int i = 0; i < 100; ++i) {
    if (fReturnLog[i] > 0) cout << Form("fReturnLog[%3d] = %6d   %s", i, fReturnLog[i], (fReturnString[i]).Data()) << endl;
  }
  cout << "----------------------------------------------------------------------" << endl;

  cout << "Writing " << fHistFile->GetName() << endl;
  fHistFile->cd();
  fHistFile->Write();
  fHistFile->Close();
  delete fHistFile;

}


