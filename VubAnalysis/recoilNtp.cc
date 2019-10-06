#include "VubAnalysis/recoilNtp.hh"

#include <fstream.h>
#include <iostream.h>
#include <iomanip.h>
#include "TH2.h"
#include "TH1.h"
#include "TProfile.h"
#include <TMap.h>
#include <TExMap.h>
#include <TObjString.h>
#include "TLorentzVector.h"
#include "TRandom.h"


//-------------------------------------------------
//FORTRAN STUFF
//-------------------------------------------------
#if 1
extern "C" 
{
  // The fortran generator:
  int ILTYP;
  float CVAL[4];
  float P_REC[16];
  float P_FIT[16];
  float CHI2T, PROBCHI2;
  int   IERR, ISMEAR, ISV;
  int abcfit_interface_vub_(int *ISMEAR, int *ILTYP, float *CVAL, float *P_REC, float *P_FIT, float *CHI2T,float *PROBCHI2,int *IERR, int *ISV);
}
#endif


Double_t p_energy_loss_corrected(Double_t pin, Double_t dip, Int_t itype);

// ----------------------------------------------------------------------
recoilNtp::recoilNtp(TTree *tree, int isMC, int isNR, int newFormat) {
  fNewFormat = newFormat; 
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/u/ec/ursl/d/output/breco-110201/mc-b0cocktail-20.root");
    if (!f) {
      f = new TFile("/u/ec/ursl/d/output/breco-110201/mc-b0cocktail-20.root");
    }
    tree = (TTree*)gDirectory->Get("h9");
    
  }
  Init(tree,isMC,isNR);
  fToBeCopied = new TEventList("toBeCopied", "Events to be copied", 1000);
  //ADDED CB
  initRest();
}

// ----------------------------------------------------------------------
recoilNtp::~recoilNtp() {
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}


// ----------------------------------------------------------------------
recoilNtp::recoilNtp(TString filename, TString treename,int isMC, int isNR, int newFormat) {
  fNewFormat = newFormat;
  TFile *f = new TFile(filename);
  TTree *tree = (TTree*)f->Get(treename);
  if (!tree) { 
    cout << "Did not find " << treename << " in file " << filename << endl;
    f->ls();
  } else {
    Init(tree,isMC,isNR);
  }
  initRest();
}
#include "init.icc"
// ----------------------------------------------------------------------
Bool_t recoilNtp::isAncestor(int ancestor, int cand) {
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
Int_t recoilNtp::isRecoed(int imc) {
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
void  recoilNtp::printLorentz(const TLorentzVector &p) {
  char line[200]; 
  sprintf(line, "p = (%7.5f, %7.5f, %7.5f, %7.5f), m = %7.5f",  p.X(), p.Y(), p.Z(), p.E(), p.Mag()); 
  cout << line;
}


// ----------------------------------------------------------------------
void recoilNtp::testJpsi() {
  static Bool_t first(kTRUE);
  if (first) {
    first = kFALSE; 
    fHistFile->cd();
    fHistFile->mkdir("jpsi", "jpsi");
    fHistFile->cd("jpsi");
    TH1D *h;
    char name[100], title[100];
    sprintf(name, "mjpsi");  sprintf(title, "mass all Jpsi");  h = new TH1D(name, title, 100, 2.0, 4.0); 
  }    

  fHistFile->cd("jpsi");
  for (int i = 0; i < nJpsi; ++i) {
    double mass = massJpsi[i];
    int d1 = d1JpsiIndex[i]-1;
    int d2 = d2JpsiIndex[i]-1;
    ((TH1D*)gDirectory->Get("mjpsi"))->Fill(mass, 1.);
    cout << "Jpsi " << i << " mass: " << mass 
	 << " d1 = " << d1 << " p = " << momentumTrk[d1] << " e(" << (isRecEl(d1)? 1:0 ) << ")" << " m(" << (isRecMu(d1)? 1:0 ) << ")"
	 << " d2 = " << d2 << " p = " << momentumTrk[d2] << " e(" << (isRecEl(d2)? 1:0 ) << ")" << " m(" << (isRecMu(d2)? 1:0 ) << ")"
	 << endl;
  }
  return;
}

// ----------------------------------------------------------------------
void recoilNtp::testDalitz() {
  static Bool_t first(kTRUE);
  if (first) {
    first = kFALSE; 
    fHistFile->cd();
    fHistFile->mkdir("dalitz", "dalitz");
    fHistFile->cd("dalitz");
    TH1D *h;
    char name[100], title[100];
    sprintf(name, "mdalitz");  sprintf(title, "mass all Dalitz");  h = new TH1D(name, title, 100, 0., 0.1); 
  }    

  fHistFile->cd("dalitz");
  for (int i = 0; i < nDalitz; ++i) {
    double mass = massDalitz[i];
    int d1 = d1DalitzIndex[i]-1;
    int d2 = d2DalitzIndex[i]-1;
    ((TH1D*)gDirectory->Get("mdalitz"))->Fill(mass, 1.);
    cout << "Dalitz " << i << " mass: " << mass 
	 << " d1 = " << d1 << " p = " << momentumTrk[d1] << " e(" << (isRecEl(d1)? 1:0 ) << ")" << " m(" << (isRecMu(d1)? 1:0 ) << ")"
	 << " d2 = " << d2 << " p = " << momentumTrk[d2] << " e(" << (isRecEl(d2)? 1:0 ) << ")" << " m(" << (isRecMu(d2)? 1:0 ) << ")"
	 << endl;
  }
  return;
}
  
// ----------------------------------------------------------------------
void recoilNtp::testConversion() {
  static Bool_t first(kTRUE);
  if (first) {
    first = kFALSE; 
    fHistFile->cd();
    fHistFile->mkdir("conversion", "conversion");
    fHistFile->cd("conversion");
    TH1D *h;
    char name[100], title[100];
    sprintf(name, "mconversion");  sprintf(title, "mass all Conversions");  h = new TH1D(name, title, 100, 0., 0.1); 
  }    

  fHistFile->cd("conversion");
  for (int i = 0; i < nGConv; ++i) {
    double mass = massGConv[i];
    int d1 = d1GConvIndex[i]-1;
    int d2 = d2GConvIndex[i]-1;
    ((TH1D*)gDirectory->Get("mconversion"))->Fill(mass, 1.);
    cout << "Conversion " << i << " mass: " << mass 
	 << " d1 = " << d1 << " p = " << momentumTrk[d1] << " e(" << (isRecEl(d1)? 1:0 ) << ")" << " m(" << (isRecMu(d1)? 1:0 ) << ")"
	 << " d2 = " << d2 << " p = " << momentumTrk[d2] << " e(" << (isRecEl(d2)? 1:0 ) << ")" << " m(" << (isRecMu(d2)? 1:0 ) << ")"
	 << endl;
  }
  return;
}

// ----------------------------------------------------------------------
void recoilNtp::maskPi0(int modes) {
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
#include "util.icc"

void recoilNtp::timestamp(const char * lun) {

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
    temp = (Int_t) (pow(2,32) + partition);
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

  //  cout<<name<<" "<<name2<<endl;

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
Int_t recoilNtp::bestKsIndex(Int_t isVerbose)
{
  fm0ks = 0;
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
	fm0ks = massKs[i];
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

#include "gamStudy.icc"
#include "trackStudy.icc"
// ----------------------------------------------------------------------
void recoilNtp::recoil(int chbcand) {

  if (fVerbose) cout << "Starting with recoil()" << endl;

  fHistFile->cd("recoil");

  Int_t Brectrktmp, Brecgamtmp;
  // -- Find lepton with highest p*
  TLorentzVector p4l; 
  Double_t mass(0.), pmax(0.), plab(0.), pcms(0.), lmass(0.), mNuSq(0.), mNuSqNC(0.);
  Int_t lmax(-1),  nB(0);
  double tmppurB(0.), tmpIpurB(0.);
  totweightfRecoilTrkMult = totweightfRecoilNutMult = 1;
  totweight = 1;
    
  //Assign the correct Lund to the B reco
  fLeptonCharge = 0;
  Int_t tmpblund;
  tmpblund = B0LUND;
  nB = nB0;
//   tmpIpurB = intpurB0[indexbestB]; //old definition
//   tmppurB = purB0[indexbestB];
  int modeB = modeB0[indexbestB];
  if(bestB0==0) {
    modeB = modeChB[indexbestB]; 
    tmpblund = CHBLUND;
    nB = nChB;
//     tmpIpurB = intpurChB[indexbestB]; //old definition
//     tmppurB = purChB[indexbestB];
  }  
  tmpIpurB = brecointpur[modeB-10000];
  tmppurB = brecosig[modeB-10000]/(brecosig[modeB-10000]+brecobkg[modeB-10000]);
  
  fPurity = tmppurB;
  fIntPurity = tmpIpurB;

  //  cout << tmpIpurB << "  " << tmppurB << endl;

  if (tmppurB < 0.1) {
    if (fReturnLog[1] == 0) fReturnString[1] = TString(Form("purity too low"));
    //fReturnLog[1]++;
    //    return;
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

#ifndef FAST
  for (int itrk = 0 ; itrk < nTrk; ++itrk) {
    weightTrk[itrk] = getweightchg(momentumTrk[itrk]);
  }

  for (int i = 0; i < nGam; ++i) {
    weightNeu[i] = getweightneu(energyGam[i]);
  }
#endif

  fHistFile->cd("recoil");
  
    if(fIsMC){
      for (int i = 0; i < nTrk; ++i) {
    	Brectrktmp = B0RecTrk[i];
	if(chbcand) Brectrktmp = chBRecTrk[i];
	if ((Brectrktmp&brecoOverlap)&&isAncestor(fBVxb,IndexTrk[i]-1)){
 	  fBadReco++;
  	}
      }      
    }


  // -- Find leading lepton
  Int_t i, jbit;
  fNLepton = fNKshort = fNKp = fNlep = 0;
  for (i = 0; i < nTrk; ++i) {
    Brectrktmp = B0RecTrk[i];
    if(chbcand) Brectrktmp = chBRecTrk[i];
    if ((Brectrktmp&brecoOverlap))      continue;
    if(fNewFormat == 1&& Brectrktmp==2)continue;// Apr 02 format

    //    cout << "recoil track " << i << "  " << momentumTrk[i] << " " << goodTrack[i] << endl; 

    if (momentumTrk[i] < PLABLO) continue;
    if (goodTrack[i] == 0) continue;

#ifndef FAST
    // -- PID selector bits
    for (jbit = 0; jbit < 10; ++jbit) {
      if (elecIdTrk[i] & (0x1<<jbit)) ((TH1D*)gDirectory->Get("elSelBits"))->Fill(jbit);
      if (muonIdTrk[i] & (0x1<<jbit)) ((TH1D*)gDirectory->Get("muSelBits"))->Fill(jbit);
      if (kaonIdTrk[i] & (0x1<<jbit)) ((TH1D*)gDirectory->Get("kaSelBits"))->Fill(jbit);
    }
#endif

    if (isRecLepton(i)) {  // do PID
      mass = MUMASS; 
      if (isRecEl(i)) {
	mass = ELMASS; // electrons override muons
	if(DOBREMRECOVERY) doBrem(i);
      }
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
  fNEl = fNMu = flepmap = 0;
  p4LeptonLab.SetXYZM(0., 0., 0., 0.);
  p4LeptonCms.SetXYZM(0., 0., 0., 0.);
  fIsPrompt = fIsCascade = kFALSE; 
  fPlab = fTlab = fTlabDR = fFlabDR = fPcms = fTcms = fFcms = fEcms = fGammaMax = -99.;
  if (lmax > -1) {
    if (fVerbose) cout << "Found a leading lepton with p = " << pcms << " at index " << lmax << endl;
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
    ftruelep    = 0;
    if (fIsMC&&fPlab>0&&TMath::Abs(fTlab-thetaMc[IndexTrk[lmax]-1])<5./180.*3.14
        &&(TMath::Abs(fFlab-phiMc[IndexTrk[lmax]-1])<5./180.*3.14
	||TMath::Abs(TMath::Abs(fFlab-phiMc[IndexTrk[lmax]-1])-2*3.1415)<5./180.*3.14)){
      ftruelep    = idMc[IndexTrk[lmax]-1];
    }
    fmothlep=0;
    if (TMath::Abs(ftruelep)==11||TMath::Abs(ftruelep)==13){
      fmothlep=idMc[mothMc[IndexTrk[lmax]-1]-1];
    // if(TMath::Abs(fmothlep)==511||TMath::Abs(fmothlep)==521)fmothlep=1;
      // else fmothlep=2;
    }
    if (fIsMC&&isAncestor(brecoI,IndexTrk[lmax]-1))fws=1;
    if (isRecEl(lmax)) {
      fElectron = kTRUE;
      fNEl = 1;
      flepmap|=1;
    }
    if (isRecMu(lmax)) {
      fMuon = kTRUE;
      fNMu = 1;
      flepmap|=4;
    }
    // -- remove overlap: misid'ed electrons (PidKilling) are automatically also misid'ed muons
    if (fElectron && fMuon) {
      fNMu = 0; 
      fMuon = kFALSE;
      //      cout << "!@#%$^& overlap between el and mu" << endl;
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
    // this is ok
  } else {
    fPcms = -96.; // reset fPcms so that the event is not dumped into 'events'
    if (fReturnLog[13] == 0) fReturnString[13] = TString("No Lepton in event");
    fReturnLog[13]++;
    return;
  }


  // -- Recoil calculation
  // ---------------------

  fEffCat = 0;
  fRecoilPi0Mult=fRecoilTrkMult = fRecoilNutMult = fRecoilNutMult80_160 = fRecoilNutMult160_320 = fRecoilNutMultfromB = fRecoilNutMultfromB80_160 = fRecoilNutMultfromB160_320=fRecoilCharge = 0; 

  vubDepleted = kFALSE;
  p4Xhad.SetXYZM(0., 0., 0., 0.);
  p4Recoil.SetXYZM(0., 0., 0., 0.);
  TLorentzVector recoillep;
  recoillep.SetXYZM(0., 0., 0., 0.);    
  int chargelep = 0;
  TLorentzVector p4t(0., 0., 0., 0.), p4BRecoilNC(0., 0., 0., 0.);

  if (fVerbose) {  
    cout << "----------------------------------------------------------------------" << endl;
    cout << "nKs = " << nKs << endl;
    cout << "-boost = "; printLorentz(TLorentzVector(-cmsBoost, 0.)); cout << endl;
    cout << "+boost = "; printLorentz(TLorentzVector(cmsBoost, 0.)); cout << endl;
  }

  // -- Kshorts
  for (i = 0; i < nKs; ++i) {
    if (goodKshort[i] != 1) continue;   // skip suboptimal kshorts 
    ++fNKshort;
    vubDepleted = kTRUE;
  }
  fcountChg=0;
  frealNKshort = fNKshort; 

  // -- calculate recoil and Xhad (= recoil - hardest lepton)
  // --------------------------------------------------------

  // -- Tracks
  if (fVerbose) cout << "== Start track list in recoil() ==" << endl;
  nGoodPi = 0;
  fMinKMom = 10000.;
  fMaxKMom = -1.;
  for (i = 0; i < nTrk; ++i) {
    int ispion(0), iskaon(0);
    if ((goodTrack[i] == 0) && (kshortLockTrk[i] == 0)) continue; // skip tracks only if they are not part of a kshort
    Brectrktmp = B0RecTrk[i];
    if(chbcand) Brectrktmp = chBRecTrk[i];
    if ((Brectrktmp&brecoOverlap)) {
      goodTrack[i] = goodHadron[i] = goodChargedKaon[i] = goodPion[i] = 0; 
      continue;
    }
    if (fNewFormat == 1&& Brectrktmp==2) {
      goodTrack[i] = goodHadron[i] = goodChargedKaon[i] = goodPion[i] = 0; 
      continue; // Apr 02 format
    }

    fRecoilCharge += chargeTrk[i];
    ++fRecoilTrkMult;
    totweightfRecoilTrkMult *= weightTrk[i];

    if (isRecEl(i)) {  
      goodHadron[i] = goodChargedKaon[i] = goodPion[i] = 0; 
      if (i==lmax) {
        mk4Vector(p4t, momentumTrk[i], thetaTrk[i], phiTrk[i], ELMASS);
      }else{
        mk4Vector(p4t, momentumTrk[i], thetaTrk[i], phiTrk[i], PIPMASS);
	ispion = 1;
	flepmap|=2;
      }
      if ((i == lmax) && (convLockTrk[i] == 1)) { 
	cout << "... -> Leading electron flagged as converted photon" << endl; 
      }
    } else if (isRecKaon(i) && (momentumTrk[i] > KAMOMLO) && i != lmax) {  
      //test    } else if (isRecKaon(i) && (momentumTrk[i] > KAMOMLO) && (chargeTrk[i]*fLeptonCharge > 0)) {  
      goodPion[i] = 0; 
      goodHadron[i] *= 1;      // Also needs to fulfill track requirements, set in selectTracks()
      goodChargedKaon[i] *= 1; // Also needs to fulfill track requirements, set in selectTracks()
      vubDepleted = kTRUE;
      ++fNKp;
      // -- correct kaon tracks for wrong tracking hypothesis
      double momentum = p_energy_loss_corrected(momentumTrk[i], 0.5*TMath::Pi() - thetaTrk[i], 1);
      ((TH2D*)fHistFile->Get("kMomCorr"))->Fill(momentumTrk[i], momentum);  
      mk4Vector(p4t, momentum, thetaTrk[i], phiTrk[i], KAPMASS);
      momentumTrk[i] = momentum;
      if(momentum<fMinKMom)fMinKMom=momentum;
      if(momentum>fMaxKMom)fMaxKMom=momentum;
      iskaon = 1;
    } else if (isRecMu(i)) {
      goodHadron[i] = goodChargedKaon[i] = goodPion[i] = 0; 
      if(i==lmax) {
        mk4Vector(p4t, momentumTrk[i], thetaTrk[i], phiTrk[i], MUMASS);
      }else{
        mk4Vector(p4t, momentumTrk[i], thetaTrk[i], phiTrk[i], PIPMASS);
	ispion = 1;
	flepmap|=8;
      }
    } else {
      goodChargedKaon[i] = 0; 
      goodHadron[i] *= 1;  // Also needs to fulfill track requirements, set in selectTracks()
      goodPion[i] *= 1;    // Also needs to fulfill track requirements, set in selectTracks()
      mk4Vector(p4t, momentumTrk[i], thetaTrk[i], phiTrk[i], PIPMASS);
      ispion = 1;
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
    if (i == lmax) {
      recoillep += p4t;    
      chargelep = chargeTrk[i];
    }

    // -- Skip tracks which are leading leptons for the rest
    if (i == lmax) continue;
    p4Xhad += p4t;
    fcountChg++;
    p4ChargPart = p4Xhad;

    if (ispion == 1 || iskaon == 1) {
      momentumGoodPi[nGoodPi] = momentumTrk[i];
      thetaGoodPi[nGoodPi] = thetaTrk[i];
      phiGoodPi[nGoodPi] = phiTrk[i];
      chargeGoodPi[nGoodPi] = chargeTrk[i];     
      isKaonGoodPi[nGoodPi] = iskaon;
      nGoodPi ++;   
    }
    fNpi = nGoodPi;


  }
 // define mxhadchg: mxhad reconstructed with charged tracks only
  // (without the neutral ones):




  fPxhadchg   = p4Xhad.Vect().Mag(); 
  fTxhadchg   = p4Xhad.Theta(); 
  fFxhadchg   = p4Xhad.Phi(); 
  fExhadchg   = p4Xhad.E();
  fMxhadchg   = p4Xhad.Mag();


  if (TMath::Abs(fMxhadchg - 0.5) < 0.05) {
    cout << "###################################################################### " << endl;
    cout << " THIS IS SUCH AN EVENT" << endl;
    cout << "###################################################################### " << endl;
  }


  if (fVerbose) cout << "== End track list in recoil() ==" << endl;
  frealNKp = fNKp; 

#ifndef FAST
  //Gammas study
  if (fOptGammas) {
    mcGam();
    GamStudy(p4Xhad);
  }
#endif

  // -- Photons
  int nEffGam(0); 

#ifndef FAST
  // count pi0s on the recoil
  for (i = 0; i < nPi0; ++i) {
    int ig=d1Pi0Index[i];
    if(ig>0){
      Brecgamtmp = B0RecGam[ig];
      if(chbcand) Brecgamtmp = chBRecGam[ig];
      if ((Brecgamtmp&brecoOverlap)) continue;
      if(fNewFormat == 1&& Brecgamtmp==2)continue;// Apr 02 format
    }
    ig=d2Pi0Index[i];
    if(ig>0){
      Brecgamtmp = B0RecGam[ig];
      if(chbcand) Brecgamtmp = chBRecGam[ig];
      if ((Brecgamtmp&brecoOverlap)) continue;
      if(fNewFormat == 1&& Brecgamtmp==2)continue;// Apr 02 format
    }
    fRecoilPi0Mult++;
  }
#endif

  p4NuPart.SetXYZM(0., 0., 0., 0.);
  fExhadnu = fFxhadnu = fTxhadnu = fPxhadnu = fMxhadnu = -99.;

  fENeu=fEPiz=0;
  nGoodGam = fp4XminPi = 0;
  fphotdeltaM = -11111;
  fp4XminPi= -11111;
  static int klprint(0);
  if ((fOptScaleKlongs == 1) && (klprint == 0)) {
    klprint = 1; 
    cout << "smearing KLONG energies with s = " << 0.274351/0.225545 << endl;
  }

  for (i = 0; i < nGam; ++i) {
    Brecgamtmp = B0RecGam[i];
    if(chbcand) Brecgamtmp = chBRecGam[i];
    if (Brecgamtmp&brecoOverlap)  goodPhoton[i] = 0; 
    if (fNewFormat == 1 && Brecgamtmp==2) goodPhoton[i] = 0;  // Apr 02 format
    if (ecalGam[i] < 0) goodPhoton[i] = 0; 
    if (goodPhoton[i] == 0) continue;
    //Reweighting of KL energy.
    //    energyGam[i] *= correctionFactor;
    if(KLlikeEMC(i)>fphotdeltaM){
      fphotdeltaM=KLlikeEMC(i);
      fp4XminPi = energyGam[i];
      fBestKlind = idGam[i];
    }
    
    // test test rescale photon energy to get agreement ('the function from E. Maly)
    //    energyGam[i]=clusterCorrection(energyGam[i], thetaGam[i],true);
    //    if(fRunnumber>40000)energyGam[i]=1.09*energyGam[i];

    // -- Alessio's magic wand for KL systematics is replaced by 
    // -- The extreme approach
    if (fOptScaleKlongs == 1) {
      if( idGam[i] == 130 ){
        energyGam[i] = 0.;
      }
    }

    mk4Vector(p4t, energyGam[i], thetaGam[i], phiGam[i], 0.);
    if (energyGam[i] < 0.16) ++nEffGam;
 
    ++fRecoilNutMult;
    totweightfRecoilNutMult *= weightNeu[i];

    if (energyGam[i] <= 0.160) ++fRecoilNutMult80_160;
    if (energyGam[i] >= 0.160 && energyGam[i] >= 0.320) ++fRecoilNutMult160_320;

    if (ifromBGam[i] == 1) ++fRecoilNutMultfromB; 
    if (energyGam[i] <= 0.160 && ifromBGam[i] == 1) ++fRecoilNutMultfromB80_160;
    if (energyGam[i] >= 0.160 && energyGam[i] >= 0.320 && ifromBGam[i] == 1) ++fRecoilNutMultfromB160_320;  
    bool pi0Dau(false);
    for (int ip=0 ; ip< nPi0 ; ip++) { 
      int indg1=d1Pi0Index[ip]-1;
      int indg2=d2Pi0Index[ip]-1;
      if(indg1>0){
	if(i == indg1 || i == indg2 ) {
	  pi0Dau = true;
	  break;
	}
	//      } else {
	//if(
      }
    }
    if(pi0Dau)fEPiz+=p4t.T();

    fENeu+=p4t.T();
    p4NuPart += p4t;
    p4Xhad += p4t;

    if (fVerbose) { cout << "  -> gam   = " << " igam = " << i << " ";  printLorentz(p4t); cout  << endl;}

    momentumGoodGam[nGoodGam] = energyGam[i];
    thetaGoodGam[nGoodGam] = thetaGam[i];
    phiGoodGam[nGoodGam] = phiGam[i];
    nGoodGam ++;   

    p4Recoil += p4t;
    p4t.Boost(-cmsBoost);
    if (p4t.Vect().Mag() > fGammaMax) fGammaMax = p4t.Vect().Mag();
  }

  // alternative photon definition
  fnaltphot=0;fEneualt=0;
  for (i = 0; i < nGam; ++i) {
    Brecgamtmp = B0RecGam[i];
    if(chbcand) Brecgamtmp = chBRecGam[i];
    if (Brecgamtmp&brecoOverlap)  altPhoton[i] = 0; 
    if (fNewFormat == 1 && Brecgamtmp==2) altPhoton[i] = 0;  // Apr 02 format
    if (ecalGam[i] < 0) altPhoton[i] = 0; 
    if (altPhoton[i] == 0) continue;
    fnaltphot++;
    fEneualt+=energyGam[i];
  }
  // =======================================================================
  // EXCLUSIVE STUFF

  fHistFile->cd();

  // pi+ calculation
  TLorentzVector bestpi;
  bestpi.SetXYZM(0., 0., 0., 0.);     
  double bestmm2pi = 1000.;
  double bestmtrk=0;
  double bestpnupi=0;
  for (i = 0; i < nGoodPi; ++i) {
    if(chargeGoodPi[i]==chargelep) continue;
    TLorentzVector tmpX = recoillep;
    TLorentzVector tmppi;
    mk4Vector(tmppi, momentumGoodPi[i], thetaGoodPi[i], phiGoodPi[i], .13957); 
    tmpX += tmppi;
    TLorentzVector tmpNeutrino = p4Upsilon - p4Breco - tmpX;
    double tmpMM2pi = tmpNeutrino.Mag2();
    if (TMath::Abs(tmpMM2pi)<TMath::Abs(bestmm2pi)) {
      bestmm2pi = tmpMM2pi;
      bestpi = tmppi;
      bestpnupi = tmpNeutrino.Vect().Mag();
      if(fNLepton>1){
        if(fNEl>0){
          mk4Vector(tmppi, momentumGoodPi[i], thetaGoodPi[i], phiGoodPi[i],ELMASS);
	}else{
	  mk4Vector(tmppi, momentumGoodPi[i], thetaGoodPi[i], phiGoodPi[i],MUMASS);
	}
      }
      bestmtrk= (recoillep+tmppi).Mag();
      }
  }
  fMM2pi = bestmm2pi;
  fmompi = bestpi.Rho();
  fthpi =  bestpi.Theta();
  fphipi = bestpi.Phi();

  fmtrkpi = bestmtrk;
  fpnupi = bestpnupi;

  // pi0 and eta -> gg calculation

  TLorentzVector bestpi0;
  bestpi0.SetXYZM(0., 0., 0., 0.);     
  TLorentzVector bestetagg;
  bestetagg.SetXYZM(0., 0., 0., 0.);     
  TLorentzVector best1phpi0;
  best1phpi0.SetXYZM(0., 0., 0., 0.);     
  TLorentzVector best2phpi0;
  best2phpi0.SetXYZM(0., 0., 0., 0.);     
  int countpi0 = 0;
  int countetagg = 0;
  double bestmm2 = 1000.;
  double bestmm2etagg = 1000.;
  int idxBestEta = -999;
  nGoodPi0 = 0;
  nGoodEta = 0;
  nGoodRho0 = 0;
  for (i = 0; i < nGoodGam; ++i) {
    for (int j = 0; j < nGoodGam; ++j) {
      if(j < i){
	TLorentzVector tmpph1;
	TLorentzVector tmpph2;
	TLorentzVector tmppi0;
	TLorentzVector tmpX = recoillep;
	tmppi0.SetXYZM(0., 0., 0., 0.);     
	mk4Vector(tmpph1, momentumGoodGam[i], thetaGoodGam[i], phiGoodGam[i], 0.);      
	tmppi0 += tmpph1;
	mk4Vector(tmpph2, momentumGoodGam[j], thetaGoodGam[j], phiGoodGam[j], 0.);              
	tmppi0 += tmpph2;
	tmpX += tmppi0;
	double pi0mass = tmppi0.Mag();
	((TH1D*)gDirectory->Get("pi0massnocutall"))->Fill(pi0mass);	   
	if(fBVxbTyp==-11) ((TH1D*)gDirectory->Get("pi0massnocut"))->Fill(pi0mass);	   
	double tmpMM2pi0;
	double tmpMM2eta;
	if(pi0mass>.11&&pi0mass<.16 && nGoodPi0<100) {
	  countpi0 ++;
	  TLorentzVector tmpNeutrino = p4Upsilon - p4Breco - tmpX;
	  tmpMM2pi0   = tmpNeutrino.Mag2();
	  gammafrompi0[i] = 1;
	  gammafrompi0[j] = 1;

	  if (TMath::Abs(tmpMM2pi0)<TMath::Abs(bestmm2)) {
	    bestmm2 = tmpMM2pi0;
	    bestpi0 = tmppi0;
	    best1phpi0 = tmpph1;
	    best2phpi0 = tmpph2;
	    best1phpi0.Boost(-cmsBoost);
	    fmom1phpi0 = best1phpi0.Vect().Mag();
	    best2phpi0.Boost(-cmsBoost);
	    fmom2phpi0 = best2phpi0.Vect().Mag();
	    if(fmom2phpi0>fmom1phpi0) {
	      double tmpmom = fmom1phpi0;
	      fmom1phpi0 = fmom2phpi0;
	      fmom2phpi0 = tmpmom;
	    }
	  }
	  momentumGoodPi0[nGoodPi0] = tmppi0.Rho();
	  thetaGoodPi0[nGoodPi0] = tmppi0.Theta();
	  phiGoodPi0[nGoodPi0] = tmppi0.Phi();
	  massGoodPi0[nGoodPi0] = tmppi0.Mag();
	  pi0dau1index[nGoodPi0] = i;
	  pi0dau2index[nGoodPi0] = j;
	  nGoodPi0 ++;   
	}
	if(pi0mass>.45&&pi0mass<.65&&nGoodEta<80) {
	  countetagg ++;
	  TLorentzVector tmpNeutrino = p4Upsilon - p4Breco - tmpX;
	  tmpMM2eta   = tmpNeutrino.Mag2();	  
	  if (TMath::Abs(tmpMM2eta)<TMath::Abs(bestmm2etagg)) {
	    bestmm2etagg = tmpMM2eta;
	    bestetagg = tmppi0;
            idxBestEta = nGoodEta;
	  }

	  momentumGoodEta[nGoodEta] = tmppi0.Rho();
	  thetaGoodEta[nGoodEta] = tmppi0.Theta();
	  phiGoodEta[nGoodEta] = tmppi0.Phi();
	  massGoodEta[nGoodEta] = tmppi0.Mag();
	  etadau1index[nGoodEta] = i;
	  etadau2index[nGoodEta] = j;	  
	  etadau1lund[nGoodEta] = 22;
	  etadau2lund[nGoodEta] = 22;	  
	  nGoodEta ++;   
	}	
      }
    }
  }

  fNrecopi0 = countpi0;
  double bestpi0mass = bestpi0.Mag();
  fmompi0 = bestpi0.Rho();
  fthpi0 =  bestpi0.Theta();
  fphipi0 = bestpi0.Phi();
  double bestetaggmass = bestetagg.Mag();
  if (countpi0){
     ((TH1D*)gDirectory->Get("pi0massall"))->Fill(bestpi0mass);
     if(bestpi0mass>.11&&bestpi0mass<.16) ((TH1D*)gDirectory->Get("mm2pi0all"))->Fill(bestmm2);	  
     if(fBVxbTyp==-11){
       ((TH1D*)gDirectory->Get("pi0mass"))->Fill(bestpi0mass);
       ((TH1D*)gDirectory->Get("mm2pi0"))->Fill(bestmm2);	  
     }
     fMM2pi0 = bestmm2;
     fMpi0 = bestpi0mass;
  }
  if (countetagg){
     fMM2etagg = bestmm2etagg;
     fMetagg = bestetaggmass;
  }
  
  // if only one photon is present, merged pi0 is assumed

  if (nGoodGam == 1){
    TLorentzVector tmppi0;
    mk4Vector(tmppi0, momentumGoodGam[0], thetaGoodGam[0], phiGoodGam[0], .13497);      
    TLorentzVector tmpX = recoillep;
    tmpX += tmppi0;
    TLorentzVector tmpNeutrino = p4Upsilon - p4Breco - tmpX;
    double tmpMM2   = tmpNeutrino.Mag2();
    ((TH1D*)gDirectory->Get("mm2pi0_0pi0all"))->Fill(tmpMM2);	  
    if(fBVxbTyp==-11) ((TH1D*)gDirectory->Get("mm2pi0_0pi0"))->Fill(tmpMM2);	 
    fMM2gamma = tmpMM2;
    bestpi0 = tmppi0;
    countpi0 ++;
    TLorentzVector tmpph;
    mk4Vector(tmpph, momentumGoodGam[0], thetaGoodGam[0], phiGoodGam[0], 0.);      
    best1phpi0 = tmpph;
    best1phpi0.Boost(-cmsBoost);
    fmom1phpi0 = best1phpi0.Vect().Mag();    
    momentumGoodPi0[nGoodPi0] = tmppi0.Rho();
    thetaGoodPi0[nGoodPi0] = tmppi0.Theta();
    phiGoodPi0[nGoodPi0] =  tmppi0.Phi();
    massGoodPi0[nGoodPi0] = tmppi0.Mag();
    pi0dau1index[nGoodPi0] = 0;
    nGoodPi0 ++;        
  }  
  
  for (int k = 0; k < nGoodPi0; ++k) {
    massGoodPi0[k] = .13497;
  }

  // rho0, omega and eta -> pipipi0 calculation

  int countrho0(0);
  int countomega(0);
  int countetapppi0(0);
  TLorentzVector best1pirho0;
  best1pirho0.SetXYZM(0., 0., 0., 0.);     
  TLorentzVector best2pirho0;
  best2pirho0.SetXYZM(0., 0., 0., 0.);     
  TLorentzVector best1piome;
  best1piome.SetXYZM(0., 0., 0., 0.);     
  TLorentzVector best2piome;
  best2piome.SetXYZM(0., 0., 0., 0.);     
  TLorentzVector bestpi0ome;
  bestpi0ome.SetXYZM(0., 0., 0., 0.);     
  TLorentzVector bestrho0;
  bestrho0.SetXYZM(0., 0., 0., 0.);     
  TLorentzVector bestomega;
  bestomega.SetXYZM(0., 0., 0., 0.);     
  TLorentzVector bestetapppi0;
  bestetapppi0.SetXYZM(0., 0., 0., 0.);     
  double bestmm2rho0 = 1000.;
  double bestmm2etapppi0 = 1000.;
  double bestmm2omega = 1000.;
  for (i = 0; i < nGoodPi; ++i) {
    for (int j = 0; j < nGoodPi; ++j) {
      if(j < i){
	if(chargeGoodPi[i]!=chargeGoodPi[j]){
	  TLorentzVector tmpph1;
	  TLorentzVector tmpph2;
	  TLorentzVector tmprho0;
	  TLorentzVector tmpX = recoillep;
	  tmprho0.SetXYZM(0., 0., 0., 0.);     
	  mk4Vector(tmpph1, momentumGoodPi[i], thetaGoodPi[i], phiGoodPi[i], .13957);      
	  tmprho0 += tmpph1;
	  mk4Vector(tmpph2, momentumGoodPi[j], thetaGoodPi[j], phiGoodPi[j], .13957);   
	  tmprho0 += tmpph2;
	  tmpX += tmprho0;
	  double rho0mass = tmprho0.Mag();
	  ((TH1D*)gDirectory->Get("rho0massnocutall"))->Fill(rho0mass);	   
	  if(fBVxbTyp==-13) ((TH1D*)gDirectory->Get("rho0massnocut"))->Fill(rho0mass);    
	  TLorentzVector tmpNeutrino = p4Upsilon - p4Breco - tmpX;
	  double tmpMM2rho0  = tmpNeutrino.Mag2();
	  countrho0 ++;
	  if (rho0mass<1.5&&TMath::Abs(tmpMM2rho0)<TMath::Abs(bestmm2rho0)) {
	    bestmm2rho0 = tmpMM2rho0;
	    fMrho0= rho0mass;
	    fMM2rho0 = bestmm2rho0;
	    bestrho0 = tmprho0;
	    best1pirho0 = tmpph1;
	    best2pirho0 = tmpph2;	      
	    best1pirho0.Boost(-cmsBoost);
	    fmom1pirho0 = best1pirho0.Vect().Mag();
	    best2pirho0.Boost(-cmsBoost);
	    fmom2pirho0 = best2pirho0.Vect().Mag();
	    if(fmom2pirho0>fmom1pirho0) {
	      double tmpmom = fmom1pirho0;
	      fmom1pirho0 = fmom2pirho0;
	      fmom2pirho0 = tmpmom;
	    }
	  }
	  if(nGoodRho0<80){
	    momentumGoodRho0[nGoodRho0] = tmprho0.Rho();
	    thetaGoodRho0[nGoodRho0] = tmprho0.Theta();
	    phiGoodRho0[nGoodRho0] = tmprho0.Phi();
	    massGoodRho0[nGoodRho0] = tmprho0.Mag();
	    rho0dau1index[nGoodRho0] = i;
	    rho0dau2index[nGoodRho0] = j;
	    nGoodRho0 ++;   
	  }

	  for (int k = 0; k < nGoodPi0; ++k) {
	    TLorentzVector tmpomega = tmprho0;
	    TLorentzVector tmpph0;
	    mk4Vector(tmpph0, momentumGoodPi0[k], thetaGoodPi0[k], phiGoodPi0[k], massGoodPi0[k]);      
	    tmpomega += tmpph0;
	    tmpX = recoillep + tmpomega;
	    double omegamass = tmpomega.Mag();
	    ((TH1D*)gDirectory->Get("omegamassnocutall"))->Fill(omegamass);	   
	    if(fBVxbTyp==-14) ((TH1D*)gDirectory->Get("omegamassnocut"))->Fill(omegamass);    
	    tmpNeutrino = p4Upsilon - p4Breco - tmpX;
	    double tmpMM2omega;
	    double tmpMM2eta;
	    if(omegamass>.68&&omegamass<.88){
	      tmpMM2omega = tmpNeutrino.Mag2();
	      countomega ++;
	      if (TMath::Abs(tmpMM2omega)<TMath::Abs(bestmm2omega)) {
		bestmm2omega = tmpMM2omega;
		fMomega= omegamass;
		fMM2omega = bestmm2omega;
		bestomega = tmpomega;
		best1piome = tmpph1;
		best2piome = tmpph2;	      
		bestpi0ome = tmpph0;
		TLorentzVector dalipi1pi2 =  best1piome + best2piome;
		TLorentzVector dalipi1pi0 =  best1piome + bestpi0ome;
		fdalitzpi1pi2ome = pow(dalipi1pi2 .Mag(),2);  
		fdalitzpi1pi0ome = pow(dalipi1pi0 .Mag(),2);  
     	
		// cos theta hel    
		TVector3  omeBoost = tmpomega.Vect();
		omeBoost.SetMag(omeBoost.Mag()/tmpomega.E());
		dalipi1pi2.Boost(-omeBoost);
		TLorentzVector pi1omeboost = best1piome;
		pi1omeboost.Boost(-omeBoost);
		TLorentzVector pi2omeboost = best2piome;
		pi2omeboost.Boost(-omeBoost);
		TLorentzVector pi0omeboost = bestpi0ome;
		pi0omeboost.Boost(-omeBoost);
		TVector3 verpi1pi2 = dalipi1pi2.Vect();
		TVector3 verpi1 = pi1omeboost.Vect();
		TVector3 verpi2 = pi2omeboost.Vect();
		TVector3 verpi0 = pi0omeboost.Vect();    
		TVector3 orthog = verpi1.Cross(verpi2);
     	
		TVector3 dalipi1pi2Boost  = dalipi1pi2.Vect();
		dalipi1pi2Boost.SetMag(dalipi1pi2Boost.Mag()/dalipi1pi2.E());
		TLorentzVector pi1daliboost = pi1omeboost;
		pi1daliboost.Boost(-dalipi1pi2Boost);
		TLorentzVector pi2daliboost = pi2omeboost;
		pi2daliboost.Boost(-dalipi1pi2Boost);
		TVector3 verpi1dali = pi1daliboost.Vect();
		TVector3 verpi2dali = pi2daliboost.Vect();
     	
		verpi1pi2.SetMag(1); 
		verpi0.SetMag(1);
		verpi1dali.SetMag(1);
		verpi2dali.SetMag(1);
		omeBoost.SetMag(1);
		orthog.SetMag(1);
     	
		double costh =  verpi1dali*verpi0;
		fcosthome = TMath::Abs(costh);
     	
		best1piome.Boost(-cmsBoost);
		fmom1piome = best1piome.Vect().Mag();
		best2piome.Boost(-cmsBoost);
		fmom2piome = best2piome.Vect().Mag();
		bestpi0ome.Boost(-cmsBoost);
		fmompi0ome = bestpi0ome.Vect().Mag();
		if(fmom2piome>fmom1piome) {
		  double tmpmom = fmom1piome;
		  fmom1piome = fmom2piome;
		  fmom2piome = tmpmom;
		}
	      }
	    }
	    if(omegamass>.45&&omegamass<.65 && nGoodEta<80){
	      countetapppi0 ++;
	      tmpMM2eta = tmpNeutrino.Mag2();
	      if (TMath::Abs(tmpMM2eta)<TMath::Abs(bestmm2etapppi0)) {
		bestmm2etapppi0 = tmpMM2eta;
		fMetapppi0 = omegamass;
		fMM2etapppi0 = bestmm2etapppi0;
		bestetapppi0 = tmpomega;
	      }
	      momentumGoodEta[nGoodEta] = tmpomega.Rho();
	      thetaGoodEta[nGoodEta] = tmpomega.Theta();
	      phiGoodEta[nGoodEta] = tmpomega.Phi();
	      massGoodEta[nGoodEta] = tmpomega.Mag();
	      etadau1index[nGoodEta] = i;
	      etadau2index[nGoodEta] = j;	  
	      etadau3index[nGoodEta] = k;	  
	      if(chargeGoodPi[i]>0){
		etadau1lund[nGoodEta] = 211;
		etadau2lund[nGoodEta] = -211;	  
	      }else{
		etadau1lund[nGoodEta] = -211;
		etadau2lund[nGoodEta] = 211;	  
	      }		
	      etadau3lund[nGoodEta] = 111;
	      nGoodEta ++;   
	    }
	  }
	}
      }
    }
  }

  fNrecorho0 = countrho0;
  fmomrho0 = bestrho0.Rho();
  fthrho0 =  bestrho0.Theta();
  fphirho0 = bestrho0.Phi();
  fNrecoomega = countomega;
  fmomomega = bestomega.Rho();
  fthomega =  bestomega.Theta();
  fphiomega = bestomega.Phi();

  // rho+- calculation

  int countrho(0);
  TLorentzVector bestrho;
  bestrho.SetXYZM(0., 0., 0., 0.);     
  TLorentzVector bestpirho;
  bestpirho.SetXYZM(0., 0., 0., 0.);     
  TLorentzVector bestpi0rho;
  bestpi0rho.SetXYZM(0., 0., 0., 0.);     
  double bestmm2rho = 1000.;
  for (i = 0; i < nGoodPi; ++i) {
    if(chargeGoodPi[i]==chargelep) continue;
    for (int k = 0; k < nGoodPi0; ++k) {
      TLorentzVector tmpph;
      TLorentzVector tmpph0;
      TLorentzVector tmprho;
      TLorentzVector tmpX = recoillep;
      tmprho.SetXYZM(0., 0., 0., 0.);     
      mk4Vector(tmpph, momentumGoodPi[i], thetaGoodPi[i], phiGoodPi[i], .13957);      
      tmprho += tmpph;
      mk4Vector(tmpph0, momentumGoodPi0[k], thetaGoodPi0[k], phiGoodPi0[k], massGoodPi0[k]);      
      tmprho += tmpph0;
      tmpX += tmprho;
      double rhomass = tmprho.Mag();
      ((TH1D*)gDirectory->Get("rhomassnocutall"))->Fill(rhomass);	   
      if(fBVxbTyp==13) ((TH1D*)gDirectory->Get("rhomassnocut"))->Fill(rhomass);    
      TLorentzVector tmpNeutrino = p4Upsilon - p4Breco - tmpX;
      double tmpMM2   = tmpNeutrino.Mag2();
      countrho ++;
      if (rhomass<1.5&&TMath::Abs(tmpMM2)<TMath::Abs(bestmm2rho)) {
	bestmm2rho = tmpMM2;
	fMrho= rhomass;
	fMM2rho = bestmm2rho;
	bestrho = tmprho;
	bestpirho = tmpph;
	bestpi0rho = tmpph0;	      
	bestpirho.Boost(-cmsBoost);
	fmompirho = bestpirho.Vect().Mag();
	bestpi0rho.Boost(-cmsBoost);
	fmompi0rho = bestpi0rho.Vect().Mag();
      }
    }
  }

  fNrecorho = countrho;
  fmomrho = bestrho.Rho();
  fthrho =  bestrho.Theta();
  fphirho = bestrho.Phi();
  flMommin = 1000.;

  // eta -> pi0pi0pi0 calculation 

  int countetapi0pi0pi0(0);
  double bestmm2etapi0pi0pi0 = 1000.;
  double bestetapi0pi0pi0mass = 1000.;
  TLorentzVector bestetapi0pi0pi0;
  bestetapi0pi0pi0.SetXYZM(0., 0., 0., 0.);     
  for (i = 0; i < nGoodPi0; ++i) {
    for (int k = 0; k < nGoodPi0; ++k) {
      if(k <= i) continue;
      if((pi0dau1index[i]==pi0dau1index[k])||
	 (pi0dau1index[i]==pi0dau2index[k])||
	 (pi0dau2index[i]==pi0dau1index[k])||
	 (pi0dau2index[i]==pi0dau2index[k])) continue;            
      for (int j = 0; j < nGoodPi0; ++j) {
	if(j <= k) continue;
	if((pi0dau1index[i]==pi0dau1index[j])||
	   (pi0dau1index[i]==pi0dau2index[j])||
	   (pi0dau2index[i]==pi0dau1index[j])||
	   (pi0dau2index[i]==pi0dau2index[j])) continue;            
	if((pi0dau1index[k]==pi0dau1index[j])||
	   (pi0dau1index[k]==pi0dau2index[j])||
	   (pi0dau2index[k]==pi0dau1index[j])||
	   (pi0dau2index[k]==pi0dau2index[j])) continue;            
	TLorentzVector tmpph1;
	TLorentzVector tmpph2;
	TLorentzVector tmpph3;
	TLorentzVector tmpetapi0pi0pi0;
	TLorentzVector tmpX = recoillep;
	tmpetapi0pi0pi0.SetXYZM(0., 0., 0., 0.);     
	mk4Vector(tmpph1, momentumGoodPi0[i], thetaGoodPi0[i], phiGoodPi0[i], massGoodPi0[i]);      
	tmpetapi0pi0pi0 += tmpph1;
	mk4Vector(tmpph2, momentumGoodPi0[j], thetaGoodPi0[j], phiGoodPi0[j], massGoodPi0[j]);              
	tmpetapi0pi0pi0 += tmpph2;
	mk4Vector(tmpph3, momentumGoodPi0[k], thetaGoodPi0[k], phiGoodPi0[k], massGoodPi0[k]);              
	tmpetapi0pi0pi0 += tmpph3;
	double etapi0pi0pi0mass = tmpetapi0pi0pi0.Mag();
	if(etapi0pi0pi0mass>.45&&etapi0pi0pi0mass<.65 && nGoodEta<80) {
	  countetapi0pi0pi0 ++;
	  tmpX += tmpetapi0pi0pi0;
	  TLorentzVector tmpNeutrino = p4Upsilon - p4Breco - tmpX;
	  double tmpMM2etapi0pi0pi0   = tmpNeutrino.Mag2();
	  if (TMath::Abs(tmpMM2etapi0pi0pi0)<TMath::Abs(bestmm2etapi0pi0pi0)) {
	    bestmm2etapi0pi0pi0 = tmpMM2etapi0pi0pi0;
	    bestetapi0pi0pi0mass = etapi0pi0pi0mass;
	    bestetapi0pi0pi0 = tmpetapi0pi0pi0;
	  }
	  momentumGoodEta[nGoodEta] = tmpetapi0pi0pi0.Rho();
	  thetaGoodEta[nGoodEta] = tmpetapi0pi0pi0.Theta();
	  phiGoodEta[nGoodEta] = tmpetapi0pi0pi0.Phi();
	  massGoodEta[nGoodEta] = tmpetapi0pi0pi0.Mag();
	  etadau1index[nGoodEta] = i;
	  etadau2index[nGoodEta] = k;	  
	  etadau3index[nGoodEta] = j;	  
	  etadau1lund[nGoodEta] = 111;
	  etadau2lund[nGoodEta] = 111;	  
	  etadau3lund[nGoodEta] = 111;
	  nGoodEta ++;   
	}
      }
    }
  }

  if (countetapi0pi0pi0){
     fMM2etapi0pi0pi0 = bestmm2etapi0pi0pi0;
     fMetapi0pi0pi0 = bestetapi0pi0pi0mass;
  }
 
  //  best eta selection
  fMeta = fMetagg;
  fMM2eta = fMM2etagg;
  fModeeta = 1;
  fmometa = bestetagg.Rho();
  ftheta =  bestetagg.Theta();
  fphieta = bestetagg.Phi();
  if(TMath::Abs(fMM2etapppi0)<TMath::Abs(fMM2eta)){
    fMeta = fMetapppi0;
    fMM2eta = fMM2etapppi0;
    fModeeta = 2;
    fmometa = bestetapppi0.Rho();
    ftheta =  bestetapppi0.Theta();
    fphieta = bestetapppi0.Phi();
  }    
  if(TMath::Abs(fMM2etapi0pi0pi0)<TMath::Abs(fMM2eta)){
    fMeta = fMetapi0pi0pi0;
    fMM2eta = fMM2etapi0pi0pi0;
    fModeeta = 3;
    fmometa = bestetapi0pi0pi0.Rho();
    ftheta =  bestetapi0pi0pi0.Theta();
    fphieta = bestetapi0pi0pi0.Phi();
  }    

  if (nGoodEta > 0 && fModeeta == 1) {
    
    int idxg1tmp = etadau1index[idxBestEta];
    int idxg2tmp = etadau2index[idxBestEta];
    fEtaflag = 0;    
    if (gammafrompi0[idxg1tmp] == 1 || gammafrompi0[idxg2tmp] == 1) {
      fEtaflag = 1;
    }
    
  }

  /*
  for (int k = 0; k < nGoodEta; ++k) {
    massGoodEta[k] = 0.54775;
  }
  */

  // eta' -> rho0g calculation

  int countetaprho0g(0);
  double bestmm2etaprho0g = 1000.;
  double bestchisqetaprhog = 1000.;
  TLorentzVector bestetaprho0g;
  bestetaprho0g.SetXYZM(0., 0., 0., 0.);     
  for (i = 0; i < nGoodRho0; ++i) {
    for (int k = 0; k < nGoodGam; ++k) {
      TLorentzVector tmprho0;
      TLorentzVector tmpph;
      TLorentzVector tmpetap;
      TLorentzVector tmpX = recoillep;
      tmpetap.SetXYZM(0., 0., 0., 0.);
      TLorentzVector tmprho0ph;
      tmprho0ph.SetXYZM(0., 0., 0., 0.);
     
      mk4Vector(tmprho0, momentumGoodRho0[i], thetaGoodRho0[i], phiGoodRho0[i], massGoodRho0[i]);      
      tmpetap += tmprho0;
      mk4Vector(tmpph, momentumGoodGam[k], thetaGoodGam[k], phiGoodGam[k], 0.);      
      tmpetap += tmpph;
      tmpX += tmpetap;
      double etapmass = tmpetap.Mag();
      TLorentzVector tmpNeutrino = p4Upsilon - p4Breco - tmpX;
      double tmpMM2   = tmpNeutrino.Mag2();
      countetaprho0g ++;
      double chisqetaprhog = pow(tmpMM2/.3,2) + pow((massGoodRho0[i]-.775)/0.150,2);
      //      if (TMath::Abs(tmpMM2)<TMath::Abs(bestmm2etaprho0g) && etapmass>.86 && etapmass<1.06) {
      if (chisqetaprhog<bestchisqetaprhog&&etapmass>.86 && etapmass<1.06) {
	bestchisqetaprhog = chisqetaprhog;
	bestmm2etaprho0g = tmpMM2;
	fMetaprho0g = etapmass;
	fMM2etaprho0g = bestmm2etaprho0g;
	bestetaprho0g = tmpetap;
        mk4Vector(tmprho0ph, momentumGoodGam[k], thetaGoodGam[k], phiGoodGam[k], 0.);
	tmprho0ph.Boost(-cmsBoost);
	fMomrho0ph = tmprho0ph.Vect().Mag();
	fEtapflag = 0;         	
	if (gammafrompi0[k] == 1) {
	  fEtapflag = 1;
	}

      }
    }
  }


  // eta' -> etapipi calculation

  int countetapetapp(0);
  double bestmm2etapetapp = 1000.;
  double bestchisqetapetapp = 1000.;
  TLorentzVector bestetapetapp;
  bestetapetapp.SetXYZM(0., 0., 0., 0.);     
  int modeetapetapp;
  for (i = 0; i < nGoodRho0; ++i) {
    for (int k = 0; k < nGoodEta; ++k) {
      if((TMath::Abs(etadau1lund[k])==211 && etadau1index[k]==rho0dau1index[i])||
	 (TMath::Abs(etadau1lund[k])==211 && etadau1index[k]==rho0dau2index[i])||
	 (TMath::Abs(etadau2lund[k])==211 && etadau2index[k]==rho0dau1index[i])||
	 (TMath::Abs(etadau2lund[k])==211 && etadau2index[k]==rho0dau2index[i])) continue;            
      TLorentzVector tmprho0;
      TLorentzVector tmpeta;
      TLorentzVector tmpetap;
      TLorentzVector tmpX = recoillep;
      tmpetap.SetXYZM(0., 0., 0., 0.);     
      mk4Vector(tmprho0, momentumGoodRho0[i], thetaGoodRho0[i], phiGoodRho0[i], massGoodRho0[i]);      
      tmpetap += tmprho0;
      mk4Vector(tmpeta, momentumGoodEta[k], thetaGoodEta[k], phiGoodEta[k], 0.54775);      
      tmpetap += tmpeta;
      tmpX += tmpetap;
      double etapmass = tmpetap.Mag();
      TLorentzVector tmpNeutrino = p4Upsilon - p4Breco - tmpX;
      double tmpMM2   = tmpNeutrino.Mag2();
      countetapetapp ++;

      double chisqetapetapp = pow(tmpMM2/.3,2) + pow((massGoodEta[i]-.547)/0.03,2);

      //      if (TMath::Abs(tmpMM2)<TMath::Abs(bestmm2etapetapp) && etapmass>.86 && etapmass<1.06) {
      if (chisqetapetapp<bestchisqetapetapp&&etapmass>.86 && etapmass<1.06) {
	if(TMath::Abs(etadau1lund[k])==22) modeetapetapp = 2;
	if(TMath::Abs(etadau1lund[k])==211) modeetapetapp = 3;
	if(TMath::Abs(etadau1lund[k])==111) modeetapetapp = 4;
	bestchisqetapetapp = chisqetapetapp;
	bestmm2etapetapp = tmpMM2;
	fMetapetapp = etapmass;
	fMM2etapetapp = bestmm2etapetapp;
	fEtapmassetadau = massGoodEta[k];
	bestetapetapp = tmpetap;
      }
    }
  }

  //  best etap selection
  fMetap = fMetaprho0g;
  fMM2etap = fMM2etaprho0g;
  fModeetap = 1;
  fmometap = bestetaprho0g.Rho();
  fthetap =  bestetaprho0g.Theta();
  fphietap = bestetaprho0g.Phi();
  //  if(TMath::Abs(fMM2etapetapp)<TMath::Abs(fMM2eta)){
  if(TMath::Abs(bestchisqetapetapp)<TMath::Abs(bestchisqetaprhog)){
    fMetap = fMetapetapp;
    fMM2etap = fMM2etapetapp;
    fModeetap = modeetapetapp;
    fEtapflag = 0;             
    fmometap = bestetapetapp.Rho();
    fthetap =  bestetapetapp.Theta();
    fphietap = bestetapetapp.Phi();
  }    

  // a_0 calculation

  int counta0(0);
  TLorentzVector besta0;
  besta0.SetXYZM(0., 0., 0., 0.);     
  TLorentzVector bestetaa0;
  bestetaa0.SetXYZM(0., 0., 0., 0.);     
  TLorentzVector bestpi0a0;
  bestpi0a0.SetXYZM(0., 0., 0., 0.);     
  double bestmm2a0 = 1000.;
  double bestchisqa0 = 1000.;
  for (i = 0; i < nGoodEta; ++i) {
    for (int k = 0; k < nGoodPi0; ++k) {
      if((TMath::Abs(etadau1lund[i])==22 && etadau1index[i]==pi0dau1index[k])||
	 (TMath::Abs(etadau1lund[i])==22 && etadau1index[i]==pi0dau2index[k])||
	 (TMath::Abs(etadau2lund[i])==22 && etadau2index[i]==pi0dau1index[k])||
	 (TMath::Abs(etadau2lund[i])==22 && etadau2index[i]==pi0dau2index[k])) continue; 
      if((TMath::Abs(etadau1lund[i])==111 && pi0dau1index[etadau1index[i]]==pi0dau1index[k])||
	 (TMath::Abs(etadau1lund[i])==111 && pi0dau1index[etadau1index[i]]==pi0dau2index[k])||
	 (TMath::Abs(etadau1lund[i])==111 && pi0dau2index[etadau1index[i]]==pi0dau1index[k])||
	 (TMath::Abs(etadau1lund[i])==111 && pi0dau2index[etadau1index[i]]==pi0dau2index[k])) continue;      
      if((TMath::Abs(etadau2lund[i])==111 && pi0dau1index[etadau2index[i]]==pi0dau1index[k])||
	 (TMath::Abs(etadau2lund[i])==111 && pi0dau1index[etadau2index[i]]==pi0dau2index[k])||
	 (TMath::Abs(etadau2lund[i])==111 && pi0dau2index[etadau2index[i]]==pi0dau1index[k])||
	 (TMath::Abs(etadau2lund[i])==111 && pi0dau2index[etadau2index[i]]==pi0dau2index[k])) continue;      
      if((TMath::Abs(etadau3lund[i])==111 && pi0dau1index[etadau3index[i]]==pi0dau1index[k])||
	 (TMath::Abs(etadau3lund[i])==111 && pi0dau1index[etadau3index[i]]==pi0dau2index[k])||
	 (TMath::Abs(etadau3lund[i])==111 && pi0dau2index[etadau3index[i]]==pi0dau1index[k])||
	 (TMath::Abs(etadau3lund[i])==111 && pi0dau2index[etadau3index[i]]==pi0dau2index[k])) continue; 
      TLorentzVector tmpeta;
      TLorentzVector tmppi0;
      TLorentzVector tmpa0;
      TLorentzVector tmpX = recoillep;
      tmpa0.SetXYZM(0., 0., 0., 0.);     
      mk4Vector(tmpeta, momentumGoodEta[i], thetaGoodEta[i], phiGoodEta[i], 0.54775);      
      tmpa0 += tmpeta;
      mk4Vector(tmppi0, momentumGoodPi0[k], thetaGoodPi0[k], phiGoodPi0[k], massGoodPi0[k]);      
      tmpa0 += tmppi0;
      tmpX += tmpa0;
      double a0mass = tmpa0.Mag();
      TLorentzVector tmpNeutrino = p4Upsilon - p4Breco - tmpX;
      double tmpMM2   = tmpNeutrino.Mag2();
      counta0 ++;
      
      double chisqa0 = pow(tmpMM2/.3,2) + pow((massGoodEta[i]-.547)/0.03,2);

      //      if (TMath::Abs(tmpMM2)<TMath::Abs(bestmm2a0)&&a0mass<1.2&&a0mass>.75) {
      if (chisqa0<bestchisqa0&&a0mass<1.2&&a0mass>.75) {
	bestchisqa0 = chisqa0;
 	bestmm2a0 = tmpMM2;
	fMa0= a0mass;
	fMM2a0 = bestmm2a0;
 	besta0 = tmpa0;
	if(TMath::Abs(etadau1lund[i])==22) fModea0 = 1;
	if(TMath::Abs(etadau1lund[i])==211) fModea0 = 2;
	if(TMath::Abs(etadau1lund[i])==111) fModea0 = 3;
	fa0massetadau = massGoodEta[i];
      }
    }
  }

  fNrecoa0 = counta0;
  fmoma0 = besta0.Rho();
  ftha0 =  besta0.Theta();
  fphia0 = besta0.Phi();

  // a_0p calculation

  int counta0p(0);
  TLorentzVector besta0p;
  besta0p.SetXYZM(0., 0., 0., 0.);     
  TLorentzVector bestetaa0p;
  bestetaa0p.SetXYZM(0., 0., 0., 0.);     
  TLorentzVector bestpi0a0p;
  bestpi0a0p.SetXYZM(0., 0., 0., 0.);     
  double bestmm2a0p = 1000.;
  double bestchisqa0p = 1000.;
  for (i = 0; i < nGoodEta; ++i) {
    for (int k = 0; k < nGoodPi; ++k) {
      if(chargeGoodPi[k]==chargelep) continue;
      if((TMath::Abs(etadau1lund[i])==211 && etadau1index[i]==k)||
	 (TMath::Abs(etadau2lund[i])==211 && etadau2index[i]==k)) continue;      
      TLorentzVector tmpeta;
      TLorentzVector tmppi;
      TLorentzVector tmpa0p;
      TLorentzVector tmpX = recoillep;
      tmpa0p.SetXYZM(0., 0., 0., 0.);     
      mk4Vector(tmpeta, momentumGoodEta[i], thetaGoodEta[i], phiGoodEta[i], 0.54775);     
      tmpa0p += tmpeta;
      mk4Vector(tmppi, momentumGoodPi[k], thetaGoodPi[k], phiGoodPi[k], .13957);      
      tmpa0p += tmppi;
      tmpX += tmpa0p;
      double a0pmass = tmpa0p.Mag();
      TLorentzVector tmpNeutrino = p4Upsilon - p4Breco - tmpX;
      double tmpMM2   = tmpNeutrino.Mag2();
      counta0p ++;
      double chisqa0p = pow(tmpMM2/.3,2) + pow((massGoodEta[i]-.547)/0.03,2);

      //      if (TMath::Abs(tmpMM2)<TMath::Abs(bestmm2a0p)&&a0pmass<1.2&&a0pmass>.75) {
      if (chisqa0p<bestchisqa0p&&a0pmass<1.2&&a0pmass>.75) {
	bestchisqa0p = chisqa0p;
	bestmm2a0p = tmpMM2;
	fMa0p= a0pmass;
	fMM2a0p = bestmm2a0p;
	besta0p = tmpa0p;
	if(TMath::Abs(etadau1lund[i])==22) fModea0p = 1;
	if(TMath::Abs(etadau1lund[i])==211) fModea0p = 2;
	if(TMath::Abs(etadau1lund[i])==111) fModea0p = 3;
	fa0pmassetadau = massGoodEta[i];
      }
    }
  }

  fNrecoa0p = counta0p;
  fmoma0p = besta0p.Rho();
  ftha0p =  besta0p.Theta();
  fphia0p = besta0p.Phi();

  // ======================================================================
 

#ifndef FAST
  for (i = 0; i < nGam; ++i) {
    Brecgamtmp = B0RecGam[i];
    if(chbcand) Brecgamtmp = chBRecGam[i];
    if (!(Brecgamtmp&brecoOverlap)) continue;
    if(fNewFormat == 1&& Brecgamtmp==2)continue;// Apr 02 format
    if(lMomGam[i]<flMommin && lMomGam[i]>0) flMommin = lMomGam[i];
  }  
#endif

  totweight = totweightfRecoilTrkMult * totweightfRecoilNutMult; 

  if ((nEffGam == 0) && (fRecoilTrkMult == 1)) fEffCat = 1; 
  if ((nEffGam == 0) && (1 < fRecoilTrkMult) && (fRecoilTrkMult <= 3) ) fEffCat = 2; 
  if ((nEffGam == 0) && (3 < fRecoilTrkMult)) fEffCat = 3; 
  if ((nEffGam > 0)  && (fRecoilTrkMult == 1)) fEffCat = 4; 
  if ((nEffGam > 0)  && (1 < fRecoilTrkMult) && (fRecoilTrkMult <= 3) ) fEffCat = 5; 
  if ((nEffGam > 0)  && (3 < fRecoilTrkMult)) fEffCat = 6; 

  fPxhad   = p4Xhad.Vect().Mag(); 
  fTxhad   = p4Xhad.Theta(); 
  fFxhad   = p4Xhad.Phi(); 
  fExhad   = p4Xhad.E();
  fMxhad   = p4Xhad.Mag();

  fPxhadnu   = p4NuPart.Vect().Mag(); 
  fTxhadnu   = p4NuPart.Theta(); 
  fFxhadnu   = p4NuPart.Phi(); 
  fExhadnu   = p4NuPart.E();
  fMxhadnu   = p4NuPart.Mag();

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


  // deltaM calculation
  TLorentzVector p4XminPi(0., 0., 0., 0.);
  fwp4XminPi = 0;
  fphotp4Xmin = 0;
  fdeltaM = 11111;
  // RF 18 Jan 03 wdeltam  changed meaning: mm2 computed a' la D*lnu
  fwdeltaM = 11111;
  double pPiMin(10000.);

  for(i = 0; i < nGoodPi; ++i) {
   if(fLeptonCharge * chargeGoodPi[i] == -1 && isKaonGoodPi[i] == 0) {
     mk4Vector(p4t, momentumGoodPi[i], thetaGoodPi[i], phiGoodPi[i], PIPMASS);  
     TLorentzVector p4XminPi(p4Xhad-p4t);
     deltaMGoodPi[i] = fMxhad - p4XminPi.Mag();
     if(fdeltaM>deltaMGoodPi[i]) {
       fdeltaM = deltaMGoodPi[i];
     }   
     if(p4t.P()<pPiMin){
       // start calculation of m_neutrino^2 a' la partial reco
       // gamma_D*=E_pi/E*_pi with E*_pi=145.0 MeV from kinematics
       // the direction of D* is assumed to be the same as the soft pion
       double DstarMass=2.010;
       double EstarPi=0.145;
       if(p4t.E()>EstarPi){
	 double EDstar=p4t.E()*DstarMass/EstarPi;
	 TLorentzVector pDstar;
	 mk4Vector(pDstar, sqrt(EDstar*EDstar-DstarMass*DstarMass), thetaGoodPi[i], phiGoodPi[i],DstarMass ); 
	 // build the partial reconstruction neutrino missing mass
	 TLorentzVector neutPR=p4Upsilon - p4Breco-p4LeptonLab-pDstar;
	 fwdeltaM=neutPR.Mag2();
	 pPiMin=p4t.P();
      } else {
	 fwdeltaM=11111;
       }
     }
   }    
  }
 


  p4Brecoil   = p4Upsilon - p4Breco;
  p4BRecoilNC = p4Upsilon - p4BrecoNC;
  TLorentzVector p4Neutrino = p4Upsilon - p4Breco - p4Recoil;
  TLorentzVector pfourneu[15]; 
  //  double tmpfMM2[11];

  fMCharpart   = p4ChargPart.Mag();
  fTCharpart   = p4ChargPart.Theta(); 
  fPCharpart   = p4ChargPart.Phi(); 
  fECharpart   = p4ChargPart.E();

#ifndef FAST
  if (fOptGammas) {
    for(int jk=0;jk<15;jk++) {
      pfourneu[jk] = p4BRecoilNC - pGS[jk] - p4LeptonLab;
      tmpfMM2[jk]= pfourneu[jk].Mag2();
      tmpfNuT[jk]= pfourneu[jk].Theta();
      p4NeuPart[jk] = pGS[jk] - p4ChargPart;
      fMNeupart[jk] = p4NeuPart[jk].Mag();
      fTNeupart[jk] = p4NeuPart[jk].Theta();
      fPNeupart[jk] = p4NeuPart[jk].Phi();
      fENeupart[jk] = p4NeuPart[jk].E();
    }
  }
#endif

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

  // -- Calculate cuts
  // -----------------
#ifndef FAST
  if (fOptGammas) {
    for(int kw =0; kw < 10 ; kw++) {
      fGoodNoHole[kw] =   kFALSE;
    }
  }
#endif

  if ((fPlab > 0.5) && (fPcms > 0.5) && (TLABLO < fTlab*DR) && (fTlab*DR < TLABHI)) fGoodAccLepton = kTRUE;
  if (fGoodAccLepton) {
    if (fElectron) fGoodAccElectron = kTRUE; 
    if (fMuon)     fGoodAccMuon     = kTRUE; 
  }
  if (fGoodAccLepton && (PCMSLO < fPcms) && (PLABLO < fPlab))  fGoodLepton = kTRUE;
  if (fGoodLepton) {
    if (fElectron) fGoodElectron = kTRUE; 
    if (fMuon)     fGoodMuon     = kTRUE; 
  }
  if (fNLepton == 1) fOneLepton = kTRUE; 

  if ((fwdeltaM < PRMM2) && (fBrecoCharge == 0)) fGoodPRMM2 = kTRUE;
  if (TMath::Abs(fBrecoCharge) > 0)              fGoodPRMM2 = kTRUE;

  if ((MM2LO < fMM2) && (fMM2 < MM2HI)) fGoodMM2 = kTRUE;
  if (TMath::Abs(fRecoilCharge + fBrecoCharge) < REQTOTALCHARGE) fGoodChargeCons = kTRUE;
  if (fLeptonCharge == -fBrecoFlavor) fGoodChargeCorr = kTRUE;
  if (fGoodLepton && fOneLepton && fGoodPRMM2 && fGoodMM2 && fGoodChargeCons && fGoodChargeCorr) fGoodEvent = kTRUE;

  // -- Charged B: right- and wrong-sign spectra   RS and WS
  // -- Neutral B: right- and wrong-flavor spectra RF and WF
  if (0 != fBrecoCharge) {
    if (fBrecoCharge == -fLeptonCharge) {
      fGoodRS = kTRUE;
      if (fElectron) fGoodERS = kTRUE; 
      if (fMuon)     fGoodMRS = kTRUE; 
    } else {
      fGoodWS = kTRUE;
      if (fElectron) fGoodEWS = kTRUE; 
      if (fMuon)     fGoodMWS = kTRUE; 
    } 
  } else {
    if (fBrecoFlavor == -fLeptonCharge) {
      fGoodRF = kTRUE;
      if (fElectron) fGoodERF = kTRUE; 
      if (fMuon)     fGoodMRF = kTRUE; 
    } else {
      fGoodWF = kTRUE;
      if (fElectron) fGoodEWF = kTRUE; 
      if (fMuon)     fGoodMWF = kTRUE; 
    }
  }

  // -- "All other" cuts
  if (fGoodAccLepton                        && fGoodChargeCorr && fGoodMM2 && fGoodPRMM2 && fGoodChargeCons) faoLepton     = kTRUE;
  if (fGoodAccElectron && !fGoodAccMuon     && fGoodChargeCorr && fGoodMM2 && fGoodPRMM2 && fGoodChargeCons) faoElectron   = kTRUE;
  if (fGoodAccMuon     && !fGoodAccElectron && fGoodChargeCorr && fGoodMM2 && fGoodPRMM2 && fGoodChargeCons) faoMuon       = kTRUE;
  if (fGoodLepton      && fOneLepton        && fGoodChargeCorr             && fGoodPRMM2 && fGoodChargeCons) faoMM2        = kTRUE;
  if (fGoodLepton      && fOneLepton        && fGoodChargeCorr && fGoodMM2               && fGoodChargeCons) faoPRMM2      = kTRUE;
  if (fGoodLepton      && fOneLepton        && fGoodChargeCorr && fGoodMM2 && fGoodPRMM2                   ) faoChargeCons = kTRUE;


  if (faoLepton) {
    if (fGoodRS) faoRS = kTRUE; 
    if (fGoodWS) faoWS = kTRUE; 
    if (fGoodRF) faoRF = kTRUE; 
    if (fGoodWF) faoWF = kTRUE; 
  } 
  if (faoElectron) {
    if (fGoodERS) faoERS = kTRUE; 
    if (fGoodEWS) faoEWS = kTRUE; 
    if (fGoodERF) faoERF = kTRUE; 
    if (fGoodEWF) faoEWF = kTRUE; 
  }
  if (faoMuon) {
    if (fGoodMRS) faoMRS = kTRUE; 
    if (fGoodMWS) faoMWS = kTRUE; 
    if (fGoodMRF) faoMRF = kTRUE; 
    if (fGoodMWF) faoMWF = kTRUE; 
  }

  if (faoPRMM2) {
    if (fGoodAccElectron) faoPRMM2E = kTRUE; 
    if (fGoodAccMuon)     faoPRMM2M = kTRUE; 

    if (fGoodRS) faoPRMM2RS = kTRUE; 
    if (fGoodWS) faoPRMM2WS = kTRUE; 
    if (fGoodRF) faoPRMM2RF = kTRUE; 
    if (fGoodWF) faoPRMM2WF = kTRUE; 

    if (fGoodERS) faoPRMM2ERS = kTRUE; 
    if (fGoodEWS) faoPRMM2EWS = kTRUE; 
    if (fGoodERF) faoPRMM2ERF = kTRUE; 
    if (fGoodEWF) faoPRMM2EWF = kTRUE; 

    if (fGoodMRS) faoPRMM2MRS = kTRUE; 
    if (fGoodMWS) faoPRMM2MWS = kTRUE; 
    if (fGoodMRF) faoPRMM2MRF = kTRUE; 
    if (fGoodMWF) faoPRMM2MWF = kTRUE; 
  }

  if (faoMM2) {
    if (fGoodAccElectron) faoMM2E = kTRUE; 
    if (fGoodAccMuon)     faoMM2M = kTRUE; 

    if (fGoodRS) faoMM2RS = kTRUE; 
    if (fGoodWS) faoMM2WS = kTRUE; 
    if (fGoodRF) faoMM2RF = kTRUE; 
    if (fGoodWF) faoMM2WF = kTRUE; 

    if (fGoodERS) faoMM2ERS = kTRUE; 
    if (fGoodEWS) faoMM2EWS = kTRUE; 
    if (fGoodERF) faoMM2ERF = kTRUE; 
    if (fGoodEWF) faoMM2EWF = kTRUE; 

    if (fGoodMRS) faoMM2MRS = kTRUE; 
    if (fGoodMWS) faoMM2MWS = kTRUE; 
    if (fGoodMRF) faoMM2MRF = kTRUE; 
    if (fGoodMWF) faoMM2MWF = kTRUE; 
  }

  if (faoChargeCons) {
    if (fGoodRS) faoChargeConsRS = kTRUE; 
    if (fGoodWS) faoChargeConsWS = kTRUE; 
    if (fGoodRF) faoChargeConsRF = kTRUE; 
    if (fGoodWF) faoChargeConsWF = kTRUE; 

    if (fGoodERS) faoChargeConsERS = kTRUE; 
    if (fGoodEWS) faoChargeConsEWS = kTRUE; 
    if (fGoodERF) faoChargeConsERF = kTRUE; 
    if (fGoodEWF) faoChargeConsEWF = kTRUE; 

    if (fGoodMRS) faoChargeConsMRS = kTRUE; 
    if (fGoodMWS) faoChargeConsMWS = kTRUE; 
    if (fGoodMRF) faoChargeConsMRF = kTRUE; 
    if (fGoodMWF) faoChargeConsMWF = kTRUE; 
  }


#ifndef FAST
  if (fOptGammas) {
    for(int ky =0; ky <= 10 ; ky++) {
      if(TMath::Cos(fTNu) < (0.5 + (0.1* ky)/2)) fGoodNoHole[ky] =  kTRUE;
    }
    if ((TLABLO < fTNu*DR) && (fTNu*DR < TLABHI)) fGoodAccNu = kTRUE;  
  }

  if (fOptGammas) {
    if(fNlep == 0 && fGoodAccNu && fGoodNoHole[8]) {
      ((TH1D*)gDirectory->Get("PREcoil"))->Fill(fPrecoil);
      ((TH1D*)gDirectory->Get("TREcoil"))->Fill(fTrecoil);
      ((TH1D*)gDirectory->Get("FREcoil"))->Fill(fFrecoil);
      ((TH1D*)gDirectory->Get("EREcoil"))->Fill(fErecoil);
      ((TH1D*)gDirectory->Get("MREcoil"))->Fill(fMrecoil);
    }
    if(fNlep == 0 && fGoodAccNu && fGoodNoHole[8]) {
      ((TH1D*)gDirectory->Get("PAllev"))->Fill(fPAllev);
      ((TH1D*)gDirectory->Get("TAllev"))->Fill(fTAllev);
      ((TH1D*)gDirectory->Get("FAllev"))->Fill(fFAllev);
      ((TH1D*)gDirectory->Get("EAllev"))->Fill(fEAllev);
      ((TH1D*)gDirectory->Get("MAllev"))->Fill(fMAllev);
    }
    if(fNlep == 0 && signalBox && fGoodChargeCons && fGoodNoHole[8] && fGoodAccNu) {
      ((TH1D*)gDirectory->Get("PREcoilsig"))->Fill(fPrecoil);
      ((TH1D*)gDirectory->Get("TREcoilsig"))->Fill(fTrecoil);
      ((TH1D*)gDirectory->Get("FREcoilsig"))->Fill(fFrecoil);
      ((TH1D*)gDirectory->Get("EREcoilsig"))->Fill(fErecoil);
      ((TH1D*)gDirectory->Get("MREcoilsig"))->Fill(fMrecoil);
    }
    if(fNlep == 0 && signalBox && fGoodChargeCons && fGoodNoHole[8] && fGoodAccNu) {
      ((TH1D*)gDirectory->Get("PAllevsig"))->Fill(fPAllev);
      ((TH1D*)gDirectory->Get("TAllevsig"))->Fill(fTAllev);
      ((TH1D*)gDirectory->Get("FAllevsig"))->Fill(fFAllev);
      ((TH1D*)gDirectory->Get("EAllevsig"))->Fill(fEAllev);
      ((TH1D*)gDirectory->Get("MAllevsig"))->Fill(fMAllev);
    }
  }
#endif

  fBmassfit = fMxhadfit = fMM2fit = fQ2Fit = 0.;

#if 1
  if ((fPcms>0.) && (-100. < fMM2) && (fMM2 < 100.)) {
    
    int ILTYP; 
    if (fElectron) {
      ILTYP = 2; 
    } else if (fMuon) {
      ILTYP = 1; 
    } else {
      cout << " ??? fGoodLepton, but neither electron nor muon" << endl;
    }
    
    float CVAL[4] = {pxUps, pyUps, pzUps, eUps};


    float P_REC[28] = {
     p4BrecoNC.Px(), p4BrecoNC.Py(), p4BrecoNC.Pz(), p4BrecoNC.E(), 
     p4LeptonLab.Px(), p4LeptonLab.Py(), p4LeptonLab.Pz(), p4LeptonLab.E(), 
     p4Xhad.Px(), p4Xhad.Py(), p4Xhad.Pz(), p4Xhad.E(), 
     p4NeutrinoNC.Px(), p4NeutrinoNC.Py(), p4NeutrinoNC.Pz(), p4NeutrinoNC.E()};

    float P_FIT[28] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
		       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
    
    float CHI2T, PROBCHI2; 
    
    int   ISMEAR(fParametrization);
    int   IERR;
    int ISV = fVerbose; 
    if (fVerbose) {
      cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
      cout << "--> Event  " << fEvent << " ILTYP = " << ILTYP << endl;
      cout << "  breco   = " ;   printLorentz(p4BrecoNC);    cout  << endl;
      cout << "  xhad    = " ;   printLorentz(p4Xhad);       cout << endl;
      cout << "  neutrino= " ;   printLorentz(p4NeutrinoNC); cout << " mm2 = " << p4NeutrinoNC.Mag2() << endl;
      cout << "  p4lep   = " ;   printLorentz(p4LeptonLab);  cout << endl;
      cout << "  brecoil = " ;   printLorentz(p4BRecoilNC);  cout  << endl;
      cout << "  p4lepcms= " ;   printLorentz(p4LeptonCms);  cout << endl;
    }
    
    /*
      Added variables to study pW+EW
    */
    TLorentzVector dov(0., 0., 0., 0.);
    TLorentzVector lepdov = p4LeptonLab;
    TLorentzVector nudov = p4NeutrinoNC;
    lepdov.Boost(-cmsBoost);
    nudov.Boost(-cmsBoost);
    dov = lepdov + nudov;
    fEwPw = dov.E() + dov.P(); 
    /*
      Finished adding new variables
    */

    if (!fIsMakeParam) {
      int i1 = abcfit_interface_vub_(&ISMEAR,&ILTYP,CVAL,P_REC,P_FIT,&CHI2T,&PROBCHI2,&IERR, &ISV);
      i1 = 0;
    } else {
      if (fGoodLepton && fOneLepton && fGoodChargeCorr && (fMM2NC < 2.0) && !vubDepleted) {
      //      if (fGoodLepton && fOneLepton && fGoodChargeCorr && (fMM2NC < 2.0)) {
	float P_GEN[28] = {
	  p4BrecoGen.Px(), p4BrecoGen.Py(), p4BrecoGen.Pz(), p4BrecoGen.E(), 
	  p4LeptonGen.Px(), p4LeptonGen.Py(), p4LeptonGen.Pz(), p4LeptonGen.E(), 
	  p4XhadGen.Px(), p4XhadGen.Py(), p4XhadGen.Pz(), p4XhadGen.E(), 
	  p4MissGen.Px(), p4MissGen.Py(), p4MissGen.Pz(), p4MissGen.E()
	};
	
	int mpType(-1);
	if (fVcb==1) {
	  mpType=1;
	} else {
	  if (fVub==1) {
	    mpType=2;
	  } else {
	    mpType=0;
	  }
	}
	
	cout << ILTYP << " " << fMM2NC << " " << fMes << " " << fRecoilTrkMult << " " << fRecoilNutMult << " ";
	
	ofstream d4Pf(fDump4ParamFile[mpType], ios::app);
	d4Pf << ILTYP << " " << fMM2NC << " " << fMes << " " << fRecoilTrkMult << " " << fRecoilNutMult << " " << endl;
	dump4Param(mpType,P_REC);
	dump4Param(mpType,P_GEN);

	cout << " dumped " << endl; 
	return;
      }
    }

    if (fVerbose) {
      for (int t = 0; t < 28; t++) {
	cout << P_FIT[t] << "  "; 
	if (t%4 == 3) cout << endl;
      }
    }

#ifndef FAST
    //GamStudy
    if (fOptGammas) {
      float GP_FIT[28];
      float GCHI2T, GPROBCHI2; 
      int   GIERR;
      int GISV = fVerbose; 
      int   GISMEAR(0);
      
      for (int yt = 0; yt < 15; yt++) {
	float GP_REC[28] = {
	  p4BrecoNC.Px(), p4BrecoNC.Py(), p4BrecoNC.Pz(), p4BrecoNC.E(), 
	  p4LeptonLab.Px(), p4LeptonLab.Py(), p4LeptonLab.Pz(), p4LeptonLab.E(), 
	  pGS[yt].Px(), pGS[yt].Py(), pGS[yt].Pz(), pGS[yt].E(), 
	  p4NeutrinoNC.Px(), p4NeutrinoNC.Py(), p4NeutrinoNC.Pz(), p4NeutrinoNC.E()};
	
	abcfit_interface_vub_(&GISMEAR,&ILTYP,CVAL,GP_REC,GP_FIT,&GCHI2T,&GPROBCHI2,&GIERR, &GISV);
	TLorentzVector tmppGSfit(GP_FIT[8], GP_FIT[9], GP_FIT[10], GP_FIT[11]); 
	pGSfit[yt] = tmppGSfit;
      }
      for(int jl=0; jl<15; jl++) {
	tmpfMxhadfit[jl] = pGSfit[jl].Mag();
	tmpfTxhadfit[jl] = pGSfit[jl].Theta(); 
	tmpfFxhadfit[jl] = pGSfit[jl].Phi(); 
	tmpfExhadfit[jl] = pGSfit[jl].E();
      }
    }
#endif

    fBmass = p4BrecoNC.Mag();
    fProbChi2 = PROBCHI2;
    fChi2     = CHI2T; 
    TLorentzVector pbrecofit(P_FIT[0], P_FIT[1], P_FIT[2], P_FIT[3]); 
    fBmassfit = p4BrecoNC.P();
    TLorentzVector pelfit(P_FIT[4], P_FIT[5], P_FIT[6], P_FIT[7]); 
    TLorentzVector pxhadfit(P_FIT[8], P_FIT[9], P_FIT[10], P_FIT[11]); 
    fMxhadfit = pxhadfit.Mag();
    fExhadfit = pxhadfit.E();
    TLorentzVector pnufit(P_FIT[12], P_FIT[13], P_FIT[14], P_FIT[15]); 
    fMM2fit = pnufit.Mag2();
    fQ2Fit  = 2.*pnufit*p4LeptonLab;
    /*
      Added variables to study pW+EW
    */
    TLorentzVector dovfit(0., 0., 0., 0.);
    TLorentzVector lepdovfit = pelfit;
    ftLepFit = pelfit.E();

    TLorentzVector nudovfit = pnufit;
    lepdovfit.Boost(-cmsBoost); 
    nudovfit.Boost(-cmsBoost);
    dovfit = lepdovfit + nudovfit;
    fEwPwfit = dovfit.E() + dovfit.P(); 
    /*
      Finished adding new variables
    */

    if (!((fMxhadfit>0) || (fMxhadfit<0) || (fMxhadfit == 0))) {
      cout << "MEZZEGA: NAN in MXHADFIT " << IERR << endl;
      fMxhadfit = -999.;
      fPcms = -97.; // reset fPcms so that the event is not dumped into 'events'
      if (fReturnLog[10] == 0) fReturnString[10] = TString("MEZZEGA");
      fReturnLog[10]++;
      return;
    }   

    if (fGoodAccLepton 
	&& fGoodChargeCorr
	&& ((fGoodAccElectron && !fGoodAccMuon) || (fGoodAccMuon && !fGoodAccElectron))
	) {
      if (fMxhadfit < 1.6) {
	fLowMx = kTRUE; 
      } else {
	fLowMx = kFALSE; 
      }
    }

    if (IERR < 0) {
      /* R.F. too much printout
	 cout  << "returning, kFit " << IERR 
	 << " mxhadfit=" << fMxhadfit 
	 << " mm2= " << fMM2 
	 << " p= " << fPcms
	 << " nlep= " << fNLepton
	 << endl;
	 
      */
      if(!DOEXCLUSIVE) fPcms = -98.; // reset fPcms so that the event is not dumped into 'events'
      
      if (fReturnLog[11] == 0) fReturnString[11] = TString("IERR < 0");
      fReturnLog[11]++;
      return;
    }

    if (fVerbose) {
      cout << "fBmassfit= "; printLorentz(pbrecofit);    cout  << endl;
      cout << "fmxhadfit= "; printLorentz(pxhadfit);     cout  << endl;
      cout << "fMM2fit=   ";   printLorentz(pnufit);      cout << " mm2 = " << pnufit.Mag2()  << endl;
      cout << "fElfit=    ";   printLorentz(pelfit);      cout << " mm2 = " << pelfit.Mag2()  << endl;
      
      cout << "ierr = " << IERR << " chi2t = " <<CHI2T << endl;
      cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
    }
    
    ((TH1D*)fHistFile->Get("p100"))->Fill(PROBCHI2);  
    ((TH1D*)fHistFile->Get("p101"))->Fill(CHI2T);  
#endif
  }

#ifndef FAST
  // vertex quantities
  if(chbcand==0){
    fDx=VtxXXB0[indexbestB]-VtxXLepB0[indexbestB];
    fDy=VtxYXB0[indexbestB]-VtxYLepB0[indexbestB];
    fDz=VtxZXB0[indexbestB]-VtxZLepB0[indexbestB];
    fS2Dxx=VtxCovXXXB0[indexbestB]+VtxCovXXLepB0[indexbestB];
    fS2Dyy=VtxCovYYXB0[indexbestB]+VtxCovYYLepB0[indexbestB];
    fS2Dzz=VtxCovZZXB0[indexbestB]+VtxCovZZLepB0[indexbestB];
    fS2Dxy=VtxCovXYXB0[indexbestB]+VtxCovXYLepB0[indexbestB];
    fS2Dyz=VtxCovYZXB0[indexbestB]+VtxCovYZLepB0[indexbestB];
    fS2Dxz=VtxCovXZXB0[indexbestB]+VtxCovXZLepB0[indexbestB];
  } else {
    fDx=VtxXXChB[indexbestB]-VtxXLepChB[indexbestB];
    //    if(TMath::Abs(fDx)>100 && fNLepton>0)cout <<"urgle "<< fDx <<" index "<< indexbestB<<" X " <<VtxXXChB[indexbestB]<<" lep "<< VtxXLepChB[indexbestB]<<endl; 
    fDy=VtxYXChB[indexbestB]-VtxYLepChB[indexbestB];
    fDz=VtxZXChB[indexbestB]-VtxZLepChB[indexbestB];
    fS2Dxx=VtxCovXXXChB[indexbestB]+VtxCovXXLepChB[indexbestB];
    fS2Dyy=VtxCovYYXChB[indexbestB]+VtxCovYYLepChB[indexbestB];
    fS2Dzz=VtxCovZZXChB[indexbestB]+VtxCovZZLepChB[indexbestB];
    fS2Dxy=VtxCovXYXChB[indexbestB]+VtxCovXYLepChB[indexbestB];
    fS2Dyz=VtxCovYZXChB[indexbestB]+VtxCovYZLepChB[indexbestB];
    fS2Dxz=VtxCovXZXChB[indexbestB]+VtxCovXZLepChB[indexbestB];
  }
#endif

#ifndef FAST
  if (fOptGammas) {
    for(int jp=0;jp<15;jp++) {
      tmpfGoodMM2[jp] = fGoodEventPh[jp] =  fGoodEventPhNMNC[jp] = kFALSE;
      if ((MM2LO < tmpfMM2[jp]) && (tmpfMM2[jp] < MM2HI)) tmpfGoodMM2[jp] = kTRUE;
      if (fGoodLepton && tmpfGoodMM2[jp] && fGoodChargeCons && fGoodChargeCorr) fGoodEventPh[jp] = kTRUE;
      if (fGoodLepton && fGoodChargeCons && fGoodChargeCorr) fGoodEventPhNMNC[jp] = kTRUE;
    }
  }
  // - Low MX montecarlo breakdown 
  if (fOptCategories)  MxStudy(chbcand);
#endif
    
  // -- Fill histograms
  // ------------------
  char name[100];

  fMxhadRes    = fMxhad - fMxhadGen; 
  fMxhadfitRes = fMxhadfit - fMxhadGen; 
  fQ2Res       = fQ2 - fQ2Gen; 

  //Added Ciuchini Variables
  fwCiuc=2*(fExhadfit/BMASS);
  fxCiuc=2*(ftLepFit/BMASS);
  double tmCiuc;
  tmCiuc = pow((1-pow((fMxhadfit/fExhadfit),2)),0.5);
  fcsiCiuc= 2*tmCiuc/(1+tmCiuc); 
  

#ifndef FAST
  // -- Blind? 
  if (fOptBlind && ((fMxhadfit < 1.5)  || (fMxhad < 1.5))) {
    if (fReturnLog[12] == 0) fReturnString[12] = TString("blinding ....");
    fReturnLog[12]++;
    return;
  }

  if (fGoodLepton) {
    fHistFile->cd("breco");  sprintf(name, "ae%d", fBmode); ((TH1D*)gDirectory->Get(name))->Fill(fMes);
  }
#endif

  fHistFile->cd(); 
  if (faoLepton)     ((TH1D*)gDirectory->Get("mesalleventsA"))->Fill(fMes);
  if (faoMM2)        ((TH1D*)gDirectory->Get("mesalleventsB"))->Fill(fMes);
  if (faoChargeCons) ((TH1D*)gDirectory->Get("mesalleventsC"))->Fill(fMes);
  
#ifndef FAST
  for (int iB0=0; iB0<nB0; iB0++){
    if(bestB0 && iB0 == indexbestB) continue;
    if(TMath::Abs(d1B0Lund[iB0])==411){
      sprintf(name, "mesdccross");}
    else if(TMath::Abs(d1B0Lund[iB0])==413){ 
      sprintf(name, "mesdstarcross");
    }else{
      continue;
    }
    ((TH1D*)gDirectory->Get(name))->Fill(mseB0[iB0]);
  }
#endif

#ifndef FAST
  for (int iChB=0; iChB<nChB; iChB++){
    if((!bestB0) && iChB == indexbestB) continue;
    if(TMath::Abs(d1ChBLund[iChB])==421){
      sprintf(name, "mesd0cross");
    }else if(TMath::Abs(d1ChBLund[iChB])==423){
      sprintf(name, "mesdstar0cross");
    }else{
      continue;
    }
    ((TH1D*)gDirectory->Get(name))->Fill(mseChB[iChB]);
  }
#endif

#ifndef FAST
  if (fGoodEvent) {
    fHistFile->cd("breco"); sprintf(name, "ac%d", fBmode); ((TH1D*)gDirectory->Get(name))->Fill(fMes);
    fHistFile->cd();    
  }    
#endif

#ifndef FAST
  if (fGoodEvent&&!(vubDepleted)&&fMxhadfit < 1.55) {
    fHistFile->cd("breco"); sprintf(name, "acm%d", fBmode); ((TH1D*)gDirectory->Get(name))->Fill(fMes);
  }    
#endif

  //    if (1) {
  //      char line[200];
  //      cout << "################################################################################" << endl;
  //      sprintf(line, "i%d:  g=%d l=%d o=%d mm=%d q=%d c=%d dep=%d m=%3.2f p=%3.2f mxh=%4.3f mxhfit=%4.3f mes=%4.3f sgbox=%d side=%d ", 
  //  	    fEvent, fGoodEvent, fGoodLepton, fOneLepton, fGoodMM2, fGoodChargeCons, fGoodChargeCorr, vubDepleted, 
  //  	    fMM2, fPcms, fMxhad, fMxhadfit, fMes, signalBox, mesSideband);
  //      cout << line << endl;
  //      cout << "################################################################################" << endl;
  //    }
  
  fHistFile->cd();

  // -- Determine Event weight
  // -------------------------
  if (fIsMC) calculateEvtW8(); 
  //  cout << "event weight: " << fEvtW8 << " lumi weight: " << fLumiW8 << endl;



#ifdef FAST
  static int first(1); 
  if (first) {
    first = 0; 
    cout << "  ... booking a_dep" << endl;   fastBookHist("a_dep");     
    cout << "  ... booking a_enh" << endl;   fastBookHist("a_enh");     
    cout << "  ... booking e_dep" << endl;   fastBookHist("e_dep");     
    cout << "  ... booking e_enh" << endl;   fastBookHist("e_enh");     
    cout << "  ... booking m_dep" << endl;   fastBookHist("m_dep");     
    cout << "  ... booking m_enh" << endl;   fastBookHist("m_enh");     

    cout << "  ... booking ars_dep" << endl; fastBookHist("ars_dep"); 
    cout << "  ... booking aws_dep" << endl; fastBookHist("aws_dep"); 
    cout << "  ... booking arf_dep" << endl; fastBookHist("arf_dep"); 
    cout << "  ... booking awf_dep" << endl; fastBookHist("awf_dep"); 

    cout << "  ... booking ars_enh" << endl; fastBookHist("ars_enh"); 
    cout << "  ... booking aws_enh" << endl; fastBookHist("aws_enh"); 
    cout << "  ... booking arf_enh" << endl; fastBookHist("arf_enh"); 
    cout << "  ... booking awf_enh" << endl; fastBookHist("awf_enh"); 

    cout << "  ... booking ers_dep" << endl; fastBookHist("ers_dep"); 
    cout << "  ... booking ews_dep" << endl; fastBookHist("ews_dep"); 
    cout << "  ... booking erf_dep" << endl; fastBookHist("erf_dep"); 
    cout << "  ... booking ewf_dep" << endl; fastBookHist("ewf_dep"); 

    cout << "  ... booking ers_enh" << endl; fastBookHist("ers_enh"); 
    cout << "  ... booking ews_enh" << endl; fastBookHist("ews_enh"); 
    cout << "  ... booking erf_enh" << endl; fastBookHist("erf_enh"); 
    cout << "  ... booking ewf_enh" << endl; fastBookHist("ewf_enh"); 

    cout << "  ... booking mrs_dep" << endl; fastBookHist("mrs_dep"); 
    cout << "  ... booking mws_dep" << endl; fastBookHist("mws_dep"); 
    cout << "  ... booking mrf_dep" << endl; fastBookHist("mrf_dep"); 
    cout << "  ... booking mwf_dep" << endl; fastBookHist("mwf_dep"); 
    
    cout << "  ... booking mrs_enh" << endl; fastBookHist("mrs_enh"); 
    cout << "  ... booking mws_enh" << endl; fastBookHist("mws_enh"); 
    cout << "  ... booking mrf_enh" << endl; fastBookHist("mrf_enh"); 
    cout << "  ... booking mwf_enh" << endl; fastBookHist("mwf_enh"); 
    
    fHistFile->cd(); 
    TH1D *hl = new TH1D("l100", "looper counter", 300, 0., 300.); 

  }

  fHistFile->cd();
  if (fGoodAccLepton) {
    ((TH1D*)gDirectory->Get("derecoil"))->Fill(fDeltaE);
    ((TH1D*)gDirectory->Get("mesrecoil"))->Fill(fMes);

    if (fGoodLepton) {
      ((TH1D*)gDirectory->Get("deall"))->Fill(fDeltaE);
      ((TH1D*)gDirectory->Get("mesall"))->Fill(fMes);
      if (0 != fBrecoCharge) {
	((TH1D*)gDirectory->Get("mesallBch"))->Fill(fMes);
      } else {
	((TH1D*)gDirectory->Get("mesallBnu"))->Fill(fMes);
      }

      if (fSeedMode == 0) ((TH1D*)gDirectory->Get("mesallS0"))->Fill(fMes);
      else if (fSeedMode == 1) ((TH1D*)gDirectory->Get("mesallS1"))->Fill(fMes);
      else if (fSeedMode == 2) ((TH1D*)gDirectory->Get("mesallS2"))->Fill(fMes);
      else if (fSeedMode == 3) ((TH1D*)gDirectory->Get("mesallS3"))->Fill(fMes);

      if (fGoodEvent) {
	
	if (0 != fBrecoCharge) {
	  ((TH1D*)gDirectory->Get("mesallcutsBch"))->Fill(fMes);
	} else {
	  ((TH1D*)gDirectory->Get("mesallcutsBnu"))->Fill(fMes);
	}
	((TH1D*)gDirectory->Get("deallcuts"))->Fill(fDeltaE);
	((TH1D*)gDirectory->Get("mesallcuts"))->Fill(fMes);

	if (fSeedMode == 0) ((TH1D*)gDirectory->Get("mesallcutsS0"))->Fill(fMes);
	else if (fSeedMode == 1) ((TH1D*)gDirectory->Get("mesallcutsS1"))->Fill(fMes);
	else if (fSeedMode == 2) ((TH1D*)gDirectory->Get("mesallcutsS2"))->Fill(fMes);
	else if (fSeedMode == 3) ((TH1D*)gDirectory->Get("mesallcutsS3"))->Fill(fMes);
	
      }
    }

    if (vubDepleted) {
      fastFillHist("a_dep");
      if (fElectron) fastFillHist("e_dep");
      if (fMuon) fastFillHist("m_dep");
      if (fGoodRS) {
	fastFillHist("ars_dep");
	if (fElectron) fastFillHist("ers_dep");
	if (fMuon)     fastFillHist("mrs_dep");
      }
      if (fGoodWS) {
	fastFillHist("aws_dep");
	if (fElectron) fastFillHist("ews_dep");
	if (fMuon)     fastFillHist("mws_dep");
      }
      if (fGoodRF) {
	fastFillHist("arf_dep");
	if (fElectron) fastFillHist("erf_dep");
	if (fMuon)     fastFillHist("mrf_dep");
      }
      if (fGoodWF) {
	fastFillHist("awf_dep");
	if (fElectron) fastFillHist("ewf_dep");
	if (fMuon)     fastFillHist("mwf_dep");
      }
    } else {
      fastFillHist("a_enh");
      if (fElectron) fastFillHist("e_enh");
      if (fMuon) fastFillHist("m_enh");
      if (fGoodRS) {
	fastFillHist("ars_enh");
	if (fElectron) fastFillHist("ers_enh");
	if (fMuon)     fastFillHist("mrs_enh");
      }
      if (fGoodWS) {
	fastFillHist("aws_enh");
	if (fElectron) fastFillHist("ews_enh");
	if (fMuon)     fastFillHist("mws_enh");
      }
      if (fGoodRF) {
	fastFillHist("arf_enh");
	if (fElectron) fastFillHist("erf_enh");
	if (fMuon)     fastFillHist("mrf_enh");
      }
      if (fGoodWF) {
	fastFillHist("awf_enh");
	if (fElectron) fastFillHist("ewf_enh");
	if (fMuon)     fastFillHist("mwf_enh");
      }
    }
  } else {
    cout << fPcms << endl;
  }


  // -- mes histograms for efficiency tables
  if (fGoodAccLepton) {
    fillMesHist("recoil", "a"); 
    if (fMuon) fillMesHist("recoil", "m"); 
    if (fElectron) fillMesHist("recoil", "e"); 
    if (fBrecoCharge == 0) {
      fillMesHist("recoil", "bnu"); 
    } else {
      fillMesHist("recoil", "bch"); 
    }

    if (fGoodLepton) {
      fillMesHist("sgall", "a"); 
      if (fMuon) fillMesHist("sgall", "m"); 
      if (fElectron) fillMesHist("sgall", "e"); 
      if (fBrecoCharge == 0) {
	fillMesHist("sgall", "bnu");
      } else {
	fillMesHist("sgall", "bch");
      }
    }
    
    if (vubDepleted) {
      fillMesHist("sgvcb", "a"); 
      if (fMuon) fillMesHist("sgvcb", "m"); 
      if (fElectron) fillMesHist("sgvcb", "e"); 
      if (fBrecoCharge == 0) {
	fillMesHist("sgvcb", "bnu"); 
      } else {
	fillMesHist("sgvcb", "bch"); 
      }
    } else {
      fillMesHist("sgvub", "a"); 
      if (fMuon) fillMesHist("sgvub", "m"); 
      if (fElectron) fillMesHist("sgvub", "e"); 
      if (fBrecoCharge == 0) {
	fillMesHist("sgvub", "bnu"); 
      } else {
	fillMesHist("sgvub", "bch"); 
      }
    }
  }
      

//    TH1D *l100 = (TH1D*)fHistFile->Get("l100"); 
//    if (fGoodAccLepton) {
//      l100->Fill(0.);
//      if (fSeedMode == 0) l100->Fill(100.); 
//      else if (fSeedMode == 1) l100->Fill(101.); 
//      else if (fSeedMode == 2) l100->Fill(102.); 
//      else if (fSeedMode == 3) l100->Fill(103.); 
//    }

//    if (fGoodEvent) {
//      l100->Fill(1.);
//      if (fSeedMode == 0) l100->Fill(200.); 
//      else if (fSeedMode == 1) l100->Fill(201.); 
//      else if (fSeedMode == 2) l100->Fill(202.); 
//      else if (fSeedMode == 3) l100->Fill(203.); 
//    }


//    for (int itrk = 0; itrk < nTrk; ++itrk) {
//      if ((goodTrack[itrk] == 0) && (loopTrack[i] == 0)) continue;
//      if (fGoodAccLepton) {
//        if (loopTrack[i] == 0) l100->Fill(10.);
//        if (loopTrack[i] == 1) l100->Fill(11.);
//        if (loopTrack[i] == 2) l100->Fill(12.);
//        if (loopTrack[i] == 3) l100->Fill(13.);
      
//        if (fSeedMode == 0) {
//  	if (loopTrack[i] == 0) l100->Fill(110.);
//  	if (loopTrack[i] == 1) l100->Fill(111.);
//  	if (loopTrack[i] == 2) l100->Fill(112.);
//  	if (loopTrack[i] == 3) l100->Fill(113.);
//        } else if (fSeedMode == 1) {
//  	if (loopTrack[i] == 0) l100->Fill(120.);
//  	if (loopTrack[i] == 1) l100->Fill(121.);
//  	if (loopTrack[i] == 2) l100->Fill(122.);
//  	if (loopTrack[i] == 3) l100->Fill(123.);
//        } else if (fSeedMode == 2) {
//  	if (loopTrack[i] == 0) l100->Fill(130.);
//  	if (loopTrack[i] == 1) l100->Fill(131.);
//  	if (loopTrack[i] == 2) l100->Fill(132.);
//  	if (loopTrack[i] == 3) l100->Fill(133.);
//        } else if (fSeedMode == 3) {
//  	if (loopTrack[i] == 0) l100->Fill(140.);
//  	if (loopTrack[i] == 1) l100->Fill(141.);
//  	if (loopTrack[i] == 2) l100->Fill(142.);
//  	if (loopTrack[i] == 3) l100->Fill(143.);
//        }

//      }			     

//      if (fGoodEvent) {	     
//        if (loopTrack[i] == 0) l100->Fill(20.);
//        if (loopTrack[i] == 1) l100->Fill(21.);
//        if (loopTrack[i] == 2) l100->Fill(22.);
//        if (loopTrack[i] == 3) l100->Fill(23.);

//        if (fSeedMode == 0) {
//  	if (loopTrack[i] == 0) l100->Fill(210.);
//  	if (loopTrack[i] == 1) l100->Fill(211.);
//  	if (loopTrack[i] == 2) l100->Fill(212.);
//  	if (loopTrack[i] == 3) l100->Fill(213.);
//        } else if (fSeedMode == 1) {
//  	if (loopTrack[i] == 0) l100->Fill(220.);
//  	if (loopTrack[i] == 1) l100->Fill(221.);
//  	if (loopTrack[i] == 2) l100->Fill(222.);
//  	if (loopTrack[i] == 3) l100->Fill(223.);
//        } else if (fSeedMode == 2) {
//  	if (loopTrack[i] == 0) l100->Fill(230.);
//  	if (loopTrack[i] == 1) l100->Fill(231.);
//  	if (loopTrack[i] == 2) l100->Fill(232.);
//  	if (loopTrack[i] == 3) l100->Fill(233.);
//        } else if (fSeedMode == 3) {
//  	if (loopTrack[i] == 0) l100->Fill(240.);
//  	if (loopTrack[i] == 1) l100->Fill(241.);
//  	if (loopTrack[i] == 2) l100->Fill(242.);
//  	if (loopTrack[i] == 3) l100->Fill(243.);
//        }
//      }
//    }

#else

  if (fGoodAccLepton) {
    ((TH1D*)gDirectory->Get("derecoil"))->Fill(fDeltaE);
    ((TH1D*)gDirectory->Get("mesrecoil"))->Fill(fMes);
    if (faoLepton)     ((TH1D*)gDirectory->Get("mesrecoilA"))->Fill(fMes);
    if (faoMM2)        ((TH1D*)gDirectory->Get("mesrecoilB"))->Fill(fMes);
    if (faoChargeCons) ((TH1D*)gDirectory->Get("mesrecoilC"))->Fill(fMes);
    fillRecoilHist("recoil",chbcand);
    if (mesSideband) fillRecoilHist("bgrecoil",chbcand);
    if (signalBox) fillRecoilHist("sgrecoil",chbcand);
    fillMesHist("recoil", "a");

    if (fGoodLepton) {
      ((TH1D*)gDirectory->Get("deall"))->Fill(fDeltaE);
      ((TH1D*)gDirectory->Get("mesall"))->Fill(fMes);
      if (faoLepton)     ((TH1D*)gDirectory->Get("mesallA"))->Fill(fMes);
      if (faoMM2)        ((TH1D*)gDirectory->Get("mesallB"))->Fill(fMes);
      if (faoChargeCons) ((TH1D*)gDirectory->Get("mesallC"))->Fill(fMes);
      if (mesSideband) fillRecoilHist("bgall",chbcand);
      if (signalBox) fillRecoilHist("sgall",chbcand);
      fillMesHist("sgall", "a");

      
      if (fElectron) {
	((TH1D*)gDirectory->Get("deallel"))->Fill(fDeltaE);
	((TH1D*)gDirectory->Get("mesallel"))->Fill(fMes);
	if (faoLepton)     ((TH1D*)gDirectory->Get("mesallelA"))->Fill(fMes);
	if (faoMM2)        ((TH1D*)gDirectory->Get("mesallelB"))->Fill(fMes);
	if (faoChargeCons) ((TH1D*)gDirectory->Get("mesallelC"))->Fill(fMes);
      } 
      if (fMuon) {
	((TH1D*)gDirectory->Get("deallmu"))->Fill(fDeltaE);
	((TH1D*)gDirectory->Get("mesallmu"))->Fill(fMes);
	if (faoLepton)     ((TH1D*)gDirectory->Get("mesallmuA"))->Fill(fMes);
	if (faoMM2)        ((TH1D*)gDirectory->Get("mesallmuB"))->Fill(fMes);
	if (faoChargeCons) ((TH1D*)gDirectory->Get("mesallmuC"))->Fill(fMes);
      }
    }  
    if (vubDepleted) {
      ((TH1D*)gDirectory->Get("devcb"))->Fill(fDeltaE);
      ((TH1D*)gDirectory->Get("mesvcb"))->Fill(fMes);
      if (faoLepton)     ((TH1D*)gDirectory->Get("mesvcbA"))->Fill(fMes);
      if (faoMM2)        ((TH1D*)gDirectory->Get("mesvcbB"))->Fill(fMes);
      if (faoChargeCons) ((TH1D*)gDirectory->Get("mesvcbC"))->Fill(fMes);
      if (mesSideband) fillRecoilHist("bgvcb",chbcand);
      if (signalBox) fillRecoilHist("sgvcb",chbcand);

      fillMesHist("sgvcb", "a");

    } else {
      ((TH1D*)gDirectory->Get("devub"))->Fill(fDeltaE);
      ((TH1D*)gDirectory->Get("mesvub"))->Fill(fMes);
      if (faoLepton)     ((TH1D*)gDirectory->Get("mesvubA"))->Fill(fMes);
      if (faoMM2)        ((TH1D*)gDirectory->Get("mesvubB"))->Fill(fMes);
      if (faoChargeCons) ((TH1D*)gDirectory->Get("mesvubC"))->Fill(fMes);
      if (mesSideband) fillRecoilHist("bgvub",chbcand);
      if (signalBox) fillRecoilHist("sgvub",chbcand);

      fillMesHist("sgvub", "a");
    }

  }

  fHistFile->cd();
  if (fGoodEvent) {
    ((TH1D*)gDirectory->Get("deallcuts"))->Fill(fDeltaE);
    ((TH1D*)gDirectory->Get("mesallcuts"))->Fill(fMes);
    if (faoLepton)     ((TH1D*)gDirectory->Get("mesallcutsA"))->Fill(fMes);
    if (faoMM2)        ((TH1D*)gDirectory->Get("mesallcutsB"))->Fill(fMes);
    if (faoChargeCons) ((TH1D*)gDirectory->Get("mesallcutsC"))->Fill(fMes);
  }

#endif
  
  fHistFile->cd();
  
  if (fDump & 16) {
    double mks(-99.), pks(-99.), tks(-99.), fks(-99.), rks(-99.), r3ks(-99.), 
      mpi1(-99.), mpi2(-99.), epi1(-99.), epi2(-99.), mctks(-99.), mctpi1(-99.), mctpi2(-99.), 
      eg1(-99.), eg2(-99.), eg3(-99.), eg4(-99.), 
      we(-99.), wk(-99.), nr(-99.), good(-99.);
    
    static Bool_t firstKs(kTRUE);
    if (firstKs) {
      firstKs = kFALSE;
      fHistFile->cd();
      
      fKsTree = new TTree("tkshorts", "tkshorts"); 
      fKsTree->Branch("good", &good, "good/D");
      fKsTree->Branch("gevt", &fGoodEvent, "gevt/b");
      fKsTree->Branch("gccons", &fGoodChargeCons, "gccons/b");
      fKsTree->Branch("glep", &fGoodLepton, "glep/b");
      fKsTree->Branch("gmm2", &fGoodMM2, "gmm2/b");
      fKsTree->Branch("we", &we, "we/D");
      fKsTree->Branch("wk", &wk, "wk/D");
      fKsTree->Branch("nr", &nr, "nr/D");
      fKsTree->Branch("mks", &mks, "mks/D");
      fKsTree->Branch("pks", &pks, "pks/D");
      fKsTree->Branch("tks", &tks, "tks/D");
      fKsTree->Branch("fks", &fks, "fks/D");
      fKsTree->Branch("rks", &rks, "rks/D");
      fKsTree->Branch("r3ks", &r3ks, "r3ks/D");
      fKsTree->Branch("mpi1", &mpi1, "mpi1/D");
      fKsTree->Branch("mpi2", &mpi2, "mpi2/D");
      fKsTree->Branch("epi1", &epi1, "epi1/D");
      fKsTree->Branch("epi2", &epi2, "epi2/D");
      fKsTree->Branch("mcks", &mctks, "mcks/D");
      fKsTree->Branch("mcpi1", &mctpi1, "mcpi1/D");
      fKsTree->Branch("mcpi2", &mctpi2, "mcpi2/D");
      fKsTree->Branch("eg1", &eg1, "eg1/D");
      fKsTree->Branch("eg2", &eg2, "eg2/D");
      fKsTree->Branch("eg3", &eg3, "eg3/D");
      fKsTree->Branch("eg4", &eg4, "eg4/D");
    }
    
    // -- Fill Ks->pi+pi- into reduced tree
    if (fGoodAccLepton) {
      fHistFile->cd();
      for (i = 0; i < nKs; ++i) {
	if (TMath::Abs(d1KsLund[i]) != 211) continue;
	if (TMath::Abs(d2KsLund[i]) != 211) continue;
	int pi1 = d1KsIndex[i]-1;
	int pi2 = d2KsIndex[i]-1; 
	we = goodWe[i];
	wk = goodWk[i];
	nr = goodNr[i];
	mks = massKs[i];
	pks = pKs[i];
	tks = thKs[i];
	fks = phiKs[i];
	r3ks= TMath::Sqrt(xKs[i]*xKs[i] + yKs[i]*yKs[i] + zKs[i]*zKs[i]);
	r3ks= TMath::Sqrt( (xKs[i]-beamSX)*(xKs[i]-beamSX) + (yKs[i]-beamSY)*(yKs[i]-beamSY) + (zKs[i]-beamSZ)*(zKs[i]-beamSZ) );
	rks = TMath::Sqrt(xKs[i]*xKs[i] + yKs[i]*yKs[i]);
	rks = TMath::Sqrt( (xKs[i]-beamSX)*(xKs[i]-beamSX) + (yKs[i]-beamSY)*(yKs[i]-beamSY) );
	epi1 = momentumTrk[pi1];
	epi2 = momentumTrk[pi2];
	
	mpi1 = mpi2 = eg1 = eg2 = eg3 = eg4 = -99.;
	if (MCKs[i] > -1) mctks = idMc[MCKs[i]-1];
	mctpi1 = idTrk[pi1];
	mctpi2 = idTrk[pi2];
	good = -1.; if (goodKshort[i] >= 1) good = 1;
	fKsTree->Fill();
      }    
      
    }
  }
  
  fHistFile->cd();

}

#include "splitOff.icc"
#include "mcTruth.icc"

#include "makeParam.icc"
// ----------------------------------------------------------------------
Int_t recoilNtp::filterK0sPi0Pi0() {
  Int_t npi0(0), nK0s(0);    
  Int_t im(0);

  static int first(1);

  TDirectory *old = gDirectory;
  fHistFile->cd();
  TH1D *h; 
  if (first) {
    first = 0; 
    h = new TH1D("k0s0", "", 1000, 0., 1.);
    h = new TH1D("k0s1", "", 200, 0., 1.);
    h = new TH1D("k0s2", "", 200, 0., 1.);
  }

  for (im = 0; im < nMc; ++im) {
    if (idMc[im] != 310) continue;
    npi0 = 0;
    for (int id = 0; id < nMc; ++id) {
      if (mothMc[id]-1 == im) {
	if (idMc[id] == 111) npi0++;
      }
    }
    if ((nDauMc[im] == 2) && (npi0 == 2)) nK0s++;
  }

  
  if (nK0s > 0) {
    for (im = 0; im < nGam; ++im) {
      if (ecalGam[im] < 0) continue;
      ((TH1D*)gDirectory->Get("k0s0"))->Fill(ecalGam[im]);
      ((TH1D*)gDirectory->Get("k0s1"))->Fill(s1s9Gam[im]);
      if (ecalGam[im] > 1.8) ((TH1D*)gDirectory->Get("k0s2"))->Fill(s1s9Gam[im]);
    }
  }
  old->cd();
  return nK0s;
}


// ----------------------------------------------------------------------
Double_t recoilNtp::kPlus() {
  double mB(5.279), mb(4.800);

  TDirectory *old = gDirectory;
  fHistFile->cd();
  static Bool_t first(kTRUE);
  if (first) {
    first = kFALSE;
    TH1D *h;
    char name[100], title[100];
    sprintf(name, "qPlus");  sprintf(title, "qPlus");  h = new TH1D(name, title, 100, 0., 2.5); 
    sprintf(name, "kPlus");  sprintf(title, "kPlus");  h = new TH1D(name, title, 100, -2., 1.); 

    sprintf(name, "kPlus1");  sprintf(title, "kPlus1");  h = new TH1D(name, title, 100, -2., 1.); 
    sprintf(name, "kPlus2");  sprintf(title, "kPlus2");  h = new TH1D(name, title, 100, -2., 1.); 
    sprintf(name, "kPlus3");  sprintf(title, "kPlus3");  h = new TH1D(name, title, 100, -2., 1.); 
    sprintf(name, "kPlus4");  sprintf(title, "kPlus4");  h = new TH1D(name, title, 100, -2., 1.); 
    sprintf(name, "kPlus5");  sprintf(title, "kPlus5");  h = new TH1D(name, title, 100, -2., 1.); 

  }    

  char line[200];
  int xu(-99), xum(-99), xud(-99); 
  for (int i = 0; i < nMc; ++i) {
    if ((TMath::Abs(idMc[i]) == 41) || (TMath::Abs(idMc[i]) == 42)) {
      xu = i;
      xum= mothMc[i]-1;
    }
    if (mothMc[i]-1 == xu) {
      xud = i;
      break;
    }
  }

  if (xu == -99) {
    old->cd();
    return -99.;
  }

  double ds = TMath::Sqrt(TMath::Power(xMc[xu]-xMc[xum],2) + TMath::Power(yMc[xu]-yMc[xum],2) + TMath::Power(zMc[xu]-zMc[xum],2));
  double ct = 10.*ds*massMc[xu]/pMc[xu];
  double qp = 10000.*ct;
  double kp = mB - mb - qp;

  ((TH1D*)gDirectory->Get("qPlus"))->Fill(qp, 1.);
  ((TH1D*)gDirectory->Get("kPlus"))->Fill(kp, 1.);

  double w8(0.);
  w8 = fw8(kp, 4.80, 1.29); 
  ((TH1D*)gDirectory->Get("kPlus1"))->Fill(kp, w8);
  w8 = fw8(kp, 4.95, 1.29); 
  ((TH1D*)gDirectory->Get("kPlus2"))->Fill(kp, w8);
  w8 = fw8(kp, 4.65, 1.29); 
  ((TH1D*)gDirectory->Get("kPlus3"))->Fill(kp, w8);
  w8 = fw8(kp, 4.80, 0.38); 
  ((TH1D*)gDirectory->Get("kPlus4"))->Fill(kp, w8);
  w8 = fw8(kp, 4.80, 3.60); 
  ((TH1D*)gDirectory->Get("kPlus5"))->Fill(kp, w8);

  if (fVerbose) cout << "xu = " << xu << " mass = " << massMc[xu] << " qPlus = " << qp << " kPlus = " << kp
		     << " -> daughter = " << xud << " mother = " << xum << endl;
  old->cd();
  return kp;
}


// ----------------------------------------------------------------------
void recoilNtp::dumpOneB(int b1) {
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
void recoilNtp::dumpGeneratorBlock(int b1, int b2) {
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

  // flavor variables: 1 for B0, -1 for B+, 2 for B0bar, -2 for B-
int recoilNtp::getBflavor(int Bid) {
  int flavor=0;
  
  int sign=Bid/TMath::Abs(Bid);
  
  // neutral
  if(TMath::Abs(Bid)==511) {
    flavor=1+(1-sign)/2;
  }
  // charged
  if(TMath::Abs(Bid)==521) {
    flavor=-1-(1-sign)/2;
  }
    
  return flavor;
}

int recoilNtp::getBflavor(int BRECOflavor, int charge) {
  int flavor=0;
  
  // fBrecoFlavor is the sign of a prompt lepton if the BRECO decayed semileptonically
  if(BRECOflavor==1 && charge==1) flavor=-1;
  if(BRECOflavor==-1 && charge==-1) flavor=-2;
  if(BRECOflavor==1 && TMath::Abs(charge)<0.5) flavor=1;
  if(BRECOflavor==-1 && TMath::Abs(charge)<0.5) flavor=2;
  
  // cout << "(" << BRECOflavor << ", " << charge << ") -> " << flavor << endl;
  if(flavor==0) {
    cout << "ERROR: getBflavor() returned 0 for (" << BRECOflavor << ", " << charge << ")! Exiting..." << endl;
  }
    
  return flavor;
}

// get the seed id of the generated B
int recoilNtp::getSeedDbyimc(int inimc) {
  int result=-1;
  
  for (Int_t imc = 0; imc < nMc; ++imc) {
    int lid = idMc[imc];
    int lstatus=0;
    if((TMath::Abs(lid)==411)||(TMath::Abs(lid)==421)||(TMath::Abs(lid)==413)||(TMath::Abs(lid)==423)) lstatus=2;
    int mothindex=mothMc[imc] - 1;
    if(mothindex==inimc && lstatus==2) {
      // charm daughter detected
      result=TMath::Abs(lid);
      break;
    }
  }
  
  if(fVerbose) cout << "getSeedDbyimc(" << inimc << ") = " << result << endl;

  return result;
}


// get the seed id of the BRECO candidate
int recoilNtp::getSeedDbymode(int inmodeB) {
  int result=-1;

  switch(int(inmodeB/1000)) {
  case 11: {
    result=421;
    break;
  }
  case 12: {
    result=411;
    break;
  }
  case 13: {
    result=413;
    break;
  }
  case 14: {
    result=423;
    break;
  }
  case 15: {
    result=423;
    break;
  }
  default: result=-1;
  }

  if(fVerbose) cout << "getSeedDbymode(" << inmodeB << ") = " << result << endl;

  return result;
}


void recoilNtp::goodBreco(int lbrecoI) {

  if (!fIsMC) return;

  TLorentzVector p4lBGen(p4BrecoGen);
  TLorentzVector p4lBReco(p4BrecoNC);
  p4lBReco.Boost(-upsBoost);
  p4lBGen.Boost(-upsBoost);
  double pcmsBReco = p4lBReco.Vect().Mag();
  double tcmsBReco = p4lBReco.Theta();
  double fcmsBReco = p4lBReco.Phi();
  double pcmsBGen = p4lBGen.Vect().Mag();
  double tcmsBGen = p4lBGen.Theta();
  double fcmsBGen = p4lBGen.Phi();
  double distance=99999.;
  Bool_t dauBmatch(kFALSE);
  Bool_t seedBmatch(kFALSE);
  int ndauB, lmodeB, fTruBflavor;
  const double PI=3.1415935;
  int lid = idMc[lbrecoI];
  if((TMath::Abs(lid)==511)||(TMath::Abs(lid)==521)){
    fTruBflavor=getBflavor(lid);
    double deltath=tcmsBGen-tcmsBReco;
    // wrap dphi
    double deltaphi=fcmsBGen-fcmsBReco;
    while (deltaphi >  PI) deltaphi -= 2*PI;
    while (deltaphi < -PI) deltaphi += 2*PI;
    distance=sqrt(deltath*deltath+deltaphi*deltaphi);
    if (bestB0){
      ndauB = ndauB0[indexbestB];
      lmodeB = modeB0[indexbestB];
    } else {
      ndauB = ndauChB[indexbestB];
      lmodeB = modeChB[indexbestB];
    }
    dauBmatch=(nDauMc[lbrecoI]==ndauB);
    seedBmatch=(getSeedDbyimc(lbrecoI)==getSeedDbymode(lmodeB));
  }
  if(fTruBflavor == getBflavor(fBrecoFlavor, fBrecoCharge)){
    if(distance<0.4) fBrecoQual += 1;
    if(dauBmatch) fBrecoQual += 2;
    if(seedBmatch) fBrecoQual += 4;  
  }
//      for(int dcnt=0;dcnt<nMc;dcnt++)
//        if(mothMc[dcnt]-1 == lbrecoI)
//          cout << idMc[dcnt] << endl;
//      cout << "---------" << endl;
//      if(bestB0)
//        cout << d1B0Lund[indexbestB] << " " << d2B0Lund[indexbestB] << " "<< d3B0Lund[indexbestB] << " " << d4B0Lund[indexbestB] << " " << d5B0Lund[indexbestB] << " " << d6B0Lund[indexbestB] << " " << d7B0Lund[indexbestB] << endl;
//      else 
//        cout << d1ChBLund[indexbestB] << " " << d2ChBLund[indexbestB] << " " << d3ChBLund[indexbestB] << " " << d4ChBLund[indexbestB] << " " << d5ChBLund[indexbestB] << " " << d6ChBLund[indexbestB] << " " << d7ChBLund[indexbestB] << endl;
//      cout << "----------" << endl;
//      cout << lid << " " << getBflavor(fBrecoFlavor, fBrecoCharge) << " " << fTruBflavor << " " << distance << " " << nDauMc[lbrecoI] << " " << ndauB << " " << getSeedDbyimc(lbrecoI) << " " << getSeedDbymode(lmodeB) << " " << fBrecoQual << endl;
//      cout << "*************" << endl;
}


// ----------------------------------------------------------------------
void recoilNtp::Loop(Int_t maxEvent, Int_t startEvent, Int_t isVerbose, Int_t lun) {

  int step(1000);
  int ischB(0);
  findPro = findUps = 0;
  // 
  fVerbose = isVerbose; 

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

  cout << "Running over " << nentries << " events" << endl;

  Int_t nbytes = 0, nb = 0;
  int nvxbevt(0); 
  const char *pChar; 

#if FAST
  int oldrunnumber(0); 
  TString rfile(fHistFile->GetName());
  rfile.ReplaceAll(".root", ".runs"); 
  ofstream fRUN(rfile); 
#endif

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


#ifdef FAST
    if (fRunnumber != oldrunnumber) {
      fRUN << fRunnumber << endl;
      oldrunnumber = fRunnumber; 
    }
#endif

    //For the events without Breco candidate:
    if(((nB0==0)&&(nChB==0))||((modeB0[0]==-1)&&(modeChB[0]==-1)&&(nB0==1)&&(nChB==1))){
      if(fIsMC){
        //        cout << jentry << "th events, only truth info" << endl;
        mcTruth(ischB);
        //        if((fDump & 4))
        //          fTree->Fill();
        if((fDump & 32) && (fVub>0))
          fTree->Fill();
      }
      continue;
    }

    findbestB();
    
    brecoOverlap = 1; // this is flag for overlap with BRECO candidate
    if(indexbestB == 1) brecoOverlap = 2;

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
#ifndef FAST
    if (fOptMakeEventList) {
      if (jentry%2 == 0) {
	//      cout << "--> Copying " << jentry << " into fToBeCopied " << bla << endl;
	//	fToBeCopied->Enter(jentry);
      }
    }
#endif

    ischB = 0;
    if (bestB0 == 0) {
      ischB = 1;
    }

    // -- MonteCarlo Truth
    brecoI=-99;

    if (fIsMC) {
      mcTruth(ischB);
      if (fBVxb == fB1Index) {
	brecoI = fB2Index;
      } else {
	brecoI = fB1Index;
      }
      fBadReco=0;
      tmpPgen = pMc[brecoI];
      tmpThetagen = thetaMc[brecoI];
      tmpPhigen = phiMc[brecoI];
      tmpMassgen = massMc[brecoI];
      mk4Vector(p4BrecoGen, tmpPgen, tmpThetagen, tmpPhigen, tmpMassgen);
    }

    //      if (fBVxb > 0) {
    //        ++nvxbevt; 
    //        cout << "beginevent: " << nvxbevt << endl;
    //        dumpOneB(fBVxb); 
    //        cout << "endevent: " << nvxbevt << endl;
    //      }

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
    
    fdpb=0;
    if(fIsMC){
      fdpb=(p4Breco.Vect()-p4BrecoGen.Vect()).Mag();	
    }
    p4Upsilon = TLorentzVector(pxUps, pyUps, pzUps, eUps); 
    //The following 2 lines were added by me (KT) since I did not find any initialization of upsBoost anywhere...
    upsBoost = TVector3(pxUps, pyUps, pzUps);
    upsBoost.SetMag(upsBoost.Mag()/eUps);

    p4Brecoil = p4Upsilon - p4Breco; 
    cmsBoost = p4Brecoil.BoostVector();

    // track and neutrals selection must be after the computation of the recoil B, which is used
    selectTracks();
    //    trackStudy(); 

    // splitoff studies --- needed for photon selection
    doSplitOffStudy();
    selectPhotons();

    // ks bug studies
#ifndef FAST
    if(fIsMC){
      fHistFile->cd("mcTruth");
      bool done[200];
      TH1D* histOk=(TH1D*)gDirectory->Get("h77001");
      TH1D* histNo=(TH1D*)gDirectory->Get("h77002");
      TH1D* histKs=(TH1D*)gDirectory->Get("h77003");
      for(int i=0;i<nMc;i++)done[i]=false;
      for (int ig=0;ig<nGam;ig++){
	if(goodPhoton[ig]==0)continue;
	Int_t iKs=isKs2Pi0Dau(ig);
	if(iKs<0){
	  histNo->Fill(energyGam[ig]);
	  continue;
	}
	
	histOk->Fill(energyGam[ig]);
	if(!done[iKs])histKs->Fill(pMc[iKs]);
	done[iKs]=true;
      }
    }
#endif

    TLorentzVector p4t = p4Breco; 
    p4t.Boost(-upsBoost);
    fPcmsBreco = p4t.Vect().Mag(); 

    if (fVerbose) cout << "CALLING breco()" << endl;
    breco(ischB);
    if (fOptCategories) mxCategory();
    maskKshorts(1);
    //     maskPi0(1);
    //    maskConversions();
    if (fVerbose) cout << "CALLING recoil()" << endl;

    recoil(ischB);

    //For the Breco quality
    goodBreco(brecoI);

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
      if (fDump & 16 &&((fPcms > 0.) || (fVub == 1) )) {
	if (isVerbose) cout << " fill fTree 16 " << endl;
	fTree->Fill();
      }
      //fDump 32 used for the signal MC without requiring a Breco candidate, but dumping only events with a true signal decay
      if ((fDump & 32) && ((fPcms > 0.) || (fVub == 1) || (fVcb == 1))) {
	if (isVerbose) cout << " fill fTree 32" << endl;
        if(fVub>0)
          fTree->Fill();
      }

    }

    if(fOptMakeEventList == 2 && fGoodEvent){
      fToBeCopied->Enter(jentry);
      fSkimcount++;
      //      cout << fSkimcount << endl;
    }

    if(fVerbose) cout << " jentry = " << jentry << "  mes = " << fMes << " pcms = " << fPcms << " mass = " << p4LeptonLab.Mag()
		      << " xhadmass = " << fMxhad << " mm2 = " << fMM2NC << " xhadfitmass = " << fMxhadfit
		      << " vub = " << fVub << "  vcb = " << fVcb << endl;
  }

  cout << "----------------------------------------------------------------------" << endl;
  if (fOptGammas)  cout<<findPro<<" "<<findUps<<endl;
}

// ----------------------------------------------------------------------

void recoilNtp::Skim(Double_t pCut, Int_t maxEvent, Int_t startEvent, Int_t isVerbose, const char *ITSfile) {

  
  findPro = findUps = 0;
  // 
  fVerbose = isVerbose; 
  double tmpMassPB, tmpMassThetaB, tmpMassPhiB ;
  double tmpPB, tmpThetaB, tmpPhiB ;
  double tmpPgen, tmpThetagen, tmpPhigen, tmpMassgen ;
  double tmpMB, tmpBevM;
  //
  
  int step(1000);
  int ischB(0);

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
  int nk0sEvents(0); 
  Int_t nbytes = 0, nb = 0;
  Int_t Brectrktmp;
  TLorentzVector p4l; 
  Double_t mass(0.), pmax(0.), plab(0.), pcms(0.);

  for (Int_t jentry = startEvent; jentry < startEvent+maxEvent; jentry++) {
   if (fReturnLog[0] == 0) fReturnString[0] = TString("Loop event counter");
    fReturnLog[0]++;
    if (isVerbose) cout << "->  new event " << endl;
    fEvent = jentry;
    // in case of a TChain, ientry is the entry number in the current file
    tsdump = kFALSE;
    Int_t ientry = LoadTree(jentry); 
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%step == 0) cout << "Event " << jentry << endl;
    if (ientry == 0) cout << "File " << fChain->GetCurrentFile()->GetName() << endl;

    //CB cut and paste from recoil and other stuff...
    // -- Initialize event
    initVariables();

    findbestB();
    
    brecoOverlap = 1; // this is flag for overlap with BRECO candidate
    if(indexbestB == 1) brecoOverlap = 2;

    ischB = 0;
    if (bestB0 == 0) {
      ischB = 1;
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

    p4Upsilon = TLorentzVector(pxUps, pyUps, pzUps, eUps); 
    upsBoost = TVector3(pxUps, pyUps, pzUps);
    upsBoost.SetMag(upsBoost.Mag()/eUps);

//  Changed sign of p4Brecoil with respect to original
    p4Brecoil = p4Upsilon - p4Breco; 
    cmsBoost = p4Brecoil.Vect();
    cmsBoost.SetMag(cmsBoost.Mag()/p4Brecoil.E());

    // -- Find leading lepton
    Int_t i;
    pmax = 0;


    for (i = 0; i < nTrk; ++i) {

      Brectrktmp = B0RecTrk[i];
      if(ischB) Brectrktmp = chBRecTrk[i];
      if((Brectrktmp&brecoOverlap)) continue;
      if(fNewFormat == 1 && Brectrktmp==2)continue;// Apr 02 format
      //      if (momentumTrk[i] < PLABLO) continue;
      //      if (goodTrack[i] == 0) continue;
      
      bool iRecLooLept(0);
      if (TMath::Abs(elecIdTrk[i])>7 || (TMath::Abs(muonIdTrk[i]))>7) iRecLooLept = 1;

      if (iRecLooLept) {  // do PID
	mass = MUMASS; 
	if (TMath::Abs(elecIdTrk[i])>7) mass = ELMASS; // electrons override muons
	mk4Vector(p4l, momentumTrk[i], thetaTrk[i], phiTrk[i], mass);
	p4l.Boost(-cmsBoost);
	plab = momentumTrk[i];
	pcms = p4l.Vect().Mag(); 
	
	if (pcms > pmax) {
	  pmax = pcms;
	}
      }
      fisSkim = (pmax >= pCut);	
    }         
    
    if (fOptMakeEventList) {
      if (fisSkim) {
	//	cout << "--> Copying " << jentry << endl;
	fToBeCopied->Enter(jentry);	
	fSkimcount++;
      }
    }
  }
}

void recoilNtp::FastSkim(Int_t maxEvent, Int_t isVerbose, const char *ITSfile) {
  char  buffer[200];
  ifstream is(ITSfile);
  const int nmaxevt=10000;
  TExMap mm;
  int getRun[nmaxevt],getUpper[nmaxevt],getLower[nmaxevt];
  int i(0);
  cout <<" reading file "<<ITSfile<<endl;
  while (is.getline(buffer, 200, '\n')) {
    is>>getRun[i]>>getLower[i]>>getUpper[i];
    mm.Add(getLower[i],getRun[i]);
    //    cout << i << " " << getRun[i]<< " " <<getUpper[i]<<" " <<getLower[i]<<endl;
    i++;
  }

  if (fChain == 0) return;
  Int_t nentries = Int_t(fChain->GetEntries());
  if (maxEvent == 0) maxEvent = nentries;
  if (nentries < 1) {
    cout << "Found no entries in " << fChain->GetName() << endl;
  } else {
    cout << "Found " << nentries << " entries in tree " << fChain->GetName() << endl;
  }

  
  
  for (Int_t jentry = 0; jentry < maxEvent; jentry++) {
    Int_t ientry = LoadTree(jentry); 
    fChain->GetEntry(jentry); 
    fisSkim=false; 
    //    for(int j=0;j<i;j++){
      //      if(runNumber==getRun[j] && upperID==getUpper[j] && lowerID==getLower[j]){
    if(mm.GetValue(lowerID)==runNumber){
	cout <<" found event "<<runNumber<<" " <<upperID<<" " <<lowerID<<endl;
	fisSkim=true; 
    }
	//	break;
	//      } else if(getRun[j]>runNumber){
	//	break;
	//      }
    
    if (fOptMakeEventList) {
      if (fisSkim) {
	//	cout << "--> Copying " << jentry << endl;
	fToBeCopied->Enter(jentry);	
	fSkimcount++;
      }
    }
  }
}

// ----------------------------------------------------------------------
TFile* recoilNtp::openHistFile(TString name) {
  cout << name << endl;
  fHistFile = new TFile(name.Data(), "RECREATE");
  fHistFile->cd();
  cout << "Opened " << fHistFile->GetName() << endl;
  return fHistFile;
}


// ----------------------------------------------------------------------
void recoilNtp::closeHistFile() {
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

// ----------------------------------------------------------------------
void recoilNtp::dumpEventList(const char *filename) {
  TFile *f = new TFile(filename, "RECREATE");
  fToBeCopied->Print();
  fChain->SetEventList(fToBeCopied);
  TTree *small = fChain->CopyTree("");
  small->Write();
  small->Print();
}


#include "bookHist.icc"
#include "fillHist.icc"
#ifdef FAST
#include "fastHist.icc"
#endif

void recoilNtp::KinFit(TLorentzVector& p4XReco, TLorentzVector& p4Nu,
                       TLorentzVector& p4XFit, TLorentzVector& p4NuFit)
{
  int ILTYP; 
  if (fElectron) {
    ILTYP = 2; 
  } else if (fMuon) {
    ILTYP = 1; 
  } 
  
  float CVAL[4] = {pxUps, pyUps, pzUps, eUps};
  float P_REC[28] = {
    p4BrecoNC.Px(), p4BrecoNC.Py(), p4BrecoNC.Pz(), p4BrecoNC.E(), 
    p4LeptonLab.Px(), p4LeptonLab.Py(), p4LeptonLab.Pz(), p4LeptonLab.E(), 
    p4XReco.Px(), p4XReco.Py(), p4XReco.Pz(), p4XReco.E(), 
    p4Nu.Px(), p4Nu.Py(), p4Nu.Pz(), p4Nu.E()};
  float P_FIT[28];
  float CHI2T, PROBCHI2; 
  int   ISMEAR(0);
  int   IERR;
  int ISV = 0;

#if 1
  abcfit_interface_vub_(&ISMEAR,&ILTYP,CVAL,P_REC,P_FIT,&CHI2T,&PROBCHI2,&IERR, &ISV);
#endif

  p4XFit.SetXYZT(P_FIT[8], P_FIT[9], P_FIT[10], P_FIT[11]);
  p4NuFit.SetXYZT(P_FIT[12], P_FIT[13], P_FIT[14], P_FIT[15]); 
}

void recoilNtp::MxStudy(int chbcand) 
{
  fMM2Ks=100.;
  fMxFitKs=fMxKs=0.;
  fGoodEventKS=kFALSE;
  if (fNKshort>0) {
    TLorentzVector pXKS(0,0,0,0), pMissKS(0,0,0,0);
    TLorentzVector p4XKSFit(0,0,0,0),pMissKSFit(0,0,0,0);
    Bool_t ksStatus = MxhadKS(fMxKs,fMM2Ks,pXKS,pMissKS);
    Bool_t GoodMM2KS(kFALSE);
    if ((ksStatus) && (MM2LO < fMM2Ks) && (fMM2Ks < MM2HI)) GoodMM2KS = kTRUE;
    if (fGoodLepton && fOneLepton && GoodMM2KS && fGoodChargeCons && fGoodChargeCorr) fGoodEventKS = kTRUE;
    if ((fPcms>0.) && (-100. < fMM2) && (fMM2 < 100.)) {
      KinFit(pXKS,pMissKS,p4XKSFit,pMissKSFit);
      fMxFitKs=p4XKSFit.M();
    }
  }


  fMxMisK=0;
  fMM2MisK=100.0;
  fMxFitchK=fMxchK=0;
  if (fNKp>0) {
    TLorentzVector pXMischK(0,0,0,0), pMissMischK(0,0,0,0);
    TLorentzVector pXMischKFit(0,0,0,0), pMissMischKFit(0,0,0,0);
    MxMisIdchK(fMxMisK,fMM2MisK,pXMischK, pMissMischK);
    if ((fPcms>0.) && (-100. < fMM2) && (fMM2 < 100.)) {
      KinFit( pXMischK,pMissMischK, pXMischKFit, pMissMischKFit);
      fMxFitMisK= pXMischKFit.M();
    }
    TLorentzVector pXchK(0,0,0,0), pMisschK(0,0,0,0);
    TLorentzVector pXchKFit(0,0,0,0), pMisschKFit(0,0,0,0);
    fchkStatus = MxHadchK(chbcand, fMxchK, fMM2chK, pXchK, pMisschK);
    if ((fPcms>0.) && (-100. < fMM2) && (fMM2 < 100.)) {
      KinFit(pXchK ,pMisschK, pXchKFit, pMisschKFit);
      fMxFitchK= pXchKFit.M();
    }
  }

  if (fIsMC) {

    fsumklemcen0=fsumklemcen22=fsumklemcen=0.0;
    fklreslen=0;
    TList * listkl = createTrueList(130);
    TList * listks = createTrueList(310,ftksdecay,20);
    TList * listchk1 = createTrueList(321);
    TList * listchk = createTrueList(-321);
    listchk->AddAll( listchk1 );
    fKLemcp4.SetXYZM(0,0,0,0);
  
    for (Int_t i=0; i<nGam; i++) {
      if (energyGam[i]>GAMMAELO) {
	TLorentzVector tmp(0,0,0,0);
	mk4Vector(tmp, energyGam[i], thetaGam[i], phiGam[i], 0.);
	Int_t match = angularMatch(*listkl,i);
	if (match>0) {
	  fKLemcp4 += tmp;
	  fsumklemcen+= energyGam[i];
	  if ( idGam[i] == 22) {
	    fsumklemcen22+=energyGam[i];
	  }
	  if ( idGam[i] == 0) {
	    fsumklemcen0+= energyGam[i];
	  }
	}
      }
    }
    
    
    // prepare ntuple variables
    fnKL = listkl->GetSize();
    fntks = listks->GetSize();
    fntchk = listchk->GetSize();
  
    TIterator*  next = listkl->MakeIterator();
    TLorentzVector *p4;
    Int_t counter = 0;
    while (( (p4 = ((TLorentzVector*) (*next)())) != 0 ) ) {
      if ( counter == 21 ) break;
      ftklp[counter] = p4->P();
      ftklth[counter] = p4->Theta();
      ftklph[counter] = p4->Phi();
      ftklisol[counter]=kTRUE;
      for (Int_t k=0; k<nMc; k++) {
	// this cannot work, we should loop only on stable particles.
	if ( TMath::Abs(p4->Phi()-phiMc[k]) < .200 ) 
	  if ( TMath::Abs(p4->Theta()-thetaMc[k]) < .200 ) 
	    ftklisol[counter]=kFALSE;
      }
      counter++;
    }
    
    next = listks->MakeIterator();
    counter = 0;
    while (( (p4 = ((TLorentzVector*) (*next)())) != 0 ) ) {
      if ( counter == 21 ) break;
      ftksp[counter] = p4->P();
      ftksth[counter] = p4->Theta();
      ftksph[counter] = p4->Phi();
      //    note : ftksdecay already filled
      counter++;

    }
    
    next = listchk->MakeIterator();
    
    counter = 0;
    while (( (p4 = ((TLorentzVector*) (*next)())) != 0 ) ) {
      if ( counter == 21 ) break;
      ftchkp[counter] = p4->P();
      ftchkth[counter] = p4->Theta();
      ftchkph[counter] = p4->Phi();
      counter++;
    }
    
    fksmatchp = 0;
    fpks = 0;
    Int_t iks = bestKsIndex(0);
    if (iks > -1) {
      fpks = pKs[iks];
      Int_t ipar = isMcKs(iks);
      if (ipar > -1) 
	fksmatchp=pMc[ipar];
      else
	fksmatchp=0;
    }
    
    Int_t j(0);
    for (Int_t i=0; i<nKs; i++) {
      if (goodKshort[i] == 1 && j<20) {
	fallKsm0[j]=massKs[i];
	fallKsp[j]=pKs[i];
	fallKsMc[j] = ( isMcKs(i) > -1 ? kTRUE : kFALSE );
	j++;
      }
    }
    
  j=0;
  for (Int_t i=0; i<nTrk; i++) {
    if (isRecKaon(i) && (! isRecEl(i) ) && j<20) {
      fallchkp[j]=momentumTrk[i];
      fallchkMc[j] = isMcKaon(i);
      j++;
    }
  }
  
  delete listkl;
  delete listks;
  delete listchk1;
  delete listchk;
  }
}


Int_t recoilNtp::isMcKs( int iks ) 
{
  Int_t i1 = IndexTrk[d1KsIndex[iks]-1];
  Int_t i2 = IndexTrk[d2KsIndex[iks]-1];
  if (i1>0 && i2>0 && i1<201 && i2<201) {
    Int_t ipar = mothMc[i1-1];
    if ( ipar == mothMc[i2-1] )
      if ( idMc[ipar-1] == 310 )
        return (ipar-1);
  }
  return -1;
}

