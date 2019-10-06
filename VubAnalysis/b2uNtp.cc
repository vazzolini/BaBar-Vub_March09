#include <iostream.h>
#include <fstream.h>
#include "TRandom.h"

#include "TH1.h"

#include "VubAnalysis/b2uNtp.hh"

Double_t p_energy_loss_corrected(Double_t pin, Double_t dip, Int_t itype);


#include "VubAnalysis/b2uStudies.icc"
#include "VubAnalysis/b2uMcTruth.icc"



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


// ----------------------------------------------------------------------
void b2uNtp::Loop(int maxEvent, int startEvent) {

  int step(1000);
  double tmpPgen, tmpThetagen, tmpPhigen, tmpMassgen ;

  if (fChain == 0) return;
  int nentries = int(fChain->GetEntries());
  if (maxEvent == 0) maxEvent = nentries;
  if (nentries < 1) {
    cout << "xx> b2uNtp::Loop> Found no entries in " << fChain->GetName() << endl;
  } else {
    cout << "==> b2uNtp::Loop> Found " << nentries << " entries in tree " << fChain->GetName() << endl;
  }

  if (startEvent > 0) {
    cout << "==> b2uNtp::Loop> Will start at event " << startEvent << endl;
    if (startEvent+maxEvent >  nentries) {
      cout << "xx> b2uNtp::Loop> Requested " << maxEvent << " events, but will run only to end of chain"  << endl;
      maxEvent = nentries - startEvent; 
    }
  }

  Int_t nbytes(0), nb(0), ientry(0);
  const char *pChar; 

  for (int jentry = startEvent; jentry < startEvent+maxEvent; jentry++) {
    fEvent = jentry;
    // -- in case of a TChain, ientry is the entry number in the current file
    if (jentry%step == 0) cout << "==> b2uNtp::Loop> Event " << jentry << endl;
    ientry = LoadTree(jentry);  nb = fChain->GetEntry(jentry);   nbytes += nb;
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
    // -- that's it for ROOT 


    // -- ANALYSIS LOOP 
    // ----------------
    if (fVerbose) cout << "----------------------------------------------------------------------" << endl;
    if (fVerbose) cout << "==> b2uNtp::Loop>  jentry = " << jentry << endl;

    if (fIsMC && fOptSmearNeut > 0) {
      smearNeut();
    }

    mcTruth(); // -- Protection against data and var initialization in mcTruth() itself
    if (fIsMC) {
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
      mk4Vector(f4BrecoGen, tmpPgen, tmpThetagen, tmpPhigen, tmpMassgen);
    }

    // -- Breco, \FourS, and boost vectors
    doBreco(); 
    if (fIndexBestB < 0) {
      if (fVerbose) cout << "==> b2uNtp::Loop> No good BRECO found, skipping " << endl;
      continue; 
    }
    
    // -- Prepare the recoil
    cleanupTracks(); 
    selectTracks();
    doSplitOffStudy();
    selectPhotons();
    selectKlongs();
    selectKshorts();
    selectKshortsZZ(); 

    fillPidMaps(); 

    // -- The analysis
    doBremRecovery(); 
    findLeadingLepton();
    if (fPcms < 0.) {
      if (fVerbose) cout << "==> b2uNtp::Loop> No leading lepton found, skipping " << endl;
      //      continue;
    }
    recoil(); 

    //For the Breco quality
    goodBreco(brecoI);

    // -- partial reconstruction, etc. These may need the full X system, so keep them AFTER recoil()!
    //    NOTE: The check for correct Breco charge is in the functions!
    doPartialB0();   
    doPartialB0bis();
    doPartialBp();   

    doDeltaM();

    // -- Some Studies
    //     doKshorts(); 
    //     doKszz(); 
    //     doKlongs(); 
    //    doTracks(); 

    // -- Dump tree
    if (fDump > 0) {
      //      if ((fDump & 4) && ((fPcms > 0.) || (fVub == 1) || (fVcb == 1))) {
      if (fDump & 4) {
        if (fVerbose) cout << "==> b2uNtp::Loop> fill fTree at dumpLevel 4" << endl;
        fTree->Fill();
      }
      if (fDump & 8) {
        if (fVerbose) cout << "==> b2uNtp::Loop> fill fTree at dumpLevel 8" << endl;
        fTree->Fill();
      }

    }



  }


}


// ----------------------------------------------------------------------
void b2uNtp::doBremRecovery() {

  TLorentzVector p4ele(0., 0., 0., 0.); 
  TLorentzVector p4t(0., 0., 0., 0.), p4g(0., 0., 0., 0.);
  double a(99.); 

  for (int itrk = 0; itrk < nTrk; ++itrk) {
    if (fBrecoTrk[itrk])         continue; 
    if (fGoodPionTrk[itrk] == 0) continue;
    if (isRecEl(itrk) && (fGoodElectronTrk[itrk] == 1)) {  
      // got a good electron
      mk4Vector(p4ele, momentumTrk[itrk], thetaTrk[itrk], phiTrk[itrk], ELMASS);
      // initialize the 4 vectors to 0
      mk4Vector(p4t, 0., 0., 0., 0.);
      mk4Vector(p4g, 0., 0., 0., 0.);

      // -- Collect all photons within Ed's cone of 0.08
      for (int ig = 0; ig < nGam; ++ig) {
	if (ecalGam[ig] < 0.) continue; 
	mk4Vector(p4t, energyGam[ig], thetaGam[ig], phiGam[ig], 0.);
	a = p4t.Angle(p4ele.Vect()); 
	if (a < ELBREMCONE) {
	  fGoodPhoton[ig] = 0; 
	  //cout << "--> b2uNtp::doBremRecovery> Reset goodPhoton because it's from bremsstrahlung " << ecalGam[ig] << endl;
	  p4g += p4t; 
	}
      }

      // -- update electron momentum
      p4ele += p4g; 
      momentumTrk[itrk] = p4ele.Vect().Mag(); 
      thetaTrk[itrk] = p4ele.Vect().Theta(); 
      phiTrk[itrk] = p4ele.Vect().Phi(); 

    }
  }
  
}

// ----------------------------------------------------------------------
void b2uNtp::findLeadingLepton() {

  // -- Reset data members
  // ---------------------
  fLeptonIndex = -99; 
  fLeptonCharge = fNlepton = fElectron = fMuon = fNEl = fNMu = 0; 
  fPlab = fTlab = fTlabDR = fFlab = fFlabDR = fPcms = fTcms = fFcms = fEcms = -99.;
  f4Lepton.SetXYZM(0., 0., 0., 0.); 
  f4LeptonCMS.SetXYZM(0., 0., 0., 0.); 

  // -- Find leading lepton
  TLorentzVector p4l(0., 0., 0., 0.); 
  double mass(0.), lmass(0.), pcms(-100.), pmax(-99.); 
  for (int i = 0; i < nTrk; ++i) {
    if (fBrecoTrk[i])         continue; 
    if (fGoodPionTrk[i] == 0) continue;

    if (isRecLepton(i)) {  // do PID
      mass = MUMASS; 
      if (isRecEl(i)) mass = ELMASS; // electrons override muons
      mk4Vector(p4l, momentumTrk[i], thetaTrk[i], phiTrk[i], mass);
      p4l.Boost(-f3CmsBoost);
      pcms = p4l.Vect().Mag(); 

      fNlepton++;

      if (pcms > pmax) {
        pmax = pcms;
        lmass = mass;
        fLeptonIndex = i;
      }
    }
  }

  // -- Set up leading lepton
  if (fLeptonIndex  > -1) {
    setupLeadingLepton(lmass); 

    if (lmass < 0.1) {
      fElectron = 1;
      fNEl = 1;
    }
    if (lmass > 0.1) {
      fMuon = kTRUE;
      fNMu = 1;
    }
    // -- remove overlap: misid'ed electrons (PidKilling) are automatically also misid'ed muons
    if (fElectron && fMuon) {
      fNMu = 0; 
      fMuon = kFALSE;
    }

    if (fVerbose) cout << "==> b2uNtp::findLeadingLepton> Found a leading lepton at " << fLeptonIndex 
		       << " with pcms = " << fPcms 
		       << " and mass = " << lmass 
		       << endl;

  }

}


// ----------------------------------------------------------------------
void b2uNtp::setupLeadingLepton(double lmass) {
  TLorentzVector p4l(0., 0., 0., 0.); 
  mk4Vector(p4l, momentumTrk[fLeptonIndex], thetaTrk[fLeptonIndex], phiTrk[fLeptonIndex], lmass);
  f4Lepton      = p4l;
  p4l.Boost(-f3CmsBoost);
  f4LeptonCMS   = p4l;
  fPlab         = momentumTrk[fLeptonIndex];
  fTlab         = thetaTrk[fLeptonIndex];
  fTlabDR       = thetaTrk[fLeptonIndex]*DR;
  fFlab         = phiTrk[fLeptonIndex];
  fFlabDR       = phiTrk[fLeptonIndex]*DR;
  fPcms         = p4l.Vect().Mag();
  fTcms         = p4l.Theta();
  fFcms         = p4l.Phi();
  fEcms         = p4l.E();
  fLeptonCharge = chargeTrk[fLeptonIndex];
}


// ----------------------------------------------------------------------
void b2uNtp::recoil() {

  // -- Reset data members
  // ---------------------
  fB2uDepleted = 0; 
  fNKs = fNKz = fNKl = fNKp = fRecoilTrkMult = fRecoilGamMult = 0; 
  fQtot = fRecoilCharge = 0; 
  fPxhad = fTxhad = fFxhad = fExhad = fMxhad = fMxtrk = fMxnut = -99.;
  fPNu = fTNu = fCosTNu = fFNu = fMM2 = fEmiss = fQ2 = fMM2NC = -99.;
  fPcmsTrkLo = 99.;

  f4Recoil.SetXYZM(0., 0., 0., 0.); 
  f4Xhad.SetXYZM(0., 0., 0., 0.); 
  f4Xtrk.SetXYZM(0., 0., 0., 0.); 
  f4Xnut.SetXYZM(0., 0., 0., 0.); 
  f4Neutrino.SetXYZM(0., 0., 0., 0.); 
  f4NeutrinoNC.SetXYZM(0., 0., 0., 0.); 


  fHistFile->cd();
  TH1D *h1; 
  static int first(1); 
  if (first) {
    first = 0; 
    h1 = new TH1D("lepton", "", 30, 0., 3.); 
    h1 = new TH1D("photon", "", 30, 0., 3.); 
    h1 = new TH1D("nphoton", "", 30, 0., 30.); 
  }

  // -- Kshorts
  for (int i = 0; i < nKs; ++i) {
    if (fGoodKs[i] == 0) continue;   // skip suboptimal kshorts 
    if (TMath::Abs(d1KsLund[i]) == 211) {
      ++fNKs;
    }
    if (TMath::Abs(d1KsLund[i]) == 111) {
      ++fNKz;
    }
    //b2u    fB2uDepleted = 1;
  }

  // -- Tracks 
  TLorentzVector p4t(0., 0., 0., 0.); 
  for (int i = 0; i < nTrk; ++i) {
    if (fBrecoTrk[i])         continue; 
    if (fGoodPionTrk[i] == 0) continue;
    fRecoilCharge += chargeTrk[i];
    ++fRecoilTrkMult;

    if (isRecEl(i) && (fGoodElectronTrk[i] == 1)) {  
      if (i == fLeptonIndex) {
        mk4Vector(p4t, momentumTrk[i], thetaTrk[i], phiTrk[i], ELMASS);
      } else{
        mk4Vector(p4t, momentumTrk[i], thetaTrk[i], phiTrk[i], PIPMASS);
      }
    } else if (isRecKaon(i)  && (fGoodKaonTrk[i] == 1)) {  
      fB2uDepleted = 1;
      ++fNKp;
      double momentum = p_energy_loss_corrected(momentumTrk[i], 0.5*TMath::Pi() - thetaTrk[i], 1);
      mk4Vector(p4t, momentum, thetaTrk[i], phiTrk[i], KAPMASS);
      momentumTrk[i] = momentum;
    } else if (isRecMu(i)  && (fGoodMuonTrk[i] == 1)) {
      if (i == fLeptonIndex) {
        mk4Vector(p4t, momentumTrk[i], thetaTrk[i], phiTrk[i], MUMASS);
      } else{
        mk4Vector(p4t, momentumTrk[i], thetaTrk[i], phiTrk[i], PIPMASS);
      }
    } else {
      mk4Vector(p4t, momentumTrk[i], thetaTrk[i], phiTrk[i], PIPMASS);
    }

    if (fVerbose) cout << "Adding track " << i << " p = " << momentumTrk[i] << " and mass = " << p4t.Mag() << endl;

    f4Recoil += p4t;

    if (i == fLeptonIndex) continue;

    f4Xhad += p4t;
    f4Xtrk += p4t;

    p4t.Boost(-f3CmsBoost);
    if (p4t.Vect().Mag() < fPcmsTrkLo) {
      fPcmsTrkLo = p4t.Vect().Mag();
    }     

  }



  // -- Photons
  int nphoton(0); 
  double emax(-99.); 
  for (int i = 0; i < nGam; ++i) {
    if (fBrecoPhoton[i])     continue; 
    if (fGoodPhoton[i] == 0) continue;

    if (fGoodKl[i] == 1) {
      fNKl++; 
    }

    ++fRecoilGamMult; 

    if (fVerbose) cout << "Adding photon " << i << " E = " << energyGam[i] << endl;
    
    mk4Vector(p4t, energyGam[i], thetaGam[i], phiGam[i], 0.);
    f4Recoil += p4t;
    f4Xhad += p4t;
    f4Xnut += p4t;

    ++nphoton; 
    p4t.Boost(-f3CmsBoost);
    if (p4t.Vect().Mag() > emax) {
      emax = p4t.Vect().Mag();
    }
  }

  ((TH1D*)gFile->Get("lepton"))->Fill(fPcms);
  ((TH1D*)gFile->Get("photon"))->Fill(emax);
  ((TH1D*)gFile->Get("nphoton"))->Fill(nphoton);

  // -- Y system
  TLorentzVector p4y = f4Xhad + f4Lepton; 
  p4y.Boost(-f3CmsBoost);



  // -- Neutrino, etc
  fPxhad   = f4Xhad.Vect().Mag(); 
  fTxhad   = f4Xhad.Theta(); 
  fFxhad   = f4Xhad.Phi(); 
  fExhad   = f4Xhad.E();
  fMxhad   = f4Xhad.Mag();
  fMxtrk   = f4Xtrk.Mag();
  fMxnut   = f4Xnut.Mag();


  fQtot    = fRecoilCharge + fBrecoCharge;  

  f4Neutrino   = f4Upsilon - f4Breco - f4Recoil;
  f4NeutrinoNC = f4Upsilon - f4BrecoNC - f4Recoil;

  fPNu   = f4Neutrino.Vect().Mag(); 
  fTNu   = f4Neutrino.Theta(); 
  fCosTNu= f4Neutrino.CosTheta(); 
  fFNu   = f4Neutrino.Phi(); 
  fMM2   = f4Neutrino.Mag2();
  fEmiss = f4Neutrino.E();
  fQ2    = 2.*f4Neutrino*f4Lepton;
  fMM2NC = f4NeutrinoNC.Mag2();
  fCosBY = 0.; 


  // -- Kinematic fitting
#if 1
  float CVAL[4] = {pxUps, pyUps, pzUps, eUps};
  float P_REC[28] = {
    f4BrecoNC.Px(), f4BrecoNC.Py(), f4BrecoNC.Pz(), f4BrecoNC.E(), 
    f4Lepton.Px(), f4Lepton.Py(), f4Lepton.Pz(), f4Lepton.E(), 
    f4Xhad.Px(), f4Xhad.Py(), f4Xhad.Pz(), f4Xhad.E(), 
    f4NeutrinoNC.Px(), f4NeutrinoNC.Py(), f4NeutrinoNC.Pz(), f4NeutrinoNC.E()};
  
  float P_FIT[28] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
		     0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
  
  float CHI2T, PROBCHI2; 
  
  int ISMEAR = 102;
  int IERR;
  int ISV = fVerbose; 

  int ILTYP(0); 
  if (fElectron) {
    ILTYP = 2; 
  } else if (fMuon) {
    ILTYP = 1; 
  }

  if (ILTYP > 0) {

    int i1 = abcfit_interface_vub_(&ISMEAR,&ILTYP,CVAL,P_REC,P_FIT,&CHI2T,&PROBCHI2,&IERR, &ISV);
    i1 = 0; 
    
    fProbChi2 = PROBCHI2;
    fChi2     = CHI2T; 
    TLorentzVector pxhadfit(P_FIT[8], P_FIT[9], P_FIT[10], P_FIT[11]); 
    fMxhadfit = pxhadfit.Mag();
    TLorentzVector pnufit(P_FIT[12], P_FIT[13], P_FIT[14], P_FIT[15]); 
    fMM2fit = pnufit.Mag2();
    fQ2fit  = 2.*pnufit*f4Lepton;
    
    if (!((fMxhadfit>0) || (fMxhadfit<0) || (fMxhadfit == 0))) {
      cout << "MEZZEGA: NAN in MXHADFIT " << IERR << endl;
      fMxhadfit = -999.;
      fPcms = -97.; // reset fPcms so that the event is not dumped into 'events'
    }   
  }

#endif


  if (fVerbose)  cout << "====================> " << fMxhad <<  " mm2 = " << fMM2 << "  " << fMM2NC << endl;

}




// ----------------------------------------------------------------------
void b2uNtp::setup(const char *histFileName, int dumpLevel) {

  // -- Define constants
  BZLUND = 511;
  BPLUND = 521;

  ELMASS  = 0.000511;
  MUMASS  = 0.10567;

  PIPMASS = 0.13957; 
  PIZMASS = 0.13498;
  KAPMASS = 0.49368;
  KAZMASS = 0.49767;
  BPMASS  = 5.2790;
  BZMASS  = 5.2794;
  BMASS   = 5.2792;

  BQMASS  = 4.800;
  A0      = 1.29;

  DR      = 57.29577951;

  // -- Open ROOT output file 
  fHistFile = new TFile(histFileName, "RECREATE");
  fHistFile->cd();
  cout << "==> b2uNtp::setup> Opened " << fHistFile->GetName() << endl;
       
  // -- Book Trees
  if (dumpLevel > 0) {
    cout << "==> b2uNtp::setup> Booking events tree" << endl;
    fDump = dumpLevel;
    fTree = new TTree("events", "events"); 
    fTree->Branch("run",    &runNumber, "run/I");

    // -- BRECO
    fTree->Branch("mes",    &fMes, "mes/D");
    fTree->Branch("intpur",     &fIntPurity, "intpur/D");
    fTree->Branch("brecoflav",  &fBrecoFlavor, "brecoflav/I");
    fTree->Branch("brecocharge",&fBrecoCharge , "brecocharge/I");

    // -- recoil
    fTree->Branch("xcharge", &fRecoilCharge, "xcharge/I"); // note that this is the entire recoil, not xhad!
    fTree->Branch("mxhad",   &fMxhad, "mxhad/D");
    fTree->Branch("mxtrk",   &fMxtrk, "mxtrk/D");
    fTree->Branch("mxnut",   &fMxnut, "mxnut/D");
    fTree->Branch("mxhadfit",&fMxhadfit, "mxhadfit/D");
    fTree->Branch("chi2",    &fChi2, "chi2/D");
    fTree->Branch("probchi2",&fProbChi2,"probchi2/D");
    fTree->Branch("pcmstrklo",&fpcmsTrklo,"pcmstrklo/D");

    // -- lepton
    fTree->Branch("lcharge", &fLeptonCharge ,"lcharge/I");
    fTree->Branch("pcms",    &fPcms, "pcms/D");
    fTree->Branch("tcms",    &fTcms, "tcms/D");
    fTree->Branch("tlab",    &fTlab, "tlab/D");
    fTree->Branch("plab",    &fPlab, "plab/D");

    // -- missing momentum (aka neutrino)
    fTree->Branch("mm2",      &fMM2,   "mm2/D");
    fTree->Branch("mm2fit",   &fMM2fit,"mm2fit/D");
    fTree->Branch("q2fit",    &fQ2,    "q2fit/D");
    fTree->Branch("pmiss",    &fPNu,   "pmiss/D");
    fTree->Branch("tmiss",    &fTNu,   "tmiss/D");
    fTree->Branch("emiss",    &fEmiss, "emiss/D");
    fTree->Branch("costmiss", &fCosTNu,"costmiss/D");

    // -- Event
    fTree->Branch("nle",    &fNlepton, "nle/I");
    fTree->Branch("nel",    &fNEl, "nel/I");
    fTree->Branch("nmu",    &fNMu, "nmu/I");

    fTree->Branch("nchg",   &fRecoilTrkMult, "nchg/I");
    fTree->Branch("nneu",   &fRecoilGamMult, "nneu/I");
    fTree->Branch("nks",    &fNKs, "nks/I");
    fTree->Branch("nkz",    &fNKz, "nkz/I");
    fTree->Branch("nkl",    &fNKl, "nkl/I");
    fTree->Branch("nkp",    &fNKp, "nkp/I");

    fTree->Branch("deltam",&fDeltaM, "wdeltam/D");
    fTree->Branch("mm1pr",  &fMM1pr, "mm1pr/D");
    fTree->Branch("mm2pr",  &fMM2pr, "mm2pr/D");
    fTree->Branch("mm3pr",  &fMM3pr, "mm3pr/D");

    fTree->Branch("oa1",    &fOA1, "oa1/D");
    fTree->Branch("oa2",    &fOA2, "oa2/D");
    fTree->Branch("oa3",    &fOA3, "oa3/D");

    fTree->Branch("pcmstrklo", &fPcmsTrkLo, "pcmstrklo/D");

    // -- Generator level
    fTree->Branch("vub",    &fVub, "vub/I");
    fTree->Branch("vcb",    &fVcb, "vcb/I");
    fTree->Branch("other",  &fOther, "other/I");
    fTree->Branch("qb",     &fQb,    "qb/I");
    fTree->Branch("kplus",  &fKplus, "fkplus/D");
    
    fTree->Branch("ecmsgen", &fEcmsGen, "ecmsgen/D");
    fTree->Branch("tcmsgen", &fTcmsGen, "tcmsgen/D");
    fTree->Branch("fcmsgen", &fFcmsGen, "fcmsgen/D");
    fTree->Branch("mxhadgen",&fMxhadGen, "mxhadgen/D");
    fTree->Branch("q2Gen",   &fQ2Gen, "q2Gen/D");

    fTree->Branch("ctvgen",  &fctvGen, "ctvgen/D");
    fTree->Branch("ctlgen",  &fctlGen, "ctlgen/D");
    fTree->Branch("chigen",  &fchiGen, "chigen/D");
    fTree->Branch("pxhadgen",&fPxhadGen, "pxhadgen/D");
    fTree->Branch("wKK",     &fWithKKbar, "wKK/I");
    fTree->Branch("mxhadgenwoph", &fMxhadGenwoPh, "mxhadgenwoph/D");
    fTree->Branch("brecoqual", &fBrecoQual, "brecoqual/I");
    fTree->Branch("chgdaugen", &fchgDau, "chgdaugen/I");
    fTree->Branch("neudaugen", &fneuDau, "neudaugen/I");

    // -- neutrino
    fTree->Branch("vpgen", &fPvGen, "vpgen/D");
    fTree->Branch("vtgen", &fTvGen, "vtgen/D");
    fTree->Branch("vfgen", &fFvGen, "vfgen/D");

    fTree->Branch("vpcmsgen", &fPvcmsGen, "vpcmsgen/D");
    fTree->Branch("vtcmsgen", &fTvcmsGen, "vtcmsgen/D");
    fTree->Branch("vfcmsgen", &fFvcmsGen, "vfcmsgen/D");
    
    fTree->Branch("Gvxbtyp", &fBVxbTyp, "Gvxbtyp/I");
    fTree->Branch("GfDpi", &fDpi, "GfDpi/I");
    fTree->Branch("GfDk", &fDk, "GfDk/I");
    fTree->Branch("GfDks", &fDks, "GfDks/I");
    fTree->Branch("GfDpiz", &fDpiz, "GfDpiz/I");
    fTree->Branch("GfDlep", &fDlep, "GfDlep/I");
 

  }

  // -- Book Histograms
  TH1D *h1; 
  h1 = new TH1D("h1", "mes", 40, 5.2, 5.3); 

}




// ----------------------------------------------------------------------
void b2uNtp::doBreco() {
  // This is from findbestB() and parts of the old recoil()
  // and skipBadBreco() 


  // -- Reset data members
  // ---------------------
  fChB = 0; // This IS the default!
  fSeedMode = fBmode = -99; 
  fIndexBestB = fBrecoOverlap = -99; 
  fIntPurity = fPurity = -99.;
  fMes = fPcmsBreco = -99.;
  fBrecoFlavor = fBrecoCharge = -99; 

  f4Breco.SetXYZM(0., 0., 0., 0.); 
  f4BrecoNC.SetXYZM(0., 0., 0., 0.); 
  f4Upsilon.SetXYZM(0., 0., 0., 0.); 
  f4Recoil.SetXYZM(0., 0., 0., 0.); 

  f3CmsBoost.SetXYZ(0., 0., 0.); 
  f3UpsBoost.SetXYZ(0., 0., 0.); 


  // -- Create histograms and read in a priori purity tables
  // -------------------------------------------------------
  static int first(1); 
  static double brecosig[6000], brecobkg[6000], brecointpur[6000];  
  
  TH1D *h; 
  if (first) {
    first = 0;
    cout << "FIXME!!!!! b2uNtp:doBreco: skipBadBreco() not rejecting anything!" << endl;
    
    char buffer[200];
    float bmode, dmode, sig, bkg, pur, sb;
    ifstream is("tables/tablepurity.dat");
    int mode;
    cout << "==> b2uNtp::doBreco> reading purity tables from " << " tables/tablepurity.dat"; 
    int count(0); 
    while (is.getline(buffer, 200, '\n')) {
      if (buffer[0] == '#') {continue;}
      sscanf(buffer, "%f %f %f %f %f %f", &bmode, &dmode, &sig, &bkg, &pur, &sb);
      mode = (dmode+100) * 100 + bmode-10000;
      brecosig[mode] = sig;      
      brecobkg[mode] = bkg;
      brecointpur[mode] = pur; 
      ++count; 
    }
    cout << " ... a total of " << count << " lines " << endl;

    
    if (fHistFile != 0) {
      char title[200];
      char name[200];
      fHistFile->cd();
      fHistFile->mkdir("doBreco", "doBreco");
      fHistFile->cd("doBreco");
      sprintf(name, "meszero");  sprintf(title, "mes zero");  h = new TH1D(name, title, 40, 5.2, 5.3);
      sprintf(name, "mesbest");  sprintf(title, "mes best");  h = new TH1D(name, title, 40, 5.2, 5.3);
      
      sprintf(name, "h101");  sprintf(title, "int purity B0");  h = new TH1D(name, title, 500, 0., 1.); 
      sprintf(name, "h111");  sprintf(title, "int purity B+");  h = new TH1D(name, title, 500, 0., 1.); 

      sprintf(name, "mesBz");  sprintf(title, "B0 mes");  h = new TH1D(name, title, 40, 5.2, 5.3);
      sprintf(name, "mesBp");  sprintf(title, "Bp mes");  h = new TH1D(name, title, 40, 5.2, 5.3);

      sprintf(name, "mesBzFailed");  sprintf(title, "B0 mes (failed purity cuts)");  h = new TH1D(name, title, 40, 5.2, 5.3);
      sprintf(name, "mesBpFailed");  sprintf(title, "Bp mes (failed purity cuts)");  h = new TH1D(name, title, 40, 5.2, 5.3);

    }
  }


  // -- Determine mode with highest purity
  // -------------------------------------
  double tmpintpur(-999);
  for (int iB0=0; iB0<nB0; iB0++){
    int mode = modeB0[iB0]-10000;
    if(brecointpur[mode]>tmpintpur) {
      fIndexBestB = iB0;    
      tmpintpur = brecointpur[mode];
    }
  }
  
  for (int iChB=0; iChB<nChB; iChB++){
    int mode = modeChB[iChB]-10000;
    if(brecointpur[mode]>tmpintpur) {
      fIndexBestB = iChB;    
      tmpintpur = brecointpur[mode];
      fChB = 1; 
    }
  }
  
  // -- BRECO characteristics
  int modeB; 
  if (fChB) {
    fMes = mseChB[fIndexBestB];
    modeB = modeChB[fIndexBestB]; 
  } else {
    fMes = mseB0[fIndexBestB];
    modeB = modeB0[fIndexBestB]; 
  }

  fIntPurity = brecointpur[modeB-10000];
  fPurity    = brecosig[modeB-10000]/(brecosig[modeB-10000]+brecobkg[modeB-10000]);
  fBmode     = modeB; 

  // -- Yuck
  fBrecoCharge = 0;
  fBrecoFlavor = 1;  

  if (fChB) {
    if(d2ChBLund[fIndexBestB]!=0&&d2ChBLund[fIndexBestB]!=111&&d2ChBLund[fIndexBestB]!=310){ fBrecoCharge=fBrecoCharge+(TMath::Abs(d2ChBLund[fIndexBestB])/d2ChBLund[fIndexBestB]);}        
    if(d3ChBLund[fIndexBestB]!=0&&d3ChBLund[fIndexBestB]!=111&&d3ChBLund[fIndexBestB]!=310){ fBrecoCharge=fBrecoCharge+(TMath::Abs(d3ChBLund[fIndexBestB])/d3ChBLund[fIndexBestB]);}
    if(d4ChBLund[fIndexBestB]!=0&&d4ChBLund[fIndexBestB]!=111&&d4ChBLund[fIndexBestB]!=310){ fBrecoCharge=fBrecoCharge+(TMath::Abs(d4ChBLund[fIndexBestB])/d4ChBLund[fIndexBestB]);}
    if(d5ChBLund[fIndexBestB]!=0&&d5ChBLund[fIndexBestB]!=111&&d5ChBLund[fIndexBestB]!=310){ fBrecoCharge=fBrecoCharge+(TMath::Abs(d5ChBLund[fIndexBestB])/d5ChBLund[fIndexBestB]);}
    if(d6ChBLund[fIndexBestB]!=0&&d6ChBLund[fIndexBestB]!=111&&d6ChBLund[fIndexBestB]!=310){ fBrecoCharge=fBrecoCharge+(TMath::Abs(d6ChBLund[fIndexBestB])/d6ChBLund[fIndexBestB]);}
    if(d7ChBLund[fIndexBestB]!=0&&d7ChBLund[fIndexBestB]!=111&&d7ChBLund[fIndexBestB]!=310){ fBrecoCharge=fBrecoCharge+(TMath::Abs(d7ChBLund[fIndexBestB])/d7ChBLund[fIndexBestB]);}        
    fBrecoFlavor=fBrecoCharge; 
  }else{
    fBrecoCharge=0;
    fBrecoFlavor=-1*(TMath::Abs(d1B0Lund[fIndexBestB])/d1B0Lund[fIndexBestB]);
  }


  // -- Define seed modes
  if ((11000 <= modeB) &&  (modeB < 12000)) {
    fSeedMode = 2;
  } else if  ((12000 <= modeB) &&  (modeB < 13000)) {
    fSeedMode = 0;
  } else if  ((13000 <= modeB) &&  (modeB < 14000)) {
    fSeedMode = 1;
  } else if  ((14000 <= modeB) &&  (modeB < 16000)) {
    fSeedMode = 3;
  } else {
    fSeedMode = -1;
  }


  if (fHistFile != 0) {
    fHistFile->cd("doBreco");
    ((TH1D*)gDirectory->Get("meszero"))->Fill(mseB0[0]);
    ((TH1D*)gDirectory->Get("meszero"))->Fill(mseChB[0]);
    if (!fChB) { 
      ((TH1D*)gDirectory->Get("mesbest"))->Fill(mseB0[fIndexBestB]);
    } else{
      ((TH1D*)gDirectory->Get("mesbest"))->Fill(mseChB[fIndexBestB]);
    }

    if (!fChB) {
      ((TH1D*)gDirectory->Get("h101"))->Fill(tmpintpur);
    } else {
      ((TH1D*)gDirectory->Get("h111"))->Fill(tmpintpur);
    }

    if (fChB) {
      ((TH1D*)gDirectory->Get("mesBp"))->Fill(fMes);
    } else {
      ((TH1D*)gDirectory->Get("mesBz"))->Fill(fMes);
    }
  }
  
  if (fVerbose) cout 
    << "==> b2uNtp::doBreco> "
    << (fChB? "B+" : "B0") << " (seed mode " << fSeedMode << ") "
    << " mes = " << fMes    << " at " << fIndexBestB
    << endl;

  // -- And now the cleanup: Purity cuts and wrong B-D flavor correlation 
  // --------------------------------------------------------------------
  if ((fSeedMode == 0) && (fPurity < IPURDC)) {
    fIndexBestB = -1; 
    ((TH1D*)gDirectory->Get("mesBzFailed"))->Fill(fMes);
    if (fVerbose) cout << "xx> b2uNtp::doBreco> Failed D+ purity cut " << fPurity << "  " << fIntPurity << endl;
    return;
  }
  if ((fSeedMode == 1) && (fPurity < IPURDSTAR)) {
    fIndexBestB = -1; 
    ((TH1D*)gDirectory->Get("mesBzFailed"))->Fill(fMes);
    if (fVerbose) cout << "xx> b2uNtp::doBreco> Failed D*+ purity cut " << fPurity << "  " << fIntPurity << endl;
    return;
  }
  if ((fSeedMode == 2) && (fPurity < IPURD0)) {
    fIndexBestB = -1; 
    ((TH1D*)gDirectory->Get("mesBpFailed"))->Fill(fMes);
    if (fVerbose) cout << "xx> b2uNtp::doBreco> Failed D0 purity cut " << fPurity << "  " << fIntPurity << endl;
    return;
  }
  if ((fSeedMode == 3) && (fPurity < IPURDSTAR0)) {
    fIndexBestB = -1; 
    ((TH1D*)gDirectory->Get("mesBpFailed"))->Fill(fMes);
    if (fVerbose) cout << "xx> b2uNtp::doBreco> Failed D*0 purity cut " << fPurity << "  " << fIntPurity << endl;
    return;
  }

  // --jump the events (reco Bch) with wrong B-D flavor correlation)
  //  int recoilNtp::skipBadBreco() {
  // -- FIXME: THIS IS NOT FLAGGING ANY BRECO CHANNELS AT THE MOMENT. INVESTIGATE! 


  // -- Setup overlap masks with BRECO for recoil analysis
  fBrecoOverlap = 1; 
  if (fIndexBestB == 1) fBrecoOverlap = 2;

  int mask(0); 
  for (int i = 0; i < nTrk; ++i) {
    mask = (fChB? chBRecTrk[i]: B0RecTrk[i]); 
    if (mask & fBrecoOverlap) {
      fBrecoTrk[i] = 1; 
    } else {
      fBrecoTrk[i] = 0; 
    }
  }

  for (int i = 0; i < nGam; ++i) {
    mask = (fChB? chBRecGam[i]: B0RecGam[i]); 
    if (mask & fBrecoOverlap) {
      fBrecoPhoton[i] = 1; 
    } else {
      fBrecoPhoton[i] = 0; 
    }
  }


  // -- Setup Upsilon(4S) and boost vectors 
  double tmpMassPB, tmpMassThetaB, tmpMassPhiB, 
    tmpPB, tmpThetaB, tmpPhiB, 
    tmpMB, tmpBevM;
  if (!fChB) {
    tmpMassPB =MassPB0[fIndexBestB];
    tmpMassThetaB =MassThetaB0[fIndexBestB];
    tmpMassPhiB =MassPhiB0[fIndexBestB];
    tmpPB = pB0[fIndexBestB];
    tmpBevM = massB0[fIndexBestB];
    tmpThetaB =thB0[fIndexBestB];
    tmpPhiB =phiB0[fIndexBestB];
    tmpMB = BZMASS;
  } else {
    tmpMassPB =MassPChB[fIndexBestB];
    tmpMassThetaB =MassThetaChB[fIndexBestB];
    tmpMassPhiB =MassPhiChB[fIndexBestB];
    tmpPB = pChB[fIndexBestB];
    tmpBevM = massChB[fIndexBestB];
    tmpThetaB =thChB[fIndexBestB];
    tmpPhiB =phiChB[fIndexBestB];
    tmpMB = BPMASS;
  }
  mk4Vector(f4Breco, tmpMassPB , tmpMassThetaB, tmpMassPhiB, tmpMB); 
  mk4Vector(f4BrecoNC, tmpPB , tmpThetaB, tmpPhiB, tmpBevM); 
  
  f4Upsilon  = TLorentzVector(pxUps, pyUps, pzUps, eUps); 
  f3UpsBoost = TVector3(pxUps, pyUps, pzUps);
  f3UpsBoost.SetMag(f3UpsBoost.Mag()/eUps);

  f4Brecoil  = f4Upsilon - f4Breco; 
  f3CmsBoost = f4Brecoil.BoostVector();

  TLorentzVector p4t = f4Breco; 
  p4t.Boost(-f3UpsBoost);
  fPcmsBreco = p4t.Vect().Mag(); 

  if(fHistFile!=0){
    fHistFile->cd("doBreco");
    ((TH1D*)gDirectory->Get("meszero"))->Fill(mseB0[0]);
  }

}

// ----------------------------------------------------------------------
// -- Partial reconstruction of B0 -> D*+ (->pi^+ X) \ell
//    (this is the original)
void b2uNtp::doPartialB0() {

  fMM2pr = -11111.;
  if (fBrecoCharge != 0) return;

  TLorentzVector p4t(0., 0., 0., 0.); 
  TLorentzVector p4Dstar(0., 0., 0., 0.);
  double pPiMin(10000.); 

  double DstarMass(2.010);
  double EstarPi(0.145);
  double pionP(0.), pionE(0.), oa(0.);


  for (int i = 0; i < nTrk; ++i) {
    if (fBrecoTrk[i])         continue; 
    if (isRecKaon(i)  && (fGoodKaonTrk[i] == 1)) continue;
    if (isRecEl(i)    && (fGoodElectronTrk[i] == 1)) continue; 
    if (isRecMu(i)    && (fGoodMuonTrk[i] == 1)) continue;    

    if (fLeptonCharge*chargeTrk[i] == -1) {
      mk4Vector(p4t, momentumTrk[i], thetaTrk[i], phiTrk[i], PIPMASS);  

      pionE = p4t.E();

      p4t.Boost(-f3CmsBoost);
      oa = p4t.Angle(f4LeptonCMS.Vect());
      pionP = p4t.P();

      if ( (pionP < 0.2) && (pionP < pPiMin) && (pionE > EstarPi)) {
	double EDstar = pionE*DstarMass/EstarPi;
	mk4Vector(p4Dstar, sqrt(EDstar*EDstar-DstarMass*DstarMass), thetaTrk[i], phiTrk[i], DstarMass); 
	TLorentzVector neutPR = f4Upsilon - f4Breco - f4Lepton - p4Dstar;
	fMM2pr = neutPR.Mag2();
	pPiMin = pionP;
	fOA2   = oa; 
      } 
    }    
  }
}

// ----------------------------------------------------------------------
// -- Partial reconstruction of B+ -> D*0 (->pi^0 X) \ell
//    (the decay D*0 -> D0 gamma may have too high Q-value)
void b2uNtp::doPartialBp() {

  fMM1pr = -11111.;
  if (fBrecoCharge == 0) return;

  TLorentzVector p4t(0., 0., 0., 0.); 
  TLorentzVector p4Dstar(0., 0., 0., 0.);
  double oaMax(0.);

  double DstarMass(2.007);
  double EstarPi(0.142);
  double pionP(0.), pionE(0.), oa(0.), pPiMin(10000.);

  //   for (int i = 0; i < nMc; ++i) {
  //     if ((idMc[i] == 111) && (TMath::Abs(idMc[mothMc[i]-1]) == 423)) {
  //       cout << "Truth pi0 with p = " << pMc[i] << endl;
  //     }
  //   }

  for (int i = 0; i < nPi0; ++i) {
    // -- skip pi0 whos daughters are part of the BRECO
    if (fBrecoPhoton[d1Pi0Index[i]-1]) {
      //      cout << "This photon of a pi0 is part of the BRECO" << endl;
      continue; 
    }
    if (fBrecoPhoton[d2Pi0Index[i]-1]) {
      //      cout << "This photon of a pi0 is part of the BRECO" << endl;
      continue; 
    }

    mk4Vector(p4t, pPi0[i], thPi0[i], phiPi0[i], PIZMASS);  

    pionE = p4t.E();
    
    p4t.Boost(-f3CmsBoost);
    oa = p4t.Angle(f4LeptonCMS.Vect());
    pionP = p4t.P();

    if ((fRecoilGamMult>4) && (oa > 2.2) && (pionP < 0.2) && (pionP < pPiMin) && (pionE > EstarPi)) {
      double EDstar = pionE*DstarMass/EstarPi;
      mk4Vector(p4Dstar, sqrt(EDstar*EDstar-DstarMass*DstarMass), thPi0[i], phiPi0[i], DstarMass); 
      TLorentzVector neutPR = f4Upsilon - f4Breco - f4Lepton - p4Dstar;
      fMM1pr = neutPR.Mag2();
      oaMax  = oa;
      fOA1   = oa; 
      pPiMin = pionP;
      //      cout << i << " Bp: " << pionE << endl;
    } 
  }    
}


// ----------------------------------------------------------------------
// -- Partial reconstruction of B0 -> D*+ (->pi^0 X) \ell
//    (to recover the other 30% of D*+ decays)
void b2uNtp::doPartialB0bis() {

  fMM3pr = -11111.;
  if (fBrecoCharge != 0) return;

  TLorentzVector p4t(0., 0., 0., 0.); 
  TLorentzVector p4Dstar(0., 0., 0., 0.);
  double oaMax(0.);

  double DstarMass(2.010);
  double EstarPi(0.140);
  double pionP(0.), pionE(0.), oa(0.), pPiMin(10000.);

  for (int i = 0; i < nPi0; ++i) {
    // -- skip pi0 whos daughters are part of the BRECO
    if (fBrecoPhoton[d1Pi0Index[i]-1]) {
      //      cout << "This photon of a pi0 is part of the BRECO" << endl;
      continue; 
    }
    if (fBrecoPhoton[d2Pi0Index[i]-1]) {
      //      cout << "This photon of a pi0 is part of the BRECO" << endl;
      continue; 
    }

    mk4Vector(p4t, pPi0[i], thPi0[i], phiPi0[i], PIZMASS);  

    pionE = p4t.E();
    
    p4t.Boost(-f3CmsBoost);
    oa = p4t.Angle(f4LeptonCMS.Vect());
    pionP = p4t.P();

    if ((fRecoilGamMult>4) && (oa > 2.2) && (pionP < 0.2) && (pionP < pPiMin) && (pionE > EstarPi)) {
      double EDstar = pionE*DstarMass/EstarPi;
      mk4Vector(p4Dstar, sqrt(EDstar*EDstar-DstarMass*DstarMass), thPi0[i], phiPi0[i], DstarMass); 
      TLorentzVector neutPR = f4Upsilon - f4Breco - f4Lepton - p4Dstar;
      fMM3pr = neutPR.Mag2();
      oaMax  = oa;
      fOA3   = oa; 
      pPiMin = pionP;
      //      cout << "B0bis: " << pionE << endl;
    } 
  }    
}


// ----------------------------------------------------------------------
// -- This is not very effective as missing particles will lead to a large 
//    delta(M). It may be better to do also a partial reco with D*+ -> D+ pi0
void b2uNtp::doDeltaM() {

  fDeltaM = 11111.;
  TLorentzVector p4t(0., 0., 0., 0.); 
  TLorentzVector p4XminPi(0., 0., 0., 0.);
  double deltaM(0.);

  int mcmom(-1), mcgmom(-1), index(-1);

  //  cout << "new event " << endl;

  for (int i = 0; i < nTrk; ++i) {
    if (fBrecoTrk[i])         continue; 
    if (isRecKaon(i)  && (fGoodKaonTrk[i] == 1)) continue;
    if (isRecEl(i)    && (fGoodElectronTrk[i] == 1)) continue; 
    if (isRecMu(i)    && (fGoodMuonTrk[i] == 1)) continue;    

    if (fLeptonCharge*chargeTrk[i] == -1) {
      mk4Vector(p4t, momentumTrk[i], thetaTrk[i], phiTrk[i], PIPMASS);  
      p4XminPi = f4Xhad - p4t;
      deltaM = fMxhad - p4XminPi.Mag();
      
      if (fDeltaM > deltaM) {
	fDeltaM = deltaM;
	// 	int mom = mothMc[IndexTrk[i]-1];
	// 	mcmom   = idMc[mom]; 
	// 	mcgmom  = idMc[mothMc[mom]-1];
	index = i; 
      }   

      //       // -- Check for lost slow pions
      //       int mom = mothMc[IndexTrk[i]-1];
      //       if (TMath::Abs(idMc[mothMc[mom]-1]) == 413) {
      // 	cout << i << " Track from D*+ with ID = " << idTrk[i] << " and p = " 
      //             << momentumTrk[i] << " deltaM = " << deltaM << endl;
      //       }

    }
  }

  //   if (fDeltaM < 0.15) {
  //     cout << index << " slow pion from mother: " << mcmom << " and grandmother: " << mcgmom << endl;
  //   }
  
}



// ----------------------------------------------------------------------
void b2uNtp::fillPidMaps() {
  int i(0);
  static int first(1), firstK(1), firstE(1), firstM(1);

  // -- Protection against running PidKilling on Data
  if ((runNumber < 100000) && (fOptPidKillingKa || fOptPidKillingEl || fOptPidKillingMu)) {
    if (first) cout << "==> b2uNtp::FillPidMaps> Resetting PidKilling options, since this is data" << endl;
    fOptPidKillingKa = 0;
    fOptPidKillingEl = 0;
    fOptPidKillingMu = 0;
  }

  for (i = 0; i < nTrk; ++i) {
    // -- Kaons
    if (TMath::Abs(kaonIdTrk[i]) & IDKA) {
      fRecKa[i] = 1; 
    } else {
      fRecKa[i] = 0; 
    }

    if (fOptPidKillingKa > 0){
      if (firstK) { 
	firstK = 0; 
	cout << "==> b2uNtp::FillPidMaps> Running with Kaon pidkilling" << endl;    
      }
      if (isPidKillKaon(i)) {
	fRecKa[i] = 1; 
      } else {
	fRecKa[i] = 0;
      }
    } 
    if (fGoodKaonTrk[i] == 0) fRecKa[i] = 0; 

    // -- Electrons
    if (TMath::Abs(elecIdTrk[i]) & IDEL) {
      fRecEl[i] = 1; 
    } else {
      fRecEl[i] = 0;
    }

    if (fOptPidKillingEl > 0) {
      if (firstE) { 
	firstE = kFALSE; 
	cout << "==> b2uNtp::FillPidMaps> Running with Electron pidkilling" << endl;        
      }
      if (isPidKillEl(i)) {
	fRecEl[i] = 1;
      } else {
	fRecEl[i] = 0;
      }
    } 
    if (fGoodElectronTrk[i] == 0) fRecEl[i] = 0; 
    
    // -- Muons
    if ((TMath::Abs(muonIdTrk[i]) & IDMU) && ((TMath::Abs(kaonIdTrk[i]) & IDKA) == 0))  {
      fRecMu[i] = 1; 
    } else {
      fRecMu[i] = 0;
    }

    if (fOptPidKillingMu > 0) {
      if (firstM) { 
	firstM = kFALSE; 
	cout << "==> b2uNtp::FillPidMaps> Running with Muon pidkilling" << endl;    
      }
      if (isPidKillMu(i)) {
	fRecMu[i] = 1;
      } else {
	fRecMu[i] = 0;
      }
    } 
    if (fGoodMuonTrk[i] == 0) fRecMu[i] = 0; 
    
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

// ----------------------------------------------------------------------
Bool_t b2uNtp::isPidKillEl(int i) {
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
Bool_t b2uNtp::isPidKillMu(int i) {
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
Bool_t b2uNtp::isPidKillKaon(int i) {
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
// 0 acceptance and energy cut
// 1 fklm: Fortin, Kowalewski, Lacker, and Muheim
void b2uNtp::selectPhotons() {

  static int killedNeuts(0);
  static int first(1);
  if (first) {
    first = 0;
    cout << "==> b2uNtp::selectPhotons> Selecting photons at level " << PHOTONSELECTION << endl;
    if (DONEUTKILLING &1) {
      cout << "==> b2uNtp::selectPhotons> Killing neutrals with flat probability " << NEUTKILL << endl;
    } else {
      cout << "==> b2uNtp::selectPhotons> No neutral killing " << endl;
    }
  }

  Bool_t gpl(kTRUE), gpd(kFALSE), acc(kFALSE), superric(kFALSE), klsel(kFALSE), ric(kFALSE), ricTight(kFALSE), gg(kFALSE);
  Bool_t fklm(kFALSE); // fortin, kowalewski, lacker, muheim
  Bool_t rsd(kFALSE);  // Rolf Dubitzky

  static double magInEndCap = sqrt( 91.*91. + 180.*180. );
  double docaMin  = 1000;
  double alphaMin = 1000;
  TVector3 p3Gam( 0., 0., 1000.);

  TLorentzVector p4Gam(0., 0., 0., 0.);
  int it(0); 
  Double_t rand(0.);

  TLorentzVector p4Trk(0., 0., 0., 0.);

  for (int i = 0; i < nGam; ++i) {
    fGoodPhoton[i] = 0;
    gpl = gpd = acc = ric = kFALSE;
    acc = ((DR*thetaGam[i] > GAMMATHETALO) && (DR*thetaGam[i] < GAMMATHETAHI));
    gpl = ((energyGam[i] >= 0.030)
           && (nCryGam[i] >= 1.)
           && (lMomGam[i] <= 0.8));
    gpd = ((energyGam[i] >= 0.100)
           && (nCryGam[i] >= 1.)
           && (lMomGam[i] <= 0.8));


    mk4Vector(p4Gam, energyGam[i], thetaGam[i], phiGam[i], 0.);

    // -- Determine minimum opening angle to all non-electrons (see BAD 633)
    double alpha(0.), alphaMin(99.); 
    if (acc && (energyGam[i] > 0))  {

      /* RD doca calculation  */
      if ( thetaGam[i] > 0.469 )
	{
	  /* in barrel */
	  // p3Gam.SetPtThetaPhi( 91.0, thetaGam[i], phiGam[i] );
	  double pt = 91.0;
	  double theta = thetaGam[i];
	  double phi =  phiGam[i];
	  double xX = pt * TMath::Cos(phi);
	  double xY = pt * TMath::Sin(phi);
	  Double_t tanTheta = TMath::Tan(theta);
	  double xZ = tanTheta ? pt / tanTheta : 0;
	  p3Gam.SetXYZ( xX, xY, xZ );
	  
	} else {
	  /* in endcap */
	  // p3Gam.SetMapThetaPhi( magInEndCap, thetaGam[i], phiGam[i] );
	  Double_t amag = TMath::Abs( magInEndCap );
	  double theta = thetaGam[i];
	  double phi = phiGam[i];
	  double xX = amag * TMath::Sin(theta) * TMath::Cos(phi);
	  double xY = amag * TMath::Sin(theta) * TMath::Sin(phi);
	  double xZ = amag * TMath::Cos(theta);
	  p3Gam.SetXYZ( xX, xY, xZ );
	}


      for (it = 0; it < nTrk; ++it) {
	if (isRecEl(it)) continue; // FIXME DO WE WANT THIS? 
	alpha = TMath::Cos(thetaGam[i])*TMath::Cos(thetaAtEMCTrk[it])
	  +  TMath::Sin(thetaGam[i])*TMath::Sin(thetaAtEMCTrk[it])*TMath::Cos(phiGam[i] - phiAtEMCTrk[it]); 
	alpha = TMath::ACos(alpha); 
	if (alpha < alphaMin) alphaMin = alpha; 	


	/* RD doca calculation  */
	TVector3 p3Trk( 0., 0., 1000. );
	if ( thetaAtEMCTrk[it] > 0.469 )
	  {
	    /* in barrel */
	    // p3Gam.SetPtThetaPhi( 91.0, thetaGam[i], phiGam[i] );
	    double pt = 91.0;
	    double theta = thetaAtEMCTrk[it];
	    double phi =  phiAtEMCTrk[it];
	    double xX = pt * TMath::Cos(phi);
	    double xY = pt * TMath::Sin(phi);
	    Double_t tanTheta = TMath::Tan(theta);
	    double xZ = tanTheta ? pt / tanTheta : 0;
	    p3Trk.SetXYZ( xX, xY, xZ );
	    
	  } else {
	    /* in endcap */
	    // p3Gam.SetMapThetaPhi( magInEndCap, thetaGamAtEMCTrkAtEMCTrk[i], phiGamAtEMCTrk[i] );
	    Double_t amag = TMath::Abs( magInEndCap );
	    double theta = thetaAtEMCTrk[it];
	    double phi = phiAtEMCTrk[it];
	    double xX = amag * TMath::Sin(theta) * TMath::Cos(phi);
	    double xY = amag * TMath::Sin(theta) * TMath::Sin(phi);
	    double xZ = amag * TMath::Cos(theta);
	    p3Trk.SetXYZ( xX, xY, xZ );
	  }
	
	TVector3 dist = p3Trk - p3Gam;
	double doca = dist.Mag();
	if ( doca<docaMin )
	  {
	    docaMin         = doca;
	  }
	if (  p3Trk.Angle( p3Gam ) < alphaMin ) {
	  alphaMin =  p3Trk.Angle( p3Gam );
	}


      }
    }


    fklm = (
	    acc
	    && (ecalGam[i] > 0.050)
	    && (nCryGam[i] > 2)
	    && (lMomGam[i] < 0.6)
	    && (alpha > 0.08)
	    ); 


    rsd = (
	   acc
	   && (ecalGam[i] > GAMMAELO)
	   && (nCryGam[i] > 2)
	   && (lMomGam[i] < 0.6)
	   && (docaMin > 20)
	   );
    
    p4Gam.Boost(-f3CmsBoost); 

    superric = energyGam[i] >= GAMMAELO && p4Gam.T()<2.8 && lMomGam[i]>0.05 && lMomGam[i]<0.5 && s9s25Gam[i]>0.9;
    //2bu    klsel = energyGam[i] >= GAMMAELO && p4Gam.T()<2.8 && KLlikeEMC(i)<0;
        
    ric = energyGam[i] >= GAMMAELO && energyGam[i]<4. &&  lMomGam[i]>0.05;
    ricTight = energyGam[i] >= GAMMAELO&& energyGam[i]<4. && lMomGam[i]>0.05 && lMomGam[i]<0.5 && s9s25Gam[i]>0.9;

    gg = ((energyGam[i] >= 0.030)
          && (nCryGam[i] >= 1.)
          && (lMomGam[i] <= 0.8)
          && (p4Gam.T() < 2.8)
	  //b2u          && (splitOffGam[i] == 0)
          );
    
    if (PHOTONSELECTION == 0) { 
      if (acc && (energyGam[i] > GAMMAELO))  fGoodPhoton[i] = 1; 
    } 

    if (PHOTONSELECTION == 1) { 
      if (fklm) fGoodPhoton[i] = 1; 
    } 

    if (PHOTONSELECTION == 2) { 
      if (rsd) {
	//	cout << "taken " << docaMin << endl;
	fGoodPhoton[i] = 1; 
      } else {
	//	cout << "not taken "  << docaMin << endl;
      }
    } 

    if (PHOTONSELECTION == 3) { 
      if (superric && acc &&  splitOffGam[i]==0 ) fGoodPhoton[i] = 1;
    } 

    if (PHOTONSELECTION == 4) { 
      if (rsd && (p4Gam.T() < 2.8)) fGoodPhoton[i] = 1;
    } 

    if (PHOTONSELECTION == 5) { 
      if (rsd && (splitOffGam[i]==0)) fGoodPhoton[i] = 1;
    } 

    if (PHOTONSELECTION == 6) { 
      if (rsd && (s9s25Gam[i] > 0.9)) fGoodPhoton[i] = 1;
    } 

    if (PHOTONSELECTION == 7) { 
      if (rsd && (lMomGam[i] < 0.5)) fGoodPhoton[i] = 1;
    } 

    if (PHOTONSELECTION == 8) { 
      if (
	  acc
	  && (ecalGam[i] > GAMMAELO)
	  && (nCryGam[i] > 2)
	  && (lMomGam[i] < 0.6)
	  && (splitOffGam[i]==0)
	  )  fGoodPhoton[i] = 1;
    }

    if (PHOTONSELECTION == 9) { 
      if (
	  acc
	  && (ecalGam[i] > GAMMAELO)
	  && (nCryGam[i] > 2)
	  && (lMomGam[i] < 0.5)
	  && (splitOffGam[i]==0)
	  )  fGoodPhoton[i] = 1;
    }

    if (PHOTONSELECTION == 10) { 
      if (
	  acc
	  && (ecalGam[i] > GAMMAELO)
	  && (nCryGam[i] > 2)
	  && (lMomGam[i] < 0.5)
	  && (s9s25Gam[i] > 0.9)
	  && (splitOffGam[i]==0)
	  )  fGoodPhoton[i] = 1;
    }


    // -- Neutral killing in MC: 
    //    DONEUTKILLING&1: All neutrals with flat probability 1.8%
    if (runNumber > 100000) {
      if (DONEUTKILLING & 1) {
	rand = gRandom->Rndm(); 
	if (rand < NEUTKILL) {
	  ++killedNeuts; 
	  fGoodPhoton[i] = 0; 

	  if (0 == killedNeuts%1000) cout << "killed " << killedNeuts << " neutrals so far (in the entire event)" << endl;
	}      
      }
    }




  }


}

// ----------------------------------------------------------------------
// -- stripped down very much. Must read Jon/Bill's stuff and implement!?
void b2uNtp::selectTracks() {

  Bool_t ct(kTRUE), gtvl(kFALSE), gtl(kFALSE), 
    electron(kFALSE), muon(kFALSE), kaon(kFALSE), 
    acc(kFALSE), noSvtGhost(kFALSE);
  Bool_t sel1(kFALSE);
  Double_t pt(0.), dca(0.), dcaz(0.), rand(0.);

  static int killedTracks(0);
  static int first(1); 
  TH1D *h, *hall, *hgl; 
  fHistFile->cd();
  if (first) {
    first = 0; 
    cout << "==> b2uNtp::selectTracks> Selecting tracks  at level " << TRACKSELECTION << endl;
    if ((DOTRACKKILLING > 0) && (runNumber > 100000)) {
      if (DOTRACKKILLING &1) cout << "==> b2uNtp::selectTracks> Killing tracks with flat probability " << TRACKKILL << endl;
    } else {
      cout << "==> b2uNtp::selectTracks> No track killing " << endl;
    }
    h = new TH1D("gltheta", "theta for good low pT tracks", 100, 0., 3.15); 
    h = new TH1D("gtheta", "theta for good low pT tracks", 100, 0., 3.15); 
    h = new TH1D("atheta", "theta for all low pT tracks", 100, 0., 3.15); 
  }

  h = (TH1D*)gDirectory->Get("gtheta");    
  hall = (TH1D*)gDirectory->Get("atheta");    
  hgl = (TH1D*)gDirectory->Get("gltheta");    

  TLorentzVector p4Trk(0., 0., 0., 0.);
  double pcms(0.), thetadr(0.);
  double  pcmslo(9999.);  int  ind_pcmslo(-99);
  for (int i = 0; i < nTrk; ++i) {
    ct = gtvl = gtl = acc = 0;
    electron = muon = kaon = 0; 

    pt      = momentumTrk[i]*sin(thetaTrk[i]);
    dcaz    = zPocaTrk[i] - beamSZ;
    dca     = TMath::Sqrt((xPocaTrk[i]-beamSX)*(xPocaTrk[i]-beamSX) + (yPocaTrk[i]-beamSY)*(yPocaTrk[i]-beamSY));
    thetadr = DR*thetaTrk[i]; 

    noSvtGhost = kTRUE; 
    if (pt > 0.2) {
      if (ndchTrk[i] > 0) {
	noSvtGhost = kTRUE;   // this is the good case
      } else {
	noSvtGhost = kFALSE;  // and this the bad 
      }
    }

    ct   = kTRUE;
    acc  = ((thetadr > PITHETALO) && (thetadr < PITHETAHI));
    gtvl = ((pt > 0.0) 
            && (momentumTrk[i] < 10.0)
            && (tproTrk[i] >= 0.) 
            && (dca  <= 1.5)
            && (TMath::Abs(dcaz) <= 10.0));
    gtl  = ((pt > 0.1) 
            && (momentumTrk[i] <= 10.0)
            && (ndchTrk[i] >= 12)
            && (tproTrk[i] >= 0.) 
            && (dca  <= 1.5)
            && (TMath::Abs(dcaz) <= 10.0));

    mk4Vector(p4Trk, momentumTrk[i], thetaTrk[i], phiTrk[i], ELMASS);
    
    p4Trk.Boost(-f3CmsBoost); 
    pcms = p4Trk.Rho();

    sel1 =  (pt > 0.06) 
      && (TMath::Abs(dcaz) <= 5.0)
      && (pcms < 2.7)
      ; 



    electron = ((momentumTrk[i] > 0.5) && (pcms > ELMOMLO) && (thetadr > ELTHETALO) && (thetadr < ELTHETAHI));
    muon     = ((momentumTrk[i] > 0.5) && (pcms > MUMOMLO) && (thetadr > MUTHETALO) && (thetadr < MUTHETAHI));
    kaon     = ((momentumTrk[i] > KAMOMLO) && (thetadr > KATHETALO) && (thetadr < KATHETAHI));

    fGoodElectronTrk[i] = fGoodMuonTrk[i] = fGoodPionTrk[i] = fGoodKaonTrk[i] = 0; 

    if (pt < 0.2) hall->Fill(thetaTrk[i]); 


    if (TRACKSELECTION == 0) { 
      if (gtvl 
	  && noSvtGhost
	  && acc
	  ) {

	if (pcms<pcmslo) {
          pcmslo = pcms;
          ind_pcmslo = i;
        }

	fGoodPionTrk[i] = 1; 
	if (electron && gtl) fGoodElectronTrk[i] = 1; 
	if (muon && gtl) fGoodMuonTrk[i] = 1; 
	if (kaon) fGoodKaonTrk[i] = 1; 
      }
    } 

    if (TRACKSELECTION == 1) { 
      if (gtvl 
	  && noSvtGhost
	  && acc
	  && (fCleanGoodTrack[i] == 1)
	  ) {
	
        if (pcms<pcmslo) {
          pcmslo = pcms;
          ind_pcmslo = i;
        }

	fGoodPionTrk[i] = 1; 
	if (electron && gtl) fGoodElectronTrk[i] = 1; 
	if (muon && gtl) fGoodMuonTrk[i] = 1; 
	if (kaon) fGoodKaonTrk[i] = 1; 

	if (pt < 0.2) h->Fill(thetaTrk[i]); 

      }

      if (gtvl 
	  && noSvtGhost
	  && acc) {
	if (pt < 0.2) hgl->Fill(thetaTrk[i]); 
      }

    }


    // -- Track killing in MC: 
    //    DOTRACKKILLING&1: All tracks with flat probability 1.3%
    if (runNumber > 100000) {
      if (DOTRACKKILLING & 1) {
	rand = gRandom->Rndm(); 
	if (rand < TRACKKILL) {
	  ++killedTracks; 
	  fGoodPionTrk[i] = 0; 
	  fGoodElectronTrk[i] = 0; 
	  fGoodMuonTrk[i] = 0; 
	  fGoodKaonTrk[i] = 0; 

	  if (0 == killedTracks%1000) cout << "killed " << killedTracks << " tracks so far (in the entire event)" << endl;
	}      
      }
    }


 
  }

  fpcmsTrklo = pcmslo;

}



// ----------------------------------------------------------------------
void b2uNtp::cleanupTracks() {

  const double phiSSmatch(0.200); 
  const double phiOSmatch(0.200); 
  const double thetaSSmatch(0.200);
  const double thetaOSmatch(0.200); 
  const double ptmatch(0.1); 

  const double paraphimatch(0.100); 
  const double parathetamatch(0.100); 
  const double paraptmatch(0.050); 

  int kill(0); 
  double tlo(1.4), thi(1.7); 
  double pi(3.1415); 

  for (int i = 0; i < 100; ++i) {
    fCleanGoodTrack[i] = 1; 
  }

  double pt, dz; 
  double ptj, dzj; 
  double dphi, dtheta; 

  for (int i = 0; i < nTrk; ++i) {

    pt = momentumTrk[i]*TMath::Sin(thetaTrk[i]); 
    dz = TMath::Abs(zPocaTrk[i] - primVtxZ);

    for (int j = i+1; j < nTrk; ++j) {
      
      ptj    = momentumTrk[j]*TMath::Sin(thetaTrk[j]); 
      dzj    = TMath::Abs(zPocaTrk[j] - primVtxZ);

      // -- opposite sign loopers at opposite phi (central theta)
      dphi   = unwrapped(3.1415 + phiTrk[i], phiTrk[j]); 
      dtheta = TMath::Abs(3.1415 - thetaTrk[i] - thetaTrk[j]); 
      if ((chargeTrk[i]*chargeTrk[j] < 0.)
	  && (thetaTrk[i] > tlo) && (thetaTrk[i] < thi) 
	  && (thetaTrk[j] > tlo) && (thetaTrk[j] < thi) 
	  && (pt  < 0.18)
	  && (ptj < 0.18)
	  && (dphi < phiOSmatch)
	  && (dtheta < thetaOSmatch)
	  ) {
          
	if (dzj < dz) {
	  fCleanGoodTrack[i] = 0; 
	  //	  kill = i; 
	} else {
	  fCleanGoodTrack[j] = 0; 
	  //	  kill = j; 
	}
      }

      // -- same sign loopers at same phi (central theta)
      dtheta   = TMath::Abs(thetaTrk[i] - thetaTrk[j]);
      dphi     = unwrapped(phiTrk[i], phiTrk[j]); 
      if ((chargeTrk[i]*chargeTrk[j] > 0.)
	  && (thetaTrk[i] > tlo) && (thetaTrk[i] < thi) 
	  && (thetaTrk[j] > tlo) && (thetaTrk[j] < thi) 
	  && (pt  < 0.18)
	  && (ptj < 0.18)
	  && (dphi < phiSSmatch)
	  && (dtheta < thetaSSmatch)
	  ) {
	
	if (dzj < dz) {
	  fCleanGoodTrack[i] = 0; 
	  //	  kill = i; 
	} else {
	  fCleanGoodTrack[j] = 0; 
	  //	  kill = j; 
	}
      }

      
      // -- ghost tracks "sharing" DCH hits
      dtheta   = TMath::Abs(thetaTrk[i] - thetaTrk[j]);
      dphi     = unwrapped(phiTrk[i], phiTrk[j]); 
      if ((chargeTrk[i]*chargeTrk[j] > 0.)
	  && (pt > 0.)   && (pt < 0.35)   
	  && (dtheta < parathetamatch)
	  && (dphi < paraphimatch)
	  && (TMath::Abs(pt - ptj) < paraptmatch)
	  ) {
	if (ndchTrk[i] < ndchTrk[j]) {
	  fCleanGoodTrack[i] = 0; 
	  //	  kill = i; 
	} else {
	  fCleanGoodTrack[j] = 0; 
	  kill = j; 
	}
	
      }
      
    }

  }

}


// ----------------------------------------------------------------------
void b2uNtp::selectKlongs() {

  static int first(1); 
  if (first) {
    first = 0; 
    cout << "==> b2uNtp::selectKlongs> Selecting Klongs  at level " << KLSELECTION << endl;
  }


  double lh0, lh1, lh2; 
  double lat, ecal; 
  for (int i = 0; i < nGam; ++i) {
    fGoodKl[i] = 0;
    lh0   = lhKlong00(i); 
    //    lh1   = lhKlong01(i); 
    //    lh2   = lhKlong02(i); 
   lat   = lMomGam[i]; 
    ecal  = ecalGam[i]; 

    // -- do not touch fGoodPhoton! 
    if (KLSELECTION == 0) {
      if ((lat > 0.5) && (ecal>0.2) && lh0 > -0.5) {
	fGoodKl[i] = 1; 
      }
    }
  }

}

// ----------------------------------------------------------------------
void b2uNtp::selectKshorts() {

  static int first(1); 
  if (first) {
    first = 0; 
    cout << "==> b2uNtp::selectKshorts> Selecting Kshorts  at level " << KSSELECTION << endl;
  }


  int noel, pi1El, pi2El, noka, pi1Ka, pi2Ka; 
  TVector3 p3K(0., 0., 0.), v3K(0., 0., 0.); 
  double chi2, alpha, mks, pks, rks; 

  //   rks   > 0.325
  //   alpha > -0.82
  //   pks   > 0.15 
  //   chi2  < 10.5

  for (int i = 0; i < nKs; ++i) {

    // -- skip K0S->pi pi
    if (TMath::Abs(d1KsLund[i]) != 211) {
      continue;
    }

    fGoodKs[i] = 0;
   
    int pi1 = d1KsIndex[i]-1;
    int pi2 = d2KsIndex[i]-1; 

    pi1El = (isRecEl(pi1)? 1:0);
    pi2El = (isRecEl(pi2)? 1:0);
    pi1Ka = (isRecKaon(pi1)? 1:0);
    pi2Ka = (isRecKaon(pi2)? 1:0);
    
    noel = pi1El + pi2El; 
    noka = pi1Ka + pi2Ka; 

    chi2 = chi2Ks[i];
    mks = massKs[i];
    pks = pKs[i];
    rks = TMath::Sqrt( (xKs[i]-beamSX)*(xKs[i]-beamSX) + (yKs[i]-beamSY)*(yKs[i]-beamSY) );

    mk3Vector(p3K, pKs[i], thKs[i], phiKs[i]); 
    v3K.SetXYZ(xKs[i], yKs[i], zKs[i]); 
    alpha = p3K.Angle(v3K); 
    
    if (KSSELECTION == 0) {
      if (
	  (noel == 0)  && (noka == 0) 
	  && (0.486 < mks) && (mks < 0.510)
	  ) {
	fGoodKs[i] = 1; 
      }
    }

    if (KSSELECTION == 1) {
      if (
	  (noel == 0)  
	  //	  && (noka == 0) 
	  && (0.490 < mks) && (mks < 0.505)
	  && (rks > 0.325) 
	  && (TMath::Cos(alpha) >-0.82)
	  && (pks > 0.15)
	  && (chi2 < 10.0)
	  ) {
	fGoodKs[i] = 1; 
	// 	cout << "Selected as good KS " << idMc[MCKs[i]-1] 
	// 	     << " with mass " << mks 
	// 	     << " chi2: " << chi2
	// 	     << endl;
      } else {
	// 	if (idMc[MCKs[i]-1] == 310)
	// 	  cout << "Failed as good KS: " << idMc[MCKs[i]-1] 
	// 	       << " with mass: " << mks 
	// 	       << " noel: " << noel 
	// 	       << " noka: " << noka
	// 	       << " pks: " << pks
	// 	       << " rks: " << rks
	// 	       << " chi2: " << chi2
	// 	       << " cos(a): " << TMath::Cos(alpha)
	// 	       << endl;
      }	
    }

  }
}


// ----------------------------------------------------------------------
void b2uNtp::selectKshortsZZ() {

  static int first(1); 
  if (first) {
    first = 0; 
    cout << "==> b2uNtp::selectKshortsZZ> Selecting KshortsZZ  at level " << KSZZSELECTION << endl;
  }

  for (int i = 0; i < nKs; ++i) {

    // -- skip K0S->pi0 pi0
    if (TMath::Abs(d1KsLund[i]) != 111) {
      continue;
    }

    fGoodKs[i] = 0;

  }
}



// ----------------------------------------------------------------------
void b2uNtp::readCuts(TString filename) {

  cout << "======================================="  << endl;
  cout << "==> b2uNtp::readCuts> reading cut file: " << filename.Data() << endl;
  cout << "---------------------------------------"  << endl;

  char  buffer[200];
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
    if (!strcmp(CutName, "purity"))     {
      PURITY = CutValue;     ok = 1;
      cout << "purity:            " << PURITY << endl;
    }
    if (!strcmp(CutName, "intPurity"))  {
      INTPURITY = CutValue;  ok = 1;
      cout << "intPurity:         " << INTPURITY << endl;
    }
    if (!strcmp(CutName, "ipurDstar"))  {
      IPURDSTAR = CutValue;  ok = 1;
      cout << "ipurDstar:         " << IPURDSTAR << endl;
    }
    if (!strcmp(CutName, "ipurDc"))     {
      IPURDC = CutValue;     ok = 1;
      cout << "ipurDc:            " << IPURDC << endl;
    }
    if (!strcmp(CutName, "ipurDstar0")) {
      IPURDSTAR0 = CutValue; ok = 1;
      cout << "ipurDstar0:        " << IPURDSTAR0 << endl;
    }
    if (!strcmp(CutName, "ipurD0"))     {
      IPURD0 = CutValue;     ok = 1;
      cout << "ipurD0:            " << IPURD0 << endl;
    }

    // -- Particle ID: selector level
    if (!strcmp(CutName, "idEl")) {
      IDEL = int(CutValue); ok = 1;
      cout << "idEl:              " << IDEL << endl;
    }
    if (!strcmp(CutName, "idMu")) {
      IDMU = int(CutValue); ok = 1;
      cout << "idMu:              " << IDMU << endl;
    }
    if (!strcmp(CutName, "idKa")) {
      IDKA = int(CutValue); ok = 1;
      cout << "idKa:              " << IDKA << endl;
    }

    if (!strcmp(CutName, "TrackSelection")) {
      TRACKSELECTION = int(CutValue); ok = 1;
      cout << "TrackSelection:    " << TRACKSELECTION << endl;
    }
    if (!strcmp(CutName, "trackKilling"))   {
      DOTRACKKILLING = (int)CutValue; ok = 1;
      cout << "trackKilling:    "   << DOTRACKKILLING << endl;
    }
    if (!strcmp(CutName, "trackKill")   )   {
      TRACKKILL = CutValue; ok = 1;
      cout << "trackKill:    "   << TRACKKILL << endl;
    }

    if (!strcmp(CutName, "PhotonSelection")) {
      PHOTONSELECTION = int(CutValue); ok = 1;
      cout << "PhotonSelection:   " << PHOTONSELECTION << endl;
    }
    if (!strcmp(CutName, "neutKilling"))   {
      DONEUTKILLING = (int)CutValue; ok = 1;
      cout << "neutKilling:    "   << DONEUTKILLING << endl;
    }
    if (!strcmp(CutName, "neutKill")   )   {
      NEUTKILL = CutValue; ok = 1;
      cout << "neutKill:    "   << NEUTKILL << endl;
    }

    if (!strcmp(CutName, "KlSelection")) {
      KLSELECTION = int(CutValue); ok = 1;
      cout << "KlSelection:    " << KLSELECTION << endl;
    }

    if (!strcmp(CutName, "KsSelection")) {
      KSSELECTION = int(CutValue); ok = 1;
      cout << "KsSelection:    " << KSSELECTION << endl;
    }

    if (!strcmp(CutName, "KszzSelection")) {
      KSZZSELECTION = int(CutValue); ok = 1;
      cout << "KszzSelection:    " << KSZZSELECTION << endl;
    }


    // -- Particles: acceptance
    if (!strcmp(CutName, "elMomLo")) {
      ELMOMLO = CutValue; ok = 1;
      cout << "elMomLo:         " << ELMOMLO << endl;
    }
    if (!strcmp(CutName, "elThetaLo")) {
      ELTHETALO = CutValue; ok = 1;
      cout << "elThetaLo:         " << ELTHETALO << endl;
    }
    if (!strcmp(CutName, "elThetaHi")) {
      ELTHETAHI = CutValue; ok = 1;
      cout << "elThetaHi:         " << ELTHETAHI << endl;
    }

    if (!strcmp(CutName, "elBremCone")) {
      ELBREMCONE = CutValue; ok = 1;
      cout << "elBremCone:        " << ELBREMCONE << endl;
    }


    if (!strcmp(CutName, "muMomLo")) {
      MUMOMLO = CutValue; ok = 1;
      cout << "muMomLo:           " << MUMOMLO << endl;
    }
    if (!strcmp(CutName, "muThetaLo")) {
      MUTHETALO = CutValue; ok = 1;
      cout << "muThetaLo:           " << MUTHETALO << endl;
    }
    if (!strcmp(CutName, "muThetaHi")) {
      MUTHETAHI = CutValue; ok = 1;
      cout << "muThetaHi:           " << MUTHETAHI << endl;
    }

    if (!strcmp(CutName, "kaMomLo")) {
      KAMOMLO = CutValue; ok = 1;
      cout << "kaMomLo:           " << KAMOMLO << endl;
    }
    if (!strcmp(CutName, "kaThetaLo")) {
      KATHETALO = CutValue; ok = 1;
      cout << "kaThetaLo:           " << KATHETALO << endl;
    }
    if (!strcmp(CutName, "kaThetaHi")) {
      KATHETAHI = CutValue; ok = 1;
      cout << "kaThetaHi:           " << KATHETAHI << endl;
    }

    if (!strcmp(CutName, "piMomLo")) {
      PIMOMLO = CutValue; ok = 1;
      cout << "piMomLo:           " << PIMOMLO << endl;
    }
    if (!strcmp(CutName, "piThetaLo")) {
      PITHETALO = CutValue; ok = 1;
      cout << "piThetaLo:           " << PITHETALO << endl;
    }
    if (!strcmp(CutName, "piThetaHi")) {
      PITHETAHI = CutValue; ok = 1;
      cout << "piThetaHi:           " << PITHETAHI << endl;
    }

    if (!strcmp(CutName, "gammaELo")) {
      GAMMAELO = CutValue; ok = 1;
      cout << "gammaELo:           " << GAMMAELO << endl;
    }
    if (!strcmp(CutName, "gammaThetaLo")) {
      GAMMATHETALO = CutValue; ok = 1;
      cout << "gammaThetaLo:           " << GAMMATHETALO << endl;
    }
    if (!strcmp(CutName, "gammaThetaHi")) {
      GAMMATHETAHI = CutValue; ok = 1;
      cout << "gammaThetaHi:           " << GAMMATHETAHI << endl;
    }

    // -- Pidmaps 
    if (!strcmp(CutName, "shiftElPidTables")) {
      SHIFTELPIDTABLES = CutValue; ok = 1;
    }
    if (!strcmp(CutName, "shiftMuPidTables")) {
      SHIFTMUPIDTABLES = CutValue; ok = 1;
    }
    if (!strcmp(CutName, "shiftKaPidTables")) {
      SHIFTKAPIDTABLES = CutValue; ok = 1;
    }
    if (!strcmp(CutName, "shiftElMisTables")) {
      SHIFTELMISTABLES = CutValue; ok = 1;
    }
    if (!strcmp(CutName, "shiftMuMisTables")) {
      SHIFTMUMISTABLES = CutValue; ok = 1;
    }
    if (!strcmp(CutName, "shiftKaMisTables")) {
      SHIFTKAMISTABLES = CutValue; ok = 1;
    }

    if (!strcmp(CutName, "PidTables")) {
      char tablefile[1000];
      sscanf(buffer, "%s %s", CutName, tablefile);
      sprintf(PIDTABLES, "%s", tablefile); ok = 1;
      getPidTables();
    } 


    if (ok == 0)  cout << "==> b2uNtp::readCuts() Error: Don't know about variable " << CutName << endl;
  }

  cout << "---------------------------------------"  << endl;
  
}


// ----------------------------------------------------------------------
void b2uNtp::getPidTables() {
  cout << "==> b2uNtp::getPidTables> Reading PidTables from " 
       << PIDTABLES << endl;

  ifstream is(PIDTABLES);
  char tableName[1000], selector[100], buffer[1000], fname[200];
  int el(0), mu(0), ka(0), source(0), sink(0);
  while (is.getline(buffer, 200, '\n')) {
    if (buffer[0] == '#') {continue;}
    sscanf(buffer, "%s %d %s %d", tableName, &source, selector, &sink);
    if (TMath::Abs(sink) == 11) {
      cout << "Electron Table " << el << "  ->  " << tableName << endl;
      sprintf(fname, "%s%d_%d", selector, sink, source); 
      fPTel[el] = new PIDTable(tableName, fname);
      if (TMath::Abs(source) == TMath::Abs(sink)) {
        cout << "shifting el pidtables relative by " << SHIFTELPIDTABLES << endl;
        if (TMath::Abs(SHIFTELPIDTABLES) > 0.001) fPTel[el]->shiftRel(SHIFTELPIDTABLES);
      } else {
        cout << "shifting el misdtables relative by " << SHIFTELMISTABLES << endl;
        if (TMath::Abs(SHIFTELMISTABLES) > 0.001) fPTel[el]->shiftRel(SHIFTELMISTABLES);
      }
      el++;
    }
    if (TMath::Abs(sink) == 13) {
      cout << "Muon Table " << mu << "  ->  " << tableName << endl;
      sprintf(fname, "%s%d_%d", selector, sink, source); 
      fPTmu[mu] = new PIDTable(tableName, fname);
      if (TMath::Abs(source) == TMath::Abs(sink)) {
        cout << "shifting mu pidtables relative by " << SHIFTMUPIDTABLES << endl;
        if (TMath::Abs(SHIFTMUPIDTABLES) > 0.001) fPTmu[mu]->shiftRel(SHIFTMUPIDTABLES);
      } else {
        cout << "shifting mu mistables relative by " << SHIFTMUMISTABLES << endl;
        if (TMath::Abs(SHIFTMUMISTABLES) > 0.001) fPTmu[mu]->shiftRel(SHIFTMUMISTABLES);
      }
      mu++;
    }
    if (TMath::Abs(sink) == 321) {
      cout << "Kaon Table " << ka << "  ->  " << tableName << endl;
      sprintf(fname, "%s%d_%d", selector, sink, source); 
      fPTka[ka] = new PIDTable(tableName, fname);
      if (TMath::Abs(source) == TMath::Abs(sink)) {
        cout << "shifting ka pidtables relative by " << SHIFTKAPIDTABLES << endl;
        if (TMath::Abs(SHIFTKAPIDTABLES) > 0.001) fPTka[ka]->shiftRel(SHIFTKAPIDTABLES);
      } else {
        cout << "shifting ka mistables relative by " << SHIFTKAMISTABLES << endl;
        if (TMath::Abs(SHIFTKAMISTABLES) > 0.001) fPTka[ka]->shiftRel(SHIFTKAMISTABLES);
      }
      ka++;
    }
  }
}    


// ----------------------------------------------------------------------
Bool_t b2uNtp::isAncestor(int ancestor, int cand) {
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
Bool_t b2uNtp::isStable(int id) {
  int aid = TMath::Abs(id);
  if ((aid == 22) || (aid == 11) || (aid == 13) 
      || (aid == 211) 
      || (aid == 321) || (aid == 130)
      || (aid == 2212) || (aid == 2112)
      ) {
    return true;
  }
  return false;
}


// ----------------------------------------------------------------------
void b2uNtp::getMcDaughters(int imc, int &d1, int &d2, int &d3, int &d4, int &d5) {
  int cnt(1);
  for (int i = 0; i <nMc; ++i) {
    if (imc == (mothMc[i]-1)) {
      if (cnt == 1) {
	d1 = i; ++cnt; 
      } else if (cnt == 2) {
	d2 = i; ++cnt; 
      } else if (cnt == 3) {
	d3 = i; ++cnt; 
      } else if (cnt == 4) {
	d4 = i; ++cnt; 
      } else if (cnt == 5) {
	d5 = i; ++cnt; 
      }
    }
  }
}


// ----------------------------------------------------------------------
// prints overlap information if indices to the two B are passed
void b2uNtp::dumpGeneratorBlock(int b1, int b2) {
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
int b2uNtp::isRecoed(int imc) {
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
void b2uNtp::closeHistFile() {
  cout << "==> b2uNtp::closeHistFile> Writing " << fHistFile->GetName() << endl;
  fHistFile->cd();
  fHistFile->Write();
  fHistFile->Close();
  delete fHistFile;
}



// ----------------------------------------------------------------------
double b2uNtp::kPlus() {
  double mB(5.279), mb(4.800);

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
    return -99.;
  }

  double ds = TMath::Sqrt(TMath::Power(xMc[xu]-xMc[xum],2) + TMath::Power(yMc[xu]-yMc[xum],2) + TMath::Power(zMc[xu]-zMc[xum],2));
  double ct = 10.*ds*massMc[xu]/pMc[xu];
  double qp = 10000.*ct;
  double kp = mB - mb - qp;

  return kp;
}





// ======================================================================
// -- Anything below is candidate for the attic.icc file
// ======================================================================

// ----------------------------------------------------------------------
// -- The rest is provided by ROOT (more or less). No need to edit. 
// ----------------------------------------------------------------------
b2uNtp::b2uNtp(TTree *tree, int isMC) {
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    cout << "xx> b2uNtp::b2uNtp> Error, no tree given" << endl;
  }

  fIsMC = isMC; 
  Init(tree);
}

// ----------------------------------------------------------------------
b2uNtp::~b2uNtp() {
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

// ----------------------------------------------------------------------
Int_t b2uNtp::GetEntry(Int_t entry) {
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

// ----------------------------------------------------------------------
Int_t b2uNtp::LoadTree(Int_t entry) {
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
void b2uNtp::Init(TTree *tree) {
  //   Set branch addresses
  if (tree == 0) return;
  fChain    = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
  
  fChain->SetBranchAddress("event",&event);
  fChain->SetBranchAddress("runNumber",&runNumber);
  fChain->SetBranchAddress("platform",&platform);
  fChain->SetBranchAddress("partition",&partition);
  fChain->SetBranchAddress("upperID",&upperID);
  fChain->SetBranchAddress("lowerID",&lowerID);
  fChain->SetBranchAddress("primVtxX",&primVtxX);
  fChain->SetBranchAddress("primVtxY",&primVtxY);
  fChain->SetBranchAddress("primVtxZ",&primVtxZ);
  fChain->SetBranchAddress("primVtxCovXX",&primVtxCovXX);
  fChain->SetBranchAddress("primVtxCovYY",&primVtxCovYY);
  fChain->SetBranchAddress("primVtxCovZZ",&primVtxCovZZ);
  fChain->SetBranchAddress("primVtxCovXY",&primVtxCovXY);
  fChain->SetBranchAddress("primVtxCovYZ",&primVtxCovYZ);
  fChain->SetBranchAddress("primVtxCovXZ",&primVtxCovXZ);
  fChain->SetBranchAddress("primVtxChi2",&primVtxChi2);
  fChain->SetBranchAddress("primVtxNdof",&primVtxNdof);
  fChain->SetBranchAddress("BCountFilter",&BCountFilter);
  fChain->SetBranchAddress("DchTrig",&DchTrig);
  fChain->SetBranchAddress("EmcTrig",&EmcTrig);
  fChain->SetBranchAddress("R2All",&R2All);
  fChain->SetBranchAddress("nGTLFid",&nGTLFid);
  fChain->SetBranchAddress("nChgFid",&nChgFid);
  fChain->SetBranchAddress("eTotFid",&eTotFid);
  fChain->SetBranchAddress("PrimVtxdr",&PrimVtxdr);
  fChain->SetBranchAddress("PrimVtxdz",&PrimVtxdz);
  fChain->SetBranchAddress("VtxProb",&VtxProb);
  fChain->SetBranchAddress("beamSX",&beamSX);
  fChain->SetBranchAddress("beamSY",&beamSY);
  fChain->SetBranchAddress("beamSZ",&beamSZ);
  fChain->SetBranchAddress("beamSCovXX",&beamSCovXX);
  fChain->SetBranchAddress("beamSCovYY",&beamSCovYY);
  fChain->SetBranchAddress("beamSCovZZ",&beamSCovZZ);
  fChain->SetBranchAddress("beamSCovXZ",&beamSCovXZ);
  fChain->SetBranchAddress("pxUps",&pxUps);
  fChain->SetBranchAddress("pyUps",&pyUps);
  fChain->SetBranchAddress("pzUps",&pzUps);
  fChain->SetBranchAddress("eUps",&eUps);
  fChain->SetBranchAddress("nTrkTot",&nTrkTot);
  fChain->SetBranchAddress("W2",&W2);
  fChain->SetBranchAddress("FoxWol2",&FoxWol2);
  fChain->SetBranchAddress("FoxWol2Neu",&FoxWol2Neu);
  fChain->SetBranchAddress("thrust",&thrust);
  fChain->SetBranchAddress("thrustNeu",&thrustNeu);
  if (fIsMC) {
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
  if (fIsMC) fChain->SetBranchAddress("MCB0",MCB0);
  fChain->SetBranchAddress("mseB0",mseB0);
  fChain->SetBranchAddress("mHatB0",mHatB0);
  fChain->SetBranchAddress("deltaeB0",deltaeB0);
  fChain->SetBranchAddress("ThruB0",ThruB0);
  fChain->SetBranchAddress("thThruB0",thThruB0);
  fChain->SetBranchAddress("phiThruB0",phiThruB0);
  fChain->SetBranchAddress("cosTBB0",cosTBB0);
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
  fChain->SetBranchAddress("MassPB0",MassPB0);
  fChain->SetBranchAddress("MassPhiB0",MassPhiB0);
  fChain->SetBranchAddress("MassThetaB0",MassThetaB0);
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
  fChain->SetBranchAddress("nChB",&nChB);
  fChain->SetBranchAddress("massChB",massChB);
  fChain->SetBranchAddress("pChB",pChB);
  fChain->SetBranchAddress("thChB",thChB);
  fChain->SetBranchAddress("phiChB",phiChB);
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
  if (fIsMC) fChain->SetBranchAddress("MCChB",MCChB);
  fChain->SetBranchAddress("mseChB",mseChB);
  fChain->SetBranchAddress("mHatChB",mHatChB);
  fChain->SetBranchAddress("deltaeChB",deltaeChB);
  fChain->SetBranchAddress("ThruChB",ThruChB);
  fChain->SetBranchAddress("thThruChB",thThruChB);
  fChain->SetBranchAddress("phiThruChB",phiThruChB);
  fChain->SetBranchAddress("cosTBChB",cosTBChB);
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
  fChain->SetBranchAddress("MassPChB",MassPChB);
  fChain->SetBranchAddress("MassPhiChB",MassPhiChB);
  fChain->SetBranchAddress("MassThetaChB",MassThetaChB);
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
  if (fIsMC) fChain->SetBranchAddress("MCDstar",MCDstar);
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
  if (fIsMC) fChain->SetBranchAddress("MCDstar0",MCDstar0);
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
  if (fIsMC) fChain->SetBranchAddress("MCD0",MCD0);
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
  if (fIsMC) fChain->SetBranchAddress("MCChD",MCChD);
  fChain->SetBranchAddress("d1ChDIndex",d1ChDIndex);
  fChain->SetBranchAddress("d1ChDLund",d1ChDLund);
  fChain->SetBranchAddress("d2ChDIndex",d2ChDIndex);
  fChain->SetBranchAddress("d2ChDLund",d2ChDLund);
  fChain->SetBranchAddress("d3ChDIndex",d3ChDIndex);
  fChain->SetBranchAddress("d3ChDLund",d3ChDLund);
  fChain->SetBranchAddress("d4ChDIndex",d4ChDIndex);
  fChain->SetBranchAddress("d4ChDLund",d4ChDLund);
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
  if (fIsMC) fChain->SetBranchAddress("MCJpsi",MCJpsi);
  fChain->SetBranchAddress("d1JpsiIndex",d1JpsiIndex);
  fChain->SetBranchAddress("d1JpsiLund",d1JpsiLund);
  fChain->SetBranchAddress("d1JpsiGamIndex",d1JpsiGamIndex);
  fChain->SetBranchAddress("d1JpsiGamBrIndex",d1JpsiGamBrIndex);
  fChain->SetBranchAddress("d1JpsiGamNumBr",d1JpsiGamNumBr);
  fChain->SetBranchAddress("d2JpsiIndex",d2JpsiIndex);
  fChain->SetBranchAddress("d2JpsiLund",d2JpsiLund);
  fChain->SetBranchAddress("d2JpsiGamIndex",d2JpsiGamIndex);
  fChain->SetBranchAddress("d2JpsiGamBrIndex",d2JpsiGamBrIndex);
  fChain->SetBranchAddress("d2JpsiGamNumBr",d2JpsiGamNumBr);
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
  if (fIsMC) fChain->SetBranchAddress("MCKs",MCKs);
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
  if (fIsMC) fChain->SetBranchAddress("MCPi0",MCPi0);
  fChain->SetBranchAddress("d1Pi0Index",d1Pi0Index);
  fChain->SetBranchAddress("d1Pi0Lund",d1Pi0Lund);
  fChain->SetBranchAddress("d2Pi0Index",d2Pi0Index);
  fChain->SetBranchAddress("d2Pi0Lund",d2Pi0Lund);
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
  if (fIsMC) fChain->SetBranchAddress("MCGConv",MCGConv);
  fChain->SetBranchAddress("d1GConvIndex",d1GConvIndex);
  fChain->SetBranchAddress("d1GConvLund",d1GConvLund);
  fChain->SetBranchAddress("d2GConvIndex",d2GConvIndex);
  fChain->SetBranchAddress("d2GConvLund",d2GConvLund);
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
  if (fIsMC) fChain->SetBranchAddress("MCDalitz",MCDalitz);
  fChain->SetBranchAddress("d1DalitzIndex",d1DalitzIndex);
  fChain->SetBranchAddress("d1DalitzLund",d1DalitzLund);
  fChain->SetBranchAddress("d2DalitzIndex",d2DalitzIndex);
  fChain->SetBranchAddress("d2DalitzLund",d2DalitzLund);
  fChain->SetBranchAddress("nTrk",&nTrk);
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
  fChain->SetBranchAddress("chargeTrk",chargeTrk);
  fChain->SetBranchAddress("momentumTrk",momentumTrk);
  fChain->SetBranchAddress("ppcov00",ppcov00);
  fChain->SetBranchAddress("ppcov10",ppcov10);
  fChain->SetBranchAddress("ppcov11",ppcov11);
  fChain->SetBranchAddress("ppcov20",ppcov20);
  fChain->SetBranchAddress("ppcov21",ppcov21);
  fChain->SetBranchAddress("ppcov22",ppcov22);
  fChain->SetBranchAddress("xPocaTrk",xPocaTrk);
  fChain->SetBranchAddress("yPocaTrk",yPocaTrk);
  fChain->SetBranchAddress("zPocaTrk",zPocaTrk);
  fChain->SetBranchAddress("thetaTrk",thetaTrk);
  fChain->SetBranchAddress("phiTrk",phiTrk);
  fChain->SetBranchAddress("muonIdTrk",muonIdTrk);
  fChain->SetBranchAddress("elecIdTrk",elecIdTrk);
  fChain->SetBranchAddress("kaonIdTrk",kaonIdTrk);
  fChain->SetBranchAddress("pionIdTrk",pionIdTrk);
  fChain->SetBranchAddress("protonIdTrk",protonIdTrk);
  if (fIsMC) fChain->SetBranchAddress("idTrk",idTrk);
  if (fIsMC) fChain->SetBranchAddress("IndexTrk",IndexTrk);
  if (fIsMC) fChain->SetBranchAddress("IndexNtTrk",IndexNtTrk);
  fChain->SetBranchAddress("B0RecTrk",B0RecTrk);
  fChain->SetBranchAddress("chBRecTrk",chBRecTrk);
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
  if (fIsMC) fChain->SetBranchAddress("idGam",idGam);
  if (fIsMC) fChain->SetBranchAddress("IndexGam",IndexGam);
  if (fIsMC) fChain->SetBranchAddress("IndexNtGam",IndexNtGam);
  fChain->SetBranchAddress("B0RecGam",B0RecGam);
  fChain->SetBranchAddress("chBRecGam",chBRecGam);
  Notify();
}

// ----------------------------------------------------------------------
Bool_t b2uNtp::Notify() {
  // Called when loading a new file.
  // Get branch pointers.
  b_event = fChain->GetBranch("event");
  b_runNumber = fChain->GetBranch("runNumber");
  b_platform = fChain->GetBranch("platform");
  b_partition = fChain->GetBranch("partition");
  b_upperID = fChain->GetBranch("upperID");
  b_lowerID = fChain->GetBranch("lowerID");
  b_primVtxX = fChain->GetBranch("primVtxX");
  b_primVtxY = fChain->GetBranch("primVtxY");
  b_primVtxZ = fChain->GetBranch("primVtxZ");
  b_primVtxCovXX = fChain->GetBranch("primVtxCovXX");
  b_primVtxCovYY = fChain->GetBranch("primVtxCovYY");
  b_primVtxCovZZ = fChain->GetBranch("primVtxCovZZ");
  b_primVtxCovXY = fChain->GetBranch("primVtxCovXY");
  b_primVtxCovYZ = fChain->GetBranch("primVtxCovYZ");
  b_primVtxCovXZ = fChain->GetBranch("primVtxCovXZ");
  b_primVtxChi2 = fChain->GetBranch("primVtxChi2");
  b_primVtxNdof = fChain->GetBranch("primVtxNdof");
  b_BCountFilter = fChain->GetBranch("BCountFilter");
  b_DchTrig = fChain->GetBranch("DchTrig");
  b_EmcTrig = fChain->GetBranch("EmcTrig");
  b_R2All = fChain->GetBranch("R2All");
  b_nGTLFid = fChain->GetBranch("nGTLFid");
  b_nChgFid = fChain->GetBranch("nChgFid");
  b_eTotFid = fChain->GetBranch("eTotFid");
  b_PrimVtxdr = fChain->GetBranch("PrimVtxdr");
  b_PrimVtxdz = fChain->GetBranch("PrimVtxdz");
  b_VtxProb = fChain->GetBranch("VtxProb");
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
  b_d1JpsiGamIndex = fChain->GetBranch("d1JpsiGamIndex");
  b_d1JpsiGamBrIndex = fChain->GetBranch("d1JpsiGamBrIndex");
  b_d1JpsiGamNumBr = fChain->GetBranch("d1JpsiGamNumBr");
  b_d2JpsiIndex = fChain->GetBranch("d2JpsiIndex");
  b_d2JpsiLund = fChain->GetBranch("d2JpsiLund");
  b_d2JpsiGamIndex = fChain->GetBranch("d2JpsiGamIndex");
  b_d2JpsiGamBrIndex = fChain->GetBranch("d2JpsiGamBrIndex");
  b_d2JpsiGamNumBr = fChain->GetBranch("d2JpsiGamNumBr");
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
  b_protonIdTrk = fChain->GetBranch("protonIdTrk");
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
void b2uNtp::Show(Int_t entry) {
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

// ----------------------------------------------------------------------
Int_t b2uNtp::Cut(Int_t entry) {
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}

// ----------------------------------------------------------------------
void  b2uNtp::printLorentz(const TLorentzVector &p) {
  char line[200];
  sprintf(line, "p = (%7.5f, %7.5f, %7.5f, %7.5f), m = %7.5f",  p.X(), p.Y(), p.Z(), p.E
(), p.Mag());
  cout << line;
}

// ----------------------------------------------------------------------
void b2uNtp::smearNeut() { 
  static Bool_t first(kTRUE);
  if (0 == fOptSmearNeut) return;
  if (runNumber < 100000) return; // no smearing for data
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
    if (first) {      
      first = kFALSE;
    }
    energyGam[i] = tempene;
  }  
  first = kFALSE;

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
b2uNtp::clusterCorrection(const double rawEnergy,
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
b2uNtp::clusterReCorrection(const double oldCorrectedEnergy, 
                    const double clusterPositionTheta) const
{

  double eCorr=oldCorrectedEnergy;
  const double eRaw=eCorr*eCorr/clusterCorrection(eCorr,clusterPositionTheta,false);
  eCorr=clusterCorrection(eRaw,clusterPositionTheta,true);

  return eCorr;

}

void b2uNtp::doSplitOffStudy(){

  double magFac=0.1; // Pt corresponding to a radius on 22 cm (the edge of the DCH) - GeV
  for (Int_t i1=0; i1<nTrk; i1++) {
    int ndch1=ndchTrk[i1];
    int nsvt1=nsvtTrk[i1];
    double phiMinTrk(1000.),thetaMinTrk(1000.);
    double thetaAtPhiMinTrk(1000.),phiAtThetaMinTrk(1000.);
    // track track splitoffs
    for (Int_t i2=i1+1; i2<nTrk; i2++) {
      int ndch2=ndchTrk[i2];
      int nsvt2=nsvtTrk[i2];
      if(
         (ndch1==0 && nsvt2==0) || (ndch2==0 && nsvt1==0)
         ){
        if(chargeTrk[i1] != chargeTrk[i2])continue;
        if(momentumTrk[i1]<0.1 || momentumTrk[i2]<0.1 )continue;

        double dph=phiTrk[i1]-phiTrk[i2];
        double corr=2*(asin(0.1/momentumTrk[i1])-asin(0.1/momentumTrk[i2]));// angular correction at 22 cm (DCH - SVT boundary)
        if(chargeTrk[i1]>0)corr=-corr;
        dph+=corr;

        double dth=thetaTrk[i1]-thetaTrk[i2];

        if(TMath::Abs(dth) < thetaMinTrk){
          thetaMinTrk = TMath::Abs(dth);
          phiAtThetaMinTrk = dph;
        }
        if(TMath::Abs(dph) < phiMinTrk){
          phiMinTrk = TMath::Abs(dph);
          thetaAtPhiMinTrk = dth;
        }
      }
    }
  }

  // track cluster splitoffs
  for (int ig=0;ig<nGam;ig++){
    double phiMinNeu(1000.),thetaMinNeu(1000.);
    double thetaAtPhiMinNeu(1000.),phiAtThetaMinNeu(1000.);
    int imin = -1;
    bool dumptree(false);
    splitOffGam[ig]=0;
    if (ecalGam[ig] < 0) continue; // no IFR only stuff
    for (Int_t i1=0; i1<nTrk; i1++) {
      double dph=phiAtEMCTrk[i1]-phiGam[ig];
      double dth=thetaAtEMCTrk[i1]-thetaGam[ig];

      if(TMath::Abs(dth) < thetaMinNeu){
        thetaMinNeu = TMath::Abs(dth);
        phiAtThetaMinNeu = dph;
      }
      if(TMath::Abs(dph) < phiMinNeu){
        phiMinNeu = TMath::Abs(dph);
        thetaAtPhiMinNeu = dth;

        imin=i1;
      }
      // select splitOff photons
      bool phicut = (chargeTrk[i1]>0 && dph>-0.03 && dph<0.07) || (chargeTrk[i1]<0 && dph>-0.07 && dph<0.03);
      bool acc = thetaGam[ig]> 0.41 && ecalTrk[i1] < 0; //track must not be matched
      if(phicut && acc && (TMath::Abs(dth)<0.03)) splitOffGam[ig]=1;
    }
  }

}

void b2uNtp::goodBreco(int lbrecoI) {

  fBrecoQual = 0;
  if (!fIsMC) return;

  TLorentzVector p4lBGen(f4BrecoGen);
  TLorentzVector p4lBReco(f4BrecoNC);
  p4lBReco.Boost(-f3UpsBoost);
  p4lBGen.Boost(-f3UpsBoost);
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
    if (!fChB){
      ndauB = ndauB0[fIndexBestB];
      lmodeB = modeB0[fIndexBestB];
    } else {
      ndauB = ndauChB[fIndexBestB];
      lmodeB = modeChB[fIndexBestB];
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
//      if(!fChB)
//        cout << d1B0Lund[fIndexBestB] << " " << d2B0Lund[fIndexBestB] << " "<< d3B0Lund[fIndexBestB] << " " << d4B0Lund[fIndexBestB] << " " << d5B0Lund[fIndexBestB] << " " << d6B0Lund[fIndexBestB] << " " << d7B0Lund[fIndexBestB] << endl;
//      else 
//        cout << d1ChBLund[fIndexBestB] << " " << d2ChBLund[fIndexBestB] << " " << d3ChBLund[fIndexBestB] << " " << d4ChBLund[fIndexBestB] << " " << d5ChBLund[fIndexBestB] << " " << d6ChBLund[fIndexBestB] << " " << d7ChBLund[fIndexBestB] << endl;
//      cout << "----------" << endl;
//      cout << lid << " " << getBflavor(fBrecoFlavor, fBrecoCharge) << " " << fTruBflavor << " " << distance << " " << nDauMc[lbrecoI] << " " << ndauB << " " << getSeedDbyimc(lbrecoI) << " " << getSeedDbymode(lmodeB) << " " << fBrecoQual << endl;
//      cout << "*************" << endl;
}

// get the seed id of the generated B
int b2uNtp::getSeedDbyimc(int inimc) {
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
int b2uNtp::getSeedDbymode(int inmodeB) {
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

  // flavor variables: 1 for B0, -1 for B+, 2 for B0bar, -2 for B-
int b2uNtp::getBflavor(int Bid) {
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

int b2uNtp::getBflavor(int BRECOflavor, int charge) {
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
