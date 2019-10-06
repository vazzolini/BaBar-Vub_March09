#include "VubAnalysisCode.hh"

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
  int abcfit_interface_vub_(int *ISMEAR, int *ILTYP, float *CVAL, float *P_REC, float *P_FIT, float
 *CHI2T,float *PROBCHI2,int *IERR, int *ISV);
}
#endif

Double_t p_energy_loss_corrected(Double_t pin, Double_t dip, Int_t itype);


VubAnalysisCode::VubAnalysisCode(TTree *tree, int isMC, int newFormat):baseClass(tree, isMC, newFormat) {
  cout << "VubAnalysisCode has been called" << endl;
  Init(tree,isMC); 
  fToBeCopied = new TEventList("toBeCopied", "Events to be copied", 1000);
  //ADDED CB
  initRest();
}

void VubAnalysisCode::makeParam(char ifile[80])  {
  fIsMakeParam=kTRUE;
  sprintf(fDump4ParamFile[0],"%s%s",ifile,"_vno");
  sprintf(fDump4ParamFile[1],"%s%s",ifile,"_vcb");
  sprintf(fDump4ParamFile[2],"%s%s",ifile,"_vub");
};

// ----------------------------------------------------------------------
void VubAnalysisCode::Skim(Double_t pCut, Int_t maxEvent, Int_t startEvent, Int_t isVerbose, const char *ITSfile) {

  
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

// ----------------------------------------------------------------------
void VubAnalysisCode::dumpEventList(const char *filename) {
  TFile *f = new TFile(filename, "RECREATE");
  fToBeCopied->Print();
  fChain->SetEventList(fToBeCopied);
  TTree *small = fChain->CopyTree("");
  small->Write();
  small->Print();
}

// ----------------------------------------------------------------------
void VubAnalysisCode::dump4Param(int type, float P[28]) {

//type 0 Vno
//type 1 Vcb
//type 2 Vub

 if ((type==0) || (type==1) || (type=2)) {
  ofstream d4Pf(fDump4ParamFile[type],ios::app);
  for (int j=0; j<4; j++) {
   for (int i=0; i<4; i++) d4Pf << P[j*4+i] << " ";
   d4Pf << endl;
  }
 } 
}


// ----------------------------------------------------------------------
void VubAnalysisCode::recoil() {

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

  fHistFile->cd("recoil");
  

  // -- Find leading lepton
  Int_t i, jbit;
  fNLepton = fNKshort = fNKp = fNlep = 0;
  for (i = 0; i < nTrk; ++i) {
    Brectrktmp = B0RecTrk[i];
    if(fChB) Brectrktmp = chBRecTrk[i];
    if ((Brectrktmp&brecoOverlap))      continue;
    if(fNewFormat == 1&& Brectrktmp==2)continue;// Apr 02 format
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
  int ispion = nGoodPi = 0;
  fMinKMom = 10000.;
  fMaxKMom = -1.;
  for (i = 0; i < nTrk; ++i) {
    if ((goodTrack[i] == 0) && (kshortLockTrk[i] == 0)) continue; // skip tracks only if they are not part of a kshort
    Brectrktmp = B0RecTrk[i];
    if(fChB) Brectrktmp = chBRecTrk[i];
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
      }
      if ((i == lmax) && (convLockTrk[i] == 1)) { 
	cout << "... -> Leading electron flagged as converted photon" << endl; 
      }
    } else if (isRecKaon(i) && (momentumTrk[i] > KAMOMLO)) {  
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
    } else if (isRecMu(i)) {
      goodHadron[i] = goodChargedKaon[i] = goodPion[i] = 0; 
      if(i==lmax) {
        mk4Vector(p4t, momentumTrk[i], thetaTrk[i], phiTrk[i], MUMASS);
      }else{
        mk4Vector(p4t, momentumTrk[i], thetaTrk[i], phiTrk[i], PIPMASS);
	ispion = 1;
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

    // -- Skip tracks which are leading leptons for the rest
    if (i == lmax) continue;
    p4Xhad += p4t;
    fcountChg++;
    p4ChargPart = p4Xhad;

    if (ispion == 1) {
      momentumGoodPi[nGoodPi] = momentumTrk[i];
      thetaGoodPi[nGoodPi] = thetaTrk[i];
      phiGoodPi[nGoodPi] = phiTrk[i];
      chargeGoodPi[nGoodPi] = chargeTrk[i];     
      nGoodPi ++;   
    }


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
    if(fChB) Brecgamtmp = chBRecGam[i];
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

  flMommin = 1000.;

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
   if(fLeptonCharge * chargeGoodPi[i] == -1 ) {
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
    fQ2Fit  = 2*pnufit*p4LeptonLab;
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
      fPcms = -98.; // reset fPcms so that the event is not dumped into 'events'
      
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

  fHistFile->cd(); 
  if (faoLepton)     ((TH1D*)gDirectory->Get("mesalleventsA"))->Fill(fMes);
  if (faoMM2)        ((TH1D*)gDirectory->Get("mesalleventsB"))->Fill(fMes);
  if (faoChargeCons) ((TH1D*)gDirectory->Get("mesalleventsC"))->Fill(fMes);
  
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

  /*  EJH EJH EJH
  if (fGoodAccLepton) {
    ((TH1D*)gDirectory->Get("derecoil"))->Fill(fDeltaE);
    ((TH1D*)gDirectory->Get("mesrecoil"))->Fill(fMes);
    if (faoLepton)     ((TH1D*)gDirectory->Get("mesrecoilA"))->Fill(fMes);
    if (faoMM2)        ((TH1D*)gDirectory->Get("mesrecoilB"))->Fill(fMes);
    if (faoChargeCons) ((TH1D*)gDirectory->Get("mesrecoilC"))->Fill(fMes);
    fillRecoilHist("recoil");
    if (mesSideband) fillRecoilHist("bgrecoil");
    if (signalBox) fillRecoilHist("sgrecoil");
    fillMesHist("recoil", "a");
    if (fGoodLepton) {
      ((TH1D*)gDirectory->Get("deall"))->Fill(fDeltaE);
      ((TH1D*)gDirectory->Get("mesall"))->Fill(fMes);
      if (faoLepton)     ((TH1D*)gDirectory->Get("mesallA"))->Fill(fMes);
      if (faoMM2)        ((TH1D*)gDirectory->Get("mesallB"))->Fill(fMes);
      if (faoChargeCons) ((TH1D*)gDirectory->Get("mesallC"))->Fill(fMes);
      if (mesSideband) fillRecoilHist("bgall");
      if (signalBox) fillRecoilHist("sgall");
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
      if (mesSideband) fillRecoilHist("bgvcb");
      if (signalBox) fillRecoilHist("sgvcb");

      fillMesHist("sgvcb", "a");

    } else {
      ((TH1D*)gDirectory->Get("devub"))->Fill(fDeltaE);
      ((TH1D*)gDirectory->Get("mesvub"))->Fill(fMes);
      if (faoLepton)     ((TH1D*)gDirectory->Get("mesvubA"))->Fill(fMes);
      if (faoMM2)        ((TH1D*)gDirectory->Get("mesvubB"))->Fill(fMes);
      if (faoChargeCons) ((TH1D*)gDirectory->Get("mesvubC"))->Fill(fMes);
      if (mesSideband) fillRecoilHist("bgvub");
      if (signalBox) fillRecoilHist("sgvub");

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
  EJH EJH EJH */

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

// ----------------------------------------------------------------------

void VubAnalysisCode::bookHist(int dump) {
  TH1 *h;
  char name[100], title[100];
  char title1[100], variable[100];
  if (!fHistFile) {
    cout << "Call openHistFile(...) before booking histograms" << endl;
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
    fTree->Branch("evtw8",  &fEvtW8, "evtw8/D");
    fTree->Branch("lw8",    &fLumiW8,"lw8/D");

    // -- breco
    fTree->Branch("bmass",      &fBmass, "bmass/D");
    fTree->Branch("bmassfit",   &fBmassfit, "bmassfit/D");
    fTree->Branch("sbox",       &signalBox, "sbox/b");
    fTree->Branch("mes",        &fMes, "mes/D");
    fTree->Branch("de",         &fDeltaE, "de/D");
    fTree->Branch("pur",        &fPurity, "pur/D");
    fTree->Branch("Gvxbtyp", &fBVxbTyp, "Gvxbtyp/I");
    fTree->Branch("GSem", &fBVSTyp, "GSem/I");
    fTree->Branch("GfDpi", &fDpi, "GfDpi/I");
    fTree->Branch("GfDpiz", &fDpiz, "GfDpiz/I");
    fTree->Branch("GfDk", &fDk, "GfDk/I");
    fTree->Branch("GfDks", &fDks, "GfDks/I");
    fTree->Branch("GfDkmiss", &fDkmiss, "GfDkmiss/I");
    fTree->Branch("GfDlep", &fDlep, "GfDlep/I");
    fTree->Branch("GfDgam", &fDgam, "GfDgam/I");
    fTree->Branch("GfDnu", &fDnu, "GfDnu/I");
    fTree->Branch("GfD0Ds", &fD0CfDs, "GfD0Ds/I");
    fTree->Branch("GfDDs", &fDCfDs, "GfDDs/I");      
    fTree->Branch("GfDkl", &fDkl, "GfDkl/I");
    fTree->Branch("GfDkspiopio", &fDkspiopio, "GfDkspiopio/I");
    fTree->Branch("GfDkspipi", &fDkspipi, "GfDkspipi/I");
    fTree->Branch("intpur",     &fIntPurity, "intpur/D");
    fTree->Branch("brecoflav",  &fBrecoFlavor, "brecoflav/I");
    fTree->Branch("brecocharge",&fBrecoCharge , "brecocharge/I");
    fTree->Branch("brecomc",    &fBrecoMc,  "brecomc/I");
    fTree->Branch("mode",       &fBmode   ,  "mode/I");    
    fTree->Branch("nnpi0",      &nnpi0   ,  "nnpi0/I");
    fTree->Branch("nnks",       &nnks   ,  "nnks/I");
    fTree->Branch("nnpar",      &nnpar   ,  "nnpar/I");

    // -- generator block
    fTree->Branch("fBchgen", &fBbchgen, "fBchgen/I");

    fTree->Branch("pcmsgen",  &fPcmsGen, "pcmsgen/D");
    fTree->Branch("tcmsgen",  &fTcmsGen, "tcmsgen/D");
    fTree->Branch("fcmsgen",  &fFcmsGen, "fcmsgen/D");
    fTree->Branch("ecmsgen",  &fEcmsGen, "ecmsgen/D");

    fTree->Branch("mxhadgen", &fMxhadGen, "mxhadgen/D");
    fTree->Branch("pxhadgen", &fPxhadGen, "pxhadgen/D");
    fTree->Branch("txhadgen", &fTxhadGen, "txhadgen/D");
    fTree->Branch("fxhadgen", &fFxhadGen, "fxhadgen/D");
    fTree->Branch("exhadgen", &fExhadGen, "exhadgen/D");

    /*
      Added Ciuchini variables 
    */
    fTree->Branch("GxCiuc", &fGxCiuc, "GxCiuc/D");
    fTree->Branch("GwCiuc", &fGwCiuc, "GwCiuc/D");
    fTree->Branch("GcsiCiuc", &fGcsiCiuc, "GcsiCiuc/D");
    fTree->Branch("xCiuc", &fxCiuc, "xCiuc/D");
    fTree->Branch("wCiuc", &fwCiuc, "wCiuc/D");
    fTree->Branch("csiCiuc", &fcsiCiuc, "csiCiuc/D");
    /*
      Added variables to study pW+ EW
    */
    fTree->Branch("EwPwfit", &fEwPwfit, "EwPwfit/D");
    fTree->Branch("EwPw"   , &fEwPw,    "EwPw/D");
    fTree->Branch("EwPwG"  , &fEwPwG,   "EwPwG/D");
    /*
      Finished
    */


    fTree->Branch("kplus",    &fKplus, "fkplus/D");
    fTree->Branch("GoodEvent", &fGoodEvent, "GoodEvent/b");

    fTree->Branch("vub",    &fVub, "vub/I");
    fTree->Branch("vcb",    &fVcb, "vcb/I");
    fTree->Branch("vxbtyp", &fBVxbTyp, "vxbtyp/I");
    fTree->Branch("other",  &fOther, "other/I");
    fTree->Branch("wKK",    &fWithKKbar, "wKK/I");


    // -- recoil
    fTree->Branch("xcharge", &fRecoilCharge, "xcharge/I"); // note that this is the entire recoil, not xhad!
    fTree->Branch("qtot",    &fQtot, "qtot/D");            // qtot is qtot
    fTree->Branch("pxhad",   &fPxhad, "pxhad/D");
    fTree->Branch("txhad",   &fTxhad, "txhad/D");
    fTree->Branch("fxhad",   &fFxhad, "fxhad/D");
    fTree->Branch("exhad",   &fExhad, "exhad/D");
    fTree->Branch("mxhad",   &fMxhad, "mxhad/D");
    fTree->Branch("gmax",    &fGammaMax,"gmax/D");
    fTree->Branch("mxhadfit",&fMxhadfit, "mxhadfit/D");
    fTree->Branch("mxminpi",&fp4XminPi, "mxminpi/D");
    fTree->Branch("bestklind",&fBestKlind, "bestklind/I");
    fTree->Branch("deltam",&fdeltaM, "deltam/D");
    fTree->Branch("mxwminpi",&fwp4XminPi, "mxwminpi/D");
    fTree->Branch("wdeltam",&fwdeltaM, "wdeltam/D");
    fTree->Branch("mxphotmin",&fphotp4Xmin, "mxphotmin/D");
    fTree->Branch("photdeltam",&fphotdeltaM, "photdeltam/D");
    fTree->Branch("minlat",&flMommin, "minlat/D");

  //mxhadchg properties in new flags 

   fTree->Branch("pxhadchg",&fPxhadchg, "fPxhadchg/D");
   fTree->Branch("txhadchg",&fTxhadchg, "fTxhadchg/D");
   fTree->Branch("fxhadchg",&fFxhadchg, "fFxhadchg/D");
   fTree->Branch("exhadchg",&fExhadchg, "fExhadchg/D");
   fTree->Branch("mxhadchg",&fMxhadchg, "fMxhadchg/D");

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
    fTree->Branch("npi0",  &fRecoilPi0Mult, "npi0/I");
    fTree->Branch("nneu",  &fRecoilNutMult, "nneu/I");
    fTree->Branch("nneu80_160",  &fRecoilNutMult80_160, "nneu80_160/I");
    fTree->Branch("nneu160_320",  &fRecoilNutMult160_320, "nneu160_320/I");
    fTree->Branch("nneufromB",  &fRecoilNutMultfromB, "nneufromB/I");
    fTree->Branch("nneufromB80_160",  &fRecoilNutMultfromB80_160, "nneufromB80_160/I");
    fTree->Branch("nneufromB160_320",  &fRecoilNutMultfromB160_320, "nneufromB160_320/I");
    fTree->Branch("nkp",  &fNKp, "nkp/I");
    fTree->Branch("nks",  &fNKshort, "nks/I");
    fTree->Branch("totweight",  &totweight, "totweight/D");
    fTree->Branch("totweightNutMult",  &totweightfRecoilNutMult, "totweightfRecoilNutMult/D");
    fTree->Branch("totweightTrkMult",  &totweightfRecoilTrkMult, "totweightfRecoilTrkMult/D");

    // -- neutrino
    fTree->Branch("enu",    &fEmiss, "enu/D");
    fTree->Branch("pnu",    &fPNu, "pnu/D");
    fTree->Branch("tnu",    &fTNu, "tnu/D");
    fTree->Branch("fnu",    &fFNu, "fnu/D");
    fTree->Branch("mm2",    &fMM2, "mm2/D");
    fTree->Branch("mm2nc",  &fMM2NC, "mm2nc/D");
    fTree->Branch("mm2fit", &fMM2fit, "mm2fit/D");
    fTree->Branch("ENeu",    &fENeu, "Eneu/D");
    fTree->Branch("EPiz",    &fEPiz, "EPiz/D");
    fTree->Branch("MinKMom",    &fMinKMom, "MinKMom/D");
    fTree->Branch("MaxKMom",    &fMaxKMom, "MaxKMom/D");

    /*#ifndef FAST
    cout << "Booking vertex stuff, FAST ist not enabled" << endl;
    // vertex info
    fTree->Branch("dx",&fDx,"dx/D");
    fTree->Branch("dy",&fDy,"dy/D");
    fTree->Branch("dz",&fDz,"dz/D");
    fTree->Branch("s2dxx",&fS2Dxx,"s2dxx/D");
    fTree->Branch("s2dyy",&fS2Dyy,"s2dyy/D");
    fTree->Branch("s2dzz",&fS2Dzz,"s2dzz/D");
    fTree->Branch("s2dxy",&fS2Dxy,"s2dxy/D");
    fTree->Branch("s2dyz",&fS2Dyz,"s2dyz/D");
    fTree->Branch("s2dxz",&fS2Dxz,"s2dxz/D");
    #endif*/

    fTree->Branch("q2",     &fQ2,    "q2/D");
    fTree->Branch("q2nc",   &fQ2NC,  "q2nc/D");
    fTree->Branch("q2fit",  &fQ2Fit, "q2fit/D");
    fTree->Branch("q2Gen",  &fQ2Gen, "q2Gen/D");

    if (fOptGammas) {  
      fGTree = new TTree("photons", "photons"); 
      for(int ju =0; ju<15; ju++) {
	sprintf(title1, "G_NuM%d", ju);    sprintf(variable, "G_NuM%d%s", ju, "/D");
	fGTree->Branch(title1,    &tmpfMM2[ju], variable);
	sprintf(title1, "G_NuT%d", ju);    sprintf(variable, "G_NuT%d%s", ju, "/D"); 
	fGTree->Branch(title1,    &tmpfNuT[ju], variable);
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
      fGTree->Branch("ContKs", &fcontKs, "ContKs/b");
      fGTree->Branch("ChargeCorr", &fGoodChargeCorr, "ChargeCorr/b");
      fGTree->Branch("ChargeCons", &fGoodChargeCons, "ChargeCons/b");
      fGTree->Branch("GoodLep", &fGoodLepton, "GoodLep/b");
      fGTree->Branch("Masschpa", &fMCharpart, "Masschpa/D");
      fGTree->Branch("Enchpa", &fECharpart, "Enchpa/D");
      fGTree->Branch("Thchpa", &fTCharpart, "Thchpa/D");
      fGTree->Branch("Phichpa", &fPCharpart, "Phichpa/D");
      fGTree->Branch("sbox",    &signalBox,   "sbox/b");
      fGTree->Branch("sideb",   &mesSideband, "sideb/b");
      fGTree->Branch("Gvxbtyp", &fBVxbTyp, "Gvxbtyp/I");
      fGTree->Branch("GfDpi", &fDpi, "GfDpi/I");
      fGTree->Branch("GfDpiz", &fDpiz, "GfDpiz/I");
      fGTree->Branch("GfDk", &fDk, "GfDk/I");
      fGTree->Branch("GfDkmiss", &fDkmiss, "GfDkmiss/I");
      fGTree->Branch("GfDks", &fDks, "GfDks/I");
      fGTree->Branch("GfDkl", &fDkl, "GfDkl/I");
      fGTree->Branch("GfDkspiopio", &fDkspiopio, "GfDkspiopio/I");
      fGTree->Branch("GfDkspipi", &fDkspipi, "GfDkspipi/I");
      fGTree->Branch("GfDlep", &fDlep, "GfDlep/I");
      fGTree->Branch("GfD0Ds", &fD0CfDs, "GfD0Ds/I");
      fGTree->Branch("GfDDs", &fDCfDs, "GfDDs/I");      
    }
    
    /*#ifndef FAST
    // MX study variables
    fTree->Branch("allksm0",fallKsm0 , "allksm0[nks]/D");
    fTree->Branch("allksp",fallKsp , "allksp[nks]/D");
    fTree->Branch("allksmc",fallKsMc , "allksmc[nks]/b");
    fTree->Branch("allchkp",fallchkp , "allchkp[nkp]/D");
    fTree->Branch("allchkmc",fallchkMc , "allchkmc[nkp]/b");
    fTree->Branch("m0ks", &fm0ks, "m0ks/D");
    fTree->Branch("pks", &fpks, "pks/D");
    fTree->Branch("pksmc", &fksmatchp, "pksmc/D");
 
    // MC TRUTH KAONS
    if (fIsMC ) {
      fTree->Branch("ntkl",&fnKL,"ntkl/I");
      fTree->Branch("tklp",ftklp,"tklp[ntkl]/D");
      fTree->Branch("tklth",ftklth,"tklth[ntkl]/D");
      fTree->Branch("tklph",ftklph,"tklph[ntkl]/D");
      fTree->Branch("tklisol",ftklisol,"tklisol[ntkl]/b");
      fTree->Branch("ntks",&fntks,"ntks/I");
      fTree->Branch("tksp",ftksp,"tksp[ntks]/D");
      fTree->Branch("tksth",ftksth,"tksth[ntks]/D");
      fTree->Branch("tksph",ftksph,"tksph[ntks]/D");
      fTree->Branch("tksdec",ftksdecay,"tksdec[ntks]/I");
      fTree->Branch("ntchk",&fntchk,"ntchk/I");
      fTree->Branch("tchkp",ftchkp,"tchkp[ntchk]/D");
      fTree->Branch("tchkth",ftchkth,"tchkth[ntchk]/D");
      fTree->Branch("tchkph",ftchkph,"tchkph[ntchk]/D");

      // KL match resolution && energy deposit from MC
      fTree->Branch("nklres",&fklreslen,"nklres/I");
      fTree->Branch("klresth",fklresth,"klresth[nklres]/D");
      fTree->Branch("klresph",fklresph,"klresph[nklres]/D");
      fTree->Branch("klid",fidCone,"klid[nklres]/I");
      fTree->Branch("klcone",finCone,"klcone[nklres]/b");
      fTree->Branch("emckl",&fsumklemcen , "emckl/D");
      fTree->Branch("emckl0",&fsumklemcen0 , "emckl0/D");
      fTree->Branch("emckl22",&fsumklemcen22 , "emckl22/D");
 
    }
    // MX reco 
    if (fOptCategories) {
      
      fTree->Branch("bgcat",  &fMxCategory, "bgcat/I");
      fTree->Branch("mxks",  &fMxKs, "mxks/D");
      fTree->Branch("mm2ks",         &fMM2Ks, "mm2ks/D");
      fTree->Branch("mxksfit",         &fMxFitKs, "mxksfit/D");
      
      fTree->Branch("mm2misk",&fMM2MisK,"mm2misk/D");
      fTree->Branch("mxmisk",&fMxMisK,"mxmisk/D");
      fTree->Branch("mxmiskfit",&fMxFitMisK,"mxmiskfit/D");
    
      fTree->Branch("mm2chk",&fMM2chK,"mm2mchk/D");
      fTree->Branch("mxchk",&fMxchK,"mxchk/D");
      fTree->Branch("mxchkfit",&fMxFitchK,"mxchkfit/D");
    }
    #endif*/

  }

  fHistFile->mkdir("mcTruth", "mcTruth");
  fHistFile->cd("mcTruth");

  sprintf(name, "TrueKlEn");  sprintf(title, "E(KL)");  h = new TH1D(name, title, 50, 0., 5.); 
  sprintf(name, "TrueKlEnNz");  sprintf(title, "E(KL)");  h = new TH1D(name, title, 50, 0., 5.); 
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
  /*#ifndef FAST
  // ks bug studies
  sprintf(name, "h77001");  sprintf(title, "ksDauEnergy");  h = new TH1D(name, title, 50, 0., 2.5); 
  sprintf(name, "h77002");  sprintf(title, "nonksDauEnergy");  h = new TH1D(name, title, 100, 0., 3.); 
  sprintf(name, "h77003");  sprintf(title, "ksMomentum");  h = new TH1D(name, title, 100, 0., 3.); 
  #endif*/
  
  sprintf(name, "h121000");  sprintf(title, "generator pcms e/mu");  h = new TH1D(name, title, PCBIN, 0., PCMAX); 
  sprintf(name, "h121005");  sprintf(title, "generator pcms e/mu Vcb"); h = new TH1D(name, title, PCBIN, 0., PCMAX); 
  sprintf(name, "h121006");  sprintf(title, "generator pcms e/mu Vub"); h = new TH1D(name, title, PCBIN, 0., PCMAX); 

  sprintf(name, "h123000");  sprintf(title, "generator Mx one lepton");  h = new TH1D(name, title, XBIN, 0., XMAX); 
  sprintf(name, "h123005");  sprintf(title, "generator Mx one lepton Vcb"); h = new TH1D(name, title, XBIN, 0., XMAX); 
  sprintf(name, "h123006");  sprintf(title, "generator Mx one lepton Vub"); h = new TH1D(name, title, XBIN, 0., XMAX); 

  sprintf(name, "h124000");  sprintf(title, "generator N hadrons with one lepton"); h = new TH1D(name, title, 20, 0., 20.); 
  sprintf(name, "h124005");  sprintf(title, "generator N hadrons with one lepton Vcb"); h = new TH1D(name, title, 20, 0., 20.); 
  sprintf(name, "h124006");  sprintf(title, "generator N hadrons with one lepton Vub"); h = new TH1D(name, title, 20, 0., 20.); 
  

  fHistFile->cd();
  fHistFile->mkdir("SplitOff", "SplitOff");     fHistFile->cd("SplitOff");
  doSplitOffStudy(true);

  if (fOptGammas) {
    fHistFile->cd();
    cout<<"making the directory"<<endl;
    fHistFile->mkdir("Photon", "Photon");     fHistFile->cd("Photon");
    
    //Photon block
    sprintf(name, "ThetaGam");  sprintf(title, "Cos(Theta) (E<80 MeV)");  h = new TH1D(name, title, 256, -1., 1.0); h->Sumw2();
    sprintf(name, "ThetaGam2");  sprintf(title, "Cos(Theta) (E<80 MeV)");  h = new TH1D(name, title, 256, 0., 3.); h->Sumw2();
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
    sprintf(name, "EnGam");  sprintf(title, "Energy");  h = new TH1D(name, title, 100, 0., 3.); h->Sumw2();
    sprintf(name, "MothGam");  sprintf(title, "MothGam");  h = new TH1D(name, title, 20000, -10000., 10000.); h->Sumw2();
    sprintf(name, "EnGamBo");  sprintf(title, "Energy");  h = new TH1D(name, title, 100, 0., 3.); h->Sumw2();
    sprintf(name, "EnGamKs");  sprintf(title, "Energy");  h = new TH1D(name, title, 100, 0., 1.); h->Sumw2();
    sprintf(name, "EnGamKspi");  sprintf(title, "Energy");  h = new TH1D(name, title, 100, 0., 1.); h->Sumw2();
    sprintf(name, "EnGamKspiz");  sprintf(title, "Energy");  h = new TH1D(name, title, 100, 0., 1.); h->Sumw2();
    sprintf(name, "EnGamNoKs");  sprintf(title, "Energy");  h = new TH1D(name, title, 100, 0., 1.); h->Sumw2();
    sprintf(name, "lMomGam");  sprintf(title, "LAT");  h = new TH1D(name, title, 256, 0.002, 1.); h->Sumw2();
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
    sprintf(name, "lMomGam (MA)");  sprintf(title, "l mom ");  h = new TH1D(name, title, 256, 0., 1.); h->Sumw2();
    sprintf(name, "ZMom42Gam (MA)");  sprintf(title, "Z42 mom ");  h = new TH1D(name, title, 100, 0., 0.2); h->Sumw2();
    sprintf(name, "PhiGam (MA)");  sprintf(title, "Phi of the recoil gamma");  h = new TH1D(name, title, 100, -3.15, 3.15); h->Sumw2();
    
    sprintf(name, "EnGam (MA) cut");  sprintf(title, "Energy");  h = new TH1D(name, title, 100, 0., 0.08); h->Sumw2();
    sprintf(name, "lMomGam (MA) cut");  sprintf(title, "l mom ");  h = new TH1D(name, title, 256, 0., 1.); h->Sumw2();
    sprintf(name, "ZMom42Gam (MA) cut");  sprintf(title, "Z42 mom ");  h = new TH1D(name, title, 100, 0., 0.2); h->Sumw2();
    sprintf(name, "PhiGam (MA) cut");  sprintf(title, "Phi of the recoil gamma");  h = new TH1D(name, title, 100, -3.15, 3.15); h->Sumw2();
    sprintf(name, "E vs theta (MA)");  sprintf(title, "Theta vs E");  h = new TH2D(name, title, 100, -1., 1., 100, 0., 2.);
    
    sprintf(name, "EnGam (MB)");  sprintf(title, "Energy");  h = new TH1D(name, title, 100, 0.08, 1.); h->Sumw2();
    sprintf(name, "lMomGam (MB)");  sprintf(title, "l mom ");  h = new TH1D(name, title, 256, 0., 1.); h->Sumw2();
    sprintf(name, "ZMom42Gam (MB)");  sprintf(title, "Z42 mom ");  h = new TH1D(name, title, 100, 0., 0.2); h->Sumw2();
    sprintf(name, "PhiGam (MB)");  sprintf(title, "Phi of the recoil gamma");  h = new TH1D(name, title, 100, -3.15, 3.15); h->Sumw2();

    sprintf(name, "EnGam (MB) cut");  sprintf(title, "Energy");  h = new TH1D(name, title, 100, 0., 0.08); h->Sumw2();
    sprintf(name, "lMomGam (MB) cut");  sprintf(title, "l mom ");  h = new TH1D(name, title, 256, 0., 1.); h->Sumw2();
    sprintf(name, "ZMom42Gam (MB) cut");  sprintf(title, "Z42 mom ");  h = new TH1D(name, title, 100, 0., 0.2); h->Sumw2();
    sprintf(name, "PhiGam (MB) cut");  sprintf(title, "Phi of the recoil gamma");  h = new TH1D(name, title, 100, -3.15, 3.15); h->Sumw2();
    sprintf(name, "E vs theta (MB)");  sprintf(title, "Theta vs E");  h = new TH2D(name, title, 100, -1., 1., 100, 0., 2.);

    sprintf(name, "EnGam (NMP)");  sprintf(title, "Energy");  h = new TH1D(name, title, 100, 0.08, 1.); h->Sumw2();
    sprintf(name, "lMomGam (NMP)");  sprintf(title, "l mom ");  h = new TH1D(name, title, 256, 0., 1.); h->Sumw2();
    sprintf(name, "ZMom42Gam (NMP)");  sprintf(title, "Z42 mom ");  h = new TH1D(name, title, 100, 0., 0.2); h->Sumw2();
    sprintf(name, "PhiGam (NMP)");  sprintf(title, "Phi of the recoil gamma");  h = new TH1D(name, title, 100, -3.15, 3.15); h->Sumw2();
    
    sprintf(name, "EnGam (NMP) cut");  sprintf(title, "Energy");  h = new TH1D(name, title, 100, 0., 0.08); h->Sumw2();
    sprintf(name, "lMomGam (NMP) cut");  sprintf(title, "l mom ");  h = new TH1D(name, title, 256, 0., 0.4); h->Sumw2();
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

  sprintf(name, "alldeltam");  sprintf(title, "DeltaM all events");  h = new TH1D(name, title, 50, 0., 1.0); 
  sprintf(name, "alldeltam2");  sprintf(title, "DeltaM all events zoom");  h = new TH1D(name, title, 50, .13, .250);
  sprintf(name, "alldeltamw");  sprintf(title, "DeltaM all events wrong sign");  h = new TH1D(name, title, 50, 0., 1.0); 
  sprintf(name, "alldeltamw2");  sprintf(title, "DeltaM all events wrong sign zoom");  h = new TH1D(name, title, 50, .13, .250); 
  sprintf(name, "alldeltamph");  sprintf(title, "DeltaM phot all events");  h = new TH1D(name, title, 50, 0., 1.0); h->Sumw2();
  sprintf(name, "alldeltamph2");  sprintf(title, "DeltaM phot all events zoom");  h = new TH1D(name, title, 50, .13, .250);

  sprintf(name, "bw8");  sprintf(title, "B decay weights");  h = new TH1D(name, title, 200, 0., 2.0); h->Sumw2();
  sprintf(name, "dw8");  sprintf(title, "B decay weights");  h = new TH1D(name, title, 200, 0., 2.0); h->Sumw2();
  sprintf(name, "ew8");  sprintf(title, "Evt weights");  h = new TH1D(name, title, 200, 0., 2.0); h->Sumw2();

  sprintf(name, "dca");  sprintf(title, "dca");  h = new TH1D(name, title, 300, 0., 6.0); h->Sumw2();
  sprintf(name, "dcaz");  sprintf(title, "dcaz");  h = new TH1D(name, title, 220, -11., 11.0); h->Sumw2();

  sprintf(name, "adca");  sprintf(title, "after cuts dca");  h = new TH1D(name, title, 300, 0., 6.0); h->Sumw2();
  sprintf(name, "adcaz");  sprintf(title, "after cuts dcaz");  h = new TH1D(name, title, 220, -11., 11.0); h->Sumw2();

  sprintf(name, "edca");  sprintf(title, "el dca");  h = new TH1D(name, title, 300, 0., 6.0); h->Sumw2();
  sprintf(name, "edcaz");  sprintf(title, "el dcaz");  h = new TH1D(name, title, 220, -11., 11.0); h->Sumw2();

  sprintf(name, "mdca");  sprintf(title, "mu dca");  h = new TH1D(name, title, 300, 0., 6.0); h->Sumw2();
  sprintf(name, "mdcaz");  sprintf(title, "mu dcaz");  h = new TH1D(name, title, 220, -11., 11.0); h->Sumw2();
  
  sprintf(name, "mesalleventsS0");  sprintf(title, "mes All seed 0 ");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesalleventsS1");  sprintf(title, "mes All seed 1 ");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesalleventsS2");  sprintf(title, "mes All seed 2 ");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesalleventsS3");  sprintf(title, "mes All seed 3 ");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesalleventsBnu");  sprintf(title, "mes All Bnu");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesalleventsBch");  sprintf(title, "mes All Bch");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesallevents");  sprintf(title, "mes All");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "deallevents");  sprintf(title, "delta E All");  h = new TH1D(name, title, 40, -0.1, 0.1); 
  sprintf(name, "mesalleventsA");  sprintf(title, "mes Vub enhanced");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesalleventsB");  sprintf(title, "mes Vub enhanced");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesalleventsC");  sprintf(title, "mes Vub enhanced");  h = new TH1D(name, title, 40, 5.2, 5.3); 

  sprintf(name, "meszero");  sprintf(title, "mes zero");  h = new TH1D(name, title, 40, 5.2, 5.3);
  sprintf(name, "mesbest");  sprintf(title, "mes best");  h = new TH1D(name, title, 40, 5.2, 5.3);

  //crossfeed plots
  sprintf(name, "mesdstarcross");  sprintf(title, "mes dstar cross");  h = new TH1D(name, title, 40, 5.2, 5.3);
  sprintf(name, "mesdccross");  sprintf(title, "mes dc cross");  h = new TH1D(name, title, 40, 5.2, 5.3);
  sprintf(name, "mesdstar0cross");  sprintf(title, "mes dstar0 cross");  h = new TH1D(name, title, 40, 5.2, 5.3);
  sprintf(name, "mesd0cross");  sprintf(title, "mes d0 cross");  h = new TH1D(name, title, 40, 5.2, 5.3);
    
  sprintf(name, "mesallcutsS0");  sprintf(title, "mes All Cuts seed 0");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesallcutsS1");  sprintf(title, "mes All Cuts seed 1");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesallcutsS2");  sprintf(title, "mes All Cuts seed 2");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesallcutsS3");  sprintf(title, "mes All Cuts seed 3");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesallcutsBnu");  sprintf(title, "mes All Cuts Bnu");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesallcutsBch");  sprintf(title, "mes All Cuts Bch");  h = new TH1D(name, title, 40, 5.2, 5.3); 
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
  sprintf(name, "mestest");  sprintf(title, "mes Recoil");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "derecoil");  sprintf(title, "delta E Recoil");  h = new TH1D(name, title, 40, -0.1, 0.1); 
  sprintf(name, "mesrecoilA");  sprintf(title, "mes Vub enhanced");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesrecoilB");  sprintf(title, "mes Vub enhanced");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesrecoilC");  sprintf(title, "mes Vub enhanced");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  
  sprintf(name, "mesallS0");  sprintf(title, "mes All seed 0");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesallS1");  sprintf(title, "mes All seed 1");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesallS2");  sprintf(title, "mes All seed 2");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesallS3");  sprintf(title, "mes All seed 3");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesallBnu");  sprintf(title, "mes All Bnu");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  sprintf(name, "mesallBch");  sprintf(title, "mes All Bch");  h = new TH1D(name, title, 40, 5.2, 5.3); 
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

  sprintf(name,"kMomCorr"); sprintf(title, "kaon momentum correction"); h=new TH2D(name, title, 100, 0., .5, 100, 0., .5); h->Sumw2();

  int no(0), le(0), i(0);
  char lc[10];
  for (i = 0; i < 11; ++i) {
    if (i==0) {fHistFile->mkdir("recoil","recoil");  fHistFile->cd("recoil"); gDirectory->mkdir("CutPlots","CutPlots"); }
    if (i==1) {fHistFile->mkdir("bgall","Sideband all");  fHistFile->cd("bgall"); gDirectory->mkdir("CutPlots","CutPlots"); }
    if (i==2) {fHistFile->mkdir("sgall","Signalbox all");  fHistFile->cd("sgall"); gDirectory->mkdir("CutPlots","CutPlots"); }
    if (i==3) {fHistFile->mkdir("bgvcb","Sideband Vub depleted");fHistFile->cd("bgvcb"); gDirectory->mkdir("CutPlots","CutPlots"); }
    if (i==4) {fHistFile->mkdir("sgvcb","Signalbox Vub depleted");fHistFile->cd("sgvcb"); gDirectory->mkdir("CutPlots","CutPlots"); }
    if (i==5) {fHistFile->mkdir("bgvub","Sideband Vub enhanced");fHistFile->cd("bgvub"); gDirectory->mkdir("CutPlots","CutPlots");}
    if (i==6) {fHistFile->mkdir("sgvub","Signalbox Vub enhanced");fHistFile->cd("sgvub"); gDirectory->mkdir("CutPlots","CutPlots");}
    if (i==7) {fHistFile->mkdir("bgallevents","Sideband");fHistFile->cd("bgallevents");gDirectory->mkdir("CutPlots","CutPlots"); }
    if (i==8) {fHistFile->mkdir("sgallevents","Signalbox");fHistFile->cd("sgallevents"); gDirectory->mkdir("CutPlots","CutPlots");}
    if (i==9) {fHistFile->mkdir("bgrecoil","Sideband recoil");  fHistFile->cd("bgrecoil"); gDirectory->mkdir("CutPlots","CutPlots"); }
    if (i==10){fHistFile->mkdir("sgrecoil","Signalbox recoil");  fHistFile->cd("sgrecoil"); gDirectory->mkdir("CutPlots","CutPlots"); }
    
    
    sprintf(name,"ks100");  sprintf(title, "mass kshorts");  h = new TH1D(name, title, 50, 0.475, 0.525); h->Sumw2();
    sprintf(name,"ks101");  sprintf(title, "mass kshorts goodLepton");  h = new TH1D(name, title, 50, 0.475, 0.525); h->Sumw2();
    sprintf(name,"ks102");  sprintf(title, "mass kshorts goodMM2");  h = new TH1D(name, title, 50, 0.475, 0.525); h->Sumw2();
    sprintf(name,"ks103");  sprintf(title, "mass kshorts chargeCorr");  h = new TH1D(name, title, 50, 0.475, 0.525); h->Sumw2();
    sprintf(name,"ks104");  sprintf(title, "mass kshorts chargeCons");  h = new TH1D(name, title, 50, 0.475, 0.525); h->Sumw2();
    sprintf(name,"ks105");  sprintf(title, "mass kshorts goodEvent");  h = new TH1D(name, title, 50, 0.475, 0.525); h->Sumw2();
  
    sprintf(name,"ks200");  sprintf(title, "plab k");  h = new TH1D(name, title, 50, 0., 5.); h->Sumw2();
    sprintf(name,"ks210");  sprintf(title, "tlab k");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
    sprintf(name,"ks220");  sprintf(title, "plab k");  h = new TH1D(name, title, 50, 0., 5.); h->Sumw2();
    sprintf(name,"ks201");  sprintf(title, "plab k");  h = new TH1D(name, title, 50, 0., 5.); h->Sumw2();
    sprintf(name,"ks211");  sprintf(title, "tlab k");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
    sprintf(name,"ks221");  sprintf(title, "plab k");  h = new TH1D(name, title, 50, 0., 5.); h->Sumw2();
    sprintf(name,"ks202");  sprintf(title, "plab k");  h = new TH1D(name, title, 50, 0., 5.); h->Sumw2();
    sprintf(name,"ks212");  sprintf(title, "tlab k");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
    sprintf(name,"ks222");  sprintf(title, "plab k");  h = new TH1D(name, title, 50, 0., 5.); h->Sumw2();
    sprintf(name,"ks203");  sprintf(title, "plab k");  h = new TH1D(name, title, 50, 0., 5.); h->Sumw2();
    sprintf(name,"ks213");  sprintf(title, "tlab k");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
    sprintf(name,"ks223");  sprintf(title, "plab k");  h = new TH1D(name, title, 50, 0., 5.); h->Sumw2();
    sprintf(name,"ks204");  sprintf(title, "plab k");  h = new TH1D(name, title, 50, 0., 5.); h->Sumw2();
    sprintf(name,"ks214");  sprintf(title, "tlab k");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
    sprintf(name,"ks224");  sprintf(title, "plab k");  h = new TH1D(name, title, 50, 0., 5.); h->Sumw2();
    sprintf(name,"ks205");  sprintf(title, "plab k");  h = new TH1D(name, title, 50, 0., 5.); h->Sumw2();
    sprintf(name,"ks215");  sprintf(title, "tlab k");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
    sprintf(name,"ks225");  sprintf(title, "plab k");  h = new TH1D(name, title, 50, 0., 5.); h->Sumw2();

    sprintf(name,"kp200");  sprintf(title, "plab k");  h = new TH1D(name, title, 50, 0., 5.); h->Sumw2();
    sprintf(name,"kp210");  sprintf(title, "tlab k");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
    sprintf(name,"kp220");  sprintf(title, "plab k");  h = new TH1D(name, title, 50, 0., 5.); h->Sumw2();
    sprintf(name,"kp201");  sprintf(title, "plab k");  h = new TH1D(name, title, 50, 0., 5.); h->Sumw2();
    sprintf(name,"kp211");  sprintf(title, "tlab k");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
    sprintf(name,"kp221");  sprintf(title, "plab k");  h = new TH1D(name, title, 50, 0., 5.); h->Sumw2();
    sprintf(name,"kp202");  sprintf(title, "plab k");  h = new TH1D(name, title, 50, 0., 5.); h->Sumw2();
    sprintf(name,"kp212");  sprintf(title, "tlab k");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
    sprintf(name,"kp222");  sprintf(title, "plab k");  h = new TH1D(name, title, 50, 0., 5.); h->Sumw2();
    sprintf(name,"kp203");  sprintf(title, "plab k");  h = new TH1D(name, title, 50, 0., 5.); h->Sumw2();
    sprintf(name,"kp213");  sprintf(title, "tlab k");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
    sprintf(name,"kp223");  sprintf(title, "plab k");  h = new TH1D(name, title, 50, 0., 5.); h->Sumw2();
    sprintf(name,"kp204");  sprintf(title, "plab k");  h = new TH1D(name, title, 50, 0., 5.); h->Sumw2();
    sprintf(name,"kp214");  sprintf(title, "tlab k");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
    sprintf(name,"kp224");  sprintf(title, "plab k");  h = new TH1D(name, title, 50, 0., 5.); h->Sumw2();
    sprintf(name,"kp205");  sprintf(title, "plab k");  h = new TH1D(name, title, 50, 0., 5.); h->Sumw2();
    sprintf(name,"kp215");  sprintf(title, "tlab k");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
    sprintf(name,"kp225");  sprintf(title, "plab k");  h = new TH1D(name, title, 50, 0., 5.); h->Sumw2();

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

    /*#ifndef FAST    
    for (int ik=1; ik<21 ; ik++){
      no = 1600;		     
      sprintf(name,"%s%d%s%d","a",no,"bin",ik);   sprintf(title, "mes plab all had nocut");  h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      no = 1610;		     
      sprintf(name,"%s%d%s%d","a",no,"bin",ik);   sprintf(title, "mes plab all had all cuts");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
      no = 9900;		     
      sprintf(name,"%s%d%s%d","a",no,"bin",ik);   sprintf(title, "mes plab all neu nocut");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
      no = 9910;		     
      sprintf(name,"%s%d%s%d","a",no,"bin",ik);   sprintf(title, "mes plab all neu all cuts");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    }
    
    for (int ik=1; ik<21 ; ik++){
      no = 1600;		     
      sprintf(name,"%s%d%s%d","pur60a",no,"bin",ik);   sprintf(title, "mes plab all had nocut");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
      no = 1610;		     
      sprintf(name,"%s%d%s%d","pur60a",no,"bin",ik);   sprintf(title, "mes plab all had all cuts");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
      no = 9900;		     
      sprintf(name,"%s%d%s%d","pur60a",no,"bin",ik);   sprintf(title, "mes plab all neu nocut");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
      no = 9910;		     
      sprintf(name,"%s%d%s%d","pur60a",no,"bin",ik);   sprintf(title, "mes plab all neu all cuts");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    }
    
    
    for (int ik=1; ik<21 ; ik++){
      no = 1600;		     
      sprintf(name,"%s%d%s%d","pur80a",no,"bin",ik);   sprintf(title, "mes plab all had nocut");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
      no = 1610;		     
      sprintf(name,"%s%d%s%d","pur80a",no,"bin",ik);   sprintf(title, "mes plab all had all cuts");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
      no = 9900;		     
      sprintf(name,"%s%d%s%d","pur80a",no,"bin",ik);   sprintf(title, "mes plab all neu nocut");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
      no = 9910;		     
      sprintf(name,"%s%d%s%d","pur80a",no,"bin",ik);   sprintf(title, "mes plab all neu all cuts");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    }
    #endif*/


    for (le = 0; le < 5; ++le) {
      if (le == 0) sprintf(lc, "a");
      if (le == 1) sprintf(lc, "e");
      if (le == 2) sprintf(lc, "m");
      if (le == 3) sprintf(lc, "bch");
      if (le == 4) sprintf(lc, "bnu");
      
      // -- mes histograms for cut studies
      no = 1; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      no = 2; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      no = 3; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      no = 4; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      no = 5; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      no = 6; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      no = 7; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      no = 8; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      
      no =11; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      no =12; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      no =13; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      no =14; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      no =15; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      no =16; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      no =17; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      no =18; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();

      no =21; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      no =22; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      no =23; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      no =24; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      no =25; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      no =26; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      no =27; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      no =28; sprintf(name,"mes%s%d",lc,no); sprintf(title, "mes %d", no); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
      

      /*#ifndef FAST
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

      no = 1500;           
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "theta lep");      h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "theta lep Vcb");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "theta lep Vub");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "theta lep D    ");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "theta lep D*   ");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "theta lep D(*)x");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "theta lep other");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
						    
      no = 1510;           			    
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "theta lep");      h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "theta lep Vcb");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "theta lep Vub");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "theta lep D    ");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "theta lep D*   ");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "theta lep D(*)x");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "theta lep other");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
						    
      no = 1520;           			    
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "theta lep");      h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "theta lep Vcb");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "theta lep Vub");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "theta lep D    ");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "theta lep D*   ");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "theta lep D(*)x");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "theta lep other");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();


      no = 1600;		     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "plab all had");      h = new TH1D(name, title, 2*PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "plab all had Vcb");  h = new TH1D(name, title, 2*PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "plab all had Vub");  h = new TH1D(name, title, 2*PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "plab all had D    ");  h = new TH1D(name, title, 2*PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "plab all had D*   ");  h = new TH1D(name, title, 2*PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "plab all had D(*)x");  h = new TH1D(name, title, 2*PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "plab all had other");  h = new TH1D(name, title, 2*PCBIN, 0., PCMAX); h->Sumw2();
				     
      no = 1610;		     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "plab all had");      h = new TH1D(name, title, 2*PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "plab all had Vcb");  h = new TH1D(name, title, 2*PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "plab all had Vub");  h = new TH1D(name, title, 2*PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "plab all had D");  h = new TH1D(name, title, 2*PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "plab all had D*");  h = new TH1D(name, title, 2*PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "plab all had D(*)x");  h = new TH1D(name, title, 2*PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "plab all had other");  h = new TH1D(name, title, 2*PCBIN, 0., PCMAX); h->Sumw2();

      no = 1620;		     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "plab all had");      h = new TH1D(name, title, 2*PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "plab all had Vcb");  h = new TH1D(name, title, 2*PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "plab all had Vub");  h = new TH1D(name, title, 2*PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "plab all had D");  h = new TH1D(name, title, 2*PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "plab all had D*");  h = new TH1D(name, title, 2*PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "plab all had D(*)x");  h = new TH1D(name, title, 2*PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "plab all had other");  h = new TH1D(name, title, 2*PCBIN, 0., PCMAX); h->Sumw2();


      no = 1700;           
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "theta had");      h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "theta had Vcb");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "theta had Vub");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "theta had D    ");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "theta had D*   ");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "theta had D(*)x");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "theta had other");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
						    
      no = 1710;           			    
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "theta had");      h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "theta had Vcb");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "theta had Vub");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "theta had D    ");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "theta had D*   ");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "theta had D(*)x");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "theta had other");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
						    
      no = 1720;           			    
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "theta had");      h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "theta had Vcb");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "theta had Vub");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "theta had D    ");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "theta had D*   ");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "theta had D(*)x");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "theta had other");  h = new TH1D(name, title, 36, 0., 180.); h->Sumw2();

      no = 1800;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "plab all photons");      h = new TH1D(name, title, 4*PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "plab all photons Vcb");  h = new TH1D(name, title, 4*PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "plab all photons Vub");  h = new TH1D(name, title, 4*PCBIN, 0., PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "plab all photons D    ");  h = new TH1D(name, title, 4*PCBIN, 0.,PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "plab all photons D*   ");  h = new TH1D(name, title, 4*PCBIN, 0.,PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "plab all photons D(*)x");  h = new TH1D(name, title, 4*PCBIN, 0.,PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "plab all photons other");  h = new TH1D(name, title, 4*PCBIN, 0.,PCMAX); h->Sumw2();
      
      no = 1810;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "plab all photons");      h = new TH1D(name, title, 4*PCBIN, 0.,PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "plab all photons Vcb");  h = new TH1D(name, title, 4*PCBIN, 0.,PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "plab all photons Vub");  h = new TH1D(name, title, 4*PCBIN, 0.,PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "plab all photons D");  h = new TH1D(name, title, 4*PCBIN, 0.,PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "plab all photons D*");  h = new TH1D(name, title, 4*PCBIN, 0.,PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "plab all photons D(*)x");  h = new TH1D(name, title, 4*PCBIN, 0.,PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "plab all photons other");  h = new TH1D(name, title, 4*PCBIN, 0.,PCMAX); h->Sumw2();

      no = 1820;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "plab all photons");      h = new TH1D(name, title, 4*PCBIN, 0.,PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "plab all photons Vcb");  h = new TH1D(name, title, 4*PCBIN, 0.,PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "plab all photons Vub");  h = new TH1D(name, title, 4*PCBIN, 0.,PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "plab all photons D");  h = new TH1D(name, title, 4*PCBIN, 0.,PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "plab all photons D*");  h = new TH1D(name, title, 4*PCBIN, 0.,PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "plab all photons D(*)x");  h = new TH1D(name, title, 4*PCBIN, 0.,PCMAX); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "plab all photons other");  h = new TH1D(name, title, 4*PCBIN, 0.,PCMAX); h->Sumw2();

      no = 1900;           
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "theta photons");      h = new TH1D(name, title, 60, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "theta photons Vcb");  h = new TH1D(name, title, 60, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "theta photons Vub");  h = new TH1D(name, title, 60, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "theta photons D    ");  h = new TH1D(name, title, 60, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "theta photons D*   ");  h = new TH1D(name, title, 60, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "theta photons D(*)x");  h = new TH1D(name, title, 60, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "theta photons other");  h = new TH1D(name, title, 60, 0., 180.); h->Sumw2();
                                                  
      no = 1910;                                  
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "theta photons");      h = new TH1D(name, title, 60, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "theta photons Vcb");  h = new TH1D(name, title, 60, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "theta photons Vub");  h = new TH1D(name, title, 60, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "theta photons D    ");  h = new TH1D(name, title, 60, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "theta photons D*   ");  h = new TH1D(name, title, 60, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "theta photons D(*)x");  h = new TH1D(name, title, 60, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "theta photons other");  h = new TH1D(name, title, 60, 0., 180.); h->Sumw2();
                                                  
      no = 1920;                                  
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "theta photons");      h = new TH1D(name, title, 60, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "theta photons Vcb");  h = new TH1D(name, title, 60, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "theta photons Vub");  h = new TH1D(name, title, 60, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "theta photons D    ");  h = new TH1D(name, title, 60, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "theta photons D*   ");  h = new TH1D(name, title, 60, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "theta photons D(*)x");  h = new TH1D(name, title, 60, 0., 180.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "theta photons other");  h = new TH1D(name, title, 60, 0., 180.); h->Sumw2();


				     
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

      
      no = 2600;           	     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "cos#theta pmiss");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "cos#theta pmiss");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "cos#theta pmiss");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "cos#theta pmiss");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "cos#theta pmiss");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "cos#theta pmiss");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "cos#theta pmiss");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
				     
      no = 2610;           	     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "cos#theta pmiss");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "cos#theta pmiss");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "cos#theta pmiss");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "cos#theta pmiss");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "cos#theta pmiss");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "cos#theta pmiss");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "cos#theta pmiss");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();

      no = 2620;           	     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "cos#theta pmiss");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "cos#theta pmiss");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "cos#theta pmiss");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "cos#theta pmiss");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "cos#theta pmiss");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "cos#theta pmiss");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "cos#theta pmiss");  h = new TH1D(name, title, 20, -1., 1.); h->Sumw2();

      no = 2700;           	     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();

      no = 2710;           	     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();

      no = 2720;           	     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "pmiss");  h = new TH1D(name, title, 25, 0., 5.); h->Sumw2();

      
      Float_t MXCUT = 1.6;
      Int_t nbins = 13;
      Int_t nbins2 = 10;
      Float_t ranges[14] = {0.,.3,.6,1.1,MXCUT,1.9,2.2,2.5,2.8,3.1,3.4,3.7,4.2,5.};
      Float_t ranges2[11] = {0.,MXCUT,1.9,2.2,2.5,2.8,3.1,3.4,3.7,4.2,5.};

      no = 2800;		     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "xhad mass");     h = new TH1D(name, title, nbins, ranges); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "xhad mass Vcb"); h = new TH1D(name, title, nbins, ranges); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "xhad mass Vub"); h = new TH1D(name, title, nbins, ranges); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "xhad mass D    "); h = new TH1D(name, title, nbins, ranges); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "xhad mass D*   "); h = new TH1D(name, title, nbins, ranges); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "xhad mass D(*)x"); h = new TH1D(name, title, nbins, ranges); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "xhad mass other"); h = new TH1D(name, title, nbins, ranges); h->Sumw2();
      				     
      no = 2810;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "xhad mass");     h = new TH1D(name, title, nbins, ranges); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "xhad mass Vcb"); h = new TH1D(name, title, nbins, ranges); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "xhad mass Vub"); h = new TH1D(name, title, nbins, ranges); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "xhad mass D    "); h = new TH1D(name, title, nbins, ranges); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "xhad mass D*   "); h = new TH1D(name, title, nbins, ranges); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "xhad mass D(*)x"); h = new TH1D(name, title, nbins, ranges); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "xhad mass other"); h = new TH1D(name, title, nbins, ranges); h->Sumw2();

      no = 2820;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "xhad mass");     h = new TH1D(name, title, nbins, ranges); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "xhad mass Vcb"); h = new TH1D(name, title, nbins, ranges); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "xhad mass Vub"); h = new TH1D(name, title, nbins, ranges); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "xhad mass D    "); h = new TH1D(name, title, nbins, ranges); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "xhad mass D*   "); h = new TH1D(name, title, nbins, ranges); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "xhad mass D(*)x"); h = new TH1D(name, title, nbins, ranges); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "xhad mass other"); h = new TH1D(name, title, nbins, ranges); h->Sumw2();

      no = 2900;		     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "xhad mass");     h = new TH1D(name, title, nbins2, ranges2); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "xhad mass Vcb"); h = new TH1D(name, title, nbins2, ranges2); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "xhad mass Vub"); h = new TH1D(name, title, nbins2, ranges2); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "xhad mass D    "); h = new TH1D(name, title, nbins2, ranges2); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "xhad mass D*   "); h = new TH1D(name, title, nbins2, ranges2); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "xhad mass D(*)x"); h = new TH1D(name, title, nbins2, ranges2); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "xhad mass other"); h = new TH1D(name, title, nbins2, ranges2); h->Sumw2();
      				     
      no = 2910;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "xhad mass");     h = new TH1D(name, title, nbins2, ranges2); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "xhad mass Vcb"); h = new TH1D(name, title, nbins2, ranges2); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "xhad mass Vub"); h = new TH1D(name, title, nbins2, ranges2); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "xhad mass D    "); h = new TH1D(name, title, nbins2, ranges2); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "xhad mass D*   "); h = new TH1D(name, title, nbins2, ranges2); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "xhad mass D(*)x"); h = new TH1D(name, title, nbins2, ranges2); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "xhad mass other"); h = new TH1D(name, title, nbins2, ranges2); h->Sumw2();

      no = 2920;                   
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "xhad mass");     h = new TH1D(name, title, nbins2, ranges2); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "xhad mass Vcb"); h = new TH1D(name, title, nbins2, ranges2); h->Sumw2();
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "xhad mass Vub"); h = new TH1D(name, title, nbins2, ranges2); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "xhad mass D    "); h = new TH1D(name, title, nbins2, ranges2); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "xhad mass D*   "); h = new TH1D(name, title, nbins2, ranges2); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "xhad mass D(*)x"); h = new TH1D(name, title, nbins2, ranges2); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "xhad mass other"); h = new TH1D(name, title, nbins2, ranges2); h->Sumw2();
				     



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

      no = 4400;           	     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "number of kaons");      h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "number of kaons Vcb");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2(); 
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "number of kaons Vub");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "number of kaons D    ");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "number of kaons D*   ");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "number of kaons D(*)x");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "number of kaons other");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
				     			    
      no = 4410;           	     			    
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "number of kaons");      h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "number of kaons Vcb");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2(); 
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "number of kaons Vub");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "number of kaons D    ");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "number of kaons D*   ");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "number of kaons D(*)x");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "number of kaons other");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
							    
      no = 4420;           	     			    
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "number of kaons");      h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "number of kaons Vcb");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2(); 
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "number of kaons Vub");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "number of kaons D    ");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "number of kaons D*   ");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "number of kaons D(*)x");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "number of kaons other");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();

      no = 4500;           	     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "number of kshorts");      h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "number of kshorts Vcb");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2(); 
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "number of kshorts Vub");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "number of kshorts D    ");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "number of kshorts D*   ");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "number of kshorts D(*)x");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "number of kshorts other");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
				     			    
      no = 4510;           	     			    
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "number of kshorts");      h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "number of kshorts Vcb");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2(); 
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "number of kshorts Vub");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "number of kshorts D    ");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "number of kshorts D*   ");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "number of kshorts D(*)x");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "number of kshorts other");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
							    
      no = 4520;           	     			    
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "number of kshorts");      h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "number of kshorts Vcb");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2(); 
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "number of kshorts Vub");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "number of kshorts D    ");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "number of kshorts D*   ");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "number of kshorts D(*)x");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "number of kshorts other");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();

      no = 4600;           	     
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "number of leptons");      h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "number of leptons Vcb");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2(); 
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "number of leptons Vub");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "number of leptons D    ");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "number of leptons D*   ");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "number of leptons D(*)x");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "number of leptons other");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
				     			    
      no = 4610;           	     			    
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "number of leptons");      h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "number of leptons Vcb");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2(); 
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "number of leptons Vub");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "number of leptons D    ");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "number of leptons D*   ");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "number of leptons D(*)x");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "number of leptons other");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
							    
      no = 4620;           	     			    
      sprintf(name,"%s%d",lc,no+0); sprintf(title, "number of leptons");      h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+1); sprintf(title, "number of leptons Vcb");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2(); 
      sprintf(name,"%s%d",lc,no+2); sprintf(title, "number of leptons Vub");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+3); sprintf(title, "number of leptons D    ");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+4); sprintf(title, "number of leptons D*   ");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+5); sprintf(title, "number of leptons D(*)x");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();
      sprintf(name,"%s%d",lc,no+6); sprintf(title, "number of leptons other");  h = new TH1D(name, title, 10, 0., 10.); h->Sumw2();

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
				     
      //#endif
      */
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
}

// ----------------------------------------------------------------------
void VubAnalysisCode::maskKshorts(int modes) {
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
    sprintf(name, "r100");  sprintf(title, "Kshort mass with r < x cm");  h = new TH1D(name, title, 100, 0.46, 0.54); 
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
	    if ((Brectrktmp&brecoOverlap))  ++cntBreco;
	    if (!(Brectrktmp&brecoOverlap))  ++cntRecoil;
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
    goodKshort[i] = 0;
    goodWe[i] = 0;
    goodWk[i] = 0;
    goodNr[i] = 0;
  }
  if (fVerbose) cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  for (i = 0; i < nKs; ++i) {
    if (TMath::Abs(d1KsLund[i]) != 211) continue;
    if (TMath::Abs(d2KsLund[i]) != 211) continue;
    if (fVerbose) cout << "Ks[" << i << "] -> (" << d1KsIndex[i]-1 << "," << d2KsIndex[i]-1 << ") " 
		  << "with daughters " << d1KsLund[i] << " " << d2KsLund[i]
		  << " and mass = " << massKs[i] << endl;
    if ((d1KsIndex[i]-1 > nTrk) || (d2KsIndex[i]-1 > nTrk)) { 
      cout << "Daughter indices messed up: d1KsIndex = " << d1KsIndex[i]-1 << " d2KsIndex = " << d2KsIndex[i]-1 
	   << "  for nTrk = " << nTrk << endl;
      continue;
    }
    partOfKs[d1KsIndex[i]-1] += 1; 
    partOfKs[d2KsIndex[i]-1] += 1; 
  }
  
  // =============
  // -- KS->pi+pi-
  // =============
  // -- Two passes over KS block: (1) Reject trivial cases: Ks with kaons and electron, overlaps with BRECO
  for (i = 0; i < nKs; ++i) {
    if (TMath::Abs(d1KsLund[i]) != 211) {
      goodKshort[i] = 0;
      continue;
    }
    if (TMath::Abs(d2KsLund[i]) != 211) {
      goodKshort[i] = 0;
      continue;
    }
    int pi1 = d1KsIndex[i]-1;
    int pi2 = d2KsIndex[i]-1; 
    double mass = massKs[i];
    int mc = MCKs[i]-1; 
    int matched(0);
    if ((mc > -1) && (mc < nMc)) {
      if (TMath::Abs(idMc[mc]) == 310) matched = 1;
    }

    if (isRecMu(pi1) || isRecMu(pi2)) ((TH1D*)gDirectory->Get("m100"))->Fill(mass);	
    // -- Veto KS candidates with electrons or kaons (note: KS candidates with muons peak at m(KS))
    if (isRecEl(pi1)) {
      if (fVerbose) cout << "Rejecting Ks[" << i <<"] due to identified electron " << pi1 << endl;
      //      if (matched == 1) cout << "Rejecting matched KS with el1 cut" << endl;
      goodKshort[i] = 0;
      goodWe[i] += 1;
      ((TH1D*)gDirectory->Get("e100"))->Fill(mass);	
      continue;
    }
    if (isRecEl(pi2)) {
      if (fVerbose) cout << "Rejecting Ks[" << i <<"] due to identified electron " << pi2 << endl;
      //      if (matched == 1) cout << "Rejecting matched KS with el2 cut" << endl;
      goodKshort[i] = 0;
      goodWe[i] += 1;
      ((TH1D*)gDirectory->Get("e100"))->Fill(mass);	
      continue;
    }
    if (isRecKaon(pi1)) {
      if (fVerbose) cout << "Rejecting Ks[" << i <<"] due to identified kaon " << pi1 << endl;
      //      if (matched == 1) cout << "Rejecting matched KS with ka1 cut" << endl;
      goodKshort[i] = 0;
      goodWk[i] += 1;
      ((TH1D*)gDirectory->Get("k100"))->Fill(mass);	
      continue;
    }
    if (isRecKaon(pi2)) {
      if (fVerbose) cout << "Rejecting Ks[" << i <<"] due to identified kaon " << pi2 << endl;
      //      if (matched == 1) cout << "Rejecting matched KS with ka2 cut" << endl;
      goodKshort[i] = 0;
      goodWk[i] += 1;
      ((TH1D*)gDirectory->Get("k100"))->Fill(mass);	
      continue;
    }
    // -- Check against overlap with BRECO candidate
    int overlap(0);
    Brectrktmp = B0RecTrk[pi1];
    if(fBrecoCharge != 0) Brectrktmp = chBRecTrk[pi1];
    if ((Brectrktmp&brecoOverlap)) {
      if(fVerbose) cout << "Rejecting Ks[" << i << "] due to pi1 at " << pi1 << " overlaps with BRECO " << endl;
      overlap++;
    }
    Brectrktmp = B0RecTrk[pi2];
    if(fBrecoCharge != 0) Brectrktmp = chBRecTrk[pi2];
    if ((Brectrktmp&brecoOverlap)) {
      if(fVerbose) cout << "Rejecting Ks[" << i << "] due to pi2 at " << pi2 << " overlaps with BRECO " << endl;
      overlap++;
    }
    if (overlap > 0) {
      goodKshort[i] = 0;
      goodNr[i] = overlap;
      //      if (overlap == 1 && matched == 1) cout << "Rejecting matched KS with 1 pi breco-overlap cut" << endl;
      if (overlap == 1) ((TH1D*)gDirectory->Get("h101"))->Fill(mass);
      if (overlap == 2) ((TH1D*)gDirectory->Get("h102"))->Fill(mass);
      continue;
    } else {
      ((TH1D*)gDirectory->Get("h103"))->Fill(mass);
    }
    // -- Check for vertex separation from origin
    double rks = TMath::Sqrt( (xKs[i]-beamSX)*(xKs[i]-beamSX) + (yKs[i]-beamSY)*(yKs[i]-beamSY) );
    if (rks < KSPIPRLO) {
      goodKshort[i] = 0;
      ((TH1D*)gDirectory->Get("r100"))->Fill(mass);	
      continue;
    }

    ((TH1D*)gDirectory->Get("h100"))->Fill(mass);
    if ((mass < KSPIPLO) || (mass > KSPIPHI)) {
      //      if (matched == 1) cout << "Rejecting matched KS with mass cut: mass = " << mass << endl;
      goodKshort[i] = 2;
      continue;
    }
      
    goodKshort[i] = 1;
    kshortLockTrk[pi1] = 1;
    kshortLockTrk[pi2] = 1;      

  } // first pass over KS block

  // -- Second Pass: Reject Ks whose daughters are part of a 'better' KS
  //    for (i = 0; i < nKs; ++i) {
  //      if (goodKshort[i] == 0) continue;
  //      if (TMath::Abs(d1KsLund[i]) != 211) continue;
  //      if (TMath::Abs(d2KsLund[i]) != 211) continue;
  //      int takeThis(1); 
  //      int pi1 = d1KsIndex[i]-1;
  //      int pi2 = d2KsIndex[i]-1; 
  //      double mass = massKs[i];
  //      double residual = TMath::Abs(KAZMASS - mass); 
  //      int mc = MCKs[i]-1; 
  //      int matched(0);
  //      if ((mc > -1) && (mc < nMc)) {
  //        if (TMath::Abs(idMc[mc]) == 310) matched = 1;
  //      }
  //      if (partOfKs[pi1] > 1) {
  //        for (j = 0; j < nKs; ++j) {
  //  	if (i == j) continue;
  //  	if (goodKshort[j] == 0) continue;
  //  	// -- Compare residuals
  //  	if ((d1KsIndex[j]-1 == pi1) || (d2KsIndex[j]-1 == pi1)) {
  //  	  if (TMath::Abs(KAZMASS - massKs[j]) < residual) {
  //  	    if (matched == 1) cout << "Rejecting matched KS with best-selection cut" << endl;
  //  	    if (fVerbose) cout << "Rejecting KS[" << i << "] -> (" << pi1 << "," << pi2 << ")"
  //  			  << " with mass = " << mass 
  //  			  << ", since KS[" << j << "] -> (" << d1KsIndex[j]-1 << "," << d2KsIndex[j]-1 << ")"
  //  			  << " has better mass: " << massKs[j] << endl;
  //  	    takeThis = 0; 
  //  	    goodKshort[i] = 0;
  //  	    break;
  //  	  }
  //  	}
  //        }
  //      }
  //      if (partOfKs[pi2] > 1) {
  //        for (j = 0; j < nKs; ++j) {
  //  	if (i == j) continue;
  //  	// -- do not compare against KS overlapping with BRECO
  //  	if (goodKshort[j] == 0) continue;
  //  	// -- Compare residuals
  //  	if ((d1KsIndex[j]-1 == pi2) || (d2KsIndex[j]-1 == pi2)) {
  //  	  if (TMath::Abs(KAZMASS - massKs[j]) < residual) {
  //  	    if (fVerbose) cout << "Rejecting KS[" << i << "] -> (" << pi1 << "," << pi2 << ")"
  //  			  << " mass = " << mass 
  //  			  << ", since KS[" << j << "] -> (" << d1KsIndex[j]-1 << "," << d2KsIndex[j]-1 << ")"
  //  			  << " has better mass: " << massKs[j] << endl;
  //  	    takeThis = 0; 
  //  	    goodKshort[i] = 0;
  //  	    break;
  //  	  }
  //  	}
  //        }
  //      }
  //      if (takeThis == 1) {
  //        ((TH1D*)gDirectory->Get("h104"))->Fill(mass);
  //        if (signalBox) ((TH1D*)gDirectory->Get("h105"))->Fill(mass);
  //        goodKshort[i] = 1;
  //        kshortLockTrk[pi1] = 1;
  //        kshortLockTrk[pi2] = 1;      
  //        if (fIsMC && (TMath::Abs(idMc[MCKs[i]-1])==310)) ((TH1D*)gDirectory->Get("h106"))->Fill(mass);
  //        if(fVerbose) cout << "Locking KS[" << i << "] -> ("  << pi1 << "," << pi2 << ") mass = " << mass << endl;
  //        ++goodKcharged;
  //      }
  //    } // second pass over KS block


  return;

  // =============
  // -- KS->pi0pi0
  // =============
  for (i = 0; i < nKs; ++i) {
    if (TMath::Abs(d1KsLund[i]) != 111) continue;
    if (TMath::Abs(d2KsLund[i]) != 111) continue;

    goodKshort[i] = 0;

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
void  VubAnalysisCode::mxCategory() {

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
  Bool_t cascade(kFALSE), klong(kFALSE), kshort(kFALSE), kspippim(kFALSE), kspi0pi0(kFALSE), neutron(kFALSE), lostk(kFALSE), kplus(kFALSE), kmiss(kFALSE) ;
  
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
     // -- charged Kaons
    if ((overlap == 2) && (TMath::Abs(idMc[i]) == 321)) {
      if (isRecoed(i) == -1) {
        kmiss= kTRUE;
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
      if (nDauMc[i] == 2) {
	Int_t id[2]={0,0};
	Int_t counter=0;
	for (Int_t j=0; j<nMc; j++) {
	  Int_t imother = mothMc[j]-1;
	  if (imother == i) {
	    id[counter]=idMc[j];
	    counter++;
	  }
	  if (counter == 2) break;
	}
	
	if ((TMath::Abs(id[0])==211) 
	    && (TMath::Abs(id[1])==211)) 
	  kspippim=kTRUE;
	if ((id[0]==111) && (id[1]==111)) 
	  kspi0pi0=kTRUE;
      }
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
  if (klong)    fMxCategory += MCklong;    // 1
  if (kspippim) fMxCategory += MCkshortpippim;   // 2
  if (kspi0pi0) fMxCategory += MCkshortpi0pi0;   // 4
  if (kplus)    fMxCategory += MCkplus;    // 8
  if (lostk)    fMxCategory += MCkineff;    // 16
  if (kmiss)    fMxCategory += MCkmiss;    // 32
  if (cascade)  fMxCategory += MCcascade;  // 64
    

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
void VubAnalysisCode::Loop(Int_t maxEvent, Int_t startEvent, Int_t isVerbose, Int_t lun) {

  int step(1000);
  fChB=0;
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
  int nvxbevt(0); 
  const char *pChar; 

  int oldrunnumber(0); 
  TString rfile(fHistFile->GetName());
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

    fChB = 0;
    if (bestB0 == 0) {
      fChB = 1;
    }

    // -- MonteCarlo Truth
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

    //      if (fBVxb > 0) {
    //        ++nvxbevt; 
    //        cout << "beginevent: " << nvxbevt << endl;
    //        dumpOneB(fBVxb); 
    //        cout << "endevent: " << nvxbevt << endl;
    //      }

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
    //    trackStudy(); 

    // splitoff studies --- needed for photon selection
    doSplitOffStudy();
    selectPhotons();

    TLorentzVector p4t = p4Breco; 
    p4t.Boost(-upsBoost);
    fPcmsBreco = p4t.Vect().Mag(); 


    if (fVerbose) cout << "CALLING breco()" << endl;
    breco();
    if (fOptCategories) mxCategory();
    maskKshorts(1);
    //     maskPi0(1);
    //    maskConversions();
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
void VubAnalysisCode::maskPi0(int modes) {
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
void VubAnalysisCode::maskConversions(int modes)  {
  Int_t i(0),ip1,ip2;
  Double_t mass;
  
  for (i = 0; i < 100; ++i) {
    goodConv[i] = 0;
    convLockTrk[i] = 0;
  }
  
  for (i = 0; i < nGConv; ++i) {
    mass = massGConv[i];
    ip1 = d1GConvIndex[i]-1; 
    ip2 = d2GConvIndex[i]-1; 
    if ((isRecEl(ip1) || isRecEl(ip2)) && (mass<0.010)) {
      goodConv[i] = 1;
      convLockTrk[ip1] = 1;
      convLockTrk[ip2] = 1;
      if (momentumTrk[ip1] > 1.0 || momentumTrk[ip2] > 1.0) {
        if (fVerbose) cout << " -> Masking conversion with daugher track lab momenta: p1 = " 
                      << momentumTrk[ip1] << " and p2 = " << momentumTrk[ip2] << endl;
      }
    }
  }
}


// ----------------------------------------------------------------------
void VubAnalysisCode::fastBookHist(const char *name) {

  fHistFile->cd(); 
  fHistFile->mkdir(name, name); 
  fHistFile->cd(name);

  int i(-1); 
  int nbins(1);
  int bins[20]; 

  Bool_t *nocuts, *sigcuts, *aocuts; 
  Bool_t *aoprmm2, *aomm2, *aocc; 

  //  allcuts = &fGoodEvent; 
  sigcuts = &fLowMx; 
      
  if      (!strcmp(name, "a_dep"))  {
    i =  1; 
    nocuts = &fGoodAccLepton; 
    aocuts = &faoLepton; 
    aoprmm2= &faoPRMM2;
    aomm2  = &faoMM2;
    aocc   = &faoChargeCons;
  }  else if (!strcmp(name, "a_enh"))  {
    i =  2; 
    nocuts = &fGoodAccLepton; 
    aocuts = &faoLepton; 
    aoprmm2= &faoPRMM2;
    aomm2  = &faoMM2;
    aocc   = &faoChargeCons;
  }  else if (!strcmp(name, "e_dep"))  {
    i =  3; 
    nocuts = &fGoodAccElectron; 
    aocuts = &faoElectron; 
    aoprmm2= &faoPRMM2E;
    aomm2  = &faoMM2E;
    aocc   = &faoChargeConsE;
  }  else if (!strcmp(name, "e_enh"))  {
    i =  4; 
    nocuts = &fGoodAccElectron; 
    aocuts = &faoElectron; 
    aoprmm2= &faoPRMM2E;
    aomm2  = &faoMM2E;
    aocc   = &faoChargeConsE;
  }  else if (!strcmp(name, "m_dep"))  {
    i =  5; 
    nocuts = &fGoodAccMuon; 
    aocuts = &faoMuon; 
    aoprmm2= &faoPRMM2M;
    aomm2  = &faoMM2M;
    aocc   = &faoChargeConsM;
  }  else if (!strcmp(name, "m_enh"))  {
    i =  6; 
    nocuts = &fGoodAccMuon; 
    aocuts = &faoMuon; 
    aoprmm2= &faoPRMM2M;
    aomm2  = &faoMM2M;
    aocc   = &faoChargeConsM;

  }  else if (!strcmp(name, "ars_dep"))  {
    i = 10; 
    nocuts = &fGoodAccLepton; 
    aocuts = &faoRS; 
    aoprmm2= &faoPRMM2RS;
    aomm2  = &faoMM2RS;
    aocc   = &faoChargeConsRS;
  }  else if (!strcmp(name, "ars_enh"))  {
    i = 11; 
    nocuts = &fGoodAccLepton; 
    aocuts = &faoRS; 
    aoprmm2= &faoPRMM2RS;
    aomm2  = &faoMM2RS;
    aocc   = &faoChargeConsRS;
  }  else if (!strcmp(name, "aws_dep"))  {
    i = 12; 
    nocuts = &fGoodAccLepton; 
    aocuts = &faoWS; 
    aoprmm2= &faoPRMM2WS;
    aomm2  = &faoMM2WS;
    aocc   = &faoChargeConsWS;
  }  else if (!strcmp(name, "aws_enh"))  {
    i = 13; 
    nocuts = &fGoodAccLepton; 
    aocuts = &faoWS; 
    aoprmm2= &faoPRMM2WS;
    aomm2  = &faoMM2WS;
    aocc   = &faoChargeConsWS;

  }  else if (!strcmp(name, "arf_dep"))  {
    i = 14; 
    nocuts = &fGoodAccLepton; 
    aocuts = &faoRF; 
    aoprmm2= &faoPRMM2RF;
    aomm2  = &faoMM2RF;
    aocc   = &faoChargeConsRF;
  }  else if (!strcmp(name, "arf_enh"))  {
    i = 15; 
    nocuts = &fGoodAccLepton; 
    aocuts = &faoRF; 
    aoprmm2= &faoPRMM2RF;
    aomm2  = &faoMM2RF;
    aocc   = &faoChargeConsRF;
  }  else if (!strcmp(name, "awf_dep"))  {
    i = 16; 
    nocuts = &fGoodAccLepton; 
    aocuts = &faoWF; 
    aoprmm2= &faoPRMM2WF;
    aomm2  = &faoMM2WF;
    aocc   = &faoChargeConsWF;
  }  else if (!strcmp(name, "awf_enh"))  {
    i = 17;
    nocuts = &fGoodAccLepton; 
    aocuts = &faoWF; 
    aoprmm2= &faoPRMM2WF;
    aomm2  = &faoMM2WF;
    aocc   = &faoChargeConsWF;

  }  else if (!strcmp(name, "ers_dep"))  {
    i = 20; 
    nocuts = &fGoodAccElectron; 
    aocuts = &faoERS; 
    aoprmm2= &faoPRMM2ERS;
    aomm2  = &faoMM2ERS;
    aocc   = &faoChargeConsERS;
  }  else if (!strcmp(name, "ers_enh"))  {
    i = 21; 
    nocuts = &fGoodAccElectron; 
    aocuts = &faoERS; 
    aoprmm2= &faoPRMM2ERS;
    aomm2  = &faoMM2ERS;
    aocc   = &faoChargeConsERS;
  }  else if (!strcmp(name, "ews_dep"))  {
    i = 22; 
    nocuts = &fGoodAccElectron; 
    aocuts = &faoEWS; 
    aoprmm2= &faoPRMM2EWS;
    aomm2  = &faoMM2EWS;
    aocc   = &faoChargeConsEWS;
  }  else if (!strcmp(name, "ews_enh"))  {
    i = 23; 
    nocuts = &fGoodAccElectron; 
    aocuts = &faoEWS; 
    aoprmm2= &faoPRMM2EWS;
    aomm2  = &faoMM2EWS;
    aocc   = &faoChargeConsEWS;

  }  else if (!strcmp(name, "erf_dep"))  {
    i = 24; 
    nocuts = &fGoodAccElectron; 
    aocuts = &faoERF; 
    aoprmm2= &faoPRMM2ERF;
    aomm2  = &faoMM2ERF;
    aocc   = &faoChargeConsERF;
  }  else if (!strcmp(name, "erf_enh"))  {
    i = 25; 
    nocuts = &fGoodAccElectron; 
    aocuts = &faoERF; 
    aoprmm2= &faoPRMM2ERF;
    aomm2  = &faoMM2ERF;
    aocc   = &faoChargeConsERF;
  }  else if (!strcmp(name, "ewf_dep"))  {
    i = 26; 
    nocuts = &fGoodAccElectron; 
    aocuts = &faoEWF; 
    aoprmm2= &faoPRMM2EWF;
    aomm2  = &faoMM2EWF;
    aocc   = &faoChargeConsEWF;
  }  else if (!strcmp(name, "ewf_enh"))  {
    i = 27;
    nocuts = &fGoodAccElectron; 
    aocuts = &faoEWF; 
    aoprmm2= &faoPRMM2EWF;
    aomm2  = &faoMM2EWF;
    aocc   = &faoChargeConsEWF;

  }  else if (!strcmp(name, "mrs_dep"))  {
    i = 30; 
    nocuts = &fGoodAccMuon; 
    aocuts = &faoMRS; 
    aoprmm2= &faoPRMM2MRS;
    aomm2  = &faoMM2MRS;
    aocc   = &faoChargeConsMRS;
  }  else if (!strcmp(name, "mrs_enh"))  {
    i = 31; 
    nocuts = &fGoodAccMuon; 
    aocuts = &faoMRS; 
    aoprmm2= &faoPRMM2MRS;
    aomm2  = &faoMM2MRS;
    aocc   = &faoChargeConsMRS;
  }  else if (!strcmp(name, "mws_dep"))  {
    i = 32; 
    nocuts = &fGoodAccMuon; 
    aocuts = &faoMWS; 
    aoprmm2= &faoPRMM2MWS;
    aomm2  = &faoMM2MWS;
    aocc   = &faoChargeConsMWS;
  }  else if (!strcmp(name, "mws_enh"))  {
    i = 33; 
    nocuts = &fGoodAccMuon; 
    aocuts = &faoMWS; 
    aoprmm2= &faoPRMM2MWS;
    aomm2  = &faoMM2MWS;
    aocc   = &faoChargeConsMWS;

  }  else if (!strcmp(name, "mrf_dep"))  {
    i = 34; 
    nocuts = &fGoodAccMuon; 
    aocuts = &faoMRF; 
    aoprmm2= &faoPRMM2MRF;
    aomm2  = &faoMM2MRF;
    aocc   = &faoChargeConsMRF;
  }  else if (!strcmp(name, "mrf_enh"))  {
    i = 35; 
    nocuts = &fGoodAccMuon; 
    aocuts = &faoMRF; 
    aoprmm2= &faoPRMM2MRF;
    aomm2  = &faoMM2MRF;
    aocc   = &faoChargeConsMRF;
  }  else if (!strcmp(name, "mwf_dep"))  {
    i = 36; 
    nocuts = &fGoodAccMuon; 
    aocuts = &faoMWF; 
    aoprmm2= &faoPRMM2MWF;
    aomm2  = &faoMM2MWF;
    aocc   = &faoChargeConsMWF;
  }  else if (!strcmp(name, "mwf_enh"))  {
    i = 37;
    nocuts = &fGoodAccMuon; 
    aocuts = &faoMWF; 
    aoprmm2= &faoPRMM2MWF;
    aomm2  = &faoMM2MWF;
    aocc   = &faoChargeConsMWF;
 
  }  else if (!strcmp(name, "all"))     {
    i = 0; 
    nocuts = &fGoodAccLepton; 
    aocuts = &faoLepton; 
    aoprmm2= &faoPRMM2;
    aomm2  = &faoMM2;
    aocc   = &faoChargeCons;

  }  else {
    cout << "Unknown directory name " << endl;
  }

  // -- Lepton (NOTE: aocuts does NOT include fOneLepton!)
  d1050[i] = new valHist(1050, "pcms ", 18, 0., 3.0);  
  d1050[i]->setup(&fPcms, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, aocuts); 

  d1550[i] = new valHist(1550, "plab ", 20, 0., 4.);  
  d1550[i]->setup(&fPlab, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, aocuts); 

  d1650[i] = new valHist(1650, "theta", 18, 0., 180.);  
  d1650[i]->setup(&fTlabDR, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, aocuts); 

  d1750[i] = new valHist(1750, "phi", 18, -180., 180.);  
  d1750[i]->setup(&fFlabDR, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, aocuts); 


  // -- Mxhad 
  d2050[i] = new valHist(2050, "MxHad ", 15, 0., 3.);  
  d2050[i]->setup(&fMxhad, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, &fGoodEvent); 

  d2150[i] = new valHist(2150, "Q^{2} ", 20, 0., 20.);  
  d2150[i]->setup(&fQ2, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, &fGoodEvent); 

  d2250[i] = new valHist(2250, "#Delta MxHad ", 10, -1., 1.);  
  d2250[i]->setup(&fMxhadRes, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, &fGoodEvent); 

  d2350[i] = new valHist(2350, "#Delta Q^{2} ", 10, -3., 3.);  
  d2350[i]->setup(&fQ2Res, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, &fGoodEvent); 


  // -- Mxhadfit 
  d2550[i] = new valHist(2550, "MxHadFit ", 15, 0., 3.);  
  d2550[i]->setup(&fMxhadfit, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, &fGoodEvent); 
  d2650[i] = new valHist(2650, "#Delta MxHadFit ", 10, -1., 1.);  
  d2650[i]->setup(&fMxhadfitRes, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, &fGoodEvent); 

  d2750[i] = new valHist(2750, "Prob(chi2)", 10, 0., 1.);  
  d2750[i]->setup(&fProbChi2, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, &fGoodEvent); 
  d2850[i] = new valHist(2850, "chi2", 10, 0., 40.);  
  d2850[i]->setup(&fChi2, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, &fGoodEvent); 


  // -- Neutrino
  d3050[i] = new valHist(3050, "m_{miss}^{2} ", 30, -5., 10.);  
  d3050[i]->setup(&fMM2, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, aomm2); 

  d3150[i] = new valHist(3150, "m_{miss}^{2} (partial reco'ed D*)", 30, -50., 10.);  
  d3150[i]->setup(&fwdeltaM, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, aoprmm2); 

  d3550[i] = new valHist(3550, "|pmiss| ", 16, 0., 4.);  
  d3550[i]->setup(&fPNu, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, aomm2); 

  d3650[i] = new valHist(3650, "Emiss ", 16, 0., 4.);  
  d3650[i]->setup(&fEmiss, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, aomm2); 

  d3750[i] = new valHist(3750, "cos #theta_{miss} ", 10, -1., 1.);  
  d3750[i]->setup(&fCosTNu, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, aomm2); 


  // -- Event quantities
  d4050[i] = new valHist(4050, "Q_{total} ", 8, -4., 4.);  
  d4050[i]->setup(&fQtot, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, aocc); 

  d4150[i] = new valHist(4150, "Q_{recoil} ", 8, -4., 4.);  
  d4150[i]->setup(&fQrecoil, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, aocc); 

  d4250[i] = new valHist(4250, "N_{Tracks}", 10, 0., 10.);  
  d4250[i]->setup(&fNtracks, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, &fGoodEvent); 

  d4350[i] = new valHist(4350, "N_{Neutrals} ", 14, 0., 14.);  
  d4350[i]->setup(&fNneutrals, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, &fGoodEvent); 

  d4450[i] = new valHist(4450, "N_{Leptons} ", 4, 0., 4., nbins, bins);  
  d4450[i]->setup(&frealNLepton, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, &fGoodEvent); 

  d4550[i] = new valHist(4550, "N_{Kaons} ", 4, 0., 4., nbins, bins);  
  d4550[i]->setup(&frealNKp, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, &fGoodEvent); 

  d4650[i] = new valHist(4650, "N_{KS} ", 4, 0., 4.);  
  d4650[i]->setup(&frealNKp, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, &fGoodEvent); 


  // -- plab for miscellaneous particles
  d5050[i] = new valHist(5050, "plab_{Tracks} ", 20, 0., 3.0);  
  d5050[i]->setup(&fPlab, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, &fGoodEvent); 

  d5150[i] = new valHist(5150, "Elab_{Neutrals} ", 10, 0., 1.5);  
  d5150[i]->setup(&fPlab, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, &fGoodEvent); 

  d5250[i] = new valHist(5250, "plab_{K+} ", 20, 0., 3.0);  
  d5250[i]->setup(&fPlab, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, &fGoodEvent); 

  d5350[i] = new valHist(5350, "plab_{Ks} ", 20, 0., 3.0);  
  d5350[i]->setup(&fPlab, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, &fGoodEvent); 

  d5450[i] = new valHist(5450, "plab_{Pions} ", 20, 0., 3.0);  
  d5450[i]->setup(&fPlab, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, &fGoodEvent); 


  // -- theta for miscellaneous particles
  d6050[i] = new valHist(6050, "#theta_{tracks}", 12, 0., 3.15);  
  d6050[i]->setup(&fTlabDR, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, &fGoodEvent); 

  d6150[i] = new valHist(6150, "#theta_{neutrals}", 12, 0., 3.15);  
  d6150[i]->setup(&fTlabDR, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, &fGoodEvent); 

  d7000[i] = new valHist(7000, "track length ", 20, 0., 200.0);  
  d7000[i]->setup(&fPlab, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, &fGoodEvent); 

  d7050[i] = new valHist(7050, "chi2/ndof ", 12, 0., 4.0);  
  d7050[i]->setup(&fPlab, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, &fGoodEvent); 

  d7100[i] = new valHist(7100, "log(dca) ", 20, -4., .5);  
  d7100[i]->setup(&fPlab, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, &fGoodEvent); 

  d7150[i] = new valHist(7150, "dca z ", 15, -3., 3.0);  
  d7150[i]->setup(&fPlab, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, &fGoodEvent); 

  d7200[i] = new valHist(7200, "NDch ", 20, 10., 50.0);  
  d7200[i]->setup(&fPlab, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, &fGoodEvent); 

  d7250[i] = new valHist(7250, "NSvt ", 13, 0., 13.0);  
  d7250[i]->setup(&fPlab, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, &fGoodEvent); 

  d8000[i] = new valHist(8000, "LAT ", 20, 0., 1.0);  
  d8000[i]->setup(&fPlab, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, &fGoodEvent); 

  d8050[i] = new valHist(8050, "secmom ", 20, 0., 0.005);  
  d8050[i]->setup(&fPlab, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, &fGoodEvent); 

  d8100[i] = new valHist(8100, "A20 ", 20, 0.5, 1.0);  
  d8100[i]->setup(&fPlab, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, &fGoodEvent); 

  d8150[i] = new valHist(8150, "A42 ", 20, 0., 0.3);  
  d8150[i]->setup(&fPlab, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, &fGoodEvent); 

  d8200[i] = new valHist(8200, "s9s25 ", 20, 0.5, 1.0);  
  d8200[i]->setup(&fPlab, &fMes, &signalBox, &mesSideband, nocuts, sigcuts, &fGoodEvent); 

 
}

// ----------------------------------------------------------------------
void VubAnalysisCode::fastFillHist(const char *name) {
  int i(-1);


  if      (!strcmp(name, "a_dep"))    i =  1; 
  else if (!strcmp(name, "a_enh"))    i =  2; 
  else if (!strcmp(name, "e_dep"))    i =  3; 
  else if (!strcmp(name, "e_enh"))    i =  4; 
  else if (!strcmp(name, "m_dep"))    i =  5; 
  else if (!strcmp(name, "m_enh"))    i =  6; 

  else if (!strcmp(name, "ars_dep"))  i = 10; 
  else if (!strcmp(name, "ars_enh"))  i = 11; 
  else if (!strcmp(name, "aws_dep"))  i = 12; 
  else if (!strcmp(name, "aws_enh"))  i = 13; 
  else if (!strcmp(name, "arf_dep"))  i = 14; 
  else if (!strcmp(name, "arf_enh"))  i = 15; 
  else if (!strcmp(name, "awf_dep"))  i = 16; 
  else if (!strcmp(name, "awf_enh"))  i = 17;

  else if (!strcmp(name, "ers_dep"))  i = 20; 
  else if (!strcmp(name, "ers_enh"))  i = 21; 
  else if (!strcmp(name, "ews_dep"))  i = 22; 
  else if (!strcmp(name, "ews_enh"))  i = 23; 
  else if (!strcmp(name, "erf_dep"))  i = 24; 
  else if (!strcmp(name, "erf_enh"))  i = 25; 
  else if (!strcmp(name, "ewf_dep"))  i = 26; 
  else if (!strcmp(name, "ewf_enh"))  i = 27;

  else if (!strcmp(name, "mrs_dep"))  i = 30; 
  else if (!strcmp(name, "mrs_enh"))  i = 31; 
  else if (!strcmp(name, "mws_dep"))  i = 32; 
  else if (!strcmp(name, "mws_enh"))  i = 33; 
  else if (!strcmp(name, "mrf_dep"))  i = 34; 
  else if (!strcmp(name, "mrf_enh"))  i = 35; 
  else if (!strcmp(name, "mwf_dep"))  i = 36; 
  else if (!strcmp(name, "mwf_enh"))  i = 37;
 
  else if (!strcmp(name, "all"))     i = 0; 
  else {
    cout << "Unknown directory name " << endl;
  }

  d1050[i]->fillHist(fEvtW8);
  d1550[i]->fillHist(fEvtW8); 
  d1650[i]->fillHist(fEvtW8); 
  d1750[i]->fillHist(fEvtW8);

  d2050[i]->fillHist(fEvtW8); 
  d2150[i]->fillHist(fEvtW8);
  d2250[i]->fillHist(fEvtW8);
  d2350[i]->fillHist(fEvtW8);

  d2550[i]->fillHist(fEvtW8);
  d2650[i]->fillHist(fEvtW8);
  d2750[i]->fillHist(fEvtW8);   
  d2850[i]->fillHist(fEvtW8);

  d3050[i]->fillHist(fEvtW8);
  d3150[i]->fillHist(fEvtW8);
  d3550[i]->fillHist(fEvtW8); 
  d3650[i]->fillHist(fEvtW8); 
  d3750[i]->fillHist(fEvtW8);

  d4050[i]->fillHist(fEvtW8); 
  d4150[i]->fillHist(fEvtW8); 
  d4250[i]->fillHist(fEvtW8); 
  d4350[i]->fillHist(fEvtW8);
  d4450[i]->fillHist(fEvtW8); 
  d4550[i]->fillHist(fEvtW8); 
  d4650[i]->fillHist(fEvtW8);

  d5050[i]->fillHist(nTrk, &momentumTrk[0], &goodHadron[0], fEvtW8);
  d5150[i]->fillHist(nGam, &energyGam[0], &goodPhoton[0], fEvtW8);
  d5250[i]->fillHist(nTrk, &momentumTrk[0], &goodChargedKaon[0], fEvtW8);
  d5350[i]->fillHist(nKs, &pKs[0], &goodKshort[0], fEvtW8);
  d5450[i]->fillHist(nTrk, &momentumTrk[0], &goodPion[0], fEvtW8);

  d6050[i]->fillHist(nTrk, &thetaTrk[0], &goodHadron[0], fEvtW8);
  d6150[i]->fillHist(nGam, &thetaGam[0], &goodPhoton[0], fEvtW8);

  d7000[i]->fillHist(nTrk, &tlenTrack[0], &goodTrack[0], fEvtW8);   
  d7050[i]->fillHist(nTrk, &c2nTrack[0],  &goodTrack[0], fEvtW8);   
  d7100[i]->fillHist(nTrk, &dcaTrack[0],  &goodTrack[0], fEvtW8);   
  d7150[i]->fillHist(nTrk, &dcazTrack[0], &goodTrack[0], fEvtW8);   
  d7200[i]->fillHist(nTrk, &ndchTrack[0], &goodTrack[0], fEvtW8);   
  d7250[i]->fillHist(nTrk, &nsvtTrack[0], &goodTrack[0], fEvtW8);   
 
  d8000[i]->fillHist(nTrk, &lMomGam[0],   &goodPhoton[0], fEvtW8);   
  d8050[i]->fillHist(nTrk, &secMomGam[0], &goodPhoton[0], fEvtW8);   
  d8100[i]->fillHist(nTrk, &ZMom20Gam[0], &goodPhoton[0], fEvtW8);   
  d8150[i]->fillHist(nTrk, &ZMom42Gam[0], &goodPhoton[0], fEvtW8);   
  d8200[i]->fillHist(nTrk, &s9s25Gam[0],  &goodPhoton[0], fEvtW8);   


}


// ----------------------------------------------------------------------
void VubAnalysisCode::fillMesHist(const char *dir, const char *le) {
  fHistFile->cd(dir);
  Int_t l; 

  char name[200];
  sprintf(name, "mes%s%d", le, 1);  if (fGoodAccLepton) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
  sprintf(name, "mes%s%d", le, 2);  if (fGoodChargeCorr) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
  sprintf(name, "mes%s%d", le, 3);  if (fGoodLepton) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
  sprintf(name, "mes%s%d", le, 4);  if (fOneLepton) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
  sprintf(name, "mes%s%d", le, 5);  if (fGoodPRMM2) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
  sprintf(name, "mes%s%d", le, 6);  if (fGoodMM2) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
  sprintf(name, "mes%s%d", le, 7);  if (fGoodChargeCons) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
  sprintf(name, "mes%s%d", le, 8);  if (fGoodEvent) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
  
  sprintf(name, "mes%s%d", le, 11);      
  if (fGoodAccLepton) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
  sprintf(name, "mes%s%d", le, 12);      
  if (fGoodAccLepton&&fGoodChargeCorr) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
  sprintf(name, "mes%s%d", le, 13);      
  if (fGoodAccLepton&&fGoodChargeCorr&&fGoodLepton) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
  sprintf(name, "mes%s%d", le, 14);      
  if (fGoodAccLepton&&fGoodChargeCorr&&fGoodLepton&&fOneLepton) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
  sprintf(name, "mes%s%d", le, 15);      
  if (fGoodAccLepton&&fGoodChargeCorr&&fGoodLepton&&fOneLepton&&fGoodPRMM2) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
  sprintf(name, "mes%s%d", le, 16);      
  if (fGoodAccLepton&&fGoodChargeCorr&&fGoodLepton&&fOneLepton&&fGoodPRMM2&&fGoodMM2) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
  sprintf(name, "mes%s%d", le, 17);      
  if (fGoodAccLepton&&fGoodChargeCorr&&fGoodLepton&&fOneLepton&&fGoodPRMM2&&fGoodMM2&&fGoodChargeCons) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
  sprintf(name, "mes%s%d", le, 18);      
  if (fGoodEvent) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
  
  sprintf(name, "mes%s%d", le, 21);      
  if (fGoodChargeCorr&&fGoodPRMM2&&fGoodMM2&&fGoodChargeCons) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
  sprintf(name, "mes%s%d", le, 22);      
  if (fGoodAccLepton&&fGoodPRMM2&&fGoodMM2&&fGoodChargeCons) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
  sprintf(name, "mes%s%d", le, 23);      
  if (fGoodAccLepton&&fGoodChargeCorr&&fGoodPRMM2&&fGoodMM2&&fGoodChargeCons) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
  sprintf(name, "mes%s%d", le, 24);      
  if (fGoodLepton&&fGoodChargeCorr&&fGoodPRMM2&&fGoodMM2&&fGoodChargeCons) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
  sprintf(name, "mes%s%d", le, 25);      
  if (fOneLepton&&fGoodLepton&&fGoodChargeCorr&&fGoodMM2&&fGoodChargeCons) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
  sprintf(name, "mes%s%d", le, 26);      
  if (fOneLepton&&fGoodLepton&&fGoodChargeCorr&&fGoodPRMM2&&fGoodChargeCons) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
  sprintf(name, "mes%s%d", le, 27);      
  if (fOneLepton&&fGoodLepton&&fGoodChargeCorr&&fGoodPRMM2&&fGoodMM2) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
  sprintf(name, "mes%s%d", le, 28);      
  if (1) ((TH1D*)gDirectory->Get(name))->Fill(fMes);
  
  fHistFile->cd();
  return;
}

// ----------------------------------------------------------------------
void VubAnalysisCode::fillRecoilHist(const char *dir) {
  fHistFile->cd(dir);

  Int_t i, l, thebin; 
  double pnuRes = fPNu - p4MissGen.Vect().Mag();
  double tnuRes = TMath::Cos(fTNu) - TMath::Cos(p4MissGen.Theta());
  double q2Res  = fQ2 - fQ2Gen;
  
  double xmassRes    = fMxhad - fMxhadGen;
  double xmassfitRes = fMxhadfit - fMxhadGen;
  if (fOptGammas) {
    for(int jh =0; jh<15; jh++) {
      tmpxmassResF[jh] = tmpfMxhadfit[jh]- fMxhadGen; 
      //      cout << tmpxmassResF[6] << " " << tmpfMxhadfit[6] << endl;
      tmpxmassRes[jh] = tmpfMxhad[jh]- fMxhadGen; 
    }
  }

  // -- Kshorts
  for (Int_t iks = 0; iks < nKs; ++iks) {
    double mass = massKs[iks];
    if (goodKshort[iks] == 0) continue;   // skip kshorts whose daughters make up better combinations
    ((TH1D*)gDirectory->Get("ks100"))->Fill(mass);
    ((TH1D*)gDirectory->Get("ks200"))->Fill(pKs[iks]);
    ((TH1D*)gDirectory->Get("ks210"))->Fill(thKs[iks]*DR);
    ((TH1D*)gDirectory->Get("ks220"))->Fill(pKs[iks]);

    if (fGoodAccLepton&&fGoodChargeCorr)     {
      ((TH1D*)gDirectory->Get("ks101"))->Fill(mass);
      ((TH1D*)gDirectory->Get("ks201"))->Fill(pKs[iks]);
      ((TH1D*)gDirectory->Get("ks211"))->Fill(thKs[iks]*DR);
      ((TH1D*)gDirectory->Get("ks221"))->Fill(pKs[iks]);
    }
    if (fGoodLepton&&fOneLepton)        {
      ((TH1D*)gDirectory->Get("ks102"))->Fill(mass);
      ((TH1D*)gDirectory->Get("ks202"))->Fill(pKs[iks]);
      ((TH1D*)gDirectory->Get("ks212"))->Fill(thKs[iks]*DR);
      if (fGoodChargeCorr) ((TH1D*)gDirectory->Get("ks222"))->Fill(pKs[iks]);
    }
    if (fGoodChargeCons) {
      ((TH1D*)gDirectory->Get("ks103"))->Fill(mass);
      ((TH1D*)gDirectory->Get("ks203"))->Fill(pKs[iks]);
      ((TH1D*)gDirectory->Get("ks213"))->Fill(thKs[iks]*DR);
      if (fGoodLepton&&fOneLepton&&fGoodChargeCorr) ((TH1D*)gDirectory->Get("ks223"))->Fill(pKs[iks]);
    }
    if (fGoodMM2) {
      ((TH1D*)gDirectory->Get("ks104"))->Fill(mass);
      ((TH1D*)gDirectory->Get("ks204"))->Fill(pKs[iks]);
      ((TH1D*)gDirectory->Get("ks214"))->Fill(thKs[iks]*DR);
      if (fGoodLepton&&fOneLepton&&fGoodChargeCorr&&fGoodChargeCons) ((TH1D*)gDirectory->Get("ks224"))->Fill(pKs[iks]);
    }
    if (fGoodEvent)      {
      ((TH1D*)gDirectory->Get("ks105"))->Fill(mass);
      ((TH1D*)gDirectory->Get("ks205"))->Fill(pKs[iks]);
      ((TH1D*)gDirectory->Get("ks215"))->Fill(thKs[iks]*DR);
      if (fGoodEvent) ((TH1D*)gDirectory->Get("ks225"))->Fill(pKs[iks]);
    }
  }

  // -- charged kaons
  for (i = 0; i < nTrk; ++i) {
    if (goodTrack[i] == 0) continue;
    if (goodChargedKaon[i] == 0)  continue;
    ((TH1D*)gDirectory->Get("kp200"))->Fill(momentumTrk[i]);
    ((TH1D*)gDirectory->Get("kp210"))->Fill(thetaTrk[i]*DR);
    ((TH1D*)gDirectory->Get("kp220"))->Fill(momentumTrk[i]);
    if (fGoodAccLepton&&fGoodChargeCorr)     {
      ((TH1D*)gDirectory->Get("kp201"))->Fill(momentumTrk[i]);
      ((TH1D*)gDirectory->Get("kp211"))->Fill(thetaTrk[i]*DR);
      ((TH1D*)gDirectory->Get("kp221"))->Fill(momentumTrk[i]);
    }
    if (fGoodLepton&&fOneLepton)        {
      ((TH1D*)gDirectory->Get("kp202"))->Fill(momentumTrk[i]);
      ((TH1D*)gDirectory->Get("kp212"))->Fill(thetaTrk[i]*DR);
      if (fGoodChargeCorr) ((TH1D*)gDirectory->Get("kp222"))->Fill(momentumTrk[i]);
    }
    if (fGoodChargeCons) {
      ((TH1D*)gDirectory->Get("kp203"))->Fill(momentumTrk[i]);
      ((TH1D*)gDirectory->Get("kp213"))->Fill(thetaTrk[i]*DR);
      if (fGoodLepton&&fOneLepton&&fGoodChargeCorr) ((TH1D*)gDirectory->Get("kp223"))->Fill(momentumTrk[i]);
    }
    if (fGoodMM2) {
      ((TH1D*)gDirectory->Get("kp204"))->Fill(momentumTrk[i]);
      ((TH1D*)gDirectory->Get("kp214"))->Fill(thetaTrk[i]*DR);
      if (fGoodLepton&&fOneLepton&&fGoodChargeCorr&&fGoodChargeCons) ((TH1D*)gDirectory->Get("kp224"))->Fill(momentumTrk[i]);
    }
    if (fGoodEvent)      {
      ((TH1D*)gDirectory->Get("kp205"))->Fill(momentumTrk[i]);
      ((TH1D*)gDirectory->Get("kp215"))->Fill(thetaTrk[i]*DR);
      if (fGoodEvent) ((TH1D*)gDirectory->Get("kp225"))->Fill(momentumTrk[i]);
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
      if (fOptGammas && i == 0 && l == 2) {

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
	sprintf(name, "%s%d", le, 1500+i);   ((TH1D*)gDirectory->Get(name))->Fill(fTlab*DR, 1.);
	Int_t Brectrktmp;
	double dca, dcaz;
	for (int itrk = 0; itrk < nTrk; ++itrk) {
	  dcaz = zPocaTrk[itrk] - beamSZ;
	  dca  = TMath::Sqrt((xPocaTrk[itrk]-beamSX)*(xPocaTrk[itrk]-beamSX) + (yPocaTrk[itrk]-beamSY)*(yPocaTrk[itrk]-beamSY));
	  ((TH1D*)fHistFile->Get("dca"))->Fill(dca);
	  ((TH1D*)fHistFile->Get("dcaz"))->Fill(dcaz);
	  if (isRecEl(itrk)) {
	    ((TH1D*)fHistFile->Get("edca"))->Fill(dca);
	    ((TH1D*)fHistFile->Get("edcaz"))->Fill(dcaz);
	  }
	  if (isRecMu(itrk)) {
	    ((TH1D*)fHistFile->Get("mdca"))->Fill(dca);
	    ((TH1D*)fHistFile->Get("mdcaz"))->Fill(dcaz);
	  }

	  if (goodTrack[itrk] == 0) continue;
	  ((TH1D*)fHistFile->Get("adca"))->Fill(dca);
	  ((TH1D*)fHistFile->Get("adcaz"))->Fill(dcaz);
	  if (isRecLepton(itrk)) continue;
	  Brectrktmp = B0RecTrk[itrk];
	  if(fChB) Brectrktmp = chBRecTrk[itrk];
	  if ((Brectrktmp&brecoOverlap)) continue;
          // if (kshortLockTrk[itrk] == 1) continue; 
	  sprintf(name, "%s%d", le, 1600+i);   ((TH1D*)gDirectory->Get(name))->Fill(momentumTrk[itrk], weightTrk[itrk]);
	  sprintf(name, "%s%d", le, 1700+i);   ((TH1D*)gDirectory->Get(name))->Fill(thetaTrk[itrk]*DR, weightTrk[itrk]);
	  thebin = (momentumTrk[itrk] / (PCMAX/20)) + 1;
	  if(momentumTrk[itrk]>3.) thebin = 20;
          if(momentumTrk[itrk]<0.) thebin = 1;
	}


        for (int igam = 0; igam < nGam; ++igam) {
          if (goodPhoton[igam] == 0) continue;
          Brectrktmp = B0RecGam[igam];
          if(fChB) Brectrktmp = chBRecGam[igam];
          if ((Brectrktmp&brecoOverlap)) continue;
          sprintf(name, "%s%d", le, 1800+i);   ((TH1D*)gDirectory->Get(name))->Fill(energyGam[igam], weightNeu[igam]);
          sprintf(name, "%s%d", le, 1900+i);   ((TH1D*)gDirectory->Get(name))->Fill(thetaGam[igam]*DR, weightNeu[igam]);
	  thebin = (energyGam[igam] / (1.5/20)) + 1; 
	  if(energyGam[igam]>1.5) thebin = 20;
	  if(energyGam[igam]<0.) thebin = 1;
	  if(l == 2 && i == 0){
	    sprintf(name,"%s%d%s%d","a",9900,"bin",thebin);    ((TH1D*)gDirectory->Get(name))->Fill(fMes, weightNeu[igam]);
	    if(fIntPurity > 60)	 {sprintf(name,"%s%d%s%d","pur60a",9900,"bin",thebin);    ((TH1D*)gDirectory->Get(name))->Fill(fMes);}
	    if(fIntPurity > 80)	 {sprintf(name,"%s%d%s%d","pur80a",9900,"bin",thebin);    ((TH1D*)gDirectory->Get(name))->Fill(fMes);}
	  }

        }

	sprintf(name, "%s%d", le, 4400+i);  ((TH1D*)gDirectory->Get(name))->Fill(fNKp, 1.);
	sprintf(name, "%s%d", le, 4500+i);  ((TH1D*)gDirectory->Get(name))->Fill(fNKshort, 1.);
	sprintf(name, "%s%d", le, 4600+i);  ((TH1D*)gDirectory->Get(name))->Fill(fNLepton, 1.);

	// -- all other cuts
	if (faoLepton) {
	  sprintf(name, "%s%d", le, 1020+i);  ((TH1D*)gDirectory->Get(name))->Fill(fPcms, 1.);
	  sprintf(name, "%s%d", le, 1120+i);   if (fIsPrompt) ((TH1D*)gDirectory->Get(name))->Fill(fPcms, 1.);
	  sprintf(name, "%s%d", le, 1220+i);   if (fIsCascade) ((TH1D*)gDirectory->Get(name))->Fill(fPcms, 1.);
	  sprintf(name, "%s%d", le, 1320+i);   ((TH1D*)gDirectory->Get(name))->Fill(TMath::Cos(fTcms), 1.);
	  sprintf(name, "%s%d", le, 1420+i);   ((TH1D*)gDirectory->Get(name))->Fill(fPups, 1.);
	  sprintf(name, "%s%d", le, 1520+i);   ((TH1D*)gDirectory->Get(name))->Fill(fTlab*DR, 1.);
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
	  sprintf(name, "%s%d", le, 2420+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhadfit);
	  sprintf(name, "%s%d", le, 2520+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhadfit);
	  sprintf(name, "%s%d", le, 2620+i);   ((TH1D*)gDirectory->Get(name))->Fill(TMath::Cos(fTNu));
	  sprintf(name, "%s%d", le, 2720+i);   ((TH1D*)gDirectory->Get(name))->Fill(fPNu);
	  sprintf(name, "%s%d", le, 2820+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhadfit);
	  sprintf(name, "%s%d", le, 2920+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhadfit);
	  sprintf(name, "%s%d", le, 3920+i);   ((TH1D*)gDirectory->Get(name))->Fill(fGammaMax);
	  sprintf(name, "%s%d", le, 9020+i);   ((TH1D*)gDirectory->Get(name))->Fill(xmassRes);
	  sprintf(name, "%s%d", le, 9420+i);   ((TH1D*)gDirectory->Get(name))->Fill(xmassfitRes);
	}
      }

      if (fGoodLepton) {
	sprintf(name, "%s%d", le, 2000+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhad);
	sprintf(name, "%s%d", le, 2100+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhad);
	sprintf(name, "%s%d", le, 2200+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhadGen);
	sprintf(name, "%s%d", le, 2400+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhadfit);
	sprintf(name, "%s%d", le, 2500+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhadfit);
	sprintf(name, "%s%d", le, 2600+i);   ((TH1D*)gDirectory->Get(name))->Fill(TMath::Cos(fTNu));
	sprintf(name, "%s%d", le, 2700+i);   ((TH1D*)gDirectory->Get(name))->Fill(fPNu);
	sprintf(name, "%s%d", le, 2800+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhadfit);
	sprintf(name, "%s%d", le, 2900+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhadfit);
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
        	
      }
      // -- after cuts
      if (fGoodEvent) {
	((TH1D*)gDirectory->Get("mNuSqNC"))->Fill(fMM2NC);
        sprintf(name, "%s%d", le, 1010+i);   ((TH1D*)gDirectory->Get(name))->Fill(fPcms);
	sprintf(name, "%s%d", le, 1110+i);   if (fIsPrompt) ((TH1D*)gDirectory->Get(name))->Fill(fPcms);
	sprintf(name, "%s%d", le, 1210+i);   if (fIsCascade) ((TH1D*)gDirectory->Get(name))->Fill(fPcms);
	sprintf(name, "%s%d", le, 1310+i);   ((TH1D*)gDirectory->Get(name))->Fill(TMath::Cos(fTcms));
	sprintf(name, "%s%d", le, 1410+i);   ((TH1D*)gDirectory->Get(name))->Fill(fPups);
	sprintf(name, "%s%d", le, 1510+i);   ((TH1D*)gDirectory->Get(name))->Fill(fTlab*DR);

	Int_t Brectrktmp;
	for (int itrk = 0; itrk < nTrk; ++itrk) {
	  if (goodTrack[itrk] == 0) continue;
	  if (isRecLepton(itrk)) continue;
	  Brectrktmp = B0RecTrk[itrk];
	  if(fChB) Brectrktmp = chBRecTrk[itrk];
	  if ((Brectrktmp&brecoOverlap)) continue;
	  // if (kshortLockTrk[itrk] == 1) continue; 
	  sprintf(name, "%s%d", le, 1610+i);   ((TH1D*)gDirectory->Get(name))->Fill(momentumTrk[itrk], weightTrk[itrk]);
	  sprintf(name, "%s%d", le, 1710+i);   ((TH1D*)gDirectory->Get(name))->Fill(thetaTrk[itrk]*DR, weightTrk[itrk]);
	  sprintf(name, "%s%d", le, 1620+i);   ((TH1D*)gDirectory->Get(name))->Fill(momentumTrk[itrk], weightTrk[itrk]);
	  sprintf(name, "%s%d", le, 1720+i);   ((TH1D*)gDirectory->Get(name))->Fill(thetaTrk[itrk]*DR, weightTrk[itrk]);
	  thebin = (momentumTrk[itrk] / (PCMAX/20)) + 1;
	  if(momentumTrk[itrk]>3.) thebin = 20;
          if(momentumTrk[itrk]<0.) thebin = 1;
	  if(l == 2 && i == 0){
	    sprintf(name,"%s%d%s%d","a",1610,"bin",thebin);    ((TH1D*)gDirectory->Get(name))->Fill(fMes, weightTrk[itrk]);
	    if(fIntPurity > 60) {sprintf(name,"%s%d%s%d","pur60a",1610,"bin",thebin);    ((TH1D*)gDirectory->Get(name))->Fill(fMes);}
	    if(fIntPurity > 80) {sprintf(name,"%s%d%s%d","pur80a",1610,"bin",thebin);    ((TH1D*)gDirectory->Get(name))->Fill(fMes);}
	  }
	}

	for (int igam = 0; igam < nGam; ++igam) {
	  if (goodPhoton[igam] == 0) continue;
	  Brectrktmp = B0RecGam[igam];
	  if(fChB) Brectrktmp = chBRecGam[igam];
	  if ((Brectrktmp&brecoOverlap)) continue;
	  //      if (kshortLockTrk[itrk] == 1) continue; 
	  sprintf(name, "%s%d", le, 1810+i);   ((TH1D*)gDirectory->Get(name))->Fill(energyGam[igam], weightNeu[igam]);
	  sprintf(name, "%s%d", le, 1910+i);   ((TH1D*)gDirectory->Get(name))->Fill(thetaGam[igam]*DR, weightNeu[igam]);
	  sprintf(name, "%s%d", le, 1820+i);   ((TH1D*)gDirectory->Get(name))->Fill(energyGam[igam], weightNeu[igam]);
	  sprintf(name, "%s%d", le, 1920+i);   ((TH1D*)gDirectory->Get(name))->Fill(thetaGam[igam]*DR, weightNeu[igam]);
	  thebin = (energyGam[igam] / (1.5/20)) + 1; 
	  if(energyGam[igam]>1.5) thebin = 20;
	  if(energyGam[igam]<0.) thebin = 1;
	  if(l == 2 && i == 0){
	    sprintf(name,"%s%d%s%d","a",9910,"bin",thebin);    ((TH1D*)gDirectory->Get(name))->Fill(fMes, weightNeu[igam]);
	    if(fIntPurity > 60)	{sprintf(name,"%s%d%s%d","pur60a",9910,"bin",thebin);    ((TH1D*)gDirectory->Get(name))->Fill(fMes);}
	    if(fIntPurity > 80)	{sprintf(name,"%s%d%s%d","pur80a",9910,"bin",thebin);    ((TH1D*)gDirectory->Get(name))->Fill(fMes);}
	  }
	}


        sprintf(name, "%s%d", le, 2010+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhad);
	sprintf(name, "%s%d", le, 2110+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhad);
        sprintf(name, "%s%d", le, 2210+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhadGen);
	sprintf(name, "%s%d", le, 2410+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhadfit);
	sprintf(name, "%s%d", le, 2510+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhadfit);
	sprintf(name, "%s%d", le, 2610+i);   ((TH1D*)gDirectory->Get(name))->Fill(TMath::Cos(fTNu));
	sprintf(name, "%s%d", le, 2710+i);   ((TH1D*)gDirectory->Get(name))->Fill(fPNu);
	sprintf(name, "%s%d", le, 2810+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhadfit);
	sprintf(name, "%s%d", le, 2910+i);   ((TH1D*)gDirectory->Get(name))->Fill(fMxhadfit);

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

	sprintf(name, "%s%d", le, 4410+i);  ((TH1D*)gDirectory->Get(name))->Fill(fNKp, 1.);
	sprintf(name, "%s%d", le, 4510+i);  ((TH1D*)gDirectory->Get(name))->Fill(fNKshort, 1.);
	sprintf(name, "%s%d", le, 4610+i);  ((TH1D*)gDirectory->Get(name))->Fill(fNLepton, 1.);

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
void VubAnalysisCode::mcTruth() {

  if (!fIsMC) return;
  
  Bool_t isSemilep(kFALSE);
  Int_t ipB[] = {-1, -1};
  Int_t typB[] = {7, 7};
  Int_t cnt(-1);

  static int counter(0); 

  Int_t Imc(-1), Jmc(-1), Kmc(-1); 
  fB1Index = fB2Index = fBVxb = fBVcb = fBVub =  -99;
  fBbchgen = 0;
  // -- Find indices of B's
  for (Int_t imc = 0; imc < nMc; ++imc) {
    if ((TMath::Abs(idMc[imc]) == 511) || (TMath::Abs(idMc[imc]) == 521)) {
      if(TMath::Abs(idMc[imc]) == 511) {
	fBbchgen = 2; //B0s
      } else if (TMath::Abs(idMc[imc]) == 521) {
	fBbchgen = 1; //B+/-s
      }
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
  Int_t ib, dCount(0);
  float dLen(-1.);
  fDpi = fDpiz = fDk = fDkmiss = fDks = fDlep = fDgam = 0;
  fDnu = fDkspipi = fDkspiopio = fDkl = 0;
  Bool_t isCharm(kFALSE);
  
  for (ib = 0; ib < 2; ++ib) { 
    for (Int_t imc = 0; imc < nMc; ++imc) {
      if ((mothMc[imc]-1) == ipB[ib]) {
 

	//D systematics
	if (isTruLepton(imc)) {
	  //There's a lepton
	  for (Int_t s = 0; s < nMc; ++s) {
	    if ((TMath::Abs(idMc[s]) == 421) && isAncestor(ipB[ib],s) && (dCount == 0)) {
	      //D0
	      if (((TMath::Abs(idMc[mothMc[s]-1]) >= 400) && (TMath::Abs(idMc[mothMc[s]-1]) <= 499)) || ((TMath::Abs(idMc[mothMc[s]-1]) >= 10400) && (TMath::Abs(idMc[mothMc[s]-1]) <= 10499))) {
		//Flag it as a D0 coming from dstar
		fD0CfDs = 1;
	      } else if ((TMath::Abs(idMc[mothMc[s]-1]) == 511) || (TMath::Abs(idMc[mothMc[s]-1]) == 521)) {
		fD0CfDs = 2;
	      } else {
		cout<<"What the hell: "<< TMath::Abs(idMc[mothMc[s]-1]) <<endl;
	      }
	      dCount++;
	      for (Int_t imcD = 0; imcD < nMc; ++imcD) { //Loop on direct D dau
		if ((mothMc[imcD]-1) == s) {
		  if (idMc[imcD] == 111) fDpiz++;   
		  dLen = sqrt(pow((xMc[imcD]-xMc[s]),2)+pow((yMc[imcD]-yMc[s]),2)+pow((zMc[imcD]-zMc[s]),2));
		  if( dLen > 0.05) {
		    if (TMath::Abs(idMc[imcD]) == 211) fDpi++;   
		    if (TMath::Abs(idMc[imcD]) == 22)  fDgam++;   
		    if ((TMath::Abs(idMc[imcD]) == 12)||(TMath::Abs(idMc[imcD]) == 14)||(TMath::Abs(idMc[imcD]) == 16))  fDnu++;   
		    if (TMath::Abs(idMc[imcD]) == 321) {
		      fDk++;   
		      //How to identify the missed K+?
		      if(goodTrack[imcD] == 0) fDkmiss++;
		    }
		    if ((idMc[imcD] == 310) || (idMc[imcD] == 130)) fDks++;   
		    if (idMc[imcD] == 130) fDkl++;   
		    if (idMc[imcD] == 310) {
		      for (Int_t imcKd = 0; imcKd < nMc; ++imcKd) {
			if((mothMc[imcKd]-1) == imcD) {
			  if(idMc[imcKd] == 211) {
			    fDkspipi++;   
			  } else if (idMc[imcKd] == 111) {
			    fDkspiopio++;   
			  }
			}
		      }
		    }
		    if (isTruLepton(imcD)) fDlep++;   
		  } else if (dLen <= 0.05)   {
		    for (Int_t imcDD = 0; imcDD < nMc; ++imcDD) {
		      dLen = sqrt(pow((xMc[imcDD]-xMc[imcD]),2)+pow((yMc[imcDD]-yMc[imcD]),2)+pow((zMc[imcDD]-zMc[imcD]),2));
		      if ((mothMc[imcDD]-1) == imcD) { 
			if (idMc[imcDD] == 111) fDpiz++;   
			if(dLen > 0.05) {
			  if (TMath::Abs(idMc[imcDD]) == 211) fDpi++;   
			  if (TMath::Abs(idMc[imcDD]) == 321) {
			    fDk++;   
			    //How to identify the missed K+?
			    if(goodTrack[imcDD] == 0) fDkmiss++;
			  }
			  if ((idMc[imcDD] == 310) || (idMc[imcDD] == 130))fDks++;   
			  if (idMc[imcDD] == 130) fDkl++;   
			  if (idMc[imcDD] == 310) {
			    for (Int_t imcKd = 0; imcKd < nMc; ++imcKd) {
			      if((mothMc[imcKd]-1) == imcDD) {
				if(idMc[imcKd] == 211) {
				  fDkspipi++;   
				} else if (idMc[imcKd] == 111) {
				  fDkspiopio++;   
				}
			      }
			    }
			  }
			  if (isTruLepton(imcDD)) fDlep++;   
			} else if (dLen <= 0.05) {
			  for (Int_t imcDDD = 0; imcDDD < nMc; ++imcDDD) {
			    dLen = sqrt(pow((xMc[imcDDD]-xMc[imcDD]),2)+pow((yMc[imcDDD]-yMc[imcDD]),2)+pow((zMc[imcDDD]-zMc[imcDD]),2));
			    if ((mothMc[imcDDD]-1) == imcDD) {
			      if (idMc[imcDDD] == 111) fDpiz++;   
			      if (dLen > 0.05) {
				if (TMath::Abs(idMc[imcDDD]) == 211) fDpi++;   
				if (TMath::Abs(idMc[imcDDD]) == 321) {
				  fDk++;   
				  //How to identify the missed K+?
				  if(goodTrack[imcDDD] == 0) fDkmiss++;
				}
				if ((idMc[imcDDD] == 310) || (idMc[imcDDD] == 130)) fDks++;   
				if (idMc[imcDDD] == 130) fDkl++;   
				if (idMc[imcDDD] == 310) {
				  for (Int_t imcKd = 0; imcKd < nMc; ++imcKd) {
				    if((mothMc[imcKd]-1) == imcDDD) {
				      if(idMc[imcKd] == 211) {
					fDkspipi++;   
				      } else if (idMc[imcKd] == 111) {
					fDkspiopio++;   
				      }
				    }
				  }
				}
				if (isTruLepton(imcDDD)) fDlep++;   
			      }
			    }
			  }
			}
		      }
		    }
		  } 
		}
	      }
	       
	    } else if ((TMath::Abs(idMc[s]) == 411) && isAncestor(ipB[ib],s) && (dCount == 0)) {
	      //D+
	      if (((TMath::Abs(idMc[mothMc[s]-1]) >= 400) && (TMath::Abs(idMc[mothMc[s]-1]) <= 499)) || ((TMath::Abs(idMc[mothMc[s]-1]) >= 10400) && (TMath::Abs(idMc[mothMc[s]-1]) <= 10499))) {
		//Flag it as a D0 coming from dstar
		fDCfDs = 1;
	      } else if ((TMath::Abs(idMc[mothMc[s]-1]) == 511) || (TMath::Abs(idMc[mothMc[s]-1]) == 521)) {
		fDCfDs = 2;
	      } else {
		cout<<"What the hell: "<< TMath::Abs(idMc[mothMc[s]-1]) <<endl;
	      }
	      for (Int_t imcD = 0; imcD < nMc; ++imcD) {
		if ((mothMc[imcD]-1) == s) {
		  if (TMath::Abs(idMc[imcD]) == 111) fDpiz++;   
		  dLen = sqrt(pow((xMc[imcD]-xMc[s]),2)+pow((yMc[imcD]-yMc[s]),2)+pow((zMc[imcD]-zMc[s]),2));
		  if( dLen > 0.05) {
		    if (TMath::Abs(idMc[imcD]) == 22)  fDgam++;   
		    if (TMath::Abs(idMc[imcD]) == 211) fDpi++;   
		    if ((TMath::Abs(idMc[imcD]) == 12)||(TMath::Abs(idMc[imcD]) == 14)||(TMath::Abs(idMc[imcD]) == 16))  fDnu++;   
		    if (TMath::Abs(idMc[imcD]) == 321) {
		      fDk++;   
		      //How to identify the missed K+?
		      if(goodTrack[imcD] == 0) fDkmiss++;
		    }
		    if ((TMath::Abs(idMc[imcD]) == 310) || (TMath::Abs(idMc[imcD]) == 130)) fDks++;   
		    if (idMc[imcD] == 130) fDkl++;   
		    if (idMc[imcD] == 310) {
		      for (Int_t imcKd = 0; imcKd < nMc; ++imcKd) {
			if((mothMc[imcKd]-1) == imcD) {
			  if(idMc[imcKd] == 211) {
			    fDkspipi++;   
			  } else if (idMc[imcKd] == 111) {
			    fDkspiopio++;   
			  }
			}
		      }
		    }
		    if (isTruLepton(imcD)) fDlep++;   
		  } else if (dLen <= 0.05)   {
		    for (Int_t imcDD = 0; imcDD < nMc; ++imcDD) {
		      dLen = sqrt(pow((xMc[imcDD]-xMc[imcD]),2)+pow((yMc[imcDD]-yMc[imcD]),2)+pow((zMc[imcDD]-zMc[imcD]),2));
		      if ((mothMc[imcDD]-1) == imcD) { 
			if (TMath::Abs(idMc[imcDD]) == 111) fDpiz++;   
			if(dLen > 0.05) {
			  if (TMath::Abs(idMc[imcDD]) == 211) fDpi++;   
			  if (TMath::Abs(idMc[imcDD]) == 321) {
			    fDk++;   
			    //How to identify the missed K+?
			    if(goodTrack[imcDD] == 0) fDkmiss++;
			  }
			  if ((TMath::Abs(idMc[imcDD]) == 310) || (TMath::Abs(idMc[imcDD]) == 130)) fDks++;   
			  if (idMc[imcDD] == 130) fDkl++;   
			  if (idMc[imcDD] == 310) {
			    for (Int_t imcKd = 0; imcKd < nMc; ++imcKd) {
			      if((mothMc[imcKd]-1) == imcDD) {
				if(idMc[imcKd] == 211) {
				  fDkspipi++;   
				} else if (idMc[imcKd] == 111) {
				  fDkspiopio++;   
				} 
			      }
			    }
			  }
			  if (isTruLepton(imcDD)) fDlep++;   
			} else if (dLen <= 0.05) {
			  for (Int_t imcDDD = 0; imcDDD < nMc; ++imcDDD) {
			    dLen = sqrt(pow((xMc[imcDDD]-xMc[imcDD]),2)+pow((yMc[imcDDD]-yMc[imcDD]),2)+pow((zMc[imcDDD]-zMc[imcDD]),2));
			    if ((mothMc[imcDDD]-1) == imcDD) { 
			      if (TMath::Abs(idMc[imcDDD]) == 111) fDpiz++;   
			      if(dLen > 0.05) {
				if (TMath::Abs(idMc[imcDDD]) == 211) fDpi++;   
				if (TMath::Abs(idMc[imcDDD]) == 321) {
				  fDk++;   
				  //How to identify the missed K+?
				  if(goodTrack[imcDDD] == 0) fDkmiss++;
				}
				if ((TMath::Abs(idMc[imcDDD]) == 310) || (TMath::Abs(idMc[imcDDD]) == 130)) fDks++;   
				if (idMc[imcDDD] == 130) fDkl++;   
				if (idMc[imcDDD] == 310) {
				  for (Int_t imcKd = 0; imcKd < nMc; ++imcKd) {
				    if((mothMc[imcKd]-1) == imcDDD) {
				      if(idMc[imcKd] == 211) {
					fDkspipi++;   
				      } else if (idMc[imcKd] == 111) {
					fDkspiopio++;   
				      }
				    }
				  }
				}
				if (isTruLepton(imcDDD)) fDlep++;   
			      }
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
	
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
	      fBVcb = ipB[ib];
	      typB[ib] = 2;
	      isCharm = kTRUE;
	      break;
	    } else {
	      fBVub = ipB[ib];
	    }

	  }
	}

      }
    }
  }


  Int_t idmc(0);
  fWithKKbar = kFALSE; 
  for (Int_t s = 0; s < nMc; ++s) {
    idmc = TMath::Abs(idMc[s]);
    if (isAncestor(fBVub, s)  && ((idmc == 321) || (idmc == 311) || (idmc == 310) || (idmc == 130))) {
      //        cout << " ======================================================================" << endl;
      //        cout << "Found KKbar" << endl;
      //        cout << counter++ << "/" << fReturnLog[0] 
      //  	   << "  " << idmc << "  " << TMath::Abs(idMc[mothMc[s]-1]) << endl;
      //        cout << " ----------------------------------------------------------------------" << endl;
      //        dumpGeneratorBlock(); 
      //        cout << " ----------------------------------------------------------------------" << endl;
      fWithKKbar = kTRUE; 
      break;
    }
  }




  fHistFile->cd("mcTruth");
  fBVSTyp = -1;
  // -- Determine special modes of sl B decay : fBVxbTyp (Gvxbtyp in the tree) : 
  //    1 D l nu (gamma)
  //    2 D* l nu (gamma)
  //    4 D2*l nu (gamma)
  //    5 D1*lnu (gamma)
  // the sign of fBVxbTyp tells if B0 (+) or B+/- (-)	

  fBVxbTyp = -99; 
  int nDaureq, nDauSt;
  nDaureq = 3;
  nDauSt = 0;
  for (Int_t ss = 0; ss < nMc; ++ss) {
    int idmc = TMath::Abs(idMc[ss]);
    if (mothMc[ss]-1 == fBVcb) {
      //      cout<<"I'd enetered at all: "<<fBVcb<<endl;
      //      cout<<idmc<<endl;
      if ((idmc == 411) || (idmc == 421)) {
	//D decay
	for (Int_t ssl = 0; ssl < nMc; ++ssl) {
	  if((ssl != ss) && (mothMc[ss] == mothMc[ssl])) {
	    if( idMc[ssl] == 22) {
	      nDaureq++;
	    } else {
	      //Count whatever it is
	      nDauSt++;
	    }
	  }
	}
	if(nDauMc[fBVcb] == nDaureq)	fBVxbTyp = 1; 
	if(nDauSt >2) fBVxbTyp = 6; 
	break;
      }
      else if ((idmc == 413) || (idmc == 423)) {
	//D* decay
	//	dumpGeneratorBlock(fB1Index,fB2Index);
	//	cout<<"Event; daughters -> "<<nDauMc[fBVcb]<<endl;
	for (Int_t ssl = 0; ssl < nMc; ++ssl) {
	  if((ssl != ss) && (mothMc[ss] == mothMc[ssl])) {
	    if(idMc[ssl] == 22) {
	      nDaureq++;
	    } else {
	      nDauSt++;
	    }
	  }
	}
	if(nDauMc[fBVcb] == nDaureq) fBVxbTyp = 2; 
	if(nDauSt >2)	fBVxbTyp = 6;
	break;
      }
      else if ((idmc == 415) || (idmc == 425)) {
	//D_2* decay
	for (Int_t ssl = 0; ssl < nMc; ++ssl) {
	  if((ssl != ss) && (mothMc[ss] == mothMc[ssl])) {
	    if(idMc[ssl] == 22) {
	      nDaureq++;
	    } else {
	      nDauSt++;
	    }
	  }
	}
	if(nDauMc[fBVcb] == nDaureq) fBVxbTyp = 4; 
	if(nDauSt >2)	fBVxbTyp = 6; 
	break;
      }
      else if ((idmc == 10413) || (idmc == 10423)  ) {
	//D_1
	for (Int_t ssl = 0; ssl < nMc; ++ssl) {
	  if((ssl != ss) && (mothMc[ss] == mothMc[ssl])) {
	    if(idMc[ssl] == 22) {
	      nDaureq++;
	    } else {
	      nDauSt++;
	    }
	  }
	}
	if(nDauMc[fBVcb] == nDaureq) fBVxbTyp = 5; 
	if(nDauSt >2)	fBVxbTyp = 6; 
	break;
      }
      else if (((idmc >= 400) && (idmc <= 499)) || ((idmc >= 10400) && (idmc <= 10499)) || ((idmc >= 20400) && (idmc <= 20499)))  {
	//Other D
	for (Int_t ssl = 0; ssl < nMc; ++ssl) {
	  if((ssl != ss) && (mothMc[ss] == mothMc[ssl])) {
	    if(idMc[ssl] == 22) {
	      nDaureq++;
	    } else {
	      nDauSt++;
	    }
	  }
	}
	if(nDauMc[fBVcb] == nDaureq) fBVxbTyp = 6; 
	if(nDauSt >2)	fBVxbTyp = 6; 
	break;
      }
      /*      else if ((idmc >= 200) && (idmc <= 300) || (idmc >= 10200) && (idmc <= 10300) || (idmc >= 20200) && (idmc <= 20300) ) {
	//a_1, b_1, pi
	for (Int_t ssl = 0; ssl < nMc; ++ssl) {
	  if((ssl != ss) && (mothMc[ss] == mothMc[ssl]) && idMc[ssl] == 22) nDaureq++;
	}
	if(nDauMc[fBVcb] == nDaureq) fBVxbTyp = 8; 
	break;
	}*/
//*******************
 
    } else if (mothMc[ss]-1 == fBVub) {
      //OTHER SPECIAL MODES OF SL B DECAYS (Vub modes)

      //  b -> u l nu DECAYS splitted in exclusive channels:
      // the flag  fBVxbTyp (Gvxbtyp in the tree) describes each channel with
      // the convention:
      /*

      fBVxbTyp = 11       B -> pi l+ nu   
      fBVxbTyp = 12       B -> eta l+ nu
      fBVxbTyp = 13       B-> rho0 l+ nu
      fBVxbTyp = 14       B-> omega l+ nu 
      fBVxbTyp = 15       B -> eta_prime l+ nu
      fBVxbTyp = 16       B -> a_10+-,a_20+-,h_1...l nu
      fBVxbTyp = 17       B -> b_10+- l+ nu
      fBVxbTyp = 18       B -> f_0, f'_0, f_2...l nu
      fBVxbTyp = 19       B -> a_00+-...l nu
      fBVxbTyp = 20       B ->f'_2, h'_1...l nu
      fBVxbTyp = 21       B -> f_1 l+ nu 
      fBVxbTyp = 22       B -> f'_1 l+ n


      // the sign of fBVxbTyp tells if B0 (+) or B+/- (-)

      */

      // B -> pi0+/- l+ nu decay for charged or neutral B

      if ( (TMath::Abs(idmc) == 211 )|| (idmc == 111))  {
//	if (TMath::Abs(idMc[fBVub]) == 511)
//   	cout << "B0"<<endl;
	for (Int_t ssl = 0; ssl < nMc; ++ssl) {
	  if((ssl != ss) && (mothMc[ss] == mothMc[ssl]) && idMc[ssl] == 22) nDaureq++;
	}
	if(nDauMc[fBVub] == nDaureq) fBVxbTyp = 11; 
	break;
      }
	 
      // 2) B -> eta l+ nu decay for charged B 

	  else if (idmc == 221)  {

	    for (Int_t ssl = 0; ssl < nMc; ++ssl) {
	  if((ssl != ss) && (mothMc[ss] == mothMc[ssl]) && idMc[ssl] == 22) nDaureq++;
	}
	 if(nDauMc[fBVub] == nDaureq) fBVxbTyp = 12; 
	break;
      }
      // B-> rho0+/- l+ nu decay for charged or neutral B

	
      else if ((TMath::Abs(idmc) == 213)  || (idmc == 113)) {

	for (Int_t ssl = 0; ssl < nMc; ++ssl) {
	  if((ssl != ss) && (mothMc[ss] == mothMc[ssl]) && idMc[ssl] == 22) nDaureq++;
	}
	if(nDauMc[fBVub] == nDaureq) fBVxbTyp = 13; 
	break;
      }

      // B-> omega l+ nu decay for charged B

 
      else if (idmc == 223) {

	for (Int_t ssl = 0; ssl < nMc; ++ssl) {
	  if((ssl != ss) && (mothMc[ss] == mothMc[ssl]) && idMc[ssl] == 22) nDaureq++;
	}
	if(nDauMc[fBVub] == nDaureq) fBVxbTyp = 14;
	break;
      }
  
      //B -> eta_prime l+ nu decay for charged B

      else if (idmc == 331) {

	for (Int_t ssl = 0; ssl < nMc; ++ssl) {
	  if((ssl != ss) && (mothMc[ss] == mothMc[ssl]) && idMc[ssl] == 22) nDaureq++;
	}
	if(nDauMc[fBVub] == nDaureq) fBVxbTyp = 15; 
	break;
      }
      
      //Decays with final multiplicity rho pi (B -> a_10+-,a_20+-,h_1 etc)
      //for neutral or charged B
      


      else if ((TMath::Abs(idmc) == 20213) || (TMath::Abs(idmc) == 215) || 
	       (idmc == 10223) || (idmc == 20113) || (idmc == 115)) {
	
	for (Int_t ssl = 0; ssl < nMc; ++ssl) {
	  if((ssl != ss) && (mothMc[ss] == mothMc[ssl]) && idMc[ssl] == 22) nDaureq++;
	}
	if(nDauMc[fBVub] == nDaureq) fBVxbTyp = 16; 
	break;
      }
 
      // B -> b_10+- l+ nu decay for charged or neutral B

      else if ((TMath::Abs(idmc) == 10213) || (idmc == 10113)) {

	for (Int_t ssl = 0; ssl < nMc; ++ssl) {
	  if((ssl != ss) && (mothMc[ss] == mothMc[ssl]) && idMc[ssl] == 22) nDaureq++;
	}
	if(nDauMc[fBVub] == nDaureq) fBVxbTyp = 17; 
	break;
      }

      //Decays with final multiplicity pi pi (B -> a_10+-,a_20+-,h_1 etc)
      //for charged B
      
 
      // B -> f_0, f'_0, f_2 etc

	 else if ((idmc == 10221) || (idmc == 10331) || (idmc == 225)){

	for (Int_t ssl = 0; ssl < nMc; ++ssl) {
	  if((ssl != ss) && (mothMc[ss] == mothMc[ssl]) && idMc[ssl] == 22) nDaureq++;
	}
	if(nDauMc[fBVub] == nDaureq) fBVxbTyp = 18; 
	break;
      }

      //Decays with final multiplicity eta pi (B -> a_10+-,a_20+-,h_1 etc)
      //for neutral or charged B
      
      //B -> a_00+- etc

	 else if ((TMath::Abs(idmc) == 10211) || (idmc == 10111)){ 

	for (Int_t ssl = 0; ssl < nMc; ++ssl) {
	  if((ssl != ss) && (mothMc[ss] == mothMc[ssl]) && idMc[ssl] == 22) nDaureq++;
	}
	if(nDauMc[fBVub] == nDaureq) fBVxbTyp = 19; 
	break;
      }

      //Decays with final multiplicity k k (B -> a_10+-,a_20+-,h_1 etc)
      //for charged B
      
      //B ->f'_2, h'_1 etc
      
      else if ((TMath::Abs(idmc) == 335) || (TMath::Abs(idmc) == 10333)){
	
	for (Int_t ssl = 0; ssl < nMc; ++ssl) {
	  if((ssl != ss) && (mothMc[ss] == mothMc[ssl]) && idMc[ssl] == 22) nDaureq++;
	}
	if(nDauMc[fBVub] == nDaureq) fBVxbTyp = 20; 
	break;
      }


      // B -> f_1 l+ nu decay for charged B
      
      
      else if (TMath::Abs(idmc) == 20223){
	
	for (Int_t ssl = 0; ssl < nMc; ++ssl) {
	  if((ssl != ss) && (mothMc[ss] == mothMc[ssl]) && idMc[ssl] == 22) nDaureq++;
	}
	if(nDauMc[fBVub] == nDaureq) fBVxbTyp = 21; 
	break;
      }
      

      // B -> f'_1 l+ n decay for charged B
      
      else if (TMath::Abs(idmc) == 20333){
	
	for (Int_t ssl = 0; ssl < nMc; ++ssl) {
	  if((ssl != ss) && (mothMc[ss] == mothMc[ssl]) && idMc[ssl] == 22) nDaureq++;
	}
	if(nDauMc[fBVub] == nDaureq) fBVxbTyp = 22; 
	break;
      }
      
      
      
      //inclusive
      
      else if ((idmc == 41)||(idmc == 42)) {
	//Xu neutral
	for (Int_t ssl = 0; ssl < nMc; ++ssl) {
	  if((ssl != ss) && (mothMc[ss] == mothMc[ssl]) && idMc[ssl] == 22) nDaureq++;
	}
	if(nDauMc[fBVub] == nDaureq) fBVxbTyp = 7; 
	break;
      }	
      else if (isCharm) {
	if(fVerbose && (idmc != 10413) && (idmc != 20413) && (idmc != 415) && (idmc != 10411))	cout<<idmc<<" This is the mother ID :: charm"<<endl;
	fBVxbTyp = 3; 
      } else {
	if (fVerbose) cout<<idmc<<" This is the mother ID :: unass"<<endl;
      }
      //      else {
      //	if((idmc<10 || idmc>20) && (idmc != 42 && idmc != 211 )) cout<<idmc<< " " << nDauMc[fBVxb]<<endl;
      //}
    }
  }

// give a sign to fBVxbTyp depending on the id of the B meson (+ if B0(~), - if B+(-))


    if (fBVxb > -1 && fBVxbTyp>0) {
	if(TMath::Abs(idMc[fBVxb]) == 511){
	} else if(TMath::Abs(idMc[fBVxb]) == 521) {
	   fBVxbTyp = -fBVxbTyp; 
	} else {
	   cout <<"?!? wrong B moth? " << idMc[fBVxb] << endl;
	}	
    }
  
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
  double massBmeson = 5.2792;
  for (ib = 0; ib < 2; ++ib) {
    // -- Breco
    if (typB[ib] >= seedMode+2) {
      int ibk = ipB[ib];
      massBmeson = massMc[ibk];
      mk4Vector(pMcBreco, pMc[ibk], thetaMc[ibk], phiMc[ibk], massMc[ibk]);
    }    
    // -- recoil
    if ((typB[ib] == 1) || (typB[ib] == 2)) {
      int ibk = ipB[ib];
      mk4Vector(pMcRecoil, pMc[ibk], thetaMc[ibk], phiMc[ibk], massMc[ibk]);
    }
  }

  // check whether two semileptonic B decays are present. If vub=vcb=1, analyze 
  // the vub one, if vub=2 || vcb=2 the first one.

  int indexjump[2] = {0,0};
  if(typB[0] == 1 && typB[1] == 2) indexjump[1] = 1; 
  if(typB[1] == 1 && typB[0] == 2) indexjump[0] = 1; 
  if(typB[0] == 1 && typB[1] == 1) indexjump[1] = 1; 
  if(typB[0] == 2 && typB[1] == 2) indexjump[1] = 1; 

  // -- Recoil
  for (ib = 0; ib < 2; ++ib) {
    if (typB[ib] > 2 || indexjump[ib]) continue;
    
    //    cout<<    ipB[ib] << " <- B index ; nMC particles ->  "<< nMc << endl;
    for (int imc = 0; imc < nMc; ++imc) {

      if (mothMc[imc]-1 == ipB[ib]) {
	mk4Vector(p, pMc[imc], thetaMc[imc], phiMc[imc], massMc[imc]);
	//NO!	p.Boost(boostVector); 
	
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
      //Energy of true Klongs
      double TrueKlEnergy;
      if(TMath::Abs(idMc[imc]) == 130) {
	TrueKlEnergy = p.E();
	((TH1D*)gDirectory->Get("TrueKlEn"))->Fill(TrueKlEnergy, 1.);
	if(TrueKlEnergy>0) ((TH1D*)gDirectory->Get("TrueKlEnNz"))->Fill(TrueKlEnergy, 1.);
	//	cout<<TrueKlEnergy<<" :: Energy od true KL"<<endl;
      }
    }
  }
  
  fKplus   = kPlus();
  //Correction of 2factor missing
  fQ2Gen   = 2.*p4MissGen*p4LeptonGen;

  TVector3 boostVector(-pMcRecoil.Vect());  

  if (fVub + fVcb > 0) {
    boostVector.SetMag(boostVector.Mag()/pMcRecoil.E());
  }
  TLorentzVector lcms = p4LeptonGen;
  lcms.Boost(boostVector); 
  fPcmsGen = lcms.Vect().Mag();
  fTcmsGen = lcms.Theta();
  fFcmsGen = lcms.Phi();
  fEcmsGen = lcms.E();

  fPxhadGen = p4XhadGen.Vect().Mag();
  fTxhadGen = p4XhadGen.Theta();
  fFxhadGen = p4XhadGen.Phi();
  fExhadGen = p4XhadGen.E();
  fMxhadGen = p4XhadGen.Mag();
  ftLepG = p4LeptonGen.E();
  /*
    Added variables to study pW+EW
  */
  TLorentzVector dovG(0., 0., 0., 0.);
  TLorentzVector lepdovG = p4LeptonGen;
  TLorentzVector nudovG = p4MissGen;
  lepdovG.Boost(boostVector); 
  nudovG.Boost(boostVector);
  dovG = lepdovG + nudovG;
  fEwPwG = dovG.E() + dovG.P(); 
  /*
    Finished adding new variables
  */


  //  fQ2Gen    = 2.*p4Neutrino*p4LeptonLab;

  fGwCiuc=2*(fExhadGen/massBmeson);
  fGxCiuc=2*(ftLepG/massBmeson);
  double tmpCiuc;
  tmpCiuc = pow((1-pow((fMxhadGen/fExhadGen),2)),0.5);
  fGcsiCiuc= 2*tmpCiuc/(1+tmpCiuc); 

  if (fVub + fVcb > 0) {
    ((TH1D*)gDirectory->Get("h77000"))->Fill(fNLeptonMC, 1.);
    if (fNLeptonMC == 1) {
      ((TH1D*)gDirectory->Get("h121000"))->Fill(fPcmsGen, 1.);
      ((TH1D*)gDirectory->Get("h123000"))->Fill(fMxhadGen, 1.);

      if (fVcb == 1) {
	((TH1D*)gDirectory->Get("h121005"))->Fill(fPcmsGen, 1.);
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
Double_t VubAnalysisCode::kPlus() {
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
void VubAnalysisCode::initRest() {
  DR = 57.29577951;
  PCBIN = 30;
  PCMAX = 3.0;
  PLBIN = 40;
  PLMAX = 4.0;
  FBIN = 36;
  B0LUND = 511;
  CHBLUND = 521;

  for (int i = 0; i < 100; ++i) fReturnLog[i] = 0; 

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
  BMASS   = 5.2792;
  BQMASS  = 4.800;
  A0      = 1.29;

  fOptBlind = fOptGammas = fOptCategories = fOptKlongs = fOptMakeEventList 
    = fOptPidKilling = fOptPidKillingKaon = fOptPidKillingEl = fOptPidKillingMu 
    = fOptSmearTracks = fOptSmearNeut = fOptScaleKlongs = 0;

  fOptFilterK0sPi0Pi0 = 1010000;// run number before which we need to apply the Ks bug fix

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

  TLABLO = 23.5;  // BAD 90 
  TLABHI = 135.8;

  PCMSLO = 1.0;
  PLABLO = 0.5;
  NLEPTON = 1;
  ELMOMLO = 1.0;
  MUMOMLO = 1.0;
  KAMOMLO = 0.1;
  //  IDEL = IDMU = IDKA = 16; 
  IDEL = 32;
  IDMU = 8; 
  IDKA = 4; 

  MM2LO = -9999.;

  fParametrization = 0; // default: old parametrization for (CT,AS)

  REQTOTALCHARGE = 1.;
  REQCHARGECORR  = kTRUE; 

  KSPIPLO = KSPIZLO = 0.2;
  KSPIPHI = KSPIZHI = 0.8;
  KSPIPRLO = -1.;
  PTLO = 0.000;
 
  GAMMAELO = 0.080;
  GAMMAEHI = 100.0;
  
  fRunRange = TString("undefined"); 
  fHistFile = 0;
  fDump = 0;
  fVerbose = 0;
  fIsMakeParam = kFALSE;

  Int_t i(0);
  for (i = 0; i < 100; ++i) { 
    recEl[i] = recMu[i] = recKa[i] = 0;
    kshortLockTrk[i] = 0; 
    kshortLockGam[i] = 0; 
    goodKshort[i] = 0;
  }

  // -- Get default Trk and Pid tables
  SMEARELPIDTABLES = 0;
  SMEARMUPIDTABLES = 0;
  SMEARKAPIDTABLES = 0;
  SMEARELMISTABLES = 0;
  SMEARMUMISTABLES = 0;
  SMEARKAMISTABLES = 0;
  sprintf(TRKTABLES, "TrkTables.txt");
  //  getTrkTables();
  SMEARTRKPX = SMEARTRKPY = 1.3;
  SIGMANEUT = 0.02;
  SHIFTNEUT = -0.0075;
  
  TRACKSELECTION = PHOTONSELECTION = 0; 
  DOTRACKKILLING = 0; 
  SLOWPIONKILL = 0.016;
  TRACKKILL = 0.013;

  sprintf(PIDTABLES, "PidTables.txt");
  //  getPidTables();

  sprintf(WEIGHTNEU, "tableneu.dat");

  sprintf(WEIGHTCHG, "tablechg.dat");

  int therandom(0);
  Dvar = new recoilDSys("ddecay.table",therandom,2);
  Bsem = new recoilDSys(therandom);
  DOBDECWEIGHT = DODDECWEIGHT = 1;
 
}

// ----------------------------------------------------------------------
void VubAnalysisCode::Init(TTree *tree, int isMC) {
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
   fChain->SetBranchAddress("runNumber",&runNumber);
   fChain->SetBranchAddress("platform",&platform);
   fChain->SetBranchAddress("partition",&partition);
   fChain->SetBranchAddress("upperID",&upperID);
   fChain->SetBranchAddress("lowerID",&lowerID);
   //just for Killing
//   fChain->SetBranchAddress("lower",&lower);
//   fChain->SetBranchAddress("upper",&upper);
//   fChain->SetBranchAddress("run",&run); 
//   fChain->SetBranchAddress("pur",&pur);
//   fChain->SetBranchAddress("mode",&mode);

   fChain->SetBranchAddress("primVtxX",&primVtxX);
   fChain->SetBranchAddress("primVtxY",&primVtxY);
   fChain->SetBranchAddress("primVtxZ",&primVtxZ);
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
   //MC block
   if(isMC) { 
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
   if(isMC)   fChain->SetBranchAddress("MCB0",MCB0);
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
   if(isMC)   fChain->SetBranchAddress("MCChB",MCChB);
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
   if(isMC)   fChain->SetBranchAddress("MCDstar",MCDstar);
   fChain->SetBranchAddress("d1DstarIndex",d1DstarIndex);
   fChain->SetBranchAddress("d1DstarLund",d1DstarLund);
   fChain->SetBranchAddress("d2DstarIndex",d2DstarIndex);
   fChain->SetBranchAddress("d2DstarLund",d2DstarLund);
   fChain->SetBranchAddress("nDstarBS",&nDstarBS);
   fChain->SetBranchAddress("massDstarBS",massDstarBS);
   fChain->SetBranchAddress("chi2DstarBS",chi2DstarBS);
   fChain->SetBranchAddress("dofDstarBS",dofDstarBS);
   //     fChain->SetBranchAddress("spixDstarBS",spixDstarBS);
   //     fChain->SetBranchAddress("spiyDstarBS",spiyDstarBS);
   //     fChain->SetBranchAddress("spizDstarBS",spizDstarBS);
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
   if(isMC)   fChain->SetBranchAddress("MCDstar0",MCDstar0);
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
   if(isMC)   fChain->SetBranchAddress("MCD0",MCD0);
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
   if(isMC)   fChain->SetBranchAddress("MCChD",MCChD);
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
   if(isMC)   fChain->SetBranchAddress("MCKs",MCKs);
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
   if(isMC)   fChain->SetBranchAddress("MCPi0",MCPi0);
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
   if(isMC)   fChain->SetBranchAddress("MCGConv",MCGConv);
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
     if(isMC)   fChain->SetBranchAddress("MCDalitz",MCDalitz);
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
     if(isMC)   fChain->SetBranchAddress("MCJpsi",MCJpsi);
     fChain->SetBranchAddress("d1JpsiIndex",d1JpsiIndex);
     fChain->SetBranchAddress("d1JpsiLund",d1JpsiLund);
     fChain->SetBranchAddress("d2JpsiIndex",d2JpsiIndex);
     fChain->SetBranchAddress("d2JpsiLund",d2JpsiLund);
   }   
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
   if (fNewFormat) {
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
   if (isMC) {
     fChain->SetBranchAddress("idTrk",idTrk);
     fChain->SetBranchAddress("IndexTrk",IndexTrk);
     fChain->SetBranchAddress("IndexNtTrk",IndexNtTrk);
   }
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
   if(isMC){
     fChain->SetBranchAddress("idGam",idGam);
     fChain->SetBranchAddress("IndexGam",IndexGam);
     fChain->SetBranchAddress("IndexNtGam",IndexNtGam);
   }
   fChain->SetBranchAddress("B0RecGam",B0RecGam);
   fChain->SetBranchAddress("chBRecGam",chBRecGam);

   Notify();
}

// ----------------------------------------------------------------------
void VubAnalysisCode::initVariables() {
  TLorentzVector p4Null(0., 0., 0., 0.);
  fRunnumber = fLower = fUpper = -99;
  fEvtW8 = 1.; 
  fWithKKbar = 0; 
  
  // -- Cuts
  fOneLepton = kFALSE;
  fGoodPRMM2 = fGoodMM2 = fGoodChargeCorr = fGoodChargeCons = fGoodEvent = fLowMx = kFALSE;

  fGoodAccLepton = fGoodAccElectron = fGoodAccMuon = kFALSE; 
  // need fGoodAccERS? -> that' s "faoRS"

  fGoodLepton    = fGoodElectron    = fGoodMuon    = kFALSE;  
  fGoodRS        = fGoodERS         = fGoodMRS     = kFALSE; 
  fGoodWS        = fGoodEWS         = fGoodMWS     = kFALSE; 
  fGoodRF        = fGoodERF         = fGoodMRF     = kFALSE; 
  fGoodWF        = fGoodEWF         = fGoodMWF     = kFALSE; 

  faoLepton = faoElectron = faoMuon = kFALSE; 
  faoRS     = faoERS      = faoMRS  = kFALSE; 
  faoWS     = faoEWS      = faoMWS  = kFALSE; 
  faoRF     = faoERF      = faoMRF  = kFALSE; 
  faoWF     = faoEWF      = faoMWF  = kFALSE; 

  faoPRMM2   = faoPRMM2E   = faoPRMM2M   = kFALSE; 
  faoPRMM2RS = faoPRMM2ERS = faoPRMM2MRS = kFALSE; 
  faoPRMM2WS = faoPRMM2EWS = faoPRMM2MWS = kFALSE; 
  faoPRMM2RF = faoPRMM2ERF = faoPRMM2MRF = kFALSE; 
  faoPRMM2WF = faoPRMM2EWF = faoPRMM2MWF = kFALSE; 

  faoMM2   = faoMM2E   = faoMM2M   = kFALSE; 
  faoMM2RS = faoMM2ERS = faoMM2MRS = kFALSE; 
  faoMM2WS = faoMM2EWS = faoMM2MWS = kFALSE; 
  faoMM2RF = faoMM2ERF = faoMM2MRF = kFALSE; 
  faoMM2WF = faoMM2EWF = faoMM2MWF = kFALSE; 
  
  faoChargeCons   = faoChargeConsE   = faoChargeConsM   = kFALSE; 
  faoChargeConsRS = faoChargeConsERS = faoChargeConsMRS = kFALSE;
  faoChargeConsWS = faoChargeConsEWS = faoChargeConsMWS = kFALSE; 
  faoChargeConsRF = faoChargeConsERF = faoChargeConsMRF = kFALSE;
  faoChargeConsWF = faoChargeConsEWF = faoChargeConsMWF = kFALSE;

    

  // -- MC Truth 
  p4XhadGen = p4MissGen = p4RecoilGen = p4LeptonGen = p4Null;

  fQ2Gen = fKplus = fMxhadGen = fPxhadGen = fTxhadGen = fFxhadGen = fExhadGen = -99.;
  fMxhadRes = fMxhadfitRes = -99.;
  fPcmsGen = fTcmsGen = fFcmsGen = fEcmsGen = -99.;
  fB1Index = fB2Index = fBVxb = fBVxbTyp = -99;
  fMxCategory = 0;
  fDpi = fDpiz = fDk = fDks = fDkmiss= fDlep = fDCfDs = fD0CfDs = 0;
  fDkl = fDkspiopio = fDkspipi = 0;
  // -- Reco
  p4Brecoil = p4Breco = p4BrecoNC = p4BrecoGen = p4BrecoGenSmear = p4Null;

  lundB = fSeedMode = fBmode = -999999;
  
  mesSideband = mesSignalband = signalBox = vubDepleted = kFALSE;
  fNLeptonMC = 0; 
  genLeptonCharge = fBrecoCharge; 
  fBrecoFlavor = 0;

  fVub =  fVcb =  fOther =  0;
  fBrecoMc = -99;
  fNlep = fNLepton = fNEl = fNMu = fNKshort = fNKp = 0;

  fPurity = fIntPurity = -99.;
  fMes = fDeltaE = fBmass = -99.;
 
  fQtot = -99; 
  fNneutrals = fNtracks = -99.;
  fRecoilCharge = fLeptonCharge = fRecoilTrkMult = fRecoilNutMult = 0;
  fRecoilNutMult80_160 = fRecoilNutMult160_320 = fRecoilNutMultfromB = fRecoilNutMultfromB80_160 = fRecoilNutMultfromB160_320 = 0;
  fIsPrompt = fIsCascade = fElectron = fMuon = kFALSE;  
  fEffCat = 0;

  fGammaMax = fQ2 = fQ2NC = fQ2Fit = fQ2Res = -99.;
  fPNu = fTNu = fFNu = fEmiss = fMM2 = -99.;
  fPNuNC = fTNuNC = fFNuNC = fEmissNC = fMM2NC = -99.;
  fPlab = fTlab = fFlab = fPcms = fTcms = fFcms = fEcms = fPups = -99.;
  fPxhadchg = fMxhadchg = fTxhadchg = fExhadchg = fFxhadchg = -99.;
  fPxhad = fTxhad = fFxhad = fExhad = fMxhad = -99.;
  fPrecoil = fTrecoil = fFrecoil = fErecoil = fMrecoil = -99.;
  fBmassfit = fMxhadfit = fMM2fit = fQ2Fit = 0.;
  fPAllev = fTAllev = fFAllev = fEAllev = fMAllev =-99.; 
  fProbChi2 = fChi2 = -99.;

  for (int i=0; i<100; i++){
    weightTrk[i] = 1;
    weightNeu[i] = 1;
    
    loopTrack[i] = 0; 
  }

  for (int i=0; i<80; i++) ifromBGam[i]=0;

  // -- Some event quantities
  fRunnumber = runNumber;
  fLower = lowerID; 
  fUpper = upperID; 



}

void VubAnalysisCode::selectPhotons() {
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
// additional fix to Ks bug: kill intermediate energy photons	
// according to k(e)=1+gaus(0.68,0.72,0.3)*poli(8)
    if(fIsMC && runNumber< fOptFilterK0sPi0Pi0&& goodPhoton[i] == 1){
      double coeff[9]={ .792903,4.80598,-32.5345,95.5484,-145.436,123.012,-58.6776,14.8144,-1.54217};
      Int_t iKs=isKs2Pi0Dau(i);
      if(iKs>=0){
	double killFac=1+0.68*exp(-0.5*pow((energyGam[i]-0.72)/0.3,2));
	double corr = 0;
	for (Int_t pol=0;pol<9;pol++)// 8th prder polinomial corerction
	  corr+=coeff[pol]*pow(energyGam[i],pol);

	killFac/=corr;
	double xrand=gRandom->Rndm() ;
	if(xrand > 1./killFac)  goodPhoton[i] = 0;
	//	if(corr<0.8)cout << "low corr: "<<corr<<"for E= "<<energyGam[i]<<endl;
      }
      
    }
    if (fVerbose) cout << i << "  " << energyGam[i] <<" "<<thetaGam[i]<<" "<<phiGam[i]<<p4Gam.T()<<" "<<lMomGam[i]<<" "<<s9s25Gam[i]<<" "<<splitOffGam[i]<< "  " << goodPhoton[i] << endl;
  }
}

Int_t VubAnalysisCode::isKs2Pi0Dau(Int_t iGam){
  if(IndexGam[iGam]<=0)return -1;
  Int_t Imoth(mothMc[IndexGam[iGam]-1]-1);
  if(Imoth<0)return -1;
  if(idMc[Imoth]!=111)return -1;
  Int_t Igrand(mothMc[Imoth]-1);
  if(Igrand<0)return -1;
  if(idMc[Igrand]!=310)return -1;
  return Igrand;
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
void VubAnalysisCode::selectTracks() {
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
    tlenTrack[i] = ndchTrack[i] = nsvtTrack[i] = c2nTrack[i]  = dcaTrack[i] = dcazTrack[i] = -99.; 
    goodTrack[i] = goodHadron[i] = goodChargedKaon[i] = goodPion[i] = 0; 
    gtvl = gtl = acc = 0;
    pt = momentumTrk[i]*sin(thetaTrk[i]);

    dcaz = zPocaTrk[i] - beamSZ;
    dcazTrack[i] = dcaz; 
    dca  = TMath::Sqrt((xPocaTrk[i]-beamSX)*(xPocaTrk[i]-beamSX) + (yPocaTrk[i]-beamSY)*(yPocaTrk[i]-beamSY));
    dcaTrack[i]  = (dca > 1.e-4 ? TMath::Log10(dca) : -3.9); 
    c2nTrack[i] = (ntdofTrk[i] > 0 ? tChi2Trk[i]/ntdofTrk[i] : -1); 
    tlenTrack[i] = tLenTrk[i]; 
    ndchTrack[i] = ndchTrk[i]; 
    nsvtTrack[i] = nsvtTrk[i]; 

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
	  //  	  if (fVerbose) cout << "Killing track: " << i << " p = " << momentumTrk[i] 
	  //  			     << " rand = " << rand << " < " << TRACKKILL 
	  //  			     << " -- killed " << killedTracks << " tracks so far "
	  //  			     << endl;
	  if (0 == killedTracks%1000) cout << "killed " << killedTracks << " tracks so far (in the entire event)" << endl;
	}      
      }
      if (DOTRACKKILLING & 2) {
        if ((pt < 0.200) && (rand < SLOWPIONKILL) && (goodTrack[i] == 1)) {
          goodTrack[i] = 0; 
	  ++killedSlowPions; 
	  //            cout << "Killing slow pion: " << i << " p = " << momentumTrk[i] 
	  //  	       << " rand = " << rand << " < " << SLOWPIONKILL
	  //  	       << " -- killed " << killedSlowPions << " slow pions so far "
	  //  	       << endl;
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
