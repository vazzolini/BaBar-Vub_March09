#include <fstream.h>

#include "b2uQA.hh"
#include <TH2.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "../RecoilAnalysis/recoilDSys.hh"

b2uQA::b2uQA(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
     cout << "no tree provided or found" << endl;
     return;
   }
   Init(tree);
}

b2uQA::b2uQA(const char *filename, int Sys, TString cutfile, int reweight)
{

  fReW = reweight;

  if (fReW == 1) cout << "Applying re-weighting to MC" << endl;

  tl = new TLatex;

  TString fileString(filename);
  if (fileString.Contains(".root")) {
    cout << "Reading file " << filename << endl;
    fFilename = TString(filename);
    TFile *f = new TFile(filename);
    fChain = (TTree*)f->Get("events");
  } else{
    cout << "Reading chain " << filename << endl;
    fFilename = TString(filename);
    fFilename.ReplaceAll("*", "");
    fFilename += TString(".root");
    TChain *a = new TChain("events");
    a->Add(filename);
    fChain = a;
  }

  if (fileString.Contains("data")) {
    fMC = 0;
  } else{
    fMC = 1;
  }

  Init(fChain);

  QTOTCUT = 0.5; MM2CUTHI = 0.5; MM2CUTLO = -1.0; 
  EMPCUTHI = 10000.; EMPCUTLO = -10000.; CTMCUTLO = -10.95;
  CTMCUTHI = 10.95; PMISSCUT = -1.; PCMSTRKLO = -1.;
  PRMM1CUT = -20.; PRMM2CUT =-3.; PRMM3CUT = -20.; PCMSCUT = 1.0;

  readCuts(cutfile);

   //Dvar = new recoilDSys("ddecay.table",0,2);
   //Bsem = new recoilDSys(0);

  TRandom random(0);
  random.SetSeed(0); 
  int therandom = 0;
  dImode = 2;
  if(Sys > 0) {
    dImode = Sys;
    therandom = random.Rndm() * 1000000;
    if(Sys > 2) {
      therandom = Sys;
      Sys = 2;
    }
    if(Sys == 2) {
      //if(unfB){
      //  Dvar = new recoilDSys("ddecay.table",0,Sys);
      //  Dvar->recoilDSys3(therandom);
      //}
      //else
      Dvar = new recoilDSys("ddecay.table",therandom,Sys);
    } else {
      //if(unfB){
      //  Dvar = new recoilDSys("dIdecay.table",0,Sys);
      //  Dvar->recoilDSys3(therandom);
      //}
      //else 
      Dvar = new recoilDSys("dIdecay.table",therandom,Sys);
    }
    //if(unfB){
    //  Bsem = new recoilDSys(0); 
    //  Bsem->recoilDSys2(therandom);
    //}
    //else
    Bsem = new recoilDSys(therandom);
  } else {
    Dvar = new recoilDSys("ddecay.table",therandom,2);
    Bsem = new recoilDSys(therandom);
  }
  //  cout<<"therandom::  "<<therandom<<"   Dvar::   "<<Dvar<<"   Bsem::  "<<Bsem<<endl;

}

// ----------------------------------------------------------------------
b2uQA::~b2uQA() {
  if (fHistFile) {
      fHistFile->cd();
      fHistFile->Write();
      fHistFile->Close();
      delete fHistFile;
  }
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

void b2uQA::readCuts(TString cutfile){

  cout << "======================================="  << endl;
  cout << "==> readCuts> reading cut file: " << cutfile.Data() << endl;
  cout << "---------------------------------------"  << endl;

  char  buffer[200];
  sprintf(buffer, "%s", cutfile.Data());
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
    if (!strcmp(CutName, "qtot"))     {
      QTOTCUT = CutValue;     ok = 1;
      cout << "qtot:            " << QTOTCUT << endl;
    }
    if (!strcmp(CutName, "mm2hi"))     {
      MM2CUTHI = CutValue;     ok = 1;
      cout << "mm2hi:            " << MM2CUTHI << endl;
    }
    if (!strcmp(CutName, "mm2lo"))     {
      MM2CUTLO = CutValue;     ok = 1;
      cout << "mm2lo:            " << MM2CUTLO << endl;
    }
    if (!strcmp(CutName, "empmhi"))     {
      EMPCUTHI = CutValue;     ok = 1;
      cout << "empmhi:            " << EMPCUTHI << endl;
    }
    if (!strcmp(CutName, "empmlo"))     {
      EMPCUTLO = CutValue;     ok = 1;
      cout << "empmlo:            " << EMPCUTLO << endl;
    }
    if (!strcmp(CutName, "ctmlo"))     {
      CTMCUTLO = CutValue;     ok = 1;
      cout << "ctmlo:            " << CTMCUTLO << endl;
    }
    if (!strcmp(CutName, "ctmhi"))     {
      CTMCUTHI = CutValue;     ok = 1;
      cout << "ctmhi:            " << CTMCUTHI << endl;
    }
    if (!strcmp(CutName, "pmisscut"))     {
      PMISSCUT = CutValue;     ok = 1;
      cout << "pmisscut:            " << PMISSCUT << endl;
    }
    if (!strcmp(CutName, "pcmstrklo"))     {
      PCMSTRKLO = CutValue;     ok = 1;
      cout << "pcmstrklo:            " << PCMSTRKLO << endl;
    }
    if (!strcmp(CutName, "pcmscut"))     {
      PCMSCUT = CutValue;     ok = 1;
      cout << "pcmscut:            " << PCMSCUT << endl;
    }
    if (!strcmp(CutName, "prmm1cut"))     {
      PRMM1CUT = CutValue;     ok = 1;
      cout << "prmm1cut:            " << PRMM1CUT << endl;
    }
    if (!strcmp(CutName, "prmm2cut"))     {
      PRMM2CUT = CutValue;     ok = 1;
      cout << "prmm2cut:            " << PRMM2CUT << endl;
    }
    if (!strcmp(CutName, "prmm3cut"))     {
      PRMM3CUT = CutValue;     ok = 1;
      cout << "prmm3cut:            " << PRMM3CUT << endl;
    }
    

    if (ok == 0)  cout << "==> readCuts() Error: Don't know about variable " << CutName << endl;
  }

  cout << "---------------------------------------"  << endl;

}

Int_t b2uQA::GetEntry(Int_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Int_t b2uQA::LoadTree(Int_t entry)
{
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

void b2uQA::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses of the tree
   // will be set. It is normaly not necessary to make changes to the
   // generated code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running with PROOF.

   // Set branch addresses
   if (tree == 0) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run",&run);
   fChain->SetBranchAddress("mes",&mes);
   fChain->SetBranchAddress("intpur",&intpur);
   fChain->SetBranchAddress("brecoflav",&brecoflav);
   fChain->SetBranchAddress("brecocharge",&brecocharge);
   fChain->SetBranchAddress("xcharge",&xcharge);
   fChain->SetBranchAddress("mxhad",&mxhad);
   fChain->SetBranchAddress("mxtrk",&mxtrk);
   fChain->SetBranchAddress("mxnut",&mxnut);
   fChain->SetBranchAddress("mxhadfit",&mxhadfit);
   fChain->SetBranchAddress("chi2",&chi2);
   fChain->SetBranchAddress("probchi2",&probchi2);
   fChain->SetBranchAddress("lcharge",&lcharge);
   fChain->SetBranchAddress("pcms",&pcms);
   fChain->SetBranchAddress("tcms",&tcms);
   fChain->SetBranchAddress("tlab",&tlab);
   fChain->SetBranchAddress("plab",&plab);
   fChain->SetBranchAddress("mm2",&mm2);
   fChain->SetBranchAddress("q2fit",&q2fit);
   fChain->SetBranchAddress("pmiss",&pmiss);
   fChain->SetBranchAddress("tmiss",&tmiss);
   fChain->SetBranchAddress("emiss",&emiss);
   fChain->SetBranchAddress("nle",&nle);
   fChain->SetBranchAddress("nel",&nel);
   fChain->SetBranchAddress("nmu",&nmu);
   fChain->SetBranchAddress("nchg",&nchg);
   fChain->SetBranchAddress("nneu",&nneu);
   fChain->SetBranchAddress("nks",&nks);
   fChain->SetBranchAddress("nkz",&nkz);
   fChain->SetBranchAddress("nkl",&nkl);
   fChain->SetBranchAddress("nkp",&nkp);
   fChain->SetBranchAddress("deltam",&deltam);
   fChain->SetBranchAddress("mm1pr",&mm1pr);
   fChain->SetBranchAddress("mm2pr",&mm2pr);
   fChain->SetBranchAddress("mm3pr",&mm3pr);
   fChain->SetBranchAddress("oa1",&oa1);
   fChain->SetBranchAddress("oa2",&oa2);
   fChain->SetBranchAddress("oa3",&oa3);
   fChain->SetBranchAddress("pcmstrklo",&pcmstrklo);
   fChain->SetBranchAddress("vub",&vub);
   fChain->SetBranchAddress("vcb",&vcb);
   fChain->SetBranchAddress("other",&other);
   fChain->SetBranchAddress("qb",&qb);
   fChain->SetBranchAddress("kplus",&kplus);
   fChain->SetBranchAddress("ecmsgen",&ecmsgen);
   fChain->SetBranchAddress("tcmsgen",&tcmsgen);
   fChain->SetBranchAddress("fcmsgen",&fcmsgen);
   fChain->SetBranchAddress("mxhadgen",&mxhadgen);
   fChain->SetBranchAddress("pxhadgen",&pxhadgen);
   fChain->SetBranchAddress("q2Gen",&q2Gen);
   fChain->SetBranchAddress("ctvgen",&ctvgen);
   fChain->SetBranchAddress("ctlgen",&ctlgen);
   fChain->SetBranchAddress("chigen",&chigen);
   fChain->SetBranchAddress("vpgen",&vpgen);
   fChain->SetBranchAddress("vtgen",&vtgen);
   fChain->SetBranchAddress("vfgen",&vfgen);
   fChain->SetBranchAddress("vpcmsgen",&vpcmsgen);
   fChain->SetBranchAddress("vtcmsgen",&vtcmsgen);
   fChain->SetBranchAddress("vfcmsgen",&vfcmsgen);
   fChain->SetBranchAddress("Gvxbtyp",&Gvxbtyp);
   fChain->SetBranchAddress("GfDpi",&GfDpi);
   fChain->SetBranchAddress("GfDk",&GfDk);
   fChain->SetBranchAddress("GfDks",&GfDks);
   fChain->SetBranchAddress("GfDpiz",&GfDpiz);
   fChain->SetBranchAddress("GfDlep",&GfDlep);
   Notify();
}

Bool_t b2uQA::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. Typically here the branch pointers
   // will be retrieved. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed.

   // Get branch pointers
   b_run = fChain->GetBranch("run");
   b_mes = fChain->GetBranch("mes");
   b_intpur = fChain->GetBranch("intpur");
   b_brecoflav = fChain->GetBranch("brecoflav");
   b_brecocharge = fChain->GetBranch("brecocharge");
   b_xcharge = fChain->GetBranch("xcharge");
   b_mxhad = fChain->GetBranch("mxhad");
   b_mxtrk = fChain->GetBranch("mxtrk");
   b_mxnut = fChain->GetBranch("mxnut");
   b_mxhadfit = fChain->GetBranch("mxhadfit");
   b_chi2 = fChain->GetBranch("chi2");
   b_probchi2 = fChain->GetBranch("probchi2");
   b_lcharge = fChain->GetBranch("lcharge");
   b_pcms = fChain->GetBranch("pcms");
   b_tcms = fChain->GetBranch("tcms");
   b_tlab = fChain->GetBranch("tlab");
   b_plab = fChain->GetBranch("plab");
   b_mm2 = fChain->GetBranch("mm2");
   b_q2fit = fChain->GetBranch("q2fit");
   b_pmiss = fChain->GetBranch("pmiss");
   b_tmiss = fChain->GetBranch("tmiss");
   b_emiss = fChain->GetBranch("emiss");
   b_nle = fChain->GetBranch("nle");
   b_nel = fChain->GetBranch("nel");
   b_nmu = fChain->GetBranch("nmu");
   b_nchg = fChain->GetBranch("nchg");
   b_nneu = fChain->GetBranch("nneu");
   b_nks = fChain->GetBranch("nks");
   b_nkz = fChain->GetBranch("nkz");
   b_nkl = fChain->GetBranch("nkl");
   b_nkp = fChain->GetBranch("nkp");
   b_wdeltam = fChain->GetBranch("deltam");
   b_mm1pr = fChain->GetBranch("mm1pr");
   b_mm2pr = fChain->GetBranch("mm2pr");
   b_mm3pr = fChain->GetBranch("mm3pr");
   b_oa1 = fChain->GetBranch("oa1");
   b_oa2 = fChain->GetBranch("oa2");
   b_oa3 = fChain->GetBranch("oa3");
   b_pcmstrklo = fChain->GetBranch("pcmstrklo");
   b_vub = fChain->GetBranch("vub");
   b_vcb = fChain->GetBranch("vcb");
   b_other = fChain->GetBranch("other");
   b_qb = fChain->GetBranch("qb");
   b_fkplus = fChain->GetBranch("kplus");
   b_ecmsgen = fChain->GetBranch("ecmsgen");
   b_tcmsgen = fChain->GetBranch("tcmsgen");
   b_fcmsgen = fChain->GetBranch("fcmsgen");
   b_mxhadgen = fChain->GetBranch("mxhadgen");
   b_pxhadgen = fChain->GetBranch("pxhadgen");
   b_q2Gen = fChain->GetBranch("q2Gen");
   b_ctvgen = fChain->GetBranch("ctvgen");
   b_ctlgen = fChain->GetBranch("ctlgen");
   b_chigen = fChain->GetBranch("chigen");
   b_vpgen = fChain->GetBranch("vpgen");
   b_vtgen = fChain->GetBranch("vtgen");
   b_vfgen = fChain->GetBranch("vfgen");
   b_vpcmsgen = fChain->GetBranch("vpcmsgen");
   b_vtcmsgen = fChain->GetBranch("vtcmsgen");
   b_vfcmsgen = fChain->GetBranch("vfcmsgen");
   b_Gvxbtyp = fChain->GetBranch("Gvxbtyp");
   b_GfDpi = fChain->GetBranch("GfDpi");
   b_GfDk = fChain->GetBranch("GfDk");
   b_GfDks = fChain->GetBranch("GfDks");
   b_GfDpiz = fChain->GetBranch("GfDpiz");
   b_GfDlep = fChain->GetBranch("GfDlep");

   return kTRUE;
}


// ----------------------------------------------------------------------
void b2uQA::bookHist(int mode) {
  fFilename = TString(Form("qa-%d-", mode)) + fFilename;
  cout << "Dumping output into " <<  fFilename.Data() << endl;
  fHistFile = new TFile(fFilename.Data(), "RECREATE");
  fHistFile->cd();

  fHmes = new TH1D("mes", "mes normalization", 40, 5.2, 5.3); fHmes->Sumw2();

  fH800  = new sHist(800,  "pcms", 18, 0., 3.0);  fH800->setup(&pcms, &mes);
  fH900  = new sHist(900,  "pcms", 18, 0., 3.0);  fH900->setup(&pcms, &mes);
  fH1000 = new sHist(1000, "pcms", 18, 0., 3.0);  fH1000->setup(&pcms, &mes);
  fH1100 = new sHist(1100, "tlab", 26, 0., 2.6);  fH1100->setup(&tlab, &mes);

  fH2000 = new sHist(2000, "mxhad", 20, 0., 4.);    fH2000->setup(&mxhad, &mes);
  fH2100 = new sHist(2100, "mxhadfit", 20, 0., 4.); fH2100->setup(&mxhadfit, &mes);
  fH2200 = new sHist(2200, "mxtrk", 15, 0., 3.);    fH2200->setup(&mxtrk, &mes);
  fH2300 = new sHist(2300, "mxnut", 15, 0., 3.);    fH2300->setup(&mxnut, &mes);

  fH3000 = new sHist(3000, "mm2", 30, -5., 10.);         fH3000->setup(&mm2, &mes);
  fH3100 = new sHist(3100, "pmiss", 12, 0., 4.);         fH3100->setup(&pmiss, &mes);
  fH3200 = new sHist(3200, "tmiss", 17, 0., 3.40);       fH3200->setup(&tmiss, &mes);
  fH3300 = new sHist(3300, "Q2", 22, 0., 22.);           fH3300->setup(&q2fit, &mes);
  fH3400 = new sHist(3400, "emiss-pmiss", 35, -1., 2.5); fH3400->setup(&depmiss, &mes);
  fH3500 = new sHist(3500, "ctmiss", 20, -1., 1.0);      fH3500->setup(&dctmiss, &mes);

  fH4000 = new sHist(4000, "nchg", 12, 0., 12.); fH4000->setup(&dnchg, &mes);
  fH4100 = new sHist(4100, "nneu", 12, 0., 12.); fH4100->setup(&dnneu, &mes);
  fH4200 = new sHist(4200, "qtot", 10,-5.,  5.); fH4200->setup(&dqtot, &mes);
  fH4300 = new sHist(4300, "pcmstrklo", 16, 0.,  0.8); fH4300->setup(&pcmstrklo, &mes);

}

// ----------------------------------------------------------------------
void b2uQA::readHist(int mode) {

  gDirectory->ReadAll();

  fHmes = (TH1D*)gDirectory->Get("mes");

  fH800  = new sHist(800);
  fH900  = new sHist(900);
  fH1000 = new sHist(1000);
  fH1100 = new sHist(1100);

  fH2000 = new sHist(2000);
  fH2100 = new sHist(2100);
  fH2200 = new sHist(2200);
  fH2300 = new sHist(2300);

  fH3000 = new sHist(3000);
  fH3100 = new sHist(3100);
  fH3200 = new sHist(3200);
  fH3300 = new sHist(3300);
  fH3400 = new sHist(3400);
  fH3500 = new sHist(3500);

  fH4000 = new sHist(4000);
  fH4100 = new sHist(4100);
  fH4200 = new sHist(4200);
  fH4300 = new sHist(4300);

}


// ----------------------------------------------------------------------
void b2uQA::makeAll(int mode ) {
  bookHist(mode);
  loop(mode);
}


// ----------------------------------------------------------------------
void b2uQA::loop(int mode) {
  int takeB0(0), takeBp(0), 
    takeEl(0), takeMu(0), 
    takeDepleted(0);

  if (mode/100 == 3) {
    takeBp = 1;
    cout << " B+ ";
    mode -= 300; 
  } else if (mode/100 == 2) {
    takeB0 = 1;
    cout << " B0 ";
    mode -= 200; 
  } else {
    cout << " B ";
  }

  if (mode/30 == 1) {
    takeMu = 1; 
    cout << " muons ";
    mode -= 30; 
  } else if (mode/20 == 1) {
    takeEl = 1;
    cout << " electrons ";
    mode -= 20; 
  } else {
    cout << " leptons ";
  }

  if (mode == 1) {
    takeDepleted = 1;     
    cout << " depleted ";
  } else {
    cout << " enhanced ";
  }
  cout << endl;


  //  for (Int_t jentry = 0; jentry < 20000; jentry++) {
  double w(1.);
  int depleted(0); 
  int prVeto(1);
  int primaryLepton(0);

  //double PRMM1CUT(-20.);
  //double PRMM2CUT(-3.);
  //double PRMM3CUT(-20.);
  //double PCMSCUT(1.0);

  // -- hard mm2
  //double MM2CUTHI(0.5);
  //double MM2CUTLO(-1.0);
  //double EMPCUTHI(100000.15);
  //double EMPCUTLO(-100000.3);

  //double QTOTCUT(0.5);
  //double CTMCUTHI(10.95);
  //double PCMSTRKLO(-0.1);
  //double CTMCUTLO(-10.95);
  //double PMISSCUT(-0.1);

//   // -- loose
//   cout << "==> running LOOSE !!! " << endl;
//   double QTOTCUT(1.0);
//   double MM2CUT(1.0);
//   double CTMCUT(1.1);

  cout << "CUTS: " << endl
       << "QTOTCUT   " << QTOTCUT << endl
       << "MM2CUTHI  " << MM2CUTHI << endl
       << "MM2CUTLO  " << MM2CUTLO << endl
       << "EMPCUTHI  " << EMPCUTHI << endl
       << "EMPCUTLO  " << EMPCUTLO << endl
       << "CTMCUTLO  " << CTMCUTLO << endl
       << "CTMCUTHI  " << CTMCUTHI << endl
       << "PMISSCUT  " << PMISSCUT << endl
       << "PCMSTRKLO " << PCMSTRKLO << endl
    ;
   if (fChain == 0) return;

   Int_t nentries = Int_t(fChain->GetEntries());
   Int_t nevents(0);

   Int_t nbytes = 0, nb = 0;
   for (Int_t jentry=0; jentry<nentries;jentry++) {
      Int_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;


      if (jentry % 50000 == 0) cout << "Event " << jentry << endl;
      
      w = 1.;

      if( fMC && fReW==1){
	//Vcb and oth weights
	//w *= getBrecoWeight(intpur);                     //breco weighting  

	//cout << "    calling bsysw" << endl;
	w *= getBsysweight(Gvxbtyp,vub);                                //Bdec weighting
	//cout << "getBsysweight: " << getBsysweight(Gvxbtyp,vub) << endl;
	//cout << "       w = " << w << endl;

	//cout << "    calling dsysw" << endl;
	//cout << "       GfDpi " << GfDpi << " GfDk " << GfDk << " GfDks " << GfDks << " GfDpiz " << GfDpiz << " GfDlep " << GfDlep << " vub " << vub << endl;
	w *= getDsysweight(GfDpi,GfDk,GfDks,GfDpiz,GfDlep,dImode,vub);  //Ddec weighting
	//cout << "getDsysweight: " << getDsysweight(GfDpi,GfDk,GfDks,GfDpiz,GfDlep,dImode,vub) << endl;
	//cout << "       w = " << w << endl;
      }


      dnchg = double(nchg); 
      dnneu = double(nneu); 
      
      primaryLepton = 0;
      if (brecoflav + lcharge == 0) {
	primaryLepton = 1;
      }
      
      depleted = 0; 
      if ((nks > 0) || (nkp > 0)) {
	depleted = 1;
      }
      
      prVeto = 0; 
      prVeto = ((mm1pr > PRMM1CUT) || (mm3pr > PRMM3CUT) || (mm2pr > PRMM2CUT));   // new PR cuts EJH
      if (prVeto) {
	depleted = 1;
      }
      
      dqtot = brecocharge + xcharge; 
      depmiss = emiss - pmiss;
      dctmiss = TMath::Cos(tmiss);
      
      // -- Cuts for the normalization
      // -----------------------------
      // o PCMSCUT
      if (pcms > 1.0) {
	fHmes->Fill(mes,w);
      }

      // -- Subsample selection
      if (takeB0 && (brecocharge != 0)) continue;
      if (takeBp && (brecocharge == 0)) continue;
      
      if (takeEl && (nel == 0)) continue;
      if (takeMu && (nmu == 0)) continue;

      if (pcms < 0.5) continue;
      if (primaryLepton == 0) continue;
      if (nle > 1) continue;

      fH800->fillHist(w);      

      if (takeDepleted  && (depleted == 0)) continue;
      if ((takeDepleted == 0) && (depleted == 1)) continue;

      fH900->fillHist(w);      

      // -- NB: The cuts employed are
      // ----------------------------
      // o primary lepton (charge correlation) for B+ and B0!
      // o N(lepton) == 1
      
      // And then some more, depending on variable
      // o PCMSCUT
      // o MM2CUT
      // o QTOTCUT
      
    // N(kaons) and partial reco vetos are used to tag "depleted" sample
      
    // -- test to remove the partial reco from the "depleted" sample
      if (prVeto == 1) continue;
      
      // -- Lepton plots
      if (
	  (TMath::Abs(dqtot) <= QTOTCUT) 
	  && (CTMCUTLO < dctmiss ) && (dctmiss < CTMCUTHI)
	  && (pmiss > PMISSCUT)
	  && (pcmstrklo > PCMSTRKLO)
	  && ((EMPCUTLO < depmiss) && (depmiss < EMPCUTHI))
	  && ((MM2CUTLO < mm2 ) && (mm2 < MM2CUTHI))
	  ) {
	fH1000->fillHist(w);
	fH1100->fillHist(w);
      }
      
      
      // -- Hadron system
      if (
	  (pcms > PCMSCUT) 
	  && (CTMCUTLO < dctmiss ) && (dctmiss < CTMCUTHI)
	  && (pmiss > PMISSCUT)
	  && (pcmstrklo > PCMSTRKLO)
	  && (TMath::Abs(dqtot) <= QTOTCUT)
	  && ((EMPCUTLO < depmiss) && (depmiss < EMPCUTHI))
	  && ((MM2CUTLO < mm2 ) && (mm2 < MM2CUTHI))
	  ) {
	fH2000->fillHist(w);
	fH2100->fillHist(w);
	fH2200->fillHist(w);
	fH2300->fillHist(w);
      }
      
      
      // -- missing system
      if (
	  (pcms > PCMSCUT) 
	  && (pcmstrklo > PCMSTRKLO)
	  && (TMath::Abs(dqtot) <= QTOTCUT)
	  ) {
	
	// -- missing mass
	if (
	    (CTMCUTLO < dctmiss ) && (dctmiss < CTMCUTHI)
	    && (pmiss > PMISSCUT)
	    ) {
	  
	  fH3000->fillHist(w);
	  fH3400->fillHist(w);
	}
	
	// -- missing momentum
	if (
	    ((EMPCUTLO < depmiss) && (depmiss < EMPCUTHI))
	    && ((MM2CUTLO < mm2 ) && (mm2 < MM2CUTHI))
	    ) {
	  
	  fH3100->fillHist(w);
	  fH3200->fillHist(w);
	  fH3500->fillHist(w);
	}
      }
      
      
      // -- multiplicities
      if (
	  (pcms > PCMSCUT) 
	  && (CTMCUTLO < dctmiss ) && (dctmiss < CTMCUTHI)
	  && (pmiss > PMISSCUT)
	  && (pcmstrklo > PCMSTRKLO)
	  && ((EMPCUTLO < depmiss) && (depmiss < EMPCUTHI))
	  && ((MM2CUTLO < mm2 ) && (mm2 < MM2CUTHI))
	  ) {
	fH3300->fillHist(w);
	fH4000->fillHist(w);
	fH4100->fillHist(w);
	fH4200->fillHist(w);
      }
      
      nevents++;
   }
}


// ----------------------------------------------------------------------
void b2uQA::show(TH1D *hd, TH1D *h1, TH1D *h2) {

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLabelFont(133, "X");
  gStyle->SetLabelFont(133, "Y");

  TCanvas *c6 = new TCanvas("c6", "c6", 50.,   0., 400, 600);
  fPads[1]= new TPad("pad1", "", 0, 0.30, 0.99, 0.99);   fPads[1]->Draw(); 
  fPads[2]= new TPad("pad2", "", 0, 0.00, 0.99, 0.28);   fPads[2]->Draw(); 

  TH1D *hratio  = new TH1D(*hd);   hratio->SetName("hratio"); hratio->Reset();
  hratio->Divide(hd, h1);

  TH1D *hratio2  = new TH1D(*hd);   hratio2->SetName("hratio2"); hratio2->Reset();
  if (h2) hratio2->Divide(hd, h2);

  // -- Set up histograms
  int fontSize(20); 
  
  hd->SetMinimum(0.0); 
  hd->SetMinimum(0.0000001); 
  hd->SetLabelFont(133, "Y"); 
  hd->SetLabelSize(fontSize*1.5, "Y"); 
  hd->SetLabelOffset(10000, "X");

  hratio->SetMinimum(0.6); 
  hratio->SetMaximum(1.4); 
  hratio->SetLabelFont(133, "X");  hratio->SetLabelFont(133, "Y"); 
  hratio->SetLabelSize(fontSize, "X");  hratio->SetLabelSize(fontSize, "Y");
  hratio->SetNdivisions(505, "Y");

  hratio->SetTitleOffset(0.9, "X");               hratio->SetTitleOffset(0.5, "Y");              
  hratio->SetTitleSize(0.15, "X");                hratio->SetTitleSize(0.15, "Y");
  hratio->SetXTitle(hd->GetXaxis()->GetTitle());  hratio->SetYTitle("Data / MC");

  hd->SetXTitle("");
  hd->SetYTitle("");
  
  setFilledHist(hd, kBlack, 8, 1000);
  setFilledHist(h1, kBlack, kBlue, 3004);
  if (h2) setFilledHist(h2, kBlack, kRed, 3005);

  c6->cd(1);
  fPads[1]->cd(); 
  shrinkPad(0.0, 0.15, 0.1, 0.05); 
  hd->Draw();
  h1->Draw("samehist");
  if (h2) h2->Draw("samehist");

  if (fString1.Contains("no legend", TString::kIgnoreCase)) {
    // -- do what?
  } else {
    leg = new TLegend(0.62,0.70,0.89,0.89);
    leg->SetHeader(fString1);
    leg->SetFillStyle(0); leg->SetBorderSize(0.); leg->SetTextSize(0.05);  leg->SetFillColor(0); 
    legge = leg->AddEntry(hd, fString2, "p"); 
    legge = leg->AddEntry(h1, fString3, "f"); legge->SetTextSize(0.05);
    if (h2) {
      legge = leg->AddEntry(h2, fString4, "f"); legge->SetTextSize(0.05);
    }
    leg->Draw();
  }

  fPads[2]->cd();
  shrinkPad(0.3, 0.15, 0.1, 0.001); 
  gPad->SetGridx(1);  gPad->SetGridy(1);
  setHist(hratio, kBlue, 24, 0.8, 1); 
  if (h2) setHist(hratio2, kRed, 20, 0.4, 1); 
  hratio->Draw("e");
  if (h2) hratio2->Draw("esame");


  double Chi2(0.), dof(0.);
  fNrProbA = chi2TestErr(hd, h1, Chi2, dof, -1);
  fNrChi2A = Chi2; 
  fNrDofA  = dof;
  fKstA = hd->KolmogorovTest(h1, "NOUD");

  if (h2) {
    fNrProbB = chi2TestErr(hd, h2, Chi2, dof, -1);
    fNrChi2B = Chi2; 
    fNrDofB  = dof;
    fKstB = hd->KolmogorovTest(h2, "NOUD");
  }

  //   hratio->Fit("pol0", "0");
  //   fCfit   = hratio->GetFunction("pol0")->GetParameter(0);
  //   fCfitE  = hratio->GetFunction("pol0")->GetParError(0);
  //   fCChi2  = hratio->GetFunction("pol0")->GetChisquare();
  //   fCDof   = hratio->GetFunction("pol0")->GetNDF();
  
  //   pl->SetLineColor(kBlue);
  //   pl->DrawLine(hratio->GetBinLowEdge(1), fCfit, hratio->GetBinLowEdge(hratio->GetNbinsX()+1), fCfit);
  
  //   hratio->Fit("pol1", "0");
  //   fOfit  = hratio->GetFunction("pol1")->GetParameter(0);
  //   fOfitE = hratio->GetFunction("pol1")->GetParError(0);
  //   fSfit  = hratio->GetFunction("pol1")->GetParameter(1);
  //   fSfitE = hratio->GetFunction("pol1")->GetParError(1);
  //   fChi2  = hratio->GetFunction("pol1")->GetChisquare();
  //   fDof   = hratio->GetFunction("pol1")->GetNDF();



  tl->SetNDC(kTRUE);
  tl->SetTextColor(kBlue);
  tl->SetTextSize(0.1);
//   tl->DrawLatex(0.16, 0.92, Form("offset = %4.3f#pm%4.3f, slope = %4.3f#pm%4.3f", fOfit, fOfitE, fSfit, fSfitE)); 
//   tl->DrawLatex(0.16, 0.825, Form("const  = %4.3f#pm%4.3f (#chi^{2}/dof = %4.2f/%3d)", fCfit, fCfitE, fCChi2, int(fCDof))); 
  tl->DrawLatex(0.16, 0.92,  Form("KS=%4.3f P=%4.3f (#chi^{2}/dof = %4.2f/%2.0f)", fKstA, fNrProbA, fNrChi2A, fNrDofA)); 
  if (h2) {
    tl->SetTextColor(kRed);
    tl->DrawLatex(0.16, 0.825, Form("KS=%4.3f P=%4.3f (#chi^{2}/dof = %4.2f/%2.0f)", fKstB, fNrProbB, fNrChi2B, fNrDofB)); 
  }

  tl->SetTextColor(1);

  //  (*hr) = hratio; 

}


// ----------------------------------------------------------------------
void b2uQA::fillStrings(const char *s1, const char *s2, const char *s3, const char *s4, const char *s5) {
  fString1 = TString(s1); 
  fString2 = TString(s2); 
  fString3 = TString(s3); 
  fString4 = TString(s4); 
  fString5 = TString(s5); 
}


// ----------------------------------------------------------------------
void b2uQA::setFilledHist(TH1D *h, Int_t color, Int_t fillcolor, Int_t fillstyle, Int_t width) {
  // Note: 3004, 3005 are crosshatches
  // ----- 1000       is solid
  //       kYellow    comes out gray on bw printers
  h->SetLineColor(color);     h->SetLineWidth(width);   
  h->SetFillStyle(fillstyle); h->SetFillColor(fillcolor);
}

// ----------------------------------------------------------------------
void b2uQA::shrinkPad(double b, double l, double r, double t) {
  gPad->SetBottomMargin(b); 
  gPad->SetLeftMargin(l);
  gPad->SetRightMargin(r);
  gPad->SetTopMargin(t);
}

// ----------------------------------------------------------------------
void b2uQA::setTitles(TH1D *h, const char *sx, const char *sy, float size, 
	       float xoff, float yoff, float lsize, int font) {
  if (h == 0) {
    cout << " Histogram not defined" << endl;
  } else {
    h->SetXTitle(sx);                  h->SetYTitle(sy); 
    h->SetTitleOffset(xoff, "x");      h->SetTitleOffset(yoff, "y");
    h->SetTitleSize(size, "x");        h->SetTitleSize(size, "y");
    h->SetLabelSize(lsize, "x");       h->SetLabelSize(lsize, "y");
    h->SetLabelFont(font, "x");        h->SetLabelFont(font, "y");
    h->GetXaxis()->SetTitleFont(font); h->GetYaxis()->SetTitleFont(font);
    h->SetNdivisions(508, "X");
  }
}

// ----------------------------------------------------------------------
void b2uQA::setHist(TH1D *h, Int_t color, Int_t symbol, Double_t size, Double_t width) {
  h->SetLineColor(color);   h->SetLineWidth(width);
  h->SetMarkerColor(color); h->SetMarkerStyle(symbol);  h->SetMarkerSize(size); 
  h->SetStats(kFALSE); 
  h->SetFillStyle(0); h->SetFillColor(color);
}


// ----------------------------------------------------------------------
double b2uQA::chi2TestErr(TH1D *h1, TH1D *h2, double& Chi2, double& ndof, int constrain) {
  int nbins = h1->GetNbinsX();
  if (nbins != h2->GetNbinsX()) {
    cout << "chi2Test: Number of bins not the same" << endl;
    return -99.;
  }
  double df = nbins - 1 - constrain; 
  double chsq(0.), a1(0.), a2(0.), e1(0.), e2(0.), A1(h1->GetSumOfWeights()), A2(h2->GetSumOfWeights());
  for (int i = 1; i <= nbins; ++i) {
    a1 = h1->GetBinContent(i);
    e1 = h1->GetBinError(i) * h1->GetBinError(i);
    a2 = h2->GetBinContent(i);
    e2 = h2->GetBinError(i) *  h2->GetBinError(i);
    if ((TMath::Abs(a1) < 1.e-8) && (TMath::Abs(a2) < 1.e-8)) {
      cout << "chi2Test: Skipping bin " << i << " with < 1.e-8 entries" << endl;
      df -= 1.;
    } else if ((a1 < 0.) || (a2 < 0.)) {
      cout << "chi2Test: Skipping bin " << i << " with negative entries" << endl;
      df -= 1.;
    } else if ((a1/A1 < 0.01) || (a2/A2 < 0.01)) {
      cout << "chi2Test: Skipping bin " << i << " with less than 1% contribution" << endl;
      df -= 1.;
    } else {
      chsq += ((a1 - a2) * (a1 - a2)) / (e1 + e2);
    }
  }
  double gamma = 1. - TMath::Gamma(0.5*df, 0.5*chsq);
  Chi2 = chsq;
  ndof = df;
  return gamma;
}


// ----------------------------------------------------------------------
//double b2uQA::getBrecoWeight(double theintpur) {
//  int thebin = theintpur/(1./20);
//  if (theintpur>1.)thebin = 19; 
//  if (theintpur<0.)thebin = 0; 
//  return BrecoWeight[thebin];
//}

//Routines needed for reweighting calculations
// ----------------------------------------------------------------------
double b2uQA::getBsysweight(int decType,int thevub) {
  double theweight;
  theweight = 1.;  
  theweight = Bsem->weight(decType); 
  //if(DOFFWEIGHT && Gvxbtyp==2&&TMath::Abs(ctvgen)<2.){// B0 -> D*lnu FF
  if(Gvxbtyp==2&&TMath::Abs(ctvgen)<2.){// B0 -> D*lnu FF
    double def=dstlnuFF(1.18,0.72,0.92);    
    def=(def>0)?dstlnuFF(1.34,0.91,0.81)/def:1;

    if(def>100){
      //      cout<<"too high dstarlnu weight"<<def;
      //      cout <<" kin "<<ctvgen<<" "<<ctlgen<<" " <<chigen<<" "<<q2Gen<<endl;
      def=10.;
    }
    
    //    cout<<" wei "<<def<<endl;
    theweight*=def;
    
  }
  if(thevub) theweight = 1.;
  return theweight;
}

// ----------------------------------------------------------------------
double b2uQA::getDsysweight(int decDpi,int decDk,int decDks,int decDpiz,int decDlep,int decImode,int thevub) {
  double theweight;
  theweight = 1.;  
  int bla(1);
  //if(DODDECWEIGHT){
  theweight = Dvar->weight(decDpi,decDk,decDks,decDpiz,decDlep,decImode,bla);    
  //}
  if(thevub) theweight = 1.;
  return theweight;
}

double b2uQA::dstlnuFF(double r1,double r2,double rho2){
  // now calculate the 4D pdf
  float mb(5.279),mdst(2.01);
  float hp,h0,hm,a1,v,a2,ha1;
  float cchi=cos(chigen);

  a2=(mb+mdst)/2/sqrt(mb*mdst);
  a1=a2*(1-q2Gen/pow((mb+mdst),2));
  v=a2*r1;
  a2=a2*r2;
  // h_a1(w) has ben factorized in front of everything
  ha1=1-rho2*((mb*mb+mdst*mdst-q2Gen)/2/mb/mdst-1);

  hp=(mb+mdst)*a1-2*mb*pxhadgen*v/(mb+mdst);
  hm=(mb+mdst)*a1+2*mb*pxhadgen*v/(mb+mdst);
  h0=(mb*mb-mdst*mdst-q2Gen)*(mb+mdst)*a1-pow((2*mb*pxhadgen),2)*a2/(mb+mdst);
  h0=h0/2/mdst/sqrt(q2Gen);

  return ha1*ha1*pxhadgen*q2Gen*(
				 (pow((1+ctlgen)*hp,2)+pow((1-ctlgen)*hm,2))*(1-ctvgen*ctvgen)+
				 (1-ctlgen*ctlgen)*pow((2*h0*ctvgen),2)-
				 2*(1-ctlgen*ctlgen)*(1-ctvgen*ctvgen)*(2*cchi*cchi-1)*hp*hm+
				 4*sqrt((1-ctlgen*ctlgen)*(1-ctvgen*ctvgen))*ctvgen*cchi*h0*
				 (-ctlgen*(hp+hm)+(hm-hp))
				 );

}
