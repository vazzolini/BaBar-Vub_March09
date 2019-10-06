//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jan  3 17:41:28 2005 by ROOT version 4.00/08
// from TTree events/events
// found on file: csx-genbch-new-2002b.root
//////////////////////////////////////////////////////////

#ifndef b2uQA_h
#define b2uQA_h

#include <iostream.h>

#include <TROOT.h>
#include <TColor.h>
#include <TChain.h>
#include <TH1.h>
#include <TFile.h>
#include <TPad.h>

#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"

#include "sHist.hh"
#include "../RecoilAnalysis/recoilDSys.hh"

class b2uQA {
public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
  // Declaration of leave types
  Int_t           run;
  Double_t        mes;
  Double_t        intpur;
  Int_t           brecoflav;
  Int_t           brecocharge;
  Int_t           xcharge;
  Double_t        mxhad;
  Double_t        mxtrk;
  Double_t        mxnut;
  Double_t        mxhadfit;
  Double_t        chi2;
  Double_t        probchi2;
  Int_t           lcharge;
  Double_t        pcms;
  Double_t        tcms;
  Double_t        tlab;
  Double_t        plab;
  Double_t        mm2;
  Double_t        q2fit;
  Double_t        pmiss;
  Double_t        tmiss;
  Double_t        emiss;
  Int_t           nle;
  Int_t           nel;
  Int_t           nmu;
  Int_t           nchg;
  Int_t           nneu;
  Int_t           nks;
  Int_t           nkz;
  Int_t           nkl;
  Int_t           nkp;
  Double_t        deltam;
  Double_t        mm1pr;
  Double_t        mm2pr;
  Double_t        mm3pr;
  Double_t        oa1;
  Double_t        oa2;
  Double_t        oa3;
  Double_t        pcmstrklo;
  Int_t           vub;
  Int_t           vcb;
  Int_t           other;
  Int_t           qb;
  Double_t        kplus;
  Double_t        ecmsgen;
  Double_t        tcmsgen;
  Double_t        fcmsgen;
  Double_t        mxhadgen;
  Double_t        pxhadgen;
  Double_t        q2Gen;
  Double_t        ctvgen;
  Double_t        ctlgen;
  Double_t        chigen;
  Double_t        vpgen;
  Double_t        vtgen;
  Double_t        vfgen;
  Double_t        vpcmsgen;
  Double_t        vtcmsgen;
  Double_t        vfcmsgen;
  Int_t           Gvxbtyp;
  Int_t           GfDpi;
  Int_t           GfDk;
  Int_t           GfDks;
  Int_t           GfDpiz;
  Int_t           GfDlep;

  // List of branches
  TBranch        *b_run;   //!
  TBranch        *b_mes;   //!
  TBranch        *b_intpur;   //!
  TBranch        *b_brecoflav;   //!
  TBranch        *b_brecocharge;   //!
  TBranch        *b_xcharge;   //!
  TBranch        *b_mxhad;   //!
  TBranch        *b_mxtrk;   //!
  TBranch        *b_mxnut;   //!
  TBranch        *b_mxhadfit;   //!
  TBranch        *b_chi2;   //!
  TBranch        *b_probchi2;   //!
  TBranch        *b_lcharge;   //!
  TBranch        *b_pcms;   //!
  TBranch        *b_tcms;   //!
  TBranch        *b_tlab;   //!
  TBranch        *b_plab;   //!
  TBranch        *b_mm2;   //!
  TBranch        *b_q2fit;   //!
  TBranch        *b_pmiss;   //!
  TBranch        *b_tmiss;   //!
  TBranch        *b_emiss;   //!
  TBranch        *b_nle;   //!
  TBranch        *b_nel;   //!
  TBranch        *b_nmu;   //!
  TBranch        *b_nchg;   //!
  TBranch        *b_nneu;   //!
  TBranch        *b_nks;   //!
  TBranch        *b_nkz;   //!
  TBranch        *b_nkl;   //!
  TBranch        *b_nkp;   //!
  TBranch        *b_wdeltam;   //!
  TBranch        *b_mm1pr;   //!
  TBranch        *b_mm2pr;   //!
  TBranch        *b_mm3pr;   //!
  TBranch        *b_oa1;   //!
  TBranch        *b_oa2;   //!
  TBranch        *b_oa3;   //!
  TBranch        *b_pcmstrklo;   //!
  TBranch        *b_vub;   //!
  TBranch        *b_vcb;   //!
  TBranch        *b_other;   //!
  TBranch        *b_qb;   //!
  TBranch        *b_fkplus;   //!
  TBranch        *b_ecmsgen;   //!
  TBranch        *b_tcmsgen;   //!
  TBranch        *b_fcmsgen;   //!
  TBranch        *b_mxhadgen;   //!
  TBranch        *b_pxhadgen;   //!
  TBranch        *b_q2Gen;   //!
  TBranch        *b_ctvgen;   //!
  TBranch        *b_ctlgen;   //!
  TBranch        *b_chigen;   //!
  TBranch        *b_vpgen;   //!
  TBranch        *b_vtgen;   //!
  TBranch        *b_vfgen;   //!
  TBranch        *b_vpcmsgen;   //!
  TBranch        *b_vtcmsgen;   //!
  TBranch        *b_vfcmsgen;   //!
  TBranch        *b_Gvxbtyp;   //!
  TBranch        *b_GfDpi;   //!
  TBranch        *b_GfDk;   //!
  TBranch        *b_GfDks;   //!
  TBranch        *b_GfDpiz;   //!
  TBranch        *b_GfDlep;   //!

  b2uQA(TTree *tree=0);
  b2uQA(const char *filename, int Sys, TString cutfile, int reweight);
  ~b2uQA();
  void   bookHist(int mode = 0);
  Int_t  GetEntry(Int_t entry);
  Int_t  LoadTree(Int_t entry);
  Bool_t Notify();

  void   readCuts(TString cutfile);
  void   Init(TTree *tree);
  void   loop(int mode);   
  // -- run some mode
  void   makeAll(int mode);
  // -- For running of the saved histrogams
  void   readHist(int mode = 0);
  // -- Creates the overlay plots
  void   show(TH1D *hd, TH1D *h1, TH1D *h2);

  // -- Utilities
  void fillStrings(const char *s1 = "", const char *s2 = "", const char *s3 = "", const char *s4 = "", const char *s5 = "");
  void   shrinkPad(double b = 0.1, double l = 0.1, double r = 0.1, double t = 0.1);
  void   setFilledHist(TH1D *h, int lcol = kBlack, int fcol = kYellow, int fstyle = 1000, int width = 1);
  void   setTitles(TH1D *h, const char *sx, const char *sy, float size = 0.05, 
		   float xoff = 1.1, float yoff = 1.1, float lsize = 0.05, int font = 132);
  void   setHist(TH1D *h, int color = kBlack, int symbol = 20, double size = 1., double width = 2.);

  // -- Chi^2 test for two histograms, where the errors are taken from the histogram
  //    If they are constrained to the same area, constrain = 0
  //    If they are NOT constrained, set constrain = -1 
  //    This function is equivalent to chi2Test() if the errors are sqrt(n_bin)
  double chi2TestErr(TH1D*, TH1D*, double& chi2, double& ndof, int constrain = 0);

  //double getBrecoWeight(double theintpur);
  double getBsysweight(int decType,int thevub);
  double getDsysweight(int decDpi,int decDk,int decDks,int decDpiz,int decDlep,int decImode,int thevub);
  double dstlnuFF(double r1,double r2,double rho2); 
  
  // ----------------------------------------------------------------------
  // -- QA stuff
  double dqtot, dnchg, dnneu;
  double depmiss, dctmiss;


  int fMC, fReW;
  recoilDSys *Dvar;
  recoilDSys *Bsem;
  int dImode;

  double QTOTCUT, MM2CUTHI, MM2CUTLO, EMPCUTHI, EMPCUTLO;
  double CTMCUTLO, CTMCUTHI, PMISSCUT, PCMSTRKLO;
  double PRMM1CUT, PRMM2CUT, PRMM3CUT, PCMSCUT;


  double fKstA, fNrChi2A, fNrProbA, fNrDofA;
  double fKstB, fNrChi2B, fNrProbB, fNrDofB;
  double fCfit, fCfitE, fCChi2, fCDof;
  double fOfit, fOfitE, fSfit, fSfitE, fChi2, fDof;

  TString fFilename; // filename of output histograms
  TFile *fHistFile;  // output histograms
  sHist *fH800, *fH900, *fH1000, *fH1100
    , *fH2000, *fH2100, *fH2200, *fH2300 
    , *fH3000, *fH3100, *fH3200, *fH3300, *fH3400, *fH3500
    , *fH4000, *fH4100, *fH4200, *fH4300 
    ;

  TH1D *fHmes;


  // -- Display utilities
  TPad *fPads[50];
  TLatex *tl; 
  TLine *pl; 
  TLegend *leg;
  TLegendEntry *legge;

  TString fString1, fString2, fString3, fString4, fString5;


};


#endif
