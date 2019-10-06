//////////////////////////////////////////////////////////
//   This class has been automatically generated 
//     (Mon May 13 08:12:21 2002 by ROOT version3.01/06)
//   from TTree events/events
//   found on file: /nfs/farm/babar/AWG7/ISL/tmp/rootfitfiles/sx-allgeneric.root
//////////////////////////////////////////////////////////


#ifndef THECOMPARISON_H
#define THECOMPARISON_H

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TPad.h"
#include "TRotation.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "../../RecoilAnalysis/recoilAnalysis.hh" 
#include "../../RecoilAnalysis/recoilDSys.hh" 
//#include "../../RecoilAnalysis/mesData.hh"  //included by recoilDSys.hh

class thecomparison {
   public :
   TTree          *fChain;   //pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //current Tree number in a TChain
   recoilDSys *Dvar;
   recoilDSys *Bsem;

   Int_t           run;
   Int_t           upper;
   Int_t           lower;
   Int_t           nchg;
   Int_t           nneu;
   Int_t           nB;
   Int_t           IdB;
   Int_t           brecoflav;
   Int_t           brecocharge;
   Int_t           xcharge;
   Int_t           mode;
   Double_t        mes;
   Double_t        de;
   Double_t        pur;
   Double_t        intpur;
   Int_t           nle;
   Int_t           nel;
   Int_t           nmu;
   Int_t           nkp;
   Int_t           nks;
   Int_t           nlept500;
   Double_t        plab;
   Double_t        tlab;
   Double_t        flab;
   Double_t        elab;
   Double_t        pcms;
   Double_t        ecms;
   Double_t        tcms;
   Double_t        fcms;
   Int_t           lcharge;
   Int_t           isele;
   Double_t        wdeltam;
   Double_t        Eneualt;

   // pi l nu

    Int_t           indexbestPi; 
    Int_t           chbestPi; 
    Double_t        barembestPi; 
    Double_t        mm2bestPi; 
    Double_t        q2bestPi; 
    Int_t           nrecoPi; 

    // Pi0 l nu

    Int_t           indexbestPi0; 
    Int_t           chbestPi0; 
    Double_t        barembestPi0; 
    Double_t        mm2bestPi0; 
    Double_t        q2bestPi0; 
    Int_t           nrecoPi0; 
    Int_t           ndauPi0[100];   //[nrecoPi0] 
    Float_t         Estar1dauPi0[100];   //[nrecoPi0] 
    Float_t         Estar2dauPi0[100];   //[nrecoPi0] 
 
    //eta l nu

   Int_t           indexbestEta;
   Double_t        barembestEta;
   Double_t        mm2bestEta;
   Int_t           nrecoEta;
   Int_t           ndauEta[100];   //[nrecoEta]
   Int_t           modeEta[100];   //[nrecoEta]
   Float_t         Estar1dauEta[100];   //[nrecoEta]
   Float_t         Estar2dauEta[100];   //[nrecoEta]
   Float_t         Estar3dauEta[100];   //[nrecoEta]
   Float_t         Elab1dauEta[100];   //[nrecoEta]
   Float_t         Elab2dauEta[100];   //[nrecoEta]
   Float_t         Elab3dauEta[100];   //[nrecoEta]

   //etap l nu

   Int_t           indexbestEtap;
   Double_t        barembestEtap;
   Double_t        mm2bestEtap;
   Int_t           nrecoEtap;
   Int_t           ndauEtap[100];   //[nrecoEtap]
   Int_t           modeEtap[100];   //[nrecoEtap]
   Float_t         EtamassdauEtap[100];   //[nrecoEtap]
   Float_t         Rho0massdauEtap[100];   //[nrecoEtap]
   Float_t         GammamomdauEtap[100];   //[nrecoEtap]
   Float_t         Estar1dauEtap[100];   //[nrecoEtap]
   Float_t         Estar2dauEtap[100];   //[nrecoEtap]
   Float_t         Estar3dauEtap[100];   //[nrecoEtap]
   Float_t         Elab1dauEtap[100];   //[nrecoEtap]
   Float_t         Elab2dauEtap[100];   //[nrecoEtap]
   Float_t         Elab3dauEtap[100];   //[nrecoEtap]
   
//MC truth

/*     Int_t           isassocB;  */
/*     Int_t           ass_deltapB; */
//     Int_t           vub;  
//    Int_t           vcb;  
//     Int_t           other; 
/*     Int_t           nvubexcl;  */
/*     Double_t        mxhadgen; */
/*     Double_t        pcmsgen;  */
/*     Double_t        ecmsgen;  */
/*     Double_t        pxhadgen;  */
/*     Double_t        exhadgen;  */
/*     Double_t        q2Gen;  */
//     Int_t           Gvxbtyp; 
/*     Int_t           GSem;  */
/*     Int_t           GfDpi;  */
/*     Int_t           GfDpiz;  */
/*     Int_t           GfDk;  */
/*     Int_t           GfDks;  */
/*     Int_t           GfDkl;  */
/*     Int_t           GfDlep;  */
/*     Int_t           GfDgam;  */
/*     Int_t           GfDnu;  */
/*     Int_t           GfD0Ds;  */
/*     Int_t           GfDDs;  */
/*     Int_t           GfDkspiopio;  */

//List of branches

   TBranch        *b_run;
   TBranch        *b_lower;
   TBranch        *b_upper;
   TBranch        *b_nchg;
   TBranch        *b_nneu;
   TBranch        *b_nB;
   TBranch        *b_IdB;
   TBranch        *b_brecoflav;
   TBranch        *b_brecocharge;
   TBranch        *b_xcharge;
   TBranch        *b_mode;
   TBranch        *b_mes;
   TBranch        *b_de;
   TBranch        *b_pur;
   TBranch        *b_intpur;
   TBranch        *b_nle;
   TBranch        *b_nel;
   TBranch        *b_nmu;
   TBranch        *b_nkp;
   TBranch        *b_nks;
   TBranch        *b_nlept500;
   TBranch        *b_plab;
   TBranch        *b_tlab;
   TBranch        *b_flab;
   TBranch        *b_elab;
   TBranch        *b_pcms;
   TBranch        *b_ecms;
   TBranch        *b_tcms;
   TBranch        *b_fcms;
   TBranch        *b_lcharge;
   TBranch        *b_isele;
   TBranch        *b_wdeltam;
   TBranch        *b_Eneualt;

   //pi l nu 

   TBranch        *b_indexbestPi;    
   TBranch        *b_chbestPi; 
   TBranch        *b_barembestPi;   
   TBranch        *b_mm2bestPi; 
   TBranch        *b_q2bestPi;   
   TBranch        *b_nrecoPi; 

   // pi0 l nu
  
   TBranch        *b_indexbestPi0; 
   TBranch        *b_chbestPi0;   
   TBranch        *b_barembestPi0;   
   TBranch        *b_mm2bestPi0;   
   TBranch        *b_q2bestPi0;  
   TBranch        *b_nrecoPi0;   
   TBranch        *b_ndauPi0;  
   TBranch        *b_Estar1dauPi0;  
   TBranch        *b_Estar2dauPi0;    

   // eta l nu 

   TBranch        *b_indexbestEta;
   TBranch        *b_barembestEta;
   TBranch        *b_mm2bestEta;
   TBranch        *b_nrecoEta;
   TBranch        *b_ndauEta; 
   TBranch        *b_modeEta; 
   TBranch        *b_Estar1dauEta;  
   TBranch        *b_Estar2dauEta;    
   TBranch        *b_Estar3dauEta; 
   TBranch        *b_Elab1dauEta; 
   TBranch        *b_Elab2dauEta; 
   TBranch        *b_Elab3dauEta; 

   //etap l nu

   TBranch        *b_indexbestEtap;
   TBranch        *b_barembestEtap;
   TBranch        *b_mm2bestEtap;
   TBranch        *b_nrecoEtap;
   TBranch        *b_ndauEtap; 
   TBranch        *b_modeEtap; 
   TBranch        *b_EtamassdauEtap; 
   TBranch        *b_Rho0massdauEtap;
   TBranch        *b_GammamomdauEtap;
   TBranch        *b_Estar1dauEtap; 
   TBranch        *b_Estar2dauEtap; 
   TBranch        *b_Estar3dauEtap; 
   TBranch        *b_Elab1dauEtap;  
   TBranch        *b_Elab2dauEtap; 
   TBranch        *b_Elab3dauEtap; 

   //MC truth

/*    TBranch        *b_isassocB;  */
/*    TBranch        *b_ass_deltapB;  */
//    TBranch        *b_vub;
//    TBranch        *b_vcb;  
//    TBranch        *b_other; 
/*    TBranch        *b_nvubexcl;  */
/*    TBranch        *b_mxhadgen; */
/*    TBranch        *b_pcmsgen;  */
/*    TBranch        *b_ecmsgen;  */
/*    TBranch        *b_pxhadgen;  */
/*    TBranch        *b_exhadgen; */
/*    TBranch        *b_q2Gen; */
//    TBranch        *b_Gvxbtyp; 
/*    TBranch        *b_GSem;  */
/*    TBranch        *b_GfDpi; */
/*    TBranch        *b_GfDpiz;  */
/*    TBranch        *b_GfDk; */
/*    TBranch        *b_GfDks;  */
/*    TBranch        *b_GfDkl;  */
/*    TBranch        *b_GfDlep;  */
/*    TBranch        *b_GfDgam; */
/*    TBranch        *b_GfDnu;  */
/*    TBranch        *b_GfD0Ds;  */
/*    TBranch        *b_GfDDs;  */
/*    TBranch        *b_GfDkspiopio;  */

  // ----------------------------------------------------------------------
  // -- CUTS
  Double_t TOTALSTAT, TOTALSTATMODEL, BRRATIOGENVALUE,  BRRATIOVALUETAIL, MXCUT, Q2CUT, PSTARFACT, LEPTONPCUT, ELECTRONPCUT, MUONPCUT, PRMM2CUT;
  Double_t MNUSQLOW, MNUSQHIGH,  CHLOW, CHHIGH, DEPL, SSBAR, BTYPE, LEPTTYPE, MAXINTPUR, MININTPUR, RUN, TOYMC;
  Double_t DELTAMB, DELTAA, FERMIAPP;
  Int_t USECB, FIXMEANVALUE, FIXSIGMA, FIXARGUS1, FIXARGUS2;
  Int_t FITCATEGORY, FIXCB1, FIXCB2, FITTOTSHAPE, MIXCORR, FITMC, FITOPT, MULTIFIT, BLINDING;
  Int_t ISSMEARALL, ISSMEARBKG, CUTNNPI0;
  Double_t  SMEARALLMEANVALUE, SMEARALLSIGMA, SMEARBKGMEANVALUE, SMEARBKGSIGMA;
  Double_t  RANDOMSEED, BLINDSIZE;
  Int_t DOBRECOWEIGHT, DOBDECWEIGHT, DODDECWEIGHT, DOTRKWEIGHT, DONEUWEIGHT, DOFERMI, DOTHEO;
  Double_t MXCUTLOWEXCL, MXCUTHIGHEXCL, DELTAM;
  Double_t MXCUTLOWEXCL1, MXCUTHIGHEXCL1, MXCUTLOWEXCL2, MXCUTHIGHEXCL2, MXCUTLOWEXCL3, MXCUTHIGHEXCL3, MXCUTLOWEXCL4, MXCUTHIGHEXCL4;
  Double_t NCHGLOWEXCL, NCHGHIGHEXCL;
  Double_t NCOMBLOWEXCL, NCOMBHIGHEXCL;      
  Double_t NPI0LOWEXCL, NPI0HIGHEXCL;
  Double_t MOM1MIN, MOM2MIN;
  Double_t MAXENEU,JPSIWIN;
  Double_t DALITZCUT;
  Double_t DAUETAMASS, DAURHOMASS, DAUGAMMAMOM;
  Double_t DAUETAMASS2, DAUETAMASS3, DAUETAMASS4;
  Int_t DORIGHTNCHG;
  Double_t MNUSQPI0LOW, MNUSQPI0HIGH, MNUSQETALOW, MNUSQETAHIGH, MNUSQRHOLOW, MNUSQRHOHIGH, MNUSQRHO0LOW, MNUSQRHO0HIGH, MNUSQOMEGALOW, MNUSQOMEGAHIGH; 

   double SHIFTNEUT;
   double SIGMANEUT;   
   int multipl;
   double intdataen;
   double intdatadepl;
   double intMCen;
   double intMCdepl;     

   double themin;
   double themax;
   int thebins;
   TString thevar;
   TPad *fPads[50];
   TCanvas *c0, *c1, *c2, *c3, *c4, *c5, *c6, *c7;
   thecomparison();
   thecomparison(char *var, double mi, double ma, int b);
   ~thecomparison();
   Int_t  Cut(Int_t entry);
   Int_t  GetEntry(Int_t entry);
   Int_t  LoadTree(Int_t entry);

   double    getBsysweight(int,int);
   double    getDsysweight(int,int,int,int,int,int,int);
   double    FermiWeight(double,double,double);
   double    fermi(double,double,double);

   void   Init(TTree *tree);
   void   Loop(int nevents, int cat, double shift, double smear);
   // seed = 2 dstar
   // seed = 3 dstar0
   // seed = 4 dc
   // seed = 5 d0
   // sys = 0: nothing; 1: B&D defaults; 2 fermi reweighting (to values set by hand still...)
   void   Bookhist();
   Bool_t Notify();
   void   Show(Int_t entry = -1);
   int    hist(double mx);
   void   Fitmes(int cat, int cut);
   void sighisto(double& signal, double& signalErr,TH1D *histo, double &resmean, double &ressigma, double &resalpha, double &resn, int fixpar,double mean, double sigma, double alpha, double n, double argus);
   void   overlap(int cut, int norm, TString dir);
   void   effplots(TString dir);
   TChain* getchain(char *thechain, char *treename);
   double chisq(TH1 *h1, TH1 *h2);
   Double_t  smeargauss(double,double,double);
   void readCuts(TString filename);
   void dumpCuts();    
   
   ClassDef(thecomparison,1)

};

#endif
