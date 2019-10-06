
#ifndef theanal_h
#define theanal_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "../RecoilAnalysis/recoilAnalysis.hh" 
#include "../VubAnalysis/recoilDSys.hh" 

class theanal {
   public :
   TFile *fHistFile;
   TTree          *fChain;   //pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //current Tree number in a TChain
   recoilDSys *Dvar;
   recoilDSys *Bsem;
//Declaration of leaves types
   Double_t        mes;
   Double_t        de;
   Double_t        pur;
   Double_t        intpur;
   Int_t           mode;
   Int_t           brecoflav;
   Int_t           brecocharge;
   Double_t        mxhadgen;
   Double_t        pcmsgen;
   Double_t        tcmsgen;
   Double_t        fcmsgen;
   Double_t        ecmsgen;
   Double_t        pxhadgen;
   Double_t        txhadgen;
   Double_t        fxhadgen;
   Double_t        exhadgen;
   Int_t           vub;
   Int_t           vcb;
   Int_t           vxbtyp;
   Int_t           other;
   Int_t           xcharge;
   Double_t        mxhad;
   Double_t        mxhadfit;
   Double_t        mxhadchg;
   Double_t        q2;
   Double_t        q2fit;
   Double_t        q2gen;
   Int_t           lcharge;
   Double_t        plab;
   Double_t        tlab;
   Double_t        flab;
   Double_t        pcms;
   Double_t        tcms;
   Double_t        fcms;
   Double_t        ecms;
   Int_t           nle;
   Int_t           nlept500;
   Int_t           nel;
   Int_t           nmu;
   Int_t           isele;
   Int_t           nchg;
   Int_t           nneu;
   Int_t           nkp;
   Int_t           nks;
   Double_t        mm2;
   Double_t        wdeltam;
   Int_t           Gvxbtyp;
   Double_t        Eneualt;

  //pi stuff    
   Int_t           nrecoPi;
   Double_t        mm2bestPi;
   Float_t         baremassjpsiPi[100];
   Int_t           indexbestPi;
   Int_t           chbestPi;
  
    //pi0 stuff
   Int_t           indexbestPi0;
   Int_t           nrecoPi0;
   Double_t        barembestPi0;
   Double_t        mm2bestPi0;
   Double_t        mm2gamma;
   Double_t        truemom1phpi0;
   Double_t        truemom2phpi0;
   Double_t        truemomlab1phpi0;
   Double_t        truemomlab2phpi0;
   Double_t        trueth1phpi0;
   Double_t        trueth2phpi0;
   Float_t         Estar1dauPi0[100];
   Float_t         Estar2dauPi0[100];

    //rho stuff;
   Double_t        truemrho;
   Int_t           nrecorho;
   Double_t        mrho;
   Double_t        mm2rho;
   Double_t        truemompirho;
   Double_t        truemompi0rho;
   Double_t        mompirho;
   Double_t        mompi0rho;
      
     //rho0 stuff;
   Double_t        truemrho0;
   Int_t           nrecorho0;
   Double_t        mrho0;
   Double_t        mm2rho0;
   Double_t        truemom1pirho0;
   Double_t        truemom2pirho0;
   Double_t        mom1pirho0;
   Double_t        mom2pirho0;

    //omega stuff
   Double_t        truemomega;
   Int_t           nrecoomega;
   Double_t        momega;
   Double_t        mm2omega;
   Double_t        truemom1piome;
   Double_t        truemom2piome;
   Double_t        truemompi0ome;
   Double_t        truedalitzpi1pi2ome;
   Double_t        truedalitzpi1pi0ome;
   Double_t        truecosthome;
   Double_t        mom1piome;
   Double_t        mom2piome;
   Double_t        mompi0ome;
   Double_t        dalitzpi1pi2ome;
   Double_t        dalitzpi1pi0ome;
   Double_t        costhome;

    //eta stuff
   Int_t           nrecoEta;
   Int_t           indexbestEta;
   Double_t        barembestEta;
   Double_t        mm2bestEta;
   Int_t           modeEta[100];
   Double_t        metagg;
   Double_t        mm2etagg;
   Double_t        metapppi0;
   Double_t        mm2etapppi0;
   Double_t        metapi0pi0pi0;
   Double_t        mm2etapi0pi0pi0;
  
   //eta' stuff
   Int_t           nrecoEtap;
   Int_t           indexbestEtap;
   Double_t        barembestEtap;
   Double_t        mm2bestEtap;
   Int_t           modeEtap[100];
   Float_t         GammamomdauEtap[100];
   Float_t         EtamassdauEtap[100];
   Float_t         Rho0massdauEtap[100];

   //a0 stuff
   Double_t        ma0;
   Double_t        mm2a0;  
   Double_t        a0massetadau;  
   Int_t           modea0;  
  
  
   //a0p stuff
   Double_t        ma0p;
   Double_t        mm2a0p;  
   Double_t        a0pmassetadau;  
   Int_t           modea0p;  
  


//List of branches
   TBranch        *b_mes;
   TBranch        *b_de;
   TBranch        *b_pur;
   TBranch        *b_intpur;
   TBranch        *b_mode;
   TBranch        *b_brecoflav;
   TBranch        *b_brecocharge;
   TBranch        *b_mxhadgen;
   TBranch        *b_pcmsgen;
   TBranch        *b_tcmsgen;
   TBranch        *b_fcmsgen;
   TBranch        *b_ecmsgen;
   TBranch        *b_vub;
   TBranch        *b_vcb;
   TBranch        *b_vxbtyp;
   TBranch        *b_other;
   TBranch        *b_xcharge;
   TBranch        *b_mxhad;
   TBranch        *b_mxhadfit;
   TBranch        *b_mxhadchg;
   TBranch        *b_q2;
   TBranch        *b_q2fit;
   TBranch        *b_q2gen;
   TBranch        *b_lcharge;
   TBranch        *b_plab;
   TBranch        *b_tlab;
   TBranch        *b_flab;
   TBranch        *b_pcms;
   TBranch        *b_tcms;
   TBranch        *b_fcms;
   TBranch        *b_ecms;
   TBranch        *b_nle;
   TBranch        *b_nlept500;
   TBranch        *b_nel;
   TBranch        *b_nmu;
   TBranch        *b_isele;
   TBranch        *b_nchg;
   TBranch        *b_nneu;
   TBranch        *b_nkp;
   TBranch        *b_nks;
   TBranch        *b_mm2;
   TBranch        *b_mm2nc;
   TBranch        *b_mm2fit;
   TBranch        *b_wdeltam;
   TBranch        *b_Gvxbtyp;
   TBranch        *b_eneualt;

   //pi stuff 
   TBranch        *b_nrecopi;
   TBranch        *b_mm2bestpi;
   TBranch        *b_baremassjpsipi;
   TBranch        *b_indexbestpi;
   TBranch        *b_chbestpi;

    //pi0 stuff
   TBranch        *b_indexbestpi0;
   TBranch        *b_nrecopi0;
   TBranch        *b_barembestpi0;
   TBranch        *b_mm2bestpi0;
   TBranch        *b_mm2gamma;
   TBranch        *b_truemom1phpi0;
   TBranch        *b_truemom2phpi0;
   TBranch        *b_truemomlab1phpi0;
   TBranch        *b_truemomlab2phpi0;
   TBranch        *b_trueth1phpi0;
   TBranch        *b_trueth2phpi0;
   TBranch        *b_Estar1daupi0;
   TBranch        *b_Estar2daupi0;

    //rho stuff;
   TBranch        *b_truemrho;
   TBranch        *b_nrecorho;
   TBranch        *b_mrho;
   TBranch        *b_mm2rho;
   TBranch        *b_truemompirho;
   TBranch        *b_truemompi0rho;
   TBranch        *b_mompirho;
   TBranch        *b_mompi0rho;
      
     //rho0 stuff;
   TBranch        *b_truemrho0;
   TBranch        *b_nrecorho0;
   TBranch        *b_mrho0;
   TBranch        *b_mm2rho0;
   TBranch        *b_truemom1pirho0;
   TBranch        *b_truemom2pirho0;
   TBranch        *b_mom1pirho0;
   TBranch        *b_mom2pirho0;

    //omega stuff
   TBranch        *b_truemomega;
   TBranch        *b_nrecoomega;
   TBranch        *b_momega;
   TBranch        *b_mm2omega;
   TBranch        *b_truemom1piome;
   TBranch        *b_truemom2piome;
   TBranch        *b_truemompi0ome;
   TBranch        *b_truedalitzpi1pi2ome;
   TBranch        *b_truedalitzpi1pi0ome;
   TBranch        *b_truecosthome;
   TBranch        *b_mom1piome;
   TBranch        *b_mom2piome;
   TBranch        *b_mompi0ome;
   TBranch        *b_dalitzpi1pi2ome;
   TBranch        *b_dalitzpi1pi0ome;
   TBranch        *b_costhome;

    //eta stuff
   TBranch        *b_nrecoeta;
   TBranch        *b_indexbesteta;
   TBranch        *b_barembesteta;
   TBranch        *b_mm2besteta;
   TBranch        *b_modeeta;
   TBranch        *b_metagg;
   TBranch        *b_mm2etagg;
   TBranch        *b_metapppi0;
   TBranch        *b_mm2etapppi0;
   TBranch        *b_metapi0pi0pi0;
   TBranch        *b_mm2etapi0pi0pi0;
 

   //eta' stuff
   TBranch        *b_nrecoetap;
   TBranch        *b_indexbestetap;
   TBranch        *b_barembestetap;
   TBranch        *b_mm2bestetap;
   TBranch        *b_modeetap;
   TBranch        *b_etamassdauetap;
   TBranch        *b_rho0massdauetap;
   TBranch        *b_gammamomdauetap;
  
   //a0 stuff
   TBranch        *b_ma0;
   TBranch        *b_mm2a0;  
   TBranch        *b_a0massetadau;
   TBranch        *b_modea0;
  
   //a0p stuff
   TBranch        *b_ma0p;
   TBranch        *b_mm2a0p;  
   TBranch        *b_a0pmassetadau;
   TBranch        *b_modea0p;

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

   double themin;
   double themax;
   int thebins;
   char *thevar;   
   int multipl;
   double intdataen;
   double intdatadepl;
   double intMCen;
   double intMCdepl;   
   double SHIFTNEUT;
   double SIGMANEUT;
   char xname[100];
   TPad *fPads[50];
   TCanvas *c0, *c1, *c2, *c3, *c4, *c5, *c6, *c7;
   theanal(char *var, char *dir, double mi, double ma, int b);
   ~theanal();
   Int_t  Cut(Int_t entry);
   Int_t  GetEntry(Int_t entry);
   Int_t  LoadTree(Int_t entry);

   double    getBsysweight(int,int);
   double    getDsysweight(int,int,int,int,int,int,int);
   double    FermiWeight(double,double,double);
   double    fermi(double,double,double);

   void   Init(TTree *tree);
   void   Loop(int nevents);
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
   void   overlap(int cut, double relnorm, TString dir);
   void readCuts(TString filename);
   void dumpCuts();
   double chisq(TH1 *h1, TH1 *h2);
   Double_t  smeargauss(double,double,double);

   TH2D mydalitz;

};

#endif

#ifdef theanal_cxx


#endif // #ifdef theanal_cxx

