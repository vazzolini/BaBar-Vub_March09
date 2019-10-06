#ifndef exclfitNtp_h
#define exclfitNtp_h

#include <iostream>
#include <vector>

#include <TString.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h> 
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TLatex.h>
#include <TVector2.h>
#include "TCanvas.h"
#include "TGraph.h"
#include "TMinuit.h"
// #include "TPostScript.h"

#include "RooFitCore/RooHist.hh"

#include "RecoilAnalysis/recoilAnalysis.hh" 
#include "RecoilAnalysis/recoilDSys.hh"
#include "RecoilAnalysis/recoilBuSys.hh"
#include "RecoilAnalysis/mesData.hh"

#include "RooFitCore/RooDataSet.hh"
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooGlobalFunc.hh" 

#include "VirVubFitter/ExclHelper.hh"


//#include "VirVubFitter/VirClass.hh"

class TF1;
class TH1D;
class TGraph;
class TFile;
class TPostScript;
class TCanvas;
class TPad;

class RooDataSet;
class RooRealVar;
class RooPlot;
class RooAbsPdf;

class exclfitNtp {

public :

  //new add
  enum fileTypeEnum {VubTotal=0, VubTotalres=1, VubTotalnres=2, Vcb=3, Vcb1=4, Vcb2=5, Data=6, VubTruthres=7, VubTruthnres=8};

  TString fWeightFile;

  //end new add

  TTree *fChain;   //pointer to the analyzed TTree or TChain
  Int_t fCurrent; //current Tree number in a TChain
  TString fOptionFile; 
  TString fCutFile; 

  TString texPrefix; 
  
  // ----------------------------------------------------------------------
  // -- CUTS

  //  Int_t NOTUNBINNED; from incl serve ???

  Double_t  PSTARFACT, LEPTONPCUT, ELECTRONPCUT, MUONPCUT, PRMM2CUT, MM2, Q2LOWCUT, Q2HIGHCUT, Q2CORR;
  Double_t MNUSQLOW, MNUSQHIGH,  CHLOW, CHHIGH, DEPL, BTYPE, LEPTTYPE, MAXINTPUR, MININTPUR, MINPUR, RUN;
  Int_t USECB, GAUSSFIT, FIXMEANVALUE, FIXSIGMA, FIXARGUS1, FIXARGUS2, NLEP;
  Int_t FITCATEGORY, FIXCB1, FIXCB2, MIXCORR, FITOPT, BLINDING;
  Double_t  RANDOMSEED, BLINDSIZE;
  //  Int_t DOBRECOWEIGHT,  DOTRKWEIGHT, DONEUWEIGHT;
  Int_t   DOBDECWEIGHT, DODDECWEIGHT, DOBUDECWEIGHT;
  Double_t RATIOBR;
  Int_t   SPTYPE;
  Double_t MXCUTLOWEXCL, MXCUTHIGHEXCL, DELTAM;
  Double_t MXCUTLOWEXCL1, MXCUTHIGHEXCL1, MXCUTLOWEXCL2, MXCUTHIGHEXCL2, MXCUTLOWEXCL3, MXCUTHIGHEXCL3, MXCUTLOWEXCL4, MXCUTHIGHEXCL4;
  Double_t NCHGLOWEXCL, NCHGHIGHEXCL;
  Double_t NCOMBLOWEXCL, NCOMBHIGHEXCL;
  Double_t NPI0LOWEXCL, NPI0HIGHEXCL;
  Double_t MOM1MIN, MOM2MIN;
  Double_t KSELE;
  Double_t MAXENEU,JPSIWIN;
  Double_t DALITZCUT;
  Double_t DAUETAMASS, DAURHOMASS, DAUGAMMAMOM;
  Double_t DAUETAMASS2, DAUETAMASS3, DAUETAMASS4;
  Int_t DORIGHTNCHG;
  Bool_t DOTHEO;
  Double_t MNUSQPI0LOW, MNUSQPI0HIGH, MNUSQETALOW, MNUSQETAHIGH, MNUSQRHOLOW, MNUSQRHOHIGH, MNUSQRHO0LOW, MNUSQRHO0HIGH, MNUSQOMEGALOW, MNUSQOMEGAHIGH;
  
  // ----------------------------------------------------------------------
  // -- FILES
  TString FILEVUBTOTAL, FILEVUBTOTALRES, FILEVUBTOTALNRES;
  TString FILEVCB, FILEDATA, PREFIXOUT, DIRNAME;

  
  //Added
  //Declaration of leaves types
  Int_t           run;
  Int_t           lower;
  Int_t           upper;
  Double_t        bmass;
  Double_t        bmassfit;
  UChar_t         sbox;
  Double_t        mes;
  Double_t        de;
  Double_t        pur;
  Double_t        intpur;
  Int_t           mode;
  Int_t           brecoflav;
  Int_t           brecocharge;
  Int_t           brecomc;
  Double_t        mxhadgen;
  Double_t        pcmsgen;
  Double_t        tcmsgen;
  Double_t        fcmsgen;
  Double_t        ecmsgen;
  Double_t        pxhadgen;
  Double_t        txhadgen;
  Double_t        fxhadgen;
  Double_t        exhadgen;
  // Stuff for MC reweighting
  Double_t        mcreweight;
  Double_t        eBrecoilgen;
  Double_t        pxBrecoilgen;
  Double_t        pyBrecoilgen;
  Double_t        pzBrecoilgen;
  Double_t        mBrecoilgen;
  Double_t        eLeptongen;
  Double_t        pxLeptongen;
  Double_t        pyLeptongen;
  Double_t        pzLeptongen;
  Double_t        mLeptongen;
  Double_t        eMesongen;
  Double_t        pxMesongen;
  Double_t        pyMesongen;
  Double_t        pzMesongen;
  Double_t        mMesongen;
  Int_t           vub;
  Int_t           vcb;
  Int_t           vxbtyp;
  Int_t           other;
  Int_t           xcharge;
  Double_t        pxhad;
  Double_t        txhad;
  Double_t        fxhad;
  Double_t        exhad;
  Double_t        mxhad;
  Double_t        mxhadfit;
  Double_t        mxhadchg;
  Double_t        q2Gen;
  Int_t           lcharge;
  Int_t           isele;
  Double_t        plab;
  Double_t        tlab;
  Double_t        flab;
  Double_t        pcms;
  Double_t        tcms;
  Double_t        fcms;
  Double_t        ecms;
  Int_t           nle;
  Int_t           nel;
  Int_t           nmu;
  Int_t           nlept500;
  Int_t           nelec500;
  Int_t           nmu500;
  Int_t           nchg;
  Int_t           nneu;
  Int_t           nkp;
  Int_t           nks;
  Int_t           npi0;
  Int_t           GSem;
  Int_t           GfDpi;
  Int_t           GfDpiz;
  Int_t           GfDk;
  Int_t           GfDks;
  Int_t           GfDlep;
  Int_t           GfDgam;
  Double_t        pnu;
  Double_t        tnu;
  Double_t        fnu;
  Double_t        eneu;
  Double_t        Eneualt;
  Double_t        epiz;
  Double_t        mm2;
  Double_t        mm2nc;
  Double_t        mm2fit;
  Double_t        deltam;
  Double_t        wdeltam;
  Double_t        totweight;
  Double_t        totweightNutMult;
  Double_t        totweightTrkMult;
  Int_t           Gvxbtyp;
  Double_t        ass_deltapB;           


  //pi stuff    
  Int_t           nrecoPi;  
  Int_t           indexbestPi; 
  Int_t           chbestPi;  
  Double_t        mm2bestPi;
  Double_t        q2bestPi;
  Float_t         baremassjpsiPi[100];
  Float_t         mm2Pi[100];
  Float_t         q2Pi[100];
  Int_t           chPi[100];  

   //pi0 stuff
  Int_t           nrecoPi0;   
  Int_t           indexbestPi0;  
  Double_t        barembestPi0;
  Double_t        mm2bestPi0;
  Double_t        q2bestPi0;
  Double_t        mm2gamma;
  Double_t        truemom1phpi0;
  Double_t        truemom2phpi0;
  Double_t        truemomlab1phpi0;
  Double_t        truemomlab2phpi0;
  Double_t        trueth1phpi0;
  Double_t        trueth2phpi0;
  Float_t         Estar1dauPi0[100];
  Float_t         Estar2dauPi0[100];
  Float_t         mm2Pi0[100];
  Float_t         q2Pi0[100];
  Float_t         baremPi0[100];

   //rho stuff;
  Int_t           nrecoRho; 
  Int_t           indexbestRho;
  Double_t        barembestRho; 
  Double_t        mm2bestRho; 
  Double_t        q2bestRho; 
  Float_t         Estar1dauRho[100]; 
  Float_t         Estar2dauRho[100]; 
  Float_t         mm2Rho[100]; 
  Float_t         q2Rho[100]; 
  Float_t         baremRho[100];
  Float_t         PimomdauRho[100]; 
  Float_t         Pi0momdauRho[100];
  Double_t        truemrho;
  Double_t        truemompirho;
  Double_t        truemompi0rho;
     
    //rho0 stuff;
  Int_t           nrecoRho0; 
  Int_t           indexbestRho0;
  Double_t        barembestRho0; 
  Double_t        mm2bestRho0; 
  Double_t        q2bestRho0; 
  Float_t         Estar1dauRho0[100];  
  Float_t         Estar2dauRho0[100];  
  Float_t         mm2Rho0[100];  
  Float_t         q2Rho0[100];  
  Float_t         baremRho0[100]; 
  Double_t        truemrho0;
  Double_t        truemom1pirho0;
  Double_t        truemom2pirho0;

   //omega stuff
  Int_t           nrecoOmega; 
  Int_t           indexbestOmega;
  Double_t        barembestOmega;  
  Double_t        mm2bestOmega;  
  Double_t        q2bestOmega;
  Float_t         Estar1dauOmega[100];   
  Float_t         Estar2dauOmega[100];
  Float_t         Estar3dauOmega[100];   
  Float_t         mm2Omega[100];   
  Float_t         q2Omega[100];   
  Float_t         baremOmega[100];     
  Float_t         Pi1momdauOmega[100];
  Float_t         Pi2momdauOmega[100];  
  Float_t         Pi0momdauOmega[100]; 
  Double_t        truemomega;
  Double_t        truemom1piome;
  Double_t        truemom2piome;
  Double_t        truemompi0ome;
  Double_t        truedalitzpi1pi2ome;
  Double_t        truedalitzpi1pi0ome;
  Double_t        truecosthome;
  Double_t        dalitzpi1pi2ome;
  Double_t        dalitzpi1pi0ome;
  Double_t        costhome;

   //eta stuff
  Int_t           nrecoEta; 
  Int_t           indexbestEta;
  Double_t        q2bestEta;   
  Double_t        metagg;
  Double_t        mm2etagg;
  Double_t        metapppi0;
  Double_t        mm2etapppi0;
  Double_t        metapi0pi0pi0;
  Double_t        mm2etapi0pi0pi0;
  Double_t        barembestEta;
  Double_t        mm2bestEta;
  Int_t           modeEta[100];    
  Float_t         mm2Eta[100]; 
  Float_t         q2Eta[100]; 
  Float_t         baremEta[100];

   //eta' stuff
  Int_t           nrecoEtap;
  Int_t           indexbestEtap;
  Double_t        q2bestEtap;   
  Double_t        barembestEtap;
  Double_t        mm2bestEtap;
  Int_t           modeEtap[100];
  Float_t         GammamomdauEtap[100];
  Float_t         EtamassdauEtap[100];
  Float_t         Rho0massdauEtap[100];    
  Float_t         mm2Etap[100]; 
  Float_t         q2Etap[100]; 
  Float_t         baremEtap[100];

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
  TBranch        *b_run;
  TBranch        *b_lower;
  TBranch        *b_upper;
  TBranch        *b_bmass;
  TBranch        *b_bmassfit;
  TBranch        *b_sbox;
  TBranch        *b_mes;
  TBranch        *b_de;
  TBranch        *b_pur;
  TBranch        *b_intpur;
  TBranch        *b_mode;
  TBranch        *b_brecoflav;
  TBranch        *b_brecocharge;
  TBranch        *b_brecomc;
  TBranch        *b_mxhadgen;
  TBranch        *b_pcmsgen;
  TBranch        *b_tcmsgen;
  TBranch        *b_fcmsgen;
  TBranch        *b_ecmsgen;
  TBranch        *b_pxhadgen;
  TBranch        *b_txhadgen;
  TBranch        *b_fxhadgen;
  TBranch        *b_exhadgen;
  TBranch        *b_mcreweight;
  TBranch        *b_ebrecoilgen;
  TBranch        *b_pxbrecoilgen;
  TBranch        *b_pybrecoilgen;
  TBranch        *b_pzbrecoilgen;
  TBranch        *b_mbrecoilgen;
  TBranch        *b_eleptongen;
  TBranch        *b_pxleptongen;
  TBranch        *b_pyleptongen;
  TBranch        *b_pzleptongen;
  TBranch        *b_mleptongen;
  TBranch        *b_emesongen;
  TBranch        *b_pxmesongen;
  TBranch        *b_pymesongen;
  TBranch        *b_pzmesongen;
  TBranch        *b_mmesongen;
  TBranch        *b_vub;
  TBranch        *b_vcb;
  TBranch        *b_vxbtyp;
  TBranch        *b_other;
  TBranch        *b_xcharge;
  TBranch        *b_pxhad;
  TBranch        *b_txhad;
  TBranch        *b_fxhad;
  TBranch        *b_exhad;
  TBranch        *b_mxhad;
  TBranch        *b_mxhadfit;
  TBranch        *b_mxhadchg;
  TBranch        *b_q2Gen;
  TBranch        *b_lcharge;
  TBranch        *b_isele;
  TBranch        *b_plab;
  TBranch        *b_tlab;
  TBranch        *b_flab;
  TBranch        *b_pcms;
  TBranch        *b_tcms;
  TBranch        *b_fcms;
  TBranch        *b_ecms;
  TBranch        *b_nle;
  TBranch        *b_nel;
  TBranch        *b_nmu;
  TBranch        *b_nlept500;
  TBranch        *b_nelec500;
  TBranch        *b_nmu500;
  TBranch        *b_nchg;
  TBranch        *b_nneu;
  TBranch        *b_nkp;
  TBranch        *b_nks;
  TBranch        *b_npi0;
  TBranch        *b_GSem;    
  TBranch        *b_GfDpi;   
  TBranch        *b_GfDpiz;  
  TBranch        *b_GfDk;    
  TBranch        *b_GfDks;   
  TBranch        *b_GfDlep;  
  TBranch        *b_GfDgam;  
  TBranch        *b_pnu;
  TBranch        *b_tnu;
  TBranch        *b_fnu;
  TBranch        *b_eneu;
  TBranch        *b_eneualt;
  TBranch        *b_epiz;
  TBranch        *b_kminmom;
  TBranch        *b_kmaxmom;
  TBranch        *b_mm2;
  TBranch        *b_mm2nc;
  TBranch        *b_mm2fit;
  TBranch        *b_deltam;
  TBranch        *b_wdeltam;
  TBranch        *btotweight;
  TBranch        *btotweightNutMult;
  TBranch        *btotweightTrkMult;
  TBranch        *b_Gvxbtyp;
  TBranch        *b_ass_deltapB;     

   //pi stuff    
  TBranch        *b_nrecopi;
  TBranch        *b_indexbestpi;
  TBranch        *b_chbestpi;
  TBranch        *b_mm2bestpi;
  TBranch        *b_baremassjpsipi;
  TBranch        *b_mm2pi;
  TBranch        *b_q2pi;
  TBranch        *b_chpi;
  TBranch        *b_q2bestpi;
  //TBranch        *b_lepmap;

   //pi0 stuff
  TBranch        *b_nrecopi0;
  TBranch        *b_indexbestpi0;
  TBranch        *b_barembestpi0;
  TBranch        *b_mm2bestpi0;
  TBranch        *b_mm2gamma;
  TBranch        *b_truemom1phpi0;
  TBranch        *b_truemom2phpi0;
  TBranch        *b_truemomlab1phpi0;
  TBranch        *b_truemomlab2phpi0;
  TBranch        *b_trueth1phpi0;
  TBranch        *b_trueth2phpi0;
  TBranch        *b_estar1daupi0;
  TBranch        *b_estar2daupi0;
  TBranch        *b_mm2pi0;
  TBranch        *b_q2pi0;
  TBranch        *b_barempi0;
  TBranch        *b_q2bestpi0;

   //rho stuff;
  TBranch        *b_nrecorho;
  TBranch        *b_indexbestrho;
  TBranch        *b_barembestrho;
  TBranch        *b_mm2bestrho;
  TBranch        *b_q2bestrho;
  TBranch        *b_estar1daurho; 
  TBranch        *b_estar2daurho;
  TBranch        *b_mm2rho;
  TBranch        *b_q2rho;    
  TBranch        *b_baremrho;
  TBranch        *b_pimomdaurho; 
  TBranch        *b_pi0momdaurho; 
  TBranch        *b_truemrho;
  TBranch        *b_truemompirho;
  TBranch        *b_truemompi0rho;
     
    //rho0 stuff;
  TBranch        *b_nrecorho0; 
  TBranch        *b_indexbestrho0; 
  TBranch        *b_barembestrho0; 
  TBranch        *b_mm2bestrho0; 
  TBranch        *b_q2bestrho0; 
  TBranch        *b_estar1daurho0;  
  TBranch        *b_estar2daurho0; 
  TBranch        *b_mm2rho0; 
  TBranch        *b_q2rho0;     
  TBranch        *b_baremrho0; 
  TBranch        *b_truemrho0; 
  TBranch        *b_truemom1pirho0;
  TBranch        *b_truemom2pirho0;

   //omega stuff
  TBranch        *b_nrecoomega; 
  TBranch        *b_indexbestomega; 
  TBranch        *b_barembestomega; 
  TBranch        *b_mm2bestomega; 
  TBranch        *b_q2bestomega; 
  TBranch        *b_estar1dauomega;  
  TBranch        *b_estar2dauomega;
  TBranch        *b_estar3dauomega; 
  TBranch        *b_mm2omega; 
  TBranch        *b_q2omega;     
  TBranch        *b_baremomega; 
  TBranch        *b_pi1momdauomega;
  TBranch        *b_pi2momdauomega;   
  TBranch        *b_pi0momdauomega;  
  TBranch        *b_truemomega;
  TBranch        *b_truemom1piome;
  TBranch        *b_truemom2piome;
  TBranch        *b_truemompi0ome;
  TBranch        *b_truedalitzpi1pi2ome;
  TBranch        *b_truedalitzpi1pi0ome;
  TBranch        *b_truecosthome;
  TBranch        *b_dalitzpi1pi2ome;
  TBranch        *b_dalitzpi1pi0ome;
  TBranch        *b_costhome;

   //eta stuff
  TBranch        *b_nrecoeta;
  TBranch        *b_indexbesteta;
  TBranch        *b_metagg;
  TBranch        *b_mm2etagg;
  TBranch        *b_metapppi0;
  TBranch        *b_mm2etapppi0;
  TBranch        *b_metapi0pi0pi0;
  TBranch        *b_mm2etapi0pi0pi0;
  TBranch        *b_baremeta; 
  TBranch        *b_barembesteta;
  TBranch        *b_mm2besteta;
  TBranch        *b_modeeta;
  TBranch        *b_q2eta;      
  TBranch        *b_q2besteta;
  TBranch        *b_mm2eta;

   //eta' stuff
  TBranch        *b_nrecoetap;
  TBranch        *b_indexbestetap;
  TBranch        *b_baremetap; 
  TBranch        *b_barembestetap;
  TBranch        *b_mm2bestetap;
  TBranch        *b_modeetap;
  TBranch        *b_etamassdauetap;
  TBranch        *b_rho0massdauetap;
  TBranch        *b_gammamomdauetap;
  TBranch        *b_q2etap;  
  TBranch        *b_q2bestetap;   
  TBranch        *b_mm2etap;

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

  // end variables


  RooPlot *xframe;
  RooPlot *pullframe;
  TPad *pPad;
  TPad *mPad;

  //Plots and canvases

  TLatex tl;

  //Mes Fit
  // Mes fit for exclusive Vub 
  mesData* exclVubMesUnb(RooDataSet *data, RooRealVar *x,  double &resmean, double &ressigma, double &resalpha, double &resn, double &resargpar,int print, int func, double mean, double sigma, double alpha, double n, double argus, double nSIG, double nBKG, char * simply); 
 
  TVector2 exclsighistounb(RooDataSet *Adata, RooRealVar *Ax, double &Aresmean, double &Aressigma, double &Aresalpha, double &Aresn, double &Aresargpar, double Amean, double Asigma, double Aalpha,  double An, double Aargus, double AnSIG, double AnBKG, char * Asimply, int fixpar); 

  //mes fit
  double dBinomial(double, double); 
  TString fmesFile;
  //double dstlnuFF(double r1,double r2,double rho2);
   TString funfFile;

  //  int nB; ???

  //mes parameters
  double mesNsl[4], mesdatacuts[4], mesvubcuts[4], mesvcbcuts[4], mesothcuts[4], mesNslMC[4], mesvcbMC[4], mesvubMC[4],messigleptcuts[4],messigcuts[4];


  exclfitNtp(TTree *tree, int Sys, TString filename);
  ~exclfitNtp();
  TFile *fHistFile;
  TPostScript *fPostScriptFile;
  Int_t  Cut(Int_t entry);
  Int_t  GetEntry(Int_t entry);
  Int_t  LoadTree(Int_t entry);
  void openEpsFile(TString name);
  void closeEpsFile();
  TFile* openHistFile(TString name);
  void closeHistFile();
  void Init(TTree *tree);
  void Bookhist();
  void Loop(int isdata, int icat, int nevents, int isMC, int nres);
  void theFit();
  TVector2 mixCorr(RooDataSet *data, char* cuts, double* mespar, bool iseps = 0, char *nameeps = "");
  TVector2 nomixCorr(RooDataSet *data, char* cuts, double* mespar, bool iseps = 0, char *nameeps = "");
  void makesubplot(RooDataSet *data, char* var, char* cuts, double* mespar, int nbins, double xmin, double xmax, char* outputplot, bool iseps = 0, char *nameeps = "");
  void makefinalplots(char *var, int nbins);
  Bool_t Notify();
  void Show(Int_t entry = -1);
  void initRest(TString filename);
  void readOptions(TString filename, int dump = 1);
  void dumpOptions();
  void readmesParam(TString filename, int dump = 1);
  void dumpmesParam();
  void readCuts(TString filename, int dump = 1);
  void dumpCuts();
  void readpstarfactor();

//read datasets
  //int readDataFile(RooDataSet**,TString);

//write datasets
//  void writeDataFile(const std::vector<RooDataSet*>&);

  //Dataset file handlers
  //void openDataSetFile(TString);
  //void closeDataSetFile();

//Functions for Fit parameters and files initializations
  TString getfileVubTotal();
  TString getfileVubTotalres();
  TString getfileVubTotalnres();
  TString getfileVcb();
  TString getfileData();
  bool * getfilechain();
  TChain * getchain(char * thechain, const TString* treeName = 0);
  void setPrefix(TString);
  void setTexPrefix(TString s) {texPrefix = s;}
  void setDirectory(TString);
  Int_t isBlind();
  double getblindfact();
  double getBsysweight(int,int);
  double getDsysweight(int,int,int,int,int,int,int);
  double getBusysweight(int);
//   double getTrackingWeight();
//   double getNeutralWeight();
//   double getBrecoWeight(double);
  double getPstarFactor(double);
 
  //Computation of reweighting factors (BR)
  //from inclusive
  //recoilDSys *Dvar;
  //recoilDSys *Bsem;
  //double getBsysweight(int decType,int thevub);
  //double getDsysweight(int decDpi,int decDk,int decDks,int decDpiz,int decDlep,int decImode,int thevub);

  //Computation of reweighting factors (detector)

  //double getTrackingWeight();
  //double getNeutralWeight();
  //double getBrecoWeight(double theintpur);
  //double FermiWeight(double kp, double deltamb, double deltaa);
  //int rHistq2(double amx);
  //int rHistel(double bmx);
  //int rHistmx(double cmx);
  //int newrHistmx(double cmx);
  //int newbinrHistmx(double cmx);
  //double MatrixW[896];
  //double newMatrixW[1024];
  //double fermi(double kp, double m, double a);


  //Old Reweighting
  double TrueMxWeight[42];
  double truehistbins[21];
  double getTrueMxWeight(double thetrumx, int index);
  int    TrueHist(double mxt);
  void   doTheo();

  RooRealVar *varMes;
  RooRealVar *varVar;  
  RooRealVar *varMm2;  
  RooRealVar *varLepYes;  
  RooRealVar *varFlavB;  
  RooRealVar *varAllnovar;  
  RooRealVar *varAllnomm2;  
  RooRealVar *varAllcuts;  
  RooRealVar *varWe;
  RooDataSet *datadata;
  RooDataSet *datamcsig;
  RooDataSet *datamcvcb;
  RooDataSet *datamcvub;
  RooDataSet *datamcoth;
  RooDataSet *datamc;
  RooDataSet *datamcvcbvub;
  RooDataSet *datamcsigforpstar;
  
  bool ischain[4];

  TCanvas* c1; 

  double correctionratiovub, correctionratiovcb,correctionratiooth;
  int dImode;
  recoilDSys *Dvar;
  recoilDSys *Bsem;
  recoilBuSys *Busem;
  double pstarfactor[16];
  double vcbmeanweight;
  double vubmeanweight;  
  double othmeanweight;
  double vcbmeanw[9], vubmeanw[9], othmeanw[9];
  int countvcb;
  int countvub;
  int countoth;   
  int covcb[9], covub[9], cooth[9];
  int vcbcounter;
  
  // fit variables
  double calcpstarfact, fact, epssig, errepssig, epsvub, errepsvub, epsvcb, errepsvcb, epsoth, errepsoth, nsl, nslmc, nslvub, nsigdata, errsigdata, nsig, nerrsig; 
  double ndatamm2,nvubdmm2,nvcbdmm2,nothdmm2,scalefact;
  
  // ranges for plots
  char VAR[100];
  int BINSMASS, BINSMM2, BINSPCMS, BINSVAR;
  double  MINMASS, MAXMASS, MINMM2, MAXMM2, MINPCMS, MAXPCMS, MINVAR, MAXVAR;
  
  // dalitz plot
  TH2D mydalitz;

  //perform the fit
  //  void theFit(int cmb,int multc,int su, int ck);
  //void compChisq(int cmb, int su); double chisq;   int NDOF;
  //  // dependence on VirFit for the unb mes fit
  //  exclfitNtp forfit;
};

#endif

#ifdef exclfitNtp_cxx
exclfitNtp::exclfitNtp(TTree *tree)
{
}

exclfitNtp::~exclfitNtp()
{
}

Int_t exclfitNtp::GetEntry(Int_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Int_t exclfitNtp::LoadTree(Int_t entry)
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



void exclfitNtp::Show(Int_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t exclfitNtp::Cut(Int_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef exclfitNtp_cxx

