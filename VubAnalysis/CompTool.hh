//////////////////////////////////////////////////////////
//   This class has been automatically generated 
//     (Mon May 13 08:12:21 2002 by ROOT version3.01/06)
//   from TTree events/events
//   found on file: /nfs/farm/babar/AWG7/ISL/tmp/rootfitfiles/sx-allgeneric.root
//////////////////////////////////////////////////////////


#ifndef CompTool_hh
#define CompTool_hh

#include <TROOT.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include "util.hh"
#include "../RecoilAnalysis/recoilAnalysis.hh" 
#include "fitNtp.hh" 

class CompTool {
   public :
   TTree          *fChain;   //pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //current Tree number in a TChain
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
   Int_t           nnpi0;
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
   Double_t        fkplus;
   Int_t           GfD0Ds;
   Int_t           GfDDs;
   Int_t           Gvxbtyp;
   Int_t           GSem;
   Int_t           GfDpi;
   Int_t           GfDpiz;
   Int_t           GfDk;
   Int_t           GfDks;
   Int_t           GfDlep;
   Int_t           GfDgam;
   UChar_t         GoodEvent;
   UChar_t         isDupli;
   Int_t           ValMap;
   Int_t           vub;
   Int_t           vcb;
   Int_t           vxbtyp;
   Int_t           other;
   Int_t           bgcat;
   Int_t           xcharge;
   Double_t        pxhad;
   Double_t        txhad;
   Double_t        fxhad;
   Double_t        exhad;
   Double_t        mxhad;
   Double_t        gmax;
   Double_t        mxhadfit;
   Int_t           lcharge;
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
   Int_t           nchg;
   Int_t           nneu;
   Int_t           nneu80_160;
   Int_t           nneu160_320;
   Int_t           nneufromB;
   Int_t           nneufromB80_160;
   Int_t           nneufromB160_320;
   Int_t           nkp;
   Int_t           nks;
   Int_t           npi0;
   Double_t        dx;
   Double_t        dy;
   Double_t        dz;
   Double_t        s2dxx;
   Double_t        s2dyy;
   Double_t        s2dzz;
   Double_t        s2dxy;
   Double_t        s2dyz;
   Double_t        s2dxz;
   Double_t        pnu;
   Double_t        tnu;
   Double_t        fnu;
   Double_t        eneu;
   Double_t        epiz;
   Double_t        kminmom;
   Double_t        kmaxmom;
   Double_t        mm2;
   Double_t        mm2nc;
   Double_t        mm2fit;
   Double_t        allksm0[18];
   Double_t        allksp[18];
   UChar_t         allksmc[18];
   Double_t        allchkp[7];
   UChar_t         allchkmc[7];
   Double_t        m0ks;
   Double_t        pks;
   Double_t        pksmc;
   Int_t           ntkl;
   Double_t        tklp[0];
   Double_t        tklth[0];
   Double_t        tklph[0];
   UChar_t         tklisol[0];
   Int_t           ntks;
   Double_t        tksp[0];
   Double_t        tksth[0];
   Double_t        tksph[0];
   Int_t           tksdec[0];
   Int_t           ntchk;
   Double_t        tchkp[0];
   Double_t        tchkth[0];
   Double_t        tchkph[0];
   Int_t           nklres;
   Double_t        klresth[0];
   Double_t        klresph[0];
   Int_t           klid[0];
   UChar_t         klcone[0];
   Double_t        emckl;
   Double_t        emckl0;
   Double_t        emckl22;
   Double_t        mxks;
   Double_t        mm2ks;
   Double_t        mxksfit;
   Double_t        mm2misk;
   Double_t        mxmisk;
   Double_t        mxmiskfit;
   Double_t        mm2mchk;
   Double_t        mxchk;
   Double_t        mxchkfit;
   Double_t        totweight;
   Double_t        totweightNutMult;
   Double_t        totweightTrkMult;

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
   TBranch        *b_nnpi0;
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
   TBranch        *b_kplus;
   TBranch        *b_Gvxbtyp; 
   TBranch        *b_GfDDs;  
   TBranch        *b_GfD0Ds;  
   TBranch        *b_GSem;    
   TBranch        *b_GfDpi;   
   TBranch        *b_GfDpiz;  
   TBranch        *b_GfDk;    
   TBranch        *b_GfDks;   
   TBranch        *b_GfDlep;  
   TBranch        *b_GfDgam;  
   TBranch        *b_GoodEvent;
   TBranch        *b_isDupli;
   TBranch        *b_ValMap;
   TBranch        *b_vub;
   TBranch        *b_vcb;
   TBranch        *b_vxbtyp;
   TBranch        *b_other;
   TBranch        *b_bgcat;
   TBranch        *b_xcharge;
   TBranch        *b_pxhad;
   TBranch        *b_txhad;
   TBranch        *b_fxhad;
   TBranch        *b_exhad;
   TBranch        *b_mxhad;
   TBranch        *b_gmax;
   TBranch        *b_mxhadfit;
   TBranch        *b_lcharge;
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
   TBranch        *b_nchg;
   TBranch        *b_nneu;
   TBranch        *b_nneu80_160;
   TBranch        *b_nneu160_320;
   TBranch        *b_nneufromB;
   TBranch        *b_nneufromB80_160;
   TBranch        *b_nneufromB160_320;
   TBranch        *b_nkp;
   TBranch        *b_nks;
   TBranch        *b_npi0;
   TBranch        *b_dx;
   TBranch        *b_dy;
   TBranch        *b_dz;
   TBranch        *b_s2dxx;
   TBranch        *b_s2dyy;
   TBranch        *b_s2dzz;
   TBranch        *b_s2dxy;
   TBranch        *b_s2dyz;
   TBranch        *b_s2dxz;
   TBranch        *b_pnu;
   TBranch        *b_tnu;
   TBranch        *b_fnu;
   TBranch        *b_eneu;
   TBranch        *b_epiz;
   TBranch        *b_kminmom;
   TBranch        *b_kmaxmom;
   TBranch        *b_mm2;
   TBranch        *b_mm2nc;
   TBranch        *b_mm2fit;
   TBranch        *b_allksm0;
   TBranch        *b_allksp;
   TBranch        *b_allksmc;
   TBranch        *b_allchkp;
   TBranch        *b_allchkmc;
   TBranch        *b_m0ks;
   TBranch        *b_pks;
   TBranch        *b_pksmc;
   TBranch        *b_ntkl;
   TBranch        *b_tklp;
   TBranch        *b_tklth;
   TBranch        *b_tklph;
   TBranch        *b_tklisol;
   TBranch        *b_ntks;
   TBranch        *b_tksp;
   TBranch        *b_tksth;
   TBranch        *b_tksph;
   TBranch        *b_tksdec;
   TBranch        *b_ntchk;
   TBranch        *b_tchkp;
   TBranch        *b_tchkth;
   TBranch        *b_tchkph;
   TBranch        *b_nklres;
   TBranch        *b_klresth;
   TBranch        *b_klresph;
   TBranch        *b_klid;
   TBranch        *b_klcone;
   TBranch        *b_emckl;
   TBranch        *b_emckl0;
   TBranch        *b_emckl22;
   TBranch        *b_mxks;
   TBranch        *b_mm2ks;
   TBranch        *b_mxksfit;
   TBranch        *b_mm2misk;
   TBranch        *b_mxmisk;
   TBranch        *b_mxmiskfit;
   TBranch        *b_mm2chk;
   TBranch        *b_mxchk;
   TBranch        *b_mxchkfit;
   TBranch        *btotweight;
   TBranch        *btotweightNutMult;
   TBranch        *btotweightTrkMult;

  TLatex tl;
  TPad *fPads[5];
  recoilAnalysis b;
  double themin;
  double themax;
  int thebins;
  const char *thevar;   
  int multipl;
  double intdataen;
  double intdatadepl;
  double intMCen;
  double intMCdepl;   
  double SHIFTNEUT;
  double SIGMANEUT;
  TCanvas *c0, *c1, *c2, *c3, *c4, *c5, *c6, *c7, *c8;
  //   CompTool(TTree *tree=0);
  ~CompTool();
  CompTool(const char *var, double mi, double ma, int b, int Sys);
  Int_t  Cut(Int_t entry);
  Int_t  GetEntry(Int_t entry);
  Int_t  LoadTree(Int_t entry);
  void   Init(TTree *tree);
  void   Loop(int nevents, int cat, double shift = 0, double smear = 0, int isbch = 2, int multcat = 7, int seed = 1);
  int dImode;
  recoilDSys *Dvar;
  recoilDSys *Bsem;
  void sameP(TString dir, TString flag);
  // seed = 2 dstar
  // seed = 3 dstar0
  // seed = 4 dc
  // seed = 5 d0
  void   Bookhist(); TFile *fHistFile;
  Bool_t Notify();
  void   Show(Int_t entry = -1);
  int    hist(double mx);
  void   Fitmes(int cat, int cut, TString flag);
  void   sighisto(double& signal, double& signalErr,TH1D *histo, double &resmean, double &ressigma, double &resalpha, double &resn, int fixpar,double mean, double sigma, double alpha, double n, double argus);
  void   overlap(int cut, int norm, TString dir, TString flag);
  void   effplots(TString dir, TString flag);
  double chisq(TH1 *h1, TH1 *h2, TString outfile);
  double chisqerr(TH1 *h1, TH1 *h2, char * outfile);
  Double_t  smeargauss(double,double,double);
  void      dumpCuts();
  void      readCuts(TString filename, int dump = 1);
  void      readintpur();
  void      initRest();
  double    getBsysweight(int,int);
  double    getDsysweight(int,int,int,int,int,int,int);
  Double_t theweighttrk[20], theweightneu[20];
  Double_t brecosig[6000], brecobkg[6000], brecointpur[6000];  
  TString fCutFile; 
  Double_t DESIGNALLO,  DESIGNALHI, MESSIGNALBANDLO, MESSIGNALBANDHI;
  Double_t MXMAXCUT, MM2CUT,
    DFROMD; //0 --> all 1-->D fromD 2--> DfromB
  Double_t INTPURITY, CUTPCMS;
  Int_t    DOBDECWEIGHT, DODDECWEIGHT; 
};

#endif

