//////////////////////////////////////////////////////////
//   This class has been automatically generated 
//     (Mon May 13 08:12:21 2002 by ROOT version3.01/06)
//   from TTree events/events
//   found on file: /nfs/farm/babar/AWG7/ISL/tmp/rootfitfiles/sx-allgeneric.root
//////////////////////////////////////////////////////////


#ifndef TheorTool_hh
#define TheorTool_hh

#include <TROOT.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TRandom.h>
#include <TCut.h>
#include <TH1.h>
#include <TLegend.h>
#include "util.hh"
#include "../RecoilAnalysis/recoilAnalysis.hh" 

class TheorTool {
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
   Int_t           fBchgen;
   Int_t           GfDpi;
   Int_t           GfDpiz;
   Int_t           GfDk;
   Int_t           GfDkmiss;
   Int_t           GfDks;
   Int_t           GfDkl;
   Int_t           GfDkspipi;
   Int_t           GfDkspiopio;
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
   TBranch        *b_fBchgen;   
   TBranch        *b_GfDpi;   
   TBranch        *b_GfDpiz;  
   TBranch        *b_GfDk;    
   TBranch        *b_GfDkmiss;   
   TBranch        *b_GfDks;   
   TBranch        *b_GfDkl;   
   TBranch        *b_GfDkspipi;   
   TBranch        *b_GfDkspiopio;   
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

  recoilAnalysis b;
  TCanvas *c1, *c2, *c3, *c5, *c4, *c6, *c7, *c0;
  double histbins[21]; double fittedBR; double nresMCfac;
  ~TheorTool();
  TheorTool(const char *var, double mi, int ma, const char * b, double fit, int seed);
  TH1 *sig;TH1 *bkg ;
  TLatex tl;
  Bool_t Notify();
  void Init(TTree *tree);
  void dumpCuts();
  void readCuts(TString filename, int dump = 1);
  void initRest();int dImode;
  void Loop(int nevents, int cat, int isbch, int doMbaRew);
  Double_t INTPURITY; Int_t CUTPCMS; TString fCutFile; 
  Double_t DESIGNALLO,  DESIGNALHI, MXMAXCUT, MM2CUT,
    DFROMD; //0 --> all 1-->D fromD 2--> DfromB
  void Bookhist();
  void overlap(int cut, TString dir, TString flag, int scal);
  void Fitmes(int cat, int cut, TString flag);
  const char *thevar;    double themin;  int themax;  int thebins;
  int theseed; int therew;
  int hist(double mx);
  void setHist(TH1 *h, Int_t color, Int_t symbol, Double_t size, Double_t width, Int_t style);
  void setLegend(TLegend *legend, Int_t style, Double_t bsize, Double_t tsize, Double_t fcolor, TH1 *fHi, TH1 *sHi, char * fTit, char * sTit, char * type);
  void   sighisto(double& signal, double& signalErr,TH1D *histo, double &resmean, double &ressigma, double &resalpha, double &resn, int fixpar,double mean, double sigma, double alpha, double n, double argus);
  //Mb and a reweighting
  double mb1;  double mb2;
  //Now the BRexcl becomes something that we can 'change'
  double BRexclBch; double BRexclB0; double TrueMxWeight[42];
  Double_t BMASS, BQMASS, A0;
  double a1;  double a2; double rnd;
  double Nnew1, Nnew2, Nnew3, Nnew4, Nold, Nsys;
  double fermi(double kp, double m, double a);
  void initMba(double dmb);
  float randomized(float mean,float sig);
  double reweightMba(double fKp, double mxHG, double mxHF, int sysWh);
  void computeWMba(double dmb);
  TRandom rndm;
  //Study on integrals
  void intStudy(int Avub, double Amxhadgen, int AGvxbtyp, int AfBchgen, int Acat);
  void PintStudy();
  double getTrueMxWeight(double thetrumx, int index);
};

#endif


