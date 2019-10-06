#ifndef fitNtp_h
#define fitNtp_h

#include <iostream.h>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h> 
#include <TF1.h>
#include <TH1.h>
#include <TTree.h>
#include <TLatex.h>
#include <TVector2.h>
#include "TCanvas.h"
#include "TGraph.h"
#include "TMinuit.h"
#include "TPostScript.h"

#include "RecoilAnalysis/recoilAnalysis.hh" 
#include "VubAnalysis/recoilDSys.hh"

class fitNtp {

   public :
   TTree          *fChain;   //pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //current Tree number in a TChain
   TString fOptionFile; 
   TString fCutFile; 
   TString fmesFile; 
   TString texPrefix; 

  // ----------------------------------------------------------------------
  // -- CUTS
  Double_t TOTALSTAT, TOTALSTATMODEL, BRRATIOGENVALUE,  BRRATIOVALUETAIL, MXCUT, Q2CUT, PSTARFACT, LEPTONPCUT, PRMM2CUT, CSILOCUT, CSIHICUT, Q2LOCUT, Q2HICUT, XLOCUT, XHICUT, WLOCUT, WHICUT, EWPWLOCUT, EWPWHICUT;
  Double_t MNUSQLOW, MNUSQHIGH,  CHLOW, CHHIGH, DEPL, BTYPE, LEPTTYPE, MAXINTPUR, MININTPUR, RUN, TOYMC;
  Double_t DELTAMB, MULFAC, DELTAA, FERMIAPP;
  Int_t USECB, FIXMEANVALUE, FIXSIGMA, FIXARGUS1, FIXARGUS2;
  Int_t FIXCB1, FIXCB2, FITTOTSHAPE, MIXCORR, FITMC, FITOPT, FITQ2, MULTIFIT, BLINDING;
  Int_t ISSMEARALL, ISSMEARBKG, CUTNNPI0;
  Double_t  VCBCOMP, OTHCOMP;
  Double_t  SMEARALLMEANVALUE, SMEARALLSIGMA, SMEARBKGMEANVALUE, SMEARBKGSIGMA;
  Double_t  RANDOMSEED, BLINDSIZE;
  Int_t DOBRECOWEIGHT, DOBDECWEIGHT, DODDECWEIGHT, DOTRKWEIGHT, DONEUWEIGHT, DOFERMI, DOTHEO;
  Int_t DOVARSTU;
  // fit results
  double dataFirstBin, dataErrFirstBin;
  double  vubcomp,errvubcomp,vcbcomp,errvcbcomp,othcomp,errothcomp;
  double  vubcompNOMC,errvubcompNOMC,vcbcompNOMC,errvcbcompNOMC,othcompNOMC,errothcompNOMC;
  double errS, S;double blindfactor;
  double vcbaftercutsbin1, errvcbaftercutsbin1, errfitotheraftercutsbin1;
  double errfitvcbaftercutsbin1, otheraftercutsbin1, errotheraftercutsbin1;
  double epsu,errepsu,epsmx ,errepsmx,epstot,errepstot;
  double tot, errtot, totmc; double fact, chisq;
  Double_t calcpstarfact, errcalcpstarfact, calcpstarfacttemp;
  double BRBR, errBRBR, errBRBRtheo, errtotalBRBR, errBRBRMCstat;
  double vubmcselected, vcbmcselected, othermcselected;
  double vubmc, errvubmc,vubmc2, vcbmc, errvcbmc;
  double correctionratiovub, correctionratiovcb;
  double vcbFirstBin, vcbErrFirstBin, othFirstBin, othErrFirstBin;
  Double_t vubmcaftercuts, vubmcaftercutsbin1,errvubmcaftercutsbin1,areavubmcaftercutsbin1,areavcbaftercutsbin1,areaothaftercutsbin1,vubmcaftertemp,vubmcleptforeff,vubmcallforeff,vubmcnocut,vcbmcnocut,othmcnocut;
  int NDOF;
  // ----------------------------------------------------------------------
  // -- FILES
  TString FILEVUBTOTAL, FILEVUBTOTALRES, FILEVUBTOTALNRES, FILEVUBDSTAR, FILEVUBDC, FILEVUBDSTAR0;
  TString FILEVUBD0, FILEVCB, FILEVCB2, FILEDATA, PREFIXOUT, DIRNAME;

  // mes parameters

  double mesNsl[4], mesdatacuts[4], mesvubcuts[4], mesvcbcuts[4], mesothcuts[4], mesNslMC[4], mesvcbMC[4], mesvubMC[4];
  double  mesvubMCmx[4], mesvubMCall[4], mesvcbMCall[4],mesvubMClepteff[4], mesvubMCalleff[4];

  //Study from Ric
  double MatrixW[896];

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
   Int_t           Gvxbtyp;
   Int_t           GSem;
   Int_t           GfDpi;
   Int_t           GfDpiz;
   Int_t           GfDk;
   Int_t           GfDks;
   Int_t           GfDlep;
   Int_t           GfDgam;
   Double_t        intpur;
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
   Double_t        csiCiuc;
   Double_t        wCiuc;
   Double_t        xCiuc;
   UChar_t         GoodEvent;
   UChar_t         isDupli;
   Int_t           ValMap;
   Int_t           vcb;
   Int_t           vub;
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
   Int_t           npi0;
   Int_t           nneu;
   Int_t           nneu80_160;
   Int_t           nneu160_320;
   Int_t           nkp;
   Int_t           nks;
   Double_t        totweight;
   Double_t        totweightNutMult;
   Double_t        totweightTrkMult;
   Double_t        enu;
   Double_t        pnu;
   Double_t        tnu;
   Double_t        fnu;
   Double_t        mm2;
   Double_t        mm2nc;
   Double_t        mm2fit;
   Double_t        deltam;
   Double_t        wdeltam;
   Double_t        q2;
   Double_t        q2nc;
   Double_t        q2fit;
   Double_t        q2Gen;
   Double_t        allksm0[14];
   Double_t        allksp[14];
   UChar_t         allksmc[14];
   Double_t        allchkp[8];
   UChar_t         allchkmc[8];
   Double_t        m0ks;
   Double_t        pks;
   Double_t        pksmc;
   Int_t           ntkl;
  /*
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
  */
   Double_t        xX;
   Double_t        yX;
   Double_t        zX;
   Double_t        v2xxX;
   Double_t        v2yyX;
   Double_t        v2zzX;
   Double_t        v2xyX;
   Double_t        v2yzX;
   Double_t        v2xzX;
   Double_t        EwPwfit;
   Double_t        xLep;
   Double_t        yLep;
   Double_t        zLep;
   Double_t        v2xxLep;
   Double_t        v2yyLep;
   Double_t        v2zzLep;
   Double_t        v2xyLep;
   Double_t        v2yzLep;
   Double_t        v2xzLep;
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
   Int_t           nnpi0;

   TF1* f1all;
   TF1* f1bkg;
   int dImode;
   recoilDSys *Dvar;
   recoilDSys *Bsem;
   double TrackingWeight[20];
   double NeutralWeight[20];
   double BrecoWeight[20];
   double TrueMxWeight[42];
  //   double TrueExMxWeight[2];
   double truehistbins[21];
   double mxBinning[10];
   double pstarfactor[16];

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
   TBranch        *b_EwPwfit;
   TBranch        *b_xLep;
   TBranch        *b_yLep;
   TBranch        *b_zLep;
   TBranch        *b_Gvxbtyp; 
   TBranch        *b_GSem;    
   TBranch        *b_GfDpi;   
   TBranch        *b_GfDpiz;  
   TBranch        *b_GfDk;    
   TBranch        *b_GfDks;   
   TBranch        *b_GfDlep;  
   TBranch        *b_GfDgam;  
   TBranch        *b_intpur;
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
   TBranch        *b_csiCiuc;
   TBranch        *b_wCiuc;
   TBranch        *b_xCiuc;
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
   TBranch        *b_npi0;
   TBranch        *b_nneu;
   TBranch        *b_nneu80_160;
   TBranch        *b_nneu160_320;
   TBranch        *b_nneufromB;
   TBranch        *b_nneufromB80_160;
   TBranch        *b_nneufromB160_320;
   TBranch        *b_nkp;
   TBranch        *b_nks;
   TBranch        *b_totweight;
   TBranch        *b_totweightNutMult;
   TBranch        *b_totweightTrkMult;
   TBranch        *b_enu;
   TBranch        *b_pnu;
   TBranch        *b_tnu;
   TBranch        *b_fnu;
   TBranch        *b_mm2;
   TBranch        *b_deltam;
   TBranch        *b_wdeltam;
   TBranch        *b_mm2nc;
   TBranch        *b_mm2fit;
   TBranch        *b_q2;
   TBranch        *b_q2nc;
   TBranch        *b_q2fit;
   TBranch        *b_q2Gen;
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
   TBranch        *b_nnpi0;
   TBranch        *b_v2xxLep;
   TBranch        *b_v2yyLep;
   TBranch        *b_v2zzLep;
   TBranch        *b_v2xyLep;
   TBranch        *b_v2yzLep;
   TBranch        *b_v2xzLep;
   TBranch        *b_xX;
   TBranch        *b_yX;
   TBranch        *b_zX;
   TBranch        *b_v2xxX;
   TBranch        *b_v2yyX;
   TBranch        *b_v2zzX;
   TBranch        *b_v2xyX;
   TBranch        *b_v2yzX;
   TBranch        *b_v2xzX;

  //End added

  TCanvas* c1; int fprlRew; 
   fitNtp(TTree *tree, int Sys, TString filename, int prl);
   ~fitNtp();
   TFile *fHistFile;
   TPostScript *fPostScriptFile;
   Int_t  Cut(Int_t entry);
   Int_t  GetEntry(Int_t entry);
   Int_t  LoadTree(Int_t entry);
   void    openEpsFile(TString name);
   void      closeEpsFile();
   TFile*    openHistFile(TString name);
   void      closeHistFile();
   void   Init(TTree *tree);
   void   Bookhist();
   void   Loop(int isdata, int icat, int nevents, int isMC, int nres);
   void   Fitmessub(const char* comp);
   void   theFit();
   void   theq2fit();
   void   effmult();
   int    hist(double mx);
   int    TrueHist(double mxt);
   int    histlept(double pstar);
   int    histq2(double theq2hi);
   int    histmqb(double mqb);

  //Ric Study/binning
  int rHistq2(double amx);
  int rHistel(double bmx);
  int rHistmx(double cmx);

   void   varStudy(char * varNa);
  //   void   mqbStudy(double vubcomp, double vcbcomp, double othcomp);
  //   void   csiStudy(double vubcomp, double vcbcomp, double othcomp);
   int    histcsi(double csi);
   int    histw(double Mew);
   int    histx(double Mex);
   TVector2 sighisto(TH1D *histo, double &resmean, double &ressigma, double &resalpha, double &resn, int fixpar,double mean, double sigma, double alpha, double n, double argus);
   Bool_t Notify();
   void   Show(Int_t entry = -1);
   void      initRest(TString filename);
   void      readOptions(TString filename, int dump = 1);
   void      dumpOptions();
   void      readmesParam(TString filename, int dump = 1);
   void      dumpmesParam();
   void      readCuts(TString filename, int dump = 1);
   void      dumpCuts();
   void      readWeights();
   void      readpstarfactor();
   TString   getfileVubTotal();
   TString   getfileVubTotalres();
   TString   getfileVubTotalnres();
   TString   getfileVubDstar();
   TString   getfileVubDc();
   TString   getfileVubDstar0();
   TString   getfileVubD0();
   TString   getfileVcb();
   TString   getfileVcb2();
   TString   getfileData();
   Int_t     isfitMC();
   Int_t     isfitOpt();
   Int_t     isfitQ2();
   Double_t  isToyMC();
   TH1D*     makeMesHist(const char * name, bool doToy);
   TH1D*     fitWithErrors(int ich=-99,int ineu=-99,char * typ="mx");
   TH1D*     fitWithoutErrors(int ich=-99,int ineu=-99);
   void      setPrefix(TString);
   void      setTexPrefix(TString s) {texPrefix = s;}
   void      setDirectory(TString);
   Double_t  totalStat();
   Double_t  totalStatModel();
   Double_t  BRRatioGenValue();
   Int_t     isBlind();
   Double_t  smeargauss(double,double,double);
   double    getblindfact();
   double    getBsysweight(int,int);
   double    getDsysweight(int,int,int,int,int,int,int);
   double    getTrackingWeight();
   double    getNeutralWeight();
   double    getBrecoWeight(double);
   double    getTrueMxWeight(double thetrumx, int index);
  //   double    getTrueExMxWeight(int ind);
   double    getPstarFactor(double);
   int       getMixcorr();
   void      doFermi();
   void      doTheo();
   void      doVarStu();
   double    FermiWeight(double,double,double);
   double    fermi(double,double,double);
  TGraph*    scanParameter(int parnum, int nsig, TMinuit &a, void (*func)(int &, double *, double &, double *, int)); 
  void DumpOFile(char *name, char* namet);
  void fillMesStu(double Mepcms,double Metheq2,double MeEwPwfit,double MecsiCiuc, double MexCiuc,double MewCiuc, char * Melgroup, char * Meflavcat, double Memes, double wei, char * Medata);
  void mesUtility(char * nameps,char * namevar, char * Xtit, double thecorrectionratio);
};

#endif

#ifdef fitNtp_cxx
fitNtp::fitNtp(TTree *tree)
{
}

fitNtp::~fitNtp()
{
}

Int_t fitNtp::GetEntry(Int_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Int_t fitNtp::LoadTree(Int_t entry)
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



void fitNtp::Show(Int_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t fitNtp::Cut(Int_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef fitNtp_cxx

