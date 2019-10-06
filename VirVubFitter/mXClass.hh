#ifndef mXClass_h
#define mXClass_h

#include <iostream>
#include <string>
#include <vector>
#include <utility>

#include <stdlib.h>

#include <TString.h>
#include <TTree.h>
#include <TChain.h>

#include <TROOT.h>
#include <TLatex.h>
#include <TVector2.h>
#include <TPaveText.h>
#include <TMinuit.h>

#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooFormulaVar.hh"

#include "RecoilAnalysis/recoilDSys.hh"
#include "RecoilAnalysis/mesData.hh"
#include "RecoilAnalysis/recoilAnalysis.hh"

#include "VirVubFitter/CMClass.hh"

class TF1;
class TH1D;
class TGraph;
class TFile;
class TPostScript;
class TCanvas;
class TPad;

class RooDataSet;
class RooPlot;
class RooAbsPdf;

class XSLBToDstrlnu_DstrToDpi_CLN;
// ################################################

class mXClass {

public :
  enum fileTypeEnum  { VubTotal=0, VubTotalres=1, VubTotalnres=2, Vcb=3, Vcb1=4, Vcb2=5, Data=6, VubTruthres=7, VubTruthnres=8 };

  enum parFitStdEnum { iMean=0, iSigma=1 , iAlpha=2, iN=3, iArgus=4 };

  enum yieldsEnum    { iNsig = 0, iNbkg, iNbpk };
  
  enum { GaussFit = 0, ArgusAndCB = 1};
  

  TTree          *fChain;   //pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //current Tree number in a TChain
  Int_t   fnewbin;
  TString fWeightFile; 
  TString fOptionFile; 
  TString texPrefix;
  TString FILEVUBTOTAL, FILEVUBTOTALRES, FILEVUBTOTALNRES;
  TString FILEVUBTRUTHRES, FILEVUBTRUTHNRES;
  TString FILEVCB, FILEDATA, PREFIXOUT, DIRNAME;
  TString FILEVCB1, FILEVCB2;
  TString DSETDIRNAME;
  Double_t MXCUT, MXBIN, CHOPCUT, CHOPBIN, LEPTONPCUT, PRMM2CUT, Q2CUT, NEUCUT; 
  Double_t TRUELEPTONPCUT;
  Double_t MNUSQLOW, MNUSQHIGH,  CHLOW, CHHIGH, DEPL, BTYPE, LEPTTYPE;
  Double_t EMPMLOW, EMPMHIGH;
  Double_t MAXINTPUR, MININTPUR, RUN;
  Double_t VCBCOMP, OTHCOMP;
  Int_t USECB, FIXMEANVALUE, FIXSIGMA, FIXARGUS1, FIXARGUS2;
  Int_t FIXCB1, FIXCB2; 
  Int_t MIXCORR, FITOPT, UNFBINNING, UNFMX2, BRECOQUAL, NUCUT, NEUCORR, UNFFIT, COMP, RECOMPUTEVARS;
  Int_t REL;
  Int_t nDATA, nVCB, nVUB, DOVARSTU;
  std::string CMDLINE;
  Int_t nVCBDATA, nVUBDATA;
  Int_t fprlRew;
  //  Int_t GAUSSFIT, TMODEL;
  Int_t MESFITMODEL;
  Int_t DOBRECOWEIGHT, DOBDECWEIGHT, DODDECWEIGHT, DOFFWEIGHT, DOSSBARWEIGHT;
  Int_t DOMESMEANCORR;
  Int_t MU;
  // fit results
  double vubcomp,errvubcomp,vcbcomp,errvcbcomp,othcomp,errothcomp;
  double vubcompNOMC,errvubcompNOMC,vcbcompNOMC,errvcbcompNOMC;
  double othcompNOMC,errothcompNOMC;
  double vcbothcomp,vubincomp,vuboutcomp;
  double errvcbothcomp,errvubincomp,errvuboutcomp;
  double vubincompNOMC,errvubincompNOMC,vcbothcompNOMC,errvcbothcompNOMC;
  double vuboutcompNOMC,errvuboutcompNOMC;
  // for multplicity categories
  double vcbcompcat[5],errvcbcompcat[5],othcompcat[5],errothcompcat[5];
  
  bool isWriteData;
  //

  CMClass compmod; //CMClass Data Member
  bool IsCM2() const;
  void SetCM2(bool);
  void SetVarfit(bool);
  char* GetEv() const;
  bool GetVarfit() const;

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
  Int_t           Gvcbtyp;
  Int_t           brecoidtrue;
  Int_t           GSem;
  Int_t           GfDpi;
  Int_t           GfDpiz;
  Int_t           GfDk;
  Int_t           GfDks;
  Int_t           GfDkl;
  Int_t           GfDlep;
  Int_t           GfDgam;
  Int_t           GfK;
  Double_t        intpur;
  Int_t           brecoflav;
  Int_t           brecocharge;
  Int_t           brecomc;
  Int_t           modeB;
  Int_t           truemodeB;
  Int_t           brecoqual;
  Double_t        mxhadgen;
  Double_t        mxhadgenwoph;
  Double_t        ctvgen;
  Double_t        ctlgen;
  Double_t        chigen;
  Double_t        pcmsgen;
  Double_t        pcmsgenwph;
  Double_t        tcmsgen;
  Double_t        fcmsgen;
  Double_t        ecmsgen;
  Double_t        pxhadgen;
  Double_t        txhadgen;
  Double_t        fxhadgen;
  Double_t        exhadgen;
  Double_t        fkplus;
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
  Double_t        mxhadneucor;
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
  Int_t           chgdaugen;
  Int_t           npi0;
  Int_t           nneu;
  Int_t           nneu80_160;
  Int_t           nneu160_320;
  Int_t           nneufromB;
  Int_t           nneufromB80_160;
  Int_t           nneufromB160_320;
  Int_t           neudaugen; 
  Int_t           nkp;
  Int_t           nks;
  Int_t           nks_T;
  Int_t           nks_VT;
  Int_t           nks_TVT;
  Double_t        totweight;
  Double_t        totweightNutMult;
  Double_t        totweightTrkMult;
  Double_t        enu;
  Double_t        pnu;
  Double_t        tnu;
  Double_t        fnu;
  Double_t        pplus;
  Double_t        pplusgen;
  Double_t        pplusfit;
  Double_t        pminus;
  Double_t        pminusgen;
  Double_t        pminusfit;
  Double_t        mm2;
  Double_t        mm2neucor;
  Double_t        mm2nc;
  Double_t        mm2fit;
  Double_t        deltam;
  Double_t        wdeltam;
  Double_t        wdeltampiz;
  Double_t        emiss;
  Double_t        pmiss;
  Double_t        emissneucor;
  Double_t        pmissneucor;
  Int_t           neu1B;
  Int_t           ch1B;
  Double_t        KSdecaylenSig[10];
  Double_t        KSmass[10];
  Double_t        KScosp[10];

  //for the different track and neutral lists
  Double_t eUps;
  Double_t pUps;
  Double_t thetaUps;
  Double_t phiUps;
  Double_t eB;
  Double_t pB;
  Double_t thetaB;
  Double_t phiB;
  Double_t esigBcms;
  Double_t psigBcms;
  Double_t thetasigBcms;
  Double_t phisigBcms;
  Double_t elab;
  Int_t numchtrk;
  Double_t trke[30];
  Double_t trkp[30];
  Double_t trktheta[30];
  Double_t trkphi[30];
  Double_t trkecms[30];
  Double_t trkpcms[30];
  Double_t trkthetacms[30];
  Double_t trkphicms[30];
  Int_t trkK[30];
  Int_t trkmu[30];
  Int_t trkel[30];
  Int_t numtotneutrk;
  Double_t neue[200];
  Double_t neup[200];
  Double_t neutheta[200];
  Double_t neuphi[200];
  Double_t neuecms[200];
  Double_t neupcms[200];
  Double_t neuthetacms[200];
  Double_t neuphicms[200];
  Int_t isneumx[200];
  Int_t isneuphloose[200];
  Int_t isneuphdefault[200];
   //the new variables from my Aug08 production and Robertos FF vars
   Double_t Xulabpx;
   Double_t Xulabpy;
   Double_t Xulabpz;
   Double_t XulabpE;
   Double_t Xudaulabpx;
   Double_t Xudaulabpy;
   Double_t Xudaulabpz;
   Double_t XudaulabpE;
   Double_t pxBlab;
   Double_t pyBlab;
   Double_t pzBlab;
   Double_t pEBlab;
   Int_t KSprodratecorr[30];
   Double_t KSp[30];
   Double_t KStheta[30];
   Double_t KSphi[30];
   Double_t KSe[30];
   Int_t KLnum;
   Int_t KLprodratecorr[30];
   Double_t KLp[30];
   Double_t KLtheta[30];
   Double_t KLphi[30];
   Double_t KLe[30];
   Double_t KLptrue[30];
   Double_t KLthetatrue[30];
   Double_t KLphitrue[30];
   Double_t KLetrue[30];
   Int_t tnumKS;
   Double_t tpupsKS[30];
   Int_t tnumKL;
   Double_t tpupsKL[30];
   Double_t dEdxDCHPullpi;
   Double_t dEdxSVTPullpi;
   Double_t dEdxDCHLLRatio;
   Double_t dEdxSVTLLRatio;

  /*
  Double_t        mm1pr;
  Double_t        mm2pr;
  Double_t        mm3pr;
  Double_t        costmiss;
  Double_t        tmiss;
  Double_t        pcmstrklo;
  */
  Double_t        q2;
  Double_t        q2neucor;
  Double_t        q2Gen;
  Double_t        q2nc;
  Double_t        q2fit;
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

  Int_t           nnpi0;

  double vubmc, errvubmc,vubmc2;
  double vuboutmc, errvuboutmc;
  double vcbmc, errvcbmc;
  
  double truehistbins[21];
  std::vector<double> chopBinning;
  double chopBinningUnf[26];

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
  TBranch        *b_Gvxbtyp; 
  TBranch        *b_Gvcbtyp; 
  TBranch        *b_brecoidtrue; 
  TBranch        *b_GSem;    
  TBranch        *b_GfDpi;   
  TBranch        *b_GfDpiz;  
  TBranch        *b_GfDk;    
  TBranch        *b_GfDks;   
  TBranch        *b_GfDkl;   
  TBranch        *b_GfDlep;  
  TBranch        *b_GfDgam;  
  TBranch        *b_GfK;  
  TBranch        *b_intpur;
  TBranch        *b_brecoflav;
  TBranch        *b_brecocharge;
  TBranch        *b_modeB;
  TBranch        *b_truemodeB;
  TBranch        *b_brecomc;
  TBranch        *b_brecoqual;
  TBranch        *b_mxhadgen;
  TBranch        *b_mxhadgenwoph;
  TBranch        *b_ctvgen;
  TBranch        *b_ctlgen;
  TBranch        *b_chigen;
  TBranch        *b_pcmsgen;
  TBranch        *b_pcmsgenwph;
  TBranch        *b_tcmsgen;
  TBranch        *b_fcmsgen;
  TBranch        *b_ecmsgen;
  TBranch        *b_pxhadgen;
  TBranch        *b_txhadgen;
  TBranch        *b_fxhadgen;
  TBranch        *b_exhadgen;
  TBranch        *b_kplus;
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
  TBranch        *b_mxhadneucor;
  TBranch        *b_gmax;
  TBranch        *b_mxhadfit;
  TBranch        *b_lcharge;
  TBranch        *b_pcms;
  TBranch        *b_tcms;
  TBranch        *b_fcms;
  TBranch        *b_ecms;
  TBranch        *b_nle;
  TBranch        *b_nel;
  TBranch        *b_nmu;
  TBranch        *b_nchg;
  TBranch        *b_chgdaugen;
  TBranch        *b_npi0;
  TBranch        *b_nneu;
  TBranch        *b_nneu80_160;
  TBranch        *b_nneu160_320;
  TBranch        *b_nneufromB;
  TBranch        *b_nneufromB80_160;
  TBranch        *b_nneufromB160_320;
  TBranch        *b_neudaugen;
  TBranch        *b_nkp;
  TBranch        *b_nks;
  TBranch        *b_nks_T;
  TBranch        *b_nks_VT;
  TBranch        *b_nks_TVT;
  TBranch        *b_totweight;
  TBranch        *b_totweightNutMult;
  TBranch        *b_totweightTrkMult;
  TBranch        *b_enu;
  TBranch        *b_pnu;
  TBranch        *b_tnu;
  TBranch        *b_fnu;
  TBranch        *b_pplus;
  TBranch        *b_pplusgen;
  TBranch        *b_pplusfit;
  TBranch        *b_pminus;
  TBranch        *b_pminusgen;
  TBranch        *b_pminusfit;
  TBranch        *b_mm2;
  TBranch        *b_mm2neucor;
  TBranch        *b_deltam;
  TBranch        *b_wdeltam;
  TBranch        *b_wdeltampiz;
  TBranch        *b_emiss;
  TBranch        *b_pmiss;
  TBranch        *b_emissneucor;
  TBranch        *b_pmissneucor;
  TBranch        *b_neu1B;   //!
  TBranch        *b_ch1B;   //!

  /*
  TBranch        *b_mm1pr;
  TBranch        *b_mm2pr;
  TBranch        *b_mm3pr;
  TBranch        *b_costmiss;
  TBranch        *b_tmiss;
  TBranch        *b_pcmstrklo;
  */
  TBranch        *b_mm2nc;
  TBranch        *b_mm2fit;
  TBranch        *b_q2;
  TBranch        *b_q2neucor;
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
  //---------------------------
  TBranch        *b_nnpi0;
  TBranch        *b_xLep;
  TBranch        *b_yLep;
  TBranch        *b_zLep;
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
  TBranch        *b_KSdecaylenSig;
  TBranch        *b_KSmass;
  TBranch        *b_KScosp;
  //for the different track and neutral lists
  TBranch *b_eUps;
  TBranch *b_pUps;
  TBranch *b_thetaUps;
  TBranch *b_phiUps;
  TBranch *b_eB;
  TBranch *b_pB;
  TBranch *b_thetaB;
  TBranch *b_phiB;
  TBranch *b_esigBcms;
  TBranch *b_psigBcms;
  TBranch *b_thetasigBcms;
  TBranch *b_phisigBcms;
  TBranch *b_elab;
  TBranch *b_plab;
  TBranch *b_tlab;
  TBranch *b_flab;
  TBranch *b_numchtrk;
  TBranch *b_trke;
  TBranch *b_trkp;
  TBranch *b_trktheta;
  TBranch *b_trkphi;
  TBranch *b_trkecms;
  TBranch *b_trkpcms;
  TBranch *b_trkthetacms;
  TBranch *b_trkphicms;
  TBranch *b_trkK;
  TBranch *b_trkel;
  TBranch *b_trkmu;
  TBranch *b_numtotneutrk;
  TBranch *b_neue;
  TBranch *b_neup;
  TBranch *b_neutheta;
  TBranch *b_neuphi;
  TBranch *b_neuecms;
  TBranch *b_neupcms;
  TBranch *b_neuthetacms;
  TBranch *b_neuphicms;
  TBranch *b_isneumx;
  TBranch *b_isneuphloose;
  TBranch *b_isneuphdefault;
   //the new variables from my Aug08 production and Robertos FF vars
   TBranch *b_Xulabpx;
   TBranch *b_Xulabpy;
   TBranch *b_Xulabpz;
   TBranch *b_XulabpE;
   TBranch *b_Xudaulabpx;
   TBranch *b_Xudaulabpy;
   TBranch *b_Xudaulabpz;
   TBranch *b_XudaulabpE;
   TBranch *b_pxBlab;
   TBranch *b_pyBlab;
   TBranch *b_pzBlab;
   TBranch *b_pEBlab;
   TBranch *b_KSprodratecorr;
   TBranch *b_KSp;
   TBranch *b_KStheta;
   TBranch *b_KSphi;
   TBranch *b_KSe;
   TBranch *b_KLnum;
   TBranch *b_KLprodratecorr;
   TBranch *b_KLp;
   TBranch *b_KLtheta;
   TBranch *b_KLphi;
   TBranch *b_KLe;
   TBranch *b_KLptrue;
   TBranch *b_KLthetatrue;
   TBranch *b_KLphitrue;
   TBranch *b_KLetrue;
   TBranch *b_tnumKS;
   TBranch *b_tpupsKS;
   TBranch *b_tnumKL;
   TBranch *b_tpupsKL;
   TBranch *b_dEdxDCHPullpi;
   TBranch *b_dEdxDVTPullpi;
   TBranch *b_dEdxDCHLLRatio;
   TBranch *b_dEdxSVTLLRatio;

  //End added

  //_______________________
  //Dataset Business


  RooRealVar *Vmes, *Vchop, *Vwe, *VlepYes, *VlepVub, *VlepVcb;
  RooRealVar *VflavB, *VlepYaSe;
  RooRealVar *Vallmes, *Vmxgenwoph, *Vmultcat, *Vmultcatgen, *Vvxbtyp, *Velmom;
  RooRealVar *Vtrumtch;
  RooRealVar *Vksele;
  RooRealVar *Vintpur;
  RooRealVar *Vch,*Vde;
  //CB
  RooRealVar *VmodeB;
  RooRealVar *VtruemodeB;
  RooRealVar *Vhaspi0;
  RooRealVar *Vch1B;
  RooRealVar *Vneu1B;

  RooDataSet *datadata;
  RooDataSet *datamcvub;
  RooDataSet *datamcvcb;
  RooDataSet *datamcoth;
  RooDataSet *unfmcvub;  //RDS for the unfolding histos, KT
  RooDataSet *unftmcvub;  //RDS for the unfolding histos, KT


  bool ischain[9];
 
  RooPlot *xframe;
  RooPlot *xframec;

  double dssR;
  int getMxBin(double mxHad);
  void SetMesFitModel(Int_t);

  //Plots and canvases
  TCanvas* c1; 
  TLatex tl;

  // create pdfs
  RooAbsPdf* createCB(RooRealVar& mes);
  RooAbsPdf* createArgus(RooRealVar& mes, const bool endpointCorrection = false);

  //Mes fit

  mesData* vubMesUnb(RooDataSet *data, RooRealVar *x, std::vector<double>& results, int print, int func, const std::vector<double>& inputPar, const double nSIG, const double nBKG, const char * simply);

  TVector2 sighistounb(RooDataSet *Adata, RooRealVar *Ax, std::vector<double>& results, const std::vector<double>& inputPar, double AnSIG, double AnBKG, const char * Asimply, int fixpar);

  // mes fits
  double dBinomial(double, double); 
  double dstlnuFF(double r1,double r2,double rho2); 
  void FitMes(const char* comp, int mcat);   // 1d fits
  int Bmode(int mode);


  void readmesParam(const TString filename, const int dump = 1);
  void dumpmesParam(const TString filename);

  void mxBinning1d(const TString filename);
  void readUnfParam(TString filename, int dump = 1);
  void dumpUnfParam(const TString filename);
  std::vector<double> mesNsl, mesdatacuts, mesvubcuts, mesvcbcuts, mesothcuts, mesNslMC, mesvcbMC, mesvubMC;
  std::vector<double> mesvubMCchop, mesvubMCall, mesvcbMCall, mesvubMClepteff, mesvubMCalleff;
  double choplowB, chophighB;
  //______________________
  TFile*    openHistFile(TString name);
  TFile*    reopenHistFile(TString name);
  mXClass();
  mXClass(TTree *tree, TString filename, int Sys, int SysD, int unf, double hiunf,int mx2u,int mu,bool iscm2,bool varfit, const std::vector<float>&, int rel, int comp);
  virtual ~mXClass();
  TFile *fHistFile, *fDatasetRootFile;
  TPostScript *fPostScriptFile;
  Int_t  Cut(Int_t entry);
  Int_t  GetEntry(Int_t entry);
  Int_t  LoadTree(Int_t entry);
  void   Init(TTree *tree);
  Bool_t Notify();
  void   InitTruth(TTree *tree);
  Bool_t NotifyTruth();
  
  void   Bookhist();
  virtual Int_t Loop(int isdata, int icat, int isMC, int nres, int truthonly );
  void readCuts(TString filename, int dump, const int runFlag);
  void readOptions(TString filename, int dump);
  void dumpOptions();
  void dumpWFermiFile();
  void dumpCuts(const TString filename);
  void Show(Int_t entry = -1);

  //Functions for Fit parameters and files iitializations
  void initRest(TString filename);
  void InitBinning();
  TString   getfile(const int i);
  TString   getfileVubTotal();
  TString   getfileVubTotalres();
  TString   getfileVubTotalnres();
  TString   getfileVcb();
  TString   getfileVcb1();
  TString   getfileVcb2();
  TString   getfileData();
  TString   getfileVubTruthres();
  TString   getfileVubTruthnres();

  bool*     getfilechain();
  TTree*    getchain(const char* thechain);

  void setTexPrefix(TString s) {texPrefix = s;}
  void setPrefix(TString);
  void setDirectory(TString);
  void setDataDir(TString);

  //Computation of reweighting factors (BR)
  int dImode;
  recoilDSys *Dvar;
  recoilDSys *Bsem;
  double getBsysweight(int decType,int thevub);
  double getFFDstarlnuWeight(const int decType);
  double dstlnuFFSP8(double R1,double R2, double rho2);
  double getDsysweight(int decDpi,int decDk,int decDks,int decDpiz,int decDlep,int decDgam,int decImode,int thevub);

  //Computation of reweighting factors (MC statistic)
  double getGenericSPWeight(int brecoid);
  double getSignalSPWeight(int brecoid);

  //charged multiplicity reweighting
  double getchMultWeight(int n);

  //neutral multiplicity reweighting
  double getneuMultWeight(int n);

  //2d-like charged and neutral multiplicity reweighting
  double getchneuMultWeight(int n);

  //K momentum reweighting
  double getpKWeight(double pk);

  //mm2 reweighting
  double getmm2Weight(double mmsq);

  //Computation of reweighting factors (detector)
  int rHistq2(double amx);
  int rHistel(double bmx);
  int rHistmx(double cmx);
  int newrHistmx(double cmx);
  int newbinrHistmx(double cmx);
  double MatrixW[896];
  double newMatrixW[1024];
  //  std::vector<double> genericMatrixW;
  void readWeights(const int runFlag);
  //  void readpstarfactor();

  //Corrects mx or q2 histos for statistics and mixing
  // double chargeCorr();
  void mixingCorr(int mult);
  void chargeCorr();
  double correctionratiovub;
  double correctionratiovcb;

  //Performs the fit
  void theFit(int multc);
  //Chisquare of Fit
  void compChisq(); 
  double chisq;   
  int NDOF;

  void fitWithErrors(char *typ);
  TGraph* scanParameter(int parnum, int nsig, TMinuit &a, void (*func)(int &, double *, double &, double *, int));

  void doBkgSub(int mult);
  void makeBkgSubPlot();
  void resultDumping();

  void  openEpsFile(TString name, Int_t style=113);
  void  closeEpsFile();
  void  closeHistFile();
  void  shrinkPad(double b, double l, double r, double t);

  //Produces the unfolding histos
  void UnfHistos(int mult);
  //Correction factors for multiplicity categories
  void MultCorr();

//   //Efficiency busines
  double vubmcaftercuts;
//   double epsu, errepsu, epschop, errepschop, epstot, errepstot;
//   double epsu_mixcorr, errepsu_mixcorr, summ, errsumm;
//   double Pepstot, errPepstot;
//   double epsphsp, errepsphsp;
//   double vubmcallforeff, vubmcleptforeff, vubmcallforeff_mixcorr, vubmcaftercuts;
//   double tot, errtot, totmc, errtotmc;
//   double vcbmcselected, vubmcselected, othermc;
//   double vubmcnocut,vcbmcorig, vubmcorig, othermcorig;
//   double vubmccut, vcbmcnocut, othmcnocut;
//   double avubmcSB, aerrvubmcSB; 
//   double vubmcSB, errvubmcSB; 

  //binning
  int nQ2B, nMxB, nB;
  std::vector<double> mxB1, chopB1;

//   //Old Reweighting
//   double TrueMxWeight[42];
//   double getTrueMxWeight(double thetrumx, int index);
//   int    TrueHist(double mxt);

};

//ReadwFermiFile
void ReadwFermiFile(const TString &, std::vector<float>&,int);
//Function for truth-matching
 void GetMultiplicity(int,int&,int&,int&,int&);
#endif
