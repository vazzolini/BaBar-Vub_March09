#ifndef VirClass_h
#define VirClass_h

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
#include <TMatrixD.h>
#include <TRandom3.h>

#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooFormulaVar.hh"

#include "RecoilAnalysis/recoilDSys.hh"
#include "RecoilAnalysis/mesData.hh"
#include "RecoilAnalysis/recoilAnalysis.hh"

#include "VirVubFitter/XSLEvtFFWeight.hh"
#include "VirVubFitter/XSLBToDstrlnu_DstrToDpi_CLN.hh"
#include "VirVubFitter/XSLBall04_pilnu.hh"
#include "VirVubFitter/XSLBall04_etalnu.hh"
#include "VirVubFitter/XSLBall05.hh"

#include "VirVubFitter/CMClass.hh"

class TF1;
class TH1;
//class TH1D;
//class TH2D;
class TGraph;
class TFile;
class TPostScript;
class TCanvas;
class TPad;

class RooDataSet;
class RooPlot;
class RooAbsPdf;

class TMesCor;
class CFDump;
class XSLEvtFFWeight;
class XSLBToDstrlnu_DstrToDpi_CLN;
// ################################################

class VirClass {

public :
  enum fileTypeEnum  { VubTotal=0, VubTotalres=1, VubTotalnres=2, Vcb=3, Vcb1=4, Vcb2=5, Data=6, VubTruthres=7, VubTruthnres=8 };

  enum var1DEnum     { iVarMx=0, iVarQ2=1, iVarPplus=2 };
  enum shapeFuncEnum { eDefault=0, eBelle04=1, eBabar05=2, eHFAGComb06=3, summer07=4 };

  enum parFitStdEnum { iMean=0, iSigma=1 , iAlpha=2, iN=3, iArgus=4 };
  enum parFitThoEnum { iCutOff=5, iEndpoint, iThoSigR, iSigma_r1, iThoSigXc, iSigma_r2, iSigma_l, iThoSigN, iThoSigAlpha };

  enum yieldsEnum    { iNsig = 0, iNbkg, iNbpk };
  
  enum { GaussFit = 0, ArgusAndCB = 1, ArgusAndThosig = 2, ThreePDFs = 3};
  

  TTree          *fChain, *fOutTree;   //pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //current Tree number in a TChain
  Int_t   fnewbin, mesIdx;
  TString fWeightFile; 
  TString fOptionFile; 
  TString texPrefix;
  TString FILEVUBTOTAL, FILEVUBTOTALRES, FILEVUBTOTALNRES;
  TString FILEVUBTRUTHRES, FILEVUBTRUTHNRES;
  TString FILEVCB, FILEDATA, PREFIXOUT, DIRNAME;
  TString DSETDIRNAME;
  TString FILEVCB1, FILEVCB2;
  Double_t TOTALSTAT, TOTALSTATMODEL, BRRATIOGENVALUE,  BRRATIOVALUETAIL_U, BRRATIOVALUETAIL_C, LUMI_GENERIC, LUMI_SIGNAL, LUMI_DATA;
  Double_t BRRATIOVALUETAIL_UCHAR, BRRATIOVALUETAIL_CCHAR, BRRATIOVALUETAIL_UNEUT, BRRATIOVALUETAIL_CNEUT;
  Double_t COSTMISSCUTLO, COSTMISSCUTHI, PCMSTRKLOCUT, PMISSCUTLO;
  Double_t MXCUT, MXBIN, PPLUSBIN, Q2CUT, Q2BIN, PPLUSCUT, CHOPCUT, CHOPBIN, PSTARFACT, LEPTONPCUT, PRMM1CUT, PRMM2CUT, PRMM3CUT;
  Double_t TRUELEPTONPCUT;
  Double_t CSILOCUT, CSIHICUT, Q2LOCUT, Q2HICUT, XLOCUT, XHICUT, WLOCUT;
  Double_t WHICUT, EWPWLOCUT, EWPWHICUT;
  Double_t MNUSQLOW, MNUSQHIGH,  CHLOW, CHHIGH, DEPL, BTYPE;
  Double_t EMPMLOW, EMPMHIGH;
  Double_t MAXINTPUR, MININTPUR, RUN, TOYMC;
  Double_t DELTAMB, DELTAA, FERMIAPP, MULFAC;
  Double_t VCBCOMP, OTHCOMP;
  Double_t SMEARALLMEANVALUE, SMEARALLSIGMA, SMEARBKGMEANVALUE, SMEARBKGSIGMA;
  Double_t magic_k_factor_B0, magic_k_factor_BCh; 
  Int_t SEME,SERVICE;
  Bool_t THECOMPARISON, SAVEPDFTREE, READPDFTREE, TOYHISTOGRAMES;
  Double_t THECASCADEWEIGHT;
  Double_t RANDOMSEED, BLINDSIZE;
  Int_t USECB, FIXMEANVALUE, FIXSIGMA, FIXARGUS1, FIXARGUS2;
  Int_t FIXCB1, FIXCB2, FITTOTSHAPE, MIXCORR, FITMC, FITOPT, FITQ2, Sun, UNFBINNING, UNFMX2, BRECOQUAL, SMALLSTATCORR, NUCUT;
  Int_t REL, SUBTRACTPEAKING, FITALLMESRANGE,  LEPTTYPE;
  Int_t DONTWEIGHTVUB;
  Int_t MULTIFIT, BLINDING;
  Float_t B0B0bar_gene,B0B0bar_skimmed,BplusBminus_gene,BplusBminus_skimmed,Vub_excl_gene,Vub_excl_skimmed,Vub_incl_gene,Vub_incl_skimmed;
  //  int_t NORE, RE;
  Int_t nDATA, nVCB, nVUB, DOVARSTU;
  std::string CMDLINE;
  Int_t nVCBDATA, nVUBDATA;
  Int_t fprlRew;
  Int_t BINNED, NOTUNBINNED;
  //  Int_t GAUSSFIT, TMODEL;
  Int_t MESFITMODEL;
  bool  EPOINTCOR,  FIXCORRRATIO;
  Int_t FITDSS;
  Double_t CUTBIN, EFFDFN;
  Double_t TOTDPESO;
  Int_t TOTD;
  Int_t ISSMEARALL, ISSMEARBKG, CUTNNPI0;
  Int_t DOBRECOWEIGHT, DOBDECWEIGHT, DODDECWEIGHT, DOTRKWEIGHT,DOFFWEIGHT, DOSSBARWEIGHT, DOMESMEANCORR, DOCASCADEWEIGHT, DOEXCLFFWEIGHT;
  Int_t ME,MU;
  Int_t DONEUWEIGHT, DOFERMI;
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

  TRandom * rnd;
  TRandom3 *randz;
  TMesCor* mesCor; // class for handling mes corrections on data

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
  Int_t           GfK;
  Double_t        intpur;
  Int_t           brecoflav;
  Int_t           brecocharge;
  Int_t           brecoid;
  Int_t           brecoidtrue;
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
  Double_t        mm2nc;
  Double_t        mm2fit;
  Double_t        deltam;
  Double_t        wdeltam;
  Double_t        emiss;
  Double_t        pmiss;
  Int_t           neu1B;
  Int_t           ch1B;
  Int_t           cascade;
  Double_t        wdeltampiz;
  Int_t           nkpL;
  Int_t	  	  nkpT;
  Int_t		  nkpVT;
  Int_t		  nkpNNL;
  Int_t		  nkpNNT;
  Int_t		  nkpNNVT;
  Double_t	   pchi2xlbs;
  Double_t        Xulabpx;
  Double_t        Xulabpy;
  Double_t        Xulabpz;
  Double_t        XulabpE;
  Double_t        Xudaulabpx;
  Double_t        Xudaulabpy;
  Double_t        Xudaulabpz;
  Double_t        XudaulabpE;
  Double_t        pxBlab;
  Double_t        pyBlab;
  Double_t        pzBlab;
  Double_t        pEBlab;
  Double_t        pxleptgen;;
  Double_t        pyleptgen;;
  Double_t        pzleptgen;;
  Double_t        elabgen;;



  /*
  Double_t        mm1pr;
  Double_t        mm2pr;
  Double_t        mm3pr;
  Double_t        costmiss;
  Double_t        tmiss;
  Double_t        pcmstrklo;
  */
  Double_t        q2;
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
  TF1* f1all;
  TF1* f1bkg;
  
  double TrackingWeight[20];
  double NeutralWeight[20];
  double BrecoWeight[20];
  double truehistbins[21];
  std::vector<double> chopBinning;
  std::vector<float> sigtopeakratiocorrmx, sigtopeakratioerrcorrmx;
  std::vector<float> sigtopeakratiovub, sigtopeakratioerrvub;
  std::vector<float> sigtopeakratiovcb, sigtopeakratioerrvcb;
  std::vector<float> sigtopeakratiooth, sigtopeakratioerroth;
  std::vector<float> sigtopeakratiovubin, sigtopeakratioerrvubin;
  std::vector<float> sigtopeakratiovubout, sigtopeakratioerrvubout;
  std::vector<float> sigtopeakratiovcboth, sigtopeakratioerrvcboth;
  double chopBinningUnf[26];
  //  double mxBinning[10];
  //  double q2Binning[13];
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
  TBranch        *b_Gvxbtyp; 
  TBranch        *b_GSem;    
  TBranch        *b_GfDpi;   
  TBranch        *b_GfDpiz;  
  TBranch        *b_GfDk;    
  TBranch        *b_GfDks;   
  TBranch        *b_GfDlep;  
  TBranch        *b_GfDgam;  
  TBranch        *b_GfK;  
  TBranch        *b_intpur;
  TBranch        *b_brecoflav;
  TBranch        *b_brecocharge;
  TBranch        *b_brecoid;
  TBranch        *b_brecoidtrue;
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
  TBranch        *b_deltam;
  TBranch        *b_wdeltam;
  TBranch        *b_emiss;
  TBranch        *b_pmiss;
  TBranch        *b_neu1B;   //!
  TBranch        *b_ch1B;   //!
  TBranch        *b_cascade;   //!
  TBranch        *b_wdeltampiz;
  TBranch        *b_nkpL;
  TBranch        *b_nkpT;
  TBranch        *b_nkpVT;
  TBranch        *b_nkpNNL;
  TBranch        *b_nkpNNT;
  TBranch        *b_nkpNNVT;
  TBranch        *b_pchi2xlbs;
  TBranch        *b_Xulabpx;
  TBranch        *b_Xulabpy;
  TBranch        *b_Xulabpz;
  TBranch        *b_XulabpE;
  TBranch        *b_Xudaulabpx;
  TBranch        *b_Xudaulabpy;
  TBranch        *b_Xudaulabpz;
  TBranch        *b_XudaulabpE;
  TBranch        *b_pxBlab;
  TBranch        *b_pyBlab;
  TBranch        *b_pzBlab;
  TBranch        *b_pEBlab;
  TBranch        *b_pxleptgen;
  TBranch        *b_pyleptgen;
  TBranch        *b_pzleptgen;
  TBranch        *b_elabgen;
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

  //End added

  //_______________________
  //Dataset Business


  RooRealVar *Vmes, *Vchop, *Vwe, *VlepYes, *VlepVub, *VlepVcb, *VlepVubSB;
  RooRealVar *VflavB, *VlepYaSe;
  RooRealVar *Vallmes, *Vmxgenwoph, *Vmultcat, *Vmultcatgen, *Vvxbtyp, *Velmom;
  RooRealVar *Vtrumtch;
  RooRealVar *Vksele;
  RooRealVar *Vpplus;
  RooRealVar *Vintpur;
  RooRealVar *Vch,*Vde;
  //CB
  RooRealVar *VmodeB;
  RooRealVar *VtruemodeB;
  RooRealVar *Vhaspi0;
  RooRealVar *Vch1B;
  RooRealVar *Vneu1B;
  RooRealVar *Vbrecoid;
  RooRealVar *Vbrecoidtrue;
  // --> for Thecomparison
  RooRealVar *Vpcms, *Vmm2, *Vnchg, *Vnkp, *Vnks, *Vnneu, *Vqtot, *Vwdeltam, *Vwdeltampiz, *Vemiss, *Vpmiss;
  RooRealVar *VlepYesBitmask, *VAllCutBitmask;

  RooDataSet *pstarsample;
  RooDataSet *pstarsamplesum;
  RooDataSet *datadata;
  RooDataSet *datamcvub;
  RooDataSet *datamcvcb;
  RooDataSet *datamcoth;
  RooDataSet *datavcboth;
  RooDataSet *datavubin;
  RooDataSet *datavubout;
  RooDataSet *unfmcvub;  //RDS for the unfolding histos, KT
  RooDataSet *unftmcvub;  //RDS for the unfolding histos, KT

  //Bidimensional fit 
  RooRealVar *Vmx, *Vq2;

  bool ischain[9];
 
  RooPlot *xframe;
  RooPlot *xframec;

  bool countMC;
  double dssR;
  void setCountingFlag(bool apply);
  void setDssRatio(double dssRatio);
  void setDssRatio(const TString dssFile);
  int getMxBin(double mxHad);
  void SetMesFitModel(Int_t);

  //Plots and canvases
  TCanvas* c1; 
  TLatex tl;

  // create pdfs
  RooAbsPdf* createThorsten(RooRealVar& mes);
  RooAbsPdf* createCCB(RooRealVar& mes, const bool endpointCorrection = false, RooRealVar* Rendpoint = 0);
  RooAbsPdf* createCB(RooRealVar& mes);
  RooAbsPdf* createArgus(RooRealVar& mes, const bool endpointCorrection = false);

  //Mes fit
  mesData* vubMesUnb(RooDataSet *data, RooRealVar *x, std::vector<double>& results, int print, int func, const std::vector<double>& inputPar, const double nSIG, const double nBKG, const char * simply,bool fixsigtopeakratio=false, int jmatch = 0);

  TVector2 sighistounb(RooDataSet *Adata, RooRealVar *Ax, std::vector<double>& results, const std::vector<double>& inputPar, double AnSIG, double AnBKG, const char * Asimply, int fixpar, bool fixsigpeakratio=false, int jmatch = 0);

  //plots for BAD
  void FitPlots(const char* comp, int su, int cut);

  // mes fits
  double dBinomial(double, double); 
  double dstlnuFF(double r1,double r2,double rho2); 
  double dstlnuFFSP8(double,double,double);
  void FitMes(const char* comp, int mcat, int su,bool fixsigpeakratio=false);   // 1d fits
  void FitMes2D(const char* comp, int su, bool fixsigpeakratio=false);          // 2d fits
  //CB
  void getSOverPk(const char* comp, int mcat, int su); 
  void setActiveCorrection(const char* comp, int signalUnfolding, int iBin);
  float SPrandomized(float mean,float lsig, float hsig, TRandom& rndm, int seed);
  int Bmode(int mode);


  void readmesParam(const TString filename, const int dump = 1);
  void dumpmesParam(const TString filename);

  void q2Binning(const TString filename);
  void mxBinning(const TString filename);
  void mxBinning1d(const TString filename);
  void pplusBinning1d(const TString filename);
  void readUnfParam(TString filename, int dump = 1);
  void dumpUnfParam(const TString filename);
  std::vector<double> mesNsl, mesdatacuts, mesvubcuts, mesvcbcuts, mesothcuts, mesNslMC, mesvcbMC, mesvubMC, mesvuboutMC, mespstarMC, mespstarcuts;
  std::vector<double> mesvubMCchop, mesvubMCall, mesvcbMCall, mesvubMClepteff, mesvubMCalleff;
  double choplowB, chophighB;
  //______________________
  TFile*    openHistFile(TString name);
  TFile*    reopenHistFile(TString name);
  void   openPdfToyFile(const TString&);
  
  VirClass();
  VirClass(TTree *tree, TString filename, int Sys, int q2, int comb, int un, int unf, double hiunf,int mx2u,int me,int mu,bool iscm2,bool varfit, const std::vector<float>&, int newbin, int SPseed, bool isWD, int rel, int bsys, bool thecmp);
  virtual ~VirClass();
  TFile *fHistFile, *fDatasetRootFile, *fPdfToyFile;
  TPostScript *fPostScriptFile;
  Int_t  Cut(Int_t entry);
  Int_t  GetEntry(Int_t entry);
  Int_t  LoadTree(Int_t entry);
  void   Init(TTree *tree);
  Bool_t Notify();
  void   InitTruth(TTree *tree);
  Bool_t NotifyTruth();
  
  void   Bookhist();
  virtual Int_t Loop(int isdata, int icat, int isMC, int nres, int comb, int truthonly );
  void readCuts(TString filename, int dump, double wsys, const int runFlag, const int Rel);
  void readOptions(TString filename, int dump);
  void dumpOptions();
  void dumpWFermiFile();
  void dumpCuts(const TString filename);
  void Show(Int_t entry = -1);

  //Functions for Fit parameters and files iitializations
  void initRest(TString filename);
  void InitBinning(const int comb);
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
  Int_t     isfitMC();
  Int_t     dontWeightVub();

  void applyEndpointCor(const bool apply, const int runFlag);
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
  double getFFXulnuWeight(const int decType);
  double getDsysweight(int decDpi,int decDk,int decDks,int decDpiz,int decDlep,int decImode,int thevub);

  vector<TObject*> grfmngr;

  //Computation of reweighting factors (MC statistic)
  double getGenericSPWeight();
  double getSignalSPWeight();

  //Computation of reweighting factors (detector)
  double getTrackingWeight();
  double getNeutralWeight();
  double getBrecoWeight(double theintpur);
  double getCascadeDecWeight();
  double FermiWeight(double kp, double deltamb, double deltaa);
  int genericHistmx(double cmx);
  int rHistq2(double amx);
  int rHistel(double bmx);
  int rHistmx(double cmx);
  int newrHistmx(double cmx);
  int newbinrHistmx(double cmx);
  double MatrixW[896];
  double newMatrixW[1024];
  //  Double_t genericWeightsB[1024]; 
  std::vector<double> genericMatrixW;
  double fermi(double kp, double m, double a);
  void readWeights(double wsys, const int runFlag);
  void readpstarfactor();

  //Corrects mx or q2 histos for statistics and mixing
  // double chargeCorr();
  void mixingCorr(int comb,int mult, int su);
  void chargeCorr(int su,int comb); //deprecated
  void computeNSLandEfficiencies(int , int);
  double correctionratiovub,correctionratiovubout,correctionratiovcb;
  
  //Computes quantities for BR extraction
  double getPstarFactor(double); // returns the precomputed pstarfactor
  void calcPstarFact(const int comb, const int su, const int ck, const int runperiod); //Gets or computes the pstar factor
  TVector2 unfolPstarFact(Double_t &mx, Double_t &q2);
  
  //Performs the fit
  void theFit(int cmb,int multc,int su, int ck);
  //Performs only the final fit (takes mx-q2 shapes from a root file)
  void theFitOnly(int cmb,int multc,int su, int ck);
  //Chisquare of Fit
  void compChisq(int cmb, int su); double chisq;   int NDOF;

  void fitWithErrors(int cmb,int su,char * typ);
  TGraph* scanParameter(int parnum, int nsig, TMinuit &a, void (*func)(int &, double *, double &, double *, int));

  //Computes, plots & writes down results
  void ruslExtr(int cmb, int su, int ck);
  TVector2 getEffFromBauer(const Double_t &mx, const Double_t &q2);
  TVector2 getEffFromDFN2D(const Double_t &mx, const Double_t &q2);
  TVector2 getEffFromDFN(const Double_t var, const int iType); // wrapper function
  TVector2 getEffFromDFN_Q2(const Double_t Q2);                // for Q2
  TVector2 getEffFromDFN_Pp(const Double_t Pp);                // for Pplus
  double getblindfact();
  //TVector2 getEffFromDFN(Double_t &mx);
  
  void doBkgSub(int mult, int su);
  void makeBkgSubPlot(int su);
  void resultDumping(int cmb, int su);

  void  openEpsFile(TString name, Int_t style=113);
  void  closeEpsFile();
  void  closeHistFile();
  //  void  closePdfToyFile();
  void  shrinkPad(double b, double l, double r, double t);

  //Produces the unfolding histos
  void UnfHistos(int mult);
  //Correction factors for multiplicity categories
  void MultCorr();

  //Global variables were to store fit output
  double blindfactor;
  double vubSB, errvubSB, MCerrvubSB;
  double vcbSB, errvcbSB, errfitvcbSB;
  double othSB, errothSB, errfitothSB;
  double dataFirstBin, dataErrFirstBin, othErrFirstBin, othFirstBin, vcbFirstBin, vcbErrFirstBin;
  double vcbothSB, vuboutSB, evcbothSB2, evuboutSB2, efitvcbothSB2, efitvuboutSB2,vuboutFirstBin, vuboutErFirBin ;

  double areavcbSB, areaothSB, areavubmcSB;
  double errvcbSBtemp, errvcbSB2, errfitvcbSBtemp,errfitvcbSB2;
  double errothSBtemp, errothSB2, errfitothSBtemp,errfitothSB2;
  double errvubSBtemp, errvubSB2, errvubSBNOFITtemp, errvubSBNOFIT2, errvubmcSBtemp,errvubmcSB2,aerrvubmcSB2;
  double dataErrFirstBintemp,dataErrFirstBin2, vcbErrFirstBintemp, vcbErrFirstBin2, othErrFirstBintemp, othErrFirstBin2;
 
  double BRBR, errBRBR, errBRBRMCstat, S, errS;
  double PBRBR, errPBRBR, errPBRBRMCstat;
  //Efficiency business
  double epsu, errepsu, epschop, errepschop, epstot, errepstot;
  double epsu_mixcorr, errepsu_mixcorr, summ, errsumm;
  double Pepstot, errPepstot;
  double epsphsp, errepsphsp;
  double vubmcallforeff, vubmcleptforeff, vubmcallforeff_mixcorr, vubmcaftercuts;
  double tot, errtot, totmc, errtotmc;
  double vcbmcselected, vubmcselected, othermc;
  double vubmcnocut,vcbmcorig, vubmcorig, othermcorig;
  double vubmccut, vcbmcnocut, othmcnocut;
  double avubmcSB, aerrvubmcSB; 
  double vubmcSB, errvubmcSB; 

  //Correction factors
  TVector2 ApplyCorrFactors(const TVector2 & fl3,const TVector2 & fl4, const TVector2& fl5, bool hasmixingcorr, Float_t chargecorr=1);
  TVector2 ApplyCorrFactors(const TVector2 & fl3,const TVector2 & fl4, const TVector2& fl5, bool hasmixingcorr, const TString&, const TString&, Float_t chargecorr=1 );
  double calcpstarfact, errcalcpstarfact, calcpstarfacttemp;
  double fact;

  //binning
  int nQ2B, nMxB, nB;
  std::vector<double> qB1, xB1, mxB1, pplusB1, chopB1;
  std::vector<double> dssRew;

  //Old Reweighting
  double TrueMxWeight[42];
  double getTrueMxWeight(double thetrumx, int index);
  int    TrueHist(double mxt);

  //2D analysis
  //  void doBkgSub2D();
  void doBkgSub2D(int su);
  void makeBkgSubPlot2D(int su);
  void computeChargeCorr(int su, int comb);

  //Debugging
  void   Debug();

//read datasets
  int readDataFile(RooDataSet**,TString);

  //read Correction ratio for mx analysis
  void readmxcorrratiosigpeak(TString);
  float* activecorrection, *activecorrectionerror;

  //Peaking background correction factors;
  void evaluatePeakingBackground(int,bool);
  //  void evaluatePeakingBackgroundMX();
  vector <TVector2> pkgbkgcorr_AC, pkgbkgcorr_SL;

  //write datasets
  void writeDataFile(const std::vector<RooDataSet*>&);

  //Dataset file handlers
  void openDataSetFile(TString);
  void closeDataSetFile();

  //MC BRratios
  void SetBRValues(Int_t,Int_t);
  void GetGeneratedEvents(Int_t,Int_t);
  void GetMCTrueNumbers(int);
  void CreateNewPstarSample(bool);
  //void SetGenericHybWeights();
  ofstream debugout;

  //mes Fit Parameters issue
  void FillPDFTree(const RooAbsPdf*, Int_t, Double_t,Int_t);
  void SavePDFTree();

  RooDataSet* RetrieveDataset(Int_t);

  //Crossfeed utils
  void CalculateCFMatrix(int,int);
  vector<CFDump*> cfdump;
  void ApplyCrossFeedCalculations(const TString&, const TString&, const TVector2&, const TVector2&, TVector2&);
  void ApplyCrossFeedCalculations(const TString&, const TString&, TH1*, TH1*);
  //void ApplyCrossFeedCalculations(const TString&, const TString&, TH2D*, TH2D*);
  //DEBUGGGONE
  Float_t mySL,mylpYesSig, mynucut, mychcut, mydepl, myWdeltaM, myAllcut;

private:

  Int_t _debug;

  Int_t _shapeFunction;
  
};

Double_t SetFlav(Int_t);
Double_t CFerrorProp(bool, bool, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t);
void getExclusiveBFRatios(Float_t&,Float_t&);
//ReadwFermiFile
void ReadwFermiFile(const TString &, std::vector<float>&,int);
//q2 background subtr
void q2bkgsub(const TH1D*, const TH1D*, const TH1D*,const int, TH1D*);
//Function for truth-matching
void GetMultiplicity(int,int&,int&,int&,int&);

class CFUtil
{
public:
  CFUtil();
  CFUtil(RooDataSet*, vector<double>&);
  ~CFUtil();
  void SetMesInitPar(vector<double>&);
  RooDataSet * GetDataset() const;
  vector<double> GetMesParams() const;
  
private:
  RooDataSet* rds;
  vector<double>* mesinitpar;
  
};

class CFDump 
{
public:
  CFDump();
  CFDump(const TString&, const TString&, const TMatrixD&, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t);
  ~CFDump();
  
public:
  TString nm;
  TString ct;
  TMatrixD probinv;
  Double_t gnn;
  Double_t sigmagnn;
  Double_t gn;
  Double_t sigmagn;
  Double_t gcc;
  Double_t sigmagcc;
  Double_t gc;
  Double_t sigmagc;
  Double_t meastotneu;
  Double_t meastotcha;
};



#endif
