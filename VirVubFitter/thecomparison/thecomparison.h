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
#include "RecoilAnalysis/recoilAnalysis.hh" 
#include "RecoilAnalysis/recoilDSys.hh" 
//#include "../../RecoilAnalysis/mesData.hh"  //included by recoilDSys.hh

#include "RooFitCore/RooAbsPdf.hh"
#include "RooFitCore/RooDataHist.hh"
#include "RooFitCore/RooRealVar.hh"

#include "RooCCB.hh"
#include "RooThorstenSig.hh"
#include "RooFitModels/RooArgusBG.hh"

#include <vector>

using namespace std;
using namespace RooFit;

class thecomparison {
   public :

     enum parFitStdEnum { iMean=0, iSigma=1 , iAlpha=2, iN=3, iArgus=4, 
			  iCutOff=5, iEndpoint, iThoSigR, iSigma_r1, iThoSigXc, iSigma_r2, iSigma_l, iThoSigN, iThoSigAlpha };

   TTree          *fChain;   //pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //current Tree number in a TChain
   recoilDSys *Dvar;
   recoilDSys *Bsem;

    // Declaration of leave types
   Int_t           run;
   Int_t           upper;
   Int_t           lower;
   Int_t           Gvxbtyp;
   Int_t           GSem;
   Int_t           GfDpi;
   Int_t           GfDpiz;
   Int_t           GfDk;
   Int_t           GfDks;
   Int_t           GfDkl;
   Int_t           GfDlep;
   Int_t           GfDgam;
   Int_t           GfDnu;
   Int_t           GfD0Ds;
   Int_t           GfDDs;
   Int_t           GfDkspiopio;
   Int_t           GfK;
   Int_t           isassocB;
   Int_t           isassocB_GHIT;
   Float_t         ass_deltapB;
   Int_t           isGoodMatch;
   Int_t           ch1B;
   Int_t           ch2B;
   Int_t           chunm;
   Int_t           neu1B;
   Int_t           neu2B;
   Int_t           neuunm;
   Int_t           brecoqual;
   Int_t           brqual;
   Float_t         brecoqualangle;
   Int_t           chgdaugen;
   Int_t           neudaugen;
   Int_t           nchg;
   Int_t           nneu;
   Int_t           xcharge;
   Int_t           nB;
   Int_t           brecoid;
   Int_t           brecoidtrue;
   Int_t           brecoflav;
   Int_t           brecocharge;
   Int_t           modeB;
   Int_t           truemodeB;
   Int_t           isdoubleD;
   Double_t        mes;
   Double_t        mesendpoint;
   Double_t        de;
   Double_t        pB;
   Double_t        eB;
   Double_t        eUps;
   Double_t        pUps;
   Double_t        thetaUps;
   Double_t        phiUps;
   Double_t        thetaB;
   Double_t        phiB;
   Double_t        pBtrue;
   Double_t        eBtrue;
   Double_t        tBtrue;
   Double_t        fBtrue;
   Double_t        pur;
   Double_t        intpur;
   Int_t           nle_nopcut;
   Int_t           nle;
   Int_t           nel;
   Int_t           nmu;
   Int_t           nkp;
   Int_t           nks;
   Int_t           nlept500;
   Int_t           nelec500;
   Int_t           nmu500;
   Int_t           nlept1000;
   Int_t           nelec1000;
   Int_t           nmu1000;
   Double_t        deltam;
   Double_t        MM1pr;
   Double_t        MM2pr;
   Double_t        MM3pr;
   Double_t        OA1;
   Double_t        OA2;
   Double_t        OA3;
   Double_t        PiMin1;
   Double_t        PiMin2;
   Double_t        PiMin3;
   Double_t        plab;
   Double_t        elab;
   Double_t        tlab;
   Double_t        flab;
   Double_t        plabgen;
   Double_t        elabgen;
   Double_t        tlabgen;
   Double_t        flabgen;
   Int_t           lchargegen;
   Double_t        pcms;
   Double_t        ecms;
   Double_t        tcms;
   Double_t        fcms;
   Int_t           lcharge;
   Int_t           leptidgen;
   Int_t           leptorg;
   Int_t           isele;
   Int_t           vub;
   Int_t           vcb;
   Int_t           other;
   Int_t           nvubexcl;
   Int_t           nvubnres;
   Int_t           ntau;
   Double_t        mxhadgen;
   Double_t        mxhadgenwoph;
   Int_t           xchargegen;
   Double_t        pcmsgen;
   Double_t        ecmsgen;
   Double_t        tcmsgen;
   Double_t        fcmsgen;
   Double_t        pxhadgen;
   Double_t        exhadgen;
   Double_t        exhadgencms;
   Double_t        txhadgen;
   Double_t        fxhadgen;
   Double_t        q2Gen;
   Double_t        ctvgen;
   Double_t        ctlgen;
   Double_t        chigen;
   Double_t        enugencms;
   Double_t        pnugencms;
   Double_t        mxhad;
   Double_t        emiss;
   Double_t        pmiss;
   Double_t        mm2;
   Double_t        q2;
   Double_t        q2new;
   Double_t        exhadcms;
   Double_t        pxhad;
   Double_t        exhad;
   Double_t        txhad;
   Double_t        fxhad;
   Double_t        pnu;
   Double_t        tnu;
   Double_t        fnu;
   Double_t        pplus;
   Double_t        pminus;
   Double_t        pplusgen;
   Double_t        pminusgen;
   Double_t        pplusfit;
   Double_t        pminusfit;
   Double_t        wdeltam;
   Double_t        mxhadfit;
   Double_t        mm2fit;
   Double_t        q2fit;
   Double_t        chisq;
   Double_t        globchisq;
   Double_t        probchisq;
   Int_t           ndof;
   Int_t           fitstatus;
   Int_t           nnpi0;
   Int_t           nneu80_160;
   Int_t           nneu160_320;
   Double_t        pstarfitlept;
   Double_t        pfitX;
   Double_t        thetafitX;
   Double_t        phifitX;
   Double_t        pfitlept;
   Double_t        thetafitlept;
   Double_t        phifitlept;
   Double_t        pfitB;
   Double_t        thetafitB;
   Double_t        phifitB;
   Double_t        totweight;
   Double_t        totweightNutMult;
   Double_t        totweightTrkMult;
   Double_t        kplus;
   Int_t           nBrems;
   Float_t         eBrems[3];   //[nBrems]
   Double_t        mxhadfit0;
   Double_t        q2fit0;
   Double_t        pplusfit0;
   Double_t        chisqfit0;
   Int_t           fitstatusfit0;
   Double_t        chisqfit;
   Int_t           fitstatusfit;
   Double_t        mxhadfit1;
   Double_t        q2fit1;
   Double_t        pplusfit1;
   Double_t        chisqfit1;
   Int_t           fitstatusfit1;
   Double_t        mxhadfit2;
   Double_t        q2fit2;
   Double_t        pplusfit2;
   Double_t        chisqfit2;
   Int_t           fitstatusfit2;
   Double_t        mxhadfit3;
   Double_t        q2fit3;
   Double_t        pplusfit3;
   Double_t        chisqfit3;
   Int_t           fitstatusfit3;
   Double_t        mxhadfit4;
   Double_t        q2fit4;
   Double_t        pplusfit4;
   Double_t        chisqfit4;
   Int_t           fitstatusfit4;
   Double_t        mxhadfit5;
   Double_t        q2fit5;
   Double_t        pplusfit5;
   Double_t        chisqfit5;
   Int_t           fitstatusfit5;
   Double_t        mxhadfit6;
   Double_t        q2fit6;
   Double_t        pplusfit6;
   Double_t        chisqfit6;
   Int_t           fitstatusfit6;
   Double_t        mxhadfit7;
   Double_t        q2fit7;
   Double_t        pplusfit7;
   Double_t        chisqfit7;
   Int_t           fitstatusfit7;
   Double_t        mxhadfit8;
   Double_t        q2fit8;
   Double_t        pplusfit8;
   Double_t        chisqfit8;
   Int_t           fitstatusfit8;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_upper;   //!
   TBranch        *b_lower;   //!
   TBranch        *b_Gvxbtyp;   //!
   TBranch        *b_GSem;   //!
   TBranch        *b_GfDpi;   //!
   TBranch        *b_GfDpiz;   //!
   TBranch        *b_GfDk;   //!
   TBranch        *b_GfDks;   //!
   TBranch        *b_GfDkl;   //!
   TBranch        *b_GfDlep;   //!
   TBranch        *b_GfDgam;   //!
   TBranch        *b_GfDnu;   //!
   TBranch        *b_GfD0Ds;   //!
   TBranch        *b_GfDDs;   //!
   TBranch        *b_GfDkspiopio;   //!
   TBranch        *b_GfK;   //!
   TBranch        *b_isassocB;   //!
   TBranch        *b_isassocB_GHIT;   //!
   TBranch        *b_ass_deltapB;   //!
   TBranch        *b_isGoodMatch;   //!
   TBranch        *b_ch1B;   //!
   TBranch        *b_ch2B;   //!
   TBranch        *b_chunm;   //!
   TBranch        *b_neu1B;   //!
   TBranch        *b_neu2B;   //!
   TBranch        *b_neuunm;   //!
   TBranch        *b_brecoqual;   //!
   TBranch        *b_brqual;   //!
   TBranch        *b_brecoqualangle;   //!
   TBranch        *b_chgdaugen;   //!
   TBranch        *b_neudaugen;   //!
   TBranch        *b_nchg;   //!
   TBranch        *b_nneu;   //!
   TBranch        *b_xcharge;   //!
   TBranch        *b_nB;   //!
   TBranch        *b_brecoid;   //!
   TBranch        *b_brecoidtrue;   //!
   TBranch        *b_brecoflav;   //!
   TBranch        *b_brecocharge;   //!
   TBranch        *b_modeB;   //!
   TBranch        *b_truemodeB;   //!
   TBranch        *b_isdoubleD;   //!
   TBranch        *b_mes;   //!
   TBranch        *b_mesendpoint;   //!
   TBranch        *b_de;   //!
   TBranch        *b_pB;   //!
   TBranch        *b_eB;   //!
   TBranch        *b_eUps;   //!
   TBranch        *b_pUps;   //!
   TBranch        *b_thetaUps;   //!
   TBranch        *b_phiUps;   //!
   TBranch        *b_thetaB;   //!
   TBranch        *b_phiB;   //!
   TBranch        *b_pBtrue;   //!
   TBranch        *b_eBtrue;   //!
   TBranch        *b_tBtrue;   //!
   TBranch        *b_fBtrue;   //!
   TBranch        *b_pur;   //!
   TBranch        *b_intpur;   //!
   TBranch        *b_nle_nopcut;   //!
   TBranch        *b_nle;   //!
   TBranch        *b_nel;   //!
   TBranch        *b_nmu;   //!
   TBranch        *b_nkp;   //!
   TBranch        *b_nks;   //!
   TBranch        *b_nlept500;   //!
   TBranch        *b_nelec500;   //!
   TBranch        *b_nmu500;   //!
   TBranch        *b_nlept1000;   //!
   TBranch        *b_nelec1000;   //!
   TBranch        *b_nmu1000;   //!
   TBranch        *b_deltam;   //!
   TBranch        *b_MM1pr;   //!
   TBranch        *b_MM2pr;   //!
   TBranch        *b_MM3pr;   //!
   TBranch        *b_OA1;   //!
   TBranch        *b_OA2;   //!
   TBranch        *b_OA3;   //!
   TBranch        *b_PiMin1;   //!
   TBranch        *b_PiMin2;   //!
   TBranch        *b_PiMin3;   //!
   TBranch        *b_plab;   //!
   TBranch        *b_elab;   //!
   TBranch        *b_tlab;   //!
   TBranch        *b_flab;   //!
   TBranch        *b_plabgen;   //!
   TBranch        *b_elabgen;   //!
   TBranch        *b_tlabgen;   //!
   TBranch        *b_flabgen;   //!
   TBranch        *b_lchargegen;   //!
   TBranch        *b_pcms;   //!
   TBranch        *b_ecms;   //!
   TBranch        *b_tcms;   //!
   TBranch        *b_fcms;   //!
   TBranch        *b_lcharge;   //!
   TBranch        *b_leptidgen;   //!
   TBranch        *b_leptorg;   //!
   TBranch        *b_isele;   //!
   TBranch        *b_vub;   //!
   TBranch        *b_vcb;   //!
   TBranch        *b_other;   //!
   TBranch        *b_nvubexcl;   //!
   TBranch        *b_nvubnres;   //!
   TBranch        *b_ntau;   //!
   TBranch        *b_mxhadgen;   //!
   TBranch        *b_mxhadgenwoph;   //!
   TBranch        *b_xchargegen;   //!
   TBranch        *b_pcmsgen;   //!
   TBranch        *b_ecmsgen;   //!
   TBranch        *b_tcmsgen;   //!
   TBranch        *b_fcmsgen;   //!
   TBranch        *b_pxhadgen;   //!
   TBranch        *b_exhadgen;   //!
   TBranch        *b_exhadgencms;   //!
   TBranch        *b_txhadgen;   //!
   TBranch        *b_fxhadgen;   //!
   TBranch        *b_q2Gen;   //!
   TBranch        *b_ctvgen;   //!
   TBranch        *b_ctlgen;   //!
   TBranch        *b_chigen;   //!
   TBranch        *b_enugencms;   //!
   TBranch        *b_pnugencms;   //!
   TBranch        *b_mxhad;   //!
   TBranch        *b_emiss;   //!
   TBranch        *b_pmiss;   //!
   TBranch        *b_mm2;   //!
   TBranch        *b_q2;   //!
   TBranch        *b_q2new;   //!
   TBranch        *b_exhadcms;   //!
   TBranch        *b_pxhad;   //!
   TBranch        *b_exhad;   //!
   TBranch        *b_txhad;   //!
   TBranch        *b_fxhad;   //!
   TBranch        *b_pnu;   //!
   TBranch        *b_tnu;   //!
   TBranch        *b_fnu;   //!
   TBranch        *b_pplus;   //!
   TBranch        *b_pminus;   //!
   TBranch        *b_pplusgen;   //!
   TBranch        *b_pminusgen;   //!
   TBranch        *b_pplusfit;   //!
   TBranch        *b_pminusfit;   //!
   TBranch        *b_wdeltam;   //!
   TBranch        *b_mxhadfit;   //!
   TBranch        *b_mm2fit;   //!
   TBranch        *b_q2fit;   //!
   TBranch        *b_chisq;   //!
   TBranch        *b_globchisq;   //!
   TBranch        *b_probchisq;   //!
   TBranch        *b_ndof;   //!
   TBranch        *b_fitstatus;   //!
   TBranch        *b_nnpi0;   //!
   TBranch        *b_nneu80_160;   //!
   TBranch        *b_nneu160_320;   //!
   TBranch        *b_pstarfitlept;   //!
   TBranch        *b_pfitX;   //!
   TBranch        *b_thetafitX;   //!
   TBranch        *b_phifitX;   //!
   TBranch        *b_pfitlept;   //!
   TBranch        *b_thetafitlept;   //!
   TBranch        *b_phifitlept;   //!
   TBranch        *b_pfitB;   //!
   TBranch        *b_thetafitB;   //!
   TBranch        *b_phifitB;   //!
   TBranch        *b_totweight;   //!
   TBranch        *b_totweightNutMult;   //!
   TBranch        *b_totweightTrkMult;   //!
   TBranch        *b_kplus;   //!
   TBranch        *b_nBrems;   //!
   TBranch        *b_eBrems;   //!
   TBranch        *b_mxhadfit0;   //!
   TBranch        *b_q2fit0;   //!
   TBranch        *b_pplusfit0;   //!
   TBranch        *b_chisqfit0;   //!
   TBranch        *b_fitstatusfit0;   //!
   TBranch        *b_chisqfit;   //!
   TBranch        *b_fitstatusfit;   //!
   TBranch        *b_mxhadfit1;   //!
   TBranch        *b_q2fit1;   //!
   TBranch        *b_pplusfit1;   //!
   TBranch        *b_chisqfit1;   //!
   TBranch        *b_fitstatusfit1;   //!
   TBranch        *b_mxhadfit2;   //!
   TBranch        *b_q2fit2;   //!
   TBranch        *b_pplusfit2;   //!
   TBranch        *b_chisqfit2;   //!
   TBranch        *b_fitstatusfit2;   //!
   TBranch        *b_mxhadfit3;   //!
   TBranch        *b_q2fit3;   //!
   TBranch        *b_pplusfit3;   //!
   TBranch        *b_chisqfit3;   //!
   TBranch        *b_fitstatusfit3;   //!
   TBranch        *b_mxhadfit4;   //!
   TBranch        *b_q2fit4;   //!
   TBranch        *b_pplusfit4;   //!
   TBranch        *b_chisqfit4;   //!
   TBranch        *b_fitstatusfit4;   //!
   TBranch        *b_mxhadfit5;   //!
   TBranch        *b_q2fit5;   //!
   TBranch        *b_pplusfit5;   //!
   TBranch        *b_chisqfit5;   //!
   TBranch        *b_fitstatusfit5;   //!
   TBranch        *b_mxhadfit6;   //!
   TBranch        *b_q2fit6;   //!
   TBranch        *b_pplusfit6;   //!
   TBranch        *b_chisqfit6;   //!
   TBranch        *b_fitstatusfit6;   //!
   TBranch        *b_mxhadfit7;   //!
   TBranch        *b_q2fit7;   //!
   TBranch        *b_pplusfit7;   //!
   TBranch        *b_chisqfit7;   //!
   TBranch        *b_fitstatusfit7;   //!
   TBranch        *b_mxhadfit8;   //!
   TBranch        *b_q2fit8;   //!
   TBranch        *b_pplusfit8;   //!
   TBranch        *b_chisqfit8;   //!
   TBranch        *b_fitstatusfit8;   //!

   double themin;
   double themax;
   int thebins;
   TString thevar;
   int multipl;
   double intdataen;
   double intdatadepl;
   double intMCen;
   double intMCdepl;   
   double SHIFTNEUT;
   double SIGMANEUT;

   RooRealVar* Thor;
   RooRealVar* Thosigma_r1;
   RooRealVar* Thoxc;
   RooRealVar* Thosigma_r2;
   RooRealVar* Thosigma_l;
   RooRealVar* Thon;
   RooRealVar* Thoalpha;
   RooRealVar* Rm;
   RooRealVar* Rs;
   RooRealVar* Ra;
   RooRealVar* Rn;
   RooRealVar* Rendpoint;
   RooRealVar* pArgPar;
   RooRealVar* pCutOff;

   TPad *fPads[50];
   TCanvas *c0, *c1, *c2, *c3, *c4, *c5, *c6, *c7;
   thecomparison();
   thecomparison(char *var, double mi, double ma, int b);
   virtual ~thecomparison();
   Int_t  Cut(Int_t entry);
   Int_t  GetEntry(Int_t entry);
   Int_t  LoadTree(Int_t entry);

   double    getBsysweight(int,int);
   double    getDsysweight(int,int,int,int,int,int,int);
   double    FermiWeight(double,double,double);
   double    fermi(double,double,double);

   void   Init(TTree *tree);
   void   Loop(int nevents, int cat, double shift = 0, double smear = 0, int isbch = 2, int multcat = 7, int seed = 1,int sys=0);
   // seed = 2 dstar
   // seed = 3 dstar0
   // seed = 4 dc
   // seed = 5 d0
   // sys = 0: nothing; 1: B&D defaults; 2 fermi reweighting (to values set by hand still...)
   void   Bookhist();
   Bool_t Notify();
   void   Show(Int_t entry = -1);
   int    hist(double mx);
   void   Fitmes(int cat, int cut,bool istmodel);
   void sighisto(double& signal, double& signalErr,TH1D *histo, double &resmean, double &ressigma, double &resalpha, double &resn, int fixpar,double mean, double sigma, double alpha, double n, double argus);
   void sighisto_newmodel(double&, double&, TH1D* ,vector<double>&, int, const vector<double>&,bool);
   void   overlap(int cut, int norm, TString dir);
   void   effplots(TString dir);
   TChain* getchain(char *thechain, char *treename);
   double chisquared(TH1 *h1, TH1 *h2);
   Double_t  smeargauss(double,double,double);
   mesData* vubMesUnb(RooDataHist *, RooRealVar *, vector<double>&, int , const vector<double>&,bool);
   
   RooAbsPdf* createThorsten(RooRealVar& mes);
   RooAbsPdf* createCCB(RooRealVar& mes);
   RooAbsPdf* createArgus(RooRealVar& mes);





   ClassDef(thecomparison,1)

};
Double_t getVal(RooAbsPdf* pdf, const char* name);
RooRealVar* getPointer(RooAbsPdf* pdf, const char* name);

#endif



// HERE THERE IS THE VARIABLES LIST FOR NTUPLES BEFORE SUMMER 06 PRODUCTION

//Declaration of leaves types
/* /\* */
/*   OLD */

/* Double_t        bmass; */
/* Double_t        bmassfit; */
/* UChar*_t         sbox; */
/* Int_t           brecomc; */
/* UChar_t         GoodEvent; */
/* UChar_t         isDupli; */
/* Int_t           ValMap; */
/* Int_t           bgcat; */
/* Double_t        gmax; */
/* Int_t           nneufromB; */
/* Int_t           nneufromB80_160; */
/* Int_t           nneufromB160_320; */
/* Double_t        mm2nc; */
/* Double_t        allksm0[18]; */
/* Double_t        allksp[18]; */
/* UChar_t         allksmc[18]; */
/* Double_t        allchkp[7]; */
/* UChar_t         allchkmc[7]; */
/* Double_t        m0ks; */
/* Double_t        pks; */
/* Double_t        pksmc; */
/* Int_t           ntkl; */
/* Double_t        tklp[0]; */
/* Double_t        tklth[0]; */
/* Double_t        tklph[0]; */
/* UChar_t         tklisol[0]; */
/* Int_t           ntks; */
/* Double_t        tksp[0]; */
/* Double_t        tksth[0]; */
/* Double_t        tksph[0]; */
/* Int_t           tksdec[0]; */
/* Int_t           ntchk; */
/* Double_t        tchkp[0]; */
/* Double_t        tchkth[0]; */
/* Double_t        tchkph[0]; */
/* Int_t           nklres; */
/* Double_t        klresth[0]; */
/* Double_t        klresph[0]; */
/* Int_t           klid[0]; */
/* UChar_t         klcone[0]; */
/* Double_t        emckl; */
/* Double_t        emckl0; */
/* Double_t        emckl22; */
/* Double_t        mxks; */
/* Double_t        mm2ks; */
/* Double_t        mxksfit; */
/* Double_t        mm2misk; */
/* Double_t        mxmisk; */
/* Double_t        mxmiskfit; */
/* Double_t        mm2mchk; */
/* Double_t        mxchk; */
/* Double_t        mxchkfit; */
/* *\/ */
/*    Int_t           run;    */
/*    Int_t           lower; */
/*    Int_t           upper; */
/*    Double_t        mes; */
/*    Double_t        de; */
/*    Double_t        pur; */
/*    Double_t        intpur; */
/*    Int_t           modeB; */
/*    Int_t           nnpi0; */
/*    Int_t           brecoflav; */
/*    Int_t           brecocharge; */
/*    Double_t        mxhadgen; */
/*    Double_t        pcmsgen; */
/*    Double_t        tcmsgen; */
/*    Double_t        fcmsgen; */
/*    Double_t        ecmsgen; */
/*    Double_t        pxhadgen; */
/*    Double_t        txhadgen; */
/*    Double_t        fxhadgen; */
/*    Double_t        exhadgen; */
/*    Double_t        fkplus; */
/*    Int_t           vub; */
/*    Int_t           vcb;  */
/*    Int_t           Gvxbtyp; */
/*    Int_t           other; */
/*    Int_t           xcharge; */
/*    Double_t        pxhad; */
/*    Double_t        txhad; */
/*    Double_t        fxhad; */
/*    Double_t        exhad; */
/*    Double_t        mxhad; */
/*    Double_t        mxhadfit; */
   
/* /\*   Double_t        csiCiuc; */
/*    Double_t        xCiuc; */
/*    Double_t        wCiuc; */
/*    Double_t        EwPwfit; */
/*    Double_t        GcsiCiuc; */
/*    Double_t        GxCiuc; */
/*    Double_t        GwCiuc; */
/*    Double_t        EwPwG; */
/* *\/ */

/*    Double_t        q2fit; */
/*    Double_t        q2; */
/*    Double_t        q2Gen; */
/*    Int_t           lcharge; */
/*    Double_t        plab; */
/*    Double_t        tlab; */
/*    Double_t        flab; */
/*    Double_t        pcms; */
/*    Double_t        tcms; */
/*    Double_t        fcms; */
/*    Double_t        ecms; */
/*    Int_t           nle; */
/*    Int_t           nel; */
/*    Int_t           nmu; */
/*    Int_t           nchg; */
/*    Int_t           nneu; */
/*    Int_t           nneu80_160; */
/*    Int_t           nneu160_320; */
/*    Int_t           nkp; */
/*    Int_t           nks; */
/*    //   Int_t           npi0; */
/*    Int_t           GSem; */
/*    Int_t           GfDpi; */
/*    Int_t           GfDpiz; */
/*    Int_t           GfDk; */
/*    Int_t           GfDks; */
/*    Int_t           GfDlep; */
/*    Int_t           GfDgam; */
 
/*   /\*   Double_t        dx; */
/*    Double_t        dy; */
/*    Double_t        dz;  */
/*    Double_t        s2dxx; */
/*    Double_t        s2dyy; */
/*    Double_t        s2dzz; */
/*    Double_t        s2dxy; */
/*    Double_t        s2dyz; */
/*    Double_t        s2dxz; */
/*    Double_t        eneu; */
/*    Double_t        epiz; */
/*    Double_t        kminmom; */
/*    Double_t        kmaxmom; */
/* *\/ */
/*    Double_t        pnu; */
/*    Double_t        tnu; */
/*    Double_t        fnu; */
/*    Double_t        mm2; */
/*    Double_t        mm2fit; */
/*    Double_t        totweight; */
/*    Double_t        totweightNutMult; */
/*    Double_t        totweightTrkMult; */

/* //List of branches */
/*    TBranch        *b_run; */
/*    TBranch        *b_lower; */
/*    TBranch        *b_upper; */
/*    /\* */
/*    TBranch        *b_bmass; */
/*    TBranch        *b_bmassfit; */
/*    TBranch        *b_sbox;*\/ */
/*    TBranch        *b_mes; */
/*    TBranch        *b_de; */
/*    TBranch        *b_pur; */
/*    TBranch        *b_intpur; */
/*    TBranch        *b_modeB; */
/*    TBranch        *b_nnpi0; */
/*    TBranch        *b_brecoflav; */
/*    TBranch        *b_brecocharge; */
/*    //   TBranch        *b_brecomc; */
/*    TBranch        *b_mxhadgen; */
/*    TBranch        *b_pcmsgen; */
/*    TBranch        *b_tcmsgen; */
/*    TBranch        *b_fcmsgen; */
/*    TBranch        *b_ecmsgen; */
/*    TBranch        *b_pxhadgen; */
/*    TBranch        *b_txhadgen; */
/*    TBranch        *b_fxhadgen; */
/*    TBranch        *b_exhadgen; */
/*    TBranch        *b_kplus; */
/*    /\*   TBranch        *b_GoodEvent; */
/*    TBranch        *b_isDupli; */
/*    TBranch        *b_ValMap;*\/ */
/*    TBranch        *b_vub; */
/*    TBranch        *b_vcb; */
/*    TBranch        *b_Gvxbtyp; */
/*    TBranch        *b_other; */
/*    //   TBranch        *b_bgcat; */
/*    TBranch        *b_xcharge; */
/*    TBranch        *b_pxhad; */
/*    TBranch        *b_txhad; */
/*    TBranch        *b_fxhad; */
/*    TBranch        *b_exhad; */
/*    TBranch        *b_mxhad; */
/*    //   TBranch        *b_gmax; */
/*    TBranch        *b_mxhadfit; */
/*    TBranch        *b_q2fit; */
/*    TBranch        *b_q2; */
/*    /\*TBranch        *b_EwPwfit; */
/*      TBranch        *b_csiCiuc; */
/*    TBranch        *b_xCiuc; */
/*    TBranch        *b_wCiuc; */
/*    TBranch        *b_GcsiCiuc;  */
/*    TBranch        *b_EwPwG; */
/*    TBranch        *b_GxCiuc; */
/*    TBranch        *b_GwCiuc;*\/ */
/*    TBranch        *b_q2Gen; */
/*    TBranch        *b_lcharge; */
/*    TBranch        *b_plab; */
/*    TBranch        *b_tlab; */
/*    TBranch        *b_flab; */
/*    TBranch        *b_pcms; */
/*    TBranch        *b_tcms; */
/*    TBranch        *b_fcms; */
/*    TBranch        *b_ecms; */
/*    TBranch        *b_nle; */
/*    TBranch        *b_nel; */
/*    TBranch        *b_nmu; */
/*    TBranch        *b_nchg; */
/*    TBranch        *b_nneu; */
/*    TBranch        *b_nneu80_160; */
/*    TBranch        *b_nneu160_320; */
/*    //TBranch        *b_nneufromB; */
/*    //TBranch        *b_nneufromB80_160; */
/*    //TBranch        *b_nneufromB160_320; */
/*    TBranch        *b_nkp; */
/*    TBranch        *b_nks; */
/*    //TBranch        *b_npi0; */
/*    TBranch        *b_GSem;     */
/*    TBranch        *b_GfDpi;    */
/*    TBranch        *b_GfDpiz;   */
/*    TBranch        *b_GfDk;     */
/*    TBranch        *b_GfDks;    */
/*    TBranch        *b_GfDlep;   */
/*    TBranch        *b_GfDgam;   */
/*    /\*TBranch        *b_dx; */
/*    TBranch        *b_dy; */
/*    TBranch        *b_dz; */
/*    TBranch        *b_s2dxx; */
/*    TBranch        *b_s2dyy; */
/*    TBranch        *b_s2dzz; */
/*    TBranch        *b_s2dxy; */
/*    TBranch        *b_s2dyz; */
/*    TBranch        *b_s2dxz; */
/*    TBranch        *b_epiz; */
/*    TBranch        *b_kminmom; */
/*    TBranch        *b_kmaxmom; */
/*    TBranch        *b_mm2nc; */
/* *\/ */
/*    TBranch        *b_pnu; */
/*    TBranch        *b_tnu; */
/*    TBranch        *b_fnu; */
/*    //TBranch        *b_eneu; */

/*    TBranch        *b_mm2; */

/*    TBranch        *b_mm2fit; */
/*    /\*TBranch        *b_allksm0; */
/*    TBranch        *b_allksp; */
/*    TBranch        *b_allksmc; */
/*    TBranch        *b_allchkp; */
/*    TBranch        *b_allchkmc; */
/*    TBranch        *b_m0ks; */
/*    TBranch        *b_pks; */
/*    TBranch        *b_pksmc; */
/*    TBranch        *b_ntkl; */
/*    TBranch        *b_tklp; */
/*    TBranch        *b_tklth; */
/*    TBranch        *b_tklph; */
/*    TBranch        *b_tklisol; */
/*    TBranch        *b_ntks; */
/*    TBranch        *b_tksp; */
/*    TBranch        *b_tksth; */
/*    TBranch        *b_tksph; */
/*    TBranch        *b_tksdec; */
/*    TBranch        *b_ntchk; */
/*    TBranch        *b_tchkp; */
/*    TBranch        *b_tchkth; */
/*    TBranch        *b_tchkph; */
/*    TBranch        *b_nklres; */
/*    TBranch        *b_klresth; */
/*    TBranch        *b_klresph; */
/*    TBranch        *b_klid; */
/*    TBranch        *b_klcone; */
/*    TBranch        *b_emckl; */
/*    TBranch        *b_emckl0; */
/*    TBranch        *b_emckl22; */
/*    TBranch        *b_mxks; */
/*    TBranch        *b_mm2ks; */
/*    TBranch        *b_mxksfit; */
/*    TBranch        *b_mm2misk; */
/*    TBranch        *b_mxmisk; */
/*    TBranch        *b_mxmiskfit; */
/*    TBranch        *b_mm2chk; */
/*    TBranch        *b_mxchk; */
/*    TBranch        *b_mxchkfit;*\/ */
/*    TBranch        *b_totweight; */
/*    TBranch        *b_totweightNutMult; */
/*    TBranch        *b_totweightTrkMult; */
