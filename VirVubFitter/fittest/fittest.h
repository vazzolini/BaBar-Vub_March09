/////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jan 19 03:22:16 2006 by ROOT version 4.01/02
// from TChain ntp1/events chain
//////////////////////////////////////////////////////////

#ifndef fittest_h
#define fittest_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

gROOT->SetStyle("Plain");

class fittest {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

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
   Float_t         eBrems[2];   //[nBrems]
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

   // Other stuff

  RooRealVar*  Vmes;   
  RooRealVar*  VlepYes; 
  RooRealVar*  VlepVub;
  RooRealVar*  VlepVcb; 
  RooRealVar*  VlepVubSB;
  RooRealVar*  VflavB;   
  RooRealVar*  VlepYaSe; 
  RooRealVar*  Vwe;      
  RooRealVar*  Vchop;
  RooRealVar*  Vmx;
  RooRealVar*  Vq2; 
  RooRealVar*  Vmultcat;
  RooRealVar*  Vksele;
  RooRealVar*  Vpplus;
  RooRealVar*  Vintpur;
  //  RooRealVar*  Vtrumtch;
  RooRealVar*  Visch;
  RooRealVar*  Vmm2;
  RooRealVar*  Vch;
  RooRealVar*  VmodeB;
  RooRealVar*  VtruemodeB;
  RooRealVar*  Vch1B;
  RooRealVar*  Vneu1B;
  RooRealVar*  Vwrongcharge;

  RooExtendPdf *ae,*se;

  RooDataSet* datadata; 

  //MODEL STUFF IN HERE
  
  RooRealVar *argpar(NULL), *cutoff(NULL), *argusfrac(NULL);
  RooRealVar *Rm(NULL), *Rs(NULL), *Ra(NULL), *Rn(NULL), *Rendpoint(NULL);
  RooRealVar *Thosigma_L(NULL), *Thosigma_r1(NULL), *Thosigma_r2(NULL), *Thor(NULL), *Thon(NULL), *Thoalpha(NULL), *Thoxc(NULL), *sigpdffraction(NULL);
  RooRealVar *CB_mean(NULL), *CB_sigma(NULL), *CB_alpha(NULL), *CB_n(NULL);
  RooRealVar *x(NULL);
  RooArgusBG *a(NULL);
  //  ccb *ccb(NULL);
  thosig *thosig(NULL);
  RooAddPdf* tmodel(NULL);
  RooCBShape *crystalball(NULL);

  //Utils here
  TString *fname, *somecut;
  bool dumppar;
  FILE *fout;
  TString parfile;
  double scala;
  RooRealVar  *tempargy(NULL), *tempsigy(NULL), *tempcrystalbally(NULL);
  
  enum PdfModel {oldmodel = 0, newmodel = 1};
  enum PdfType  {sigonly=1, bkgonly=2, all=3};
  enum fileType {ccbar=0, udsudsbar=1, BBbar=2, Data=3};

  fittest(TTree *tree=0);
  virtual ~fittest();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     WriteDataSet(Int_t ,Int_t);
  virtual void     Loop(Int_t);
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
  virtual void     Test(bool,int,int,bool,bool,bool dumpparams=false,const char*,float intp = 0.,int chb = 2, float mx_low = -300, float mx_high = 300,float q2_low=-999, float q2_high=300, float ppl_low=-100, float ppl_high=100);
  virtual void     BuildModel(bool,int,int,bool);
  virtual void     ResetArgus(bool,int,bool);
  //  virtual void     ResetCCB(bool,int,bool,float);
  virtual void     ResetThosig(bool,int,bool);
  virtual void     ResetCrystalBall(bool,int,bool);
  virtual void     NewTree(Int_t,Int_t);
  virtual RooRealVar* readParameters(const char*, const char*);
};

#endif

#ifdef fittest_cxx
fittest::fittest(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
    if (!f) {
      f = new TFile("Memory Directory");
      f->cd("Rint:/");
    }
    tree = (TTree*)gDirectory->Get("ntp1");

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
    TChain * chain = new TChain("ntp1","events chain");
    
    tree=chain;
      
#endif // SINGLE_TREE
    
  }
  Init(tree);
}

fittest::NewTree(Int_t select, Int_t runperiod)
{
  if(fChain != NULL)
    delete fChain;

  TChain *chain = new TChain("ntp1","events chain");
  
  switch(select){
  case ccbar: //ccbar
    cout<<"NOW LOOPING ON CCBAR CHAINS"<<endl;
    /*     chain->Add("/nfs/farm/babar/AWGsemilep01/menges/summer06/store.2/Vubrec-def/reduced/ccbarRedRun1.root"); */
    /*     chain->Add("/nfs/farm/babar/AWGsemilep01/menges/summer06/store.2/Vubrec-def/reduced/ccbarRedRun2.root"); */
    /*     chain->Add("/nfs/farm/babar/AWGsemilep01/menges/summer06/store.2/Vubrec-def/reduced/ccbarRedRun3.root"); */
    /*     chain->Add("/nfs/farm/babar/AWGsemilep01/menges/summer06/store.2/Vubrec-def/reduced/ccbarRedRun4.root"); */
    break;
  case udsudsbar: //uds
    cout<<"Now Looping on  UDS CHAIN"<<endl;
    /* chain->Add("/nfs/farm/babar/AWGsemilep01/menges/summer06/store.2/Vubrec-def/reduced/udsRedRun1.root"); */ 
    /* chain->Add("/nfs/farm/babar/AWGsemilep01/menges/summer06/store.2/Vubrec-def/reduced/udsRedRun2.root"); */
    /* chain->Add("/nfs/farm/babar/AWGsemilep01/menges/summer06/store.2/Vubrec-def/reduced/udsRedRun3.root"); */
    /* chain->Add("/nfs/farm/babar/AWGsemilep01/menges/summer06/store.2/Vubrec-def/reduced/udsRedRun4.root"); */
    break;
  case BBbar: //BBbar
    cout<<"Now Looping on BBBAR CHAIN"<<endl;
    switch(runperiod){
    case 12: 
      cout<<" Run 1-2"<<endl;
      chain->Add("/nfs/farm/babar/AWGsemilep01/petrella/Rootples/1111-R18/reduced/genRedChainR18-Run1-2.root");
      break;
    case 3: 
      cout<<" Run3"<<endl;
      chain->Add("/nfs/farm/babar/AWGsemilep01/petrella/Rootples/1111-R18/reduced/genRedChainR18-Run3.root"); 
      break;
    case 4: 
      cout<<" Run4"<<endl;
      chain->Add("/nfs/farm/babar/AWGsemilep01/petrella/Rootples/1111-R18/reduced/genRedChainR18-Run4.root");
      break;
    case 5:
      cout<<" Run5"<<endl;
      chain->Add("/nfs/farm/babar/AWGsemilep01/petrella/Rootples/1111-R18/reduced/genRedChainR18-Run5.root");
      break;
    case 14:
      cout<<" Run1-4"<<endl;
      chain->Add("/nfs/farm/babar/AWGsemilep01/petrella/Rootples/1111-R18/reduced/genbchRedChainR18-Run1-4.root");
      chain->Add("/nfs/farm/babar/AWGsemilep01/petrella/Rootples/1111-R18/reduced/genbnuRedChainR18-Run1-4.root");
      break;
    case 15:
      cout<< "Run1-5"<<endl;
      chain->Add("/nfs/farm/babar/AWGsemilep01/petrella/Rootples/1111-R18/reduced/genRedChainR18-Run1-2.root");
      chain->Add("/nfs/farm/babar/AWGsemilep01/petrella/Rootples/1111-R18/reduced/genRedChainR18-Run3.root");
      chain->Add("/nfs/farm/babar/AWGsemilep01/petrella/Rootples/1111-R18/reduced/genRedChainR18-Run4.root");
      chain->Add("/nfs/farm/babar/AWGsemilep01/petrella/Rootples/1111-R18/reduced/genRedChainR18-Run5.root");
      break;
    }  break; //End of BBbar
  case Data : //Data
    cout<<"Now Looping On Data"<<endl;
    switch(runperiod){
    case 12: 
      cout<<" Run 1-2 Data"<<endl;
      chain->Add("/nfs/farm/babar/AWGsemilep01/petrella/Rootples/1111-R18/reduced/OnPeakRedChainR18-Run1-2.root");
      break;
    case 3: 
      cout<<" Run3 Data"<<endl;
      chain->Add("/nfs/farm/babar/AWGsemilep01/petrella/Rootples/1111-R18/reduced/OnPeakRedChainR18-Run3.root");
      break;
    case 4: 
      cout<<" Run4 Data"<<endl;
      chain->Add("/nfs/farm/babar/AWGsemilep01/petrella/Rootples/1111-R18/reduced/OnPeakRedChainR18-Run4.root");
      break;
    case 5: 
      cout<<" Run5 Data"<<endl;
      chain->Add("/nfs/farm/babar/AWGsemilep01/petrella/Rootples/1111-R18/reduced/OnPeakRedChainR18-Run5.root");
      break;
    case 14: 
      cout<<" Run1-4"<<endl;
      chain->Add("/nfs/farm/babar/AWGsemilep01/petrella/Rootples/1111-R18/reduced/OnPeakRedChainR18-Run1-4.root");
      break;
    case 15: 
      cout<<" Run1-5"<<endl;
      chain->Add("/nfs/farm/babar/AWGsemilep01/petrella/Rootples/1111-R18/reduced/OnPeakRedChainR18-Run1-2.root");
      chain->Add("/nfs/farm/babar/AWGsemilep01/petrella/Rootples/1111-R18/reduced/OnPeakRedChainR18-Run3.root");
      chain->Add("/nfs/farm/babar/AWGsemilep01/petrella/Rootples/1111-R18/reduced/OnPeakRedChainR18-Run4.root");
      chain->Add("/nfs/farm/babar/AWGsemilep01/petrella/Rootples/1111-R18/reduced/OnPeakRedChainR18-Run5.root");
      break;
    default:
      cout<<" Run Period Not Valid: "<<runperiod<<endl;
      chain=NULL;
      break;
    } break; 
  default:
    cout<<" MonteCarlo Not Valid!"<<endl;
    chain=NULL;
    break;
  }
  Init(chain);
}
fittest::~fittest()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t fittest::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t fittest::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->IsA() != TChain::Class()) return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void fittest::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses of the tree
   // will be set. It is normaly not necessary to make changes to the
   // generated code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running with PROOF.

   // Set branch addresses
   if (tree == 0) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run",&run);
   fChain->SetBranchAddress("upper",&upper);
   fChain->SetBranchAddress("lower",&lower);
   fChain->SetBranchAddress("Gvxbtyp",&Gvxbtyp);
   fChain->SetBranchAddress("GSem",&GSem);
   fChain->SetBranchAddress("GfDpi",&GfDpi);
   fChain->SetBranchAddress("GfDpiz",&GfDpiz);
   fChain->SetBranchAddress("GfDk",&GfDk);
   fChain->SetBranchAddress("GfDks",&GfDks);
   fChain->SetBranchAddress("GfDkl",&GfDkl);
   fChain->SetBranchAddress("GfDlep",&GfDlep);
   fChain->SetBranchAddress("GfDgam",&GfDgam);
   fChain->SetBranchAddress("GfDnu",&GfDnu);
   fChain->SetBranchAddress("GfD0Ds",&GfD0Ds);
   fChain->SetBranchAddress("GfDDs",&GfDDs);
   fChain->SetBranchAddress("GfDkspiopio",&GfDkspiopio);
   fChain->SetBranchAddress("GfK",&GfK);
   fChain->SetBranchAddress("isassocB",&isassocB);
   fChain->SetBranchAddress("isassocB_GHIT",&isassocB_GHIT);
   fChain->SetBranchAddress("ass_deltapB",&ass_deltapB);
   fChain->SetBranchAddress("isGoodMatch",&isGoodMatch);
   fChain->SetBranchAddress("ch1B",&ch1B);
   fChain->SetBranchAddress("ch2B",&ch2B);
   fChain->SetBranchAddress("chunm",&chunm);
   fChain->SetBranchAddress("neu1B",&neu1B);
   fChain->SetBranchAddress("neu2B",&neu2B);
   fChain->SetBranchAddress("neuunm",&neuunm);
   fChain->SetBranchAddress("brecoqual",&brecoqual);
   fChain->SetBranchAddress("brqual",&brqual);
   fChain->SetBranchAddress("brecoqualangle",&brecoqualangle);
   fChain->SetBranchAddress("chgdaugen",&chgdaugen);
   fChain->SetBranchAddress("neudaugen",&neudaugen);
   fChain->SetBranchAddress("nchg",&nchg);
   fChain->SetBranchAddress("nneu",&nneu);
   fChain->SetBranchAddress("xcharge",&xcharge);
   fChain->SetBranchAddress("nB",&nB);
   fChain->SetBranchAddress("brecoid",&brecoid);
   fChain->SetBranchAddress("brecoidtrue",&brecoidtrue);
   fChain->SetBranchAddress("brecoflav",&brecoflav);
   fChain->SetBranchAddress("brecocharge",&brecocharge);
   fChain->SetBranchAddress("modeB",&modeB);
   fChain->SetBranchAddress("truemodeB",&truemodeB);
   fChain->SetBranchAddress("isdoubleD",&isdoubleD);
   fChain->SetBranchAddress("mes",&mes);
   fChain->SetBranchAddress("mesendpoint",&mesendpoint);
   fChain->SetBranchAddress("de",&de);
   fChain->SetBranchAddress("pB",&pB);
   fChain->SetBranchAddress("eB",&eB);
   fChain->SetBranchAddress("eUps",&eUps);
   fChain->SetBranchAddress("pUps",&pUps);
   fChain->SetBranchAddress("thetaUps",&thetaUps);
   fChain->SetBranchAddress("phiUps",&phiUps);
   fChain->SetBranchAddress("thetaB",&thetaB);
   fChain->SetBranchAddress("phiB",&phiB);
   fChain->SetBranchAddress("pBtrue",&pBtrue);
   fChain->SetBranchAddress("eBtrue",&eBtrue);
   fChain->SetBranchAddress("tBtrue",&tBtrue);
   fChain->SetBranchAddress("fBtrue",&fBtrue);
   fChain->SetBranchAddress("pur",&pur);
   fChain->SetBranchAddress("intpur",&intpur);
   fChain->SetBranchAddress("nle_nopcut",&nle_nopcut);
   fChain->SetBranchAddress("nle",&nle);
   fChain->SetBranchAddress("nel",&nel);
   fChain->SetBranchAddress("nmu",&nmu);
   fChain->SetBranchAddress("nkp",&nkp);
   fChain->SetBranchAddress("nks",&nks);
   fChain->SetBranchAddress("nlept500",&nlept500);
   fChain->SetBranchAddress("nelec500",&nelec500);
   fChain->SetBranchAddress("nmu500",&nmu500);
   fChain->SetBranchAddress("nlept1000",&nlept1000);
   fChain->SetBranchAddress("nelec1000",&nelec1000);
   fChain->SetBranchAddress("nmu1000",&nmu1000);
   fChain->SetBranchAddress("deltam",&deltam);
   fChain->SetBranchAddress("MM1pr",&MM1pr);
   fChain->SetBranchAddress("MM2pr",&MM2pr);
   fChain->SetBranchAddress("MM3pr",&MM3pr);
   fChain->SetBranchAddress("OA1",&OA1);
   fChain->SetBranchAddress("OA2",&OA2);
   fChain->SetBranchAddress("OA3",&OA3);
   fChain->SetBranchAddress("PiMin1",&PiMin1);
   fChain->SetBranchAddress("PiMin2",&PiMin2);
   fChain->SetBranchAddress("PiMin3",&PiMin3);
   fChain->SetBranchAddress("plab",&plab);
   fChain->SetBranchAddress("elab",&elab);
   fChain->SetBranchAddress("tlab",&tlab);
   fChain->SetBranchAddress("flab",&flab);
   fChain->SetBranchAddress("plabgen",&plabgen);
   fChain->SetBranchAddress("elabgen",&elabgen);
   fChain->SetBranchAddress("tlabgen",&tlabgen);
   fChain->SetBranchAddress("flabgen",&flabgen);
   fChain->SetBranchAddress("lchargegen",&lchargegen);
   fChain->SetBranchAddress("pcms",&pcms);
   fChain->SetBranchAddress("ecms",&ecms);
   fChain->SetBranchAddress("tcms",&tcms);
   fChain->SetBranchAddress("fcms",&fcms);
   fChain->SetBranchAddress("lcharge",&lcharge);
   fChain->SetBranchAddress("leptidgen",&leptidgen);
   fChain->SetBranchAddress("leptorg",&leptorg);
   fChain->SetBranchAddress("isele",&isele);
   fChain->SetBranchAddress("vub",&vub);
   fChain->SetBranchAddress("vcb",&vcb);
   fChain->SetBranchAddress("other",&other);
   fChain->SetBranchAddress("nvubexcl",&nvubexcl);
   fChain->SetBranchAddress("nvubnres",&nvubnres);
   fChain->SetBranchAddress("ntau",&ntau);
   fChain->SetBranchAddress("mxhadgen",&mxhadgen);
   fChain->SetBranchAddress("mxhadgenwoph",&mxhadgenwoph);
   fChain->SetBranchAddress("xchargegen",&xchargegen);
   fChain->SetBranchAddress("pcmsgen",&pcmsgen);
   fChain->SetBranchAddress("ecmsgen",&ecmsgen);
   fChain->SetBranchAddress("tcmsgen",&tcmsgen);
   fChain->SetBranchAddress("fcmsgen",&fcmsgen);
   fChain->SetBranchAddress("pxhadgen",&pxhadgen);
   fChain->SetBranchAddress("exhadgen",&exhadgen);
   fChain->SetBranchAddress("exhadgencms",&exhadgencms);
   fChain->SetBranchAddress("txhadgen",&txhadgen);
   fChain->SetBranchAddress("fxhadgen",&fxhadgen);
   fChain->SetBranchAddress("q2Gen",&q2Gen);
   fChain->SetBranchAddress("ctvgen",&ctvgen);
   fChain->SetBranchAddress("ctlgen",&ctlgen);
   fChain->SetBranchAddress("chigen",&chigen);
   fChain->SetBranchAddress("enugencms",&enugencms);
   fChain->SetBranchAddress("pnugencms",&pnugencms);
   fChain->SetBranchAddress("mxhad",&mxhad);
   fChain->SetBranchAddress("emiss",&emiss);
   fChain->SetBranchAddress("pmiss",&pmiss);
   fChain->SetBranchAddress("mm2",&mm2);
   fChain->SetBranchAddress("q2",&q2);
   fChain->SetBranchAddress("q2new",&q2new);
   fChain->SetBranchAddress("exhadcms",&exhadcms);
   fChain->SetBranchAddress("pxhad",&pxhad);
   fChain->SetBranchAddress("exhad",&exhad);
   fChain->SetBranchAddress("txhad",&txhad);
   fChain->SetBranchAddress("fxhad",&fxhad);
   fChain->SetBranchAddress("pnu",&pnu);
   fChain->SetBranchAddress("tnu",&tnu);
   fChain->SetBranchAddress("fnu",&fnu);
   fChain->SetBranchAddress("pplus",&pplus);
   fChain->SetBranchAddress("pminus",&pminus);
   fChain->SetBranchAddress("pplusgen",&pplusgen);
   fChain->SetBranchAddress("pminusgen",&pminusgen);
   fChain->SetBranchAddress("pplusfit",&pplusfit);
   fChain->SetBranchAddress("pminusfit",&pminusfit);
   fChain->SetBranchAddress("wdeltam",&wdeltam);
   fChain->SetBranchAddress("mxhadfit",&mxhadfit);
   fChain->SetBranchAddress("mm2fit",&mm2fit);
   fChain->SetBranchAddress("q2fit",&q2fit);
   fChain->SetBranchAddress("chisq",&chisq);
   fChain->SetBranchAddress("globchisq",&globchisq);
   fChain->SetBranchAddress("probchisq",&probchisq);
   fChain->SetBranchAddress("ndof",&ndof);
   fChain->SetBranchAddress("fitstatus",&fitstatus);
   fChain->SetBranchAddress("nnpi0",&nnpi0);
   fChain->SetBranchAddress("nneu80_160",&nneu80_160);
   fChain->SetBranchAddress("nneu160_320",&nneu160_320);
   fChain->SetBranchAddress("pstarfitlept",&pstarfitlept);
   fChain->SetBranchAddress("pfitX",&pfitX);
   fChain->SetBranchAddress("thetafitX",&thetafitX);
   fChain->SetBranchAddress("phifitX",&phifitX);
   fChain->SetBranchAddress("pfitlept",&pfitlept);
   fChain->SetBranchAddress("thetafitlept",&thetafitlept);
   fChain->SetBranchAddress("phifitlept",&phifitlept);
   fChain->SetBranchAddress("pfitB",&pfitB);
   fChain->SetBranchAddress("thetafitB",&thetafitB);
   fChain->SetBranchAddress("phifitB",&phifitB);
   fChain->SetBranchAddress("totweight",&totweight);
   fChain->SetBranchAddress("totweightNutMult",&totweightNutMult);
   fChain->SetBranchAddress("totweightTrkMult",&totweightTrkMult);
   fChain->SetBranchAddress("kplus",&kplus);
   fChain->SetBranchAddress("nBrems",&nBrems);
   fChain->SetBranchAddress("eBrems",eBrems);
   fChain->SetBranchAddress("mxhadfit0",&mxhadfit0);
   fChain->SetBranchAddress("q2fit0",&q2fit0);
   fChain->SetBranchAddress("pplusfit0",&pplusfit0);
   fChain->SetBranchAddress("chisqfit0",&chisqfit0);
   fChain->SetBranchAddress("fitstatusfit0",&fitstatusfit0);
   fChain->SetBranchAddress("chisqfit",&chisqfit);
   fChain->SetBranchAddress("fitstatusfit",&fitstatusfit);
   fChain->SetBranchAddress("mxhadfit1",&mxhadfit1);
   fChain->SetBranchAddress("q2fit1",&q2fit1);
   fChain->SetBranchAddress("pplusfit1",&pplusfit1);
   fChain->SetBranchAddress("chisqfit1",&chisqfit1);
   fChain->SetBranchAddress("fitstatusfit1",&fitstatusfit1);
   fChain->SetBranchAddress("mxhadfit2",&mxhadfit2);
   fChain->SetBranchAddress("q2fit2",&q2fit2);
   fChain->SetBranchAddress("pplusfit2",&pplusfit2);
   fChain->SetBranchAddress("chisqfit2",&chisqfit2);
   fChain->SetBranchAddress("fitstatusfit2",&fitstatusfit2);
   fChain->SetBranchAddress("mxhadfit3",&mxhadfit3);
   fChain->SetBranchAddress("q2fit3",&q2fit3);
   fChain->SetBranchAddress("pplusfit3",&pplusfit3);
   fChain->SetBranchAddress("chisqfit3",&chisqfit3);
   fChain->SetBranchAddress("fitstatusfit3",&fitstatusfit3);
   fChain->SetBranchAddress("mxhadfit4",&mxhadfit4);
   fChain->SetBranchAddress("q2fit4",&q2fit4);
   fChain->SetBranchAddress("pplusfit4",&pplusfit4);
   fChain->SetBranchAddress("chisqfit4",&chisqfit4);
   fChain->SetBranchAddress("fitstatusfit4",&fitstatusfit4);
   fChain->SetBranchAddress("mxhadfit5",&mxhadfit5);
   fChain->SetBranchAddress("q2fit5",&q2fit5);
   fChain->SetBranchAddress("pplusfit5",&pplusfit5);
   fChain->SetBranchAddress("chisqfit5",&chisqfit5);
   fChain->SetBranchAddress("fitstatusfit5",&fitstatusfit5);
   fChain->SetBranchAddress("mxhadfit6",&mxhadfit6);
   fChain->SetBranchAddress("q2fit6",&q2fit6);
   fChain->SetBranchAddress("pplusfit6",&pplusfit6);
   fChain->SetBranchAddress("chisqfit6",&chisqfit6);
   fChain->SetBranchAddress("fitstatusfit6",&fitstatusfit6);
   fChain->SetBranchAddress("mxhadfit7",&mxhadfit7);
   fChain->SetBranchAddress("q2fit7",&q2fit7);
   fChain->SetBranchAddress("pplusfit7",&pplusfit7);
   fChain->SetBranchAddress("chisqfit7",&chisqfit7);
   fChain->SetBranchAddress("fitstatusfit7",&fitstatusfit7);
   fChain->SetBranchAddress("mxhadfit8",&mxhadfit8);
   fChain->SetBranchAddress("q2fit8",&q2fit8);
   fChain->SetBranchAddress("pplusfit8",&pplusfit8);
   fChain->SetBranchAddress("chisqfit8",&chisqfit8);
   fChain->SetBranchAddress("fitstatusfit8",&fitstatusfit8);
   Notify();
}

Bool_t fittest::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. Typically here the branch pointers
   // will be retrieved. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed.

   // Get branch pointers
   b_run = fChain->GetBranch("run");
   b_upper = fChain->GetBranch("upper");
   b_lower = fChain->GetBranch("lower");
   b_Gvxbtyp = fChain->GetBranch("Gvxbtyp");
   b_GSem = fChain->GetBranch("GSem");
   b_GfDpi = fChain->GetBranch("GfDpi");
   b_GfDpiz = fChain->GetBranch("GfDpiz");
   b_GfDk = fChain->GetBranch("GfDk");
   b_GfDks = fChain->GetBranch("GfDks");
   b_GfDkl = fChain->GetBranch("GfDkl");
   b_GfDlep = fChain->GetBranch("GfDlep");
   b_GfDgam = fChain->GetBranch("GfDgam");
   b_GfDnu = fChain->GetBranch("GfDnu");
   b_GfD0Ds = fChain->GetBranch("GfD0Ds");
   b_GfDDs = fChain->GetBranch("GfDDs");
   b_GfDkspiopio = fChain->GetBranch("GfDkspiopio");
   b_GfK = fChain->GetBranch("GfK");
   b_isassocB = fChain->GetBranch("isassocB");
   b_isassocB_GHIT = fChain->GetBranch("isassocB_GHIT");
   b_ass_deltapB = fChain->GetBranch("ass_deltapB");
   b_isGoodMatch = fChain->GetBranch("isGoodMatch");
   b_ch1B = fChain->GetBranch("ch1B");
   b_ch2B = fChain->GetBranch("ch2B");
   b_chunm = fChain->GetBranch("chunm");
   b_neu1B = fChain->GetBranch("neu1B");
   b_neu2B = fChain->GetBranch("neu2B");
   b_neuunm = fChain->GetBranch("neuunm");
   b_brecoqual = fChain->GetBranch("brecoqual");
   b_brqual = fChain->GetBranch("brqual");
   b_brecoqualangle = fChain->GetBranch("brecoqualangle");
   b_chgdaugen = fChain->GetBranch("chgdaugen");
   b_neudaugen = fChain->GetBranch("neudaugen");
   b_nchg = fChain->GetBranch("nchg");
   b_nneu = fChain->GetBranch("nneu");
   b_xcharge = fChain->GetBranch("xcharge");
   b_nB = fChain->GetBranch("nB");
   b_brecoid = fChain->GetBranch("brecoid");
   b_brecoidtrue = fChain->GetBranch("brecoidtrue");
   b_brecoflav = fChain->GetBranch("brecoflav");
   b_brecocharge = fChain->GetBranch("brecocharge");
   b_modeB = fChain->GetBranch("modeB");
   b_truemodeB = fChain->GetBranch("truemodeB");
   b_isdoubleD = fChain->GetBranch("isdoubleD");
   b_mes = fChain->GetBranch("mes");
   b_mesendpoint = fChain->GetBranch("mesendpoint");
   b_de = fChain->GetBranch("de");
   b_pB = fChain->GetBranch("pB");
   b_eB = fChain->GetBranch("eB");
   b_eUps = fChain->GetBranch("eUps");
   b_pUps = fChain->GetBranch("pUps");
   b_thetaUps = fChain->GetBranch("thetaUps");
   b_phiUps = fChain->GetBranch("phiUps");
   b_thetaB = fChain->GetBranch("thetaB");
   b_phiB = fChain->GetBranch("phiB");
   b_pBtrue = fChain->GetBranch("pBtrue");
   b_eBtrue = fChain->GetBranch("eBtrue");
   b_tBtrue = fChain->GetBranch("tBtrue");
   b_fBtrue = fChain->GetBranch("fBtrue");
   b_pur = fChain->GetBranch("pur");
   b_intpur = fChain->GetBranch("intpur");
   b_nle_nopcut = fChain->GetBranch("nle_nopcut");
   b_nle = fChain->GetBranch("nle");
   b_nel = fChain->GetBranch("nel");
   b_nmu = fChain->GetBranch("nmu");
   b_nkp = fChain->GetBranch("nkp");
   b_nks = fChain->GetBranch("nks");
   b_nlept500 = fChain->GetBranch("nlept500");
   b_nelec500 = fChain->GetBranch("nelec500");
   b_nmu500 = fChain->GetBranch("nmu500");
   b_nlept1000 = fChain->GetBranch("nlept1000");
   b_nelec1000 = fChain->GetBranch("nelec1000");
   b_nmu1000 = fChain->GetBranch("nmu1000");
   b_deltam = fChain->GetBranch("deltam");
   b_MM1pr = fChain->GetBranch("MM1pr");
   b_MM2pr = fChain->GetBranch("MM2pr");
   b_MM3pr = fChain->GetBranch("MM3pr");
   b_OA1 = fChain->GetBranch("OA1");
   b_OA2 = fChain->GetBranch("OA2");
   b_OA3 = fChain->GetBranch("OA3");
   b_PiMin1 = fChain->GetBranch("PiMin1");
   b_PiMin2 = fChain->GetBranch("PiMin2");
   b_PiMin3 = fChain->GetBranch("PiMin3");
   b_plab = fChain->GetBranch("plab");
   b_elab = fChain->GetBranch("elab");
   b_tlab = fChain->GetBranch("tlab");
   b_flab = fChain->GetBranch("flab");
   b_plabgen = fChain->GetBranch("plabgen");
   b_elabgen = fChain->GetBranch("elabgen");
   b_tlabgen = fChain->GetBranch("tlabgen");
   b_flabgen = fChain->GetBranch("flabgen");
   b_lchargegen = fChain->GetBranch("lchargegen");
   b_pcms = fChain->GetBranch("pcms");
   b_ecms = fChain->GetBranch("ecms");
   b_tcms = fChain->GetBranch("tcms");
   b_fcms = fChain->GetBranch("fcms");
   b_lcharge = fChain->GetBranch("lcharge");
   b_leptidgen = fChain->GetBranch("leptidgen");
   b_leptorg = fChain->GetBranch("leptorg");
   b_isele = fChain->GetBranch("isele");
   b_vub = fChain->GetBranch("vub");
   b_vcb = fChain->GetBranch("vcb");
   b_other = fChain->GetBranch("other");
   b_nvubexcl = fChain->GetBranch("nvubexcl");
   b_nvubnres = fChain->GetBranch("nvubnres");
   b_ntau = fChain->GetBranch("ntau");
   b_mxhadgen = fChain->GetBranch("mxhadgen");
   b_mxhadgenwoph = fChain->GetBranch("mxhadgenwoph");
   b_xchargegen = fChain->GetBranch("xchargegen");
   b_pcmsgen = fChain->GetBranch("pcmsgen");
   b_ecmsgen = fChain->GetBranch("ecmsgen");
   b_tcmsgen = fChain->GetBranch("tcmsgen");
   b_fcmsgen = fChain->GetBranch("fcmsgen");
   b_pxhadgen = fChain->GetBranch("pxhadgen");
   b_exhadgen = fChain->GetBranch("exhadgen");
   b_exhadgencms = fChain->GetBranch("exhadgencms");
   b_txhadgen = fChain->GetBranch("txhadgen");
   b_fxhadgen = fChain->GetBranch("fxhadgen");
   b_q2Gen = fChain->GetBranch("q2Gen");
   b_ctvgen = fChain->GetBranch("ctvgen");
   b_ctlgen = fChain->GetBranch("ctlgen");
   b_chigen = fChain->GetBranch("chigen");
   b_enugencms = fChain->GetBranch("enugencms");
   b_pnugencms = fChain->GetBranch("pnugencms");
   b_mxhad = fChain->GetBranch("mxhad");
   b_emiss = fChain->GetBranch("emiss");
   b_pmiss = fChain->GetBranch("pmiss");
   b_mm2 = fChain->GetBranch("mm2");
   b_q2 = fChain->GetBranch("q2");
   b_q2new = fChain->GetBranch("q2new");
   b_exhadcms = fChain->GetBranch("exhadcms");
   b_pxhad = fChain->GetBranch("pxhad");
   b_exhad = fChain->GetBranch("exhad");
   b_txhad = fChain->GetBranch("txhad");
   b_fxhad = fChain->GetBranch("fxhad");
   b_pnu = fChain->GetBranch("pnu");
   b_tnu = fChain->GetBranch("tnu");
   b_fnu = fChain->GetBranch("fnu");
   b_pplus = fChain->GetBranch("pplus");
   b_pminus = fChain->GetBranch("pminus");
   b_pplusgen = fChain->GetBranch("pplusgen");
   b_pminusgen = fChain->GetBranch("pminusgen");
   b_pplusfit = fChain->GetBranch("pplusfit");
   b_pminusfit = fChain->GetBranch("pminusfit");
   b_wdeltam = fChain->GetBranch("wdeltam");
   b_mxhadfit = fChain->GetBranch("mxhadfit");
   b_mm2fit = fChain->GetBranch("mm2fit");
   b_q2fit = fChain->GetBranch("q2fit");
   b_chisq = fChain->GetBranch("chisq");
   b_globchisq = fChain->GetBranch("globchisq");
   b_probchisq = fChain->GetBranch("probchisq");
   b_ndof = fChain->GetBranch("ndof");
   b_fitstatus = fChain->GetBranch("fitstatus");
   b_nnpi0 = fChain->GetBranch("nnpi0");
   b_nneu80_160 = fChain->GetBranch("nneu80_160");
   b_nneu160_320 = fChain->GetBranch("nneu160_320");
   b_pstarfitlept = fChain->GetBranch("pstarfitlept");
   b_pfitX = fChain->GetBranch("pfitX");
   b_thetafitX = fChain->GetBranch("thetafitX");
   b_phifitX = fChain->GetBranch("phifitX");
   b_pfitlept = fChain->GetBranch("pfitlept");
   b_thetafitlept = fChain->GetBranch("thetafitlept");
   b_phifitlept = fChain->GetBranch("phifitlept");
   b_pfitB = fChain->GetBranch("pfitB");
   b_thetafitB = fChain->GetBranch("thetafitB");
   b_phifitB = fChain->GetBranch("phifitB");
   b_totweight = fChain->GetBranch("totweight");
   b_totweightNutMult = fChain->GetBranch("totweightNutMult");
   b_totweightTrkMult = fChain->GetBranch("totweightTrkMult");
   b_kplus = fChain->GetBranch("kplus");
   b_nBrems = fChain->GetBranch("nBrems");
   b_eBrems = fChain->GetBranch("eBrems");
   b_mxhadfit0 = fChain->GetBranch("mxhadfit0");
   b_q2fit0 = fChain->GetBranch("q2fit0");
   b_pplusfit0 = fChain->GetBranch("pplusfit0");
   b_chisqfit0 = fChain->GetBranch("chisqfit0");
   b_fitstatusfit0 = fChain->GetBranch("fitstatusfit0");
   b_chisqfit = fChain->GetBranch("chisqfit");
   b_fitstatusfit = fChain->GetBranch("fitstatusfit");
   b_mxhadfit1 = fChain->GetBranch("mxhadfit1");
   b_q2fit1 = fChain->GetBranch("q2fit1");
   b_pplusfit1 = fChain->GetBranch("pplusfit1");
   b_chisqfit1 = fChain->GetBranch("chisqfit1");
   b_fitstatusfit1 = fChain->GetBranch("fitstatusfit1");
   b_mxhadfit2 = fChain->GetBranch("mxhadfit2");
   b_q2fit2 = fChain->GetBranch("q2fit2");
   b_pplusfit2 = fChain->GetBranch("pplusfit2");
   b_chisqfit2 = fChain->GetBranch("chisqfit2");
   b_fitstatusfit2 = fChain->GetBranch("fitstatusfit2");
   b_mxhadfit3 = fChain->GetBranch("mxhadfit3");
   b_q2fit3 = fChain->GetBranch("q2fit3");
   b_pplusfit3 = fChain->GetBranch("pplusfit3");
   b_chisqfit3 = fChain->GetBranch("chisqfit3");
   b_fitstatusfit3 = fChain->GetBranch("fitstatusfit3");
   b_mxhadfit4 = fChain->GetBranch("mxhadfit4");
   b_q2fit4 = fChain->GetBranch("q2fit4");
   b_pplusfit4 = fChain->GetBranch("pplusfit4");
   b_chisqfit4 = fChain->GetBranch("chisqfit4");
   b_fitstatusfit4 = fChain->GetBranch("fitstatusfit4");
   b_mxhadfit5 = fChain->GetBranch("mxhadfit5");
   b_q2fit5 = fChain->GetBranch("q2fit5");
   b_pplusfit5 = fChain->GetBranch("pplusfit5");
   b_chisqfit5 = fChain->GetBranch("chisqfit5");
   b_fitstatusfit5 = fChain->GetBranch("fitstatusfit5");
   b_mxhadfit6 = fChain->GetBranch("mxhadfit6");
   b_q2fit6 = fChain->GetBranch("q2fit6");
   b_pplusfit6 = fChain->GetBranch("pplusfit6");
   b_chisqfit6 = fChain->GetBranch("chisqfit6");
   b_fitstatusfit6 = fChain->GetBranch("fitstatusfit6");
   b_mxhadfit7 = fChain->GetBranch("mxhadfit7");
   b_q2fit7 = fChain->GetBranch("q2fit7");
   b_pplusfit7 = fChain->GetBranch("pplusfit7");
   b_chisqfit7 = fChain->GetBranch("chisqfit7");
   b_fitstatusfit7 = fChain->GetBranch("fitstatusfit7");
   b_mxhadfit8 = fChain->GetBranch("mxhadfit8");
   b_q2fit8 = fChain->GetBranch("q2fit8");
   b_pplusfit8 = fChain->GetBranch("pplusfit8");
   b_chisqfit8 = fChain->GetBranch("chisqfit8");
   b_fitstatusfit8 = fChain->GetBranch("fitstatusfit8");

   return kTRUE;
}

void fittest::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t fittest::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef fittest_cxx
