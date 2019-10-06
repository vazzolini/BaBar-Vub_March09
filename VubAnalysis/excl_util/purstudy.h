//////////////////////////////////////////////////////////
//   This class has been automatically generated 
//     (Sat Jan  3 02:21:46 2004 by ROOT version3.01/06)
//   from TTree events/events
//   found on file: /u/ec/daniele/scra/excl-prod/data.root
//////////////////////////////////////////////////////////


#ifndef purstudy_h
#define purstudy_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <iostream.h>

class purstudy {
   public :
   TTree          *fChain;   //pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //current Tree number in a TChain
//Declaration of leaves types
   Int_t           run;
   Int_t           lower;
   Int_t           upper;
   Double_t        evtw8;
   Double_t        lw8;
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
   Int_t           GfDkmiss;
   Int_t           GfDlep;
   Int_t           GfDgam;
   Int_t           GfDnu;
   Int_t           GfD0Ds;
   Int_t           GfDDs;
   Int_t           GfDkl;
   Int_t           GfDkspiopio;
   Int_t           GfDkspipi;
   Double_t        intpur;
   Int_t           brecoflav;
   Int_t           brecocharge;
   Int_t           brecomc;
   Int_t           mode;
   Int_t           nnpi0;
   Int_t           nnks;
   Int_t           nnpar;
   Int_t           fBchgen;
   Double_t        pcmsgen;
   Double_t        tcmsgen;
   Double_t        fcmsgen;
   Double_t        ecmsgen;
   Double_t        mxhadgen;
   Double_t        pxhadgen;
   Double_t        txhadgen;
   Double_t        fxhadgen;
   Double_t        exhadgen;
   Double_t        GxCiuc;
   Double_t        GwCiuc;
   Double_t        GcsiCiuc;
   Double_t        xCiuc;
   Double_t        wCiuc;
   Double_t        csiCiuc;
   Double_t        EwPwfit;
   Double_t        EwPw;
   Double_t        EwPwG;
   Double_t        fkplus;
   UChar_t         GoodEvent;
   Int_t           vub;
   Int_t           vcb;
   Int_t           vxbtyp;
   Int_t           other;
   Int_t           wKK;
   Int_t           xcharge;
   Double_t        qtot;
   Double_t        pxhad;
   Double_t        txhad;
   Double_t        fxhad;
   Double_t        exhad;
   Double_t        mxhad;
   Double_t        gmax;
   Double_t        mxhadfit;
   Double_t        mxminpi;
   Int_t           bestklind;
   Double_t        deltam;
   Double_t        mxwminpi;
   Double_t        wdeltam;
   Double_t        mxphotmin;
   Double_t        photdeltam;
   Double_t        minlat;
   Double_t        mm2pi;
   Double_t        nrecopi0;
   Double_t        mpi0;
   Double_t        mm2pi0;
   Double_t        mm2gamma;
   Double_t        truemom1phpi0;
   Double_t        truemom2phpi0;
   Double_t        truemomlab1phpi0;
   Double_t        truemomlab2phpi0;
   Double_t        trueth1phpi0;
   Double_t        trueth2phpi0;
   Double_t        mom1phpi0;
   Double_t        mom2phpi0;
   Double_t        truemrho;
   Double_t        nrecorho;
   Double_t        mrho;
   Double_t        mm2rho;
   Double_t        truemompirho;
   Double_t        truemompi0rho;
   Double_t        mompirho;
   Double_t        mompi0rho;
   Double_t        truemrho0;
   Double_t        nrecorho0;
   Double_t        mrho0;
   Double_t        mm2rho0;
   Double_t        truemom1pirho0;
   Double_t        truemom2pirho0;
   Double_t        mom1pirho0;
   Double_t        mom2pirho0;
   Double_t        truemomega;
   Double_t        nrecoomega;
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
   Double_t        fPxhadchg;
   Double_t        fTxhadchg;
   Double_t        fFxhadchg;
   Double_t        fExhadchg;
   Double_t        fMxhadchg;
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
   Int_t           npi;
   Int_t           npi0;
   Int_t           nneu;
   Int_t           nneu80_160;
   Int_t           nneu160_320;
   Int_t           nneufromB;
   Int_t           nneufromB80_160;
   Int_t           nneufromB160_320;
   Int_t           nkp;
   Int_t           nks;
   Double_t        totweight;
   Double_t        totweightfRecoilNutMult;
   Double_t        totweightfRecoilTrkMult;
   Double_t        enu;
   Double_t        pnu;
   Double_t        tnu;
   Double_t        fnu;
   Double_t        mm2;
   Double_t        mm2nc;
   Double_t        mm2fit;
   Double_t        Eneu;
   Double_t        EPiz;
   Double_t        MinKMom;
   Double_t        MaxKMom;
   Double_t        q2;
   Double_t        q2Gen;
   Double_t        q2nc;
   Double_t        q2fit;
   //eta stuff
  Double_t        metagg;
  Double_t        mm2etagg;
  Double_t        metapppi0;
  Double_t        mm2etapppi0;
  Double_t        metapi0pi0pi0;
  Double_t        mm2etapi0pi0pi0;
  Double_t        meta;
  Double_t        mm2eta;  
  Int_t           modeeta;  

   //eta' stuff
  Double_t        metap;
  Double_t        mm2etap;  
  Int_t           modeetap;  
  Double_t        momrho0ph;  
  Double_t        etapmassetadau;  
  
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
   TBranch        *b_evtw8;
   TBranch        *b_lw8;
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
   TBranch        *b_GfDkmiss;
   TBranch        *b_GfDlep;
   TBranch        *b_GfDgam;
   TBranch        *b_GfDnu;
   TBranch        *b_GfD0Ds;
   TBranch        *b_GfDDs;
   TBranch        *b_GfDkl;
   TBranch        *b_GfDkspiopio;
   TBranch        *b_GfDkspipi;
   TBranch        *b_intpur;
   TBranch        *b_brecoflav;
   TBranch        *b_brecocharge;
   TBranch        *b_brecomc;
   TBranch        *b_mode;
   TBranch        *b_nnpi0;
   TBranch        *b_nnks;
   TBranch        *b_nnpar;
   TBranch        *b_fBchgen;
   TBranch        *b_pcmsgen;
   TBranch        *b_tcmsgen;
   TBranch        *b_fcmsgen;
   TBranch        *b_ecmsgen;
   TBranch        *b_mxhadgen;
   TBranch        *b_pxhadgen;
   TBranch        *b_txhadgen;
   TBranch        *b_fxhadgen;
   TBranch        *b_exhadgen;
   TBranch        *b_GxCiuc;
   TBranch        *b_GwCiuc;
   TBranch        *b_GcsiCiuc;
   TBranch        *b_xCiuc;
   TBranch        *b_wCiuc;
   TBranch        *b_csiCiuc;
   TBranch        *b_EwPwfit;
   TBranch        *b_EwPw;
   TBranch        *b_EwPwG;
   TBranch        *b_kplus;
   TBranch        *b_GoodEvent;
   TBranch        *b_vub;
   TBranch        *b_vcb;
   TBranch        *b_vxbtyp;
   TBranch        *b_other;
   TBranch        *b_wKK;
   TBranch        *b_xcharge;
   TBranch        *b_qtot;
   TBranch        *b_pxhad;
   TBranch        *b_txhad;
   TBranch        *b_fxhad;
   TBranch        *b_exhad;
   TBranch        *b_mxhad;
   TBranch        *b_gmax;
   TBranch        *b_mxhadfit;
   TBranch        *b_mxminpi;
   TBranch        *b_bestklind;
   TBranch        *b_deltam;
   TBranch        *b_mxwminpi;
   TBranch        *b_wdeltam;
   TBranch        *b_mxphotmin;
   TBranch        *b_photdeltam;
   TBranch        *b_minlat;
   TBranch        *b_mm2pi;
   TBranch        *b_nrecopi0;
   TBranch        *b_mpi0;
   TBranch        *b_mm2pi0;
   TBranch        *b_mm2gamma;
   TBranch        *b_truemom1phpi0;
   TBranch        *b_truemom2phpi0;
   TBranch        *b_truemomlab1phpi0;
   TBranch        *b_truemomlab2phpi0;
   TBranch        *b_trueth1phpi0;
   TBranch        *b_trueth2phpi0;
   TBranch        *b_mom1phpi0;
   TBranch        *b_mom2phpi0;
   TBranch        *b_truemrho;
   TBranch        *b_nrecorho;
   TBranch        *b_mrho;
   TBranch        *b_mm2rho;
   TBranch        *b_truemompirho;
   TBranch        *b_truemompi0rho;
   TBranch        *b_mompirho;
   TBranch        *b_mompi0rho;
   TBranch        *b_truemrho0;
   TBranch        *b_nrecorho0;
   TBranch        *b_mrho0;
   TBranch        *b_mm2rho0;
   TBranch        *b_truemom1pirho0;
   TBranch        *b_truemom2pirho0;
   TBranch        *b_mom1pirho0;
   TBranch        *b_mom2pirho0;
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
   TBranch        *b_pxhadchg;
   TBranch        *b_txhadchg;
   TBranch        *b_fxhadchg;
   TBranch        *b_exhadchg;
   TBranch        *b_mxhadchg;
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
   TBranch        *b_npi;
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
   TBranch        *b_mm2nc;
   TBranch        *b_mm2fit;
   TBranch        *b_ENeu;
   TBranch        *b_EPiz;
   TBranch        *b_MinKMom;
   TBranch        *b_MaxKMom;
   TBranch        *b_q2;
   TBranch        *b_q2Gen;
   TBranch        *b_q2nc;
   TBranch        *b_q2fit;
   //eta stuff
  TBranch        *b_metagg;
  TBranch        *b_mm2etagg;
  TBranch        *b_metapppi0;
  TBranch        *b_mm2etapppi0;
  TBranch        *b_metapi0pi0pi0;
  TBranch        *b_mm2etapi0pi0pi0;
  TBranch        *b_meta;
  TBranch        *b_mm2eta;
  TBranch        *b_modeeta;

   //eta' stuff
  TBranch        *b_metap;
  TBranch        *b_mm2etap;
  TBranch        *b_modeetap;
  TBranch        *b_etapmassetadau;
  TBranch        *b_momrho0ph;
  
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

   TFile *fHistFile;

   purstudy(TTree *tree=0);
   ~purstudy();
   Int_t  Cut(Int_t entry);
   Int_t  GetEntry(Int_t entry);
   Int_t  LoadTree(Int_t entry);
   void   Init(TTree *tree);
   void   Loop(char *pippo, char* mode);
   void   bookHist();
   Bool_t Notify();
   void   Show(Int_t entry = -1);
/*    TChain * getchain(char * thechain);  */
};


#endif

#ifdef purstudy_cxx
purstudy::purstudy(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/u/ec/daniele/scra/excl-prod/data.root");
      if (!f) {
         f = new TFile("/u/ec/daniele/scra/excl-prod/data.root");
      }
      tree = (TTree*)gDirectory->Get("events");

   }
   Init(tree);
}


purstudy::~purstudy()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t purstudy::GetEntry(Int_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Int_t purstudy::LoadTree(Int_t entry)
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

void purstudy::Init(TTree *tree)
{
//   Set branch addresses
   if (tree == 0) return;
   fChain    = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run",&run);
   fChain->SetBranchAddress("lower",&lower);
   fChain->SetBranchAddress("upper",&upper);
   fChain->SetBranchAddress("evtw8",&evtw8);
   fChain->SetBranchAddress("lw8",&lw8);
   fChain->SetBranchAddress("bmass",&bmass);
   fChain->SetBranchAddress("bmassfit",&bmassfit);
   fChain->SetBranchAddress("sbox",&sbox);
   fChain->SetBranchAddress("mes",&mes);
   fChain->SetBranchAddress("de",&de);
   fChain->SetBranchAddress("pur",&pur);
   fChain->SetBranchAddress("Gvxbtyp",&Gvxbtyp);
   fChain->SetBranchAddress("GSem",&GSem);
   fChain->SetBranchAddress("GfDpi",&GfDpi);
   fChain->SetBranchAddress("GfDpiz",&GfDpiz);
   fChain->SetBranchAddress("GfDk",&GfDk);
   fChain->SetBranchAddress("GfDks",&GfDks);
   fChain->SetBranchAddress("GfDkmiss",&GfDkmiss);
   fChain->SetBranchAddress("GfDlep",&GfDlep);
   fChain->SetBranchAddress("GfDgam",&GfDgam);
   fChain->SetBranchAddress("GfDnu",&GfDnu);
   fChain->SetBranchAddress("GfD0Ds",&GfD0Ds);
   fChain->SetBranchAddress("GfDDs",&GfDDs);
   fChain->SetBranchAddress("GfDkl",&GfDkl);
   fChain->SetBranchAddress("GfDkspiopio",&GfDkspiopio);
   fChain->SetBranchAddress("GfDkspipi",&GfDkspipi);
   fChain->SetBranchAddress("intpur",&intpur);
   fChain->SetBranchAddress("brecoflav",&brecoflav);
   fChain->SetBranchAddress("brecocharge",&brecocharge);
   fChain->SetBranchAddress("brecomc",&brecomc);
   fChain->SetBranchAddress("mode",&mode);
   fChain->SetBranchAddress("nnpi0",&nnpi0);
   fChain->SetBranchAddress("nnks",&nnks);
   fChain->SetBranchAddress("nnpar",&nnpar);
   fChain->SetBranchAddress("fBchgen",&fBchgen);
   fChain->SetBranchAddress("pcmsgen",&pcmsgen);
   fChain->SetBranchAddress("tcmsgen",&tcmsgen);
   fChain->SetBranchAddress("fcmsgen",&fcmsgen);
   fChain->SetBranchAddress("ecmsgen",&ecmsgen);
   fChain->SetBranchAddress("mxhadgen",&mxhadgen);
   fChain->SetBranchAddress("pxhadgen",&pxhadgen);
   fChain->SetBranchAddress("txhadgen",&txhadgen);
   fChain->SetBranchAddress("fxhadgen",&fxhadgen);
   fChain->SetBranchAddress("exhadgen",&exhadgen);
   fChain->SetBranchAddress("GxCiuc",&GxCiuc);
   fChain->SetBranchAddress("GwCiuc",&GwCiuc);
   fChain->SetBranchAddress("GcsiCiuc",&GcsiCiuc);
   fChain->SetBranchAddress("xCiuc",&xCiuc);
   fChain->SetBranchAddress("wCiuc",&wCiuc);
   fChain->SetBranchAddress("csiCiuc",&csiCiuc);
   fChain->SetBranchAddress("EwPwfit",&EwPwfit);
   fChain->SetBranchAddress("EwPw",&EwPw);
   fChain->SetBranchAddress("EwPwG",&EwPwG);
   fChain->SetBranchAddress("kplus",&fkplus);
   fChain->SetBranchAddress("GoodEvent",&GoodEvent);
   fChain->SetBranchAddress("vub",&vub);
   fChain->SetBranchAddress("vcb",&vcb);
   fChain->SetBranchAddress("vxbtyp",&vxbtyp);
   fChain->SetBranchAddress("other",&other);
   fChain->SetBranchAddress("wKK",&wKK);
   fChain->SetBranchAddress("xcharge",&xcharge);
   fChain->SetBranchAddress("qtot",&qtot);
   fChain->SetBranchAddress("pxhad",&pxhad);
   fChain->SetBranchAddress("txhad",&txhad);
   fChain->SetBranchAddress("fxhad",&fxhad);
   fChain->SetBranchAddress("exhad",&exhad);
   fChain->SetBranchAddress("mxhad",&mxhad);
   fChain->SetBranchAddress("gmax",&gmax);
   fChain->SetBranchAddress("mxhadfit",&mxhadfit);
   fChain->SetBranchAddress("mxminpi",&mxminpi);
   fChain->SetBranchAddress("bestklind",&bestklind);
   fChain->SetBranchAddress("deltam",&deltam);
   fChain->SetBranchAddress("mxwminpi",&mxwminpi);
   fChain->SetBranchAddress("wdeltam",&wdeltam);
   fChain->SetBranchAddress("mxphotmin",&mxphotmin);
   fChain->SetBranchAddress("photdeltam",&photdeltam);
   fChain->SetBranchAddress("minlat",&minlat);
   fChain->SetBranchAddress("mm2pi",&mm2pi);
   fChain->SetBranchAddress("nrecopi0",&nrecopi0);
   fChain->SetBranchAddress("mpi0",&mpi0);
   fChain->SetBranchAddress("mm2pi0",&mm2pi0);
   fChain->SetBranchAddress("mm2gamma",&mm2gamma);
   fChain->SetBranchAddress("truemom1phpi0",&truemom1phpi0);
   fChain->SetBranchAddress("truemom2phpi0",&truemom2phpi0);
   fChain->SetBranchAddress("truemomlab1phpi0",&truemomlab1phpi0);
   fChain->SetBranchAddress("truemomlab2phpi0",&truemomlab2phpi0);
   fChain->SetBranchAddress("trueth1phpi0",&trueth1phpi0);
   fChain->SetBranchAddress("trueth2phpi0",&trueth2phpi0);
   fChain->SetBranchAddress("mom1phpi0",&mom1phpi0);
   fChain->SetBranchAddress("mom2phpi0",&mom2phpi0);
   fChain->SetBranchAddress("truemrho",&truemrho);
   fChain->SetBranchAddress("nrecorho",&nrecorho);
   fChain->SetBranchAddress("mrho",&mrho);
   fChain->SetBranchAddress("mm2rho",&mm2rho);
   fChain->SetBranchAddress("truemompirho",&truemompirho);
   fChain->SetBranchAddress("truemompi0rho",&truemompi0rho);
   fChain->SetBranchAddress("mompirho",&mompirho);
   fChain->SetBranchAddress("mompi0rho",&mompi0rho);
   fChain->SetBranchAddress("truemrho0",&truemrho0);
   fChain->SetBranchAddress("nrecorho0",&nrecorho0);
   fChain->SetBranchAddress("mrho0",&mrho0);
   fChain->SetBranchAddress("mm2rho0",&mm2rho0);
   fChain->SetBranchAddress("truemom1pirho0",&truemom1pirho0);
   fChain->SetBranchAddress("truemom2pirho0",&truemom2pirho0);
   fChain->SetBranchAddress("mom1pirho0",&mom1pirho0);
   fChain->SetBranchAddress("mom2pirho0",&mom2pirho0);
   fChain->SetBranchAddress("truemomega",&truemomega);
   fChain->SetBranchAddress("nrecoomega",&nrecoomega);
   fChain->SetBranchAddress("momega",&momega);
   fChain->SetBranchAddress("mm2omega",&mm2omega);
   fChain->SetBranchAddress("truemom1piome",&truemom1piome);
   fChain->SetBranchAddress("truemom2piome",&truemom2piome);
   fChain->SetBranchAddress("truemompi0ome",&truemompi0ome);
   fChain->SetBranchAddress("truedalitzpi1pi2ome",&truedalitzpi1pi2ome);
   fChain->SetBranchAddress("truedalitzpi1pi0ome",&truedalitzpi1pi0ome);
   fChain->SetBranchAddress("truecosthome",&truecosthome);
   fChain->SetBranchAddress("mom1piome",&mom1piome);
   fChain->SetBranchAddress("mom2piome",&mom2piome);
   fChain->SetBranchAddress("mompi0ome",&mompi0ome);
   fChain->SetBranchAddress("dalitzpi1pi2ome",&dalitzpi1pi2ome);
   fChain->SetBranchAddress("dalitzpi1pi0ome",&dalitzpi1pi0ome);
   fChain->SetBranchAddress("costhome",&costhome);
   fChain->SetBranchAddress("pxhadchg",&fPxhadchg);
   fChain->SetBranchAddress("txhadchg",&fTxhadchg);
   fChain->SetBranchAddress("fxhadchg",&fFxhadchg);
   fChain->SetBranchAddress("exhadchg",&fExhadchg);
   fChain->SetBranchAddress("mxhadchg",&fMxhadchg);
   fChain->SetBranchAddress("lcharge",&lcharge);
   fChain->SetBranchAddress("plab",&plab);
   fChain->SetBranchAddress("tlab",&tlab);
   fChain->SetBranchAddress("flab",&flab);
   fChain->SetBranchAddress("pcms",&pcms);
   fChain->SetBranchAddress("tcms",&tcms);
   fChain->SetBranchAddress("fcms",&fcms);
   fChain->SetBranchAddress("ecms",&ecms);
   fChain->SetBranchAddress("nle",&nle);
   fChain->SetBranchAddress("nel",&nel);
   fChain->SetBranchAddress("nmu",&nmu);
   fChain->SetBranchAddress("nchg",&nchg);
   fChain->SetBranchAddress("npi",&npi);
   fChain->SetBranchAddress("npi0",&npi0);
   fChain->SetBranchAddress("nneu",&nneu);
   fChain->SetBranchAddress("nneu80_160",&nneu80_160);
   fChain->SetBranchAddress("nneu160_320",&nneu160_320);
   fChain->SetBranchAddress("nneufromB",&nneufromB);
   fChain->SetBranchAddress("nneufromB80_160",&nneufromB80_160);
   fChain->SetBranchAddress("nneufromB160_320",&nneufromB160_320);
   fChain->SetBranchAddress("nkp",&nkp);
   fChain->SetBranchAddress("nks",&nks);
   fChain->SetBranchAddress("totweight",&totweight);
   fChain->SetBranchAddress("totweightNutMult",&totweightfRecoilNutMult);
   fChain->SetBranchAddress("totweightTrkMult",&totweightfRecoilTrkMult);
   fChain->SetBranchAddress("enu",&enu);
   fChain->SetBranchAddress("pnu",&pnu);
   fChain->SetBranchAddress("tnu",&tnu);
   fChain->SetBranchAddress("fnu",&fnu);
   fChain->SetBranchAddress("mm2",&mm2);
   fChain->SetBranchAddress("mm2nc",&mm2nc);
   fChain->SetBranchAddress("mm2fit",&mm2fit);
   fChain->SetBranchAddress("ENeu",&Eneu);
   fChain->SetBranchAddress("EPiz",&EPiz);
   fChain->SetBranchAddress("MinKMom",&MinKMom);
   fChain->SetBranchAddress("MaxKMom",&MaxKMom);
   fChain->SetBranchAddress("q2",&q2);
   fChain->SetBranchAddress("q2Gen",&q2Gen);
   fChain->SetBranchAddress("q2nc",&q2nc);
   fChain->SetBranchAddress("q2fit",&q2fit);
   fChain->SetBranchAddress("meta",&meta);
  fChain->SetBranchAddress("mm2eta",&mm2eta);
  fChain->SetBranchAddress("modeeta",&modeeta);
  fChain->SetBranchAddress("metap",&metap);
  fChain->SetBranchAddress("mm2etap",&mm2etap);
  fChain->SetBranchAddress("modeetap",&modeetap);
  fChain->SetBranchAddress("etapmassetadau",&etapmassetadau);
  fChain->SetBranchAddress("momrho0ph",&momrho0ph);
  fChain->SetBranchAddress("ma0",&ma0);
  fChain->SetBranchAddress("mm2a0",&mm2a0);
  fChain->SetBranchAddress("a0massetadau",&a0massetadau);
  fChain->SetBranchAddress("modea0",&modea0);
  fChain->SetBranchAddress("ma0p",&ma0p);
  fChain->SetBranchAddress("mm2a0p",&mm2a0p);
  fChain->SetBranchAddress("a0pmassetadau",&a0pmassetadau);
  fChain->SetBranchAddress("modea0p",&modea0p);
  Notify();
}

Bool_t purstudy::Notify()
{
//   called when loading a new file
//   get branch pointers
   b_run = fChain->GetBranch("run");
   b_lower = fChain->GetBranch("lower");
   b_upper = fChain->GetBranch("upper");
   b_evtw8 = fChain->GetBranch("evtw8");
   b_lw8 = fChain->GetBranch("lw8");
   b_bmass = fChain->GetBranch("bmass");
   b_bmassfit = fChain->GetBranch("bmassfit");
   b_sbox = fChain->GetBranch("sbox");
   b_mes = fChain->GetBranch("mes");
   b_de = fChain->GetBranch("de");
   b_pur = fChain->GetBranch("pur");
   b_Gvxbtyp = fChain->GetBranch("Gvxbtyp");
   b_GSem = fChain->GetBranch("GSem");
   b_GfDpi = fChain->GetBranch("GfDpi");
   b_GfDpiz = fChain->GetBranch("GfDpiz");
   b_GfDk = fChain->GetBranch("GfDk");
   b_GfDks = fChain->GetBranch("GfDks");
   b_GfDkmiss = fChain->GetBranch("GfDkmiss");
   b_GfDlep = fChain->GetBranch("GfDlep");
   b_GfDgam = fChain->GetBranch("GfDgam");
   b_GfDnu = fChain->GetBranch("GfDnu");
   b_GfD0Ds = fChain->GetBranch("GfD0Ds");
   b_GfDDs = fChain->GetBranch("GfDDs");
   b_GfDkl = fChain->GetBranch("GfDkl");
   b_GfDkspiopio = fChain->GetBranch("GfDkspiopio");
   b_GfDkspipi = fChain->GetBranch("GfDkspipi");
   b_intpur = fChain->GetBranch("intpur");
   b_brecoflav = fChain->GetBranch("brecoflav");
   b_brecocharge = fChain->GetBranch("brecocharge");
   b_brecomc = fChain->GetBranch("brecomc");
   b_mode = fChain->GetBranch("mode");
   b_nnpi0 = fChain->GetBranch("nnpi0");
   b_nnks = fChain->GetBranch("nnks");
   b_nnpar = fChain->GetBranch("nnpar");
   b_fBchgen = fChain->GetBranch("fBchgen");
   b_pcmsgen = fChain->GetBranch("pcmsgen");
   b_tcmsgen = fChain->GetBranch("tcmsgen");
   b_fcmsgen = fChain->GetBranch("fcmsgen");
   b_ecmsgen = fChain->GetBranch("ecmsgen");
   b_mxhadgen = fChain->GetBranch("mxhadgen");
   b_pxhadgen = fChain->GetBranch("pxhadgen");
   b_txhadgen = fChain->GetBranch("txhadgen");
   b_fxhadgen = fChain->GetBranch("fxhadgen");
   b_exhadgen = fChain->GetBranch("exhadgen");
   b_GxCiuc = fChain->GetBranch("GxCiuc");
   b_GwCiuc = fChain->GetBranch("GwCiuc");
   b_GcsiCiuc = fChain->GetBranch("GcsiCiuc");
   b_xCiuc = fChain->GetBranch("xCiuc");
   b_wCiuc = fChain->GetBranch("wCiuc");
   b_csiCiuc = fChain->GetBranch("csiCiuc");
   b_EwPwfit = fChain->GetBranch("EwPwfit");
   b_EwPw = fChain->GetBranch("EwPw");
   b_EwPwG = fChain->GetBranch("EwPwG");
   b_kplus = fChain->GetBranch("kplus");
   b_GoodEvent = fChain->GetBranch("GoodEvent");
   b_vub = fChain->GetBranch("vub");
   b_vcb = fChain->GetBranch("vcb");
   b_vxbtyp = fChain->GetBranch("vxbtyp");
   b_other = fChain->GetBranch("other");
   b_wKK = fChain->GetBranch("wKK");
   b_xcharge = fChain->GetBranch("xcharge");
   b_qtot = fChain->GetBranch("qtot");
   b_pxhad = fChain->GetBranch("pxhad");
   b_txhad = fChain->GetBranch("txhad");
   b_fxhad = fChain->GetBranch("fxhad");
   b_exhad = fChain->GetBranch("exhad");
   b_mxhad = fChain->GetBranch("mxhad");
   b_gmax = fChain->GetBranch("gmax");
   b_mxhadfit = fChain->GetBranch("mxhadfit");
   b_mxminpi = fChain->GetBranch("mxminpi");
   b_bestklind = fChain->GetBranch("bestklind");
   b_deltam = fChain->GetBranch("deltam");
   b_mxwminpi = fChain->GetBranch("mxwminpi");
   b_wdeltam = fChain->GetBranch("wdeltam");
   b_mxphotmin = fChain->GetBranch("mxphotmin");
   b_photdeltam = fChain->GetBranch("photdeltam");
   b_minlat = fChain->GetBranch("minlat");
   b_mm2pi = fChain->GetBranch("mm2pi");
   b_nrecopi0 = fChain->GetBranch("nrecopi0");
   b_mpi0 = fChain->GetBranch("mpi0");
   b_mm2pi0 = fChain->GetBranch("mm2pi0");
   b_mm2gamma = fChain->GetBranch("mm2gamma");
   b_truemom1phpi0 = fChain->GetBranch("truemom1phpi0");
   b_truemom2phpi0 = fChain->GetBranch("truemom2phpi0");
   b_truemomlab1phpi0 = fChain->GetBranch("truemomlab1phpi0");
   b_truemomlab2phpi0 = fChain->GetBranch("truemomlab2phpi0");
   b_trueth1phpi0 = fChain->GetBranch("trueth1phpi0");
   b_trueth2phpi0 = fChain->GetBranch("trueth2phpi0");
   b_mom1phpi0 = fChain->GetBranch("mom1phpi0");
   b_mom2phpi0 = fChain->GetBranch("mom2phpi0");
   b_truemrho = fChain->GetBranch("truemrho");
   b_nrecorho = fChain->GetBranch("nrecorho");
   b_mrho = fChain->GetBranch("mrho");
   b_mm2rho = fChain->GetBranch("mm2rho");
   b_truemompirho = fChain->GetBranch("truemompirho");
   b_truemompi0rho = fChain->GetBranch("truemompi0rho");
   b_mompirho = fChain->GetBranch("mompirho");
   b_mompi0rho = fChain->GetBranch("mompi0rho");
   b_truemrho0 = fChain->GetBranch("truemrho0");
   b_nrecorho0 = fChain->GetBranch("nrecorho0");
   b_mrho0 = fChain->GetBranch("mrho0");
   b_mm2rho0 = fChain->GetBranch("mm2rho0");
   b_truemom1pirho0 = fChain->GetBranch("truemom1pirho0");
   b_truemom2pirho0 = fChain->GetBranch("truemom2pirho0");
   b_mom1pirho0 = fChain->GetBranch("mom1pirho0");
   b_mom2pirho0 = fChain->GetBranch("mom2pirho0");
   b_truemomega = fChain->GetBranch("truemomega");
   b_nrecoomega = fChain->GetBranch("nrecoomega");
   b_momega = fChain->GetBranch("momega");
   b_mm2omega = fChain->GetBranch("mm2omega");
   b_truemom1piome = fChain->GetBranch("truemom1piome");
   b_truemom2piome = fChain->GetBranch("truemom2piome");
   b_truemompi0ome = fChain->GetBranch("truemompi0ome");
   b_truedalitzpi1pi2ome = fChain->GetBranch("truedalitzpi1pi2ome");
   b_truedalitzpi1pi0ome = fChain->GetBranch("truedalitzpi1pi0ome");
   b_truecosthome = fChain->GetBranch("truecosthome");
   b_mom1piome = fChain->GetBranch("mom1piome");
   b_mom2piome = fChain->GetBranch("mom2piome");
   b_mompi0ome = fChain->GetBranch("mompi0ome");
   b_dalitzpi1pi2ome = fChain->GetBranch("dalitzpi1pi2ome");
   b_dalitzpi1pi0ome = fChain->GetBranch("dalitzpi1pi0ome");
   b_costhome = fChain->GetBranch("costhome");
   b_pxhadchg = fChain->GetBranch("pxhadchg");
   b_txhadchg = fChain->GetBranch("txhadchg");
   b_fxhadchg = fChain->GetBranch("fxhadchg");
   b_exhadchg = fChain->GetBranch("exhadchg");
   b_mxhadchg = fChain->GetBranch("mxhadchg");
   b_lcharge = fChain->GetBranch("lcharge");
   b_plab = fChain->GetBranch("plab");
   b_tlab = fChain->GetBranch("tlab");
   b_flab = fChain->GetBranch("flab");
   b_pcms = fChain->GetBranch("pcms");
   b_tcms = fChain->GetBranch("tcms");
   b_fcms = fChain->GetBranch("fcms");
   b_ecms = fChain->GetBranch("ecms");
   b_nle = fChain->GetBranch("nle");
   b_nel = fChain->GetBranch("nel");
   b_nmu = fChain->GetBranch("nmu");
   b_nchg = fChain->GetBranch("nchg");
   b_npi = fChain->GetBranch("npi");
   b_npi0 = fChain->GetBranch("npi0");
   b_nneu = fChain->GetBranch("nneu");
   b_nneu80_160 = fChain->GetBranch("nneu80_160");
   b_nneu160_320 = fChain->GetBranch("nneu160_320");
   b_nneufromB = fChain->GetBranch("nneufromB");
   b_nneufromB80_160 = fChain->GetBranch("nneufromB80_160");
   b_nneufromB160_320 = fChain->GetBranch("nneufromB160_320");
   b_nkp = fChain->GetBranch("nkp");
   b_nks = fChain->GetBranch("nks");
   b_totweight = fChain->GetBranch("totweight");
   b_totweightNutMult = fChain->GetBranch("totweightNutMult");
   b_totweightTrkMult = fChain->GetBranch("totweightTrkMult");
   b_enu = fChain->GetBranch("enu");
   b_pnu = fChain->GetBranch("pnu");
   b_tnu = fChain->GetBranch("tnu");
   b_fnu = fChain->GetBranch("fnu");
   b_mm2 = fChain->GetBranch("mm2");
   b_mm2nc = fChain->GetBranch("mm2nc");
   b_mm2fit = fChain->GetBranch("mm2fit");
   b_ENeu = fChain->GetBranch("ENeu");
   b_EPiz = fChain->GetBranch("EPiz");
   b_MinKMom = fChain->GetBranch("MinKMom");
   b_MaxKMom = fChain->GetBranch("MaxKMom");
   b_q2 = fChain->GetBranch("q2");
   b_q2Gen = fChain->GetBranch("q2Gen");
   b_q2nc = fChain->GetBranch("q2nc");
   b_q2fit = fChain->GetBranch("q2fit");
  b_meta = fChain->GetBranch("meta");
  b_mm2eta = fChain->GetBranch("mm2eta");
  b_metap = fChain->GetBranch("metap");
  b_mm2etap = fChain->GetBranch("mm2etap");
  b_ma0 = fChain->GetBranch("ma0");
  b_mm2a0 = fChain->GetBranch("mm2a0");
  b_ma0p = fChain->GetBranch("ma0p");
  b_mm2a0p = fChain->GetBranch("mm2a0p");
   return kTRUE;
}

void purstudy::Show(Int_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t purstudy::Cut(Int_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef purstudy_cxx

