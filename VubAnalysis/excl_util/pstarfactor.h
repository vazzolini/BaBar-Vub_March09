//////////////////////////////////////////////////////////
//   This class has been automatically generated 
//     (Thu Feb 13 17:41:22 2003 by ROOT version3.02/07)
//   from TTree events/events
//   found on file: /u/ec/daniele/scra/allcock.root
//////////////////////////////////////////////////////////


#ifndef pstarfactor_h
#define pstarfactor_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class pstarfactor {
   public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
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
   UChar_t         GoodEvent;
   Int_t           vub;
   Int_t           vcb;
   Int_t           vxbtyp;
   Int_t           other;
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
   Double_t        deltam;
   Double_t        mxwminpi;
   Double_t        wdeltam;
   Double_t        mxphotmin;
   Double_t        photdeltam;
   Double_t        minlat;
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
   Double_t        q2nc;
   Double_t        q2fit;

//List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lower;   //!
   TBranch        *b_upper;   //!
   TBranch        *b_evtw8;   //!
   TBranch        *b_lw8;   //!
   TBranch        *b_bmass;   //!
   TBranch        *b_bmassfit;   //!
   TBranch        *b_sbox;   //!
   TBranch        *b_mes;   //!
   TBranch        *b_de;   //!
   TBranch        *b_pur;   //!
   TBranch        *b_Gvxbtyp;   //!
   TBranch        *b_GSem;   //!
   TBranch        *b_GfDpi;   //!
   TBranch        *b_GfDpiz;   //!
   TBranch        *b_GfDk;   //!
   TBranch        *b_GfDks;   //!
   TBranch        *b_GfDkmiss;   //!
   TBranch        *b_GfDlep;   //!
   TBranch        *b_GfDgam;   //!
   TBranch        *b_GfD0Ds;   //!
   TBranch        *b_GfDDs;   //!
   TBranch        *b_GfDkl;   //!
   TBranch        *b_GfDkspiopio;   //!
   TBranch        *b_GfDkspipi;   //!
   TBranch        *b_intpur;   //!
   TBranch        *b_brecoflav;   //!
   TBranch        *b_brecocharge;   //!
   TBranch        *b_brecomc;   //!
   TBranch        *b_mode;   //!
   TBranch        *b_nnpi0;   //!
   TBranch        *b_nnks;   //!
   TBranch        *b_nnpar;   //!
   TBranch        *b_fBchgen;   //!
   TBranch        *b_mxhadgen;   //!
   TBranch        *b_pcmsgen;   //!
   TBranch        *b_tcmsgen;   //!
   TBranch        *b_fcmsgen;   //!
   TBranch        *b_ecmsgen;   //!
   TBranch        *b_pxhadgen;   //!
   TBranch        *b_txhadgen;   //!
   TBranch        *b_fxhadgen;   //!
   TBranch        *b_exhadgen;   //!
   TBranch        *b_fkplus;   //!
   TBranch        *b_GoodEvent;   //!
   TBranch        *b_vub;   //!
   TBranch        *b_vcb;   //!
   TBranch        *b_vxbtyp;   //!
   TBranch        *b_other;   //!
   TBranch        *b_xcharge;   //!
   TBranch        *b_qtot;   //!
   TBranch        *b_pxhad;   //!
   TBranch        *b_txhad;   //!
   TBranch        *b_fxhad;   //!
   TBranch        *b_exhad;   //!
   TBranch        *b_mxhad;   //!
   TBranch        *b_gmax;   //!
   TBranch        *b_mxhadfit;   //!
   TBranch        *b_mxminpi;   //!
   TBranch        *b_deltam;   //!
   TBranch        *b_mxwminpi;   //!
   TBranch        *b_wdeltam;   //!
   TBranch        *b_mxphotmin;   //!
   TBranch        *b_photdeltam;   //!
   TBranch        *b_minlat;   //!
   TBranch        *b_fPxhadchg;   //!
   TBranch        *b_fTxhadchg;   //!
   TBranch        *b_fFxhadchg;   //!
   TBranch        *b_fExhadchg;   //!
   TBranch        *b_fMxhadchg;   //!
   TBranch        *b_lcharge;   //!
   TBranch        *b_plab;   //!
   TBranch        *b_tlab;   //!
   TBranch        *b_flab;   //!
   TBranch        *b_pcms;   //!
   TBranch        *b_tcms;   //!
   TBranch        *b_fcms;   //!
   TBranch        *b_ecms;   //!
   TBranch        *b_nle;   //!
   TBranch        *b_nel;   //!
   TBranch        *b_nmu;   //!
   TBranch        *b_nchg;   //!
   TBranch        *b_npi0;   //!
   TBranch        *b_nneu;   //!
   TBranch        *b_nneu80_160;   //!
   TBranch        *b_nneu160_320;   //!
   TBranch        *b_nneufromB;   //!
   TBranch        *b_nneufromB80_160;   //!
   TBranch        *b_nneufromB160_320;   //!
   TBranch        *b_nkp;   //!
   TBranch        *b_nks;   //!
   TBranch        *b_totweight;   //!
   TBranch        *b_totweightfRecoilNutMult;   //!
   TBranch        *b_totweightfRecoilTrkMult;   //!
   TBranch        *b_enu;   //!
   TBranch        *b_pnu;   //!
   TBranch        *b_tnu;   //!
   TBranch        *b_fnu;   //!
   TBranch        *b_mm2;   //!
   TBranch        *b_mm2nc;   //!
   TBranch        *b_mm2fit;   //!
   TBranch        *b_Eneu;   //!
   TBranch        *b_EPiz;   //!
   TBranch        *b_MinKMom;   //!
   TBranch        *b_MaxKMom;   //!
   TBranch        *b_q2;   //!
   TBranch        *b_q2nc;   //!
   TBranch        *b_q2fit;   //!

   double nvcb, nvub;
   double nvuberr, nvcberr;
   double nvcblep[16],nvublep[16];
   double nvcbleperr[16],nvubleperr[16];

   int lepton;

   pstarfactor(TTree *tree=0);
   ~pstarfactor();
   TFile *fHistFile;
   Int_t  Cut(Int_t entry);
   Int_t  GetEntry(Int_t entry);
   Int_t  LoadTree(Int_t entry);
   void   Init(TTree *tree);
   void   Loop(int leptype, int imode, int isvub);
   void   Bookhist(int imode);
   void   FitMes();
   void   Finalize(int imode);
   Bool_t Notify();
   void   Show(Int_t entry = -1);
};

#endif

#ifdef pstarfactor_cxx
pstarfactor::pstarfactor(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/u/ec/daniele/scra/allcock.root");
      if (!f) {
         f = new TFile("/u/ec/daniele/scra/allcock.root");
      }
      tree = (TTree*)gDirectory->Get("events");

   }
   Init(tree);
}

pstarfactor::~pstarfactor()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t pstarfactor::GetEntry(Int_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Int_t pstarfactor::LoadTree(Int_t entry)
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

void pstarfactor::Init(TTree *tree)
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
   fChain->SetBranchAddress("mxhadgen",&mxhadgen);
   fChain->SetBranchAddress("pcmsgen",&pcmsgen);
   fChain->SetBranchAddress("tcmsgen",&tcmsgen);
   fChain->SetBranchAddress("fcmsgen",&fcmsgen);
   fChain->SetBranchAddress("ecmsgen",&ecmsgen);
   fChain->SetBranchAddress("pxhadgen",&pxhadgen);
   fChain->SetBranchAddress("txhadgen",&txhadgen);
   fChain->SetBranchAddress("fxhadgen",&fxhadgen);
   fChain->SetBranchAddress("exhadgen",&exhadgen);
   fChain->SetBranchAddress("kplus",&fkplus);
   fChain->SetBranchAddress("GoodEvent",&GoodEvent);
   fChain->SetBranchAddress("vub",&vub);
   fChain->SetBranchAddress("vcb",&vcb);
   fChain->SetBranchAddress("vxbtyp",&vxbtyp);
   fChain->SetBranchAddress("other",&other);
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
   fChain->SetBranchAddress("deltam",&deltam);
   fChain->SetBranchAddress("mxwminpi",&mxwminpi);
   fChain->SetBranchAddress("wdeltam",&wdeltam);
   fChain->SetBranchAddress("mxphotmin",&mxphotmin);
   fChain->SetBranchAddress("photdeltam",&photdeltam);
   fChain->SetBranchAddress("minlat",&minlat);
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
   fChain->SetBranchAddress("q2nc",&q2nc);
   fChain->SetBranchAddress("q2fit",&q2fit);
   Notify();
}

Bool_t pstarfactor::Notify()
{
   // Called when loading a new file.
   // Get branch pointers.
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
   b_mxhadgen = fChain->GetBranch("mxhadgen");
   b_pcmsgen = fChain->GetBranch("pcmsgen");
   b_tcmsgen = fChain->GetBranch("tcmsgen");
   b_fcmsgen = fChain->GetBranch("fcmsgen");
   b_ecmsgen = fChain->GetBranch("ecmsgen");
   b_pxhadgen = fChain->GetBranch("pxhadgen");
   b_txhadgen = fChain->GetBranch("txhadgen");
   b_fxhadgen = fChain->GetBranch("fxhadgen");
   b_exhadgen = fChain->GetBranch("exhadgen");
   b_GoodEvent = fChain->GetBranch("GoodEvent");
   b_vub = fChain->GetBranch("vub");
   b_vcb = fChain->GetBranch("vcb");
   b_vxbtyp = fChain->GetBranch("vxbtyp");
   b_other = fChain->GetBranch("other");
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
   b_deltam = fChain->GetBranch("deltam");
   b_mxwminpi = fChain->GetBranch("mxwminpi");
   b_wdeltam = fChain->GetBranch("wdeltam");
   b_mxphotmin = fChain->GetBranch("mxphotmin");
   b_photdeltam = fChain->GetBranch("photdeltam");
   b_minlat = fChain->GetBranch("minlat");
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
   b_enu = fChain->GetBranch("enu");
   b_pnu = fChain->GetBranch("pnu");
   b_tnu = fChain->GetBranch("tnu");
   b_fnu = fChain->GetBranch("fnu");
   b_mm2 = fChain->GetBranch("mm2");
   b_mm2nc = fChain->GetBranch("mm2nc");
   b_mm2fit = fChain->GetBranch("mm2fit");
   b_EPiz = fChain->GetBranch("EPiz");
   b_MinKMom = fChain->GetBranch("MinKMom");
   b_MaxKMom = fChain->GetBranch("MaxKMom");
   b_q2 = fChain->GetBranch("q2");
   b_q2nc = fChain->GetBranch("q2nc");
   b_q2fit = fChain->GetBranch("q2fit");
   return kTRUE;
}

void pstarfactor::Show(Int_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t pstarfactor::Cut(Int_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef pstarfactor_cxx

