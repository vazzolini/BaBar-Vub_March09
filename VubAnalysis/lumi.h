//////////////////////////////////////////////////////////
//   This class has been automatically generated 
//     (Mon Aug 19 08:24:03 2002 by ROOT version3.01/06)
//   from TTree events/events
//   found on file: scratch/alldata.root
//////////////////////////////////////////////////////////


#ifndef lumi_h
#define lumi_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class lumi {
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
   Int_t           mode;
   Int_t           nnpi0;
   Int_t           nnks;
   Int_t           nnpar;
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
   UChar_t         isDupli;
   Int_t           ValMap;
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
   Double_t        gmax;
   Double_t        mxhadfit;
   Double_t        mxminpi;
   Double_t        deltam;
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
   Double_t        MinKMom;
   Double_t        MaxKMom;
   Double_t        dx;
   Double_t        dy;
   Double_t        dz;
   Double_t        s2dxx;
   Double_t        s2dyy;
   Double_t        s2dzz;
   Double_t        s2dxy;
   Double_t        s2dyz;
   Double_t        s2dxz;
   Double_t        q2;
   Double_t        q2nc;
   Double_t        q2fit;
   Double_t        allksm0[12];
   Double_t        allksp[12];
   UChar_t         allksmc[12];
   Double_t        allchkp[4];
   UChar_t         allchkmc[4];
   Double_t        m0ks;
   Double_t        pks;
   Double_t        pksmc;

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
   TBranch        *b_intpur;
   TBranch        *b_brecoflav;
   TBranch        *b_brecocharge;
   TBranch        *b_brecomc;
   TBranch        *b_mode;
   TBranch        *b_nnpi0;
   TBranch        *b_nnks;
   TBranch        *b_nnpar;
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
   TBranch        *b_GoodEvent;
   TBranch        *b_isDupli;
   TBranch        *b_ValMap;
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
   TBranch        *b_gmax;
   TBranch        *b_mxhadfit;
   TBranch        *b_mxminpi;
   TBranch        *b_deltam;
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
   TBranch        *b_mm2nc;
   TBranch        *b_mm2fit;
   TBranch        *b_ENeu;
   TBranch        *b_MinKMom;
   TBranch        *b_MaxKMom;
   TBranch        *b_dx;
   TBranch        *b_dy;
   TBranch        *b_dz;
   TBranch        *b_s2dxx;
   TBranch        *b_s2dyy;
   TBranch        *b_s2dzz;
   TBranch        *b_s2dxy;
   TBranch        *b_s2dyz;
   TBranch        *b_s2dxz;
   TBranch        *b_q2;
   TBranch        *b_q2nc;
   TBranch        *b_q2fit;
   TBranch        *b_allksm0;
   TBranch        *b_allksp;
   TBranch        *b_allksmc;
   TBranch        *b_allchkp;
   TBranch        *b_allchkmc;
   TBranch        *b_m0ks;
   TBranch        *b_pks;
   TBranch        *b_pksmc;

   lumi(TTree *tree=0);
   ~lumi();
   Int_t  Cut(Int_t entry);
   Int_t  GetEntry(Int_t entry);
   Int_t  LoadTree(Int_t entry);
   void   Init(TTree *tree);
   void   Loop();
   Bool_t Notify();
   void   Show(Int_t entry = -1);
};

#endif

#ifdef lumi_cxx
lumi::lumi(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("scratch/alldata.root");
      if (!f) {
         f = new TFile("scratch/alldata.root");
      }
      tree = (TTree*)gDirectory->Get("events");

   }
   Init(tree);
}

lumi::~lumi()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t lumi::GetEntry(Int_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Int_t lumi::LoadTree(Int_t entry)
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

void lumi::Init(TTree *tree)
{
//   Set branch addresses
   if (tree == 0) return;
   fChain    = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run",&run);
   fChain->SetBranchAddress("lower",&lower);
   fChain->SetBranchAddress("upper",&upper);
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
   fChain->SetBranchAddress("GfDlep",&GfDlep);
   fChain->SetBranchAddress("GfDgam",&GfDgam);
   fChain->SetBranchAddress("intpur",&intpur);
   fChain->SetBranchAddress("brecoflav",&brecoflav);
   fChain->SetBranchAddress("brecocharge",&brecocharge);
   fChain->SetBranchAddress("brecomc",&brecomc);
   fChain->SetBranchAddress("mode",&mode);
   fChain->SetBranchAddress("nnpi0",&nnpi0);
   fChain->SetBranchAddress("nnks",&nnks);
   fChain->SetBranchAddress("nnpar",&nnpar);
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
   fChain->SetBranchAddress("isDupli",&isDupli);
   fChain->SetBranchAddress("ValMap",&ValMap);
   fChain->SetBranchAddress("vub",&vub);
   fChain->SetBranchAddress("vcb",&vcb);
   fChain->SetBranchAddress("vxbtyp",&vxbtyp);
   fChain->SetBranchAddress("other",&other);
   fChain->SetBranchAddress("xcharge",&xcharge);
   fChain->SetBranchAddress("pxhad",&pxhad);
   fChain->SetBranchAddress("txhad",&txhad);
   fChain->SetBranchAddress("fxhad",&fxhad);
   fChain->SetBranchAddress("exhad",&exhad);
   fChain->SetBranchAddress("mxhad",&mxhad);
   fChain->SetBranchAddress("gmax",&gmax);
   fChain->SetBranchAddress("mxhadfit",&mxhadfit);
   fChain->SetBranchAddress("mxminpi",&mxminpi);
   fChain->SetBranchAddress("deltam",&deltam);
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
   fChain->SetBranchAddress("MinKMom",&MinKMom);
   fChain->SetBranchAddress("MaxKMom",&MaxKMom);
   fChain->SetBranchAddress("dx",&dx);
   fChain->SetBranchAddress("dy",&dy);
   fChain->SetBranchAddress("dz",&dz);
   fChain->SetBranchAddress("s2dxx",&s2dxx);
   fChain->SetBranchAddress("s2dyy",&s2dyy);
   fChain->SetBranchAddress("s2dzz",&s2dzz);
   fChain->SetBranchAddress("s2dxy",&s2dxy);
   fChain->SetBranchAddress("s2dyz",&s2dyz);
   fChain->SetBranchAddress("s2dxz",&s2dxz);
   fChain->SetBranchAddress("q2",&q2);
   fChain->SetBranchAddress("q2nc",&q2nc);
   fChain->SetBranchAddress("q2fit",&q2fit);
   fChain->SetBranchAddress("allksm0",allksm0);
   fChain->SetBranchAddress("allksp",allksp);
   fChain->SetBranchAddress("allksmc",allksmc);
   fChain->SetBranchAddress("allchkp",allchkp);
   fChain->SetBranchAddress("allchkmc",allchkmc);
   fChain->SetBranchAddress("m0ks",&m0ks);
   fChain->SetBranchAddress("pks",&pks);
   fChain->SetBranchAddress("pksmc",&pksmc);
   Notify();
}

Bool_t lumi::Notify()
{
//   called when loading a new file
//   get branch pointers
   b_run = fChain->GetBranch("run");
   b_lower = fChain->GetBranch("lower");
   b_upper = fChain->GetBranch("upper");
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
   b_GfDlep = fChain->GetBranch("GfDlep");
   b_GfDgam = fChain->GetBranch("GfDgam");
   b_intpur = fChain->GetBranch("intpur");
   b_brecoflav = fChain->GetBranch("brecoflav");
   b_brecocharge = fChain->GetBranch("brecocharge");
   b_brecomc = fChain->GetBranch("brecomc");
   b_mode = fChain->GetBranch("mode");
   b_nnpi0 = fChain->GetBranch("nnpi0");
   b_nnks = fChain->GetBranch("nnks");
   b_nnpar = fChain->GetBranch("nnpar");
   b_mxhadgen = fChain->GetBranch("mxhadgen");
   b_pcmsgen = fChain->GetBranch("pcmsgen");
   b_tcmsgen = fChain->GetBranch("tcmsgen");
   b_fcmsgen = fChain->GetBranch("fcmsgen");
   b_ecmsgen = fChain->GetBranch("ecmsgen");
   b_pxhadgen = fChain->GetBranch("pxhadgen");
   b_txhadgen = fChain->GetBranch("txhadgen");
   b_fxhadgen = fChain->GetBranch("fxhadgen");
   b_exhadgen = fChain->GetBranch("exhadgen");
   b_kplus = fChain->GetBranch("kplus");
   b_GoodEvent = fChain->GetBranch("GoodEvent");
   b_isDupli = fChain->GetBranch("isDupli");
   b_ValMap = fChain->GetBranch("ValMap");
   b_vub = fChain->GetBranch("vub");
   b_vcb = fChain->GetBranch("vcb");
   b_vxbtyp = fChain->GetBranch("vxbtyp");
   b_other = fChain->GetBranch("other");
   b_xcharge = fChain->GetBranch("xcharge");
   b_pxhad = fChain->GetBranch("pxhad");
   b_txhad = fChain->GetBranch("txhad");
   b_fxhad = fChain->GetBranch("fxhad");
   b_exhad = fChain->GetBranch("exhad");
   b_mxhad = fChain->GetBranch("mxhad");
   b_gmax = fChain->GetBranch("gmax");
   b_mxhadfit = fChain->GetBranch("mxhadfit");
   b_mxminpi = fChain->GetBranch("mxminpi");
   b_deltam = fChain->GetBranch("deltam");
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
   b_MinKMom = fChain->GetBranch("MinKMom");
   b_MaxKMom = fChain->GetBranch("MaxKMom");
   b_dx = fChain->GetBranch("dx");
   b_dy = fChain->GetBranch("dy");
   b_dz = fChain->GetBranch("dz");
   b_s2dxx = fChain->GetBranch("s2dxx");
   b_s2dyy = fChain->GetBranch("s2dyy");
   b_s2dzz = fChain->GetBranch("s2dzz");
   b_s2dxy = fChain->GetBranch("s2dxy");
   b_s2dyz = fChain->GetBranch("s2dyz");
   b_s2dxz = fChain->GetBranch("s2dxz");
   b_q2 = fChain->GetBranch("q2");
   b_q2nc = fChain->GetBranch("q2nc");
   b_q2fit = fChain->GetBranch("q2fit");
   b_allksm0 = fChain->GetBranch("allksm0");
   b_allksp = fChain->GetBranch("allksp");
   b_allksmc = fChain->GetBranch("allksmc");
   b_allchkp = fChain->GetBranch("allchkp");
   b_allchkmc = fChain->GetBranch("allchkmc");
   b_m0ks = fChain->GetBranch("m0ks");
   b_pks = fChain->GetBranch("pks");
   b_pksmc = fChain->GetBranch("pksmc");
   return kTRUE;
}

void lumi::Show(Int_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t lumi::Cut(Int_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef lumi_cxx

