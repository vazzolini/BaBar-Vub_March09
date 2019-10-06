
//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jan 20 03:59:02 2009 by ROOT version 4.04/02b
// from TTree PdfTree/PdfTree
// found on file: PDFfile.root
//////////////////////////////////////////////////////////

#ifndef toyclass_h
#define toyclass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <map>
#include "RooFitCore/RooDataSet.hh"

class RooRealVar;


class toyclass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leave types
   Float_t         B;
   Float_t         B_err;
   Float_t         S;
   Float_t         S_err;
   Float_t         ThoSigAlpha;
   Float_t         ThoSigAlpha_err;
   Float_t         ThoSigN;
   Float_t         ThoSigN_err;
   Float_t         ThoSigR;
   Float_t         ThoSigR_err;
   Float_t         ThoSigXc;
   Float_t         ThoSigXc_err;
   Float_t         ar;
   Float_t         ar_err;
   Float_t         cutoff;
   Float_t         cutoff_err;
   Float_t         sigma_l;
   Float_t         sigma_l_err;
   Float_t         sigma_r1;
   Float_t         sigma_r1_err;
   Float_t         sigma_r2;
   Float_t         sigma_r2_err;
   Int_t           mesfitIndex;
   Double_t        lumi;
   Int_t           identifier;

   // List of branches
   TBranch        *b_B;   //!
   TBranch        *b_B_err;   //!
   TBranch        *b_S;   //!
   TBranch        *b_S_err;   //!
   TBranch        *b_ThoSigAlpha;   //!
   TBranch        *b_ThoSigAlpha_err;   //!
   TBranch        *b_ThoSigN;   //!
   TBranch        *b_ThoSigN_err;   //!
   TBranch        *b_ThoSigR;   //!
   TBranch        *b_ThoSigR_err;   //!
   TBranch        *b_ThoSigXc;   //!
   TBranch        *b_ThoSigXc_err;   //!
   TBranch        *b_ar;   //!
   TBranch        *b_ar_err;   //!
   TBranch        *b_cutoff;   //!
   TBranch        *b_cutoff_err;   //!
   TBranch        *b_sigma_l;   //!
   TBranch        *b_sigma_l_err;   //!
   TBranch        *b_sigma_r1;   //!
   TBranch        *b_sigma_r1_err;   //!
   TBranch        *b_sigma_r2;   //!
   TBranch        *b_sigma_r2_err;   //!
   TBranch        *b_mesfitIndex;   //!
   TBranch        *b_luminosity;   //!
   TBranch        *b_identifier;   //!

   //s   RooDataSet* rds;

   Double_t LUMI_DATA_OLD;
   Double_t LUMI_GENERIC_OLD;
   Double_t LUMI_SIGNAL_OLD;

   toyclass(const TString &);
   virtual ~toyclass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(Double_t,Double_t,Double_t,Int_t ,Int_t);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   void RetrieveInfo();
};

#endif

#ifdef toyclass_cxx
/* toyclass::toyclass(TTree *tree) */
/* { */
/* // if parameter tree is not specified (or zero), connect the file */
/* // used to generate this class and read the Tree. */
/*    if (tree == 0) { */
/*       TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("PDFfile.root"); */
/*       if (!f) { */
/*          f = new TFile("PDFfile.root"); */
/*       } */
/*       tree = (TTree*)gDirectory->Get("PdfTree"); */

/*    } */
/*    Init(tree); */
/* } */

toyclass::~toyclass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t toyclass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t toyclass::LoadTree(Long64_t entry)
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

void toyclass::Init(TTree *tree)
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

   fChain->SetBranchAddress("B",&B);
   fChain->SetBranchAddress("B_err",&B_err);
   fChain->SetBranchAddress("S",&S);
   fChain->SetBranchAddress("S_err",&S_err);
   fChain->SetBranchAddress("ThoSigAlpha",&ThoSigAlpha);
   fChain->SetBranchAddress("ThoSigAlpha_err",&ThoSigAlpha_err);
   fChain->SetBranchAddress("ThoSigN",&ThoSigN);
   fChain->SetBranchAddress("ThoSigN_err",&ThoSigN_err);
   fChain->SetBranchAddress("ThoSigR",&ThoSigR);
   fChain->SetBranchAddress("ThoSigR_err",&ThoSigR_err);
   fChain->SetBranchAddress("ThoSigXc",&ThoSigXc);
   fChain->SetBranchAddress("ThoSigXc_err",&ThoSigXc_err);
   fChain->SetBranchAddress("ar",&ar);
   fChain->SetBranchAddress("ar_err",&ar_err);
   fChain->SetBranchAddress("cutoff",&cutoff);
   fChain->SetBranchAddress("cutoff_err",&cutoff_err);
   fChain->SetBranchAddress("sigma_l",&sigma_l);
   fChain->SetBranchAddress("sigma_l_err",&sigma_l_err);
   fChain->SetBranchAddress("sigma_r1",&sigma_r1);
   fChain->SetBranchAddress("sigma_r1_err",&sigma_r1_err);
   fChain->SetBranchAddress("sigma_r2",&sigma_r2);
   fChain->SetBranchAddress("sigma_r2_err",&sigma_r2_err);
   fChain->SetBranchAddress("mesfitIndex",&mesfitIndex);
   fChain->SetBranchAddress("lumi",&lumi);
   fChain->SetBranchAddress("identifier",&identifier);
   Notify();
}

Bool_t toyclass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. Typically here the branch pointers
   // will be retrieved. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed.

   // Get branch pointers
   b_B = fChain->GetBranch("B");
   b_B_err = fChain->GetBranch("B_err");
   b_S = fChain->GetBranch("S");
   b_S_err = fChain->GetBranch("S_err");
   b_ThoSigAlpha = fChain->GetBranch("ThoSigAlpha");
   b_ThoSigAlpha_err = fChain->GetBranch("ThoSigAlpha_err");
   b_ThoSigN = fChain->GetBranch("ThoSigN");
   b_ThoSigN_err = fChain->GetBranch("ThoSigN_err");
   b_ThoSigR = fChain->GetBranch("ThoSigR");
   b_ThoSigR_err = fChain->GetBranch("ThoSigR_err");
   b_ThoSigXc = fChain->GetBranch("ThoSigXc");
   b_ThoSigXc_err = fChain->GetBranch("ThoSigXc_err");
   b_ar = fChain->GetBranch("ar");
   b_ar_err = fChain->GetBranch("ar_err");
   b_cutoff = fChain->GetBranch("cutoff");
   b_cutoff_err = fChain->GetBranch("cutoff_err");
   b_sigma_l = fChain->GetBranch("sigma_l");
   b_sigma_l_err = fChain->GetBranch("sigma_l_err");
   b_sigma_r1 = fChain->GetBranch("sigma_r1");
   b_sigma_r1_err = fChain->GetBranch("sigma_r1_err");
   b_sigma_r2 = fChain->GetBranch("sigma_r2");
   b_sigma_r2_err = fChain->GetBranch("sigma_r2_err");
   b_mesfitIndex = fChain->GetBranch("mesfitIndex");
   b_luminosity = fChain->GetBranch("lumi");
   b_identifier = fChain->GetBranch("identifier");

   return kTRUE;
}

void toyclass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}


Int_t toyclass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef toyclass_cxx
