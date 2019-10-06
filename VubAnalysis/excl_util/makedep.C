#define makedep_cxx
#include "makedep.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include <fstream.h>



void makedep::Loop(char *histfile, bool mctruth)
{
//   In a Root session, you can do:
//      Root > .L makedep.C
//      Root > makedep t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(i);  // read all branches
//by  b_branchname->GetEntry(i); //read only this branch

   char name[100];

   sprintf(name, "%s", histfile); 
   fHistFile = new TFile(name, "RECREATE");
   fHistFile->cd();
   cout << "Opened " << fHistFile->GetName() << endl;

   bookHist();
     
   if (fChain == 0) return;

   Int_t nentries = Int_t(fChain->GetEntries());
   cout << "Tot events :" << nentries << endl;

   Int_t nbytes = 0, nb = 0;
   for (Int_t jentry=0; jentry<nentries;jentry++) {
      Int_t ientry = LoadTree(jentry); //in case of a TChain, ientry is the entry number in the current file
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (jentry%10000 == 0) cout << "Event " << jentry << endl;
      // if (Cut(ientry) < 0) continue;
      if(pcms<1) continue;
//       if(vcb+vub == 0) continue;

      for(int i=1; i<101; i++){
	double cut = i/200.; 
	if(pur>cut)  {
	  if(brecocharge == 0){
 	    if(!mctruth || Gvxbtyp > 0){
	      sprintf(name, "h%d", 21000+i);  ((TH1D*)gDirectory->Get(name))->Fill(mes);
	      sprintf(name, "h%d", 11000+i);  ((TH1D*)gDirectory->Get(name))->Fill(mes);
	      if(mm2pi < 2.)  {
		sprintf(name, "h%d", 22000+i);  ((TH1D*)gDirectory->Get(name))->Fill(mes);
		sprintf(name, "h%d", 12000+i);  ((TH1D*)gDirectory->Get(name))->Fill(mes);
	      }
	      if(mm2pi < 2. && nchg <3)  {
		sprintf(name, "h%d", 23000+i);  ((TH1D*)gDirectory->Get(name))->Fill(mes);
		sprintf(name, "h%d", 13000+i);  ((TH1D*)gDirectory->Get(name))->Fill(mes);
	      }
 	    }
	  }else{
 	    if(!mctruth || Gvxbtyp < 0){
	      sprintf(name, "h%d", 31000+i);  ((TH1D*)gDirectory->Get(name))->Fill(mes);
	      sprintf(name, "h%d", 11000+i);  ((TH1D*)gDirectory->Get(name))->Fill(mes);
	      if(mm2pi < 2.)  {
		sprintf(name, "h%d", 32000+i);  ((TH1D*)gDirectory->Get(name))->Fill(mes);
		sprintf(name, "h%d", 12000+i);  ((TH1D*)gDirectory->Get(name))->Fill(mes);
	      }
	      if(mm2pi < 2. && nchg <3)  {
		sprintf(name, "h%d", 33000+i);  ((TH1D*)gDirectory->Get(name))->Fill(mes);
		sprintf(name, "h%d", 13000+i);  ((TH1D*)gDirectory->Get(name))->Fill(mes);
	      }
 	    }
	  }
	}	  			
      }
   }
   fHistFile->cd();
   fHistFile->Write();
   fHistFile->Close();
   delete fHistFile;
}

// ----------------------------------------------------------------------
TChain * getchain(char *thechain) {

  TChain *chain = new TChain("events");
  cout << "Chaining ... " << thechain << endl;
  char pName[2000]; 
  char buffer[200];
  sprintf(buffer, "%s", thechain);
  ifstream is(buffer);  
  cout << "files " << buffer <<  endl;
  while(is.getline(buffer, 200, '\n')){
    //    if (buffer[0] == '#') continue; 
    sscanf(buffer, "%s", pName); 
    cout << "   Add: " << buffer << endl; 
    chain->Add(pName); 
  }
  is.close();
  return chain;

}

void makedep::bookHist()
{
  char name[100],title[100];
  fHistFile->cd();  

  TH1 *h;

  for (int i = 1; i < 101; ++i) {
    sprintf(name, "h%d", 11000+i);  sprintf(title, "mes per pur lepton");  h = new TH1D(name, title, 40, 5.2, 5.3); 
    sprintf(name, "h%d", 12000+i);  sprintf(title, "mes per pur lepton and mm2");  h = new TH1D(name, title, 40, 5.2, 5.3); 
    sprintf(name, "h%d", 13000+i);  sprintf(title, "mes per pur allcuts");  h = new TH1D(name, title, 40, 5.2, 5.3);     

    sprintf(name, "h%d", 21000+i);  sprintf(title, "mes per pur lepton B0");  h = new TH1D(name, title, 40, 5.2, 5.3); 
    sprintf(name, "h%d", 22000+i);  sprintf(title, "mes per pur lepton and mm2 B0");  h = new TH1D(name, title, 40, 5.2, 5.3); 
    sprintf(name, "h%d", 23000+i);  sprintf(title, "mes per pur allcuts B0");  h = new TH1D(name, title, 40, 5.2, 5.3);     

    sprintf(name, "h%d", 31000+i);  sprintf(title, "mes per pur lepton Bch");  h = new TH1D(name, title, 40, 5.2, 5.3); 
    sprintf(name, "h%d", 32000+i);  sprintf(title, "mes per pur lepton and mm2 Bch");  h = new TH1D(name, title, 40, 5.2, 5.3); 
    sprintf(name, "h%d", 33000+i);  sprintf(title, "mes per pur allcuts Bch");  h = new TH1D(name, title, 40, 5.2, 5.3);     
  }

}
