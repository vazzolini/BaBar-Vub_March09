#define purstudy_cxx
#include "purstudy.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include <fstream.h>

void purstudy::Loop(char *histfile, char *mode)
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

  sprintf(name, "%s%s", mode, "_puropt.root"); 
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
    
    bool ALLCUTS(0);
    if(!strcmp(mode,"pi")) ALLCUTS = TMath::Abs(mm2pi)<2. && nchg ==2 && nkp ==0 && mm2rho>-2. && wdeltam<-1.5;
    if(!strcmp(mode,"pi0")) ALLCUTS = mpi0>0 && TMath::Abs(mm2pi0)<2.&& nkp == 0 && nchg == 1;
    if(!strcmp(mode,"rho0")) ALLCUTS = mrho0>0.5 && mrho0<1. && TMath::Abs(mm2rho0)<2. && nkp ==0 && nchg ==3;
    if(!strcmp(mode,"rho")) ALLCUTS = mrho>0.5 && mrho<1. && TMath::Abs(mm2rho)<2. && nkp ==0 && nchg ==2 && wdeltam<-1.5;
    if(!strcmp(mode,"omega")) ALLCUTS = momega>0 && TMath::Abs(mm2omega)<2. && nkp ==0 && nchg ==3;
    if(!strcmp(mode,"eta")) ALLCUTS = meta>0 && TMath::Abs(mm2eta)<2. && nkp ==0 && nchg <4;
    if(!strcmp(mode,"etap")) ALLCUTS = metap>0 && TMath::Abs(mm2etap)<2. && nkp ==0 && nchg <6;
    if(!strcmp(mode,"a0")) ALLCUTS = ma0>0.9 && ma0<1.1 && TMath::Abs(mm2a0)<2. && nkp ==0 && nchg <4;
    if(!strcmp(mode,"a0p")) ALLCUTS = ma0p>0.9 && ma0p<1.1 && TMath::Abs(mm2a0p)<2. && nkp ==0 && nchg <5 && wdeltam<-1.5;
    
    if(!(pcms>1&&nle==1&&brecocharge+xcharge==0&&brecoflav+lcharge==0)) continue;
    if(!ALLCUTS) continue;
    for(int i=1; i<101; i++){
      double cut = i/200.; 
      if(pur>cut)  {
	sprintf(name, "h%d", 11000+i);  ((TH1D*)gDirectory->Get(name))->Fill(mes);
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

void purstudy::bookHist()
{
  char name[100],title[100];
  fHistFile->cd();  

  TH1 *h;

  for (int i = 1; i < 101; ++i) {
    sprintf(name, "h%d", 11000+i);  sprintf(title, "mes per pur lepton");  h = new TH1D(name, title, 40, 5.2, 5.3); 
  }

}
