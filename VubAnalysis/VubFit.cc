#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>
#include <math.h>
#include <stdlib.h>

#include "TROOT.h"
#include "TRint.h"
#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TString.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TMinuit.h"
#include "TString.h"

#include "VubAnalysis/recoilDSys.hh"
#include "VubAnalysis/fitNtp.hh"


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int main(int argc, char *argv[]) {

  gROOT->SetBatch(kTRUE);

  TString texPrefix("all"); int prlRew = 1;
  TString progName  = argv[0]; progName.ReplaceAll("bin/", "");
  TString fileName, optionFile(""),wFile(""),cutFile(""),mesFile(""),prefixOut,dirName;
  int sys = 0, fermi=0, theosys =0;
  // -- command line arguments
  for (int i = 0; i < argc; i++){
    if(strcmp(argv[i],"-F")  == 0) {optionFile    = TString(argv[++i]); }    
    if(strcmp(argv[i],"-W")  == 0) {prlRew = 0; wFile = TString(argv[++i]); }
    if(strcmp(argv[i],"-C")  == 0) {cutFile    = TString(argv[++i]); }    
    if(strcmp(argv[i],"-Mes")  == 0) {mesFile     = TString(argv[++i]);  }
    if(strcmp(argv[i],"-P")  == 0) {prefixOut    = TString(argv[++i]); }    
    if(strcmp(argv[i],"-D")  == 0) {dirName     = TString(argv[++i]);  }
    if(strcmp(argv[i],"-s")  == 0) {texPrefix = TString(argv[++i]);  }
    if(strcmp(argv[i],"-Sys")  == 0) {sys = atoi(argv[++i]);  }
    if(strcmp(argv[i],"-Fermi")  == 0) {fermi=1;}
    if(strcmp(argv[i],"-Theosys")  == 0) {theosys=1;}
  }

  // -- This is the real analysis part
  
  TFile g("/u/ec/ursl/root/anaQA-j00/csx-vubmix.root");

  fitNtp a((TTree*)gDirectory->Get("events"),sys,wFile,prlRew);

  if (strcmp(optionFile.Data(), "")) a.readOptions(optionFile, 0);
  a.dumpOptions();
  
  if (theosys) a.doTheo();
  if (strcmp(cutFile.Data(), "")) a.readCuts(cutFile, 0);
  a.dumpCuts();

  if (strcmp(mesFile.Data(), "")) a.readmesParam(mesFile, 0);
  a.dumpmesParam();
 
  if (strcmp(prefixOut.Data(), "")) a.setPrefix(prefixOut);
  a.setTexPrefix(texPrefix); 

  if (strcmp(dirName.Data(), "")) a.setDirectory(dirName);

  if (fermi) a.doFermi();

  a.openHistFile("fitresult.root");
  
  // b->clnu and "other" component model
  TFile ga((a.getfileVcb()).Data());
  a.Init((TTree*)gDirectory->Get("events"));                             
  a.Bookhist();
  a.Loop(0,1,1111111111,1,0);

  // b->ulnu model
  
  if(!theosys){
    TFile h((a.getfileVubTotal()).Data());
    a.Init((TTree*)gDirectory->Get("events"));                             
    a.Loop(0,0,111111111,1,0);
  }else{
    TFile h1((a.getfileVubTotalres()).Data());
    a.Init((TTree*)gDirectory->Get("events"));                             
    a.Loop(0,0,111111111,1,0);
    TFile h2((a.getfileVubTotalnres()).Data());
    a.Init((TTree*)gDirectory->Get("events"));                             
    a.Loop(0,0,111111111,1,1);
  }
  
  
  // fit the data file...
  
  if(a.isfitMC()){
     // data file from MC...

     TFile i((a.getfileVcb()).Data());
     a.Init((TTree*)gDirectory->Get("events"));                             
     int vcbstat= (427000/270)*a.totalStat();
     a.Loop(1,1,vcbstat,1,0);
  
     TFile l((a.getfileVubTotal()).Data());                                 
     a.Init((TTree*)gDirectory->Get("events"));
     a.Loop(1,0,int((6400/(270*0.0017))*a.totalStat()*a.BRRatioGenValue()),1,0);
  
   }else{
     // data file from real data...

     TFile o((a.getfileData()).Data());                                 
     a.Init((TTree*)gDirectory->Get("events"));                             
     a.Loop(1,0,1111111111,0,0);
   }
  
  // fit mes distributions


  a.Fitmessub("data");                                    

  a.Fitmessub("cvcb");
  
  a.Fitmessub("uvub");                                    
  
  a.Fitmessub("other");                                    
  
  

  // fit Mx distribution and Vub extraction
   
   a.theFit();
   
   a.closeHistFile();
   
   return 0;
}

