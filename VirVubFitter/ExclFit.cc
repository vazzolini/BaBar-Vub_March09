#include <iostream>
#include <fstream>
#include <iomanip>
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

#include "RecoilAnalysis/recoilDSys.hh"
#include "VirVubFitter/exclfitNtp.hh"

class RooDataSet;
class RooRealVar;
class RooAbsPdf;

using namespace std;
using namespace RooFit;
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int main(int argc, char *argv[]) {
  
  gROOT->SetBatch(kTRUE);
  
  TString texPrefix("all"); 
  TString progName  = argv[0]; progName.ReplaceAll("bin/", "");
  TString fileName, optionFile(""),wFile(""),cutFile(""),mesFile(""),prefixOut,dirName;
  TString wFermiFile("");
  ////DEFAULT tree name
  TString treeName("events");
  int sys = 0;
  int hyb = 0;
  int ck=0, wsys=0, fermi=0, bsys=-99;
  // -- command line arguments
  for (int i = 0; i < argc; i++){
    if(strcmp(argv[i],"-F")  == 0) {optionFile    = TString(argv[++i]); } //filename for option file   
    if(strcmp(argv[i],"-C")  == 0) {cutFile    = TString(argv[++i]); }  // filename for cuts file  
    if(strcmp(argv[i],"-ntp")  ==0) {treeName = TString(argv[++i]); }  //set Tree name : events in CM1, ntp1 in CM2
    if(strcmp(argv[i],"-Mes")  == 0) {mesFile     = TString(argv[++i]);  } // filename for mes fit parameter file
    if(strcmp(argv[i],"-P")  == 0) {prefixOut    = TString(argv[++i]); }  // prefix for all output files  
    if(strcmp(argv[i],"-D")  == 0) {dirName     = TString(argv[++i]);  } // directory for output files
    if(strcmp(argv[i],"-s")  == 0) {texPrefix = TString(argv[++i]);  }// prefix for output tex file
    if(strcmp(argv[i],"-Sys")  == 0) {sys = atoi(argv[++i]);  }
    if(strcmp(argv[i],"-Hyb")  == 0) {hyb = 1;  }


  }
  
  // -- This is the real analysis part
  
  //TFile g("/u/ec/ursl/root/anaQA-j00/csx-vubmix.root");
  
  TTree* ntuple = dynamic_cast<TTree*>( gDirectory->Get( treeName.Data() ) ) ;

  exclfitNtp a( ntuple ,sys,wFile);
  
  if(hyb) a.doTheo();

  if (strcmp(optionFile.Data(), "")) a.readOptions(optionFile, 0);
  a.dumpOptions();
  
  if (strcmp(cutFile.Data(), "")) a.readCuts(cutFile, 0);
  a.dumpCuts();
  
  if (strcmp(mesFile.Data(), "")) a.readmesParam(mesFile, 0);
  a.dumpmesParam();
  
  if (strcmp(prefixOut.Data(), "")) a.setPrefix(prefixOut);
  a.setTexPrefix(texPrefix); 
  
  if (strcmp(dirName.Data(), "")) a.setDirectory(dirName);
  
  a.openHistFile("fitresult.root");
  
  bool *isc = a.getfilechain();

  // b->clnu and "other" component model
  
  cout << isc[0] << " " << isc[1] << " " << isc[2] << " " << isc[3] << endl;

  TFile ga;
  if(!isc[1]){
    ga.Open((a.getfileVcb()).Data());
    a.Init((TTree*)gDirectory->Get(treeName.Data()));
  }else{
    a.Init(a.getchain((char*)((a.getfileVcb()).Data()), &treeName ));
  }
  a.Bookhist();
  a.Loop(0,1,111111111,1,0);
  
  // b->ulnu model
  
  if(!hyb){
    TFile h;
    if(!isc[0]){
      h.Open((a.getfileVubTotal()).Data());
      a.Init((TTree*)gDirectory->Get(treeName.Data()));
    }else{
      a.Init(a.getchain((char*)((a.getfileVubTotal()).Data()),  &treeName ) );
    }
    a.Loop(0,0,111111111,1,0);
  }else{
    TFile h1;
    if(!isc[0]){
      h1.Open((a.getfileVubTotal()).Data());
      a.Init((TTree*)gDirectory->Get(treeName.Data()));
    }else{
      a.Init(a.getchain((char*)((a.getfileVubTotal()).Data()), &treeName) );
    }
    a.Loop(0,0,111111111,1,0);
    
    TFile h2;
    if(!isc[3]){
      h2.Open((a.getfileVubTotalnres()).Data());
      a.Init((TTree*)gDirectory->Get(treeName.Data()));
    }else{
      a.Init(a.getchain((char*)((a.getfileVubTotalnres()).Data()), &treeName));
    }
    a.Loop(0,0,111111111,1,1);
       
  }
  // fit the data file...
  
  TFile o;
  if(!isc[2]){
    o.Open((a.getfileData()).Data());
    a.Init((TTree*)gDirectory->Get(treeName.Data()));
  }else{
    a.Init(a.getchain((char*)((a.getfileData()).Data()), &treeName));
  }
  a.Loop(1,0,1111111111,0,0);
    
  // fit Mx distribution and Vub extraction
   
  a.theFit();
  
  a.closeHistFile();
  
  return 0;
}

