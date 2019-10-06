#include "BkgbrekTool.hh"
#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>
#include <math.h>
#include <stdlib.h>

#include "util.hh"
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
#include "TTree.h"

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int main(int argc, char *argv[]) {
  gROOT->SetBatch(kTRUE);
  gROOT->SetStyle("Plain"); 

  TString dir("/afs/slac.stanford.edu/u/ec/asarti/public_html/recoil/mx_stu/splot0");
  TString fileData(""), var(""), var2(""), cutFile("compcut.dat"), flag(""); 
  double min = -9999;  double max =  9999; int norm(0); int cut(0);
  int isbch = 0; int bins =25; int nevents = 100000000;    
  int semilep = 0; int all =0; int laststudy = 0; int mxmns(0);
  double min2 = -9999;  double max2 =  9999;int bins2 =25; 
  //Example
  //../bin/Linux2/Comp -d /afs/slac.stanford.edu/u/ec/asarti/public_html/recoil/mx_stu/splot0 -fd data/data_090902_all-3.root -fm data/gene_090902_all-3.root -v mxhadfit -m 0 -M 5 -C 0 -b 18
  // -- command line arguments
  for (int i = 0; i < argc; i++){
    if(strcmp(argv[i],"-d") == 0)   dir      = TString(argv[++i]);                   // Outputdir    
    if(strcmp(argv[i],"-F") == 0)   flag      = TString(argv[++i]);                   // Outputdir    
    if(strcmp(argv[i],"-fd") == 0)  fileData = TString(argv[++i]);                   // Data File
    if(strcmp(argv[i],"-v") == 0)   var   = TString(argv[++i]);                     // Variable under study
    if(strcmp(argv[i],"-vv") == 0)   var2   = TString(argv[++i]);                     // Variable under study
    if(strcmp(argv[i],"-NM") == 0)  norm  = atoi(argv[++i]);                         // D breakdown
    if(strcmp(argv[i],"-x") == 0)   mxmns  = atoi(argv[++i]);                        // Mx vs Mnu study
    if(strcmp(argv[i],"-a") == 0)   all  = 1;                                        // D breakdown
    if(strcmp(argv[i],"-s") == 0)   semilep = 1;                                     // Semilep
    if(strcmp(argv[i],"-ls") == 0)  laststudy = 1;                                   // Last study
    if(strcmp(argv[i],"-m") == 0)   min = atoi(argv[++i]);                           // Min (range of given var)
    if(strcmp(argv[i],"-M") == 0)   max = atoi(argv[++i]);                           // Max (range of given var)
    if(strcmp(argv[i],"-mm") == 0)   min2 = atoi(argv[++i]);                           // Min (range of given var)
    if(strcmp(argv[i],"-MM") == 0)   max2 = atoi(argv[++i]);                           // Max (range of given var)
    if(strcmp(argv[i],"-C") == 0)   isbch = atoi(argv[++i]);                         // B selection: 0=B0; 1=B+; 2=All
    if(strcmp(argv[i],"-b") == 0)   bins = atoi(argv[++i]);                          // Number of bins for the studied var
    if(strcmp(argv[i],"-bb") == 0)   bins2 = atoi(argv[++i]);                          // Number of bins for the studied var
    if(strcmp(argv[i],"-n") == 0)   nevents = atoi(argv[++i]);                       // Min (range of given var)
    if(strcmp(argv[i],"-c")  == 0) {cutFile    = TString(argv[++i]); }                           // file with cuts
    if(strcmp(argv[i],"-lc")  == 0) {cut = 1; }                                                 // lepton cut
  }

  TFile k(fileData.Data());  
  if(mxmns) {
    BkgbrekTool h(var, min, max, bins, var2, min2, max2, bins2);
    if (strcmp(cutFile.Data(), "")) h.readCuts(cutFile, 0);
    h.dumpCuts();
    k.cd(); 
    h.Init((TTree*)gDirectory->Get("events"));
    h.Bookhist(mxmns);             
    h.Loop(nevents,4,isbch,0);                
    h.overlap(2,dir,flag);

  } else {
    BkgbrekTool h(var, min, max, bins,"NoComp");
    if (strcmp(cutFile.Data(), "")) h.readCuts(cutFile, 0);
    h.dumpCuts();
    k.cd(); 
    h.Init((TTree*)gDirectory->Get("events"));
    h.Bookhist(mxmns);             
    //D decays
    if(all) {
      h.Loop(nevents,4,isbch,0);                
      h.Fitmes(4,0,flag,0);
      h.overlap(0,dir,flag);
    }
    //Semilep
    if(semilep) {
      h.Loop(nevents,4,isbch,1);                
      h.Fitmes(4,0,flag,1);
      h.overlapS(0,dir,flag);
    }
    
    //All plots
    if(laststudy) {
      h.LoopT(nevents,4,0);                
      h.LoopT(nevents,5,1);                
      h.Fitmes(4,cut,flag,2);
      h.Fitmes(5,cut,flag,2);
      //      h.sameP(0,dir,flag,norm);
      //      h.sameK(0,dir,flag,norm);
      h.Brek(cut,dir,flag,norm);
      h.BrekB(cut,dir,flag,norm,0);
      h.BrekB(cut,dir,flag,norm,1);
    }
  }
  return 0;
}


