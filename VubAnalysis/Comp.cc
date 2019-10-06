#include "CompTool.hh"
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

  TString dir("/afs/slac.stanford.edu/u/ec/asarti/public_html/recoil/mx_stu/splot0");
  TString fileData(""), fileMc(""), var(""), cutFile("compcut.dat"), flag(""); 
  double min = -9999;  double shift =0; double smear =0;  double max =  9999;
  int isbch = 0; int bins =25;  int cat =7; int nevents = 100000000;    
  int sys(-1);
  //Example
  //../bin/Linux2/Comp -d /afs/slac.stanford.edu/u/ec/asarti/public_html/recoil/mx_stu/splot0 -fd data/data_090902_all-3.root -fm data/gene_090902_all-3.root -v mxhadfit -m 0 -M 5 -C 0 -b 18
  // -- command line arguments
  for (int i = 0; i < argc; i++){
    if(strcmp(argv[i],"-d") == 0)   dir      = TString(argv[++i]);                   // Outputdir    
    if(strcmp(argv[i],"-F") == 0)   flag      = TString(argv[++i]);                   // Outputdir    
    if(strcmp(argv[i],"-fd") == 0)  fileData = TString(argv[++i]);                   // Data File
    if(strcmp(argv[i],"-fm") == 0)  fileMc   = TString(argv[++i]);                   // Mc file
    if(strcmp(argv[i],"-v") == 0)   var   = TString(argv[++i]);                     // Variable under study
    if(strcmp(argv[i],"-sh") == 0)  shift = atoi(argv[++i]);                         // Shifting
    if(strcmp(argv[i],"-sm") == 0)  smear = atoi(argv[++i]);                         // Smearing
    if(strcmp(argv[i],"-m") == 0)   min = atoi(argv[++i]);                           // Min (range of given var)
    if(strcmp(argv[i],"-M") == 0)   max = atoi(argv[++i]);                           // Max (range of given var)
    if(strcmp(argv[i],"-C") == 0)   isbch = atoi(argv[++i]);                         // B selection: 0=B0; 1=B+; 2=All
    if(strcmp(argv[i],"-b") == 0)   bins = atoi(argv[++i]);                          // Number of bins for the studied var
    if(strcmp(argv[i],"-n") == 0)   nevents = atoi(argv[++i]);                           // Min (range of given var)
    if(strcmp(argv[i],"-c")  == 0) {cutFile    = TString(argv[++i]); }               // file with cuts
    if(strcmp(argv[i],"-S") == 0)   sys = atoi(argv[++i]);                           // Sys
  }

  TFile k(fileData.Data());  TFile j(fileMc.Data());
  CompTool h(var.Data(), min, max, bins, sys);
  if(sys >= 0) {
    if (strcmp(cutFile.Data(), "")) h.readCuts(cutFile, 0);
    h.dumpCuts();
    k.cd(); 
    h.Init((TTree*)gDirectory->Get("events"));
    h.Bookhist();             
    h.Loop(nevents,4,0,0,isbch,cat);                
    h.Fitmes(4,0,flag);
    h.Fitmes(4,1,flag);
    j.cd();
    h.Init((TTree*)gDirectory->Get("events"));
    h.Loop(nevents,5,0,0,isbch,cat);                
    h.Fitmes(5,0,flag);
    h.Fitmes(5,1,flag);
    h.overlap(0,0,dir,flag);
    h.overlap(1,0,dir,flag);
    h.sameP(dir,flag);
  } else {
    if (strcmp(cutFile.Data(), "")) h.readCuts(cutFile, 0);
    h.dumpCuts();
    k.cd(); 
    h.Init((TTree*)gDirectory->Get("events"));
    h.Bookhist();             
    h.Loop(nevents,4,0,0,isbch,cat);                
    h.Fitmes(4,0,flag);
    h.Fitmes(4,1,flag);
    j.cd();
    h.Init((TTree*)gDirectory->Get("events"));
    h.Loop(nevents,5,0,0,isbch,cat);                
    h.Fitmes(5,0,flag);
    h.Fitmes(5,1,flag);
    h.overlap(0,0,dir,flag);
    h.overlap(1,0,dir,flag);
    //    h.overlap(0,1,dir,flag);
    //    h.effplots(dir,flag);
    h.sameP(dir,flag);
  }
  return 0;
}


