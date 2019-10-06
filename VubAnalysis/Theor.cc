#include "TheorTool.hh"
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
  TString var(""), var2(""), cutFile("compcut.dat"), flag(""); 
  double min = 0;  int max =  1; int plot = 0;
  int isbch = 0; TString bins("refWeights"); int nevents = 100000000;    
  int fitBR = 0; int seed(0); int doMba = 0; int gen = 0;
  // -- command line arguments
  for (int i = 0; i < argc; i++){
    if(strcmp(argv[i],"-d") == 0)   dir      = TString(argv[++i]);                   // Outputdir    
    if(strcmp(argv[i],"-F") == 0)   flag      = TString(argv[++i]);                   // Outputdir    
    if(strcmp(argv[i],"-v") == 0)   var   = TString(argv[++i]);                     // Variable under study
    if(strcmp(argv[i],"-fbr") == 0)  fitBR =  atoi(argv[++i]);                      // Fitted BR
    if(strcmp(argv[i],"-Mba") == 0)  doMba =  atoi(argv[++i]);                      // Mb and a reweighting
    if(strcmp(argv[i],"-M") == 0)    max =  atoi(argv[++i]);                      // Mb and a reweighting
    if(strcmp(argv[i],"-m") == 0)   min = 1;                                        // Smear BRexcl within PDG val
    if(strcmp(argv[i],"-C") == 0)   isbch = atoi(argv[++i]);                         // B selection: 0=B0; 1=B+; 2=All
    if(strcmp(argv[i],"-n") == 0)   nevents = atoi(argv[++i]);                       // Min (range of given var)
    if(strcmp(argv[i],"-s") == 0)   seed = atoi(argv[++i]);                       // Seed for random generation
    if(strcmp(argv[i],"-b") == 0)   bins  = TString(argv[++i]);                       // Weightfile
    if(strcmp(argv[i],"-pl") == 0)  plot = atoi(argv[++i]);                       // Seed for random generation
    if(strcmp(argv[i],"-c")  == 0) {cutFile    = TString(argv[++i]); }                           // file with cuts
  }

  TFile kda("/u/ec/ursl/root/anaQA-n00/csx-data.root");
  TFile kgen("/u/ec/ursl/root/anaQA-n00/csx-genb-new.root");
  TFile kmix("/u/ec/ursl/root/anaQA-n00/csx-vubmix.root");
  TFile knre("/u/ec/ursl/root/anaQA-n00/csx-vubnre.root");

  //TFile kcmix("/u/ec/ursl/root/anaQA-n00/csx-vubmix.root");
  //TFile kcnre("/u/ec/ursl/root/anaQA-n00/csx-vubnre.root");

  TFile kcmix("datasec/csx-brevubmix-old2001.root");
  TFile kcnre("datasec/csx-brevubnre-old2001.root");

  TheorTool h(var.Data(), min, max, bins.Data(), fitBR, seed);
  if (strcmp(cutFile.Data(), "")) h.readCuts(cutFile, 0);
  h.dumpCuts();  h.Bookhist();             
  h.initMba(90/1000);

  kcnre.cd(); 
  h.Init((TTree*)gDirectory->Get("events"));
  h.Loop(nevents,4,isbch,doMba); //Model               
  kcmix.cd();
  h.Init((TTree*)gDirectory->Get("events"));
  h.Loop(nevents,5,isbch,doMba); //MC  
  if(gen) {
    knre.cd(); 
    h.Init((TTree*)gDirectory->Get("events"));
    h.Loop(nevents,6,isbch,doMba); //Model gen              
    kmix.cd();
    h.Init((TTree*)gDirectory->Get("events"));
    h.Loop(nevents,7,isbch,doMba); //MC gen 
  }

  //Bkg subtraction  
  h.Fitmes(4,0,flag); //Model               
  h.Fitmes(5,0,flag); //MC               
  if(gen) {
    h.Fitmes(6,0,flag); //Model gen              
    h.Fitmes(7,0,flag); //MC gen              
  }
  if(doMba) {
    h.Fitmes(4,1,flag); //Model reweighted              
    h.Fitmes(5,1,flag); //MC reweighted              
  }
  //Works out weights on fermi motion
  if(doMba && plot) h.computeWMba(90/1000);

  h.overlap(doMba,dir,flag,0); //B0  weights
  h.overlap(doMba,dir,flag,1); //Bch weights
  h.PintStudy();
  return 0;
}


