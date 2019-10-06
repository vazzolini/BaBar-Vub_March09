
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>

#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TString.h"

#include "RecoilAnalysis/recoilDSys.hh"
#include "VirVubFitter/b2uClass.hh"


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int main(int argc, char *argv[]) {

  gROOT->SetBatch(kTRUE);
  int doth = 0;
  TString texPrefix("all"); 
  TString progName  = argv[0]; progName.ReplaceAll("bin/", "");
  TString fileName, optionFile(""),cutFile(""),wFile(""),b2uwFile("");
  TString q2File(""), mxFile(""), mxFile1d("");
  TString mesFile(""),prefixOut(""),dirName(""), UFile("");
  int sys = 0, fermi=0, theosys =0, q2 =0, unf=0, munf=0, comb=0, Sun=0, ck =0, me=0, mu=0, usemxfit=0, newbin=0;
  int NORE=0, RE=0;
  bool iscm2=false, varfit=false;
  vector<float> wFermivec;
  double hiunf=0., wsys=0., width=0.3;
  // -- command line arguments
  for (int i = 0; i < argc; i++){
    if(strcmp(argv[i],"-F")  == 0)  {optionFile    = TString(argv[++i]); }  
    if(strcmp(argv[i],"-W")  == 0)  {wFile = TString(argv[++i]);}
    if(strcmp(argv[i],"-b2uw")  == 0)  {b2uwFile = TString(argv[++i]);}
    if(strcmp(argv[i],"-Sys")  == 0){sys = atoi(argv[++i]);  }
    if(strcmp(argv[i],"-q2")  == 0) {q2 = 1;  }
    if(strcmp(argv[i],"-comb")  == 0) {comb = 1;  }
    if(strcmp(argv[i],"-mxfit")  == 0) {usemxfit = 1;  }
    if(strcmp(argv[i],"-Sun")  == 0) {Sun = 1;  }
    if(strcmp(argv[i],"-ckBauer")  == 0) {ck = 1;  }
    // if(strcmp(argv[i],"-nore")  == 0) {NORE = 1;  }
    //if(strcmp(argv[i],"-re")  == 0) {RE = 1;  }
    if(strcmp(argv[i],"-unf")  == 0) {unf = atoi(argv[++i]);  } // Number of bins for unfolding, KT
    if(strcmp(argv[i],"-hiunf")  == 0) {hiunf = atof(argv[++i]);} // Endpoint for unfolding, KT
    if(strcmp(argv[i],"-mult")  == 0) {munf = 1;  } // Flag unfolding for multiplicity categories, KT
    if(strcmp(argv[i],"-Fermi")  == 0) {fermi=1;}
    if(strcmp(argv[i],"-Mes")  == 0) {mesFile     = TString(argv[++i]);  }
    if(strcmp(argv[i],"-s")  == 0) {texPrefix = TString(argv[++i]);  }
    if(strcmp(argv[i],"-P")  == 0) {prefixOut    = TString(argv[++i]); } 
    if(strcmp(argv[i],"-me")  == 0) {me = atoi(argv[++i]);  }
    if(strcmp(argv[i],"-mu")  == 0) {mu = atoi(argv[++i]);  }
    if(strcmp(argv[i],"-C")  == 0) {cutFile    = TString(argv[++i]); } 
    if(strcmp(argv[i],"-D")  == 0) {dirName     = TString(argv[++i]);  }
    if(strcmp(argv[i],"-U")  == 0) {UFile       = TString(argv[++i]); }
    if(strcmp(argv[i],"-me")  == 0) {me = atoi(argv[++i]);  }
    if(strcmp(argv[i],"-mu")  == 0) {mu = atoi(argv[++i]);  }
    if(strcmp(argv[i],"-q2V")  == 0)  {q2File = TString(argv[++i]);}
    if(strcmp(argv[i],"-mxV")  == 0)  {mxFile = TString(argv[++i]);}
    if(strcmp(argv[i],"-mx1d")  == 0)  {mxFile1d = TString(argv[++i]);}
    if(strcmp(argv[i],"-wisys")  == 0) {wsys = atof(argv[++i]);  }
    if(strcmp(argv[i],"-newbin")  == 0)  {newbin = 1;}
    if(strcmp(argv[i],"-width")  == 0) {width = atof(argv[++i]);  }
  }

  b2uClass a((TTree*)gDirectory->Get("events"),wFile,sys,q2, comb,Sun,unf,hiunf,me,mu,usemxfit,iscm2,varfit,wFermivec,newbin);

  a.b2uGetWeights(b2uwFile);

  //Read files from thefiles/ dir
  if (strcmp(optionFile.Data(), "")) a.readOptions(optionFile, 0);
  a.dumpOptions();
  //Read cuts for ascii file in theset/ dir

  if (strcmp(cutFile.Data(), "")) a.readCuts(cutFile, 0, wsys); 
  cout << endl;
  cout << "now dumping the cuts" << endl;
  a.dumpCuts();

  //Set the prefix for output storage
  if (strcmp(prefixOut.Data(), "")) a.setPrefix(prefixOut);
  a.setTexPrefix(texPrefix); 

  //Read the mES parameters obtained from global fits
  //This helps reducing instability of mES fits.
  if (strcmp(mesFile.Data(), "")) a.readmesParam(mesFile, 0);
  cout << "now readmesParam" << endl;
  a.dumpmesParam();

  if (strcmp(q2File.Data(), "")) a.q2Binning(q2File); /////////////////////////////
  if (strcmp(mxFile.Data(), "")) a.mxBinning(mxFile); /////////////////////////////
  if (newbin==0) {
    if (strcmp(mxFile1d.Data(), "")) a.mxBinning1d(mxFile1d); ///////////////////////
  } else {
    if (strcmp(mxFile1d.Data(), "")) a.b2umxBinning1d(width); ///////////////////////
  }

  //Read the parameters for the unfolding histos
  if (unf){
    if (strcmp(UFile.Data(), "")) a.readUnfParam(UFile, 0);
    a.dumpUnfParam();
  }

  if (strcmp(dirName.Data(), "")) a.setDirectory(dirName);

  a.openHistFile("Fitres.root");
  a.Bookhist();

  cout<<"isfitMC: "<<a.isfitMC()<<endl;
  bool *isc = a.getfilechain();

  if(a.isfitMC()){
    // data file from MC...
    cout<<"Entering in VcB"<<endl;
    TFile f5;
    if(!isc[4]){
      f5.Open((a.getfileVcb1()).Data());
      a.Init((TTree*)gDirectory->Get("events"));
    }else{
      a.Init(a.getchain((char*)((a.getfileVcb1()).Data())));
    }
    a.Loop(1,1,1,0,comb);   //generic  as  MC Bch

    TFile f7;
    if(!isc[5]){
      f7.Open((a.getfileVcb2()).Data());
      a.Init((TTree*)gDirectory->Get("events"));
    }else{
      a.Init(a.getchain((char*)((a.getfileVcb2()).Data())));
    }
    a.Loop(1,1,1,0,comb);   //generic  as  MC BO

    TFile f6;
    if(!isc[2]){
      f6.Open((a.getfileVubTotalnres()).Data());
      a.Init((TTree*)gDirectory->Get("events"));
    }else{
      a.Init(a.getchain((char*)((a.getfileVubTotalnres()).Data())));
    }
    a.Loop(1,0,1,1,comb);   //signal non resonant  as  MC
   
  }else{
    // data file from real data...

    TFile f4;
    if(!isc[6]){
      f4.Open((a.getfileData()).Data());
      a.Init((TTree*)gDirectory->Get("events"));
    }else{
      a.Init(a.getchain((char*)((a.getfileData()).Data())));
    }
    a.Loop(1,1,0,0,comb);  //data
  }
  
  TFile f3;
  if(!isc[3]){
    f3.Open((a.getfileVcb()).Data());
    a.Init((TTree*)gDirectory->Get("events"));
  }else{
    a.Init(a.getchain((char*)((a.getfileVcb()).Data())));
  }
  a.Loop(0,1,1,0,comb);     //generic 
  
  //  if(NORE){
  //   TFile f1((a.getfileVubTotalnres()).Data()); 
  //   a.Init((TTree*)gDirectory->Get("events"));                             
  //   a.Loop(0,0,1,1,comb);     //signal non resonant
  
  //     } else if(RE){
  //       TFile f2((a.getfileVubTotalres()).Data()); 
  //       a.Init((TTree*)gDirectory->Get("events"));                             
  //       a.Loop(0,0,1,0,comb);     //signal resonant
  
  //     } else{
  
  //    }
  

  TFile f1;
  if(!isc[2]){
    f1.Open((a.getfileVubTotalnres()).Data());
    a.Init((TTree*)gDirectory->Get("events"));
  }else{
    a.Init(a.getchain((char*)((a.getfileVubTotalnres()).Data())));
  }
  a.Loop(0,0,1,1,comb);     //signal non resonant
  TFile f2;
  if(!isc[1]){
    f2.Open((a.getfileVubTotalres()).Data());
    a.Init((TTree*)gDirectory->Get("events"));
  }else{
    a.Init(a.getchain((char*)((a.getfileVubTotalres()).Data())));
  }
  a.Loop(0,0,1,0,comb);     //signal resonant
  //    }
  
  //if(!(unf)){
  //  if(Sun){
  //    a.calcPstarFact(comb,1);
  //  }else{
  //    a.calcPstarFact(comb,0);
  //  }
  //}
  //a.Debug();


  if(!comb) {
    //Fit Mes distributions for cont subtraction
    
    if(munf){ //Fit Mes for multiplicity categories
      a.FitMes("data",1,0);
    
      a.FitMes("vcb",1,0);
    
      a.FitMes("vub",1,0);        
    
      a.FitMes("other",1,0);  
    }

    if(Sun){ //Fit Mes with signal unfolding

      a.FitMes("data",0,1);

      a.FitMes("vubin",0,1);
      
      a.FitMes("vubout",0,1);
      
      a.FitMes("vcboth",0,1);
      
    } else {
      
      a.FitMes("data",0,0);
      
      a.FitMes("vcb",0,0);
      
      a.FitMes("vub",0,0);
      
      a.FitMes("other",0,0);

    }
    
    
  } else {
    //Fit Mes distributions for cont subtraction
    
    if(Sun){ //Fit Mes with signal unfolding
      a.FitMes2D("data",0);
      
      a.FitMes2D("vubin",1);
      
      a.FitMes2D("vubout",1);        
      
      a.FitMes2D("vcboth",1);
    }else{
      
      a.FitMes2D("data",0);
      
      a.FitMes2D("vcb",0);
      
      a.FitMes2D("vub",0);
      
      a.FitMes2D("other",0);
    }
  }

  if(!(unf)){
    if(Sun){
      a.calcPstarFact(comb,1);
    }else{
      a.calcPstarFact(comb,0);
    }
  }

  if(Sun){
    a.theFit(comb,munf,1,ck);// a.theFit(comb);
  } else{
    a.theFit(comb,munf,0,ck);
  }

  if(Sun){
    a.b2upstar(1);
  } else {
    a.b2upstar(0);
  }

  if(Sun){
    a.b2ueff(1);
  } else {
    a.b2ueff(0);
  }

  cout << "Now the unfolding part" << endl;
  if(munf)
    a.MultCorr();
  if(unf)
    a.UnfHistos(munf);

  a.closeHistFile();

  return 0;
}
