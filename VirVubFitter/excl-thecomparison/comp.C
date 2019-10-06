//#include "TLatex.h"
#include <iostream>
//#include "thecomparison.h"
//#include "thecomparison.C"

void comp(TString file1, char* tree1, TString file2, char* tree2, TString dir, TString cutFile, char *var, double min, double max, int bins,double shift, double smear){
  
    cout<<"Start ROOT Settings\n";
    gSystem->Load("libPhysics.so");
    gSystem->Load("libthecomparison.so");
    
    gROOT->SetStyle("Plain");
    gStyle->SetPalette(1);
    gStyle->SetOptStat(1111111);  // Show overflow, underflow + SumOfWeights 
    //  gStyle->SetStatStyle(0); // for a completely transparent stat box
    gStyle->SetOptFit(111110); 
    gStyle->SetOptFile(1); 
    gStyle->SetMarkerStyle(20);
    gStyle->SetMarkerSize(0.8);
    gStyle->SetMarkerColor(1);
    gStyle->SetTitleBorderSize(0);  // no border around histogram title (font size can't be calculated anyways ...)
    gROOT->ForceStyle();
    
    cout<<"End Root Settings\n";
    
    cout<<" var = "<<var<<" min "<<min<<" max "<<max<<" bins "<<bins<<endl;
    
    thecomparison h(var, min, max, bins);              
    //  k.cd(); 
    //  h.Init(h3);
    h.Init(h.getchain(file1.Data(),tree1));
    gROOT->cd();
    h.readCuts(cutFile);
    h.Bookhist();             
    gROOT->cd();
  
    h.Loop(100000000,4,0,0);                
    h.Fitmes(4,0);
    h.Fitmes(4,1);

    //  j.cd();
    h.Init(h.getchain(file2.Data(),tree2));
    gROOT->cd();
    h.Loop(100000000,5,0,0);                
    gROOT->cd();
    
    h.Fitmes(5,0);
    //  std::cout<<"FitMes(5,0) done"<<endl;
    h.Fitmes(5,1);
    // std::cout<<"FitMes(5,1) done"<<endl;

    h.overlap(0,0,dir); //cut=0(leptoncut) norm=0
    h.overlap(1,0,dir); //cut=1 norm=0
    h.overlap(0,1,dir); //cut=0 norm=1
    //    h.effplots(dir);
}

