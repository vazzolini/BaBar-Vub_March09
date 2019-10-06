//#include "TLatex.h"
#include <iostream>
//#include "thecomparison.h"
//#include "thecomparison.C"

void comp(TString file1, char* tree1, TString file2, char* tree2, TString dir, char *var, double min, double max, int bins, double shift, double smear, int isbch, int cat,int sys){
  
    cout<<"Start ROOT Settings\n";
    gSystem->Load("libPhysics.so");
    gSystem->Load("libthecomparison.so");
    
    gROOT->SetStyle("Plain");
    gStyle->SetStatX(0.91);  //to move the stat Window a little bit to the left
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
    //TCanvas *c = 0;
    //c = (TCanvas*)gROOT->FindObject("c0"); if (c) c->Delete(); c = 0;
    //p = (TPad*)gROOT->FindObject("p0"); if (p) p->Delete(); p = 0;
    // --- Create a new canvas.
    
    //TCanvas c0("c0","--c0--",472,0,800,900);
    //&c0->ToggleEventStatus();
    
    //TLatex tl;
    //tl.SetNDC(kTRUE);
    
    cout<<"End Root Settings\n";
    
    // TFile k(file1);
    // TFile j(file2);

    cout<<" var = "<<var<<" min "<<min<<" max "<<max<<" bins "<<bins<<endl;
    
    thecomparison h(var, min, max, bins);              
    //  k.cd(); 
    //  h.Init(h3);
    h.Init(h.getchain(file1.Data(),tree1));
    gROOT->cd();
    h.Bookhist();             
    gROOT->cd();
  
    h.Loop(100000000,4,0,0, isbch,cat,1,0);                
    h.Fitmes(4,0,1);
    h.Fitmes(4,1,1);

    //  j.cd();
    h.Init(h.getchain(file2.Data(),tree2));
    gROOT->cd();
    h.Loop(100000000,5,shift,smear,isbch,cat,1,sys);                
    gROOT->cd();
    
    h.Fitmes(5,0,1);
    //  std::cout<<"FitMes(5,0) done"<<endl;
    h.Fitmes(5,1,1);
    // std::cout<<"FitMes(5,1) done"<<endl;

    h.overlap(0,0,dir); //cut=0(leptoncut) norm=0
    h.overlap(1,0,dir); //cut=1 norm=0
    h.overlap(0,1,dir); //cut=0 norm=1
    h.effplots(dir);
}

