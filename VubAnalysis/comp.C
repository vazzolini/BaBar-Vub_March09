#include "TLatex.h"
#include "thecomparison.C";
// bSel= 0 B0, 1 B-, 2 all
//       -1 B0, -2 B-, -3 all  separately vub and vcb (file1 contains vcb, file2 vub)
void comp(TString file1, TString file2, TString dir, char *var, double min, double max, int bins, double shift, double smear, int isbch, int cat,int sys){
  TFile k(file1);
  TFile j(file2);
  thecomparison h(var, min, max, bins);              
  k.cd(); 
  //  h.Init(h3);
  h.Init(events);
  gROOT.cd();
  h.Bookhist()    ;             
  gROOT.cd();

  h.Loop(100000000,4,0,0, isbch,cat,1,0);                
  h.Fitmes(4,0);
  h.Fitmes(4,1);
  
  j.cd();
  h.Init(events);
  gROOT.cd();  
  h.Loop(100000000,5,shift,smear,isbch,cat,1,sys);                
  gROOT.cd();

  h.Fitmes(5,0);
  h.Fitmes(5,1);
  h.overlap(0,0,dir);
  h.overlap(1,0,dir);
  h.overlap(0,1,dir);
  h.effplots(dir);
}
