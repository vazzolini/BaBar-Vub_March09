#include "TFile.h"
#include "TKey.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TH1D.h"
#include <iostream>
#include <vector>
using namespace std;

void shrinkPad(double b=0.1, double l=0.1, double r=0.1, double t=0.1);
void setFilledHist(TH1 *h, int lcol = kBlack, int fcol = kYellow, int fstyle = 1000, int width = 1);
void MergeRootFiles(const TString&, const TString&);
void MakeOverlap(Int_t, const TString& , const TString&, const TString&, const TString&);
Double_t chisquared(const TH1* h1, const TH1* h2);

void Plot(TString dir) {

  gROOT->SetStyle("Plain");
  gROOT->ProcessLine("gSystem->Load\(\"../../RooFitCore/tmp/libRooFitCore.so\"\)");

  vector<TString> varname;
  vector<TString> vartitle;
  vector<TString> vardescr;
  
   varname.push_back("pcms"); vartitle.push_back("P_{lep}"); vardescr.push_back(" GeV/c");
   varname.push_back("mm2"); vartitle.push_back("m_{#nu}^{2}"); vardescr.push_back(" GeV^{2}/c^{4}");
   varname.push_back("mx"); vartitle.push_back("m_{X}"); vardescr.push_back(" GeV/c^{2}");
   varname.push_back("pplus"); vartitle.push_back("P_{+}"); vardescr.push_back(" GeV/c");
   varname.push_back("q2"); vartitle.push_back("q^{2}"); vardescr.push_back(" GeV^{2}/c^{4}");
   varname.push_back("nchg"); vartitle.push_back("N. charged"); vardescr.push_back(" tracks");  
   varname.push_back("nneu"); vartitle.push_back("N. neutral"); vardescr.push_back(" clusters");  
   varname.push_back("emiss"); vartitle.push_back("Missing Energy"); vardescr.push_back(" GeV");  
   varname.push_back("pmiss"); vartitle.push_back("Missing Momentum"); vardescr.push_back(" GeV/c");  
   varname.push_back("qtot"); vartitle.push_back("Total Charge"); vardescr.push_back("");  
   varname.push_back("wdeltam"); vartitle.push_back("M_{#nu}^{2}, PR"); vardescr.push_back(" GeV^{2}/c^{4}");  
   varname.push_back("wdeltampiz"); vartitle.push_back("M_{#nu}^{2}, PR"); vardescr.push_back(" GeV^{2}/c^{4}");  
   varname.push_back("nks"); vartitle.push_back("N. K_{S}"); vardescr.push_back("");  
   varname.push_back("nkp"); vartitle.push_back("N. K^{#pm}"); vardescr.push_back("");
  
  
  for(UInt_t i = 0; i < varname.size(); i++) {
    
    cout << "Processing " << varname[i].Data() << endl;
    
    MergeRootFiles(dir,varname[i]);
    MakeOverlap(0,dir,varname[i],vartitle[i],vardescr[i]);
    MakeOverlap(1,dir,varname[i],vartitle[i],vardescr[i]);
  }
}


Double_t chisquared(const TH1* h1, const TH1* h2) {

  Double_t chisq (0), tempchisq(0);
  Double_t inter, temperr1, temperr2, inter2;

  inter = h1->Integral();
  inter2 = h2->Integral();
  
  Int_t nbins = h1->GetNbinsX();

  for(Int_t i=1; i<nbins+1; i++) {
    tempchisq = 0;
    temperr1 = 1;
    temperr2 = 1;
    if ( h1->GetBinError(i)>.9 ) temperr1 = h1->GetBinError(i);
    if ( h2->GetBinError(i)>.9 ) temperr2 = h2->GetBinError(i);
    if (((h1->GetBinError(i))*(h1->GetBinError(i))+(h2->GetBinError(i))*(h2->GetBinError(i)))>0) tempchisq = ((h1->GetBinContent(i)) - (h2->GetBinContent(i))*inter/inter2) / sqrt((temperr1*temperr1)+(temperr2*temperr2)*(inter/inter2)*(inter/inter2));  
    chisq = chisq + tempchisq*tempchisq;
  }
  chisq = chisq / nbins;
  return chisq;
}

void MergeRootFiles(const TString& dir,const TString& varname){
  
  cout << "Merging Rootfiles " << endl;
  
  TFile *out = new TFile( TString(dir + varname + ".root"), "RECREATE" );
  TFile *dat = new TFile( TString(dir + varname + "Data.root"), "READ" );
  TFile *mc =  new TFile( TString(dir + varname + "MC.root"), "READ" );

  TDirectory *current_sourcedir;

  Int_t keysum = 0;

  for(Int_t i = 0; i < 2; i++) {
    i == 0 ? dat->cd() : mc->cd();
    cout << "Looping on " << (i == 0 ? dat->GetName() : mc->GetName());
    
    current_sourcedir  = gDirectory;
    
    TIter nextkey( current_sourcedir->GetListOfKeys() );

    TKey *key=0;
    cout << " it has " << current_sourcedir->GetListOfKeys()->GetSize() << " keys " << endl;
    keysum += current_sourcedir->GetListOfKeys()->GetSize();
    
    while ( (key = (TKey*)nextkey())) {
      TObject *obj = key->ReadObj();
      out->cd();
      obj->Write();
      i == 0 ? dat->cd() : mc->cd();
    }
  }
  cout << "Outfile has " << out->GetListOfKeys()->GetSize() << " keys.. ";
  
  if( (out->GetListOfKeys()->GetSize() - keysum) != 0 )
    cout << "... List of keys not matching. Something is badly screwed. Don't trust your plots" << endl;
  else
    cout << " ... Everything matches" << endl;

  out->Close();
}


void MakeOverlap(Int_t cut, const TString& dir, const TString& varname, const TString& vardescr, const TString& units){

  gROOT->SetStyle("Plain");
  cout << "------------> Making overlap plots" << endl;

  TFile *file = new TFile( TString(dir + varname + ".root"),"READ");
  file->cd();

  char line[100];
  
  TString name, mycut;
  TPad *fPads[50];

  TCanvas *c6 = new TCanvas("canvas", "canvas", 300, 0, 400,800);  
  // -- top left
  fPads[1]= new TPad("pad1", "", 0.00, 0.65, 0.99, 0.99);   fPads[1]->Draw(); 
  fPads[2]= new TPad("pad2", "", 0.00, 0.50, 0.99, 0.64);   fPads[2]->Draw();  
  // -- bottom left
  fPads[3]= new TPad("pad3", "", 0.00, 0.15, 0.99, 0.49);   fPads[3]->Draw(); 
  fPads[4]= new TPad("pad4", "", 0.00, 0.00, 0.99, 0.14);   fPads[4]->Draw(); 
  
  cut == 0 ? mycut = "SL" : mycut = "AC";

  name = "data" + mycut + "h";
  TH1D y1;
  ( (TH1D*)gROOT->FindObject(name) )->Copy(y1);

  name = "MC" + mycut + "h";
  TH1D y2;
  ( (TH1D*)gROOT->FindObject(name) )->Copy(y2);

  name = "data" + mycut + "d";
  TH1D y3;
  ( (TH1D*)gROOT->FindObject(name) )->Copy(y3);

  name = "MC" + mycut + "d";
  TH1D y4;
  ( (TH1D*)gROOT->FindObject(name) )->Copy(y4);
  
  fPads[1]->cd(); shrinkPad(0.001, 0.2); 
   
  setFilledHist(&y1, kBlack, kBlue, 3004);
  Double_t inter1;
  Double_t max1 = 1.4*y1.GetMaximum();
  y1.SetLabelSize(0.07, "Y");
  y1.SetMaximum(max1);
  y1.SetMinimum(0.);
  y1.SetMarkerSize(0.5);
  y1.SetName("Enriched");
  inter1 = y1.Integral();
  // y1.SetNormFactor(inter1);
  y1.SetLineColor(kRed);
  name = cut == 0 ? varname + " after Semileptonic cuts: enriched" : varname + " after All cuts: enriched";
  y1.SetTitle(name);
  y1.SetMarkerStyle(20);
  y1.Draw();
  
  Double_t inter2;
  inter2 = y2.Integral();
  y2.SetMarkerSize(0.5);
  
  /* Don't understand this.. maybe it's too old.
    
    if( !cut ) {
    intdataen = inter1;
    intMCen = inter2;
  }

  cout << intdataen << "  " << intMCen << endl;

  if (norm) 
    if( intMCen>0 )
      y2.Scale(intdataen/intMCen);
    else if( inter2>0 )
      y2.Scale(inter1/inter2);
  */

  y2.Scale(inter1/inter2);

  setFilledHist(&y2 , kBlack, kRed, 3005);
  y2.DrawCopy("histosame");
  
  double g1 = chisquared(&y1,&y2);
  sprintf(line, "#chi^{2} = %5.4f", g1);
  
  TLatex tl;
  tl.SetNDC(kTRUE);
  tl.SetTextSizePixels(50); tl.DrawLatex(0.22, 0.75, line);
  
  fPads[2]->cd(); 
  shrinkPad(0.001, 0.2); 
  shrinkPad(0.4, 0.2, 0.1, 0.001);
  
  gPad->SetGridx(1);  gPad->SetGridy(1);
  TH1D *hratio = new TH1D(y1); hratio->SetName("hratio"); hratio->Reset();
  hratio->Divide(&y1, &y2);
  hratio->SetMinimum(0.5); hratio->SetMaximum(1.5);
  hratio->SetMarkerStyle(24);
  hratio->SetNdivisions(504, "Y");
  hratio->GetXaxis()->SetTitle(TString(vardescr + " "+units));
  hratio->SetLabelSize(0.17, "X");
  hratio->SetTitleOffset(1.2,"X");
  hratio->SetTitleSize(0.12,"X");
  hratio->SetLabelSize(0.17, "Y");
  
  hratio->SetStats(0);
  hratio->SetTitle("");
  hratio->Draw();

  //depleted now 
  fPads[3]->cd(); shrinkPad(0.001, 0.2); 

  setFilledHist(&y3, kBlack, kBlue, 3004);
  double inter3;
  double max2 = 1.4*y3.GetMaximum();
  y3.SetLabelSize(0.07, "Y");
  y3.SetMaximum(max2);
  y3.SetMinimum(0.);
  y3.SetMarkerSize(0.5);
  y3.SetName("Depleted");
  inter3 = y3.Integral();
  // y3.SetNormFactor(inter1);
  y3.SetLineColor(kRed);
  name = cut == 0 ? varname + " after Semileptonic cuts: depleted" : varname + " after All cuts: depleted";
  y3.SetTitle(name);
  y3.SetMarkerStyle(20);
  y3.Draw();
 
  double inter4;

  inter4 = y4.Integral();
  y4.SetMarkerSize(0.5);
  
  /*
    if(!cut) {
    intdatadepl = inter3;
    intMCdepl = inter4;
  }
  
  if (norm) 
    y4.Scale(intdatadepl/intMCdepl);
  else
    y4.Scale(inter3/inter4);
  */
  
  y4.Scale(inter3/inter4);

  setFilledHist(&y4 , kBlack, kRed, 3005);
  y4.DrawCopy("histsame");

  double g2 = chisquared(&y3,&y4);
  sprintf(line, "#chi^{2} = %5.4f", g2); 
  tl.SetTextSizePixels(50); tl.DrawLatex(0.22, 0.75, line);
  
  fPads[4]->cd(); shrinkPad(0.001, 0.2); 
  shrinkPad(0.4, 0.2, 0.1, 0.001);

  gPad->SetGridx(1);  gPad->SetGridy(1);
  TH1D *hratio2 = new TH1D(y4); hratio2->SetName("hratio2"); hratio2->Reset();
  hratio2->Divide(&y3, &y4);
  hratio2->SetMinimum(0.5); hratio2->SetMaximum(1.5);
  hratio2->SetMarkerStyle(24);
  hratio2->SetNdivisions(504, "Y");
  hratio2->GetXaxis()->SetTitle(TString(vardescr+ " "+units));
  hratio2->SetLabelSize(0.17, "X");  hratio2->SetLabelSize(0.17, "Y");
  hratio2->SetTitleOffset(1.2,"X");
  hratio2->SetTitleSize(0.12,"X");
  hratio2->SetStats(0);
  hratio2->SetTitle("");
  hratio2->Draw();

  cut == 0 ? mycut = "LeptonCuts" : mycut = "AllCuts";

  name = dir + "/comparison_" + varname + mycut + ".ps";
  
  c6->SaveAs(name);  
}

void shrinkPad(double b, double l, double r, double t) {
  gPad->SetBottomMargin(b); 
  gPad->SetLeftMargin(l);
  gPad->SetRightMargin(r);
  gPad->SetTopMargin(t);
}
void setFilledHist(TH1 *h, Int_t color, Int_t fillcolor, Int_t fillstyle, Int_t width) {
  // Note: 3004, 3005 are crosshatches
  // ----- 1000       is solid
  //       kYellow    comes out gray on bw printers
  h->SetLineColor(color);     h->SetLineWidth(width);   
  h->SetFillStyle(fillstyle); h->SetFillColor(fillcolor);
}
