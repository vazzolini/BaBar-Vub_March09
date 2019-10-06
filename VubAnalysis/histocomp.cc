#include <fstream.h>
#include <iostream.h>
#include <string.h>
#include <iomanip.h>
#include <stdlib.h>

#include "histocomp.hh"

#include "TSystem.h"
#include "TRegexp.h"

#include "functions.hh"
#include "../RecoilAnalysis/recoilAnalysis.hh"


// ----------------------------------------------------------------------
histocomp::histocomp() {

  c0 = (TCanvas*)gROOT->FindObject("c0"); 
  if (c0 == 0) {
    cout << "TCanvas c0 not found. Creating my own version." << endl;
    c0 = new TCanvas("c0","--c0--",356,0,656,700);
  }
  tl  = new TLatex();
  tl->SetNDC(kTRUE);
  pl  = new TLine();
  pa  = new TArrow();
  box = new TBox();

  
  fNormalization1 = fNormalization2 = 0.;
  nMc = nDa = 0; 
  fNormArea = 1;

  BMASS   = 5.2792;
  BQMASS  = 4.800;
  A0      = 1.29;


  f0 = new TF1("f0", f_gauss, 5.27, 5.29, 3);
  f1 = new TF1("f1", f_aag, 5.2, 5.29, 5); // main function
  f2 = new TF1("f2", f_aag, 5.2, 5.29, 5); // for bg
  f3 = new TF1("f3", f_p1ag, 0.480, 0.520, 5); // for Ks mass peak


  fc = new TF1("fc", f_cb, 5.2, 5.29, 5); fc->SetNpx(1000);
  fg = new TF1("fg", f_gauss, 5.27, 5.29, 3); fg->SetNpx(1000);
  fa = new TF1("fa", f_argus, 5.2, 5.29, 2); fa->SetNpx(1000);
  fga = new TF1("fga", f_aag, 5.2, 5.29, 5); fga->SetNpx(1000);
  fca = new TF1("fca", f_aacb, 5.2, 5.29, 7); fca->SetNpx(1000);

  fgg = new TF1("fgg", f_2g, -5.0, 5.0, 6); fgg->SetNpx(1000);
  fggg = new TF1("fggg", f_3g, -5.0, 5.0, 9); fggg->SetNpx(1000);
  fp2g = new TF1("fp2g", f_p2ag, -5.0, 5.0, 6); fp2g->SetNpx(1000);

  //  loadFiles();

  int i(-1);
  // lepton spectrum
  i++; sprintf(fHistName[i], "a1000"); sprintf(fMesHist[i], "mesa1");  
  i++; sprintf(fHistName[i], "a1020"); sprintf(fMesHist[i], "mesa21");  
  i++; sprintf(fHistName[i], "e1000"); sprintf(fMesHist[i], "mese1");  
  i++; sprintf(fHistName[i], "e1020"); sprintf(fMesHist[i], "mese21");  
  i++; sprintf(fHistName[i], "m1000"); sprintf(fMesHist[i], "mesm1");   
  i++; sprintf(fHistName[i], "m1020"); sprintf(fMesHist[i], "mesm21");  
  i++; sprintf(fHistName[i], "m1500"); sprintf(fMesHist[i], "mesm1");   
  i++; sprintf(fHistName[i], "m1510"); sprintf(fMesHist[i], "mesm21");   
  i++; sprintf(fHistName[i], "e1500"); sprintf(fMesHist[i], "mese1");   
  i++; sprintf(fHistName[i], "e1510"); sprintf(fMesHist[i], "mese21");   
  // hadron und photon spectrum
  i++; sprintf(fHistName[i], "a1600"); sprintf(fMesHist[i], "mesa1");   
  i++; sprintf(fHistName[i], "a1610"); sprintf(fMesHist[i], "mesa7");   
  i++; sprintf(fHistName[i], "a1700"); sprintf(fMesHist[i], "mesa1");   
  i++; sprintf(fHistName[i], "a1710"); sprintf(fMesHist[i], "mesa7");    
  i++; sprintf(fHistName[i], "a1800"); sprintf(fMesHist[i], "mesa1");   
  i++; sprintf(fHistName[i], "a1810"); sprintf(fMesHist[i], "mesa7");   
  i++; sprintf(fHistName[i], "a1900"); sprintf(fMesHist[i], "mesa1");   
  i++; sprintf(fHistName[i], "a1910"); sprintf(fMesHist[i], "mesa7");    
  // mxrec
  i++; sprintf(fHistName[i], "a2000"); sprintf(fMesHist[i], "mesa1");    
  i++; sprintf(fHistName[i], "a2020"); sprintf(fMesHist[i], "mesa7");   
  // mxfit
  i++; sprintf(fHistName[i], "a2400"); sprintf(fMesHist[i], "mesa1");    
  i++; sprintf(fHistName[i], "a2420"); sprintf(fMesHist[i], "mesa7");   
  i++; sprintf(fHistName[i], "e2400"); sprintf(fMesHist[i], "mese7");   
  i++; sprintf(fHistName[i], "e2420"); sprintf(fMesHist[i], "mese7");   
  i++; sprintf(fHistName[i], "m2400"); sprintf(fMesHist[i], "mesm7");   
  i++; sprintf(fHistName[i], "m2420"); sprintf(fMesHist[i], "mesm7");   
  // theta(pmiss)
  i++; sprintf(fHistName[i], "a2600"); sprintf(fMesHist[i], "mesa1");    
  i++; sprintf(fHistName[i], "a2620"); sprintf(fMesHist[i], "mesa7");   
  // |pmiss|
  i++; sprintf(fHistName[i], "a2700"); sprintf(fMesHist[i], "mesa1");    
  i++; sprintf(fHistName[i], "a2720"); sprintf(fMesHist[i], "mesa7");   
  // mxfit
  i++; sprintf(fHistName[i], "a2800"); sprintf(fMesHist[i], "mesa1");    
  i++; sprintf(fHistName[i], "a2820"); sprintf(fMesHist[i], "mesa7");   
  i++; sprintf(fHistName[i], "a2900"); sprintf(fMesHist[i], "mesa1");    
  i++; sprintf(fHistName[i], "a2920"); sprintf(fMesHist[i], "mesa7");   
  // mm2
  i++; sprintf(fHistName[i], "a3000"); sprintf(fMesHist[i], "mesa1");    
  i++; sprintf(fHistName[i], "a3020"); sprintf(fMesHist[i], "mesa25");   
  // multiplicities
  i++; sprintf(fHistName[i], "a4000"); sprintf(fMesHist[i], "mesa1");    
  i++; sprintf(fHistName[i], "a4020"); sprintf(fMesHist[i], "mesa26");   
  i++; sprintf(fHistName[i], "a4100"); sprintf(fMesHist[i], "mesa1");   
  i++; sprintf(fHistName[i], "a4120"); sprintf(fMesHist[i], "mesa26");   
  // qtot
  i++; sprintf(fHistName[i], "a4300"); sprintf(fMesHist[i], "mesa1");    
  i++; sprintf(fHistName[i], "a4320"); sprintf(fMesHist[i], "mesa25");   
  // nKp, nKs, nLepton
  i++; sprintf(fHistName[i], "a4400"); sprintf(fMesHist[i], "mesa1");    
  i++; sprintf(fHistName[i], "a4420"); sprintf(fMesHist[i], "mesa7");   
  i++; sprintf(fHistName[i], "a4500"); sprintf(fMesHist[i], "mesa1");    
  i++; sprintf(fHistName[i], "a4520"); sprintf(fMesHist[i], "mesa7");    
  i++; sprintf(fHistName[i], "a4600"); sprintf(fMesHist[i], "mesa1");    
  i++; sprintf(fHistName[i], "a4620"); sprintf(fMesHist[i], "mesa7");    
  fX1X2 = i; 

}

// ----------------------------------------------------------------------
histocomp::~histocomp() {
}


// ----------------------------------------------------------------------
void histocomp::x1x2(const char *file1, const char *file2,  const char *dir) {

  makeCanvas(32);
  char hist[100], meshist[100];

  TFile *fFiles[10];
  char label1[200], label2[200];
  fFiles[0] = new TFile(file1);
  fFiles[1] = new TFile(file2);
  sprintf(label1, "Data");
  sprintf(label2, "MC");


  for (int i = 0; i < fX1X2+1; ++i) {
    sprintf(hist, fHistName[i]);
    sprintf(meshist, "sgvub/%s", fMesHist[i]);
    cout << "  hist: " << hist << "  meshist: " << meshist << endl;
    //fillStrings("b->u enh.", label1, label2); 
    theOverlay(fFiles[0], fFiles[1], hist, "vub", meshist, 1); 
    fPads[1]->cd(); sprintf(line, "%s:%s", fFiles[0]->GetName(), hist); tl->SetTextSizePixels(40);  tl->DrawLatex(0.2, 0.92, line);
    sprintf(meshist, "sgvcb/%s", fMesHist[i]);
    cout <<"  hist: " << hist << "  meshist: " << meshist << endl;
    //fillStrings("b->u depl.", label1, label2); 
    theOverlay(fFiles[0], fFiles[1], hist, "vcb", meshist, 3);
    fPads[3]->cd(); sprintf(line, "%s:%s", fFiles[1]->GetName(), hist); tl->SetTextSizePixels(40);  tl->DrawLatex(0.2, 0.92, line);
    sprintf(line, "%s%s%s.eps",dir,"/", hist); 
    cout << line << endl;
    c6->SaveAs(line);
  }
  
}
// ----------------------------------------------------------------------
void histocomp::theOverlay(TFile *file1, TFile *file2, const char *hist, const char *dir, const char *meshist, int nPad){
  gStyle->SetOptStat(0);
  TH1D *h1, *h2; 
  char label[200]; 
  char bla; 
  double nb1(0.), nb2(0.), xmax(-99.), xmin(-99.);
  int num(0), ndivx(0);
  sscanf(hist, "%1s%d", &bla, &num); 
  if ((num >= 1000) && (num <= 1026)) sprintf(label, "p_{cms} [GeV]");
  if ((num >= 1100) && (num <= 1126)) sprintf(label, "p_{cms} [GeV]");
  if ((num >= 1200) && (num <= 1226)) sprintf(label, "p_{cms} [GeV]");
  if ((num >= 1300) && (num <= 1326)) sprintf(label, "cos#theta_{#ell}");
  if ((num >= 1400) && (num <= 1426)) sprintf(label, "p_{#Upsilon} [GeV]");
  if ((num >= 1500) && (num <= 1526)){sprintf(label, "#theta_{lab} [deg]"); ndivx = -606;}
  if ((num >= 1600) && (num <= 1626)) sprintf(label, "p_{lab} [GeV]");
  if ((num >= 1700) && (num <= 1726)) {sprintf(label, "#theta_{lab} [deg]"); ndivx = -606;}
  if ((num >= 1800) && (num <= 1826)) {sprintf(label, "E_{lab} [GeV]"); xmin = 0.; xmax = 1.;}
  if ((num >= 1900) && (num <= 1926)) {sprintf(label, "#theta_{lab} [deg]"); ndivx = -606;}
  if ((num >= 2000) && (num <= 2026)) sprintf(label, "M_{X} [GeV]");
  if ((num >= 2100) && (num <= 2126)) sprintf(label, "M_{X} [GeV]");
  if ((num >= 2200) && (num <= 2226)) sprintf(label, "gen. M_{X} [GeV]");
  if ((num >= 2300) && (num <= 2326)) sprintf(label, "M_{X}**2 [GeV^{2}]");
  if ((num >= 2400) && (num <= 2426)) sprintf(label, "fitted M_{X} [GeV]");
  if ((num >= 2500) && (num <= 2526)) sprintf(label, "fitted M_{X} [GeV]");
  if ((num >= 2600) && (num <= 2626)) sprintf(label, "cos #theta_{miss}");
  if ((num >= 2700) && (num <= 2726)) sprintf(label, "|p_{miss}| [GeV]");
  if ((num >= 2800) && (num <= 2826)) sprintf(label, "fitted M_{X} [GeV]");
  if ((num >= 2900) && (num <= 2926)) sprintf(label, "fitted M_{X} [GeV]");
  if ((num >= 3000) && (num <= 3026)) sprintf(label, "M_{miss}^{2} [GeV^{2}]");
  if ((num >= 3100) && (num <= 3126)) sprintf(label, "Cat 1: M_{miss}^{2} [GeV^{2}]");
  if ((num >= 3200) && (num <= 3226)) sprintf(label, "Cat 2: M_{miss}^{2} [GeV^{2}]");
  if ((num >= 3300) && (num <= 3326)) sprintf(label, "Cat 3: M_{miss}^{2} [GeV^{2}]");
  if ((num >= 3400) && (num <= 3426)) sprintf(label, "Cat 4: M_{miss}^{2} [GeV^{2}]");
  if ((num >= 3500) && (num <= 3526)) sprintf(label, "Cat 5: M_{miss}^{2} [GeV^{2}]");
  if ((num >= 3600) && (num <= 3626)) sprintf(label, "Cat 6: M_{miss}^{2} [GeV^{2}]");
  if ((num >= 4000) && (num <= 4026)) {sprintf(label, "N_{Trk} "); xmin = 0.; xmax = 15.;}
  if ((num >= 4100) && (num <= 4126)) {sprintf(label, "N_{Nut} "); xmin = 0.; xmax = 15.;}
  if ((num >= 4200) && (num <= 4226)) sprintf(label, "Q_{recoil} ");
  if ((num >= 4300) && (num <= 4326)) sprintf(label, "Q_{total} ");
  if ((num >= 4400) && (num <= 4426)) {sprintf(label, "N_{K^{+}} "); xmin = 0.; xmax = 15.;}
  if ((num >= 4500) && (num <= 4526)) {sprintf(label, "N_{K_{S}} "); xmin = 0.; xmax = 15.;}
  if ((num >= 4600) && (num <= 4626)) {sprintf(label, "N_{Lepton} "); xmin = 0.; xmax = 15.;}
  if ((num >= 9000) && (num <= 9026)) sprintf(label, "resolution M_{X} [GeV] ");
  if ((num >= 9100) && (num <= 9126)) sprintf(label, "resolution p_{miss} [GeV] ");
  if ((num >= 9200) && (num <= 9226)) sprintf(label, "resolution #theta_{miss} [rad] ");
  if ((num >= 9300) && (num <= 9326)) sprintf(label, "resolution Q^{2} [GeV^{2}] ");

  if (file1) {
    file1->cd(); 
  } else {
    cout << "histocomp::theOverlay: file1 not defined" << endl;
    return;
  }
  recoilAnalysis rec;
  c0->cd();
   h1 =  rec.bgSubtracted(hist, dir, meshist, 1);  h1->SetName("f1sub");
   TH1D *m1 = (TH1D*)gFile->Get(meshist);
   mesData tempthemes1(*(rec.newMes(m1,1,1)));
   fNormalization1 = tempthemes1.theSig();
   //mesFit a1(m1, "cb");  fNormalization1 = a1.getSig();
   cout << "####> #B normalization, nB1 = " << fNormalization1 << endl;
   if (file2) {
     file2->cd(); 
   } else {
     cout << "histocomp::theOverlay: file2 not defined" << endl;
     return;
   }
   h2 =  rec.bgSubtracted(hist, dir, meshist, 1);  h2->SetName("f2sub");
   TH1D *m2 = (TH1D*)gFile->Get(meshist);
   mesData tempthemes2(*(rec.newMes(m2,1,1)));
   fNormalization2 = tempthemes2.theSig(); 
 //   mesFit a2(m2, "cb");  fNormalization2 = a2.getSig();
   cout << "####> #B normalization, nB1 = " << fNormalization2 << endl;

  if (h1==0 || h2==0) {
    cout << "histocomp::theOverlay: h1 or h2 not defined" << endl;
    return;
  }

  if (fNormArea == 1) {
    cout << "Normalizing to equal area" << endl;
    if (h2->GetSumOfWeights() > 0.) h2->Scale(h1->GetSumOfWeights()/h2->GetSumOfWeights());
  } else {
    cout << "Normalizing to #B" << endl;
    if (fNormalization2) h2->Scale(fNormalization1/fNormalization2);
  }

  setFilledHist(h1, kBlack, kBlue, 3004);
  setFilledHist(h2, kBlack, kRed, 3005);

  double max = (h1->GetMaximum() > h2->GetMaximum()) ? 1.4*h1->GetMaximum() : 1.4*h2->GetMaximum();
  h1->SetLabelSize(0.07, "Y");
  h1->SetMaximum(max);
  h1->SetMinimum(0.1);
  if (ndivx != 0)  h1->SetNdivisions(ndivx, "X");
  if (xmax != -99.) h1->SetAxisRange(xmin, xmax);

  fPads[nPad]->cd();  shrinkPad(0.001, 0.2); 
  char legStyle[10];
  if (!strcmp(fString2.Data(), "SP4") || (!strcmp(fString2.Data(), "SP3"))) {
    sprintf(legStyle, "f");
    h1->Draw("hist");
  } else {
    sprintf(legStyle, "p");
    h1->Draw("e");
  }
  setTitles(h1, "", "Entries/Bin", 0.08, 0.99, 0.75);   
  h2->Draw("histsame");
  h1->Draw("axissame");

  double chi2, ndof; 
  double kolmo = h1->KolmogorovTest(h2);
  sprintf(line, "K = %5.4f", kolmo); tl->SetTextSizePixels(50);  tl->DrawLatex(0.22, 0.83, line);
  double chi2T = chi2Test(h1, h2, chi2, ndof);
  //  sprintf(line, "#chi^{2} = %5.4f", chi2T); tl->SetTextSizePixels(50); tl->DrawLatex(0.22, 0.75, line);
  sprintf(line, "#chi^{2}/df = %3.1f/%3.0f", chi2, ndof); tl->SetTextSizePixels(50); tl->DrawLatex(0.22, 0.75, line);

  ofstream OUT(fTexFile.Data(), ios::app);
  sprintf(line, "\\vdef{%s_%c_cndof%s}{%3.1f/%3.0f}", fFil.Data(), (nPad==1?'u':'c'), hist, chi2, ndof);
  OUT << line << endl;


  leg = new TLegend(0.65,0.6,0.9,0.89);
  leg->SetHeader(fString1);
  leg->SetFillStyle(0); leg->SetBorderSize(0.); leg->SetTextSize(0.08);  leg->SetFillColor(0); 
  legge = leg->AddEntry(h1, fString2, legStyle); 
  legge = leg->AddEntry(h2, fString3, "f"); 
  leg->Draw();

  fPads[nPad+1]->cd(); 
  shrinkPad(0.4, 0.2, 0.1, 0.001);
  //  gPad->SetTopMargin(0.001);  gPad->SetBottomMargin(0.4); gPad->SetLeftMargin(0.2); 
  gPad->SetGridx(1);  gPad->SetGridy(1);
  TH1D *hratio = new TH1D(*h1); hratio->SetName("hratio"); hratio->Reset();
  hratio->Divide(h1, h2);
  hratio->SetMinimum(0.5); hratio->SetMaximum(1.5);
  setTitles(hratio, label, "ratio", 0.2, 1.2, 0.25);   
  hratio->SetMarkerStyle(24);
  hratio->SetNdivisions(504, "Y");
  if (ndivx != 0)  hratio->SetNdivisions(ndivx, "X");
  if (xmax != -99.) hratio->SetAxisRange(xmin, xmax);
  hratio->SetLabelSize(0.22, "X");  hratio->SetLabelSize(0.17, "Y");
  hratio->Draw();

}

// ----------------------------------------------------------------------
void histocomp::makeCanvas(int i) {
  if (i & 64) {
    c7 = new TCanvas("c7", "c7", 200,  20, 800, 800);
    c7->ToggleEventStatus();
  }
  if (i & 32) {
    c6 = new TCanvas("c6", "c6", 300.,   0., 400,800);
    // -- top left
    fPads[1]= new TPad("pad1", "", 0.00, 0.65, 0.99, 0.99);   fPads[1]->Draw(); 
    fPads[2]= new TPad("pad2", "", 0.00, 0.50, 0.99, 0.64);   fPads[2]->Draw(); 
    // -- bottom left
    fPads[3]= new TPad("pad3", "", 0.00, 0.15, 0.99, 0.49);   fPads[3]->Draw(); 
    fPads[4]= new TPad("pad4", "", 0.00, 0.00, 0.99, 0.14);   fPads[4]->Draw(); 
    
  }
  if (i & 16) { 
    c5 = new TCanvas("c5", "c5", 210,   0, 800, 1000);
    c5->ToggleEventStatus();
  }
  if (i & 8) { 
    c4 = new TCanvas("c4", "c4", 210,   0, 800, 600);
    c4->ToggleEventStatus();
  }
  if (i & 4) {
    c3 = new TCanvas("c3", "c3", 200,  20, 800, 800);
    c3->ToggleEventStatus();
  }
  if (i & 1) {
    c1 = new TCanvas("c1", "c1", 20,  60, 1200, 400);
    c1->ToggleEventStatus();
  }
  if (i & 2) { 
    c2 = new TCanvas("c2", "c2", 300, 200, 400, 800);
    c2->ToggleEventStatus();
  }
}


// ----------------------------------------------------------------------
void histocomp::fillStrings(const char *s1, const char *s2, const char *s3, const char *s4, const char *s5) {
  fString1 = TString(s1); 
  fString2 = TString(s2); 
  fString3 = TString(s3); 
  fString4 = TString(s4); 
  fString5 = TString(s5); 
}
