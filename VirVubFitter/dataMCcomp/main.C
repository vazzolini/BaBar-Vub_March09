#include <iostream>
#include "TFile.h"

#include "TCanvas.h"

#include "mesFit.hh"
#include "util.hh"

#include "TLatex.h"
#include "TPad.h"

using namespace std;

// ===============================================================================================  

void usage();
void Bookhist(Int_t, Double_t, Double_t, Int_t);
void ApplyMixingCorrection(Int_t);
void MakeOverlap(Int_t, const TString&, const TString&, const TString&, const TString&);
void MakeEffPlots(const TString&, const TString&);
Double_t chisquared(const TH1*, const TH1*);

// ===============================================================================================  

int main(int argc, char* argv[]){
  
  gROOT->SetStyle("Plain");
  
// ============== Command Line Arguments =========================================================

  if(argc < 10){
    usage();
    return 1;
  }

  TString rootfilename( argv[1] );
  TString outdir( argv[2] );
  Int_t bch( atoi( argv[3]) );
  Int_t isdata( atoi (argv[4] ));
  TString varname( argv[5] );
  Double_t lowbin( atof(argv[6]) );
  Double_t hibin( atof(argv[7]) );
  Int_t nbin( atoi( argv[8]) );
  TString vardescr( argv[9] );
  TString units( argv[10] );

// ===============================================================================================  

  vector<TVector2*> fitres;
  vector<Double_t> results;
  TString cut,name,ndata,ncut,ndepl,kcut;
  
  TFile *datfile = new TFile( rootfilename.Data() );  

  name = isdata == 1 ? "Data" : "MC";
  
  TFile *fout = new TFile( TString(outdir+"/"+varname+name+".root"),"RECREATE" );    
  RooRealVar *myvar = new RooRealVar( varname.Data(),vardescr.Data(),lowbin,hibin );
   
  RooDataSet *datadata = (RooDataSet*)datfile->Get("DATA");
  RooDataSet *mc = (RooDataSet*)datfile->Get("PStar");
  
  RooArgSet *arg = (RooArgSet*)datadata->get();
  RooRealVar *mes = (RooRealVar*)arg->find("mes");
  
  mesFit *themesfit;

  RooDataSet *reduced = NULL;
  RooDataSet *rds = NULL;

  static Int_t IFLAV_START, IFLAV_STOP;

  switch(bch) {
  case 0: IFLAV_START = 4; IFLAV_STOP = 6; break;
  case 1: IFLAV_START = 3; IFLAV_STOP = 4; break;
  case 2: IFLAV_START = 3; IFLAV_STOP = 6; break;
  default: IFLAV_START = 3; IFLAV_STOP = 6; break;
  }
    		
  Bookhist(nbin,lowbin,hibin,isdata);

  Int_t stop,istart;
  TH1D *h;

  istart = isdata == 1 ? 0 : 1;
  
  //----------------> LOOPS <-------------------

  // DATA/MC LOOP
  for(Int_t idata = istart; idata < istart+1; idata++){
    idata == 0 ? ndata = "data" : ndata = "MC";
    rds = idata == 0 ? datadata : mc;

    // SL/AC cut Loop
    for(Int_t icut = 0; icut < 2; icut++){
      icut == 0 ? ncut = "SL" : ncut = "AC";

      // Enriched / Depleted Loop
      for(Int_t idepl = 0; idepl < 2; idepl++){
	idepl == 0 ? ndepl="h" : ndepl ="d";

	//fix enriched/depleted sample when analyzing nks or nkp;

	if( strcmp(varname.Data(),"nkp") == 0 && idepl == 0) kcut = " && nks == 0";
	if( strcmp(varname.Data(),"nks") == 0 && idepl == 0) kcut = " && nkp == 0";
	
	if( strcmp(varname.Data(),"nkp") == 0 && idepl == 1) kcut = " && nks > 0"; 
	if( strcmp(varname.Data(),"nks") == 0 && idepl == 1) kcut = " && nkp > 0";
	
	// Flav Loop
	for(Int_t iflav = IFLAV_START; iflav < IFLAV_STOP; iflav++){
	  name = ndata + ncut + ndepl + Form("%d",iflav);
	  h = (TH1D*)gROOT->FindObject(name);
	  results.resize(0);
	  
	  for(Int_t isFitBin = 0; isFitBin < 2; isFitBin++) {
	    //Get Best fit parameters before fitting every bin
	    stop = (isFitBin == 0) ? 2 : h->GetNbinsX()+1;
	    
	    // Main loop
	    //	    for(Int_t i = 1; i < stop; i++){
	    for(Int_t i = 1; i < stop; i++){
	      
	      isFitBin == 0 ? cout << "Getting best fit parameters" : cout << "-----> Fitting " << ndata.Data() << " bin " << i;
	      cout << endl;
	      
	      cut = varname + " >= " + Form("%f",h->GetBinLowEdge(i)) + " && " + varname + " < ";
	      isFitBin == 0 ? cut = cut + Form("%f",h->GetBinLowEdge(h->GetNbinsX()) + h->GetBinWidth(h->GetNbinsX())) : 
		cut = cut + Form("%f",h->GetBinLowEdge(i) + h->GetBinWidth(1));
	      
	      cut = cut + " && flavB == " + Form("%d",iflav) + kcut;

	      Int_t mask_sl = 0x3f; //SL MASK enr //= 0111111
	      Int_t mask_ac = 0x3ff; // AC MASK enr //= 01111111111
	      Int_t mask;

	      if( icut == 0 ) { // SL CUT CASE
		mask = mask_sl;
		if( strcmp(varname.Data(),"pcms") == 0 ) mask = 0x3e; //= 0111110
		if( strcmp(varname.Data(),"nle") == 0 ) mask = 0x1f;  //= 0011111
		
		if(idepl) {
		  mask |= 0x40;
		  mask_sl |= 0x40;
		}

		//define again bitmask for Kaon
		
		if( strcmp(varname.Data(),"nks") == 0 || strcmp(varname.Data(),"nkp") == 0 ) {
		  mask = 0x3f;
		  mask_sl = 0x7f;
		}
		
		cut += Form(" && (lepYesBitmask == %d || lepYesBitmask == %d)",mask_sl,mask);
		
	      } else if (icut == 1) {  // ALLCUT CASE
		mask = mask_ac;
		if( strcmp(varname.Data(),"pcms") == 0 ) mask = 0x3fe;          //= 01111111110
		if( strcmp(varname.Data(),"nle") == 0 ) mask = 0x3df;           //= 01111011111
		if( strcmp(varname.Data(),"mm2") == 0  ) mask = 0x3bf;          //= 01110111111
		if( strcmp(varname.Data(),"qtot") == 0  ) mask = 0x37f;         //= 01101111111
		if( strcmp(varname.Data(),"wdeltam") == 0  ) mask = 0x2ff;      //= 01011111111
		if( strcmp(varname.Data(),"wdeltampiz") == 0  ) mask = 0x1ff ;  //= 00111111111
		if( strcmp(varname.Data(),"emiss") == 0 ) mask = 0x3bf;         //= 01110111111
		if( strcmp(varname.Data(),"pmiss") == 0 ) mask = 0x3bf;         //= 01110111111   

		if(idepl) {
		  mask |= 0x400;
		  mask_ac |= 0x400;
		}
		
		if( strcmp(varname.Data(),"nks") == 0 || strcmp(varname.Data(),"nkp") == 0 ) {
		  mask = 0x3ff;
		  mask_ac = 0x7ff;
		} 
		
		cut += Form(" && (AllCutBitmask == %d || AllCutBitmask == %d)",mask_ac,mask);
	      }
	      
	      cout << cut.Data() << endl;
	      
	      reduced = (RooDataSet*)rds->reduce(cut);
	      reduced->SetTitle(cut);
	      
	      themesfit = new mesFit(ArgusAndThosig, mes, fout);
	      name += Form("%s%d","bin",i);
	      if(!isFitBin) name += Form("fitPar");
	      
	      fitres.push_back( themesfit->fitModel(reduced, bool(results.size()),results, name,idata,icut) );
	      
	      if(!isFitBin)
		name.Resize(name.Length() - TString( Form("%s%d","bin",i)).Length() - TString( Form("fitPar") ).Length() );	    
	      else {
		name.Resize(name.Length() - TString( Form("%s%d","bin",i)).Length() );	    
		h->SetBinContent( i,fitres[i-1]->X() );
		h->SetBinError( i,fitres[i-1]->Y() );
	      }

	      // clean up
	      delete themesfit;
	      themesfit = NULL;
	      
	      delete reduced;
	      reduced = NULL;
	    }
	    
	    for(UInt_t jk = 0; jk < fitres.size(); jk++)
	      delete fitres[jk];
	    fitres.resize(0);
	    h->Write(outdir);
	  } //best fit parameters
	} //B flav  loop
      } // depleted loop
    } // SL/AC  loop
  } //data / MC loop


// ========== Apply Mixing Correction=============================================================
  ApplyMixingCorrection(isdata);
    
// =========== Make Plot =========================================================================

//   MakeOverlap(0,outdir,varname,vardescr,units); // SL Cuts 
//   MakeOverlap(1,outdir,varname,vardescr,units); // AC Cuts 
  

// =========== Make EffPlot ======================================================================
  
  //  MakeEffPlots(outdir,varname);

  datfile->Close(); delete datfile;
  fout->Close(); delete fout;
  delete myvar;

  cout << ">---+-------> This is the end <-------+---<" << endl;
  return 0;
}

// ===============================================================================================

void MakeOverlap(Int_t cut, const TString& dir, const TString& varname, const TString& vardescr, const TString& units){

  cout << "------------> Making overlap plots" << endl;

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

// ===============================================================================================

void ApplyMixingCorrection(Int_t isdata){

  cout << "------------> Applying Mixing correction " << endl;
  TH1D *h,*h3,*h4,*h5;
  TString data,cut,depl,name;
  Double_t tempbinchb, tempbinb0os, tempbinb0ss, temperrchb, temperrb0os, temperrb0ss;
  Double_t tempbin, temperr;
  Double_t chid = 0.181;
  
  Int_t istart = isdata == 1 ? 0 : 1;

  for(Int_t idata = istart; idata < istart + 1; idata++) {
    idata == 0 ? data = "data" : data = "MC";
    
    for(Int_t icut = 0; icut < 2; icut ++) {
      icut == 0 ? cut = "SL" : cut = "AC";
      
      for(Int_t idepl = 0; idepl < 2; idepl++) {
	idepl == 0 ? depl="h" : depl ="d";
	
	//	for(Int_t iflav = 3; iflav < 6; iflav++) {
	name = data+cut+depl+Form("%d",3);
	h3 =(TH1D*) gROOT->FindObject(name);
	name = data+cut+depl+Form("%d",4);
	h4 =(TH1D*) gROOT->FindObject(name);
	name = data+cut+depl+Form("%d",5);
	h5 =(TH1D*) gROOT->FindObject(name);
	
	name = data+cut+depl;
	h = (TH1D*) gROOT->FindObject(name);

	for(Int_t i = 1; i < h3->GetNbinsX()+1; i++){
	  tempbinchb = h3->GetBinContent(i);
	  tempbinb0os = h4->GetBinContent(i);
	  tempbinb0ss = h5->GetBinContent(i);
	  temperrchb = h3->GetBinError(i);
	  temperrb0os = h4->GetBinError(i);
	  temperrb0ss = h5->GetBinError(i);
	  
	  tempbin = tempbinchb + ((1-chid)/(1-2*chid)) * tempbinb0os - (chid/(1-2*chid)) * tempbinb0ss;
	  temperr = sqrt(temperrchb*temperrchb + ((1-chid)/(1-2*chid)) * ((1-chid)/(1-2*chid)) * tempbinb0os + (chid/(1-2*chid)) * (chid/(1-2*chid))* tempbinb0ss);
	  
	  h->SetBinContent(i,tempbin);
	  h->SetBinError(i,temperr);
	}
	h->Write();
      }
    }
  }
}

// ===============================================================================================

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

// ===============================================================================================

void MakeEffPlots(const TString& dir,const TString& varname){
  cout << "------------> Making Efficiency plots" << endl;

  TString name;
  TPad *fPads[40];

  name = "dataACh";
  TH1D *y1(((TH1D*)gDirectory->Get(name)));
  name = "dataSLh";
  TH1D *y2(((TH1D*)gDirectory->Get(name)));
  name = "dataACd";
  TH1D *y3(((TH1D*)gDirectory->Get(name)));
  name = "dataSLd";
  TH1D *y4(((TH1D*)gDirectory->Get(name)));
  name = "MCACh";
  TH1D *y5(((TH1D*)gDirectory->Get(name)));
  name = "MCSLh";
  TH1D *y6(((TH1D*)gDirectory->Get(name)));
  name = "MCACd";
  TH1D *y7(((TH1D*)gDirectory->Get(name)));
  name = "MCSLd";
  TH1D *y8(((TH1D*)gDirectory->Get(name)));

  TH1D *hratio1 = new TH1D(*y1); hratio1->SetName("hratio1"); hratio1->Reset();
  TH1D *hratio2 = new TH1D(*y5); hratio2->SetName("hratio2"); hratio2->Reset();
  TH1D *hratio3 = new TH1D(*y3); hratio3->SetName("hratio3"); hratio3->Reset();
  TH1D *hratio4 = new TH1D(*y7); hratio4->SetName("hratio4"); hratio4->Reset();
  TH1D *hratioratio1 = new TH1D(*hratio1); hratioratio1->SetName("hratioratio1"); hratioratio1->Reset();
  TH1D *hratioratio2 = new TH1D(*hratio3); hratioratio2->SetName("hratioratio2"); hratioratio2->Reset();
  
  hratio1->Divide(y1, y2, 1, 1, "B");
  hratio2->Divide(y5, y6, 1, 1, "B");
  hratio3->Divide(y3, y4, 1, 1, "B");
  hratio4->Divide(y7, y8, 1, 1, "B");
  hratioratio1->Divide(hratio1, hratio2);
  hratioratio2->Divide(hratio3, hratio4);

  TCanvas *c7 = new TCanvas("c7", "c7", 300,   0, 400,800);  
  // -- top left
  fPads[1]= new TPad("pad1", "", 0.00, 0.65, 0.99, 0.99);   fPads[1]->Draw(); 
  fPads[2]= new TPad("pad2", "", 0.00, 0.50, 0.99, 0.64);   fPads[2]->Draw();  
  // -- bottom left
  fPads[3]= new TPad("pad3", "", 0.00, 0.15, 0.99, 0.49);   fPads[3]->Draw(); 
  fPads[4]= new TPad("pad4", "", 0.00, 0.00, 0.99, 0.14);   fPads[4]->Draw(); 
  
  fPads[1]->cd(); shrinkPad(0.001, 0.2); 

  //   double max = 1.4*hratio1->GetMaximum();
  hratio1->SetLabelSize(0.07, "Y");
  //   hratio1->SetMaximum(max);
  //cout << max << endl;
  hratio1->SetMaximum(1.);
  hratio1->SetMinimum(-0.1);
  setFilledHist(hratio1 , kBlack, kBlue, 3005);
  hratio1->Draw();

  setFilledHist(hratio2 , kBlack, kRed, 3005);
  hratio2->DrawCopy("histsame");
  
  fPads[2]->cd(); shrinkPad(0.001, 0.2); 
  shrinkPad(0.4, 0.2, 0.1, 0.001);
  
  gPad->SetGridx(1);  gPad->SetGridy(1);
  
  hratioratio1->SetMinimum(0.5); hratioratio1->SetMaximum(1.5);
  hratioratio1->SetMarkerStyle(24);
  hratioratio1->SetNdivisions(504, "Y");
  hratioratio1->SetLabelSize(0.22, "X");  hratioratio1->SetLabelSize(0.17, "Y");
  hratioratio1->SetStats(0);
  hratioratio1->SetTitle("");
  hratioratio1->Draw();
  
  fPads[3]->cd(); shrinkPad(0.001, 0.2); 

  //   double max = 1.4*hratio1->GetMaximum();
  hratio3->SetLabelSize(0.07, "Y");
  //   hratio1->SetMaximum(max);
  //cout << max << endl;
  hratio3->SetMaximum(1.);
  hratio3->SetMinimum(-0.1);
  setFilledHist(hratio3 , kBlack, kBlue, 3005);
  hratio3->Draw();

  setFilledHist(hratio4 , kBlack, kRed, 3005);
  hratio4->DrawCopy("histsame");
  
  fPads[4]->cd(); shrinkPad(0.001, 0.2); 
  shrinkPad(0.4, 0.2, 0.1, 0.001);

  gPad->SetGridx(1);  gPad->SetGridy(1);
  
  hratioratio2->SetMinimum(0.5); hratioratio2->SetMaximum(1.5);
  hratioratio2->SetMarkerStyle(24);
  hratioratio2->SetNdivisions(504, "Y");
  hratioratio2->SetLabelSize(0.22, "X");  hratioratio2->SetLabelSize(0.17, "Y");
  hratioratio2->SetStats(0);
  hratioratio2->SetTitle("");
  hratioratio2->Draw();
  
  name = dir + "/comparisoneff" + varname + ".ps";
  c7->SaveAs(name);  
}

// ===============================================================================================

void Bookhist(Int_t n, Double_t low, Double_t hi, Int_t isdata) {
  
  TH1 *h;
  TString data,cut,depl,name;

  Int_t istart = isdata == 1 ? 0 : 1;

  cout << "BOOKING HISTOGRAMS ";
  cout << (isdata == 1 ? "Data only " : "MC only ");
  cout << endl ;

  for(Int_t idata = istart; idata < istart + 1; idata++) {
    idata == 0 ? data = "data" : data = "MC";
    
    for(Int_t icut = 0; icut < 2; icut ++) {
      icut == 0 ? cut = "SL" : cut = "AC";
      
      for(Int_t idepl = 0; idepl < 2; idepl++) {
	idepl == 0 ? depl = "h" : depl = "d";
	name = data + cut + depl;
	h = new TH1D(name,"",n,low,hi); h->Sumw2();
	
	for(Int_t iflav = 3; iflav < 6 ; iflav++) {
	  name = data + cut + depl + Form("%d",iflav);

	  h = new TH1D(name,"",n,low,hi); h->Sumw2();
	}
      }
    }
  }
}

// ===============================================================================================

void usage() {
  cout << "Please use the proper command line arguments\n" << endl;
  cout << " ./thecomparison rootfilename outdir bch varname lowbin hibin nbins vardescrption " << endl;
  cout << endl;
}


// |_____|___|____|______|_____|__________|____|
// |ksele|nle|flav|intpur|Btype|acceptance|pcms|

// |_____|____________|_________|___________|___|___|____|______|___|______|____|
// |ksele|wdeltacutpiz|wdeltaCut|totalcharge|mm2|nle|flav|intpur|chg|accept|pcms|
