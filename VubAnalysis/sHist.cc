#include <iostream.h>
#include <TCanvas.h>


#include "sHist.hh"
#include "mesFit.hh"


// ----------------------------------------------------------------------
sHist::sHist() {

  fHist = 0; 
  fhSum = 0;
  for (int i = 0; i < BINMAX; ++i) {
    fhHist[i] = 0;
    fVarBins[i] = 0; 
  }
  
}

// ----------------------------------------------------------------------
sHist::sHist(int base, const char *tit, int nbins, double lo, double hi) {
  int *bins = new int[nbins];
  for (int i = 0; i < nbins; ++i) bins[i] = i+1; // FIXME?? no mes histogram for overflow
  bookHist(base, tit, nbins, lo, hi, bins); 
  delete [] bins; 
}

// ----------------------------------------------------------------------
sHist::sHist(int base) {
  bookHist(base);
}


// ----------------------------------------------------------------------
void sHist::bookHist(int base, const char *tit, int nbins, double lo, double hi, int *bins) {

  fBase = base; 

  const int mesbins  = 40; 
  const double meslo = 5.20; 
  const double meshi = 5.30; 

  fMax  = hi;
  fNmes = nbins; 

  if (fNmes > BINMAX) {
    cout << "too many mes bins requested" << endl;
    return;
  }

  TString title = TString(tit) + TString(Form("; varbin = %d", fNmes)); 

  fHist  = new TH1D(Form("a%d", fBase), title, nbins, lo, hi); fHist->Sumw2(); 
  char line[200];
  char name[200];
  fVarBins[0] = 1; 
  for (int i = 0; i < fNmes; ++i) {
    if (bins) fVarBins[i] = bins[i]; 
    sprintf(line, "sH%dvarbin%d %s: %d ..  (%f .. ) ", fBase, fVarBins[i], tit, fVarBins[i], fHist->GetBinLowEdge(fVarBins[i])); 
    //    cout << line << endl;
    sprintf(name, "a%dvarbin%d", fBase, fVarBins[i]); 
    cout << name << " .. " << line << " .. " << mesbins << " .. " << meslo << " .. " << meshi << endl;
    fhHist[i] = new TH1D(name, line, mesbins, meslo, meshi); fhHist[i]->Sumw2(); 
  }

  sprintf(name, "a%dsum", fBase); 
  sprintf(line, "sH%dsum ", fBase); 
  cout << name << " .. " << line << " .. " << mesbins << " .. " << meslo << " .. " << meshi << endl;
  fhSum = new TH1D(name, line, mesbins, meslo, meshi); fhSum->Sumw2(); 
}


// ----------------------------------------------------------------------
void sHist::bookHist(int base) {

  cout << "Reading histogram h" << base << endl;
  fBase = base; 

  char name[200];
  sprintf(name, "a%d", fBase); 
  fhSum = (TH1D*)gDirectory->Get(name); // temporary abuse of hsum. 

  cout << "bookHist: " 
       << fhSum->GetName() << " " << fhSum->GetTitle() << " " << fhSum->GetNbinsX() << " " << fhSum->GetBinLowEdge(1) << endl;
  fNmes = fhSum->GetNbinsX();
  
  fVarBins[0] = 1; 
  for (int i = 0; i < fNmes; ++i) {
    sprintf(name, "a%dvarbin%d", fBase, i+1); 
    //     cout << "Reading Bin " << name << endl;
    fhHist[i] = (TH1D*)gDirectory->Get(name);
  }

  fHist  = new TH1D(Form("a%d", fBase), 
		    fhSum->GetTitle(), 
		    fNmes,
		    fhSum->GetBinLowEdge(1),
		    fhSum->GetBinLowEdge(fhSum->GetNbinsX()+1)
		    ); fHist->Sumw2();

  // -- now fill hsum with the real sum of all mes histograms
  sprintf(name, "a%dsum", fBase); 
  fhSum = (TH1D*)gDirectory->Get(name);


  //   cout << "Created histogram " << Form("a%d", fBase) 
  //        << " " << fhSum->GetTitle() 
  //        << " " << fNmes
  //        << " " << fhSum->GetBinLowEdge(1)
  //        << " " << fhSum->GetBinLowEdge(fhSum->GetNbinsX()+1)
  //        << endl;

}



// ----------------------------------------------------------------------
void sHist::setup(double *var, double *mes) {
  fpVar = var;
  fpMes = mes; 
}


// ----------------------------------------------------------------------
void sHist::fillHist(double w8) {
  int vbin(-1), tmp(-1);

  double var = *fpVar; 
  double mes = *fpMes; 

  tmp = fHist->FindBin(var);
  for (int i = 0; i < fNmes; ++i) {
    if (tmp >= fVarBins[i]) vbin = i; 
  }

  //  cout << "var = " << var << " bin = " << tmp << "  -> vbin = " << vbin << endl;

  if ((vbin > -1) && (var < fMax)) {
    fhHist[vbin]->Fill(mes, w8); 
    fhSum->Fill(mes, w8);
  }

}


// ----------------------------------------------------------------------
TH1D* sHist::signalHist(const char *func, int prefit) {

  //  mesFit *mesf[BINMAX];
  mesFit mesf[BINMAX];


  TCanvas *c0 = (TCanvas*)gROOT->FindObject("c0");
  c0->Clear();
  c0->Divide(4,3);
  c0->cd(1);

  cout <<"Fitting sum start" << endl;
  mesFit a(fhSum, func, 63, prefit); 
  cout <<"Fitting sum end" << endl;


  cout << "signalHist: " 
       << fhSum->GetName() << " " << fhSum->GetTitle() << " " << fhSum->GetNbinsX() << " " << fhSum->GetBinLowEdge(1) << endl;

  double s(0.), sE(0.);
  for (int i = 0; i < fNmes; ++i) {
    c0->cd(i+2);
    mesf[i].hist(fhHist[i]);
    mesf[i].fixSignal(a); 
    mesf[i].fitCB(63); 

    s  = mesf[i].getSig();
    sE = mesf[i].getSigE();

    // -- Not much statistics in mes peak. Two reasons: mes histogram as a whole is low statistics
    //    or the signal-constrained fit just provided no signal component. 
    if (s < 1.) {
      if (fhHist[i]->GetSumOfWeights() < 40.) {
	// -- Not enough entries in mES histogram for stable fits
	cout << "sHist: " << fBase << " Warning, low statistics, resetting to integral in bin " 
	     << i+1 << " at abscissa value " << fHist->GetBinCenter(i+1) << endl;
	s = fhHist[i]->Integral(fhHist[i]->FindBin(5.27), fhHist[i]->FindBin(5.29));
	sE = TMath::Sqrt(s);
      } else {
	// -- Signal-constrained fit has no signal ...
	sE = 1.;
	sE = mesf[i].getBgE(); // this is a fishy argument: The BG could is uncertain by bgE, 
	                       // therefore the same amount  could be signal
      }
    }

    if (sE < 0.99*TMath::Sqrt(s)) {
      cout << "sHist: " << fBase << " Warning, error smaller than sqrt(entries), resetting in bin " 
	   << i+1 << " at abscissa value " << fHist->GetBinCenter(i+1) 
	   << " to sqrt(sigE**2 + bgE**2)"
	   << endl;
      sE = TMath::Sqrt(mesf[i].getSigE()*mesf[i].getSigE() + mesf[i].getBgE()*mesf[i].getBgE());

//       s  = fhHist[i]->Integral(fhHist[i]->FindBin(5.27), fhHist[i]->FindBin(5.29)) - mesf[i].getBg();
//       cout << "sHist: " << fBase << " Warning, resetting bin content " 
// 	   << i+1 << " at abscissa value " << fHist->GetBinCenter(i+1) 
// 	   << " from " << mesf[i].getSig() 
// 	   << " to " << s
// 	   << endl;


      //      delete mesf[i];
    }

    cout << "sHist: SetBinContent(" << i+1 << ", " << s << ")  at " << fhHist[i]->GetTitle() << endl;
    fHist->SetBinContent(i+1, s); 
    fHist->SetBinError(i+1, sE); 
  }

  for (int i = 0; i <= fHist->GetNbinsX(); ++i) {
    cout << Form("%2d %3.3f: %5.3f +/- %5.3f", i, fHist->GetBinCenter(i), 
		 fHist->GetBinContent(i), fHist->GetBinError(i)) << endl;
  }


  c0->cd(9);
  fHist->Draw();

  return fHist;
}

// ----------------------------------------------------------------------
void sHist::print() {
  cout << "sHist " << fBase << endl;

  for (int i = 0; i < fNmes; ++i)  {
    if (fhHist[i]) (fhHist[i])->Print();
  }

  cout << "end printout for " << fBase << endl;
}
