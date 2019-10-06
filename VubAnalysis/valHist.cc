#include "valHist.hh"

#include <iostream.h>


// ----------------------------------------------------------------------
valHist::valHist() {
  fBase = fNmes = 0; 
  fpVar = fpMes = 0; 
  fpSignalBox = fpSideband = fpNoCuts = fpSigCuts = fpAoCuts = 0; 
  a00 = a10 = a20 = a00sg = a10sg = a20sg = a00bg = a10bg = a20bg = 0; 
  for (int i = 0; i < MESMAX; ++i) {
    a00mes[i] = a10mes[i] = a20mes[i] = 0; 
    fVarBins[i] = 0; 
  }
}

// ----------------------------------------------------------------------
valHist::valHist(int base, const char *tit, int nbins, double lo, double hi, int nmes, int *bins) {
  bookHist(base, tit, nbins, lo, hi, nmes, bins); 
}

// ----------------------------------------------------------------------
valHist::valHist(int base, const char *tit, int nbins, double lo, double hi) {
  int *bins = new int[nbins];
  for (int i = 0; i < nbins; ++i) bins[i] = i+1; // no mes histogram for underflow
  bookHist(base, tit, nbins, lo, hi, nbins, bins); 
  delete [] bins; 
}

// ----------------------------------------------------------------------
valHist::valHist(int base, TH1D *zb, int nmes, int *bins) {
  bookHist(base, zb->GetTitle(), zb->GetNbinsX(), zb->GetBinLowEdge(1), zb->GetBinLowEdge(zb->GetNbinsX()+1), nmes, bins);
}

// ----------------------------------------------------------------------
void valHist::bookHist(int base, const char *tit, int nbins, double lo, double hi, int nmes, int *bins) {
  fBase = base; 

  const int mesbins  = 40; 
  const double meslo = 5.20; 
  const double meshi = 5.30; 

  fNmes = nmes; 
  if (nmes > MESMAX) {
    cout << "too many mes bins requested" << endl;
    return;
  }

  TString title = TString(tit) + TString(Form("; varbin = %d", fNmes)); 

  a00   = new TH1D(Form("a%d",   fBase),    title, nbins, lo, hi); a00->Sumw2(); 
  a10   = new TH1D(Form("a%d",   fBase+10), title, nbins, lo, hi); a10->Sumw2(); 
  a20   = new TH1D(Form("a%d",   fBase+20), title, nbins, lo, hi); a20->Sumw2(); 
					    
  a00sg = new TH1D(Form("a%dsg", fBase),    title, nbins, lo, hi); a00sg->Sumw2(); 
  a10sg = new TH1D(Form("a%dsg", fBase+10), title, nbins, lo, hi); a10sg->Sumw2(); 
  a20sg = new TH1D(Form("a%dsg", fBase+20), title, nbins, lo, hi); a20sg->Sumw2(); 
					    	    
  a00bg = new TH1D(Form("a%dbg", fBase),    title, nbins, lo, hi); a00bg->Sumw2(); 
  a10bg = new TH1D(Form("a%dbg", fBase+10), title, nbins, lo, hi); a10bg->Sumw2(); 
  a20bg = new TH1D(Form("a%dbg", fBase+20), title, nbins, lo, hi); a20bg->Sumw2(); 

  char line[200];
  char name[200];
  fVarBins[0] = 1; 
  for (int i = 0; i < fNmes; ++i) {
    if (bins) fVarBins[i] = bins[i]; 
    sprintf(line, "a%dvarbin%d %s: %d ..  (%f .. ) ", fBase, fVarBins[i], tit, fVarBins[i], a00->GetBinLowEdge(fVarBins[i])); 
    //    cout << line << endl;
    sprintf(name, "a%dvarbin%d", fBase, fVarBins[i]);    a00mes[i] = new TH1D(name, line, mesbins, meslo, meshi); a00mes[i]->Sumw2(); 
    sprintf(name, "a%dvarbin%d", fBase+10, fVarBins[i]); a10mes[i] = new TH1D(name, line, mesbins, meslo, meshi); a10mes[i]->Sumw2(); 
    sprintf(name, "a%dvarbin%d", fBase+20, fVarBins[i]); a20mes[i] = new TH1D(name, line, mesbins, meslo, meshi); a20mes[i]->Sumw2(); 
  }
}


// ----------------------------------------------------------------------
void valHist::print() {
  cout << "Base = " << fBase << endl;
}


// ----------------------------------------------------------------------
void valHist::setup(double *var, double *mes, Bool_t *signalBox, Bool_t *sideband, 
		    Bool_t *nocuts, Bool_t *sigcuts, Bool_t *aocuts
		    ) {

  fpVar = var; 
  fpMes = mes; 
  fpSignalBox = signalBox; 
  fpSideband = sideband; 
  fpNoCuts  = nocuts; 
  fpSigCuts = sigcuts; 
  fpAoCuts  = aocuts; 
}


// ----------------------------------------------------------------------
void valHist::fillHist(double w8) {
  int vbin(0), tmp(-1);

  double var = *fpVar; 

  tmp = a00sg->FindBin(var);
  for (int i = 0; i < fNmes; ++i) {
    if (tmp >= fVarBins[i]) vbin = i; 
  }

  //  cout << "var = " << var << " bin = " << tmp << "  -> vbin = " << vbin << endl;

  //    char line[200];
  //    sprintf(line, "  ncuts %d acuts %d  aocuts%d mes %4.3f  sgbox %d  sideb %d var: %5.4f %s", 
  //  	  *fpNoCuts, *fpSigCuts, *fpAoCuts, *fpMes, *fpSignalBox,*fpSideband, var, a00sg->GetTitle()); 
  //    cout << line << endl;

  if (*fpNoCuts) {
    a00->Fill(var, w8); 
    if (*fpSignalBox) a00sg->Fill(var, w8); 
    if (*fpSideband)  a00bg->Fill(var, w8); 
    a00mes[vbin]->Fill(*fpMes, w8); 
  }

  //  if (*fpSigCuts) {
  if ((*fpSigCuts) && (*fpAoCuts)) {
    a10->Fill(var, w8); 
    if (*fpSignalBox) a10sg->Fill(var, w8); 
    if (*fpSideband)  a10bg->Fill(var, w8); 
    a10mes[vbin]->Fill(*fpMes, w8); 
  }

  if (*fpAoCuts) {
    a20->Fill(var, w8); 
    if (*fpSignalBox) a20sg->Fill(var, w8); 
    if (*fpSideband)  a20bg->Fill(var, w8); 
    a20mes[vbin]->Fill(*fpMes, w8); 
  }
}


// ----------------------------------------------------------------------
void valHist::fillHist(int ntracks, float *tracks, int *goodTracks, double w8) {

  int vbin(0), tmp(-1);

  for (int itrack = 0; itrack < ntracks; ++itrack) {
    if (goodTracks[itrack] == 0 ) {
      //      cout << "Skipping track " << itrack << endl;
      continue;
    }
    //    cout << "Taking   track " << itrack << endl;
    double var = tracks[itrack]; 

    tmp = a00sg->FindBin(var);
    for (int i = 0; i < fNmes; ++i) {
      if (tmp >= fVarBins[i]) vbin = i; 
    }

    if (*fpNoCuts) {
      a00->Fill(var, w8); 
      if (*fpSignalBox) a00sg->Fill(var, w8); 
      if (*fpSideband)  a00bg->Fill(var, w8); 
      a00mes[vbin]->Fill(*fpMes, w8); 
    }

    //    if (*fpSigCuts) {
    if ((*fpSigCuts) && (*fpAoCuts)) {
      a10->Fill(var, w8); 
      if (*fpSignalBox) a10sg->Fill(var, w8); 
      if (*fpSideband)  a10bg->Fill(var, w8); 
      a10mes[vbin]->Fill(*fpMes, w8); 
    }

    if (*fpAoCuts) {
      a20->Fill(var, w8); 
      if (*fpSignalBox) a20sg->Fill(var, w8); 
      if (*fpSideband)  a20bg->Fill(var, w8); 
      a20mes[vbin]->Fill(*fpMes, w8); 
    }
  }
}
