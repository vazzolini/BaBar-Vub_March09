#include <fstream.h>
#include <iostream.h>
#include <iomanip.h>
#include <math.h>
#include <stdlib.h>

#include "TF1.h"
#include "TStyle.h"
#include "TMinuit.h"

#include "mesFit.hh"
#include "functions.hh"
#include "util.hh"

ClassImp(mesFit)

// ----------------------------------------------------------------------
mesFit::~mesFit() {
  if (f0) {
    delete f0; 
  } else {
    cout << "f0 not present?" << endl;
  }
  if (f1) {
    delete f1; 
  } else {
    cout << "f1 not present?" << endl;
  }
  if (f2) {
    delete f2; 
  } else {
    cout << "f2 not present?" << endl;
  }

  if (f3) {
    delete f3; 
  } else {
    cout << "f3 not present?" << endl;
  }

  if (f4) {
    delete f4; 
  } else {
    cout << "f4 not present?" << endl;
  }

  if (tl) {
    delete tl; 
  } else {
    cout << "tl not present?" << endl;
  }

}


// ----------------------------------------------------------------------
mesFit::mesFit() {
  init();
}


// ----------------------------------------------------------------------
mesFit::mesFit(const char *s, const char *f, int verbose, int preFit) {
  TH1D *h = (TH1D*)gDirectory->Get(s); 
  if (h) {
    init();
    hist(h); 
    fVerbose = verbose; 
    
    if ((!strcmp(f, "gaus")) || (!strcmp(f, "Gaus")) || (!strcmp(f, "GAUS"))
        || (!strcmp(f, "gauss")) || (!strcmp(f, "Gauss")) || (!strcmp(f, "GAUSS"))) {
      fitGauss(fPrintLevel);
      print();
    } else if ((!strcmp(f, "cb")) || (!strcmp(f, "Cb")) || (!strcmp(f, "CB"))
               || (!strcmp(f, "xb")) || (!strcmp(f, "Xb")) || (!strcmp(f, "XB"))) {
      fitCB(fPrintLevel, preFit);
      print();
    }
    
  } else {
    cout << "Error: " << s << " not found in gDirectory" << endl;
  }
}


// ----------------------------------------------------------------------
mesFit::mesFit(TH1D *h, const char *f, int verbose, int preFit) {
  init();
  hist(h); 
  fVerbose = verbose; 

  if ((!strcmp(f, "gaus")) || (!strcmp(f, "Gaus")) || (!strcmp(f, "GAUS"))
      || (!strcmp(f, "gauss")) || (!strcmp(f, "Gauss")) || (!strcmp(f, "GAUSS"))) {
    fitGauss(fPrintLevel);
    print();
  } else if ((!strcmp(f, "cb")) || (!strcmp(f, "Cb")) || (!strcmp(f, "CB"))
	     || (!strcmp(f, "xb")) || (!strcmp(f, "Xb")) || (!strcmp(f, "XB"))) {
    fitCB(fPrintLevel, preFit);
    print();
  }

}

// ----------------------------------------------------------------------
mesFit::mesFit(TH1D *h, const mesFit &fix, const char *f, int verbose) {
  init();
  hist(h); 
  fVerbose = verbose; 

  fixSignal(fix); 

  if ((!strcmp(f, "gaus")) || (!strcmp(f, "Gaus")) || (!strcmp(f, "GAUS"))
      || (!strcmp(f, "gauss")) || (!strcmp(f, "Gauss")) || (!strcmp(f, "GAUSS"))) {
    fitGauss(fPrintLevel);
    print();
  } else if ((!strcmp(f, "cb")) || (!strcmp(f, "Cb")) || (!strcmp(f, "CB"))
	     || (!strcmp(f, "xb")) || (!strcmp(f, "Xb")) || (!strcmp(f, "XB"))) {
    fitCB(fPrintLevel, 0);
    print();
  }

}


// ----------------------------------------------------------------------
void mesFit::fixSignal(const mesFit &a) {
  fixPar(0, a.getPar(0));
  fixPar(1, a.getPar(1));
  fixPar(2, a.getPar(2));
  fixPar(3, a.getPar(3));
}

// ----------------------------------------------------------------------
void mesFit::fixSignal(double peak, double sigma, double xover, double tail) {
  fixPar(0, peak);
  fixPar(1, sigma);
  fixPar(2, xover);
  fixPar(3, tail);
}



// ----------------------------------------------------------------------
void mesFit::limitSignal(const mesFit &a) {
  limitPar(0, a.getPar(0) - a.getParErr(0), a.getPar(0) + a.getParErr(0));
  limitPar(1, a.getPar(1) - a.getParErr(1), a.getPar(1) + a.getParErr(1));
  limitPar(2, a.getPar(2) - a.getParErr(2), a.getPar(2) + a.getParErr(2));
  limitPar(3, a.getPar(3) - a.getParErr(3), a.getPar(3) + a.getParErr(3));
}


// ----------------------------------------------------------------------
void mesFit::limitSignal(double peak, double peakE, double sigma, double sigmaE, 
			 double xover, double xoverE, double tail, double tailE) {
  limitPar(0, peak  - peakE,   peak  + peakE);
  limitPar(1, sigma - sigmaE,  sigma + sigmaE);
  limitPar(2, xover - xoverE,  xover + xoverE);
  limitPar(3, tail  - tailE,   tail  + tailE);
}


// ----------------------------------------------------------------------
void mesFit::releaseAllPar() {
  for (int ipar = 0; ipar < 10; ++ipar) {
    fParLo[ipar] = fParHi[ipar] = fChangePar[ipar] = fFixPar[ipar] = fLimitPar[ipar] = 0;
  }
}


// ----------------------------------------------------------------------
void mesFit::init() {
  fSg =  fSgE = fBg = fBgE = fPur = fPurE = fR = fRE = 0.;
  fSigma = fSigmaE = fMean = fMeanE = 0.;

  fVerbose = 1;

  fDontChange = kFALSE;
  fPrintLevel = 63;
  fMesMinEvents = 40.;

  fMesMin  = 5.200;
  
  fLo = 5.270;
  fHi = 5.290;

  fSbLo = 5.200;
  fSbHi = 5.260;

  fPadShrinkX = 0.15; 
  fPadShrinkY = 0.15;

  fTitleX = TString("m_{ES} [GeV]"); 
  fTitleY = TString("Entries / 2.5 MeV"); 

  for (int i = 0; i < 10; ++i) {
    fPar[i] = fParErr[i] = fParLo[i] = fParHi[i] = -9999.;
    fChangePar[i] = fFixPar[i] = fLimitPar[i] = 0;
  }

  for (int i = 0; i < 2; ++i) {
    fParArgus[i] = fParArgusE[i] = -9999.;
  }

  fTxtSize = 0.08;
  fTxtX    = 0.19;
  fTxtS    = 0.08;
  fLineW   = 2;

  fLabelSize = 0.06; 
  fTitleSize = 0.07; 
  fOffsetX = 1.0; 
  fOffsetY = 1.0; 

  fNdivX = -405;
  fNdivY = 306;

  fMarkerSize = 0.6; 
  fMarkerStyle = 24; 

  tl  = new TLatex();
  tl->SetTextSize(fTxtSize);
  
  TString fname = uniqueName("mesFitf0");
  f0 = new TF1(fname.Data(), f_gauss, 5.27, 5.29, 3);

  fname = uniqueName("mesFitf1");
  f1 = new TF1(fname.Data(), f_cb, 5.27, 5.29, 5);

  fname = uniqueName("mesFitf2");
  f2 = new TF1(fname.Data(), f_aag, 5.2, 5.29, 5);

  fname = uniqueName("mesFitf3");
  f3 = new TF1(fname.Data(), f_aacb, 5.2, 5.29, 7); 

  fname = uniqueName("mesFitf4");
  f4 = new TF1(fname.Data(), f_argus, 5.2, 5.29, 2); 
}

// ----------------------------------------------------------------------
TString mesFit::uniqueName(const char *name) {
  TString fname(name);
  TF1 *f = (TF1*)gROOT->GetListOfFunctions()->FindObject(fname.Data());
  int cnt(0);
  while (f) {
    fname = TString(name) + TString(Form("v%d", cnt));
    f = (TF1*)gROOT->GetListOfFunctions()->FindObject(fname.Data());
    ++cnt;
  }
  return fname;
}

// ----------------------------------------------------------------------
void mesFit::hist(TH1D *h) {
  fHist = h;
}


// ----------------------------------------------------------------------
void mesFit::setTextSize(double x) {
  fTxtSize = x;
  tl->SetTextSize(fTxtSize);
}

// ----------------------------------------------------------------------
void mesFit::fitGauss(int printLevel) {
  fPrintLevel = printLevel;
  if (!fDontChange) {
    gStyle->SetOptTitle(0);
    shrinkPad(fPadShrinkX, fPadShrinkY);
    setHist(fHist, kBlack, fMarkerStyle, fMarkerSize, fLineW);
    setTitles(fHist, fTitleX.Data(), fTitleY.Data(), fTitleSize, fOffsetX, fOffsetY, fLabelSize);
  }
  fHist->SetNdivisions(fNdivX, "X");
  fHist->SetNdivisions(fNdivY, "Y");
  htemp = (TH1D*)fHist->DrawCopy("e");

  f2->SetParameters(0.9*fHist->GetMaximum(), 5.28, 0.003, 3.*fHist->GetBinContent(3), -10.);

  f2->SetLineStyle(1); f2->SetLineWidth(fLineW); f2->SetLineColor(kBlack);  f2->SetNpx(500);
  if (fHist->GetSumOfWeights() > fMesMinEvents) {
    cout <<"  ==> FITTING" << endl;
    fHist->Fit(f2, "RN0", "", fMesMin, fHi);
    ftemp = f2->DrawCopy("same");
  } else {
    cout << "not enough entries" << endl;
  }

  int i(0); 
  for (i = 0; i < 5; ++i) {
    fPar[i] = f2->GetParameter(i); 
    fParErr[i] = f2->GetParError(i); 
  }

  fMean  = f2->GetParameter(1);   fMeanE  = f2->GetParError(1);
  fSigma = f2->GetParameter(2);   fSigmaE = f2->GetParError(2);

  // Gaussian normalization and +/- 1 sigma 
  double a1(0.), a1E(0.), a2(0.), a2E(0.), a3(0.), a3E(0.); 
  a1 = f2->GetParameter(0); a1E = f2->GetParError(0);
  a2 = f2->GetParameter(1); a2E = f2->GetParError(1);
  a3 = f2->GetParameter(2); a3E = f2->GetParError(2);
  f0->SetParameters(a1, a2, a3);
  fSg = f0->Integral(fLo, fHi)/fHist->GetBinWidth(1);
  f0->SetParameters(a1+a1E, a2, a3);
  double max = f0->Integral(fLo, fHi)/fHist->GetBinWidth(1);
  f0->SetParameters(a1-a1E, a2, a3);
  double min = f0->Integral(fLo, fHi)/fHist->GetBinWidth(1);
  fSgE = (max - min)/2.;
  
  // -- Argus background and error estimate
  a1 = f2->GetParameter(3); a1E = f2->GetParError(3);
  a2 = f2->GetParameter(4); a2E = f2->GetParError(4);
  
  f4->SetParameters(a1, a2);
  fBg = f4->Integral(fLo, fHi)/fHist->GetBinWidth(1);
  f4->SetParameters(a1+a1E, a2);
  max = f4->Integral(fLo, fHi)/fHist->GetBinWidth(1);
  f4->SetParameters(a1-a1E, a2);
  min = f4->Integral(fLo, fHi)/fHist->GetBinWidth(1);
  fBgE = (max - min)/2.;

  // -- Sideband 
  f4->SetParameters(a1, a2);
  fSb = f4->Integral(fSbLo, fSbHi)/fHist->GetBinWidth(1);
  f4->SetParameters(a1+a1E, a2);
  max = f4->Integral(fSbLo, fSbHi)/fHist->GetBinWidth(1);
  f4->SetParameters(a1-a1E, a2);
  min = f4->Integral(fSbLo, fSbHi)/fHist->GetBinWidth(1);
  fSbE = (max - min)/2.;

  fPur  = 100.*fSg/(fSg+fBg); 
  fPurE = 100.*dRatio(fSg, fSgE, fSg+fBg, TMath::Sqrt(fSgE*fSgE + fBgE*fBgE)); 

  if (fSb > 0.) {
    fR  = fBg/fSb;
    min = (fBg - fBgE) / (fSb - fSbE);
    max = (fBg + fBgE) / (fSb + fSbE);
    fRE = 0.5 * TMath::Abs(min - max);
  }

  f4->SetParameters(a1, a2);
  f4->SetLineStyle(2); f4->SetLineWidth(fLineW); f4->SetLineColor(kRed);
  ftemp = f4->DrawCopy("same");

  if (fVerbose > 0) {
    sprintf(line, "===============    \"%s\" with title \"%s\"     ================", 
	    fHist->GetName(), fHist->GetTitle());    cout << line << endl;
    sprintf(line, "# Signal         # Bg         R       Argus Param.     Peak [MeV]   Sigma [MeV]"); cout << line << endl;
    //    sprintf(line, "-------------------------------------------------------------------------------"); cout << line << endl;
    sprintf(line, "%6.1f+/-%4.1f  %6.1f+/-%4.1f  %4.3f   %4.2f+/-%4.2f   %4.1f+/-%3.1f  %4.1f+/-%3.1f", 
	    fSg, fSgE, fBg, fBgE, fR,  
	    fArgus, fArgusE, 1000.*fMean, 1000.*fMeanE, 1000.*fSigma, 1000.*fSigmaE); cout << line << endl; 
  }
  if (fPrintLevel > 0) print();
    
}




// ----------------------------------------------------------------------
void mesFit::fitCB(int printLevel, int preFit) {
  fPrintLevel = printLevel;
  int i;
  if (!fDontChange) {
    gStyle->SetOptTitle(0);
    shrinkPad(fPadShrinkX, fPadShrinkY);
    setHist(fHist, kBlack, fMarkerStyle, fMarkerSize, fLineW);
    setTitles(fHist, fTitleX.Data(), fTitleY.Data(), fTitleSize, fOffsetX, fOffsetY, fLabelSize, 132);
  }
  fHist->SetNdivisions(fNdivX, "X");
  fHist->SetNdivisions(fNdivY, "Y");
  f3->SetLineStyle(1); f3->SetLineWidth(fLineW); f3->SetLineColor(kBlack);  
  f3->SetFillStyle(0); f3->SetFillColor(1);
  f3->SetNpx(500);
  htemp = (TH1D*)fHist->DrawCopy("e");

  // -- Set Errors of bins with zero entries to 1
  int cnt(0);
  for (i = 0; i < fHist->GetNbinsX(); ++i) {
    if (fHist->GetBinContent(i+1) < 1.e-15) {
      cnt++;
      fHist->SetBinError(i+1, 1.);
    }
  }
  if (cnt > 3) cout << "mesFit: Reset bin error to 1 in " << cnt << " bins" << endl;


  // -- Fresh start
  for (i = 0; i < 7; ++i) f3->ReleaseParameter(i); 
  // -- Parameter settings
  if (fChangePar[0] == 0) fPar[0] = 5.28;
  if (fChangePar[1] == 0) fPar[1] = 0.003;
  if (fChangePar[2] == 0) fPar[2] = 1.;
  if (fChangePar[3] == 0) fPar[3] = 3.;
  if (fChangePar[4] == 0) fPar[4] = 0.9*fHist->GetMaximum();
  //  if (fChangePar[5] == 0) fPar[5] = 3.*fHist->GetBinContent(3); // old setting: 1. * ...
  int nbin = fHist->FindBin(fSbHi) - fHist->FindBin(fSbLo);
  if (fChangePar[5] == 0) fPar[5] = 3.*fHist->Integral(fHist->FindBin(fSbLo), fHist->FindBin(fSbHi))/nbin;
  if (fChangePar[6] == 0) fPar[6] = -10.;
  f3->SetParameters(fPar[0], fPar[1], fPar[2], fPar[3], fPar[4], fPar[5], fPar[6]);
  // -- Parameter limits
  for (i = 0; i < 7; ++i) {
    if (fLimitPar[i] == 1) f3->SetParLimits(i, fParLo[i], fParHi[i]);
  }
  // -- Fixed Parameters
  for (i = 0; i < 7; ++i) {
    if (fFixPar[i] == 1) f3->FixParameter(i, fPar[i]);
  }
  // ??? REMNANT ??? 
  if ((fFixPar[3] == 0) && (preFit == 0))  f3->FixParameter(3, 5.);

  if (fVerbose >= 2) {
    cout << "PRE FIT" << endl;
    for (i = 0; i < 7; ++i) {
      cout << " par["       << i << "] = " << fPar[i] 
	   << " parErr["    << i << "] = " << fParErr[i]
	   << " parChange[" << i << "] = " << fChangePar[i]
	   << " parLimit["  << i << "] = " << fLimitPar[i]
	   << " parFix["    << i << "] = " << fFixPar[i]
	   << endl;
    }
  }

  // -- Fixing ARGUS parameters for the combined fit based on sideband
  //    fit to the pure ARGUS function.
  if (preFit == 1) {
    cout << "mesFit:  prefitting ARGUS in sideband" << endl;
    if (fHist->GetSumOfWeights() > fMesMinEvents) {
      int nbin = fHist->FindBin(fSbHi) - fHist->FindBin(fSbLo) + 1;
      f4->SetParameters(3.*fHist->Integral(fHist->FindBin(fSbLo), fHist->FindBin(fSbHi))/nbin, -10.);
      if (fVerbose < 2) {
	fHist->Fit(f4, "LRN0Q", "samee", fMesMin, fSbHi);
      }	else {
	fHist->Fit(f4, "LRN0", "samee", fMesMin, fSbHi);
      }	
      fParArgus[0]=f4->GetParameter(0); 
      fParArgus[1]=f4->GetParameter(1);
      fParArgusE[0]=f4->GetParError(0); 
      fParArgusE[1]=f4->GetParError(1);
    } else {
      if (fVerbose > 0) cout << "mesFit: not enough entries for prefitting" << endl;
    }
    f3->FixParameter(5, fParArgus[0]); f3->FixParameter(6, fParArgus[1]);
    // -- reset the x-axis limits
    f4->SetRange(fMesMin, fHi);
  }



  // -- Fit it
  if (fHist->GetSumOfWeights() > fMesMinEvents) {
    if (fVerbose < 2) {
      fHist->Fit(f3, "LRN0Q", "samee", fMesMin, fHi);
    } else {
      fHist->Fit(f3, "LRN0", "samee", fMesMin, fHi);
    }
    ftemp = f3->DrawCopy("same");
  } else {
    cout << "mesFit: not enough entries for fitting" << endl;
    return;
  }

  // -- Get parameters and errors 
  for (i = 0; i < 7; ++i) {
    fPar[i] = f3->GetParameter(i); 
    fParErr[i] = f3->GetParError(i); 
  }

  if (fVerbose >= 2) {
    cout << "POST FIT" << endl;
    for (i = 0; i < 7; ++i) {
      cout << " par["       << i << "] = " << fPar[i] 
	   << " parErr["    << i << "] = " << fParErr[i]
	   << " parChange[" << i << "] = " << fChangePar[i]
	   << " parLimit["  << i << "] = " << fLimitPar[i]
	   << " parFix["    << i << "] = " << fFixPar[i]
	   << endl;
    }
  }

  
  fMean  = f3->GetParameter(0);   fMeanE  = f3->GetParError(0);
  fSigma = f3->GetParameter(1);   fSigmaE = f3->GetParError(1);
  fArgus = f3->GetParameter(6);   fArgusE = f3->GetParError(6); 

  double a0(0.), a0E(0.), a1(0.), a1E(0.), a2(0.), a2E(0.), a3(0.), a3E(0.), a4(0.), a4E(0.); 
  a0 = f3->GetParameter(0); a0E = f3->GetParError(0);
  a1 = f3->GetParameter(1); a1E = f3->GetParError(1);
  a2 = f3->GetParameter(2); a2E = f3->GetParError(2);
  a3 = f3->GetParameter(3); a3E = f3->GetParError(3);
  a4 = f3->GetParameter(4); a4E = f3->GetParError(4);

  f1->SetParameters(a0, a1, a2, a3, a4);
  fSg = f1->Integral(fLo, fHi)/fHist->GetBinWidth(1);
  f1->SetParameters(a0, a1, a2, a3, a4+a4E);
  double max = f1->Integral(fLo, fHi)/fHist->GetBinWidth(1);
  f1->SetParameters(a0, a1, a2, a3, a4-a4E);
  double min = f1->Integral(fLo, fHi)/fHist->GetBinWidth(1);
  fSgE = (max - min)/2.;
  
  // -- Argus background and error estimate in signalbox
  a1 = f3->GetParameter(5); a1E = f3->GetParError(5);
  a2 = f3->GetParameter(6); a2E = f3->GetParError(6);

  if (preFit == 1) {
    a1E = fParArgusE[0];
    a2E = fParArgusE[1];
  } else {
    a1E = f3->GetParError(5);
    a2E = f3->GetParError(6);
  }
  
  f4->SetParameters(a1, a2);
  fBg = f4->Integral(fLo, fHi)/fHist->GetBinWidth(1);
  f4->SetParameters(a1+a1E, a2);
  max = f4->Integral(fLo, fHi)/fHist->GetBinWidth(1);
  f4->SetParameters(a1-a1E, a2);
  min = f4->Integral(fLo, fHi)/fHist->GetBinWidth(1);
  fBgE = (max - min)/2.;

  // -- Sideband 
  f4->SetParameters(a1, a2);
  fSb = f4->Integral(fSbLo, fSbHi)/fHist->GetBinWidth(1);
  f4->SetParameters(a1+a1E, a2);
  max = f4->Integral(fSbLo, fSbHi)/fHist->GetBinWidth(1);
  f4->SetParameters(a1-a1E, a2);
  min = f4->Integral(fSbLo, fSbHi)/fHist->GetBinWidth(1);
  fSbE = (max - min)/2.;


  fPur  = 100.*fSg/(fSg+fBg); 
  fPurE = 100.*dEff(fSg, fSgE, fSg+fBg, TMath::Sqrt(fSgE*fSgE + fBgE*fBgE)); 

  if (fSb > 0.) {
    fR  = fBg/fSb;
    min = (fBg - fBgE) / (fSb - fSbE);
    max = (fBg + fBgE) / (fSb + fSbE);
    fRE = 0.5 * TMath::Abs(min - max);
  }

  f4->SetParameters(a1, a2);
  f4->SetLineStyle(2); f4->SetLineWidth(fLineW); f4->SetLineColor(kRed);
  f4->SetFillStyle(0); f4->SetFillColor(1);
  ftemp = f4->DrawCopy("same");

  if (preFit == 1) {
    f4->SetLineColor(kBlue);
    f4->SetRange(fSbLo, fSbHi);
    ftemp = f4->DrawCopy("same");
  }    

  if (fVerbose > 0) {
    sprintf(line, "===============    \"%s\" with title \"%s\", Entries = %.1f     ================", 
	    fHist->GetName(), fHist->GetTitle(), fHist->GetSumOfWeights());    cout << line << endl;
    sprintf(line, "# Signal         # Bg         R       Argus Param.     Peak [MeV]   Sigma [MeV]"); cout << line << endl;
    //    sprintf(line, "-------------------------------------------------------------------------------"); cout << line << endl;
    sprintf(line, "%6.1f+/-%4.1f  %6.1f+/-%4.1f  %4.3f   %4.2f+/-%4.2f   %4.1f+/-%3.1f  %4.1f+/-%3.1f", 
	    fSg, fSgE, fBg, fBgE, fR,  
	    fArgus, fArgusE, 1000.*fMean, 1000.*fMeanE, 1000.*fSigma, 1000.*fSigmaE); cout << line << endl; 
  }
  if (fPrintLevel > 0) print();
}  

// ----------------------------------------------------------------------
void mesFit::dump() {
    sprintf(line, "===============    \"%s\" with title \"%s\"     ================", 
	    fHist->GetName(), fHist->GetTitle());    cout << line << endl;
    sprintf(line, "# Signal         # Bg         R       Argus Param.     Peak [MeV]   Sigma [MeV]"); cout << line << endl;
    //    sprintf(line, "-------------------------------------------------------------------------------"); cout << line << endl;
    sprintf(line, "%6.1f+/-%4.1f  %6.1f+/-%4.1f  %4.3f   %4.2f+/-%4.2f   %4.1f+/-%3.1f  %4.1f+/-%3.1f", 
	    fSg, fSgE, fBg, fBgE, fR,  
	    fArgus, fArgusE, 1000.*fMean, 1000.*fMeanE, 1000.*fSigma, 1000.*fSigmaE); cout << line << endl; 
}

// ----------------------------------------------------------------------
void mesFit::print() {
  tl->SetNDC(kTRUE);
  if (fPrintLevel & 1) {
    sprintf(line, "S = %6.1f #pm %5.1f", fSg, fSgE); 
    tl->DrawLatex(fTxtX, 0.8, line);
  }
  if (fPrintLevel & 2) {
    sprintf(line, "B = %6.1f #pm %5.1f", fBg, fBgE); 
    tl->DrawLatex(fTxtX, 0.8-1*1.1*fTxtSize, line);
  }
  if (fPrintLevel & 4) {
    if (fLimitPar[0] == 1) {
      sprintf(line, "m = %6.2f", 1000.*fMean); 
    } else if (fFixPar[0] == 1) {
      sprintf(line, "m #equiv %6.2f", 1000.*fMean); 
    } else {
      sprintf(line, "m = %6.2f #pm %6.2f", 1000.*fMean, 1000.*fMeanE); 
    }
    tl->DrawLatex(fTxtX, 0.8-2*1.1*fTxtSize, line);
  }
  if (fPrintLevel & 8) {
    if (fLimitPar[1] == 1) {
      sprintf(line, "s = %4.2f", 1000.*fSigma); 
    } else  if (fFixPar[1] == 1) {
      sprintf(line, "s #equiv %4.2f", 1000.*fSigma); 
    } else {
      sprintf(line, "s = %4.2f #pm %4.2f", 1000.*fSigma, 1000.*fSigmaE); 
    }
    tl->DrawLatex(fTxtX, 0.8-3*1.1*fTxtSize, line);
  }
  if (fPrintLevel & 16) {
    //    sprintf(line, "P = %4.1f #pm %4.1f", 100.*fSg/(fSg+fBg), 100.*dRatio(fSg, fSgE, fSg+fBg, TMath::Sqrt(fSgE*fSgE + fBgE*fBgE))); 
    sprintf(line, "P = %4.1f #pm %4.1f", fPur, fPurE); 
    tl->DrawLatex(fTxtX, 0.8-4*1.1*fTxtSize, line);
  }
  if (fPrintLevel & 32) {
    sprintf(line, "R = %5.4f", getR()); 
    tl->DrawLatex(fTxtX, 0.8-5*1.1*fTxtSize, line);
  }

}

