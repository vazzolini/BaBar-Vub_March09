#include "mesFIT.hh"

#include <iostream>
#include <cmath>
#include <stdlib.h>

#include "TF1.h"
#include "TStyle.h"
#include "TMinuit.h"
#include "TVirtualPad.h"  // access to gPad

using namespace std;

//#ifndef __CINT__ 
ClassImp(mesFIT)
//#endif


  double mesFIT_ebeam=-1;

// ----------------------------------------------------------------------
double F_GAUSS(double *x, double *par) {
  // par[0] -> const
  // par[1] -> mean 
  // par[2] -> sigma

  if (par[2] > 0.) {
    Double_t arg = (x[0] - par[1]) / par[2];
    Double_t fitval =  par[0]*TMath::Exp(-0.5*arg*arg);
    return fitval;
  }
  else {
    return -1.;
  }
}



// ----------------------------------------------------------------------
double F_CB(double *x, double *par) {
  // par[0]:  mean
  // par[1]:  sigma
  // par[2]:  alpha, crossover point
  // par[3]:  n, length of tail
  // par[4]:  N, normalization

  Double_t cb = 0.0;
  Double_t exponent = 0.0;

  if (x[0] > par[0] - par[2]*par[1]) {
    exponent = (x[0] - par[0])/par[1];
    cb = TMath::Exp(-exponent*exponent/2.);
  } else {
    double nenner  = TMath::Power(par[3]/par[2], par[3])*TMath::Exp(-par[2]*par[2]/2.);
    double zaehler = (par[0] - x[0])/par[1] + par[3]/par[2] - par[2];
    zaehler = TMath::Power(zaehler, par[3]);
    cb = nenner/zaehler;
  }

  if (par[4] > 0.) {
    cb *= par[4];
  }

  return cb;
}

// ----------------------------------------------------------------------
// Argus only
double F_ARGUS(double *x, double *par) {
  //   par[0] = normalization of argus
  //   par[1] = exponential factor of argus

  // --- If tail is small then Gauss
  double ebeam = mesFIT_ebeam/2;
  double ebeam2 = ebeam*ebeam;
  double background = 0.;
  double x2 = x[0]*x[0];
  if (x2/ebeam2 < 1.) {
    background = par[0]*x[0] * sqrt(1. - (x2/ebeam2)) * exp(par[1] * (1. - (x2/ebeam2))); 
  } else {
    background = 0.;
  }
  return background;
}


// ----------------------------------------------------------------------
// Argus and Gauss 
double F_AAG(double *x, double *par) {
  //   par[0] = normalization of gaussian
  //   par[1] = mean of gaussian
  //   par[2] = sigma of gaussian
  //   par[3] = normalization of argus
  //   par[4] = exponential factor of argus

//    double ebeam = 10.58/2;
//    double signal = 0.;
//    double background = 0.;
//    double result=0.;
//    if (par[2] > 0.)  signal     = par[0] * exp(-(x[0]-par[1]) * (x[0]-par[1]) / (2*par[2]*par[2]));
//    background = par[3] * x[0] * sqrt(1 - (x[0]*x[0])/(ebeam*ebeam)) * exp(par[4] * (1 - (x[0]*x[0])/(ebeam*ebeam))); 
//    result = signal + background;
//    return result;

  return  (F_ARGUS(x, &par[3]) + F_GAUSS(x, &par[0]));

}


// ----------------------------------------------------------------------
double F_AACB(double *x, double *par) {
  //   par[0] = mean of cb
  //   par[1] = sigma of cb
  //   par[2] = alpha
  //   par[3] = n
  //   par[4] = N
  //   par[5] = normalization of argus
  //   par[6] = exponential factor of argus
  return  (F_ARGUS(x, &par[5]) + F_CB(x, &par[0]));
}



// ----------------------------------------------------------------------
mesFIT::mesFIT(int histNameModifier) {
  _histNameModifier=histNameModifier;

  init(0);
}

// ----------------------------------------------------------------------
mesFIT::mesFIT(TH1D *h, const char *f, int setting, int histNameModifier, int algorithm, double en_beam, double arg_uplim) {
  init(0);
  hist(h); 

  _histNameModifier=histNameModifier;

  fVerbose = 0;
  fEbeam = en_beam;
  fSbHi = arg_uplim;

  if ((!strcmp(f, "gaus")) || (!strcmp(f, "Gaus")) || (!strcmp(f, "GAUS"))) {
    fitGauss(31, setting);
  } else if ((!strcmp(f, "gauss")) || (!strcmp(f, "Gauss")) || (!strcmp(f, "GAUSS"))) {
    fitGauss(31, setting);
  } else if ((!strcmp(f, "cb")) || (!strcmp(f, "Cb")) || (!strcmp(f, "CB"))) {
    fitCB(31, setting, algorithm);
  }

  print();
}

// ----------------------------------------------------------------------
void mesFIT::init(int pass) {
  if(pass==0) {
   fR=fRE= fSg =  fSgE = fBg = fBgE = fPur = fPurE = 0;
    fSigma = fSigmaE = fMean = fMeanE = 0;
    
    fVerbose = 0;
    
    fDontChange = kFALSE;
    fPrintLevel = 31;
    fMesMinEvents = 10.;
    
    fMesMin  = 5.20;

    fEbeam=10.5795;
    
    fLo = 5.270;
    fHi = 5.290;
    
    fSbLo = 5.21;
    fSbHi = 5.260;
    
    fPadShrinkX = 0.15; 
    fPadShrinkY = 0.12;
    
    for (int i = 0; i < 10; ++i) {
      fPar[i] = fParLo[i] = fParHi[i] = -9999.;
      fChangePar[i] = fFixPar[i] = fLimitPar[i] = 0;
    }
    
    fTxtSize = 0.08;
    fTxtX    = 0.15;
    fTxtS    = 0.08;
    fLineW   = 2;
    
    fLabelSize = 0.06; 
    fTitleSize = 0.07; 
    fOffsetX = 1.0; 
    fOffsetY = 1.0; 
    
    fMarkerSize = 0.8; 
    fMarkerStyle = 24; 
    
    tl  = new TLatex();
    tl->SetTextSize(fTxtSize);
  }

  if(pass==1) {
    mesFIT_ebeam=fEbeam;
    f0 = new TF1("mesFITf0", F_GAUSS, fLo, fHi, 3);
    f1 = new TF1("mesFITf1", F_CB, fLo, fHi, 5);
    f2 = new TF1("mesFITf2", F_AAG, fMesMin, fHi, 5);
    f3 = new TF1("mesFITf3", F_AACB, fMesMin, fHi, 7); 
    f4 = new TF1("mesFITf4", F_ARGUS, fMesMin, fHi, 2); 
  }
}

// ----------------------------------------------------------------------
void mesFIT::hist(TH1D *h) {
  char histname[255];

  fHist = h;
  fitHist=new TH1D(*h);
  sprintf(histname, "%s_fitHist_%d", h->GetName(), _histNameModifier);
  fitHist->SetName(histname);
  fitHist->Reset();
  fitHist->Sumw2();
}
//---------------------- end constructor related functions --------------------------



//------------------------------- get stuff -----------------------------------------
// get the difference between the histogram and the fit
TH1D *mesFIT::GetResidual(int mode) {
  char histname[255];

  TH1D *diff=(TH1D*)fHist->Clone();

  sprintf(histname, "%s_diff_%d", fHist->GetName(), _histNameModifier);
  diff->SetName(histname);
  diff->Sumw2();
  diff->Add(fitHist, -1.);

  if(mode==1) {
    diff->Divide(fHist);
  }

  //diff->SetMaximum(diff->GetMaximum());
  //diff->SetMinimum(diff->GetMinimum());

  return diff;
}


//------------------------------- set stuff -----------------------------------------
// ----------------------------------------------------------------------
void mesFIT::setTextSize(double x) {
  fTxtSize = x;
  tl->SetTextSize(fTxtSize);
}

// ----------------------------------------------------------------------
void mesFIT::fitGauss(int printLevel, int setting) {
  init(1);

  fPrintLevel = printLevel;
  if (!fDontChange) {
    gStyle->SetOptTitle(0);
    shrinkPad(fPadShrinkX, fPadShrinkY);
    fHist->SetNdivisions(-405, "X");
    setHist(fHist, kBlack, fMarkerStyle, fMarkerSize, fLineW);
    setTitles(fHist, "m_{ES} [GeV]", "Entries / 2.5 MeV", fTitleSize, fOffsetX, fOffsetY, fLabelSize);
  }
  htemp = (TH1D*)fHist->DrawCopy("e");

  f2->SetParameters(0.9*fHist->GetMaximum(), (fLo+fHi)/2, 0.003, 3.*fHist->GetBinContent(3), -10.);

  if (fVerbose > 0) cout << endl << 0.9*fHist->GetMaximum() << endl << endl;

  f2->SetLineStyle(1); f2->SetLineWidth(fLineW); f2->SetLineColor(kBlack);  f2->SetNpx(500);
  if (fHist->GetSumOfWeights() > fMesMinEvents) {
    if (fVerbose > 0) cout <<"  ==> FITTING" << endl;
    fHist->Fit(f2, "RN0", "", fMesMin, fHi);
    ftemp = f2->DrawCopy("same");
  } else {
    if (fVerbose > 0) cout << "not enough entries" << endl;
  }

  // fill the histogram with fit function only. No errors
  fitHist->Eval(f2);
  for(int ii=1; ii<fitHist->GetNbinsX(); ii++) fitHist->SetBinError(ii, 0);

  fMean  = f2->GetParameter(1);   fMeanE  = f2->GetParError(1);
  fSigma = f2->GetParameter(2);   fSigmaE = f2->GetParError(2);

  if (fVerbose > 0) cout << "mesFIT: fMean  = " << fMean << " +/- " << fMeanE << endl;
  if (fVerbose > 0) cout << "mesFIT: fSigma = " << fSigma << " +/- " << fSigmaE << endl;
  
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
  if (fVerbose > 0) cout << "mesFIT: fSg   = " << fSg << " +/- " << fSgE << endl;
    
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
  if (fVerbose > 0) cout << "mesFIT: fBg   = " << fBg << " +/- " << fBgE << endl;

  // -- Sideband 
  f4->SetParameters(a1, a2);
  fSb = f4->Integral(fSbLo, fSbHi)/fHist->GetBinWidth(1);
  f4->SetParameters(a1+a1E, a2);
  max = f4->Integral(fSbLo, fSbHi)/fHist->GetBinWidth(1);
  f4->SetParameters(a1-a1E, a2);
  min = f4->Integral(fSbLo, fSbHi)/fHist->GetBinWidth(1);
  fSbE = (max - min)/2.;
  if (fVerbose > 0) cout << "mesFIT: fSb   = " << fSb << " +/- " << fSbE << endl;

  // define purity
  fPur=fSg/(fBg+fSg);
  // error: 1/fPur=1+fBg/fSg;
  double fPurE_rel=0;
  if(fBg>0 && fSg>0) fPurE_rel=(fBg/fSg)*sqrt((fBgE/fBg)*(fBgE/fBg)+(fSgE/fSg)*(fSgE/fSg))/(1+fBg/fSg);
  fPurE=fPur*fPurE_rel;

  // define sideband scaling factor
  fR=0; fRE=0;
  if(fSb>0) {
    fR=(fBg/fSb);
    if(fBg>0) {
      fRE=(fBg/fSb)*TMath::Sqrt((fSbE/fSb)*(fSbE/fSb) + (fBgE/fBg)*(fBgE/fBg));
    }
  }

  f4->SetParameters(a1, a2);
  f4->SetLineStyle(2); f4->SetLineWidth(fLineW); f4->SetLineColor(kRed);
  ftemp = f4->DrawCopy("same");

  if (printLevel > 0) print();
    
}




// ----------------------------------------------------------------------
void mesFIT::fitCB(int printLevel, int setting, int newAlgorithm) {
  init(1);
  int alg1fix(1);
  // if not specified explicitly (newAlgorithm=999),
  // determine the correct algorithm for this mes distribution
  if(newAlgorithm==999) {
    if(fHist->Integral()>5000) newAlgorithm=2;
    else newAlgorithm=1;
  }

  if (newAlgorithm==3) {
    newAlgorithm = 1;
    alg1fix = 0;
  }
  _newAlgorithm=newAlgorithm;

  cout << "mesFIT::fitCB() starting: (SB=" << fSbLo << "-" << fSbHi << ", SG=" << fLo << "-" << fHi << "), newAlgorithm=" << newAlgorithm << endl;

  fPrintLevel = printLevel;
  int i;
  if (!fDontChange) {
    gStyle->SetOptTitle(0);
    shrinkPad(fPadShrinkX, fPadShrinkY);
    fHist->SetNdivisions(-405, "X");
    setHist(fHist, kBlack, fMarkerStyle, fMarkerSize, fLineW);
    setTitles(fHist, "m_{ES} [GeV]", "Entries / 2.5 MeV", fTitleSize, fOffsetX, fOffsetY, fLabelSize);
  }
  f3->SetLineStyle(1); f3->SetLineWidth(fLineW); f3->SetLineColor(kBlack);  f3->SetNpx(500);
  htemp = (TH1D*)fHist->DrawCopy("e");
  
  // -- if parameters have not been changed manually, set them to default values
  for (i = 0; i < 7; ++i) {if (fChangePar[i] == 1) f3->SetParameter(i, fPar[i]);}
  if (fVerbose > 0) cout << "Setting = " << setting << endl;
  if (setting == 0) {
    if (fChangePar[0] == 0) fPar[0] = (fLo+fHi)/2;
    if (fChangePar[1] == 0) fPar[1] = 0.003;
    if (fChangePar[2] == 0) fPar[2] = 1.;
    if (fChangePar[3] == 0) fPar[3] = 2.;
    if (fChangePar[4] == 0) fPar[4] = 0.9*fHist->GetMaximum();
    if (fChangePar[5] == 0) fPar[5] = 3.*fHist->GetBinContent(3); 
    if (fChangePar[6] == 0) fPar[6] = -10.;
  } else if (setting == 1) {
    if (fChangePar[0] == 0) fPar[0] = (fLo+fHi)/2;
    if (fChangePar[1] == 0) fPar[1] = 0.003;
    if (fChangePar[2] == 0) fPar[2] = 1.;
    if (fChangePar[3] == 0) fPar[3] = 2.;
    if (fChangePar[4] == 0) fPar[4] = 0.9*fHist->GetMaximum();
    if (fChangePar[5] == 0) fPar[5] = 1.*fHist->GetBinContent(3); 
    if (fChangePar[6] == 0) fPar[6] = -10.;
  }
  f3->SetParameters(fPar[0], fPar[1], fPar[2], fPar[3], fPar[4], fPar[5], fPar[6]);
  // -- if parameters have not been limited manually, choose default settings: 
  for (i = 0; i < 7; ++i) {if (fLimitPar[i] == 1)  f3->SetParLimits(i, fParLo[i], fParHi[i]);}
  if (setting == 0) {
    if (fLimitPar[0] == 0)  f3->SetParLimits(0, fLo, fHi); 
    if (fLimitPar[1] == 0)  f3->SetParLimits(1, 0.0024, 0.005); 
    if (fLimitPar[2] == 0)  f3->SetParLimits(2, 0.5, 1.5);
    if (fLimitPar[4] == 0)  f3->SetParLimits(4, 0., 10000000.);
    if (fLimitPar[5] == 0)  f3->SetParLimits(5, 0., 10000000.);
    if (fLimitPar[6] == 0)  f3->SetParLimits(6, -1000., 50.);
  } else if (setting == 1) {
    if (fLimitPar[1] == 0)  f3->SetParLimits(1, 0.0024, 0.008); 
    if (fLimitPar[2] == 0)  f3->SetParLimits(2, 0.5, 1.5);
    if (fLimitPar[4] == 0)  f3->SetParLimits(4, 10., 10000000.);
    if (fLimitPar[5] == 0)  f3->SetParLimits(5, 0., 10000000.);
    if (fLimitPar[6] == 0)  f3->SetParLimits(6, -1000., 50.);
  }
  // -- if parameters have not been fixed manually, choose default settings: 
  for (i = 0; i < 7; ++i) {if (fFixPar[i] == 1)  f3->FixParameter(i, fPar[i]);}
  // fix the p3 if running old algorithm - otherwise fits are VERY ustable
  if (fFixPar[3] == 0 && newAlgorithm==0)  f3->FixParameter(3, 5.);
  if (fVerbose > 0) {
    for (i = 0; i < 6; ++i) cout << "par["<< i << "] = " << fPar[i] << endl;
  }

  // -- Set Errors of bins with zero entries to 1
  for (i = 0; i < fHist->GetNbinsX(); ++i) if (fHist->GetBinContent(i+1) < 0) fHist->SetBinError(i+1, 1.);

  // new algorithm - fixing ARGUS parameters for the combined fit based on the limited-range
  // fit of the sideband to the pure ARGUS function.
  if(newAlgorithm==1 || newAlgorithm==2) {
    cout << "Fit ARGUS" << endl;
    if (fHist->GetSumOfWeights() > fMesMinEvents) {
      f4->SetParameters(3.*fHist->Integral(1, 20)/20., -10.);
      fHist->Fit(f4, "LRN0", "samee", fMesMin, fSbHi);
      f4->GetParameters(fPar_AR);
      fPar_ARE[0]=f4->GetParError(0); 
      fPar_ARE[1]=f4->GetParError(1);
      // force errors to 0 when params are 0: VALERY(050503) - WHY THE F&CK would one want that
      //if(fPar_AR[0]<=0) fPar_ARE[0]=0;
      // if(fPar_AR[1]<=0) fPar_ARE[1]=0;
      f4->Draw("same");
      for(int ipar=0; ipar<2; ipar++) {
	cout << "\tpar[" << ipar << "]=" << fPar_AR[ipar] << endl;
      }
    } else {
      if (fVerbose > 0) cout << "not enough entries" << endl;
    }
    // limit f3 parameters 
    if (alg1fix == 0) {
      f3->SetParLimits(5, fPar_AR[0]-fPar_ARE[0], fPar_AR[0]+fPar_ARE[0]); 
      f3->SetParLimits(6, fPar_AR[1]-fPar_ARE[1], fPar_AR[1]+fPar_ARE[1]); 
    } else {
      f3->FixParameter(5, fPar_AR[0]); f3->FixParameter(6, fPar_AR[1]);
    }
    // reset the f4 X limits
    f4->SetRange(fMesMin, fHi);
  }

  // fit parameters
  double a0(0.), a0E(0.), a1(0.), a1E(0.), a2(0.), a2E(0.), a3(0.), a3E(0.), a4(0.), a4E(0.); 
  // min/max for error estimations
  double min, max;

  if(newAlgorithm!=2) {
    // -- Fit it
    if (fHist->GetSumOfWeights() > fMesMinEvents) {
      fHist->Fit(f3, "LRN0", "samee", fMesMin, fHi);
      f3->GetParameters(fPar);
      f3->Draw("same");
      fchi2=f3->GetChisquare();
    } else {
      if (fVerbose > 0) cout << "not enough entries" << endl;
    }
    
    // fill the histogram with fit function only. No errors
    fitHist->Eval(f3);
    for(int ii=1; ii<fitHist->GetNbinsX(); ii++) fitHist->SetBinError(ii, 0);
    
    fMean  = f3->GetParameter(0);   fMeanE  = f3->GetParError(0);
    fSigma = f3->GetParameter(1);   fSigmaE = f3->GetParError(1);
    if (fVerbose > 0) cout << "mesFIT: fMean  = " << fMean << " +/- " << fMeanE << endl;
    if (fVerbose > 0) cout << "mesFIT: fSigma = " << fSigma << " +/- " << fSigmaE << endl;
    
    a0 = f3->GetParameter(0); a0E = f3->GetParError(0);
    a1 = f3->GetParameter(1); a1E = f3->GetParError(1);
    a2 = f3->GetParameter(2); a2E = f3->GetParError(2);
    a3 = f3->GetParameter(3); a3E = f3->GetParError(3);
    a4 = f3->GetParameter(4); a4E = f3->GetParError(4);
    // a4 sometimes is unstable
    if(a4E>a4) a4E=0.02*a4;
    
    f1->SetParameters(a0, a1, a2, a3, a4);
    fSg = f1->Integral(fLo, fHi)/fHist->GetBinWidth(1);
    f1->SetParameters(a0, a1, a2, a3, a4+a4E);
    max = f1->Integral(fLo, fHi)/fHist->GetBinWidth(1);
    f1->SetParameters(a0, a1, a2, a3, a4-a4E);
    min = f1->Integral(fLo, fHi)/fHist->GetBinWidth(1);
    fSgE = (max - min)/2.;
    if (fVerbose > 0) cout << "mesFIT: fSg   = " << fSg << " +/- " << fSgE << endl;
    if(fSg<0) {fSg=0; fSgE=0;}
  }

  
  // -- Argus background and error estimate
  if(newAlgorithm==1 || newAlgorithm==2) {
    a1 = fPar_AR[0];
    a2 = fPar_AR[1];
    a1E = fPar_ARE[0];
    a2E = fPar_ARE[1];
  }
  if(newAlgorithm==1 && alg1fix==0) {
    a1 = f3->GetParameter(5);
    a2 = f3->GetParameter(6); 
  }
  if(newAlgorithm==0) {
    a1 = f3->GetParameter(5);
    a2 = f3->GetParameter(6); 
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
  if (fVerbose > 0) cout << "mesFIT: fBg   = " << fBg << " +/- " << fBgE << endl;
  if(fBg<0) {fBg=0; fBgE=0;}

  if(newAlgorithm!=2) {
    fSgE = sqrt(fSgE*fSgE + fBgE*fBgE);
  }

  if(newAlgorithm==1||newAlgorithm==0) {
    farg1 = f3->GetParameter(5);
    farg2 = f3->GetParameter(6);
    farg1e = f3->GetParError(5);
    farg2e = f3->GetParError(6);
  }
  if(newAlgorithm==2) {
    farg1 = f4->GetParameter(0);
    farg2 = f4->GetParameter(1);
    farg1e = f4->GetParError(0);
    farg2e = f4->GetParError(1);
  }


  // -- Sideband 
  f4->SetParameters(a1, a2);
  fSb = f4->Integral(fSbLo, fSbHi)/fHist->GetBinWidth(1);
  f4->SetParameters(a1+a1E, a2);
  max = f4->Integral(fSbLo, fSbHi)/fHist->GetBinWidth(1);
  f4->SetParameters(a1-a1E, a2);
  min = f4->Integral(fSbLo, fSbHi)/fHist->GetBinWidth(1);
  fSbE = (max - min)/2.;
  if (fVerbose > 0) cout << "mesFIT: fSb   = " << fSb << " +/- " << fSbE << endl;
  if(fSb<0) {fSb=0; fSbE=0;}


  if(newAlgorithm==2) {
    // -- determine signal from the raw integral only
    int loBin=fHist->FindBin(fLo)+1;
    int hiBin=fHist->FindBin(fHi);
    cout << "mesFIT.cc: " << loBin << " (" << fHist->GetBinLowEdge(loBin) << "), " << hiBin << endl;
    double rawCount=fHist->Integral(loBin, hiBin);
    fSg=rawCount - fBg;
    // error from the histogram count and argus error
    fSgE=sqrt(rawCount + fBgE*fBgE);
  }

  // define sideband scaling factor
  fR=0; fRE=0;
  if(fSb>0) {
    fR=(fBg/fSb);
    if(fBg>0) {
      fRE=(fBg/fSb)*TMath::Sqrt((fSbE/fSb)*(fSbE/fSb) + (fBgE/fBg)*(fBgE/fBg));
    }
    // or get the error by varying parameters of argus function within errors
    // a1 is normalization - no point in varying...
    f4->SetParameters(a1, a2+a2E);
    max = f4->Integral(fLo, fHi)/f4->Integral(fSbLo, fSbHi);
    f4->SetParameters(a1, a2-a2E);
    min = f4->Integral(fLo, fHi)/f4->Integral(fSbLo, fSbHi);
    fRE = std::abs(max - min)/2.;
    if (fVerbose > 0) {
      cout << "mesFIT: fR   = " << fR << " +/- " << fRE << endl;
      cout << "\tparameters varied: a1: " << a1 << ", a2: " << a2 << " +-" << a2E << endl;
    }
  }

  // define purity
  fPur=0;
  if(fBg+fSg>0) {
    fPur=fSg/(fBg+fSg);
  }
  // error: 1/fPur=1+fBg/fSg;
  double fPurE_rel=0;
  if(fBg>0 && fSg>0) fPurE_rel=(fBg/fSg)*TMath::Sqrt((fBgE/fBg)*(fBgE/fBg)+(fSgE/fSg)*(fSgE/fSg))/(1+fBg/fSg);
  fPurE=fPur*fPurE_rel;
  
  f4->SetParameters(a1, a2);
  f4->SetLineStyle(2); f4->SetLineWidth(fLineW); f4->SetLineColor(kRed);
  f4->Draw("same");

  if (printLevel > 0) print();
}  



// ----------------------------------------------------------------------
void mesFIT::print() {
  tl->SetNDC(kTRUE);
  int position=9;
  double step=0.1;
  if (fPrintLevel & 1) {
    sprintf(line, "S = %6.1f #pm %5.1f", fSg, fSgE); 
    tl->DrawLatex(fTxtX, step*(--position), line);
  }
  if (fPrintLevel & 2) {
    sprintf(line, "B = %6.1f #pm %5.1f", fBg, fBgE); 
    tl->DrawLatex(fTxtX, step*(--position), line);
  }
  if(_newAlgorithm!=2) { 
    if (fPrintLevel & 4) {
      sprintf(line, "m = %6.2f #pm %6.2f", 1000.*fMean, 1000.*fMeanE); 
      tl->DrawLatex(fTxtX, step*(--position), line);
    }
    if (fPrintLevel & 8) {
      sprintf(line, "s = %4.2f #pm %4.2f", 1000.*fSigma, 1000.*fSigmaE); 
      tl->DrawLatex(fTxtX, step*(--position), line);
    }
  }
  if (fPrintLevel & 16) {
    sprintf(line, "P = %4.1f #pm %4.1f", 100.*fPur, 100.*fPurE); 
    tl->DrawLatex(fTxtX, step*(--position), line);
  }
  //if (fPrintLevel & 16) {
  //  sprintf(line, "R = %4.1f #pm %4.1f", 100.*fR, 100.*fRE); 
  //  tl->DrawLatex(fTxtX, step*(--position), line);
  //}
  if (fPrintLevel & 16) {
    sprintf(line, "argus = %4.1f #pm %4.1f", farg2, farg2e);
    tl->DrawLatex(fTxtX, step*(--position), line);
  }

}


// ----------------------------------------------------------------------
void mesFIT::shrinkPad(double b, double l, double r, double t) {
  gPad->SetBottomMargin(b); 
  gPad->SetLeftMargin(l);
  gPad->SetRightMargin(r);
  gPad->SetTopMargin(t);
}

// ----------------------------------------------------------------------
void mesFIT::setTitles(TH1 *h, const char *sx, const char *sy, float size, 
               float xoff, float yoff, float lsize, int font) {
  if (h == 0) {
    cout << " Histogram not defined" << endl;
  } else {
    h->SetXTitle(sx);                  h->SetYTitle(sy); 
    h->SetTitleOffset(xoff, "x");      h->SetTitleOffset(yoff, "y");
    h->SetTitleSize(size, "x");        h->SetTitleSize(size, "y");
    h->SetLabelSize(lsize, "x");       h->SetLabelSize(lsize, "y");
    h->SetLabelFont(font, "x");        h->SetLabelFont(font, "y");
    h->GetXaxis()->SetTitleFont(font); h->GetYaxis()->SetTitleFont(font);
    h->SetNdivisions(508, "X");
  }
}

// ----------------------------------------------------------------------
void mesFIT::setHist(TH1 *h, Int_t color, Int_t symbol, Double_t size, Double_t width) {
  h->SetLineColor(color);   h->SetLineWidth(width);
  h->SetMarkerColor(color); h->SetMarkerStyle(symbol);  h->SetMarkerSize(size); 
  h->SetStats(kFALSE); 
  h->SetFillStyle(0); h->SetFillColor(color);
}
