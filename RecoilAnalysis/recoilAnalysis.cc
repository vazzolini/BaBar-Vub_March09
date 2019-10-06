#include "recoilAnalysis.hh"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

ClassImp(recoilAnalysis)

// ----------------------------------------------------------------------
double gauss(double *x, double *par) {
  //   par[0] = normalization 
  //   par[1] = mean 
  //   par[2] = sigma
  double arg = 0;
  if (par[2]) arg = (x[0] - par[1])/par[2];
  return par[0]*TMath::Exp(-0.5*arg*arg);
}

// ----------------------------------------------------------------------
double argus(double *x, double *par) {
  //   par[0] = normalization of gaussian
  //   par[1] = mean of gaussian
  //   par[2] = sigma of gaussian
  //   par[3] = normalization of argus
  //   par[4] = exponential factor of argus

  double signal = 0.;
  double background = 0.;

  double ebeam = 10.58/2;
  double ebeam2 = ebeam*ebeam;
  double result = 0.;
  double x2 = x[0]*x[0];

  if (par[2] > 0.)  signal     = par[0] * exp(-(x[0]-par[1]) * (x[0]-par[1]) / (2*par[2]*par[2]));
  if (x2/ebeam2 < 1.) background = par[3] * x[0] * sqrt(1 - x2/(ebeam*ebeam)) * exp(par[4] * (1 - x2/(ebeam*ebeam))); 
  result = signal + background;
  return result;
}

// ----------------------------------------------------------------------
double argusOnly(double *x, double *par) {
  //   par[0] = normalization of argus
  //   par[1] = exponential factor of argus

  double background = 0.;

  double ebeam = 10.58/2;
  double ebeam2 = ebeam*ebeam;
  double x2 = x[0]*x[0];

  if (x2/ebeam2 < 1.) background = par[0] * x[0] * sqrt(1 - x2/(ebeam*ebeam)) * exp(par[1] * (1 - x2/(ebeam*ebeam))); 
  return background;
}

// ----------------------------------------------------------------------
double crystalball(double *x, double *par) {
  // par[0]:  mean
  // par[1]:  sigma
  // par[2]:  alpha, crossover point
  // par[3]:  n, length of tail
  // par[4]:  N, normalization

  Double_t cb = 0.0;
  Double_t exponent = 0.0;
  Double_t bla = 0.0;

  if (x[0] > par[0] - par[2]*par[1]) {
    exponent = (x[0] - par[0])/par[1];
    cb = TMath::Exp(-exponent*exponent/2.);
  } else {
    //      bla = (par[0] -x[0])/par[1] + par[3]/par[2] - par[2];
    //      exponent = TMath::Power(bla, par[3]);
    //      bla = par[3]/par[2];
    //      cb = (TMath::Power(bla, par[3]) * TMath::Exp(-par[2]*par[2]/2.)) / exponent;

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
double argusAndGauss(double *x, double *par) {
  //   par[0] = normalization of gaussian
  //   par[1] = mean of gaussian
  //   par[2] = sigma of gaussian

  //   par[3] = normalization of argus
  //   par[4] = exponential factor of argus

  double a = argusOnly(x, &par[3]);
  double g = gauss(x, &par[0]);

  return a+g; 
}

// ----------------------------------------------------------------------
double cbAndArgus(double *x, double *par) {
  //   par[0] = mean of cb
  //   par[1] = sigma of cb
  //   par[2] = alpha
  //   par[3] = n
  //   par[4] = N
  
  //   par[5] = normalization of argus
  //   par[6] = exponential factor of argus

  double a  = argusOnly(x, &par[5]);
  double cb = crystalball(x, &par[0]);

  return a+cb; 
}



// ----------------------------------------------------------------------
recoilAnalysis::recoilAnalysis() {
  
  nMc = nDa = nTs = nTsc = 0; 
  for(int jn = 0; jn <MAXSIZE; jn++) {
    for(int jo = 0; jo <105; jo++) {
      fflagC[jn][jo] = 0;
      fflag[jn][jo] = 0;
    }
  }
  fMesMinEvents = 0.;

  f0 = new TF1("f0", "gaus", 5.27, 5.29);
  f1 = new TF1("f1", argus, 5.2, 5.29, 5); // main function
  f2 = new TF1("f2", argus, 5.2, 5.29, 5); // for bg


  fc = new TF1("fc", crystalball, 5.2, 5.29, 5); fc->SetNpx(1000);
  fg = new TF1("fg", gauss, 5.27, 5.29, 3); fg->SetNpx(1000);
  fa = new TF1("fa", argusOnly, 5.2, 5.29, 2); fa->SetNpx(1000);
  fga = new TF1("fga", argusAndGauss, 5.2, 5.29, 5); fga->SetNpx(1000);
  fca = new TF1("fca", cbAndArgus, 5.2, 5.29, 7); fca->SetNpx(1000);

  //  loadFiles();

}

// ----------------------------------------------------------------------
recoilAnalysis::~recoilAnalysis() {
}


// ----------------------------------------------------------------------
void recoilAnalysis::loadMc(const char *name, double lumi) {
  if (nMc > 19) {
    cout << "Too many open MC files. Increase nMc. " << endl;
    return;
  } 
  fMc[nMc] = new TFile(name);
  fMcLumi[nMc] = lumi;
  ++nMc; 
}

// ----------------------------------------------------------------------
void recoilAnalysis::loadDa(const char *name, double lumi) {
  if (nDa > 19) {
    cout << "Too many open DATA files. Increase nDa. " << endl;
    return;
  } 
  fDa[nDa] = new TFile(name);
  fDaLumi[nDa] = lumi;
  ++nDa; 
}

// ----------------------------------------------------------------------
void recoilAnalysis::loadTs(const char *name) {
  if (nTs > 104) {
    cout << "Too many open timestamp files. Increase nTs. " << endl;
    return;
  } 
  fTs[nTs] = new TObjString(name);
  ++nTs; 
}

// ----------------------------------------------------------------------
void recoilAnalysis::loadTsC(const char *name) {
  if (nTsc > 104) {
    cout << "Too many open timestamp files. Increase nTsc. " << endl;
    return;
  } 
  fTsc[nTsc] = new TObjString(name);
  ++nTsc; 
}


// ----------------------------------------------------------------------
void recoilAnalysis::lsTs() { 
  for (Int_t i = 0; i < nTs; ++i) {
    cout << "fTs[" << i << "]: " << fTs[i]->GetString() << endl;
  }
}


// ----------------------------------------------------------------------
void recoilAnalysis::lsTsc() { 
  for (Int_t i = 0; i < nTsc; ++i) {
    cout << "fTsc[" << i << "]: " << fTsc[i]->GetString() << endl;
  }
}


// ----------------------------------------------------------------------
void recoilAnalysis::loadTimestamp(const char *filename) {
  ifstream is(filename);
  char buffer[200];
  char type[100];
  char file[100];
  while (is.getline(buffer, 200, '\n')) {
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %s", type, file);
    if (!strcmp(type, "ts")) {
      //      cout << "Loading MC file   " << file << endl;
      loadTs(file);
    }
    else if (!strcmp(type, "co")) {
      //      cout << "Loading data file " << file << endl;
      loadTsC(file);
    }
    else {
      //
    }
  }
  lsTs();
  lsTsc();
  //  lsPur();
}

// ----------------------------------------------------------------------
void recoilAnalysis::loadTimestampC(const char *filename) {
  ifstream is(filename);
  char buffer[200];
  char file[100];
  while (is.getline(buffer, 200, '\n')) {
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s", file);
    loadTsC(file);
  }
  lsTsc();
}

// ----------------------------------------------------------------------
void recoilAnalysis::loadTimestampI(const char *filename) {
  ifstream is(filename);
  char buffer[200];
  char file[100];
  while (is.getline(buffer, 200, '\n')) {
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s", file);
    loadTs(file);
  }
  lsTs();
}

// ----------------------------------------------------------------------
void recoilAnalysis::lsMc() { 
  for (Int_t i = 0; i < nMc; ++i) {
    cout << "fMc[" << i << "]: " << fMc[i]->GetName() << endl;
  }
}

// ----------------------------------------------------------------------
void recoilAnalysis::lsDa() { 
  for (Int_t i = 0; i < nDa; ++i) {
    cout << "fDa[" << i << "]: " << fDa[i]->GetName() << endl;
  }
}


// ----------------------------------------------------------------------
void recoilAnalysis::loadFiles(const char *filename) {
  ifstream is(filename);
  char buffer[200];
  char type[100];
  char file[100];
  double lumi;
  while (is.getline(buffer, 200, '\n')) {
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    lumi = -1.;
    sscanf(buffer, "%s %s %f", type, file, &lumi);
    if (!strcmp(type, "mc")) {
      //      cout << "Loading MC file   " << file << endl;
      loadMc(file);
    }
    else if (!strcmp(type, "da")) {
      //      cout << "Loading data file " << file << endl;
      loadDa(file);
    }
    else {
      //
    }
  }
  lsMc();
  lsDa();
}


// ----------------------------------------------------------------------
void recoilAnalysis::setHist(TH1 *h, int color, int symbol, float size, int width) {
  if (!h) return;
  h->SetLineColor(color);   h->SetLineWidth(width);
  h->SetMarkerColor(color); h->SetMarkerStyle(symbol);  h->SetMarkerSize(size); 
  h->SetStats(kFALSE); 
  h->SetFillStyle(0); h->SetFillColor(color);
}

// ----------------------------------------------------------------------
void recoilAnalysis::setFilledHist(TH1 *h,  int color, int fillcolor, int fillstyle, int width) {
  if (!h) return;
  h->SetLineColor(color);     
  h->SetLineWidth(width);   
  h->SetFillStyle(fillstyle); 
  h->SetFillColor(fillcolor);
}


// ----------------------------------------------------------------------
void recoilAnalysis::setTitles(TH1 *h, const char *sx, const char *sy, float size, float xoff, float yoff, float lsize, int font) {
  h->SetXTitle(sx);                  
  h->SetYTitle(sy); 
  h->SetTitleOffset(xoff, "x");      
  h->SetTitleOffset(yoff, "y");
  h->SetTitleSize(size, "x");        
  h->SetTitleSize(size, "y");
  h->SetLabelSize(lsize, "x");       
  h->SetLabelSize(lsize, "y");
}



// ----------------------------------------------------------------------
double recoilAnalysis::dBinomial(double n, double N) {
  if ((N <= 0) || (n < 0.)) return 0.;
  double w = n/N;
  return TMath::Sqrt(TMath::Abs(w*(1-w)/N));
}


// ----------------------------------------------------------------------
TF1* recoilAnalysis::getMesFunction(int i) {
  if (i == 0) {
    return fga;
  } else if (i == 1) {
    return fca;
  }
}

// ----------------------------------------------------------------------
// Displays the mes histogram, fits  argus+gaus and prints some numbers
mesData* recoilAnalysis::vubMes(TH1D *h, double &resmean, double &ressigma, double &resalpha, double &resn, int print, int func, double mean, double sigma, double alpha, double n, double argus) {


  cout << " FITTING THE HISTOGRAM " <<  h->GetName() << endl;
  double x = 0.15;
  int lsiz = 2.;
  double tsiz = 0.08;

  double p0, p1, p2, p3;
  double Dp0, Dp1 ,Dp2, Dp3 ;
  double lo, hi;
  double min, max;
  char line[200];

  double cchi2;

  setHist(h, kBlack, 20, 1., 1);
  h->SetNdivisions(-405, "X");
  setTitles(h, "m_{ES} [GeV]", "Entries / 2.5 MeV", 0.08, 2., 1.5, 0.06);
  TH1 *htemp = h->DrawCopy("e");
  
  if (print >= 1) {
    tl.SetTextSize(tsiz); 
    //    tl.DrawTextNDC(0.1, 0.95, h->GetTitle());
  }

  TF1 *farg;

  farg = fa;  
  
  farg->SetParameters(3.*h->GetBinContent(3), -10.);   
  farg->SetParLimits(0, 0., 100000000.);
  farg->SetParLimits(1, -200., 10.);
  
  h->Fit(farg, "RLM", "samee", 5.21, 5.26);

  double par1arg = farg->GetParameter(0); 
  double par2arg = farg->GetParameter(1); 

  fa->SetRange(5.21,5.3);
    
  TF1 *f;

  if (func == 0) {
    f = fga;  
    f->SetParameters(0.9*h->GetMaximum(), 5.28, 0.003, par1arg, par2arg); 
    f->SetParLimits(0, 0., 1000000.);
    f->SetParLimits(1, 5.27, 5.30);
    f->SetParLimits(2, 0.002, 0.005);
    if (mean > 0.) f->FixParameter(1, mean);
    if (sigma > 0.) f->FixParameter(2, sigma);
    if (argus > -100000.) f->FixParameter(4, argus);
  }
  else if (func == 1) {
    f = fca; 
    f->SetParameters(5.28, 0.0025, 1., 2., 0.9*h->GetMaximum(), par1arg, par2arg);
    f->SetParLimits(1, 0.002, 0.005);
    f->SetParLimits(3, -100000., 100000.);
    f->SetParLimits(4, 0., 10000000.);
    f->SetParLimits(5, 0., 10000000.);
    f->SetParLimits(6, -200., 10.);
    if (mean > 0.) f->FixParameter(0, mean);
    if (sigma > 0.) f->FixParameter(1, sigma);
    if (alpha > -100000.) f->FixParameter(2, alpha);
    if (n > -100000.) f->FixParameter(3, n);
    if (argus > -200.) f->FixParameter(6, argus);
  }  
  else if (func == 2) {
    f = fga; 
  }    
  else if (func == 3) {
    f = fca; 
  }    


  f->SetLineStyle(1); f->SetLineWidth(lsiz); f->SetLineColor(kBlack);  f->SetNpx(500);
  if (h->GetSumOfWeights() > 0) {
    h->Fit(f, "RLM", "samee", 5.2, 5.3);
  } else {
    cout << "not enough entries" << endl;
    mesData *pBla = new mesData(h->Integral(h->FindBin(5.27), h->FindBin(5.29)));
    return pBla;
  }
  f->Draw("same");

  //  double nsig(3.);
  //  double sigma = f1->GetParameter(2);
  //  double peak = f1->GetParameter(1);
  lo = 5.27; // peak - nsig*sigma;
  hi = 5.29; // peak + nsig*sigma;

  // -- default number of B's from gauss
  if (func == 0) {
    fg->SetParameters(f->GetParameter(0), f->GetParameter(1), f->GetParameter(2));
    p0 = fg->Integral(lo, hi)/h->GetBinWidth(1);
  }
  else if (func == 1) {
    fc->SetParameters(f->GetParameter(0), f->GetParameter(1), f->GetParameter(2),  f->GetParameter(3),  f->GetParameter(4));
    p0 = fc->Integral(lo, hi)/h->GetBinWidth(1);
  }

  // +1 sigma
  if (func == 0) {
    fg->SetParameters(f->GetParameter(0)+f->GetParError(0),f->GetParameter(1),f->GetParameter(2));
    max = fg->Integral(lo, hi)/h->GetBinWidth(1);
  }
  else if (func == 1) {
    fc->SetParameters(f->GetParameter(0), f->GetParameter(1), f->GetParameter(2), f->GetParameter(3), 
		      f->GetParameter(4)+f->GetParError(4));
    max = fc->Integral(lo, hi)/h->GetBinWidth(1);
  }

  // -1 sigma
  if (func == 0) {
    fg->SetParameters(f->GetParameter(0)-f->GetParError(0),f->GetParameter(1),f->GetParameter(2));
    min = fg->Integral(lo, hi)/h->GetBinWidth(1);
  }
  else if (func == 1) {
    fc->SetParameters(f->GetParameter(0), f->GetParameter(1), f->GetParameter(2), f->GetParameter(3), 
		      f->GetParameter(4)-f->GetParError(4));
    min = fc->Integral(lo, hi)/h->GetBinWidth(1);
  }

 
  Dp0 = (max - min)/2.;
  if(Dp0<1) Dp0 = sqrt(h->Integral(28,40));
  if (print >= 1) {
    sprintf(line, "S = %6.1f +/- %5.1f", p0, Dp0); 
    tl.DrawTextNDC(x, 0.8, line);
  }

  // -- argus background error estimate
  double a1(0.), a1E(0.), a2(0.), a2E(0.); 
  if (func == 0) {
    a1 = f->GetParameter(3); a1E = f->GetParError(3);
    a2 = f->GetParameter(4); a2E = f->GetParError(4);
  }
  else if (func == 1) {
    a1 = f->GetParameter(5); a1E = f->GetParError(5);
    a2 = f->GetParameter(6); a2E = f->GetParError(6);
  }    
  
  fa->SetParameters(a1, a2);
  p1 = fa->Integral(lo, hi)/h->GetBinWidth(1);

  // +1 sigma
  fa->SetParameters(a1+a1E, a2);
  max = fa->Integral(lo, hi)/h->GetBinWidth(1);
  // -1 sigma
  fa->SetParameters(a1-a1E, a2);
  min = fa->Integral(lo, hi)/h->GetBinWidth(1);
  Dp1 = (max - min)/2.;
  if (print >= 1) {
    sprintf(line, "B = %6.1f +/- %5.1f", p1, Dp1); 
    tl.DrawTextNDC(x, 0.7, line);
  }

  fa->SetParameters(a1, a2);
  fa->SetLineStyle(2); fa->SetLineWidth(3); fa->SetLineColor(kRed);
  fa->DrawCopy("same");

  if (func == 0) {
    p2 = f->GetParameter(1);
    p3 = f->GetParameter(2);
    Dp2 = f->GetParError(1);
    Dp3 = f->GetParError(2);
    resmean = f->GetParameter(1);
    ressigma = f->GetParameter(2);
    cchi2 = f->GetChisquare();
  }
  else if (func == 1) {
    p2 = f->GetParameter(0);
    p3 = f->GetParameter(1);
    Dp2 = f->GetParError(0);
    Dp3 = f->GetParError(1);
    resmean = f->GetParameter(0);
    ressigma = f->GetParameter(1);
    resalpha = f->GetParameter(2);
    resn = f->GetParameter(3);
    cchi2 = f->GetChisquare();

  }

  
  if (print >= 1) {
    if (sigma < 0.){sprintf(line, "m = %7.2f +/- %5.2f", 1000.*p3, 1000.*Dp3); }
    else {sprintf(line, "m = %7.2f", 1000.*p3);}
    tl.DrawTextNDC(x, 0.6, line);
    if (mean < 0.){sprintf(line, "s = %7.2f +/- %5.2f", 1000.*p2, 1000.*Dp2); }
    else {sprintf(line, "s = %7.2f ", 1000.*p2);}
    tl.DrawTextNDC(x, 0.5, line);
    
    sprintf(line, "chi2 = %7.2f ", cchi2);
    tl.DrawTextNDC(x, 0.4, line);
  }

 mesData *pD = new mesData(TString(h->GetTitle()),p0,Dp0,p1,Dp1,p0/(p0+p1),dBinomial(p0, p0+p1));

  return pD;
}

// ----------------------------------------------------------------------
// Displays the mes histogram, fits  argus+gaus and prints some numbers
mesData* recoilAnalysis::newMes(TH1D *h, int print, int func) {

  double x = 0.15;
  int lsiz = 5;
  double tsiz = 0.08;

  double p0, p1, p2, p3, p2E, p3E;
  double Dp0, Dp1;
  double lo, hi;
  double min, max;
  char line[200];

  setHist(h, kBlack, 20, 1.2, 1);
  h->SetNdivisions(-405, "X");
  setTitles(h, "m_{ES} [GeV]", "Entries / 2.5 MeV", 0.08, 2., 1.5, 0.06);
  TH1 *htemp = h->DrawCopy("e");
  
  if (print >= 1) {
    tl.SetTextSize(tsiz); 
    tl.DrawTextNDC(0.1, 0.95, h->GetTitle());
  }

  TF1 *f;
  if (func == 0) {
    f = fga; 
    f->SetParameters(0.9*h->GetMaximum(), 5.28, 0.003, 3.*h->GetBinContent(3), -10.);
    f->SetParLimits(0, 0., 1000000.);
    f->SetParLimits(1, 5.27, 5.30);
    f->SetParLimits(2, 0.002, 0.005);
    f->FixParameter(2, 0.00275);  
    f->FixParameter(1, 5.28 );  
  }
  else if (func == 1) {
    f = fca; 
    f->SetParameters(5.28, 0.003, 1., 2., 0.9*h->GetMaximum(), 3.*h->GetBinContent(3), -10.);
    f->SetParLimits(1, 0.002, 0.005);
    f->FixParameter(1, 0.00275);  
    f->FixParameter(0, 5.28 );  
    f->FixParameter(3, 5.);
    f->SetParLimits(4, 0., 10000000.);
    f->SetParLimits(5, 0., 10000000.);
    f->SetParLimits(6, -200., 10.);
  }  
  else if (func == 2) {
    f = fga; 
  }    
  else if (func == 3) {
    f = fca; 
  }    


  f->SetLineStyle(1); f->SetLineWidth(lsiz); f->SetLineColor(kBlack);  f->SetNpx(500); 
  if (h->GetSumOfWeights() > 0.) {
    h->Fit(f, "RL", "samee", 5.205, 5.29);
  } else {
    cout << "not enough entries" << endl;
    mesData *pBla = new mesData(h->Integral(h->FindBin(5.27), h->FindBin(5.29)));
    return pBla;
  }

  //  double nsig(3.);
  //  double sigma = f1->GetParameter(2);
  //  double peak = f1->GetParameter(1);
  lo = 5.27; // peak - nsig*sigma;
  hi = 5.29; // peak + nsig*sigma;

  // -- default number of B's from gauss
  if (func == 0) {
    fg->SetParameters(f->GetParameter(0), f->GetParameter(1), f->GetParameter(2));
    p0 = fg->Integral(lo, hi)/h->GetBinWidth(1);
  }
  else if (func == 1) {
    fc->SetParameters(f->GetParameter(0), f->GetParameter(1), f->GetParameter(2),  f->GetParameter(3),  f->GetParameter(4));
    p0 = fc->Integral(lo, hi)/h->GetBinWidth(1);
  }

  // +1 sigma
  if (func == 0) {
    fg->SetParameters(f->GetParameter(0)+f->GetParError(0),f->GetParameter(1),f->GetParameter(2));
    max = fg->Integral(lo, hi)/h->GetBinWidth(1);
  }
  else if (func == 1) {
    fc->SetParameters(f->GetParameter(0), f->GetParameter(1), f->GetParameter(2), f->GetParameter(3), 
		      f->GetParameter(4)+f->GetParError(4));
    max = fc->Integral(lo, hi)/h->GetBinWidth(1);
  }

  // -1 sigma
  if (func == 0) {
    fg->SetParameters(f->GetParameter(0)-f->GetParError(0),f->GetParameter(1),f->GetParameter(2));
    min = fg->Integral(lo, hi)/h->GetBinWidth(1);
  }
  else if (func == 1) {
    fc->SetParameters(f->GetParameter(0), f->GetParameter(1), f->GetParameter(2), f->GetParameter(3), 
		      f->GetParameter(4)-f->GetParError(4));
    min = fc->Integral(lo, hi)/h->GetBinWidth(1);
  }


  Dp0 = (max - min)/2.;
  if (print >= 1) {
    sprintf(line, "S = %6.1f +/- %5.1f", p0, Dp0); 
    tl.DrawTextNDC(x, 0.8, line);
  }

  // -- argus background error estimate
  double a1(0.), a1E(0.), a2(0.), a2E(0.); 
  if (func == 0) {
    a1 = f->GetParameter(3); a1E = f->GetParError(3);
    a2 = f->GetParameter(4); a2E = f->GetParError(4);
  }
  else if (func == 1) {
    a1 = f->GetParameter(5); a1E = f->GetParError(5);
    a2 = f->GetParameter(6); a2E = f->GetParError(6);
  }    

  fa->SetParameters(a1, a2);
  p1 = fa->Integral(lo, hi)/h->GetBinWidth(1);

  // +1 sigma
  fa->SetParameters(a1+a1E, a2);
  max = fa->Integral(lo, hi)/h->GetBinWidth(1);
  // -1 sigma
  fa->SetParameters(a1-a1E, a2);
  min = fa->Integral(lo, hi)/h->GetBinWidth(1);
  Dp1 = (max - min)/2.;
  if (print >= 1) {
    sprintf(line, "B = %6.1f +/- %5.1f", p1, Dp1); 
    tl.DrawTextNDC(x, 0.7, line);
  }

  fa->SetParameters(a1, a2);
  fa->SetLineStyle(2); fa->SetLineWidth(6); fa->SetLineColor(kRed);
  fa->DrawCopy("same");

  if (func == 0) {
    p2 = f->GetParameter(1);
    p2E= f->GetParError(1);
    p3 = f->GetParameter(2);
    p3E= f->GetParError(2);
  }
  else if (func == 1) {
    p2 = f->GetParameter(0);
    p2E= f->GetParError(0);
    p3 = f->GetParameter(1);
    p3E= f->GetParError(1);
  }

  
  if (print >= 1) {
    sprintf(line, "s = %7.2f", 1000.*p3); 
    tl.DrawTextNDC(x, 0.6, line);
    sprintf(line, "m = %7.2f", 1000.*p2); 
    tl.DrawTextNDC(x, 0.5, line);
  }


  mesData *pD = new mesData(TString(h->GetTitle()),p0,Dp0,p1,Dp1,p0/(p0+p1),dBinomial(p0, p0+p1));
  pD->setMean(p2);
  pD->setMeanE(p2E);

  pD->setSigma(p3);
  pD->setSigmaE(p3E);

  return pD;
}


// ----------------------------------------------------------------------
// Displays the mes histogram, fits  argus+gaus and prints some numbers
mesData* recoilAnalysis::mes(TH1D *h, int print) {

  double x = 0.15;
  int lsiz = 3;
  double tsiz = 0.08;

  double p0, p1, p2, p3;
  double Dp0, Dp1;
  double lo, hi;
  double min, max;
  char line[200];

  setHist(h, kBlack, 20, 1.2, 1);
  h->SetNdivisions(-405, "X");
  setTitles(h, "m_{ES} [GeV]", "Entries / 2.5 MeV", 0.08, 2., 1.5, 0.06);
  TH1 *htemp = h->DrawCopy("e");
  
  if (print >= 1) {
    tl.SetTextSize(tsiz); 
    tl.DrawTextNDC(0.1, 0.95, h->GetTitle());
  }

  f1->SetParameters(0.9*h->GetMaximum(), 5.28, 0.003, 3.*h->GetBinContent(3), -10.);
  f1->SetParLimits(0, 0., 1000000.);
  f1->SetParLimits(1, 5.27, 5.30);
  //f1->FixParameter(2, 0.00275);
  f1->FixParameter(1, 5.279);
  f1->SetParLimits(3, -0.5, 1000000.);
  f1->SetParLimits(2, 0.002, 0.005);
  f1->SetParLimits(4, -100., 100.);

  f1->SetLineStyle(1); f1->SetLineWidth(lsiz); f1->SetLineColor(kBlack);  f1->SetNpx(500);
  if (h->GetSumOfWeights() > fMesMinEvents) {
    h->Fit(f1, "LERM", "samee", 5.2, 5.29);
  } else {
    cout << "not enough entries" << endl;
    mesData *pBla = new mesData(h->Integral(h->FindBin(5.27), h->FindBin(5.29)));
    return pBla;
  }

  //  double nsig(3.);
  //  double sigma = f1->GetParameter(2);
  //  double peak = f1->GetParameter(1);
  lo = 5.27; // peak - nsig*sigma;
  hi = 5.29; // peak + nsig*sigma;

  // -- default number of B's from gauss
  f0->SetParameters(f1->GetParameter(0), f1->GetParameter(1), f1->GetParameter(2));
  p0 = f0->Integral(lo, hi)/h->GetBinWidth(1);
  // +1 sigma
  f0->SetParameters(f1->GetParameter(0)+f1->GetParError(0),f1->GetParameter(1),f1->GetParameter(2));
  max = f0->Integral(lo, hi)/h->GetBinWidth(1);
  // -1 sigma
  f0->SetParameters(f1->GetParameter(0)-f1->GetParError(0),f1->GetParameter(1),f1->GetParameter(2));
  min = f0->Integral(lo, hi)/h->GetBinWidth(1);
  Dp0 = (max - min)/2.;
  if (print >= 1) {
    sprintf(line, "S = %6.1f +/- %5.1f", p0, Dp0); 
    tl.DrawTextNDC(x, 0.8, line);
  }

  // -- argus background error estimate
  f2->SetParameters(0., 0., 0., f1->GetParameter(3), f1->GetParameter(4));
  p1 = f2->Integral(lo, hi)/h->GetBinWidth(1);
  // +1 sigma
  f2->SetParameters(0., 0., 0., f1->GetParameter(3)+f1->GetParError(3), f1->GetParameter(4)+f1->GetParError(4));
  max = f2->Integral(lo, hi)/h->GetBinWidth(1);
  // -1 sigma
  f2->SetParameters(0., 0., 0., f1->GetParameter(3)-f1->GetParError(3), f1->GetParameter(4)-f1->GetParError(4));
  min = f2->Integral(lo, hi)/h->GetBinWidth(1);
  Dp1 = (max - min)/2.;
  if (print >= 1) {
    sprintf(line, "B = %6.1f +/- %5.1f", p1, Dp1); 
    tl.DrawTextNDC(x, 0.7, line);
  }

  f2->SetParameters(0., 0., 0., f1->GetParameter(3), f1->GetParameter(4));
  f2->SetLineStyle(2); f2->SetLineWidth(6); f2->SetLineColor(kRed);
  f2->DrawCopy("same");

  p2 = f1->GetParameter(2);
  p3 = f1->GetParameter(1);
  
  if (print >= 1) {
    sprintf(line, "s = %7.2f", 1000.*p3); 
    tl.DrawTextNDC(x, 0.6, line);
    sprintf(line, "m = %7.2f", 1000.*p2); 
    tl.DrawTextNDC(x, 0.5, line);
  }

  mesData *pD = new mesData(TString(h->GetTitle()),p0,Dp0,p1,Dp1,p0/(p0+p1),dBinomial(p0, p0+p1));
  pD->setMean(f1->GetParameter(1));
  pD->setMeanE(f1->GetParError(1));

  pD->setSigma(f1->GetParameter(2));
  pD->setSigmaE(f1->GetParError(2));

  return pD;
}


// ----------------------------------------------------------------------
// Fit ARGUS+GAUS to mes distribution, determine argus(5.27 .. 5.29)/argus(5.20 .. 5.26)
double recoilAnalysis::getBgScale(TH1D *h) {
  TF1 *f1l = new TF1("f1l", argus, 5.2, 5.29, 5);
  f1l->SetNpx(5000);

  f1l->SetParameters(0.9*h->GetMaximum(), 5.279, 0.003, 3*h->GetBinContent(3), -10.);
  f1l->SetParLimits(0, 0., 1000000.);
  f1l->SetParLimits(1, 5.27, 5.30);
  f1l->SetParLimits(2, 0.002, 0.005);

  h->Draw("e");
  h->Fit("f1l", "lq", "esame", 5.2, 5.29);
  TF1 *f = h->GetFunction("f1l");

  double lo = 5.27;
  double hi = 5.29;

  TF1 *f3l = new TF1("f3l", "gaus", 5.2, 5.29);
  f3l->SetParameters(f->GetParameter(0), f->GetParameter(1), f->GetParameter(2), 0., 0.);
  fNB = f3l->Integral(5.27, 5.29)/h->GetBinWidth(1);

  TF1 *f2l = new TF1("f2l", argus, 5.2, 5.29, 5);
  f2l->SetParameters(0., 0., 0., f->GetParameter(3), f->GetParameter(4));
  f2l->SetLineStyle(2);
  f2l->Draw("same"); 

  double bg = f2l->Integral(5.200, 5.260)/h->GetBinWidth(1); 
  double pe = f2l->Integral(lo, hi)/h->GetBinWidth(1);

  double bgScale;
  if (bg > 0) {
    bgScale = pe/bg;
  }
  else {
    bgScale = 0.;
  }

  delete f1l;
  delete f2l;
  delete f3l;

  cout << "nB = " << fNB << endl;
  cout << "bgScale = " << bgScale << endl;

  return bgScale;
 
}


// ----------------------------------------------------------------------
// provides a sideband subtracted spectrum
TH1D* recoilAnalysis::bgSubtracted(const char *hist, const char *dir) {
  
  char name[100];
  
  sprintf(name, "mes%s", dir);
  TH1D *hmes = (TH1D*)gDirectory->Get(name);
  double bgScale = getBgScale(hmes);
  
  sprintf(name, "bg%s/%s", dir, hist);
  TH1D *hbg =  (TH1D*)gDirectory->Get(name);
  
  sprintf(name, "sg%s/%s", dir, hist);
  TH1D *hsig =  (TH1D*)gDirectory->Get(name);
  
  sprintf(name, "sub_%s", hist);
  TH1D *hsub = new TH1D(name, "", hsig->GetNbinsX(), hsig->GetBinLowEdge(1), hsig->GetBinLowEdge(hsig->GetNbinsX()+1)); 
  hsub->Add(hsig, hbg, 1., -1.*bgScale);
  
  return hsub;
}

// ----------------------------------------------------------------------
// provides a sideband subtracted spectrum
TH1D* recoilAnalysis::bgSubtracted(const char *hist, const char *dir, const char *subdir) {
  
  char name[100];

  sprintf(name, "mes%s", dir);
  TH1D *hmes = (TH1D*)gDirectory->Get(name);
  double bgScale = getBgScale(hmes);
  
  sprintf(name, "bg%s/%s/%s", dir, subdir, hist);
  TH1D *hbg =  (TH1D*)gDirectory->Get(name);
  
  sprintf(name, "sg%s/%s/%s", dir, subdir, hist);
  TH1D *hsig =  (TH1D*)gDirectory->Get(name);
  
  sprintf(name, "sub_%s", hist);
  TH1D *hsub = new TH1D(name, "", hsig->GetNbinsX(), hsig->GetBinLowEdge(1), hsig->GetBinLowEdge(hsig->GetNbinsX()+1)); 
  hsub->Add(hsig, hbg, 1., -1.*bgScale);
  
  return hsub;
}


// ----------------------------------------------------------------------
// provides a sideband subtracted spectrum for directories and an arbitrary mes histogram
TH1D* recoilAnalysis::bgSubtracted(const char *hist, const char *dir, const char *meshist, int func) {
  char name[100];
  

  TH1D *hmes = (TH1D*)gDirectory->Get(meshist);
  double bgScale = getBgScale(hmes);
  cout << "anaVub::bgSubtracted: bgScale = " << bgScale << endl;

  sprintf(name, "bg%s/%s", dir, hist);
  TH1D *hbg =  (TH1D*)gDirectory->Get(name);
  sprintf(name, "sg%s/%s", dir, hist);
  TH1D *hsig =  (TH1D*)gDirectory->Get(name);
  
  sprintf(name, "sub_%s", hist);
  //  TH1D *hsub = new TH1D(name, "", hsig->GetNbinsX(), hsig->GetBinLowEdge(1), hsig->GetBinLowEdge(hsig->GetNbinsX()+1)); 
  TH1D *hsub = new TH1D(*hsig); hsub->SetName(name);
  hsub->Add(hsig, hbg, 1., -1.*bgScale);

  return hsub;
}


// ----------------------------------------------------------------------
// Eliminates crossfeed events in timestamp files
void recoilAnalysis::eraser(char * fileI, int inf, char * fileC, int infC) {

  //Variables and reading of the first file

  double Ipur, pury;
  int plat, part, up, low;
  char arr;
  int lines, linef;
  double Arpury[MAXSIZE];
  Bool_t runeq, plateq, parteq, arreq, upeq, loweq; 
  char buffert[200], filets[100], run[13];
  char buffertC[200], filetsC[100];
  double ArpuryC[MAXSIZE];

  //Reset vectors
  for(int idum = 0; idum<105; idum++){
    fLines[idum] =0;
    fLinef[idum] =0;
  }
  for(int irv = 0; irv<MAXSIZE; irv++){
    fArplatC[irv] = fArpartC[irv] = fArlowC[irv] = fArupC[irv] = 0;
    fArCandC[irv] = ArpuryC[irv] = 0;
    fArarrC[irv] = 0;
    for(int irvr = 0; irvr<12; irvr++) {
      fArrunC[irv][irvr]= 0;
    }
    //    fArrunC[irv][12]= (char)"\0";
    fArrunC[irv][12]= '\0';
    fArplat[irv] = fArpart[irv] = fArlow[irv] = fArup[irv] = 0;
    fArCand[irv] = Arpury[irv] = 0;
    fArarr[irv] = 0;
    for(int irvr = 0; irvr<12; irvr++) {
      fArrun[irv][irvr]= 0;
    }
    //    fArrun[irv][12]= (char)"\0";
    fArrun[irv][12]= '\0';
  }
  
  //Starts reading the "to be compared" file

  sprintf(filetsC,"%s",fileC);
  ifstream tsC(filetsC);
  lines = 0;
  while (tsC.getline(buffertC, 200, '\n')) {
    sscanf(buffertC, "%x:%x:%x/%x:%c %s %lf %lf", &plat, &part, &up, &low,
    	   &arr, run, &pury, &Ipur);
    fArplatC[lines] = plat;
    fArpartC[lines] = part;
    fArupC[lines] = up;
    fArlowC[lines] = low;
    fArarrC[lines] = arr;
    for(int j=0;j<12;j++) {
      fArrunC[lines][j] =run[j];
    }
    //    fArrunC[lines][12] = (char)"\0";
    fArrunC[lines][12] = '\0';
    ArpuryC[lines] = pury;
    fArCandC[lines] = Ipur;
    lines++;
  }

  //Ended the reading of the first file
  //Start the reading of the different file that 
  //have to be comapred
  
  sprintf(filets,"%s",fileI);
  ifstream ts(filets);
  
  linef = 0 ;
  while (ts.getline(buffert, 200, '\n')) {
    sscanf(buffert, "%x:%x:%x/%x:%c %s %lf %lf", &plat, &part, &up, &low,
	   &arr, run, &pury, &Ipur);
    fArplat[linef] = plat;
    fArpart[linef] = part;
    fArup[linef] = up;
    fArlow[linef] = low;
    fArarr[linef] = arr;
    for(int j=0;j<12;j++) {
      fArrun[linef][j] =run[j];
    }
    //    fArrun[linef][12] = (char)"\0";
    fArrun[linef][12] = '\0';
    Arpury[linef] = pury;
    fArCand[linef] = Ipur;
    linef++;
  }

  int numeq =0;
  int numeqC =0;
  for(int ff = 0; ff < linef; ff++) {
    for (int jk = 0; jk < lines; jk++) {
      runeq = plateq = parteq = arreq = upeq = loweq = kFALSE; 
      if(fArlowC[jk] == fArlow[ff] ) {
	loweq = kTRUE;
      }
      if(fArupC[jk] == fArup[ff] ) {
	upeq = kTRUE;
      }
      if(fArarrC[jk] == fArarr[ff] ) {
	arreq = kTRUE;
      }
      if(fArpartC[jk] == fArpart[ff] ) {
	parteq = kTRUE;
      }
      if(fArplatC[jk] == fArplat[ff] ) {
	plateq = kTRUE;
      }
      int count = 0;
      for(int j=0;j<12;j++) {
	if(fArrunC[jk][j] == fArrun[ff][j] ) {
	  count++;
	  if(count == 12) runeq = kTRUE;
	}
      }
      if(runeq && plateq && parteq && arreq && upeq && loweq) {
	fflag[ff][inf] = 1;
	if(Arpury[ff] > ArpuryC[jk]) {
	  numeqC++;
	  fflagC[jk][infC] = 1;
	  fflag[ff][inf] = 0;
	} else {
	  numeq++;
	}
      }
    } 
  }
  fLines[infC] = lines;
  fLinef[inf] = linef;
  cout<<dec<<numeq<< " <-- I file lines"<<endl;
  cout<<dec<<numeqC<<" <-- C file lines"<<endl;
}



// ----------------------------------------------------------------------
// Dump to timestamp output file
void recoilAnalysis::Dump(const char * filename, int fi, int fiC) {

  char tmprun[13], tmprunC[13];
  char name[100];
  char name2[100];    
  double flun;

  sprintf(name2,"%s%s",filename,"_cand");
//    ofstream outFil2(name2,ios::app|ios::uppercase);
//    outFil2.setf(ios::uppercase);
  ofstream outFil2(name2);
  sprintf(name,"%s%s",filename,".root");
//    ofstream outFil(name,ios::app|ios::uppercase);
//    outFil.setf(ios::uppercase);
  ofstream outFil(name);

  int numeqD;
  int tmpLf, tmpLs;

  if(fi != -1) {
    tmpLf = fLinef[fi];
  } else {
    tmpLf = 0;
  }
  if(fiC != -1) {
    tmpLs = fLines[fiC];
  } else {
    tmpLs = 0;
  }

  numeqD = 0;
  if(fi != -1) {

    for(int jtc = 0;jtc<12; jtc++) {
      tmprun[jtc] = 0;
    }
    //    tmprun[12] = (char)"\0";
    tmprun[12] = '\0';
    for(int jnf = 0; jnf < tmpLf; jnf++) {
      for(int jt = 0;jt<12; jt++) {
	tmprun[jt] = fArrun[jnf][jt];
      }
      //      tmprun[12] = (char)"\0";
      tmprun[12] = '\0';
      if(fflag[jnf][fi] == 1) {
	numeqD++;
	outFil<<fArCand[jnf]<<endl;
	outFil2<<hex<<setw(8)<<setfill('0')<<fArup[jnf]<<"/"<<setw(8)<<setfill('0')<<fArlow[jnf]<<endl;
      }
    }
  }
  cout<<numeqD<<" I file lines"<<endl;
  numeqD = 0;
  if(fiC != -1) {

    for(int jtcc = 0;jtcc<12; jtcc++) {
      tmprunC[jtcc] = 0;
    }
    //    tmprunC[12] = (char)"\0";
    tmprunC[12] = '\0';
    for(int jnk = 0; jnk < tmpLs; jnk++) {
      for(int jt = 0;jt<12; jt++) {
	tmprunC[jt] = fArrunC[jnk][jt];
      }
      //      tmprunC[12] = (char)"\0";
      tmprunC[12] = '\0';
      if(fflagC[jnk][fiC] == 1 )  {
	numeqD++;
	outFil<<fArCandC[jnk]<<endl;
	outFil2<<hex<<setw(8)<<setfill('0')<<fArupC[jnk]<<"/"<<setw(8)<<setfill('0')<<fArlowC[jnk]<<endl;
      }
    }
  }
  cout<<numeqD<<" C file lines"<<endl;
}
