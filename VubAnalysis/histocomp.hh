#ifndef HISTOCOMP_h
#define HISTOCOMP_h

#include <iostream.h>

#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
#include "TArrow.h"
#include "TBox.h"


#include "TString.h"
#include "TObject.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TVirtualPad.h"  // access to gPad

// ----------------------------------------------------------------------
class number {
public:
  //    number(double x = 0.) {val = x; stat = err1 = err2 = 0.;}
  //    number(double x = 0., double e1) {val = x; stat = e1; err1 = err2 = 0.;}
  //    number(double x = 0., double e1, double e2) {val = x; stat = e1; err1 = e2; err2 = 0.;}
  number(double x = 0., double e1 = 0., double e2 = 0., double e3 = 0.) {val = x; stat = e1; err1 = e2; err2 = e3;}
  number& operator=(const number &a) {val = a.val; stat = a.stat; err1 = a.err1; err2 = a.err2; return *this;} 
  number& operator=(double a) {val = a; stat = err1 = err2 = 0.; return *this; }
  operator double() {return val;}
  void  set(double x = 0., double e1 = 0., double e2 = 0., double e3 = 0.) {val = x; stat = e1; err1 = e2; err2 = e3;}

  double val; 
  double stat, err1,  err2;
};

// ----------------------------------------------------------------------
class fitResult {
public: 
  fitResult() {r = b2u = b2c = oth = fB2u = fB2c = fOth = eff = effmx = efftot = nsl = chi2 = ssb = 0.;}
  double calc() {if (eff > 0.) return b2u/eff/effmx/(nsl*fac*fpstar); else return 0.;}
  void print() {
    cout << Form("nsl:*    %5.4f", nsl) <<  endl;
    cout << Form("b2u:*    %5.4f +/- %5.4f +/- %5.4f", b2u.val, b2u.stat, b2u.err1) <<  endl;
    cout << Form("b2c:     %5.4f +/- %5.4f +/- %5.4f", b2c.val, b2c.stat, b2c.err1) <<  endl;
    cout << Form("oth:     %5.4f +/- %5.4f +/- %5.4f", oth.val, oth.stat, oth.err1) <<  endl;
    cout << Form("efftot:  %5.4f +/- %5.4f", efftot.val, efftot.stat) <<  endl;
    cout << Form("effmx:*  %5.4f +/- %5.4f", effmx.val, effmx.stat) <<  endl;
    cout << Form("eff:*    %5.4f +/- %5.4f", eff.val, eff.stat) <<  endl;
    cout << Form("fpstar:* %5.4f", fpstar) <<  endl;
    cout << Form("fac:*    %5.4f", fac) <<  endl;
    cout << Form("chi2:    %5.4f", chi2) <<  endl;
    cout << Form("ssb:     %5.4f", ssb) <<  endl;
    cout << Form("BRBR:    %5.4f +/- %5.4f +/- %5.4f", r.val, r.stat, r.err1) <<  endl;
    //    cout << Form("->   r = %5.4f", calc()) << endl;
  }

  number r;                  // r = Gamma(b2ulnu)/Gamma(b2xlnu)
  number b2u, b2c, oth;      // first bin number of events
  number fB2u, fB2c, fOth;   // fractions 
  number eff, effmx, efftot; // efficiency
  double fpstar, fac; 
  double nsl; 
  double chi2, ssb;          // chi2 of fit and s/sqrt(s+b)

};

// ----------------------------------------------------------------------
class histocomp {

public :
  histocomp();
  ~histocomp();

  // -- Main analysis methods 
  void makeAll();
  void allYields();
  void allEffTables();
  void allKshorts(const char *what = "all"); 

  void generatorMxPlots(const char *file="genmx");
  void kaonSpectra(const char *hist = "kp", const char *seed = "all", int offset = 5);
  void playKs(const char *cuts="nr==0", const char *epsname="playKs.eps", const char *files="test");
  void yields(const char *seed = "dstar", const char *hist = "mesallevents"); 
  void mesSys(const char *seed = "mes", const char *hist = "mesall");
  // void x1x2(const char *seed = "dataMc", int mode = 0, const char *what = "all", const char *whatMesHist="mesa1"); // mode = 0 Data/MC, 2 Run1/Run2, 3 cock/gen
  void x1x2(const char *file1, const char *file2, const char *dir);
  void envelop(const char *fil="envelop", const char *mc1="smeared", const char *mc2="not smeared");
  void effTables(const char *label = "dstar");
  void bias(const char *seed = "sp3sp4"); 

  void plotCorr(const char *variable="mm2", const char *sel="vub&&pcms>1&&GoodEvent", 
		double min = -10., double max = 10.,
		const char *filename = "all");
  void plotMxCorrelations(const char *sel="vub&&pcms>1&&GoodEvent", 
			  const char *epsname="mx-vcbpcms.eps", 
			  const char *filename = "all");
  void plotMxEfficiency(const char *cuts = "vub&&pcms>1", const char *filename = "brecovub"); // change to merged breco-vub
  void sfDependency(double deltamb=0.150, int allcuts=1, const char *filename = "sx-sp3run1-brecovubnres.root");
  void resolutions(const char *seed="default");
  void mbreco(const char *file="sx-cocktail.root");
  void displayFits(const char *pattern, double min = 0., double max = 2.0, 
		   double ydef = 0.0170, double ydefE = 0.0043, 
		   double xdef = 1.55, 
		   const char *xlabel = "M_{X} [GeV]");
  void dumpTable(int choice = 1); 
  void oliver(double pcut = 0.9);
  void loopOliver(); 
  
  // -- collection of miscellaneous plots
  void smearValidation();
  void pidTableSigma();

  // -- Helper methds
  void dumpTableRow(const char *table, const char *filename, const char *prefix);
  fitResult* readFitOutput(const char *file);
  void setNorm(int normArea = 1) {fNormArea = normArea;} 
  void effTable(const char *label, const char *dir="bla"); 
  TH1D* bgSub(const char *hist, const char *dir, const char *f = "gauss", int controlPlots = 0);
  TH1D* bgSubtracted(const char *hist, const char *dir, const char *f = "gauss", const char *postfix = "A", int controlPlots = 0);
  void theOverlay(TFile *f1, TFile *f2, const char *hist = "a1000", const char *dir = "vub", const char *meshist="mesvub", int nPad = 1);
  void triOverlay(TFile *f1, TFile *f2, TFile *f3, const char *hist = "a1000", const char *dir = "vub", const char *meshist="mesvub", int nPad = 1);
  void makeCanvas(int i = 3);
  void fillStrings(const char *s1 = "", const char *s2 = "", const char *s3 = "", const char *s4 = "", const char *s5 = "");

  Double_t fermi(double kp, double m, double a);
  Double_t fw8(double kp, double mb, double a);

  // -- Files
  void loadFiles(const char *name);
  void loadMc(const char *name, double lumi = 0.);
  void loadDa(const char *name, double lumi = 0.);
  void unloadFiles();
  void unloadMc();
  void unloadDa();
  void lsMc();
  void lsDa();



  // -- public data members
  Int_t nMc, nDa, nTs, nTsc; 
  double fMcLumi[50], fDaLumi[50];   // array with corresponding luminosities
  TFile *fMc[50], *fDa[50];

  // -- Display utilities
  TPad *fPads[50];
  TCanvas *c0, *c1, *c2, *c3, *c4, *c5, *c6, *c7;
  TLatex *tl; 
  TBox *box;
  TArrow *pa; 
  TLine *pl; 
  TLegend *leg;
  TLegendEntry *legge;

  char line[200];


private: 
  TString fString1, fString2, fString3, fString4, fString5;
  TString fTexFile, fFil;

  int fX1X2;
  char fHistName[100][100];
  char fMesHist[100][100];

  TF1 *fc, *fg, *fa, *fga, *fca; 
  TF1 *f0, *f1, *f2, *f3; 
  TF1 *fgg, *fggg, *fp2g;

  Double_t BMASS, BQMASS, A0;

  double fNormalization1, fNormalization2, fNormalization3;
  int fNormArea;

  // ----------------------------------------------------------------------

};


// -- Shape function
inline Double_t histocomp::fermi(double kp, double m, double a) {
  double x = kp/(BMASS - m);
  if ((kp>-m) && (x <= 1.)) {
    return TMath::Power((1-x), a) * TMath::Exp((1+a)*x); 
  } 
  return 0.;
}
inline Double_t histocomp::fw8(double kp, double mb, double a) {
  return fermi(kp, mb, a) / fermi(kp, BQMASS, A0);
}


//  double fermimotion(double kplus,double mB,double mb,double a)
//  {
//    double result = 0;
//    double x = kplus/(mB-mb);
//    if ( kplus >= -mb && x <= 1 ) 
//      result = pow(1-x,a)*exp((1+a)*x);
//    return result;
//  }

//  double fermimotionACCMM(double pF,double mB,double mSp,double spF)
//  {
//    double result = 0;
//    if ( pF > 0 && mB*mB + mSp*mSp - 2*mB*sqrt(pF*pF+mSp*mSp) > 0)
//      result = 4*pF*pF/sqrt(M_PI)/pow(spF,3)*exp(-pF*pF/spF/spF);
//    return result;
//  }



#endif
