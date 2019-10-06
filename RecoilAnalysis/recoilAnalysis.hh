#ifndef recoilAnalysis_h
#define recoilAnalysis_h

#include <iostream>

#include <TROOT.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TChain.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <TObjString.h>
#include <TTree.h>
#include <TLatex.h>
#include "mesData.hh"


#define MAXSIZE 12

class recoilAnalysis {

private: 
  double fNB, fNBe; 

  TF1 *fc, *fg, *fa, *fga, *fca; 
  TF1 *f0, *f1, *f2; 

  TLatex tl;

  int fMesMinEvents; 

  //  double argus(double *x, double *par); 

public :
  recoilAnalysis();
  virtual ~recoilAnalysis();

  Int_t nMc, nDa, nTs, nTsc; 
  TFile *fMc[20], *fDa[20];          // array for files 
  TObjString *fTs[105], *fTsc[105];  
  double fMcLumi[20], fDaLumi[20];   // array with corresponding luminosities

  // -- Display
  void setTitles(TH1 *h, const char *sx, const char *sy,
		 float size=0.05, float xoff=1.1, float yoff=1.0, float lsize=0.05, int font=62);
  void setHist(TH1 *h, int color=1, int symbol=20, float size=1., int width=2);
  void setFilledHist(TH1 *h, int color=1, int fillcolor=1, int fillstyle=3004, int width=1);

  void loadFiles(const char *name);
  void loadTimestamp(const char *filename);
  void loadTimestampI(const char *filename);
  void loadTimestampC(const char *filename);
  void loadMc(const char *name, double lumi=-1.);
  void loadDa(const char *name, double lumi=-1.);
  void loadTs(const char *name);
  void loadTsC(const char *name);
  void lsMc();
  void lsTs();
  void lsTsc();
  void lsDa();
  void eraser(char * fileI = 0, int inf = 0, char * fileC = 0, int infC = 0);
  void Dump(const char * filename, int fi = -1, int fiC = -1);
  double dBinomial(double, double); 
  // -- fit mes distribution in a histogram
  mesData* mes(TH1D *h, int print=1);
  mesData* newMes(TH1D *h, int print=1, int func=0);
  mesData* vubMes(TH1D *h, double &resmean, double &ressigma, double &resalpha, double &resn, int print=1, int func=0, double mean=0, double sigma=0, double alpha=0, double n=0, double argus=0);
  void setMesMinimum(int i) {fMesMinEvents = i;} 
  int  getMesMinimum() {return fMesMinEvents;}
  TF1*  getMesFunction(int i);

  // -- Determine scale factor for sideband histograms. Stores number of B candidates in fNB. 
  double getBgScale(TH1D *h);

  // -- Make sideband subtracted histogram
  TH1D* bgSubtracted(const char *hist, const char *dir);
  TH1D* bgSubtracted(const char *hist, const char *dir, const char *subdir );
  TH1D* bgSubtracted(const char *hist, const char *dir, const char *postfix, int func);

  double getNB() {return fNB; }

  //Timestamp stuff

  int fArplatC[MAXSIZE], fArlowC[MAXSIZE], fArupC[MAXSIZE], fArpartC[MAXSIZE];
  char fArarrC[MAXSIZE], fArrunC[MAXSIZE][13];
  double fArCand[MAXSIZE], fArCandC[MAXSIZE];
  int fArplat[MAXSIZE], fArlow[MAXSIZE], fArup[MAXSIZE], fArpart[MAXSIZE];
  char fArarr[MAXSIZE], fArrun[MAXSIZE][13];
  //How many files to compare? 
  int fflagC[MAXSIZE][105];
  int fflag[MAXSIZE][105];
  int fLines[105], fLinef[105];

  // ----------------------------------------------------------------------
  ClassDef(recoilAnalysis,1) 

};

#endif
