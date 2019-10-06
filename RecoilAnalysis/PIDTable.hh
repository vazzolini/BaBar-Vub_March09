#ifndef PIDTABLE
#define PIDTABLE

#include "TString.h"
#include "TObject.h"
#include "TFile.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"

#include "TSortedList.h"
#include "TObjArray.h"
#include "TCollection.h"
#include "TOrdCollection.h"
#include "TList.h"

#include "PIDData.hh"

class PIDTable: public TObject {
public:

  PIDTable();
  PIDTable(const char *, const char *);
  ~PIDTable();

  // used in many places, including those outside of PIDTable
  static double getDR() {return 57.2956;} 

  // -- Input 
  void readFromFile(const char *filename = "bla.dat");
  // VALERY: instantiate from 2 3D histograms, need for hadronic fractions dump
  void readFromHist(TH3D *passed, TH3D *total);

  // -- Output
  // VALERY: restored from the previous version, need for hadronic fractions dump
  void dumpToFile(const char *filename = "bla.dat");  // vanilla dump, nothing is done to the table


  // -- Utilities
  void     setVerbosity(int t) {fVerbose = t;}
  // VALERY: extract PIDTable structure into an empty 3D histogram
  TH3D     *extractPIDTableStructure(const char *inhistname);
  void     printAll();
  void     print();
  void     clear();    // clears the data Vector
  void     flush();    // fills all cells with 0-0-0-0

  //    PIDData* begin() {return fDataVector.begin();}
  //    PIDData* end() {return fDataVector.end();}
  //    int      size() {return fDataVector.size();}

  TString  getFileName() {return fFileName;}
  PIDData* getDataR(double p, double t, double f);
  PIDData* getData(double p, double t, double f);

  void smear(int mode = 0);                         // replace cell efficiency by smeared (gaussian with sigma = err) value
  void shiftRel(double shift = 0.1);                // replace cell efficiency by shifted (oldvalue +/- shift*oldvalue) value
  void shiftSigma(double shift = 0.1);              // replace cell efficiency by shifted (oldvalue +/- shift*sigma) value

  double effD(double momentum, double theta, double phi);  
  double effR(double momentum, double theta, double phi);
  Bool_t idD(double momentum, double theta, double phi);  
  Bool_t idR(double momentum, double theta, double phi);  
  // VALERY: restored from the older version
  double errD(double momentum, double theta, double phi);  
  double errR(double momentum, double theta, double phi);

private:

  void printHeader(ofstream &);

  //  vector<PIDData> fDataVector;
  TList *fDataVector;

  TString fFileName;
  int fVerbose;
  TF1 *fRandom, *fFlat;
  
  ClassDef(PIDTable,1) //Testing PIDTable

};


#endif
