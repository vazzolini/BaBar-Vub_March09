#ifndef TRKTABLE
#define TRKTABLE

#include "TString.h"
#include "TObject.h"
#include "TF1.h"

#include <fstream>
#include "TSortedList.h"
#include "TObjArray.h"
#include "TCollection.h"
#include "TOrdCollection.h"
#include "TList.h"

#include "TRKData.hh"

class TRKTable: public TObject {

public:

  TRKTable();
  TRKTable(const char *, const char *);
  ~TRKTable();

  // -- Input 
  void readFromFile(const char *filename = "bla.dat");

  // -- Utilities
  void     setVerbosity(int t) {fVerbose = t;}
  void     printAll();
  void     print();
  void     clear();    // clears the data Vector
  void     flush();    // fills all cells with 0-0-0-0

//    TRKData* begin() {return fDataVector.begin();}
//    TRKData* end() {return fDataVector.end();}
//    int      size() {return fDataVector.size();}

  TString  getFileName() {return fFileName;}
  TRKData* getData(double p, double m, double t, double f);

  void smear(int mode = 0);   // replace eff. by smeared (gaussian with sigma = err) value
  void shiftRel(double shift = 0.1);  // replace eff. by shifted (oldvalue +/- shift*oldvalue) value

  // -- eff and err with angles in degrees and radian  
  double w8D(double momentum, double mult, double theta, double phi);  
  double w8R(double momentum, double mult, double theta, double phi);
  double w8Derr(double momentum, double mult, double theta, double phi);  
  double w8Rerr(double momentum, double mult, double theta, double phi);


private:

  void printHeader(ofstream &);

  //  vector<TRKData> fDataVector;
  TList *fDataVector;

  TString fFileName;
  int fVerbose;
  double DR; 

  TF1 *fRandom, *fFlat;
  
  ClassDef(TRKTable,1) //Testing TRKTable

};


#endif
