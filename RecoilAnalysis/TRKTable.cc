#include "TRKTable.hh"

#include <iostream>
#include <math.h>
#include <stdlib.h>

using namespace std;

ClassImp(TRKTable)

// ----------------------------------------------------------------------
TRKTable::TRKTable() {
  DR = 57.2956;
  fVerbose = 1;
  fDataVector = new TList();
  fRandom = new TF1("f1", "gaus", -10., 10.);
  fRandom->SetParameters(1., 0., 1.);
  fFlat = new TF1("fTRKTable", "pol0", 0., 1.);
  fFlat->SetParameter(0, 1.);
}

// ----------------------------------------------------------------------
TRKTable::TRKTable(const char *s, const char *f) {
  DR = 57.2956;
  fVerbose = 1;
  fDataVector = new TList();
  fRandom = new TF1("f1", "gaus", -10., 10.);
  fRandom->SetParameters(1., 0., 1.);
  fFlat = new TF1(f, "pol0", 0., 1.);
  fFlat->SetParameter(0, 1.);
  readFromFile(s);
}

// ----------------------------------------------------------------------
TRKTable::~TRKTable() {
  delete fRandom;
}


// ----------------------------------------------------------------------
void TRKTable::readFromFile(const char *filename) {
  
  TString fullName;
  if (filename[0] == '/') {
    fullName = filename;
  }
  else {
    fullName = filename;
  }

  fFileName = fullName;

  float pmin, pmax, mmin, mmax, tmin, tmax, fmin, fmax, e, se;
  if (fVerbose > 0) cout << "Read table contents from " << fullName 
			 << " and function name " << fFlat->GetName()
			 << endl;

  char  buffer[200];
  ifstream is(fullName);
  int readLines(0);
  while (is.getline(buffer, 200, '\n')) {
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%f %f %f %f %f %f %f %f %f %f", 
	   &pmin, &pmax, &mmin, &mmax, &tmin, &tmax, &fmin, &fmax, &e, &se);
    //    TRKData datum(pmin, pmax, mmin, mmax, tmin, tmax, fmin, fmax, e, se);
    //    fDataVector.push_back(datum);
    TRKData *datum = new TRKData(pmin, pmax, mmin, mmax, tmin, tmax, fmin, fmax, e, se);
    fDataVector->Add(datum);    
    readLines++;
  }
  if (readLines < 1) {
    if (fVerbose > -1) cout << "  Nothing read from " << fullName << endl; 
  }
  is.close();
}


// ----------------------------------------------------------------------
void TRKTable::printAll() {
  //  for (vector<TRKData>::iterator iTable = fDataVector.begin(); iTable != fDataVector.end(); iTable++) {
  TIter next(fDataVector); TRKData *iTable;
  while ((iTable = (TRKData*)next())) {
    cout << *iTable << endl;
  }
}


// ----------------------------------------------------------------------
void TRKTable::print() {

  int cells(0);
  double a(0.), wi(0.), w(0.), meanE(0.), meanS(0.), meanU(0.);
  double tmin(9999.), tmax(-9999.), pmin(9999.), pmax(-9999.);
  
  //    for (vector<TRKData>::iterator iTable = fDataVector.begin(); iTable != fDataVector.end(); iTable++) {
  TIter next(fDataVector); TRKData *iTable;
  while ((iTable = (TRKData*)next())) {
    cells++;
    if (iTable->getPmin() < pmin) pmin = iTable->getPmin();
    if (iTable->getPmax() > pmax) pmax = iTable->getPmax();
    if (iTable->getTmin() < tmin) tmin = iTable->getTmin();
    if (iTable->getTmax() > tmax) tmax = iTable->getTmax();
    
    a = iTable->getS();
    if (a > 1e-8) {
      wi = 1./(a*a);
      w     += wi;
      meanE += wi*iTable->getE();
      meanS += wi*iTable->getS();
    }
  }

  meanE /= w;
  meanU  = 1./TMath::Sqrt(w);  // uncertainty of mean efficiency
  meanS /= w;                  // mean of efficiency uncertainties

  char line[200];
  sprintf(line, "============");  for (int i = 0; i < getFileName().Length(); ++i) sprintf(line, "%s=", line);
  cout << line << endl;
  cout << "Summary for " << getFileName() << endl;
  sprintf(line, "------------");  for (int i = 0; i < getFileName().Length(); ++i) sprintf(line, "%s-", line);
  cout << line << endl;

  cout << " Number of cells:    " << cells << endl;
  cout << " -------------------------------" << endl;

  cout << " min momentum:       " << pmin << endl;
  cout << " max momentum:       " << pmax << endl;
  cout << " min polar angle:    " << tmin << endl;
  cout << " max polar angle:    " << tmax << endl;
  cout << " -------------------------------" << endl;

  cout << " <Efficiency>:       " << meanE << endl;
  cout << "                 +/- " << meanU << endl;
  cout << " <Delta Efficiency>: " << meanS << endl;
  cout << " -------------------------------" << endl;

}

// ----------------------------------------------------------------------
void TRKTable::smear(int mode) { 
  double bla(0.);

  //    for (vector<TRKData>::iterator iTable = fDataVector.begin(); iTable != fDataVector.end(); ++iTable) {
  TIter next(fDataVector); TRKData *iTable;
  while ((iTable = (TRKData*)next())) {
    if (mode == 0) {
      // do nothing
    }
    else if (mode == 99) {
      bla = iTable->getE() + fRandom->GetRandom()*iTable->getS();
      if (bla < 0.) {
	bla = 0.;
	cout << "Reset to 0 in cell (" << iTable->getPmin() << "," << iTable->getTmin() << "," << iTable->getFmin()<<")"<<endl;
      } else if (bla > 1.) {
	bla = 1.;
	cout << "Reset to 1 in cell (" << iTable->getPmin() << "," << iTable->getTmin() << "," << iTable->getFmin()<<")"<<endl;
      }
      iTable->setE(bla);
    } else {
      bla = iTable->getE() + mode*iTable->getS();
      if (bla < 0.) bla = 0.;
      iTable->setE(bla);
    }
  }
}
  

// ----------------------------------------------------------------------
void TRKTable::shiftRel(double shift) { 
  double bla(0.);
  if (TMath::Abs(shift) > 0.0001) {
    //      for (vector<TRKData>::iterator iTable = fDataVector.begin(); iTable != fDataVector.end(); ++iTable) {
    TIter next(fDataVector); TRKData *iTable;
    while ((iTable = (TRKData*)next())) {
      bla = iTable->getE() + shift*iTable->getE();
      if (bla < 0.) {
	bla = 0.;
	cout << "Reset to 0 in cell (" << iTable->getPmin() << "," << iTable->getTmin() << "," << iTable->getFmin()<<")"<<endl;
      } else if (bla > 1.) {
	bla = 1.;
	cout << "Reset to 1 in cell (" << iTable->getPmin() << "," << iTable->getTmin() << "," << iTable->getFmin()<<")"<<endl;
      }
      iTable->setE(bla);
    }
  }  else {
    cout << fFileName << " no shifting, because shift = " << shift << endl;
  }
}
  

// ----------------------------------------------------------------------
double TRKTable::w8D(double p, double m, double t, double f) { 
  return w8R(p, m, t*DR, f*DR);
}


// ----------------------------------------------------------------------
double TRKTable::w8R(double p, double m, double t, double f) {
  //    for (vector<TRKData>::iterator iTable = fDataVector.begin(); iTable != fDataVector.end(); ++iTable) {
  TIter next(fDataVector); TRKData *iTable;
  while ((iTable = (TRKData*)next())) {
    // cout << "TRKTable::w8R()..." << endl;
    if (iTable->isCell(p, m, t, f)) {
      return iTable->getE();
    }
  }
  return -99.;
}


double TRKTable::w8Derr(double p, double m, double t, double f) { 
  return w8Rerr(p, m, t*DR, f*DR);
}

// ----------------------------------------------------------------------
double TRKTable::w8Rerr(double p, double m, double t, double f) {
  // for (vector<TRKData>::iterator iTable = fDataVector.begin(); iTable != fDataVector.end(); ++iTable) {
  TIter next(fDataVector); TRKData *iTable;
  while ((iTable = (TRKData*)next())) {
    // cout << "TRKTable::w8R()..." << endl;
    if (iTable->isCell(p, m, t, f)) {
      return iTable->getS();
    }
  }
  return -99.;
}


// ----------------------------------------------------------------------
void TRKTable::printHeader(ofstream &EFF) {
  EFF << "# Efficiency table derived from" << endl; 
  EFF << "# Filename: " << getFileName() << endl;
  EFF << "# " << endl;
  EFF << "# Momentum bins are 125 MeV wide and go from 0. to 5 GeV." << endl;
  EFF << "# " << endl;
  EFF << "# ----------------------------------------------------------------------" << endl;
}


// ----------------------------------------------------------------------
TRKData* TRKTable::getData(double p, double m, double t, double f) {
  //    for (vector<TRKData>::iterator iTable = fDataVector.begin(); iTable != fDataVector.end(); ++iTable) {
  TIter next(fDataVector); TRKData *iTable;
  while ((iTable = (TRKData*)next())) {
    if (iTable->isCell(p, m, t, f)) return iTable;
  }
  if (fVerbose > -1) cout << "did not find TRKData for given values of p, t, f." << endl;
  return 0;
}


// ----------------------------------------------------------------------
void TRKTable::clear() {
  //    if (fVerbose > 0) cout << "Data size before clearance: " << fDataVector.size() << endl;
  //  fDataVector.erase(fDataVector.begin(), fDataVector.end());
  fDataVector->Delete();
}


// ----------------------------------------------------------------------
void TRKTable::flush() {
  //    for (vector<TRKData>::iterator iTable = begin(); iTable != end(); ++iTable) {
  TIter next(fDataVector); TRKData *iTable;
  while ((iTable = (TRKData*)next())) {
    iTable->setE(0.);
    iTable->setS(0.);
  }
}
