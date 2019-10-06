#include "PIDTable.hh"

#include "TF1.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>


ClassImp(PIDTable)

// ----------------------------------------------------------------------
PIDTable::PIDTable() {
  fVerbose = 1;
  fDataVector = new TList();
  fRandom = new TF1("f1", "gaus", -10., 10.);
  fRandom->SetParameters(1., 0., 1.);
  fFlat = new TF1("fPIDTable", "pol0", 0., 1.);
  fFlat->SetParameter(0, 1.);
}


// ----------------------------------------------------------------------
PIDTable::PIDTable(const char *s, const char *f) {
  fVerbose = 1;
  fDataVector = new TList();
  fRandom = new TF1("f1", "gaus", -10., 10.);
  fRandom->SetParameters(1., 0., 1.);
  fFlat = new TF1(f, "pol0", 0., 1.);
  fFlat->SetParameter(0, 1.);

  readFromFile(s);
}


PIDTable::~PIDTable() {
  delete fRandom;
}



// ----------------------------------------------------------------------
void PIDTable::readFromFile(const char *filename) {
  
  TString fullName;
  if (filename[0] == '/') {
    fullName = filename;
  }
  else {
    fullName = filename;
  }

  fFileName = fullName;

  float pmin, pmax, tmin, tmax, fmin, fmax, e, se;
  float pass, total;
  if (fVerbose > 0) std::cout << "Read table contents from " << fullName 
			 << " and function name " << fFlat->GetName()
			 << std::endl;

  char  buffer[200];
  ifstream is(fullName);
  int readLines(0);
  int ntottot=0;
  int npasstot=0;
  while (is.getline(buffer, 200, '\n')) {
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%f %f %f %f %f %f %f %f %f %f", 
	   &pmin, &pmax, &tmin, &tmax, &fmin, &fmax, &e, &se, &pass, &total);
    ntottot+=total;
    npasstot+=pass;
    if(fVerbose > 1) std::cout << "\tread from cell: " << pmin << ", " << pmax << ", " << tmin << ", " << tmax << ", " << fmin << ", " << fmax << ", " << e << ", " << se << ", " << pass << ", " << total << std::endl;
    if(fVerbose > 3) std::cout << "\t\t1";
    PIDData *datum = new PIDData(pmin, pmax, tmin, tmax, fmin, fmax, e, se, pass, total);
    if(fVerbose > 3) std::cout << "\t2 (" << datum << ")";
    fDataVector->Add(datum);
    if(fVerbose > 3) std::cout << "\t3" << std::endl;
    readLines++;
  }
  if (readLines < 1) {
    if (fVerbose > -1) std::cout << "  Nothing read from " << fullName << std::endl; 
  }
  if (fVerbose > -1) std::cout << "  Read from table: passed=" << npasstot << " out of total=" << ntottot << std::endl;
  is.close();
}


// ----------------------------------------------------------------------
void PIDTable::readFromHist(TH3D *hPass, TH3D *hTot) {
  char passName[255], totName[255];
  strcpy(passName, hPass->GetName());
  strcpy(totName, hTot->GetName());

  if (fVerbose > 0) std::cout << "Reading from " << passName  << " and " << totName << std::endl; 
  double n, N;
  
  for (Int_t ix = 1; ix <= hPass->GetNbinsX(); ++ix) {
    for (Int_t iy = 1; iy <= hPass->GetNbinsY(); ++iy) {
      for (Int_t iz = 1; iz <= hPass->GetNbinsZ(); ++iz) {
	n = hPass->GetBinContent(hPass->GetBin(ix, iy, iz));
	N = hTot->GetBinContent(hPass->GetBin(ix, iy, iz)); 

	double thmin, thmax, pmin, pmax, phimin, phimax;
	pmin=hPass->GetXaxis()->GetBinLowEdge(ix);
	pmax=hPass->GetXaxis()->GetBinLowEdge(ix+1);
	thmin=hPass->GetYaxis()->GetBinLowEdge(iy);
	thmax=hPass->GetYaxis()->GetBinLowEdge(iy+1);
	phimin=hPass->GetZaxis()->GetBinLowEdge(iz);
	phimax=hPass->GetZaxis()->GetBinLowEdge(iz+1);
	
	// std::cout << "Reading bin from hist: " << pmin << ", " << pmax << ", " << thmin << ", " << thmax << ", " << phimin << ", " << phimax << ", " << std::endl;

	PIDData *datum= new PIDData(pmin, pmax, thmin, thmax, phimin, phimax, 
				    -99., -99., n, N);
	fDataVector->Add(datum);
      }
    }
  }
}



// ----------------------------------------------------------------------
void PIDTable::dumpToFile(const char *psname) {
  if (fVerbose > 0) std::cout << "Dumping table to file " << psname << std::endl;
  ofstream EFF(psname);
  printHeader(EFF);
  TIter next(fDataVector); PIDData *iTable;
  while ((iTable = (PIDData*)next())) {
    EFF << *iTable << std::endl;
    
    /*
    double thmin, thmax, pmin, pmax, phimin, phimax;
    pmin = iTable->getPmin();
    pmax = iTable->getPmax();
    thmin = iTable->getTmin();
    thmax = iTable->getTmax();
    phimin = iTable->getFmin();
    phimax = iTable->getFmax();
    
    std::cout << "Dumping bin: " << pmin << ", " << pmax << ", " << thmin << ", " << thmax << ", " << phimin << ", " << phimax << ", " << std::endl;
    */
  }
  EFF.close();
}


// ----------------------------------------------------------------------
TH3D *PIDTable::extractPIDTableStructure(const char *inhistname) {
  // iterate over the thing and extract the structure in the form of an empry 3-D histogram
  float radconv=3.14159/180.;
  float xbins[100], ybins[100], zbins[100];
  float ybins_rad[100], zbins_rad[100];
  int nbinsx=0, nbinsy=0, nbinsz=0;
  bool xstop=false, ystop=false, zstop=false;
  double thmin, thmax, pmin, pmax, phimin, phimax;
  double prev_pmin=-99, prev_pmax=-99, prev_thmin=-99, prev_thmax=-99, prev_phimin=-999, prev_phimax=-999;

  TIter next(fDataVector); PIDData *iTable;
  while ((iTable = (PIDData*)next())) {
    pmin = iTable->getPmin();
    pmax = iTable->getPmax();
    thmin = iTable->getTmin();
    thmax = iTable->getTmax();
    phimin = iTable->getFmin();
    phimax = iTable->getFmax();
    
    std::cout << "Read: " << pmin << ", " << pmax << ", " << thmin << ", " << thmax << ", " << phimin << ", " << phimax << ", " << std::endl;
	  
    // check for new bin!!!
    if((pmin>=prev_pmax)) {
      prev_pmin=pmin; prev_pmax=pmax;
      std::cout << "bin X: " << nbinsx << " - (" << pmin << ", " << pmax << ")" << std::endl; 
      xbins[nbinsx]=pmin;
      xbins[nbinsx+1]=pmax;
      nbinsx++;
    }
    if((thmin>=prev_thmax)) {
      prev_thmin=thmin; prev_thmax=thmax;
      std::cout << "bin Y: " << nbinsy << " - (" << thmin << ", " << thmax << ")" << std::endl; 
      ybins[nbinsy]=thmin;
      ybins[nbinsy+1]=thmax;
      ybins_rad[nbinsy]=thmin*radconv;
      ybins_rad[nbinsy+1]=thmax*radconv;
      nbinsy++;
    }
    if((phimin>=prev_phimax)) {
      prev_phimin=phimin; prev_phimax=phimax;
      std::cout << "bin Z: " << nbinsz << " - (" << phimin << ", " << phimax << ")" << std::endl; 
      zbins[nbinsz]=phimin;
      zbins[nbinsz+1]=phimax;
      zbins_rad[nbinsz]=phimin*radconv;
      zbins_rad[nbinsz+1]=phimax*radconv;
      nbinsz++;
    }
  }
  std::cout << "\nloaded table structure with (nbinsx=" << nbinsx << ", nbinsy=" << nbinsy << ", nbinsz=" << nbinsz << "):\n\n";
  for(int i=0; i<nbinsx; i++) {
    for(int j=0; j<nbinsy; j++) {
      for(int k=0; k<nbinsz; k++) {
	// std::cout << "\t++ " << xbins[i] << "\t" << xbins[i+1] << "\t" << ybins[j] << "\t" << ybins[j+1] << "\t" << zbins[k] << "\t" << zbins[k+1] << std::endl;
	// std::cout << "\t-- " << xbins[i] << "\t" << xbins[i+1] << "\t" << ybins_rad[j] << "\t" << ybins_rad[j+1] << "\t" << zbins_rad[k] << "\t" << zbins_rad[k+1] << std::endl;
      }
    }
  }
  
  char histname[255];
  TH3D *hist=new TH3D(inhistname, inhistname, nbinsx, xbins, nbinsy, ybins, nbinsz, zbins);
  hist->Clear();

  return hist;
}


// ----------------------------------------------------------------------
void PIDTable::printAll() {
  //    for (vector<PIDData>::iterator iTable = fDataVector.begin(); iTable != fDataVector.end(); iTable++) {
  TIter next(fDataVector); PIDData *iTable;
  while ((iTable = (PIDData*)next())) {
    std::cout << *iTable << std::endl;
  }

}


// ----------------------------------------------------------------------
void PIDTable::print() {

  int cells(0);
  double pass(0.), total(0.);
  double a(0.), wi(0.), w(0.), meanE(0.), meanS(0.), meanU(0.);
  double tmin(9999.), tmax(-9999.), pmin(9999.), pmax(-9999.);
  
  //    for (vector<PIDData>::iterator iTable = fDataVector.begin(); iTable != fDataVector.end(); iTable++) {
  TIter next(fDataVector); PIDData *iTable;
  while ((iTable = (PIDData*)next())) {
    pass  += iTable->getPass();
    total += iTable->getTot();
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
  std::cout << line << std::endl;
  std::cout << "Summary for " << getFileName() << std::endl;
  sprintf(line, "------------");  for (int i = 0; i < getFileName().Length(); ++i) sprintf(line, "%s-", line);
  std::cout << line << std::endl;

  std::cout << " Number of cells:    " << cells << std::endl;
  std::cout << " Number of pass :    " << pass  << std::endl;
  std::cout << " Number of total:    " << total << std::endl;
  std::cout << " -------------------------------" << std::endl;

  std::cout << " min momentum:       " << pmin << std::endl;
  std::cout << " max momentum:       " << pmax << std::endl;
  std::cout << " min polar angle:    " << tmin << std::endl;
  std::cout << " max polar angle:    " << tmax << std::endl;
  std::cout << " -------------------------------" << std::endl;

  std::cout << " <Efficiency>:       " << meanE << std::endl;
  std::cout << "                 +/- " << meanU << std::endl;
  std::cout << " <Delta Efficiency>: " << meanS << std::endl;
  std::cout << " -------------------------------" << std::endl;

}


// ----------------------------------------------------------------------
void PIDTable::smear(int mode) { 
  double bla(0.);
  //  for (vector<PIDData>::iterator iTable = fDataVector.begin(); iTable != fDataVector.end(); ++iTable) {
  TIter next(fDataVector); PIDData *iTable;
  while ((iTable = (PIDData*)next())) {
    if (mode == 0) {
      // do nothing
    }
    else if (mode == 99) {
      bla = iTable->getE() + fRandom->GetRandom()*iTable->getS();
      if (bla < 0.) {
	bla = 0.;
	std::cout << "Reset to 0 in cell (" << iTable->getPmin() << "," << iTable->getTmin() << "," << iTable->getFmin()<<")"<<std::endl;
      } else if (bla > 1.) {
	bla = 1.;
	std::cout << "Reset to 1 in cell (" << iTable->getPmin() << "," << iTable->getTmin() << "," << iTable->getFmin()<<")"<<std::endl;
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
void PIDTable::shiftRel(double shift) { 
  double bla(0.);
  if (TMath::Abs(shift) > 0.0001) {
    //      for (vector<PIDData>::iterator iTable = fDataVector.begin(); iTable != fDataVector.end(); ++iTable) {
    TIter next(fDataVector); PIDData *iTable;
    while ((iTable = (PIDData*)next())) {
      bla = iTable->getE() + shift*iTable->getE();
      if (bla < 0.) {
	bla = 0.;
	std::cout << "Reset to 0 in cell (" << iTable->getPmin() << "," << iTable->getTmin() << "," << iTable->getFmin()<<")"<<std::endl;
      } else if (bla > 1.) {
	bla = 1.;
	std::cout << "Reset to 1 in cell (" << iTable->getPmin() << "," << iTable->getTmin() << "," << iTable->getFmin()<<")"<<std::endl;
      }
      iTable->setE(bla);
    }
  }  else {
    std::cout << fFileName << " no shifting, because shift = " << shift << std::endl;
  }
}


// ----------------------------------------------------------------------
void PIDTable::shiftSigma(double shift) { 
  double bla(0.);
  if (TMath::Abs(shift) > 0.0001) {
    TIter next(fDataVector); PIDData *iTable;
    while ((iTable = (PIDData*)next())) {
      bla = iTable->getE() + shift*iTable->getS();
      if (bla < 0.) {
	bla = 0.;
	std::cout << "Reset to 0 in cell (" << iTable->getPmin() << "," << iTable->getTmin() << "," << iTable->getFmin()<<")"<<std::endl;
      } else if (bla > 1.) {
	bla = 1.;
	std::cout << "Reset to 1 in cell (" << iTable->getPmin() << "," << iTable->getTmin() << "," << iTable->getFmin()<<")"<<std::endl;
      }
      iTable->setE(bla);
    }
  }  else {
    std::cout << fFileName << " no shifting, because shift = " << shift << std::endl;
  }
}


// ----------------------------------------------------------------------
Bool_t PIDTable::idD(double p, double t, double f) { 
  //    for (vector<PIDData>::iterator iTable = fDataVector.begin(); iTable != fDataVector.end(); ++iTable) {
  TIter next(fDataVector); PIDData *iTable;
  while ((iTable = (PIDData*)next())) {
    if (iTable->isCell(p, t, f)) {
      double bla = fFlat->GetRandom();
      // std::cout << fFlat->GetName() << "  rv = " << bla << " table value = " << iTable->getE() << " p = " << p << std::endl;
      return (bla < iTable->getE()); 
    }
  }
  return kFALSE;
}
// ----------------------------------------------------------------------
Bool_t PIDTable::idR(double p, double t, double f) { 
  return idD(p, t*getDR(), f*getDR());
}
  

// ----------------------------------------------------------------------
double PIDTable::effD(double p, double t, double f) { 
  //    for (vector<PIDData>::iterator iTable = fDataVector.begin(); iTable != fDataVector.end(); ++iTable) {
  TIter next(fDataVector); PIDData *iTable;
  while ((iTable = (PIDData*)next())) {
    if (iTable->isCell(p, t, f)) {
      return iTable->getE();
    }
  }
  return -99.;
}
// ----------------------------------------------------------------------
double PIDTable::effR(double momentum, double theta, double phi) {
  return effD(momentum, theta*getDR(), phi*getDR());
}


// ----------------------------------------------------------------------
double PIDTable::errD(double p, double t, double f) { 
  //    for (vector<PIDData>::iterator iTable = fDataVector.begin(); iTable != fDataVector.end(); ++iTable) {
  TIter next(fDataVector); PIDData *iTable;
  while ((iTable = (PIDData*)next())) {
    if (iTable->isCell(p, t, f)) {
      return iTable->getS();
    }
  }
  return -99.;
}
// ----------------------------------------------------------------------
double PIDTable::errR(double momentum, double theta, double phi) {
  return errD(momentum, theta*getDR(), phi*getDR());
}


// ----------------------------------------------------------------------
void PIDTable::printHeader(ofstream &EFF) {
  EFF << "# Efficiency table derived from" << std::endl; 
  EFF << "# Filename: " << getFileName() << std::endl;
  EFF << "# " << std::endl;
  EFF << "# Momentum bins are 125 MeV wide and go from 0. to 5 GeV." << std::endl;
  EFF << "# " << std::endl;
  EFF << "# ----------------------------------------------------------------------" << std::endl;
}


    
// ----------------------------------------------------------------------
PIDData* PIDTable::getDataR(double momentum, double theta, double phi) {
  return getData(momentum, theta*getDR(), phi*getDR());
}

PIDData* PIDTable::getData(double p, double t, double f) {
  //    for (vector<PIDData>::iterator iTable = fDataVector.begin(); iTable != fDataVector.end(); ++iTable) {
  TIter next(fDataVector); PIDData *iTable;
  while ((iTable = (PIDData*)next())) {
    if (iTable->isCell(p, t, f)) return iTable;
  }
  if (fVerbose > -1) std::cout << "did not find PIDData for given values of p, t, f." << std::endl;
  return 0;
}


// ----------------------------------------------------------------------
void PIDTable::clear() {
  //    if (fVerbose > 0) std::cout << "Data size before clearance: " << fDataVector.size() << std::endl;
  //    fDataVector.erase(fDataVector.begin(), fDataVector.end());
  fDataVector->Delete();
}


// ----------------------------------------------------------------------
void PIDTable::flush() {
  //    for (vector<PIDData>::iterator iTable = begin(); iTable != end(); ++iTable) {
  TIter next(fDataVector); PIDData *iTable;
  while ((iTable = (PIDData*)next())) {
    iTable->setPass(0.);
    iTable->setTot(0.);
    iTable->setE(0.);
    iTable->setS(0.);
  }
}
