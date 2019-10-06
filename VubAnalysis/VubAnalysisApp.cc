#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>
#include <math.h>
#include <stdlib.h>

#include "TROOT.h"
#include "TRint.h"
#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TString.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TMinuit.h"
#include "TString.h"
#include "TRandom.h"
#include "TUnixSystem.h"

#include "VubAnalysis/VubAnalysisCode.hh"


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int main(int argc, char *argv[]) {

  int processID = gSystem->GetPid();
  cout << "Running under process ID  " << processID << endl;

  TString progName  = argv[0]; progName.ReplaceAll("bin/", "");
  TString writeName, fileName, dirName, cutFile("");

  TString treeName("h1");
  Int_t dump(0), writeG(0);
  Int_t file(0), lun(0);
  Int_t isMC(0), isVerbose(0), runCategories(0);
  Int_t dirspec(0);
  Int_t nevents(0), skipevents(0);
  Int_t newFormat(2); 
  Int_t pidKilling(0), pidKillingKaon(0), pidKillingEl(0), pidKillingMu(0);
  Int_t smearTracks(0);
  Int_t smearNeut(0), scaleKlongs(0), filterK0s(1010000);
  Int_t randomSeed(processID);

  // -- command line arguments
  for (int i = 0; i < argc; i++){
    if(strcmp(argv[i],"-C")  == 0) {cutFile    = TString(argv[++i]); }           // file with cuts
    if(strcmp(argv[i],"-c")  == 0) {fileName   = TString(argv[++i]); file = 0; } // file with chain definition
    if(strcmp(argv[i],"-D")  == 0) {dirName = TString(argv[++i]);  dirspec = 1; }// where to put the output
    if(strcmp(argv[i],"-d")  == 0) dump        = atoi(argv[++i]);                // printout
    if(strcmp(argv[i],"-F")  == 0) {filterK0s  = atoi(argv[++i]);}               // require 0 or x K0S->pi0pi0 decays
    if(strcmp(argv[i],"-f")  == 0) {fileName   = TString(argv[++i]); file = 1; } // single file instead of chain
    if(strcmp(argv[i],"-l")  == 0) {runCategories = 1;}                          // run Categories
    if(strcmp(argv[i],"-MC") == 0)  isMC = 1;                                    // flag to avoid warnings
    if(strcmp(argv[i],"-mc") == 0)  isMC = 1;                                    // flag to avoid warnings
    if(strcmp(argv[i],"-n")  == 0) nevents     = atoi(argv[++i]);                // number of events to run 
    if(strcmp(argv[i],"-p")  == 0) pidKilling  = 1;                              // redo PidKilling on tuple
    if(strcmp(argv[i],"-pe")  == 0) pidKillingEl  = 1;                           // redo elPidKilling on tuple
    if(strcmp(argv[i],"-pk")  == 0) pidKillingKaon  = 1;                         // redo muPidKilling on tuple
    if(strcmp(argv[i],"-pm")  == 0) pidKillingMu  = 1;                           // redo kaPidKilling on tuple
    if(strcmp(argv[i],"-r")  == 0) randomSeed  = atoi(argv[++i]);                // set seed for random gen.
										 
    if(strcmp(argv[i],"-t")  == 0) treeName    = TString(argv[++i]);             // ...
    if(strcmp(argv[i],"-s")  == 0) skipevents  = atoi(argv[++i]);                // skip events at beginning
    if(strcmp(argv[i],"-KL") == 0) scaleKlongs = 1;                              // scale KL energies
    if(strcmp(argv[i],"-ST") == 0) smearTracks = 1;                              // smear tracks
    if(strcmp(argv[i],"-SN") == 0) smearNeut = 1;                                // smear neut
    if(strcmp(argv[i],"-v")  == 0)  isVerbose = 1;                               // be verbose
    if(strcmp(argv[i],"-y")  == 0)  newFormat = atoi(argv[++i]);                 // flag to indicate  version
                                                                                 // 0: Feb02, 1: apr 02, 2: aug02
    if(strcmp(argv[i],"-wG")  == 0) {writeName = TString(argv[++i]); writeG = 1;}// write event list into red. file (after all cuts)  

  }

  cout << "Setting gRandom seed: " << randomSeed << endl;
  gRandom->SetSeed(randomSeed);

  // -- Rootfile with output histograms
  TString  barefile(fileName), chainFile, meta, histfile;
  if (file == 0) {
    if (barefile.Contains("chains/")) {
      meta = barefile;
      //      histfile = "outp/" + barefile.ReplaceAll("chains/", "") + ".root";
      histfile = barefile.ReplaceAll("chains/", "") + ".root";
      if(dirspec) histfile = dirName + "/" + histfile;
      //      meta = "chains/" + barefile;
    } else {
      meta = barefile;
      histfile =  barefile + ".root";
      if(dirspec) histfile = dirName + "/" + histfile;
    }

    // -- The following lines strip everything from the string up to and including the last '/'
    int fl = barefile.Last('/');
    TString bla(barefile);
    bla.Replace(0, fl+1, ' '); bla.Strip(TString::kLeading, ' ');  bla.Remove(0,1);
    histfile =  bla + ".root";
    if(dirspec) histfile = dirName + "/" + histfile;
  }  else if (file == 1) {
    // -- The following lines strip everything from the string up to and including the last '/'
    int fl = barefile.Last('/');
    TString bla(barefile);
    bla.Replace(0, fl+1, ' '); bla.Strip(TString::kLeading, ' ');  bla.Remove(0,1);
    histfile =  bla;
    if(dirspec) histfile = dirName + "/" + histfile;
  }
  // use the following line to dump the output into some other place: "scratch/anaRecoil/bla.root"
  TString histFile = histfile;



  cout << "Opening " << histfile.Data() << " for output histograms" << endl;
  cout << "Opening " << fileName.Data() << " for input" << endl;

  // -- get lumi-weight of this chain
  char buffer[200]; 
  float lw8(1.); 
  TString lfile = TString("chains/lumi.") + barefile; 
  ifstream lf(lfile.Data()); 
  if (!lf) {
    cout << "Cannot open " << lfile.Data()    << " for lumi weight, assuming weight 1" << endl;
  } else {
    while (lf.getline(buffer, 200, '\n')) {
      if (buffer[0] == '#') continue; 
      sscanf(buffer, "%f", &lw8); 
    }
    lf.close();
    cout << "Opened " << lfile.Data()    << " for lumi weight " << lw8 << endl;
  }


  // -- Set up chain
  TChain *chain = new TChain(treeName);
  cout << "Chaining ... " << treeName << endl;
  char pName[2000]; 
  int nentries; 
  if (file == 0) {
    ifstream is(meta);  
    while(meta.ReadLine(is) && (!meta.IsNull())){ 
      nentries = -1;
      if (meta.Data()[0] == '#') continue; 
      sscanf(meta.Data(), "%s %d", pName, &nentries); 
      if (nentries > -1) {
	cout << pName << " -> " << nentries << " entries" << endl; 
	chain->Add(pName, nentries); 
      } else {
	cout << meta << endl; 
	chain->Add(meta); 
      }
    }
    is.close();
  }
  else if (file == 1) {
    cout << fileName << endl;
    chain->Add(fileName);
  }

  VubAnalysisCode a(chain,isMC, newFormat);
  if (pidKilling) {
    a.redoPidKilling(1);
    a.redoPidKillingEl(1);
    a.redoPidKillingMu(1);
    a.redoPidKillingKaon(1);
  }
  if (pidKillingKaon) a.redoPidKillingKaon(1);
  if (pidKillingEl)   a.redoPidKillingEl(1);
  if (pidKillingMu)   a.redoPidKillingMu(1);
  
  if (strcmp(cutFile.Data(), "")) {
    a.readCuts(cutFile, 0);
  } else {
    a.readCuts(TString("b0cuts.dat"), 0);
  }
  a.dumpCuts();
  a.openHistFile(histfile);
  
  a.doKlongScaling(scaleKlongs);
  a.doTrackSmearing(smearTracks);
  a.doNeutSmearing(smearNeut);
  a.runCategories(runCategories); 
  a.filterK0s(filterK0s); 
  a.setLumiWeight(lw8); 
  if (writeG) a.makeEventList(2); 

  a.bookHist(dump);
  cout << "====================" << endl;
  cout << "isVerbose      = " << isVerbose << endl;
  cout << "PidKilling     = " << pidKilling << endl;
  cout << "PidKillingKaon = " << pidKillingKaon << endl;
  cout << "TrkSmearing    = " << smearTracks << endl;
  cout << "NeuSmearing    = " << smearNeut << endl;
  cout << "====================" << endl;


  a.switchOffReading("platform");
  a.switchOffReading("partition");
  //    a.switchOffReading("upperID");
  //    a.switchOffReading("lowerID");

  //      a.switchOffReading("beamS*");
  a.switchOffReading("beamSCov*");

  a.switchOffReading("nTrkTot");
  a.switchOffReading("W2");
  a.switchOffReading("FoxWol2");
  a.switchOffReading("FoxWol2Neu");
  a.switchOffReading("thrust");
  a.switchOffReading("thrustNeu");

  a.switchOffReading("errMassB0");
  a.switchOffReading("m0B0");
  a.switchOffReading("xB0");
  a.switchOffReading("yB0");
  a.switchOffReading("zB0");
  a.switchOffReading("s2xB0");
  a.switchOffReading("s2yB0");
  a.switchOffReading("s2zB0");
  a.switchOffReading("chi2B0");
  a.switchOffReading("dofB0");
  a.switchOffReading("stB0");
  a.switchOffReading("ndauB0");
  a.switchOffReading("ThruB0");
  a.switchOffReading("thThruB0");
  a.switchOffReading("phiThruB0");
  //      a.switchOffReading("d1B0Index");
  //      a.switchOffReading("d1B0Lund");
  //      a.switchOffReading("d2B0Index");
  //      a.switchOffReading("d2B0Lund");
  //      a.switchOffReading("d3B0Index");
  //      a.switchOffReading("d3B0Lund");
  //      a.switchOffReading("d4B0Index");
  //      a.switchOffReading("d4B0Lund");
  //      a.switchOffReading("d5B0Index");
  //      a.switchOffReading("d5B0Lund");
  //      a.switchOffReading("d6B0Index");
  //      a.switchOffReading("d6B0Lund");
  //      a.switchOffReading("d7B0Index");
  //      a.switchOffReading("d7B0Lund");
  a.switchOffReading("VtxXLepB0");
  a.switchOffReading("VtxYLepB0");
  a.switchOffReading("VtxZLepB0");
  a.switchOffReading("VtxCovXXLepB0");
  a.switchOffReading("VtxCovYYLepB0");
  a.switchOffReading("VtxCovXYLepB0");
  a.switchOffReading("VtxCovZZLepB0");
  a.switchOffReading("VtxCovXZLepB0");
  a.switchOffReading("VtxCovYZLepB0");
  a.switchOffReading("VtxChiSqLepB0");
  a.switchOffReading("VtxNDofLepB0");
  a.switchOffReading("VtxStatLepB0");
  a.switchOffReading("VtxNUsedLepB0");
  a.switchOffReading("DocaLepB0");
  a.switchOffReading("DocaErrLepB0");
  a.switchOffReading("VtxXXB0");
  a.switchOffReading("VtxYXB0");
  a.switchOffReading("VtxZXB0");
  a.switchOffReading("VtxCovXXXB0");
  a.switchOffReading("VtxCovYYXB0");
  a.switchOffReading("VtxCovXYXB0");
  a.switchOffReading("VtxCovZZXB0");
  a.switchOffReading("VtxCovXZXB0");
  a.switchOffReading("VtxCovYZXB0");
  a.switchOffReading("VtxChiSqXB0");
  a.switchOffReading("VtxNDofXB0");
  a.switchOffReading("VtxStatXB0");
  a.switchOffReading("VtxNUsedXB0");
  a.switchOffReading("VtxPXB0");
  a.switchOffReading("VtxPhiXB0");
  a.switchOffReading("VtxThetaXB0");
  a.switchOffReading("ThrustXB0");
  a.switchOffReading("ThrustXPhiB0");
  a.switchOffReading("ThrustXThetaB0");
  //      a.switchOffReading("MassPB0");
  //      a.switchOffReading("MassPhiB0");
  //      a.switchOffReading("MassThetaB0");
  a.switchOffReading("Cov00B0");
  a.switchOffReading("Cov10B0");
  a.switchOffReading("Cov11B0");
  a.switchOffReading("Cov20B0");
  a.switchOffReading("Cov21B0");
  a.switchOffReading("Cov22B0");
  a.switchOffReading("Cov30B0");
  a.switchOffReading("Cov31B0");
  a.switchOffReading("Cov32B0");
  a.switchOffReading("Cov33B0");


  a.switchOffReading("errMassChB");
  a.switchOffReading("m0ChB");
  a.switchOffReading("xChB");
  a.switchOffReading("yChB");
  a.switchOffReading("zChB");
  a.switchOffReading("s2xChB");
  a.switchOffReading("s2yChB");
  a.switchOffReading("s2zChB");
  a.switchOffReading("chi2ChB");
  a.switchOffReading("dofChB");
  a.switchOffReading("stChB");
  a.switchOffReading("ndauChB");
  a.switchOffReading("ThruChB");
  a.switchOffReading("thThruChB");
  a.switchOffReading("phiThruChB");
  //      a.switchOffReading("d1ChBIndex");
  //      a.switchOffReading("d1ChBLund");
  //      a.switchOffReading("d2ChBIndex");
  //      a.switchOffReading("d2ChBLund");
  //      a.switchOffReading("d3ChBIndex");
  //      a.switchOffReading("d3ChBLund");
  //      a.switchOffReading("d4ChBIndex");
  //      a.switchOffReading("d4ChBLund");
  //      a.switchOffReading("d5ChBIndex");
  //      a.switchOffReading("d5ChBLund");
  //      a.switchOffReading("d6ChBIndex");
  //      a.switchOffReading("d6ChBLund");
  //      a.switchOffReading("d7ChBIndex");
  //      a.switchOffReading("d7ChBLund");
  a.switchOffReading("VtxXLepChB");
  a.switchOffReading("VtxYLepChB");
  a.switchOffReading("VtxZLepChB");
  a.switchOffReading("VtxCovXXLepChB");
  a.switchOffReading("VtxCovYYLepChB");
  a.switchOffReading("VtxCovXYLepChB");
  a.switchOffReading("VtxCovZZLepChB");
  a.switchOffReading("VtxCovXZLepChB");
  a.switchOffReading("VtxCovYZLepChB");
  a.switchOffReading("VtxChiSqLepChB");
  a.switchOffReading("VtxNDofLepChB");
  a.switchOffReading("VtxStatLepChB");
  a.switchOffReading("VtxNUsedLepChB");
  a.switchOffReading("DocaLepChB");
  a.switchOffReading("DocaErrLepChB");
  a.switchOffReading("VtxXXChB");
  a.switchOffReading("VtxYXChB");
  a.switchOffReading("VtxZXChB");
  a.switchOffReading("VtxCovXXXChB");
  a.switchOffReading("VtxCovYYXChB");
  a.switchOffReading("VtxCovXYXChB");
  a.switchOffReading("VtxCovZZXChB");
  a.switchOffReading("VtxCovXZXChB");
  a.switchOffReading("VtxCovYZXChB");
  a.switchOffReading("VtxChiSqXChB");
  a.switchOffReading("VtxNDofXChB");
  a.switchOffReading("VtxStatXChB");
  a.switchOffReading("VtxNUsedXChB");
  a.switchOffReading("VtxPXChB");
  a.switchOffReading("VtxPhiXChB");
  a.switchOffReading("VtxThetaXChB");
  a.switchOffReading("ThrustXChB");
  a.switchOffReading("ThrustXPhiChB");
  a.switchOffReading("ThrustXThetaChB");
  //      a.switchOffReading("MassPChB");
  //      a.switchOffReading("MassPhiChB");
  //      a.switchOffReading("MassThetaChB");
  a.switchOffReading("Cov00ChB");
  a.switchOffReading("Cov10ChB");
  a.switchOffReading("Cov11ChB");
  a.switchOffReading("Cov20ChB");
  a.switchOffReading("Cov21ChB");
  a.switchOffReading("Cov22ChB");
  a.switchOffReading("Cov30ChB");
  a.switchOffReading("Cov31ChB");
  a.switchOffReading("Cov32ChB");
  a.switchOffReading("Cov33ChB");

  a.switchOffReading("d1DstarIndex");
  a.switchOffReading("d1DstarLund");
  a.switchOffReading("d2DstarIndex");
  a.switchOffReading("d2DstarLund");
  a.switchOffReading("nDstarBS");
  a.switchOffReading("massDstarBS");
  a.switchOffReading("chi2DstarBS");
  a.switchOffReading("dofDstarBS");
  

  a.switchOffReading("nD0");
  a.switchOffReading("massD0");
  a.switchOffReading("pD0");
  a.switchOffReading("thD0");
  a.switchOffReading("phiD0");
  a.switchOffReading("errMassD0");
  a.switchOffReading("m0D0");
  a.switchOffReading("xD0");
  a.switchOffReading("yD0");
  a.switchOffReading("zD0");
  a.switchOffReading("s2xD0");
  a.switchOffReading("s2yD0");
  a.switchOffReading("s2zD0");
  a.switchOffReading("chi2D0");
  a.switchOffReading("dofD0");
  a.switchOffReading("stD0");
  a.switchOffReading("ndauD0");
  if(isMC)   a.switchOffReading("MCD0");
  a.switchOffReading("d1D0Index");
  a.switchOffReading("d1D0Lund");
  a.switchOffReading("d2D0Index");
  a.switchOffReading("d2D0Lund");
  a.switchOffReading("d3D0Index");
  a.switchOffReading("d3D0Lund");
  a.switchOffReading("d4D0Index");
  a.switchOffReading("d4D0Lund");
  a.switchOffReading("nChD");
  a.switchOffReading("massChD");
  a.switchOffReading("pChD");
  a.switchOffReading("thChD");
  a.switchOffReading("phiChD");
  a.switchOffReading("errMassChD");
  a.switchOffReading("m0ChD");
  a.switchOffReading("xChD");
  a.switchOffReading("yChD");
  a.switchOffReading("zChD");
  a.switchOffReading("s2xChD");
  a.switchOffReading("s2yChD");
  a.switchOffReading("s2zChD");
  a.switchOffReading("chi2ChD");
  a.switchOffReading("dofChD");
  a.switchOffReading("stChD");
  a.switchOffReading("ndauChD");
  if(isMC)   a.switchOffReading("MCChD");
  a.switchOffReading("d1ChDIndex");
  a.switchOffReading("d1ChDLund");
  a.switchOffReading("d2ChDIndex");
  a.switchOffReading("d2ChDLund");
  a.switchOffReading("d3ChDIndex");
  a.switchOffReading("d3ChDLund");
  a.switchOffReading("d4ChDIndex");
  a.switchOffReading("d4ChDLund");
  
  a.switchOffReading("nPi0");
  a.switchOffReading("massPi0");
  a.switchOffReading("pPi0");
  a.switchOffReading("thPi0");
  a.switchOffReading("phiPi0");
  a.switchOffReading("errMassPi0");
  a.switchOffReading("m0Pi0");
  a.switchOffReading("xPi0");
  a.switchOffReading("yPi0");
  a.switchOffReading("zPi0");
  a.switchOffReading("s2xPi0");
  a.switchOffReading("s2yPi0");
  a.switchOffReading("s2zPi0");
  a.switchOffReading("chi2Pi0");
  a.switchOffReading("dofPi0");
  a.switchOffReading("stPi0");
  a.switchOffReading("ndauPi0");
  if(isMC)   a.switchOffReading("MCPi0");
  //    a.switchOffReading("d1Pi0Index");
  a.switchOffReading("d1Pi0Lund");
  //    a.switchOffReading("d2Pi0Index");
  a.switchOffReading("d2Pi0Lund");
  
  
  a.switchOffReading("nDalitz");
  a.switchOffReading("massDalitz");
  a.switchOffReading("pDalitz");
  a.switchOffReading("thDalitz");
  a.switchOffReading("phiDalitz");
  a.switchOffReading("errMassDalitz");
  a.switchOffReading("m0Dalitz");
  a.switchOffReading("xDalitz");
  a.switchOffReading("yDalitz");
  a.switchOffReading("zDalitz");
  a.switchOffReading("s2xDalitz");
  a.switchOffReading("s2yDalitz");
  a.switchOffReading("s2zDalitz");
  a.switchOffReading("chi2Dalitz");
  a.switchOffReading("dofDalitz");
  a.switchOffReading("stDalitz");
  a.switchOffReading("ndauDalitz");
  if(isMC)   a.switchOffReading("MCDalitz");
  a.switchOffReading("d1DalitzIndex");
  a.switchOffReading("d1DalitzLund");
  a.switchOffReading("d2DalitzIndex");
  a.switchOffReading("d2DalitzLund");
  a.switchOffReading("nJpsi");
  a.switchOffReading("massJpsi");
  a.switchOffReading("pJpsi");
  a.switchOffReading("thJpsi");
  a.switchOffReading("phiJpsi");
  a.switchOffReading("errMassJpsi");
  a.switchOffReading("m0Jpsi");
  a.switchOffReading("xJpsi");
  a.switchOffReading("yJpsi");
  a.switchOffReading("zJpsi");
  a.switchOffReading("s2xJpsi");
  a.switchOffReading("s2yJpsi");
  a.switchOffReading("s2zJpsi");
  a.switchOffReading("chi2Jpsi");
  a.switchOffReading("dofJpsi");
  a.switchOffReading("stJpsi");
  a.switchOffReading("ndauJpsi");
  if(isMC)   a.switchOffReading("MCJpsi");
  a.switchOffReading("d1JpsiIndex");
  a.switchOffReading("d1JpsiLund");
  a.switchOffReading("d2JpsiIndex");
  a.switchOffReading("d2JpsiLund");

  a.switchOffReading("Ifr*");

  a.switchOffReading("lMomTrk");
  a.switchOffReading("ZMom42Trk");
  //  a.switchOffReading("ecalTrk");
  a.switchOffReading("ecalXTrk");
  a.switchOffReading("ecalYTrk");
  a.switchOffReading("ecalZTrk");
  a.switchOffReading("nCryTrk");
  a.switchOffReading("nBumpTrk");
  a.switchOffReading("ZMom20Trk");
  a.switchOffReading("secMomTrk");
  a.switchOffReading("s1s9Trk");
  a.switchOffReading("s9s25Trk");
  a.switchOffReading("erawTrk");
  a.switchOffReading("phiClusterTrk");
  a.switchOffReading("thetaClusterTrk");
  a.switchOffReading("phicMatTrk");
  a.switchOffReading("trkcMatTrk");
  a.switchOffReading("nPidTrk");
  a.switchOffReading("emcStatusTrk");
  //  a.switchOffReading("phiAtEMCTrk");
  a.switchOffReading("isvtTrk");
  //  a.switchOffReading("nsvtTrk");
  a.switchOffReading("fhitTrk");
  a.switchOffReading("lhitTrk");
  //  a.switchOffReading("tLenTrk");
  //  a.switchOffReading("ntdofTrk");
  //  a.switchOffReading("tproTrk");
  //  a.switchOffReading("tChi2Trk");
  a.switchOffReading("cPidTrk");
  a.switchOffReading("sfRangeTrk");
  a.switchOffReading("trkFitStatusTrk");
  a.switchOffReading("covEETrk");
  a.switchOffReading("covTTTrk");
  a.switchOffReading("covPPTrk");
  a.switchOffReading("covRRTrk");


  //  a.switchOffReading("lMomGam");
  //  a.switchOffReading("ZMom42Gam");
  a.switchOffReading("ecalXTrk");
  a.switchOffReading("ecalYTrk");
  a.switchOffReading("ecalZTrk");
  //  a.switchOffReading("nCryGam");
  a.switchOffReading("nBumpGam");
  //  a.switchOffReading("ZMom20Gam");
  //  a.switchOffReading("secMomGam");
  //  a.switchOffReading("s1s9Gam");
  //  a.switchOffReading("s9s25Gam");
  a.switchOffReading("erawGam");
  a.switchOffReading("phiClusterGam");
  a.switchOffReading("thetaClusterGam");
  a.switchOffReading("emcStatusGam");
  a.switchOffReading("covEEGam");
  //    a.switchOffReading("covTTGam");
  //    a.switchOffReading("covPPGam");
  a.switchOffReading("covRRGam");

  
  a.Loop(nevents, skipevents, isVerbose, lun);

  if (writeG) {
    a.dumpEventList(writeName.Data());
  }

  a.closeHistFile();

  return 0;
    
}
