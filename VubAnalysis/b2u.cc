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

#include "VubAnalysis/b2uNtp.hh"


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int main(int argc, char *argv[]) {

  int processID = gSystem->GetPid();
  cout << "==> b2u> Running under process ID  " << processID << endl;

  TString progName  = argv[0]; progName.ReplaceAll("bin/", "");
  TString writeName, fileName, dirName, cutFile("a00");

  TString treeName("h1");
  int dump(4);
  int file(0);
  int isMC(0), isVerbose(0);
  int dirspec(0);
  int nevents(0), skipevents(0);
  int randomSeed(processID);
  int pidKillingEl(0), pidKillingMu(0), pidKillingKa(0), smearNeut(0); 

  // -- command line arguments
  for (int i = 0; i < argc; i++){
    if(strcmp(argv[i],"-C")  == 0) {cutFile    = TString(argv[++i]); }           // file with cuts
    if(strcmp(argv[i],"-c")  == 0) {fileName   = TString(argv[++i]); file = 0; } // file with chain definition
    if(strcmp(argv[i],"-D")  == 0) {dirName = TString(argv[++i]);  dirspec = 1; }// where to put the output
    if(strcmp(argv[i],"-d")  == 0) dump        = atoi(argv[++i]);                // level of root output
    if(strcmp(argv[i],"-f")  == 0) {fileName   = TString(argv[++i]); file = 1; } // single file instead of chain
    if(strcmp(argv[i],"-MC") == 0)  isMC = 1;                                    // flag to avoid warnings
    if(strcmp(argv[i],"-mc") == 0)  isMC = 1;                                    // flag to avoid warnings
    if(strcmp(argv[i],"-n")  == 0) nevents     = atoi(argv[++i]);                // number of events to run 
    if(strcmp(argv[i],"-r")  == 0) randomSeed  = atoi(argv[++i]);                // set seed for random gen.
    if(strcmp(argv[i],"-SN") == 0) smearNeut = 1;                                // smear neut

    if(strcmp(argv[i],"-pe") == 0) pidKillingEl  = 1;                           // redo elPidKilling on tuple
    if(strcmp(argv[i],"-pm") == 0) pidKillingMu  = 1;                           // redo muPidKilling on tuple
    if(strcmp(argv[i],"-pk") == 0) pidKillingKa  = 1;                           // redo kaPidKilling on tuple
    if(strcmp(argv[i],"-p")  == 0) {                                            // redo all PidKilling on tuple
      pidKillingKa  = 1;
      pidKillingEl  = 1;
      pidKillingMu  = 1;
    }
    
    if(strcmp(argv[i],"-t")  == 0) treeName    = TString(argv[++i]);             // ...
    if(strcmp(argv[i],"-s")  == 0) skipevents  = atoi(argv[++i]);                // skip events at beginning
    if(strcmp(argv[i],"-v")  == 0)  isVerbose = 1;                               // be verbose
  }

  cout << "==> b2u> Setting gRandom seed: " << randomSeed << endl;
  gRandom->SetSeed(randomSeed);

  // -- Rootfile with output histograms
  TString  barefile(fileName), chainFile, meta, histfile("b2u.root");
  if (file == 0) {
    if (barefile.Contains("chains/")) {
      meta = barefile;
      //      histfile = "outp/" + barefile.ReplaceAll("chains/", "") + ".root";
      histfile = barefile.ReplaceAll("chains/", "") + ".root";
      if (dirspec) histfile = dirName + "/" + histfile;
      //      meta = "chains/" + barefile;
    } else {
      meta = barefile;
      histfile =  barefile + ".root";
      if (dirspec) histfile = dirName + "/" + histfile;
    }

    // -- The following lines strip everything from the string up to and including the last '/'
    int fl = barefile.Last('/');
    TString bla(barefile);
    bla.Replace(0, fl+1, ' '); bla.Strip(TString::kLeading, ' ');  bla.Remove(0,1);
    histfile =  bla + ".root";
    if (dirspec) histfile = dirName + "/" + histfile;
  }  else if (file == 1) {
    // -- The following lines strip everything from the string up to and including the last '/'
    int fl = barefile.Last('/');
    TString bla(barefile);
    bla.Replace(0, fl+1, ' '); bla.Strip(TString::kLeading, ' ');  bla.Remove(0,1);
    histfile =  bla;
    if (dirspec) histfile = dirName + "/" + histfile;
  }
  
  TString histFile = histfile;

  cout << "==> b2u> Use " << histfile.Data() << " for output histograms" << endl;
  cout << "==> b2u> Use " << fileName.Data() << " for input" << endl;

  // -- Set up chain
  TChain *chain = new TChain(treeName);
  cout << "==> b2u> Chaining ... " << treeName << endl;
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
	cout << "    " << meta << endl; 
	chain->Add(meta); 
      }
    }
    is.close();
  }
  else if (file == 1) {
    cout << fileName << endl;
    chain->Add(fileName);
  }

  // -- Analysis 
  b2uNtp a(chain, isMC);
  a.setup(histfile.Data(), dump); 
  a.readCuts(cutFile); 

  a.fOptPidKillingEl = pidKillingEl; 
  a.fOptPidKillingMu = pidKillingMu; 
  a.fOptPidKillingKa = pidKillingKa; 
  a.fOptSmearNeut    = smearNeut; 

  a.fVerbose = isVerbose; 

  a.Loop(nevents, skipevents);
  a.closeHistFile(); 

  return 0;
    
}

