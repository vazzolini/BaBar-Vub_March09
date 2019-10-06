
//  $Id: anaRecoil.cc,v 1.23 2003/03/19 02:48:34 ursl Exp $
//
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

#include "VubAnalysis/recoilNtp.hh"


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int main(int argc, char *argv[]) {

  int processID = gSystem->GetPid();
  cout << "Running under process ID  " << processID << endl;

  TString progName  = argv[0]; progName.ReplaceAll("bin/", "");
  TString writeName, fileName, dirName, cutFile("");
  char mpFileName[80];

  TString treeName("h1");
  Int_t dump(0), write(0), writeG(0);
  Int_t file(0), lun(0);
  Int_t isMC(0), isVerbose(0), switchOff(0);
  Int_t dirspec(0);
  Int_t nevents(0), skipevents(0);
  Int_t blind(0), runGammas(0), runKlongs(0), runCategories(0), filterK0s(1010000), newFormat(2); 
  Int_t pidKilling(0), pidKillingKaon(0), pidKillingEl(0), pidKillingMu(0);
  Int_t isTauNu(0), isVcb(0), isDstar(0);
  Int_t isMakeParam(0);
  Int_t doUnBrem(0), vDotagstudy(0), generatemodelist(0), trackcutstudy(0);
  Int_t smearTracks(0);
  Int_t smearNeut(0);
  Int_t scaleKlongs(0); 
  Int_t oneProng(0);
  Int_t randomSeed(processID);

  // -- command line arguments
  for (int i = 0; i < argc; i++){
    if(strcmp(argv[i],"-b")  == 0) {blind      = 1;}                                             // skip events with Mx(fit) < 1.6 GeV
    if(strcmp(argv[i],"-C")  == 0) {cutFile    = TString(argv[++i]); }                           // file with cuts
    if(strcmp(argv[i],"-c")  == 0) {fileName   = TString(argv[++i]); file = 0; }                 // file with chain definition
    if(strcmp(argv[i],"-D")  == 0) {dirName     = TString(argv[++i]);  dirspec = 1; }            // where to put the output
    if(strcmp(argv[i],"-d")  == 0) dump        = atoi(argv[++i]);                                // printout
    if(strcmp(argv[i],"-Dstar")  == 0) isDstar = 1;
    if(strcmp(argv[i],"-f")  == 0) {fileName   = TString(argv[++i]); file = 1; }                 // single file instead of chain
    if(strcmp(argv[i],"-F")  == 0) {filterK0s  = atoi(argv[++i]);}                               // require 0 or x K0S->pi0pi0 decays
    if(strcmp(argv[i],"-g")  == 0) {runGammas  = atoi(argv[++i]);}                                             
    // run gamma study - splitoff studies 
    if(strcmp(argv[i],"-K")  == 0) {runKlongs  = 1;}                                             // run Klong study
    if(strcmp(argv[i],"-l")  == 0) {runCategories = 1;}                                          // run Categories
    if(strcmp(argv[i],"-MC") == 0)  isMC = 1;                                                    // flag to avoid warnings
    if(strcmp(argv[i],"-mc") == 0)  isMC = 1;                                                    // flag to avoid warnings
    if(strcmp(argv[i],"-mp") == 0) {strcpy(mpFileName,argv[++i]); isMakeParam = 1;}           // make parameters
    if(strcmp(argv[i],"-n")  == 0) nevents     = atoi(argv[++i]);                                // there is still a small bug in here
    if(strcmp(argv[i],"-p")  == 0) pidKilling  = 1;                                              // redo PidKilling on tuple (overall)
    if(strcmp(argv[i],"-pe")  == 0) pidKillingEl  = 1;                                           // redo PidKilling on tuple for el
    if(strcmp(argv[i],"-pk")  == 0) pidKillingKaon  = 1;                                         // redo PidKilling on tuple for kaons
    if(strcmp(argv[i],"-pm")  == 0) pidKillingMu  = 1;                                           // redo PidKilling on tuple for muons
    if(strcmp(argv[i],"-r")  == 0) randomSeed  = atoi(argv[++i]);                                // set seed for random number generator
    if(strcmp(argv[i],"-TAUNU")  == 0) isTauNu = 1;

    if(strcmp(argv[i],"-vcb")  == 0) isVcb = 1;
    if(strcmp(argv[i],"-vcbtagstudy")  == 0) vDotagstudy = 1;
    if(strcmp(argv[i],"-unbrem")  == 0) doUnBrem = 1;
    if(strcmp(argv[i],"-modelist")  == 0) generatemodelist = 1;
    if(strcmp(argv[i],"-trackcutstudy") == 0) trackcutstudy = 1;

    if(strcmp(argv[i],"-t")  == 0) treeName    = TString(argv[++i]);                             // ...
    if(strcmp(argv[i],"-ts") == 0)  lun = atoi(argv[++i]);                                       // timestamp business
    if(strcmp(argv[i],"-s")  == 0) skipevents  = atoi(argv[++i]);                                // skip events at the beginning
    if(strcmp(argv[i],"-ST") == 0) smearTracks = 1;                                              // smear tracks
    if(strcmp(argv[i],"-SN") == 0) smearNeut = 1;                                                // smear neut
    if(strcmp(argv[i],"-KL") == 0) scaleKlongs = 1;                              // scale KL energies
    if(strcmp(argv[i],"-v")  == 0)  isVerbose = 1;                                               // be verbose
    if(strcmp(argv[i],"-w")  == 0) {writeName = TString(argv[++i]); write = 1;}                  // write event list into red. file
    if(strcmp(argv[i],"-wG")  == 0) {writeName = TString(argv[++i]); writeG = 1;}                // write event list into red. file (after all cuts)
    if(strcmp(argv[i],"-y")  == 0)  newFormat = atoi(argv[++i]);                                 // flag to indicate the tree version
                                                                                                 // 0: Feb02, 1: apr 02, 2: aug02
    if(strcmp(argv[i],"-z")  == 0)  switchOff = 1;                                               // flag to speed up reading of tree
    if(strcmp(argv[i],"-1P")  == 0) oneProng  = 1;                                               // flag for one prong tau decays
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
    } else {
      //      meta = "chains/" + barefile;
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
  cout << "Opening " << histFile.Data() << " for output histograms" << endl;
  cout << "Opening " << fileName.Data() << " for input" << endl;
  
  
  // -- Set up chain
  TChain *chain = new TChain(treeName);
  cout << "Chaining ... " << treeName << endl;
  char pName[2000]; 
  int nentries; 
  if (file == 0) {
    ifstream is(meta);  
    while(meta.ReadLine(is) && (!meta.IsNull())){ 
      nentries = -1;
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
  
  recoilNtp a(chain,isMC, newFormat);
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
  if (!(write)) a.openHistFile(histFile);
    
  a.runGammas(runGammas); 
  a.doKlongScaling(scaleKlongs);
  a.doTrackSmearing(smearTracks);
  a.doNeutSmearing(smearNeut);
  a.filterK0s(filterK0s); 
  
  a.runKlongs(runKlongs); 
  a.runCategories(runCategories); 
  a.runBlind(blind); 
  if (isMakeParam) a.makeParam(mpFileName);
  
  if (write) a.makeEventList(1); 
  if (writeG) a.makeEventList(2); 
    
  if (!(write)) a.bookHist(dump);
  cout << "====================" << endl;
  cout << "isVerbose      = " << isVerbose << endl;
  cout << "runGammas      = " << runGammas << endl;
  cout << "runKlongs      = " << runKlongs << endl;
  cout << "PidKilling     = " << pidKilling << endl;
  cout << "PidKillingKaon = " << pidKillingKaon << endl;
  cout << "TrkSmearing    = " << smearTracks << endl;
  cout << "NeuSmearing    = " << smearNeut << endl;
  cout << "Filter K0S     = " << filterK0s << endl;
  cout << "categories     = " << runCategories << endl;
  cout << "blind          = " << blind << endl;
  cout << "write          = " << write << endl;
  cout << "isMakeParam    = " << isMakeParam << endl;
  cout << "====================" << endl;

  if (switchOff) {
    a.switchOffReading("Vtx*");
    a.switchOffReading("Cov*");
    // leave the following in when doing smearing ...
    // a.switchOffReading("cov*");  
    // a.switchOffReading("ppcov*");
    
    a.switchOffReading("ecal*");
    
    a.switchOffReading("Ifr*");
    
    a.switchOffReading("Thrust*");
    a.switchOffReading("thrust*");
    a.switchOffReading("beamSCov*");
    
    a.switchOffReading("Doca*");
    
    a.switchOffReading("st*");
    a.switchOffReading("dof*");
    a.switchOffReading("chi2*");
    
    a.switchOffReading("d3*");
    a.switchOffReading("d4*");
    a.switchOffReading("d5*");
    a.switchOffReading("d6*");
    a.switchOffReading("d7*");
    a.switchOffReading("s2*");
    
  }
  
  if (!(write)) a.Loop(nevents, skipevents, isVerbose, lun);
  if (write) {
    a.Skim(.8,nevents, skipevents,isVerbose,"");
    a.dumpEventList(writeName.Data());
  }
  if (writeG) {
    a.dumpEventList(writeName.Data());
  }
  if (!(write)) a.closeHistFile();
  
  return 0;
  
}

