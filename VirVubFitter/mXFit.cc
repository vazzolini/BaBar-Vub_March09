
#include <iostream>
#include <vector>

#include <math.h>
#include <stdlib.h>

#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TString.h"
#include "TStopwatch.h"
#include "TStyle.h"

#include "VirVubFitter/mXClass.hh"
#include "VirVubFitter/CMClass.hh"
#include "VirVubFitter/PidCorrectMesMean.hh"

class RooDataSet;

using namespace std;

void initGraphicsStyle(TStyle* pStyle);

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int main(int argc, char *argv[]) {
  
  gROOT->SetBatch(kTRUE);

  initGraphicsStyle(gStyle);

  TString texPrefix("all"); 
  TString progName  = argv[0]; progName.ReplaceAll("bin/", "");
  TString fileName, optionFile(""),cutFile(""),wFile("");
  TString mesFile(""),prefixOut(""),dirName(""), UFile("");
  TString mxFile1d(""), wFermiFile("");

  int sys=0, sysD=0, fermi=0, theosys=0, mu=0, runflag=0, rel=18, comp=0;
  bool iscm2(true), varfit(false), endpointCor(false);
  std::vector<float> wFermivec;

  int unf=0, munf=0, mx2unf=0;
  double hiunf=0.;

  // =================== command line arguments ==============================

  for (int i = 1; i < argc; i++){

    if(strcmp(argv[i],"-Sys")     == 0) { sys = atoi(argv[++i]); }              // int flag for systematics on B BF reweighting (1=,2=,other=RND for 2)
    if(strcmp(argv[i],"-SysD")     == 0) { sysD = atoi(argv[++i]); }              // int flag for systematics on D BF reweighting (1=,2=,other=RND for 2)
    if(strcmp(argv[i],"-varfit")  == 0) { varfit = (bool)atoi(argv[++i]); }     // flag for using fitted kinematic variables (mx/q2/pplus)
    if(strcmp(argv[i],"-epCor")   == 0) { endpointCor = true; }                 // apply endpoint correction on data
    if(strcmp(argv[i],"-RunFl")   == 0) { runflag = atoi(argv[++i]); }          // specify run periods (e.g. 12, 3, 4, 14)
    if(strcmp(argv[i],"-cm")      == 0) { if(strcmp(argv[++i],"CM1")==0) iscm2 = false; } // choose ntuple style (CM1/CM2)
    if(strcmp(argv[i],"-rel")     == 0) { rel = atoi(argv[++i]); if(rel==14||rel==18||rel==22){}else{cout << "Warning: rel " << rel << " not supported, using rel 22 assumption..." << endl; rel=22;}} // choose ntuple production release for CM2 (14/18) this determines the reweightings used

    if(strcmp(argv[i],"-F")       == 0) { optionFile = TString(argv[++i]); }    // filename for options file  
    if(strcmp(argv[i],"-W")       == 0) { wFile = TString(argv[++i]); }         // filname for hybrid weights
    if(strcmp(argv[i],"-Mes")     == 0) { mesFile = TString(argv[++i]); }       // filename for mes fit parameter file
    if(strcmp(argv[i],"-s")       == 0) { texPrefix = TString(argv[++i]); }     // prefix for output tex file
    if(strcmp(argv[i],"-P")       == 0) { prefixOut = TString(argv[++i]); }     // prefix for all output files
    if(strcmp(argv[i],"-C")       == 0) { cutFile = TString(argv[++i]); }       // filename for cuts file 
    if(strcmp(argv[i],"-D")       == 0) { dirName = TString(argv[++i]); }       // directory for output files
    if(strcmp(argv[i],"-mx1d")    == 0) { mxFile1d = TString(argv[++i]); }      // file with binning for mxhad (1d: mx)	  
    if(strcmp(argv[i],"-FFile")   == 0) { wFermiFile = TString(argv[++i]); }    // filename for fermi weights

    if(strcmp(argv[i],"-U")       == 0) { UFile = TString(argv[++i]); }         // filename for unfolding parameters
    if(strcmp(argv[i],"-unf")     == 0) { unf = atoi(argv[++i]); }              // Number of bins for unfolding, KT
    if(strcmp(argv[i],"-hiunf")   == 0) { hiunf = atof(argv[++i]); }            // Endpoint for unfolding, KT
    if(strcmp(argv[i],"-mx2unf")  == 0) { mx2unf = 1;  }
    if(strcmp(argv[i],"-mu")      == 0) { mu = atoi(argv[++i]); }               // cut on generated multiplicity category
    if(strcmp(argv[i],"-mult")    == 0) { munf = 1; }                           // Flag unfolding for multiplicity categories, KT
     if(strcmp(argv[i],"-comp")   == 0) { comp = 1;}                           //do data-MC comparison
  }
    
  CMClass *cmctemp=new CMClass(iscm2,varfit,wFermivec);
  if (wFermiFile.Length() > 0) {
    ReadwFermiFile(wFermiFile,wFermivec,runflag);
    if (wFermivec.size()==0) {
      std::cout << "ERROR wfermivec not initialized. Check -FFile flag" << std::endl; 
      return 1;
    }
  }

  // ==============  Create the big fat one do-it-all-object here. ===========

  mXClass a((TTree*)gDirectory->Get(cmctemp->ev),wFile,sys,sysD,unf,hiunf,mx2unf,mu,iscm2,varfit,wFermivec,rel,comp);
  delete cmctemp;

  //  a.applyEndpointCor(endpointCor, runflag);

  // =============== Read files from thefiles/ dir ===========================

  if (optionFile.Length() > 0) a.readOptions(optionFile, 1);
  
  // ==============  Read cuts for ascii file in theset/ dir =================
  
  if (cutFile.Length() > 0) a.readCuts(cutFile, 1, runflag); 
  a.dumpWFermiFile();

  // ==============  Set the prefix for output storage =======================

  if (prefixOut.Length() > 0) a.setPrefix(prefixOut);
  a.setTexPrefix(texPrefix); 

  // ================ Read the mES parameters obtained from global fits 
  //This helps reducing instability of mES fits.      ========================

  std::cout << "now readmesParam" << std::endl;
  if (mesFile.Length() > 0) a.readmesParam(mesFile, 1);
  
  if (mxFile1d.Length() > 0) a.mxBinning1d(mxFile1d); 
  a.InitBinning();

  //================= Initialize the PidCorrectMesMean for the mes corrections 

  cout << "Initializing mES endpoint corrections....... " << endl;
  if(rel==22){
    PidCorrectMesMean::initialize("mES_means_Runs123456_final.txt");
  }else if(rel==18){
    PidCorrectMesMean::initialize("mES_means_Runs12345.txt");
  }else{
    cout << "No mES endpoint corrections available for release " << rel << endl;
  }

  
  // ================ Read the parameters for the unfolding histos 
  if (unf) {
    if (UFile.Length() > 0) a.readUnfParam(UFile, 1);
  }
  if (dirName.Length() > 0) a.setDirectory(dirName);

  // ================ Prepare histograms 
  a.openHistFile("Fitres.root");
  a.Bookhist();

  /* *************************************************************************

    now prepare datasets 
     - either read ntuples and loop over them
       (with the option to write out RooDataSets (set writeDataDir))
     - or read directly RooDataSets (set readDataDir)  

  ************************************************************************* */

  TStopwatch timer; // for timing purpose  

  bool *isc = a.getfilechain();

  // =============  Real data stuff     //
  //read reduced dataset if corresponding option is set
  //Loop o'er entries and write dataset to rootfile
  std::cout << std::endl << "Loading data" << std::endl;
  a.Init(a.getchain(a.getfile(mXClass::Data)));
  if (a.Loop(1,1,0,0,0) == 0) {  //data
    std::cout << "Error: Zero events in data chain. Stopping execution!" << std::endl;
    exit(EXIT_FAILURE);
  }

  // ================  VCB stuff   //
  //read generic datasets if corresponding option is set
    
    //Loop o'er entries and write dataset to rootfile
  std::cout << std::endl << "Loading MC: Vcb " << std::endl;
  a.Init(a.getchain(a.getfile(mXClass::Vcb)));
  a.Loop(0,1,1,0,0);     //generic 
  
  // ====================  VUBish stuff   //
  //read generic datasets if corresponding option is set

  //Loop o'er entries and ... 
  std::cout << std::endl << "Loading MC: Vub nonres" << std::endl;
  a.Init(a.getchain(a.getfile(mXClass::VubTotalnres)));
  cout << "COMP: Are we looping over the nonres signal?" << endl;
  if(comp<1){
    cout << "COMP: We are..." << endl;
    a.Loop(0,0,1,1,0);     //signal non resonant
  }
  
  //   to commment for MC study
  std::cout << std::endl << "Loading MC: Vub res" << std::endl;
  a.Init(a.getchain(a.getfile(mXClass::VubTotalres)));
  cout << "COMP: Are we looping over the res signal?" << endl;
  if(comp<1){
    cout << "COMP: We are..." << endl;
    a.Loop(0,0,1,0,0);     //signal resonant
  }


  
  //Also read the truth MC that is needed for the spectral unfolding
  if (unf) {
    //Loop o'er entries and write dataset to rootfile
    TFile f8;
    if (!isc[mXClass::VubTruthres]) {
      if (f8.Open(a.getfile(mXClass::VubTruthres)) == 0) {
	std::cout << "Error: Can't open file " << a.getfile(mXClass::VubTruthres) << ". Stopping execution!" << std::endl;
	exit(EXIT_FAILURE);
      }
      a.InitTruth((TTree*)gDirectory->Get("h6"));
    } else {
      a.InitTruth(a.getchain(a.getfile(mXClass::VubTruthres)));
    }
    std::cout << std::endl << "Loading MC: Vub res truth" << std::endl;
    a.Loop(0,0,1,0,1);     //signal resonant
      
    TFile f9;
    if (!isc[mXClass::VubTruthnres]) {
      if (f9.Open(a.getfile(mXClass::VubTruthnres)) == 0) {
	std::cout << "Error: Can't open file " << a.getfile(mXClass::VubTruthnres) << ". Stopping execution!" << std::endl;
	exit(EXIT_FAILURE);
      }
      a.InitTruth((TTree*)gDirectory->Get("h6"));
    } else {
      a.InitTruth(a.getchain(a.getfile(mXClass::VubTruthnres)));
    }
    std::cout << std::endl << "Loading MC: Vub nres truth" << std::endl;
    a.Loop(0,0,1,1,1);     //signal nonresonant
  }
    
  // print out timer
  std::cout << "Timer for fitMes: "; timer.Print();

  /* *************************************************************************

    prepare and do fitting:
     - fit yields for the different dataset component
     - do final fit

  ************************************************************************* */

  // =========== This to use the histos file 
  a.fHistFile->cd();

  // =========== Now mes fits 
  std::vector<TString> dataSetNames(4);
  dataSetNames[0] = "data";
  dataSetNames[1] = "vcb"; dataSetNames[2] = "vub"; dataSetNames[3] = "other";

  // Fit Mes 1D distributions for cont subtraction 
    
  if(munf){ //Fit Mes for multiplicity categories
    a.FitMes(dataSetNames[0], 1);
    a.FitMes(dataSetNames[1], 1);
    a.FitMes(dataSetNames[2], 1);        
    a.FitMes(dataSetNames[3], 1);  
  }

  a.FitMes(dataSetNames[0], 0);
  a.FitMes(dataSetNames[1], 0);
  a.FitMes(dataSetNames[2], 0);        
  a.FitMes(dataSetNames[3], 0); 
  
  a.theFit(munf);

  // ========== unfolding 
  if (unf) {
    std::cout << "Now the unfolding part" << std::endl;

    if(munf) a.MultCorr();
    a.UnfHistos(munf);
  }

  // ========= clean up 
  a.closeHistFile();

  std::cout << std::endl << "Final message: fits successfully completed!" << std::endl;

  return 0;
}


 void initGraphicsStyle(TStyle* pStyle)
{
  std::cout  << "Adjusting ROOT graphics environment!" << std::endl;

  // no border and color white
  pStyle->SetFrameBorderMode(0);
  pStyle->SetFrameBorderSize(0);

  pStyle->SetDrawBorder(0);

  pStyle->SetTextFont(22);

  // canvas
  pStyle->SetCanvasBorderMode(0);
  pStyle->SetCanvasBorderSize(0);
  pStyle->SetCanvasColor(0);

  // pads
  pStyle->SetPadBorderMode(0);             //                                 - default 1
  pStyle->SetPadBorderSize(0);             //                                 - default 1
  pStyle->SetPadColor(kWhite);             // set pad background color        - default 19

  // histogram title
  pStyle->SetTitleBorderSize(1);           // border size of histogram title  - default 2
  pStyle->SetTitleFillColor(kWhite);       // histogram title box fill color  - default 19
  pStyle->SetTitleFont(22);                // histogram title font            - default 62

  pStyle->SetLabelFont(22,"xyz");

  // histogram stat box
  pStyle->SetStatBorderSize(1);            // stat box border size              - default 2 (1 for no shadow)
  pStyle->SetStatColor(0);                 // stat box background color         - default 19
  pStyle->SetStatFont(22);                 // stat box font type                - default 62
  pStyle->SetStatTextColor(kBlack);        // stat box text color               - default kBlack

  return;
}

