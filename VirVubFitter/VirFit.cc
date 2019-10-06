
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

#include "VirVubFitter/VirClass.hh"
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
  TString q2File(""), mxFile(""), mxFile1d(""), wFermiFile(""), pplusFile1d(""), dssFile("");
  TString writeDataDir(""), readDataDir(""), sigpeakcorrmxfile(""), readPdfFile("");

  int sys=0, fermi=0, theosys=0, q2=0, comb=0, Sun=0, ck=0, me=0, mu=0, newbin=0, runflag=0, rel=18, bsys=-99;
  int SPseed = 0;
  double wsys=0.;
  bool iscm2(true), varfit(false), endpointCor(false), fixsigpeakmesratio(false), countMC(false), thecmp(false);
  std::vector<float> wFermivec;

  int unf=0, munf=0, mx2unf=0;
  double hiunf=0.;

  double dssRatio=1.;

  int NORE=0, RE=0;

  // =================== command line arguments ==============================

  for (int i = 1; i < argc; i++){

    if(strcmp(argv[i],"-comb")    == 0) { comb = 1; }                           // flag for 2d analysis
    if(strcmp(argv[i],"-Sun")     == 0) { Sun = 1; }                            // signal unfolding using Bauer et al.
    if(strcmp(argv[i],"-ckBauer") == 0) { ck = 1; }                             // use DFN instead of default BLL for 2D efficiency unfolding
    // if(strcmp(argv[i],"-nore")  == 0) {NORE = 1; }
    // if(strcmp(argv[i],"-re")  == 0) {RE = 1; }
    if(strcmp(argv[i],"-me")      == 0) { me = atoi(argv[++i]); }               // flag for systematics: variation of floating parameters
    if(strcmp(argv[i],"-newbin")  == 0) { newbin = 1; }                         // flag for new first bin in 1d unfolding
    if(strcmp(argv[i],"-q2")      == 0) { q2 = atoi(argv[++i]); }               // int flag for 1D variable: 0=mx, 1=q2, 2=pplus (->FITQ2)
    if(strcmp(argv[i],"-Sys")     == 0) { sys = atoi(argv[++i]); }              // int flag for systematics on B/D BF reweighting (1=,2=,other=RND for 2)
    if(strcmp(argv[i],"-SPseed")  == 0) { SPseed = atoi(argv[++i]); }           // int flag for randomizing S/P ratios in mES fits
    if(strcmp(argv[i],"-varfit")  == 0) { varfit = (bool)atoi(argv[++i]); }     // flag for using fitted kinematic variables (mx/q2/pplus)
    if(strcmp(argv[i],"-epCor")   == 0) { endpointCor = true; }                 // apply endpoint correction on data
    if(strcmp(argv[i],"-count")   == 0) { countMC = true;}                      // count MC events on truth-matched samples instead of fitting mES

    if(strcmp(argv[i],"-RunFl")   == 0) { runflag = atoi(argv[++i]); }          // specify run periods (e.g. 12, 3, 4, 14, 15)
    if(strcmp(argv[i],"-cm")      == 0) { if(strcmp(argv[++i],"CM1")==0) iscm2 = false; } // choose ntuple style (CM1/CM2)
    if(strcmp(argv[i],"-rel")      == 0) { rel = atoi(argv[++i]); if(rel==14||rel==18||rel==22){cout << "RELEASE: " << rel << endl;}else{cout << "Warning: rel " << rel << " not supported, using rel 18 assumption..." << endl; rel=18;}} // choose ntuple production release for CM2 (14/18) this determines the reweightings used

    if(strcmp(argv[i],"-F")       == 0) { optionFile = TString(argv[++i]); }    // filename for options file  
    if(strcmp(argv[i],"-W")       == 0) { wFile = TString(argv[++i]); }         // filname for hybrid weights
    if(strcmp(argv[i],"-Mes")     == 0) { mesFile = TString(argv[++i]); }       // filename for mes fit parameter file
    if(strcmp(argv[i],"-s")       == 0) { texPrefix = TString(argv[++i]); }     // prefix for output tex file
    if(strcmp(argv[i],"-P")       == 0) { prefixOut = TString(argv[++i]); }     // prefix for all output files
    if(strcmp(argv[i],"-C")       == 0) { cutFile = TString(argv[++i]); }       // filename for cuts file 
    if(strcmp(argv[i],"-wisys")   == 0) { wsys = atof(argv[++i]); }             // systematic scaling factor for hybrid weights
    if(strcmp(argv[i],"-D")       == 0) { dirName = TString(argv[++i]); }       // directory for output files
    if(strcmp(argv[i],"-q2V")     == 0) { q2File = TString(argv[++i]); }        // file with binning for q2 (2d: mx/q2)	  
    if(strcmp(argv[i],"-mxV")     == 0) { mxFile = TString(argv[++i]); }        // file with binning for mxhad (2d: mx/q2)
    if(strcmp(argv[i],"-mx1d")    == 0) { mxFile1d = TString(argv[++i]); }      // file with binning for mxhad (1d: mx)	  
    if(strcmp(argv[i],"-pplus1d") == 0) { pplusFile1d = TString(argv[++i]); }   // file with binning for pplus (1d: pplus)
    if(strcmp(argv[i],"-FFile")   == 0) { wFermiFile = TString(argv[++i]); }    // filename for fermi weights
    if(strcmp(argv[i],"-sdDir")   == 0) { writeDataDir = TString(argv[++i]); }  // directory for writing RooDataSet
    if(strcmp(argv[i],"-rdDir")   == 0) { readDataDir = TString(argv[++i]); }   // directory for reading RooDataSet
    if(strcmp(argv[i],"-readpdftree")   == 0) { readPdfFile = TString(argv[++i]); }  // directory for reading RooDataSet toy mes
    if(strcmp(argv[i],"-dssRatio")== 0) { dssRatio = atof(argv[++i]); }         // fix D**lnu+nonSL ratio wrt Dlnu+D*lnu 
    if(strcmp(argv[i],"-dssFile") == 0) { dssFile = TString(argv[++i]);}        // reweight D** bin-by-bin in kinematic variable
    if(strcmp(argv[i],"-bsys") == 0)    { bsys = atoi(argv[++i]);}              // B decays systematics

    if(strcmp(argv[i],"-U")       == 0) { UFile = TString(argv[++i]); }         // filename for unfolding parameters
    if(strcmp(argv[i],"-unf")     == 0) { unf = atoi(argv[++i]); }              // Number of bins for unfolding, KT
    if(strcmp(argv[i],"-hiunf")   == 0) { hiunf = atof(argv[++i]); }            // Endpoint for unfolding, KT
    if(strcmp(argv[i],"-mx2unf")  == 0) { mx2unf = 1;  }
    if(strcmp(argv[i],"-mu")      == 0) { mu = atoi(argv[++i]); }               // cut on generated multiplicity category
    if(strcmp(argv[i],"-mult")    == 0) { munf = 1; }                           // Flag unfolding for multiplicity categories, KT
    //    if(strcmp(argv[i],"-mesfitmodel") == 0) { mesFitModel = atoi(argv[++i]); }
    if(strcmp(argv[i],"-sigpeakcorrmx") == 0) {sigpeakcorrmxfile=TString(argv[++i]);} //File with correction ratio signal/peaking BKG
    if(strcmp(argv[i],"-thecomparison") == 0) { thecmp = true; }            //Flag for thecomparison dataset writing.
    if(strcmp(argv[i],"-Fermi")   == 0) { fermi = 1;}                           // obsolete
  }
    
  CMClass *cmctemp=new CMClass(iscm2,varfit,wFermivec);
  if (wFermiFile.Length() > 0) {
    ReadwFermiFile(wFermiFile,wFermivec,runflag);
    if (wFermivec.size()==0) {
      std::cout << "ERROR wfermivec not initialized. Check -FFile flag" << std::endl; 
      return 1;
    }
  }

  bool isWD = false;
  if (writeDataDir.Length() > 0) isWD = true;

  // ==============  Create the big fat one do-it-all-object here. ===========

  VirClass a((TTree*)gDirectory->Get(cmctemp->ev),wFile,sys,q2,comb,Sun,unf,hiunf,mx2unf,me,mu,iscm2,varfit,wFermivec,newbin,SPseed,isWD,rel,bsys,thecmp);
  delete cmctemp;

  a.setCountingFlag(countMC);

  //set the relative contribution of D**+nonSL wrt D and D*
  if(Sun != 1) dssRatio = 1;
  std::cout << "Applying D** correction ratio = " << dssRatio << std::endl;
  if(dssFile.Length() > 0) 
    a.setDssRatio(dssFile);
  else
    a.setDssRatio(dssRatio);
  
  a.applyEndpointCor(endpointCor, runflag);

  // =============== Read files from thefiles/ dir ===========================

  if (optionFile.Length() > 0) a.readOptions(optionFile, 1);
    
  // ==============  Read cuts for ascii file in theset/ dir =================
  if (cutFile.Length() > 0) a.readCuts(cutFile, 1, wsys, runflag, rel); 
  a.dumpWFermiFile();

  // ==============  Set the prefix for output storage =======================

  if (prefixOut.Length() > 0) a.setPrefix(prefixOut);
  a.setTexPrefix(texPrefix); 

  // ================ Read the mES parameters obtained from global fits 
  //This helps reducing instability of mES fits.      ========================

  std::cout << "now readmesParam" << std::endl;
  if (mesFile.Length() > 0) a.readmesParam(mesFile, 1);
  
  if (q2File.Length() > 0) a.q2Binning(q2File); 
  if (mxFile.Length() > 0) a.mxBinning(mxFile); 
  if (mxFile1d.Length() > 0) a.mxBinning1d(mxFile1d); 
  if (pplusFile1d.Length() > 0) a.pplusBinning1d(pplusFile1d); 
  a.InitBinning(comb);

  //================= Initialize the PidCorrectMesMean for the mes corrections 

  if(! a.isfitMC())
    PidCorrectMesMean::initialize("mES_means_Runs123456_final.txt");
  
  // ================ Read the correction ratio  signal/peaking BKG for mx analysis 
  if (sigpeakcorrmxfile.Length() > 0) {
    a.readmxcorrratiosigpeak(sigpeakcorrmxfile);
    fixsigpeakmesratio = true;
  }

  // ================ Read the parameters for the unfolding histos 
  if (unf) {
    if (UFile.Length() > 0) a.readUnfParam(UFile, 1);
  }
  if (dirName.Length() > 0) a.setDirectory(dirName);

  // ================ Prepare histograms 
  a.openHistFile("Fitres.root");
  a.Bookhist();

  // =============== Set readir and writedir for dataset root file 
  if (readDataDir.Length() > 0) a.setDataDir(readDataDir);
  if (writeDataDir.Length() > 0) a.setDataDir(writeDataDir);
  if (readPdfFile.Length() > 0)  a.openPdfToyFile(readPdfFile);

  /* *************************************************************************

    now prepare datasets 
     - either read ntuples and loop over them
       (with the option to write out RooDataSets (set writeDataDir))
     - or read directly RooDataSets (set readDataDir)  

  ************************************************************************* */

  TStopwatch timer; // for timing purpose  

  std::cout << "isfitMC: " << a.isfitMC() << std::endl;
  bool *isc = a.getfilechain();
  std::vector<RooDataSet*> idataset;

  if(a.isfitMC()){
    // data file from MC...
    cout<<" This section has been commented out by Antonio on 21/feb/2007 to use properly the flag isfitmc"<<endl;
    cout<<" -fitMC automatically selects Vcb chainfile as data"<<endl;

    /* COMMENT START HERE

    // first Vcb ...
    std::cout << "Entering in Vcb as data" << std::endl;
    

    a.Init(a.getchain(a.getfile(VirClass::Vcb1)));
    a.Loop(1,1,1,0,comb,0);   //generic  as  MC Bch

    a.Init(a.getchain(a.getfile(VirClass::Vcb2)));
    a.Loop(1,1,1,0,comb,0);   //generic  as  MC BO

    // ... and then Vub
    std::cout << "Entering in Vub as data" << std::endl;

    a.Init(a.getchain(a.getfile(VirClass::VubTotalnres)));
    a.Loop(1,0,1,1,comb,0);   //signal non resonant  as  MC
   
    // Why VubTotalres is missing? Roberto said: Maybe historic?
    */   // 21/Feb/2007 Comment ends here.

  } // else {
    
    // =============  Real data stuff     //
    //read reduced dataset if corresponding option is set

  if (readDataDir.Length() > 0) {
    a.readDataFile(&(a.datadata),"DATA");
  } else if ( a.READPDFTREE ) ;
  else {
    //Loop o'er entries and write dataset to rootfile
    
    std::cout << std::endl << "Loading data" << std::endl;
    a.Init(a.getchain(a.getfile(VirClass::Data)));
    if (a.Loop(1,1,0,0,comb,0) == 0) {  //data
      std::cout << "Error: Zero events in data chain. Stopping execution!" << std::endl;
      exit(EXIT_FAILURE);
    }
  }  
  cout << "SURVIVORS TO lpYes CUT " << a.mySL/a.LUMI_DATA << endl;
  cout << "SURVIVORS TO lpYesSig CUT " << a.mylpYesSig/a.LUMI_DATA << endl;
  cout << "SURVIVORS TO lpYesSig && nuCut CUT "<<a.mynucut/a.LUMI_DATA << endl;
  cout << "SURVIVORS TO lpYesSig && nuCut && ch< " << a.CHHIGH<<" && ch > " << a.CHLOW << " CUT " << a.mychcut/a.LUMI_DATA << endl;
  cout << "SURVIVORS TO lpYesSig && nuCut && ch< " << a.CHHIGH<<" && ch > " << a.CHLOW << " && ksele == " << a.DEPL << " CUT " << a.mydepl/a.LUMI_DATA << endl;
  cout << "SURVIVORS TO lPYesSig && nuCut && ch < " << a.CHHIGH<<" && ch > " << a.CHLOW << " &&  ksele == " << a.DEPL<<" && !(WdeltaCut) CUT "
       << a.myWdeltaM/a.LUMI_DATA << endl;
  cout << "SURVIVORS TO ALL CUT " << a.myAllcut/a.LUMI_DATA << endl;
  
  //Fill-in idataset with dataset data
  if (writeDataDir.Length() > 0) {
    idataset.push_back(a.datadata);
  } 
  
  
  //CB test
  //readDataDir.Resize(0);
  
  // ================  VCB stuff   //
  //read generic datasets if corresponding option is set
  
  if (readDataDir.Length() > 0) {
    a.readDataFile(&(a.datamcvcb),"Vcb");
    a.readDataFile(&(a.pstarsample),"PStar");
    if (Sun) {
      a.readDataFile(&(a.datavcboth),"VcbOth");
    } else {
      a.readDataFile(&(a.datamcoth),"Oth");
    }
  } else if ( a.READPDFTREE ) ;
  else {
    //Loop o'er entries and write dataset to rootfile
    std::cout << std::endl << "Loading MC: Vcb " << std::endl;
    a.Init(a.getchain(a.getfile(VirClass::Vcb)));
    a.Loop(0,1,1,0,comb,0);     //generic 
  }
  //Now Fill-in idataset with dataset vcb/pstar
  if (writeDataDir.Length() > 0) { 
    //Now Fill-in idataset with dataset vcb
    idataset.push_back(a.datamcvcb);
    idataset.push_back(a.pstarsample);
    if (Sun) {
      idataset.push_back(a.datavcboth);
    } else {
      idataset.push_back(a.datamcoth);
    }
  }
    
  // ====================  VUBish stuff   //
  //read generic datasets if corresponding option is set

  if (readDataDir.Length() > 0) {
    if (Sun) {
      a.readDataFile(&(a.datavubin),"VubIN");
      a.readDataFile(&(a.datavubout),"VubOUT");
    } else {
      a.readDataFile(&(a.datamcvub),"Vub");
    }      
  } else if ( a.READPDFTREE ) ;
  else { 
      //Loop o'er entries and ... 
      std::cout << std::endl << "Loading MC: Vub nres" << std::endl;
      a.Init(a.getchain(a.getfile(VirClass::VubTotalnres)));
      a.Loop(0,0,1,1,comb,0);     //signal non resonant
      
      //   to commment for MC study
      std::cout << std::endl << "Loading MC: Vub res" << std::endl;
      a.Init(a.getchain(a.getfile(VirClass::VubTotalres)));
      a.Loop(0,0,1,0,comb,0);     //signal resonant
  }
    //Now Fill-in idataset with dataset vub
  if (writeDataDir.Length() > 0) { 
    if (Sun) {
      idataset.push_back(a.datavubin);
      idataset.push_back(a.datavubout);
    } else {
      idataset.push_back(a.datamcvub);
    }
  }  
  
  
  //Also read the truth MC that is needed for the spectral unfolding
  if (unf) {
    //read generic datasets if corresponding option is set
    if (readDataDir.Length() > 0) {
      a.readDataFile(&(a.unfmcvub),"VubUnf");
      a.readDataFile(&(a.unftmcvub),"VubTru");
    } else { 
      //Loop o'er entries and write dataset to rootfile
      TFile f8;
      if (!isc[VirClass::VubTruthres]) {
	if (f8.Open(a.getfile(VirClass::VubTruthres)) == 0) {
 	  std::cout << "Error: Can't open file " << a.getfile(VirClass::VubTruthres) << ". Stopping execution!" << std::endl;
	  exit(EXIT_FAILURE);
	}
	a.InitTruth((TTree*)gDirectory->Get("h6"));
      } else {
	a.InitTruth(a.getchain(a.getfile(VirClass::VubTruthres)));
      }
      std::cout << std::endl << "Loading MC: Vub res truth" << std::endl;
      a.Loop(0,0,1,0,comb,1);     //signal resonant
      
      TFile f9;
      if (!isc[VirClass::VubTruthnres]) {
	if (f9.Open(a.getfile(VirClass::VubTruthnres)) == 0) {
	  std::cout << "Error: Can't open file " << a.getfile(VirClass::VubTruthnres) << ". Stopping execution!" << std::endl;
	  exit(EXIT_FAILURE);
	}
	a.InitTruth((TTree*)gDirectory->Get("h6"));
      } else {
	a.InitTruth(a.getchain(a.getfile(VirClass::VubTruthnres)));
      }
      std::cout << std::endl << "Loading MC: Vub nres truth" << std::endl;
      a.Loop(0,0,1,1,comb,1);     //signal nonresonant

      //Now Fill-in idataset with dataset vub
      if (writeDataDir.Length() > 0) { 
	idataset.push_back(a.unfmcvub);
	idataset.push_back(a.unftmcvub);
      }
    }
  }
  
  if(a.isfitMC()){
    a.CreateNewPstarSample(writeDataDir.Length() > 0);
    cout<<" WARNING !!! SETTING DATA DATASET AS NEW PSTARSAMPLE "<<endl;
    a.datadata=a.pstarsamplesum;
  }

  // ========= Now write all these fine datasets     
  if (writeDataDir.Length() > 0) {
    a.writeDataFile(idataset);
    std::cout << "VirFit: Writing datasets and exiting" << std::endl;
    return 0;
  }
    

  //  finish comment for MC study

  //================== CHECK MC TRUTH Numbers and write them in a file=======
  
  //  cout<<"...---... ...---... ...---... NOW ENTERING section to get Numbers for BRBR from MC for DEBUG puroposes"<<endl;
  //  a.GetMCTrueNumbers(comb);
  //  cout<<"...---... ...---... ...---... NOW EXITING section to get Numbers for BRBR from MC for DEBUG puroposes"<<endl;
  //  return 1;

  //========================================================================

  // print out timer
  std::cout << "Timer for fitMes: "; timer.Print();

  /* *************************************************************************

    prepare and do fitting:
     - fit yields for the different dataset component
     - do final fit

  ************************************************************************* */

  // ======== Calculate peaking background correction fitting Generic MC
  // PEAKING BACKGROUND IS EVALUATED FOR SL CUTS & ALL CUTS
  // For EACH Bin of Kinematic (1D) Variable. For 2D case we evaluate on mx cut only.
    
  if( a.SUBTRACTPEAKING != 0 )
    a.evaluatePeakingBackground(comb,bool(a.SUBTRACTPEAKING-1));

  if( a.BTYPE == 2 )
    // ========= This computes the correction ratio B+/B0
    a.computeChargeCorr(Sun,comb);
  else
    // ========= CALCULATE CrossFeed Matrix
    a.CalculateCFMatrix(Sun,comb);

  // ===========  Pstar Factor 
  if(!(unf)) a.calcPstarFact(comb, Sun, ck, runflag); 
  
  //CB  a.Debug();
  // =========== This to use the histos file 
  a.fHistFile->cd();

  // =========== Now mes fits 
  std::vector<TString> dataSetNames(4);
  dataSetNames[0] = "data";
  if (Sun) {
    dataSetNames[1] = "vubin"; dataSetNames[2] = "vubout"; dataSetNames[3] = "vcboth";
  } else {
    dataSetNames[1] = "vcb"; dataSetNames[2] = "vub"; dataSetNames[3] = "other";  }



  if (!comb) {
    // Fit Mes 1D distributions for cont subtraction 
    
    if(munf){ //Fit Mes for multiplicity categories
      a.FitMes(dataSetNames[0], 1, 0, fixsigpeakmesratio);
      a.FitMes(dataSetNames[1], 1, 0 );
      a.FitMes(dataSetNames[2], 1, 0 );        
      a.FitMes(dataSetNames[3], 1, 0 );  
    }

    a.FitMes(dataSetNames[0], 0, Sun, fixsigpeakmesratio);
    a.FitMes(dataSetNames[1], 0, Sun, fixsigpeakmesratio);
    a.FitMes(dataSetNames[2], 0, Sun, fixsigpeakmesratio);        
    a.FitMes(dataSetNames[3], 0, Sun, fixsigpeakmesratio); 
   
  } else {
    //Fit Mes 2D distributions for cont subtraction
    
    a.FitMes2D(dataSetNames[0], Sun);
    a.FitMes2D(dataSetNames[1], Sun);
    a.FitMes2D(dataSetNames[2], Sun);        
    a.FitMes2D(dataSetNames[3], Sun); 

  }

  //  //CB plots for BAD
  //  if(!comb){
  //    a.FitPlots("data",0,0);
  //    a.FitPlots("data",0,1);
  //  }else{
  //    a.FitPlots("data",1,0);
  //    a.FitPlots("data",1,1);
  //  }

  // Close file opened for debug

  a.theFit(comb, munf, Sun, ck);

  // ========== unfolding 
  if (unf) {
    std::cout << "Now the unfolding part" << std::endl;

    if(munf) a.MultCorr();
    a.UnfHistos(munf);
  }

  if(a.SAVEPDFTREE)
    a.SavePDFTree();

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

