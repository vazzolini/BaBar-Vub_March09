
#include "VirVubFitter/VirClass.hh"

#include <fstream>
#include <sstream>
#include <map>


#include "VirVubFitter/VirHelper.hh"

#include "RecoilAnalysis/recoilDSys.hh"

#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2.h>

#include "TStopwatch.h"

#include "TStyle.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include <TPaveText.h>
#include <TLatex.h>
#include "TCanvas.h"
#include "TPostScript.h"

#include <TRandom3.h>

#include "RooFitCore/RooGlobalFunc.hh"
#include "RooFitCore/RooArgSet.hh"

#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooDataSet.hh"
#include "RooFitCore/RooPlot.hh"

#include "VirVubFitter/TMesCor.hh"
#include "VirVubFitter/PidCorrectMesMean.hh"

using namespace std;
using namespace RooFit;

TH1D* vubHist, *vcbHist, * othHist, *dataHist; // global histograms
TH2D* vub2Hist, *vcb2Hist, * oth2Hist, *data2Hist; // global histograms
// Multiplicity category reweighting global variables
double datmult[5],datmulterr[5],datsum;
double vubmult[5],vubsum;
double multmatrix[5][5];


// ----------------------------------------------------------------------
VirClass::VirClass() 
  : compmod(), _shapeFunction(VirClass::summer07) //old _shapeFunction(VirClass::eBelle04)
{
}

// ----------------------------------------------------------------------
VirClass::VirClass(TTree* tree,TString filename, int Sys, int q2F, int comb, int un, int unfB, double hiunfB, int mx2u,int me, int mu, bool iscm2,bool varfit,const vector<float>& wFermivec, int newbin, int SPseed, bool isWD, int rel, int bsys, bool thecmp) 
  : compmod(iscm2,varfit,wFermivec), _shapeFunction(VirClass::summer07) //old _shapeFunction def: eBelle04
{

  mySL=mylpYesSig= mynucut= mychcut= mydepl= myWdeltaM= myAllcut=0;
  
  THECOMPARISON = thecmp;
  
  isWriteData = isWD; if(THECOMPARISON) isWriteData = true;

  fnewbin = newbin;

  fOutTree = NULL;
  mesIdx = 1;

  if(TOYHISTOGRAMES) randz = new TRandom3(0);

  Init(tree);
  
  if (q2F == 0) {
    choplowB = 0.; chophighB = 5.; nB = 11;
    Vchop= new RooRealVar("chop","m_X(GeV)",-10000.,10000.);
  }else if (q2F == 1) {
    choplowB = 0.; chophighB = 26.; nB = 14;
    Vchop= new RooRealVar("chop","q^2(GeV^2)",-10000.,10000.);
  } else if (q2F == 2) {
    choplowB = 0.; chophighB = 5.28; nB = 9;
    Vchop= new RooRealVar("chop","P+(GeV)",-10000.,10000.);
  } else {
    std::cout << "Error: don't know option q2 = " << q2F << "! Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }

  if(unfB!=0){ // for the unfolding binning, KT
    cout << "Equidistant mX binning: " << unfB << " bins, endpoint " << hiunfB << endl;
    choplowB = 0.; chophighB = hiunfB; nB = unfB+1;
    Vchop = new RooRealVar("chop","mx(GeV)",-10000.,10000.);
  }

  ME = me;
  MU = mu;

  FITQ2 = q2F;
  UNFBINNING = unfB;
  Sun = un;
  UNFMX2 = mx2u;
  REL = rel;
  // needs to be called after setting main flags, ie FITQ2
  initRest(filename);

  SEME = SPseed; 

  char ddecayfile[20];
  if(iscm2){
    if(REL == 22)
      sprintf(ddecayfile,"ddecay.table.CM2.R22");
    else
      sprintf(ddecayfile,"ddecay.table.CM2.R18");
  } else {
    sprintf(ddecayfile,"ddecay.table.CM1");
  }

  correctionratiovub = correctionratiovcb = correctionratiovubout = 1.;

  TRandom random(0);
  random.SetSeed(0); 
  int therandom = 0;
  dImode = 2;
  if(Sys > 0) {
    dImode = Sys;
    therandom = (int)random.Rndm() * 1000000;
    if(Sys > 2) {
      therandom = Sys;
      Sys = 2;
    }
    if(Sys == 2) {
      if(UNFBINNING){
        Dvar = new recoilDSys(ddecayfile,0,Sys);
        Dvar->recoilDSys3(therandom);
      }
      else
	Dvar = new recoilDSys(ddecayfile,therandom,Sys);
    } else {
      if(UNFBINNING){
        Dvar = new recoilDSys("dIdecay.table",0,Sys);
        Dvar->recoilDSys3(therandom);
      }
      else 
	Dvar = new recoilDSys("dIdecay.table",therandom,Sys);
    }
    if(UNFBINNING){
      Bsem = new recoilDSys(0,REL,bsys);
      Bsem->recoilDSys2(therandom);
    }
    else
      Bsem = new recoilDSys(therandom,REL,bsys);
  } else {
    Dvar = new recoilDSys(ddecayfile,therandom,2);
    Bsem = new recoilDSys(therandom,REL,bsys);
  }
  //  cout<<"FROM CTOR: bsys "<<bsys<<endl;
  //  cout<<"therandom::  "<<therandom<<"   Dvar::   "<<Dvar<<"   Bsem::  "<<Bsem<<endl;

  Vmes      = new RooRealVar("mes","m_{ES} (GeV/c^{2})",5.22,5.3);
  VlepYes   = new RooRealVar("lepYes","lepYes",0,1);
  VlepVub   = new RooRealVar("lepVub","lepVub",0,1);
  VlepVcb   = new RooRealVar("lepVcb","lepVcb",0,1);
  VlepVubSB = new RooRealVar("lepVubSB","lepVubSB",0,1);
  VflavB    = new RooRealVar("flavB","flavB",0,5);
  VlepYaSe  = new RooRealVar("lepYaSe","lepYaSe",0,1);
  Vwe       = new RooRealVar("weight","weight",0.,100.);
  Vtrumtch  = new RooRealVar("trumtch","trumtch",0,3);

  //Unfolding
  Vallmes     = new RooRealVar("allmes","allmes(GeV)",-100.,10);
  Vmxgenwoph  = new RooRealVar("mxgenwoph","m_{x} genwoph (GeV/c^{2})",0.,25.);
  Vmultcat    = new RooRealVar("multcat","multcat",0,5);
  Vmultcatgen = new RooRealVar("multcatgen","multcatgen",0,5);
  Vvxbtyp     = new RooRealVar("vxbtyp","vxbtyp",-30,30);
  Velmom      = new RooRealVar("elmom","elmom",0.,5.);

  //Bidimensional study
  Vmx = new RooRealVar("mx","m_{x} (GeV/c^{2})",-10000.,10000.);
  Vq2 = new RooRealVar("q2","q^{2} (GeV^{2})",-10000.,10000.);

  //Other useful stuff 
  Vksele = new RooRealVar("ksele","ksele",0,1);
  Vpplus = new RooRealVar("pplus","P_{+} (GeV/c^2)",-10000.,10000.);
  Vintpur = new RooRealVar("intpur","intpur",0,1);
  Vch = new RooRealVar("isBch","brecochage",0,1);
  
  //Breco study
  VmodeB = new RooRealVar("modeB","modeB",-16000,16000);
  VtruemodeB = new RooRealVar("truemodeB","truemodeB",-16000,16000);
  Vhaspi0 = new RooRealVar("haspi0","haspi0",0,1);
  Vbrecoid = new RooRealVar("brecoid","brecoid",-3,3); 
  Vbrecoidtrue = new RooRealVar("brecoidtrue","brecoidtrue", -3, 3);

  //new truth-matching 
  Vch1B = new RooRealVar("ch1B","ch1B",0,15);
  Vneu1B = new RooRealVar("neu1B","neu1B",0,15);

  //DeltaE study
  Vde = new RooRealVar("de","#delta E (GeV)",-0.12,0.12);

  //useful stuff for the comparison
  Vpcms = new RooRealVar("pcms","pcms",0,20);
  Vmm2 = new RooRealVar("mm2","mm2",-1000,1000);
  Vnchg = new RooRealVar("nchg","nchg",0,15);
  Vnneu = new RooRealVar("nneu","nneu",0,30);
  Vnkp = new RooRealVar("nkp","nkp",0,10);
  Vnks = new RooRealVar("nks","nks",0,10);
  Vqtot = new RooRealVar("qtot","qtot",-5,10);
  Vwdeltam = new RooRealVar("wdeltam","wdeltam",-10000,100000);
  Vwdeltampiz = new RooRealVar("wdeltampiz","wdeltampiz",-10000,100000);
  Vemiss = new RooRealVar("emiss","emiss",-1000,1000);
  Vpmiss = new RooRealVar("pmiss","pmiss",-100,100);
  VlepYesBitmask = new RooRealVar("lepYesBitmask","lepYesBitmask",0,127);
  VAllCutBitmask = new RooRealVar("AllCutBitmask","AllCutBitmask",0,2047);

  // and RooDataSets

  //Standard variables in the dataset
  RooArgSet mySet = RooArgSet(*Vmes,*Vchop,*Vmx,*Vq2,*Vwe,*VlepYes,*VflavB,*VlepYaSe,*Vtrumtch);
  mySet.add(*Vmultcat); mySet.add(*Vde); mySet.add(*Vbrecoidtrue); mySet.add(*Vbrecoid);
  
  RooArgSet pstarSet = RooArgSet(*Vmes, *VlepYes, *VlepYaSe, *VflavB, *VlepVcb, *VlepVub, *VlepVubSB, *Vtrumtch, *Vwe);
  pstarSet.add(*Vde); pstarSet.add(*Vchop); pstarSet.add(*Vmx); pstarSet.add(*Vq2); pstarSet.add(*Vbrecoidtrue); pstarSet.add(*Vbrecoid);

  RooArgSet SetTheComp = RooArgSet(*Vmes,*Vksele,*VflavB,*VlepYesBitmask, *VAllCutBitmask, *Vwe);
  
  //  SetTheComp.add(*Vpcms);

  SetTheComp.add(*Vmx); SetTheComp.add(*Vq2); SetTheComp.add(*Vpplus);
  SetTheComp.add(*Vmm2); SetTheComp.add(*Vqtot);
  SetTheComp.add(*Vnchg); SetTheComp.add(*Vnneu); SetTheComp.add(*Vnkp); SetTheComp.add(*Vnks);
  SetTheComp.add(*Vwdeltam); SetTheComp.add(*Vwdeltampiz); SetTheComp.add(*Vemiss); SetTheComp.add(*Vpmiss);

  //other variables in the dataset to be written on file
  if(isWriteData && !THECOMPARISON) { 
    std::cout << "Adding variables to datasets, writing datasets" << std::endl;
    
    mySet.add(*Vintpur); mySet.add(*Vch); mySet.add(*Vpplus); 
    mySet.add(*VmodeB); mySet.add(*VtruemodeB); mySet.add(*Vhaspi0);
    mySet.add(*Vch1B); mySet.add(*Vneu1B);
    
    std::cout << "Adding more variables to pstarsample datasets" << std::endl;
    pstarSet.add(*Vintpur); pstarSet.add(*Vch); 
    pstarSet.add(*Vmultcat), pstarSet.add(*Vpplus);
  }
  
  if(THECOMPARISON)
    datadata    = new RooDataSet("DATA",   "DATA",   SetTheComp, "weight");
  else
    datadata    = new RooDataSet("DATA",   "DATA",   mySet, "weight");
  
  if (!Sun) {
    datamcvub   = new RooDataSet("Vub",    "Vub",    mySet, "weight");
    datamcoth   = new RooDataSet("Oth",    "Oth",    mySet, "weight");
    if(UNFBINNING){
      RooArgSet unfRAS(*Vallmes,*Vchop,*Vwe,*VlepYaSe,*Vmxgenwoph,*Vtrumtch,*Vmultcat,*Vmultcatgen,*Vvxbtyp);
      unfRAS.add(*Velmom);
      unfmcvub  = new RooDataSet("VubUnf", "VubUnf", unfRAS, "weight");
      RooArgSet unfTRAS(*Vwe,*Vmxgenwoph,*Vmultcatgen,*Vvxbtyp,*Velmom);
      unftmcvub = new RooDataSet("VubTru", "VubTru", unfTRAS, "weight");
    }
  } else { 
    datavcboth  = new RooDataSet("VcbOth", "VcbOth", mySet, "weight");
    datavubin   = new RooDataSet("VubIN",  "VubIN",  mySet, "weight");  
    datavubout  = new RooDataSet("VubOUT", "VubOUT", mySet, "weight");  
  }
  datamcvcb     = new RooDataSet("Vcb",    "Vcb",    mySet, "weight");

  if(THECOMPARISON)
    pstarsample   = new RooDataSet("PStar", "PStar", SetTheComp, "weight");
  else
    pstarsample   = new RooDataSet("PStar", "PStar", pstarSet, "weight");
  
  pstarsamplesum = new RooDataSet("PStarsum", "PStarsum", pstarSet, "weight");
  
  fDatasetRootFile = NULL;
}

// ----------------------------------------------------------------------
VirClass::~VirClass()
{
  // delete objects allocated by this class
  if (mesCor != 0) { delete mesCor; mesCor = 0; }
  
  compmod.~CMClass();
}

bool VirClass::IsCM2() const
{
  return compmod.cm2;
}
char* VirClass::GetEv() const
{
  return compmod.ev;
}
bool VirClass::GetVarfit() const
{
  return compmod.varfit;
}
void VirClass::SetCM2(bool choice)
{
  compmod.SetCM(choice);
}
void VirClass::SetVarfit(bool choice)
{
  compmod.SetVF(choice);
}
Int_t VirClass::Loop(int isdata, int icat, int isMC, int nres,int comb, int truthonly)
{
  //---- flags legenda:
  // isdata  0 = MC to get the shapes and the efficiencies, 1 = data 
  // icat    0 = vub component, 1 = vcb and other component, 
  // nevents = number of events
  // isMC    1 = fit on fake Data(MC), 0 = fit on data
  // truthonly 0 = normal, 1 = tree contains only some truth variables needed for the spectral unfolding

  fHistFile->cd();


  int countFB = 0;
  if (fChain == 0) return 0;

  //Variables for truth-matching.
  Int_t chgreco(0), neureco(0), chgtrue(0), neutrue(0), mnchg(0), mnneu(0), mnpi0(0), mnks(0), trmtch(0);
  Int_t lepYesBitmask, AllCutBitmask;
  
  Int_t nentries = Int_t(fChain->GetEntries());
  double mesTmp;

  int nevents = 1111111111;
  if(isdata == 0 && icat == 0) {
    nevents = nVUB;
  } else if (isdata == 0 && icat == 1) {
    nevents = nVCB;
  } else {
    if (isdata == 1 && isMC == 1) {
      if(icat == 1){
	nevents = nVCBDATA;
      } else if (icat == 0) {
	nevents = nVUBDATA;
      }
    } else if (isMC ==0){
      nevents = nDATA;
    }
  }

  std::cout << "nevents:" << nevents << std::endl;

  if( nentries > nevents) nentries = nevents;

  std::cout <<  nentries << " Entries" << std::endl;
   
  f1all = new TF1("gaussall","gaus",-5*SMEARALLSIGMA,5*SMEARALLSIGMA);
  f1bkg = new TF1("gaussbck","gaus",-5*SMEARBKGSIGMA,5*SMEARBKGSIGMA);
  TRandom rand(0);


  Int_t nbytes = 0, nb = 0;

  // DEFINE VARIALBE OUTSIDE THE LOOP
  
  Int_t haspi0,flav,flavB,type,area,ys,ch,ksele,mult,multgen;
  bool selectBType, acceptance, isLp, isLpSig, isChargeCor, isHighELept, lPYes, lPYesSig, lepVub, lepVcb, lepVubSB;
  bool WdeltaCut,WdeltaCutPiz,nuCut,AllCut,ischarged;
  Double_t wtotal, w, wfermi, wfermiG, wnre, wferminok, wssbar, wNew, wOld, mxhadTmp,  q2Tmp, pplusTmp, wexcl;

  // Add a data point, with its coordinates specified in the 'mySet' argset, to the dataset (the callee).
  // Any variables present in 'mySet' but not in the dataset (the callee) is silently ignored

  RooArgSet *mySet;
  
  if(THECOMPARISON) {
    mySet = new RooArgSet(*Vmes,*Vksele,*VflavB,*VlepYesBitmask, *VAllCutBitmask, *Vwe);
    //    mySet->add(*Vpcms); 

    
    mySet->add(*Vmx); mySet->add(*Vq2); mySet->add(*Vpplus);
    mySet->add(*Vmm2); mySet->add(*Vqtot);
    mySet->add(*Vnchg); mySet->add(*Vnneu); mySet->add(*Vnkp); mySet->add(*Vnks);
    mySet->add(*Vwdeltam); mySet->add(*Vwdeltampiz); mySet->add(*Vemiss); mySet->add(*Vpmiss);

  } else {
    mySet = new RooArgSet(*Vmes,*Vchop,*VflavB,*VlepYes,*VlepYaSe,*Vmx,*Vq2,*Vtrumtch,*Vwe);
    mySet->add(*VlepVcb); mySet->add(*VlepVub); mySet->add(*VlepVubSB); mySet->add(*Vintpur);
    mySet->add(*Vch); mySet->add(*Vpplus); mySet->add(*Vhaspi0); mySet->add(*Vch1B); mySet->add(*Vneu1B);
    mySet->add(*VmodeB); mySet->add(*VtruemodeB);  mySet->add(*Vmultcat); mySet->add(*Vde); mySet->add(*Vbrecoidtrue);
    mySet->add(*Vbrecoid);
  }

  // HERE BEGINS THE LOOP

  for (Int_t jentry=0; jentry<nentries; jentry++) {     
    Int_t ientry = LoadTree(jentry); //---- in case of a TChain, ientry is the entry number in the current file
    nb = fChain->GetEntry(jentry);  
    nbytes += nb;
    lepYesBitmask = 0x0;
    AllCutBitmask = 0x0;

    // select only neutral or charged Bs or both
    // std::cout << "BTYPE = " << BTYPE << std::endl;
    // ANTONIO, 08 jan 2009 : COMMENT THIS TO ISOLATE BTYPE VARIABLE FOR CHARGE SEPARATION
    
    //    if (!(TMath::Abs(brecocharge) == BTYPE || BTYPE == 2 || truthonly==1)) continue;
    if ((truthonly==1)) continue;
    // std::cout << "brecocharge = " << brecocharge << std::endl;

    // some definitions plus the creation of the pstarfact dataset.
    //    ANTONIO, 08 jan 2009 : COMMENT THIS TO ISOLATE BTYPE VARIABLE FOR CHARGE SEPARATION
    //    selectBType = (TMath::Abs(brecocharge) == BTYPE || BTYPE == 2);
    selectBType = true;


    //Semileptonic events selection
    acceptance = (tlab<2.37) && (tlab>0.36) && (plab>0.5);
    isLp = acceptance && selectBType  && intpur > MININTPUR;  
    if (LEPTTYPE == 0) isLp = isLp && (nel > 0);
    if (LEPTTYPE == 1) isLp = isLp && (nmu > 0);
    if (LEPTTYPE == 2) isLp = isLp && (nle > 0);
    //Signal Events Selection    
    
    isLpSig = acceptance && selectBType && intpur > MININTPUR;
    if (LEPTTYPE == 0) isLpSig = isLpSig && (nel == 1);
    if (LEPTTYPE == 1) isLpSig = isLpSig && (nmu == 1);
    if (LEPTTYPE == 2) isLpSig = isLpSig && (nle == 1);

    //---- charge correlation
    flav =  lcharge + brecoflav;

    isChargeCor = !( TMath::Abs(brecocharge) != 0 && (flav != 0) );
    isHighELept = pcms > LEPTONPCUT;

    lPYes    = isHighELept && isLp && isChargeCor;
    lPYesSig = isHighELept && isLpSig && isChargeCor;
    lepVub = lPYes && vub == 1;
    lepVcb = lPYes && vcb == 1;
    lepVubSB = false;
   
    if (comb) {
      lepVubSB = lepVub && mxhadgenwoph < MXBIN && q2Gen > Q2BIN && q2Gen < Q2HICUT;
    } else {
      if (FITQ2 == 0) {
	lepVubSB = lepVub && mxhadgenwoph < MXBIN;
      } else if (FITQ2 == 1) {
	lepVubSB = lepVub && q2Gen > Q2BIN;
      } else if (FITQ2 == 2) {
	lepVubSB = lepVub && pplusgen < PPLUSBIN;
      }
    }
  
    //---- flavor category 
    //   3 = charged B,  4 = neutral B OS,  5 = neutral B SS
    
    flavB = 5;
    if (TMath::Abs(brecocharge) == 1) flavB = 3;
    if (TMath::Abs(brecocharge) == 0 && flav == 0) flavB = 4;

    // Modify flavor category in case of charge separation
    //    if ( BTYPE != 2 ) flavB = flavB == 3 ? 3 : 4;

    type = 0;
    area = 0;

    //dataset component type
    //at this point: isdata is = 0,icat !=0  
    if (isMC) {
      if (vcb && icat != 0)              type = 1; //vcb
      if (vub && icat == 0)              type = 2; //vub
      if ((vcb + vub) == 0 && icat != 0) type = 3; //other 
      
      //put D** in other
      //CB warning we screw up nsl here...
      if(FITDSS == 1){
	 if(vcb && (TMath::Abs(Gvxbtyp)==1 || TMath::Abs(Gvxbtyp)==2) && icat!=0) type=1; //vcb without D**
	 if(vub && icat ==0) type = 2;             //vub
	 if(((vcb && !(TMath::Abs(Gvxbtyp)==1 || TMath::Abs(Gvxbtyp)==2)) || (vcb+vub) ==0) && icat!=0) type =3; //other + D**
      }
      //other = D** only
      //CB warning we screw up nsl here...
      if(FITDSS == 2){
	 if(vcb && (TMath::Abs(Gvxbtyp)==1 || TMath::Abs(Gvxbtyp)==2) && icat!=0) type=1; //vcb without D**
	 if((vcb + vub) == 0 && icat!=0) type = 1; //other together with vcb 
	 if(vub && icat ==0) type = 2;             //vub
	 if(vcb && !(TMath::Abs(Gvxbtyp)==1 || TMath::Abs(Gvxbtyp)==2) && icat!=0) type =3; // D** only
      }

      // define signal box (area==1 -> vubin, area==2 -> vubout)
      if (Sun && type == 2) { 
	if (comb) {	             // 2d mx/q2
	  if (mxhadgenwoph < MXBIN && q2Gen > Q2BIN && q2Gen < Q2HICUT) { area = 1; } else { area = 2; }
	} else {
	  if (FITQ2 == 0) {          // 1d mx 
	    if (mxhadgenwoph < MXBIN)    { area = 1; } else { area = 2; }
	  } else if (FITQ2 == 1) {   // 1d q2
	    if (q2Gen > Q2BIN)       { area = 1; } else { area = 2; }
	  } else if (FITQ2 == 2) {   // 1d pplus
	    if (pplusgen < PPLUSBIN) { area = 1; } else { area = 2; }
	  }
	}
      }
    }

    //correct mes if that is what we want
    mesTmp = mes;
    if(isdata && !isfitMC()){
      mesTmp = PidCorrectMesMean::correctMes(mes,run);
    }
 
    // --- Calculate ev-by-ev reweightings --- 

    wtotal = w = wfermi = wfermiG = wnre = wferminok = wssbar = wexcl = 1.; //Weight initialization  

    if (!isdata) {
      
      // Vcb and oth weights
      w *= getTrackingWeight();                        // trk weighting
      w *= getNeutralWeight();                         // neu weighting
      w *= getBrecoWeight(intpur);                     // breco weighting  
      if(cascade!=0 && DOCASCADEWEIGHT)
	w *= getCascadeDecWeight();                    // D cascade weighting

      if (truthonly) w = 1.; // reset reweighting for truth only processing
      
      if (icat == 1) {
	w *= getBsysweight(Gvxbtyp,vub);                                // Bdec weighting
	// ANTONIO 21 Jun 2007: PATCH TO Exclude wrong caluclated Dstarlnu FF
	if(getFFDstarlnuWeight(Gvxbtyp)>2) continue; 
	else
	  w *= getFFDstarlnuWeight(Gvxbtyp);                              // B->D*lnu FF weighting
	//	if (type == 1 || type ==3)  w *= getDsysweight(GfDpi,GfDk,GfDks,GfDpiz,GfDlep,dImode,vub);  // Ddec weighting
	if (type == 1 )  w *= getDsysweight(GfDpi,GfDk,GfDks,GfDpiz,GfDlep,dImode,vub);  // Ddec weighting
	w *= getGenericSPWeight();
      }

      // --- Start vub specific part ---


      if (vub) {
	
	// ssbar popping systematics (for vub only) 
	// Warning: this breaks the total normalisation of the b->ulv BF 
	if(DOSSBARWEIGHT > 0 && GfK > 0) { 
	  if (DOSSBARWEIGHT == 1 && TMath::Abs(Gvxbtyp) == 7) wssbar *= (1.-0.3); //  -30% for non-res
	  if (DOSSBARWEIGHT == 2 && TMath::Abs(Gvxbtyp) == 7) wssbar *= (1.+0.3); //  +30% for non-res

	  if (icat == 0) ((TH1D*)gDirectory->Get("ssbar_signal"))->Fill(wssbar, 1.);
	  if (icat == 1) ((TH1D*)gDirectory->Get("ssbar_generic"))->Fill(wssbar, 1.);
	}

	//CB IMPLEMENT NEW REWEIGHTING FROM D. FORTIN
	//NEW REWEIGHTING IS APPLIED ALSO ON WIDE RESONANT EVENTS 
	//ONLY rho, pi, eta, eta', omega ARE CONSIDERED AS RESONANT...
	//rescale BRs according to the table presented in 
	//http://www.slac.stanford.edu/BFROOT/www/Physics/Analysis/AWG/InclusiveSL/common/dominique_hybrid.html
	//CB ALLOW FOR NEW MX BINNING IN DOMS FILES (-rew 11)
	//CB	if(fprlRew==10){
	if (fprlRew==10 || fprlRew==11) {
	  if (Gvxbtyp == 11) {  // B0bar ->pi+ l- nubar + CC
	    //CB am I applying this also for generic BBbar events? 
	    wfermi *= compmod.wfermivec[0];     //wfermi *= 1.33/1.8;
	  }else if (Gvxbtyp == -11) {     // B- -> pi0 l- nubar + CC
	      wfermi *= compmod.wfermivec[1];     //wfermi *= 0.72/0.9;
	    }else if (Gvxbtyp == -12) {     // B- ->eta0 l- nubar + CC
	      wfermi *= compmod.wfermivec[2];     //wfermi *= 0.84/0.3;
	    }else if (Gvxbtyp == 13) {      // B0bar -> rho+ l- nubar + CC
	      wfermi *= compmod.wfermivec[3];     //wfermi *= 2.69/2.6;
	    }else if (Gvxbtyp == -13) {     // B- -> rho0 l- nubar + CC
	      wfermi *= compmod.wfermivec[4];     //wfermi *= 1.45/1.3;
	    }else if (Gvxbtyp == -14) {     // B- -> omega0 l- nubar + CC
	      wfermi *= compmod.wfermivec[5];     //wfermi *= 1.34/1.3;
	    }else if (Gvxbtyp == -15) {     // B- -> eta_prime l- nubar + CC
	      wfermi *= compmod.wfermivec[6];     //wfermi *= 0.84/0.6;
	    }else if (Gvxbtyp == 7) {       // B0bar -> Xu+ l- nubar + CC
	    
	      // if (icat==1) continue; //changed to fix pstarsample
	    //
	    //	    wfermi *= 17.48/13.65; actually I am taking non-res from a different datafile...
	    //                             so I have to care about nonres/res normalization 
            //compute normalization factor between resonant and old hybrid...
            //RIC says the ratio between resonant (narrow+wide) and non-resonant events in root files is 0.91
            //I assume is the same for B+- and B0 (?)
            //so I have to scale 0.91 with (BR_res^NEW / BR_res^OLD) / (BR_res^NEW / BR_nres^NEW)=BR_nres^NEW /BR_res^OLD=1.71/0.735
	    
	    //wfermi *= 2.33;
	    
	    ys = newrHistmx(mxhadgen)+rHistq2(q2Gen)*8 + rHistel(ecmsgen)*8*8;
	    if(fprlRew==11){
	      ys = newbinrHistmx(mxhadgen)+rHistq2(q2Gen)*8 + rHistel(ecmsgen)*8*8;
	    }
	    
// 	    	    cout << " B0 ys mx q2 el bins:   " << ys << "  " 
// 	    		 << newrHistmx(mxhadgen) << "  "
// 	    		 << rHistq2(q2Gen) << "  "
// 	    		 << rHistel(ecmsgen) << endl;
	    
	    wfermi *= newMatrixW[ys];
	    
	    // if we are dealing with the generic sample, we have to undo that reweighting
	    
	    if (icat==1) 
	      if(REL==14){
		wfermiG = 1./genericMatrixW[genericHistmx(mxhadgen)];
		wfermi *= wfermiG;
	      } else {
		//		wfermiG = 1./genericMatrixW[ys];
		//CB we don't have to reweight the vub part in the generics 
		if(newMatrixW[ys]!=0) wfermiG = 1./newMatrixW[ys];
		wfermi *= wfermiG;
	      }
	    
	    if (FERMIAPP) wfermi *=FermiWeight(fkplus,DELTAMB,DELTAA);

	    wferminok = wfermi;

	    if(icat==0)
	      if(this->IsCM2()) {
		  wfermi *= compmod.wfermivec[7];//wfermi *= 0.46;
	      } else {
		wfermi *= compmod.wfermivec[8];//wfermi *= 2.18;
	      }
	      
		
	  } else if(Gvxbtyp == -7) {      // B- -> Xu0 l- nubar + CC
	    // if (icat==1) continue; //changed to fix pstarsample
	    //
	    //	    wfermi *= 17.90/13.70; actually I am taking non-res from a different datafile...
	    //                             so I have to care about nonres/res normalization 
	    //so I have to scale 0.91 with (BR_res^NEW / BR_res^OLD) / (BR_res^NEW / BR_nres^NEW)=BR_nres^NEW /BR_res^OLD=1.71/0.735

	    //	    wfermi *= 2.33;


	    ys = (newrHistmx(mxhadgen)+rHistq2(q2Gen)*8 + rHistel(ecmsgen)*8*8)+512;
	    if(fprlRew==11){
	      ys = (newbinrHistmx(mxhadgen)+rHistq2(q2Gen)*8 + rHistel(ecmsgen)*8*8)+512;
	    }
	    //
	    wfermi *= newMatrixW[ys];

	    // if we are dealing with the generic sample, we have to undo that reweighting
 	    if (icat==1) 
	      if(REL==14) {
		wfermiG = 1./genericMatrixW[genericHistmx(mxhadgen)];
		wfermi *= wfermiG;
	      } else {
		//		wfermiG = 1./genericMatrixW[ys];
		//CB we don't have to reweight the vub part in the generics 
		if(newMatrixW[ys]!=0) wfermiG = 1./newMatrixW[ys];
		wfermi *= wfermiG;
	      }
	    
	    if(FERMIAPP) wfermi *=FermiWeight(fkplus,DELTAMB,DELTAA);

	    wferminok = wfermi;

	    if(icat==0) //apply reweighting only on signal samples
	      if(this->IsCM2()) {
		wfermi *= compmod.wfermivec[9];// wfermi *= 0.427;
	      } else {
		wfermi *= compmod.wfermivec[10];//wfermi *= 1.92;
	      }
	    
	  } else 
	    //reject all other vubish decays
	    continue;
	  
	} else {     //end of if(fprlRew==10 || fprlRew==11){  reweighting for vub events
	  //old implementation

	  if (TMath::Abs(Gvxbtyp)==7) {
	    //CB the following line is to avoid to take nonresonant from any sample containing Vub (e.g. MC fits...)
	    if(icat==0&&nres == 0) continue;  //questo continue mi uccide il pstarsample!
	    if(FERMIAPP) {
	      wnre *= FermiWeight(fkplus,DELTAMB,DELTAA);
	      wfermi *= wnre;
	    }
	    if(fprlRew==1) {
	      wfermi *= getTrueMxWeight(mxhadgen,Gvxbtyp); 
	      //	    cout<<"Weight :: "<<getTrueMxWeight(mxhadgen,Gvxbtyp)<<" "<<mxhadgen<<" "<<Gvxbtyp<<endl;
	    } else if(fprlRew==3) {
	      int ys = rHistmx(mxhadgen)+ rHistq2(q2Gen)*14 + rHistel(ecmsgen)*14*8;
	      //int ys = rHistmx(mxhadgen);
	      wfermi *= 1.7755*MatrixW[ys];
	      if(MULFAC != 0) {
		wfermi *= 1/MULFAC;
		//cout<<"MatrixW[ys]:: "<<MatrixW[ys]<<" "<<ys<<"  "<<MULFAC<<"            nel LOOP"<<endl;
	      }
	    } else {
	      //	      cout<<"not reweighting Vub events"<<endl;
	    }
	  } else {
	    //CB this avoids double counting resonant events when doing pure MC fits... 
	    if(icat==0 && nres!=0) continue; 
	  }

	} // if (fprlRew==10 || fprlRew==11)

	// correct ratio of SP5/SP6 signal MC to reflect ratio of Run1-3/Run4 = 1.105 in data
	wNew =  wOld = 1.;

	if(REL == 14){
	  if (icat==0 && truthonly==0) {
	    wfermi *= getSignalSPWeight();
	  }
	}

	// fill histograms for weight checking
	if (truthonly==0 && icat==0) {
	  ((TH2D*)gDirectory->Get("signal_wfermi_2d"))->Fill(Gvxbtyp, wfermi, 1.);
	  ((TH2D*)gDirectory->Get("signal_wfermiMxhadgen_2d"))->Fill(Gvxbtyp, mxhadgen, wfermi);

	  ((TH1D*)gDirectory->Get("signal_wfermi"))->Fill(wfermi);
	  ((TH1D*)gDirectory->Get("signal_Mxhadgen"))->Fill(mxhadgen);
	  ((TH1D*)gDirectory->Get("signal_wfermiMxhadgen"))->Fill(mxhadgen, wfermi);
	  if (TMath::Abs(Gvxbtyp)==7) {
	    ((TH1D*)gDirectory->Get("signal_nre_Mxhadgen"))->Fill(mxhadgen);
	    ((TH1D*)gDirectory->Get("signal_nre_wfermiMxhadgen"))->Fill(mxhadgen,wfermi);
	  } else {
	    ((TH1D*)gDirectory->Get("signal_res_Mxhadgen"))->Fill(mxhadgen);
	    ((TH1D*)gDirectory->Get("signal_res_wfermiMxhadgen"))->Fill(mxhadgen,wfermi);
	  }
	}
	if (truthonly==0 && icat==1) {
	  ((TH2D*)gDirectory->Get("generic_wfermi_2d"))->Fill(Gvxbtyp, wfermi, 1.);
	  ((TH2D*)gDirectory->Get("generic_wfermiMxhadgen_2d"))->Fill(Gvxbtyp, mxhadgen, wfermi);

	  ((TH1D*)gDirectory->Get("generic_wfermi"))->Fill(wfermi, 1.);
	  ((TH1D*)gDirectory->Get("generic_Mxhadgen"))->Fill(mxhadgen);
	  ((TH1D*)gDirectory->Get("generic_wfermiGMxhadgen"))->Fill(mxhadgen, wfermiG);
	  ((TH1D*)gDirectory->Get("generic_wfermiMxhadgen"))->Fill(mxhadgen, wfermi);
	  if (TMath::Abs(Gvxbtyp)==7) {
	    ((TH1D*)gDirectory->Get("generic_nre_Mxhadgen"))->Fill(mxhadgen);
	    ((TH1D*)gDirectory->Get("generic_nre_wfermiGMxhadgen"))->Fill(mxhadgen,wfermiG);
	    ((TH1D*)gDirectory->Get("generic_nre_wfermiMxhadgen"))->Fill(mxhadgen,wfermi);
	  } else {
	    ((TH1D*)gDirectory->Get("generic_res_Mxhadgen"))->Fill(mxhadgen);
	    ((TH1D*)gDirectory->Get("generic_res_wfermiGMxhadgen"))->Fill(mxhadgen,wfermiG);
	    ((TH1D*)gDirectory->Get("generic_res_wfermiMxhadgen"))->Fill(mxhadgen,wfermi);
	  }
	}

	if (TMath::Abs(Gvxbtyp)==7) {
	  //	  cout  << mxhadgen << "  " <<q2Gen<< "  "<<ecmsgen<< "  "<<wfermi<<endl;
	  ((TH1D*)gDirectory->Get("inclwei"))->Fill(wfermi, 1.);
	} else {
	  ((TH1D*)gDirectory->Get("exclwei"))->Fill(wfermi, 1.);
	}

      } // if (vub) 

      //implement here the FF reweighting for exclusive Vub events                                                                                                                
      if(DOEXCLFFWEIGHT){
	if (TMath::Abs(Gvxbtyp)==11||TMath::Abs(Gvxbtyp)==12||TMath::Abs(Gvxbtyp)==13||TMath::Abs(Gvxbtyp)==14||TMath::Abs(Gvxbtyp)==15) wexcl *= getFFXulnuWeight(Gvxbtyp);


	if (wexcl>10) continue;

	
	// build final weight 
	wtotal     = w * wssbar * wfermi * wexcl;
	wferminok *= w * wssbar;
	
	// plots for pilnu q2 spectrum                                                                                                                                                
	if (TMath::Abs(Gvxbtyp)==11&&icat==0) ((TH1D*)gDirectory->Get("pilnuq2"))->Fill(q2Gen,wtotal);
	if (TMath::Abs(Gvxbtyp)==11&&icat==0) {
	  if (TMath::IsNaN(wexcl)) {
	    
	  } else {
	    //	    cout << "TEST TEST pilnu q2 plot, wexcl = " << wexcl << endl;
	    ((TH1D*)gDirectory->Get("pilnuq2Ball"))->Fill(q2Gen,wexcl*wtotal);
	  }
	}
	
	// plots for rholnu q2 spectrum                                                                                                                                               
	if (TMath::Abs(Gvxbtyp)==13&&icat==0) ((TH1D*)gDirectory->Get("rholnuq2"))->Fill(q2Gen,wtotal);
	if (TMath::Abs(Gvxbtyp)==13&&icat==0) {
	  if (TMath::IsNaN(wexcl)) {
	    
	  } else {
	    //	    cout << "TEST TEST rholnu q2 plot, wexcl = " << wexcl << endl;
	    ((TH1D*)gDirectory->Get("rholnuq2Ball"))->Fill(q2Gen,wexcl*wtotal);
	  }
	}
	
	// plots for omegalnu q2 spectrum 
	if (TMath::Abs(Gvxbtyp)==14&&icat==0) ((TH1D*)gDirectory->Get("omegalnuq2"))->Fill(q2Gen,wtotal);
	if (TMath::Abs(Gvxbtyp)==14&&icat==0) {
	  if (TMath::IsNaN(wexcl)) {
	    
	  } else {
	    //	    cout << "TEST TEST omegalnu q2 plot, wexcl = " << wexcl << endl;
	    ((TH1D*)gDirectory->Get("omegalnuq2Ball"))->Fill(q2Gen,wexcl*wtotal);
	  }
	}
	
	// plots for etalnu q2 spectrum                                                                                                                                               
	if (TMath::Abs(Gvxbtyp)==12&&icat==0) ((TH1D*)gDirectory->Get("etalnuq2"))->Fill(q2Gen,wtotal);
	if (TMath::Abs(Gvxbtyp)==12&&icat==0) {
	  if (TMath::IsNaN(wexcl)) {
	    
	  } else {
	    //	    cout << "TEST TEST etalnu q2 plot, wexcl = " << wexcl << endl;
	    ((TH1D*)gDirectory->Get("etalnuq2Ball"))->Fill(q2Gen,wexcl*wtotal);
	  }
	}
	
	// plots for etaplnu q2 spectrum 

	if (TMath::Abs(Gvxbtyp)==15&&icat==0) ((TH1D*)gDirectory->Get("etaplnuq2"))->Fill(q2Gen,wtotal);
	if (TMath::Abs(Gvxbtyp)==15&&icat==0) {
	  if (TMath::IsNaN(wexcl)) {
	    
	  } else {
	    //	    cout << "TEST TEST etaplnu q2 plot, wexcl = " << wexcl << endl;
	    ((TH1D*)gDirectory->Get("etaplnuq2Ball"))->Fill(q2Gen,wexcl*wtotal);
	  }
	}
	
	
      } else {
	
	// build final weight
	wtotal     = w * wssbar * wfermi;
	wferminok *= w * wssbar;
      }
      
      if (truthonly == 0) {
	if(icat==0){
	  if(vub&&(Gvxbtyp!=7&&Gvxbtyp!=-7)&&mesTmp>5.27&&intpur>.5){
	    ((TH1D*)gDirectory->Get("plotres"))->Fill(mxhadgen,wtotal);
	    ((TH1D*)gDirectory->Get("plotall"))->Fill(mxhadgen,wtotal);
	  }
	  if(vub&&!(Gvxbtyp!=7&&Gvxbtyp!=-7)&&mesTmp>5.27&&intpur>.5){
	    ((TH1D*)gDirectory->Get("plotnres"))->Fill(mxhadgen,wtotal);
	    ((TH1D*)gDirectory->Get("plotall"))->Fill(mxhadgen,wtotal);
	    ((TH1D*)gDirectory->Get("plotallnonres"))->Fill(mxhadgen,wnre);
	  }
	}
      } // if (thruthonly == 0)

    } // if (!isdata)
    

    //CB set all weights equal to 1 if we are performing the fit on MC
    if(this->dontWeightVub()) {w=1; wfermi=1; wtotal=1; wferminok=1;}

    //      int ibch = 1;
    //      AllCut  # survivors lepYes to mass nu cut  

    ch = TMath::Abs(xcharge + brecocharge);  // total charge   
    if (truthonly) ch = -99;

    ksele = (nkp + nks) > 0;  // fit on the depleted sample?

    //   Multiplicity category (X l nu side)
    mult = 0;
    
    if(nchg == 1 && nneu>0)               mult = 1;
    if((nchg == 2 || nchg==3) && nneu==0) mult = 2;
    if((nchg == 2 || nchg==3) && nneu>0)  mult = 3;
    if(nchg > 3 && nneu==0)               mult = 4;
    if(nchg > 3 && nneu>0)                mult = 5;

    //Generated multiplicity category (X l nu side)
    multgen = 0;
    if (chgdaugen==1 && neudaugen>0)                    multgen = 1;
    if ((chgdaugen==2 || chgdaugen==3) && neudaugen==0) multgen = 2;
    if ((chgdaugen==2 || chgdaugen==3) && neudaugen>0)  multgen = 3;
    if (chgdaugen>3 && neudaugen==0)                    multgen = 4;
    if (chgdaugen>3 && neudaugen>0)                     multgen = 5;

    // partial D* reconstrucction: cut on invariant mass wdeltam
    WdeltaCut = (wdeltam > PRMM2CUT && brecocharge == 0);
    WdeltaCutPiz = (wdeltampiz > -2.0 && wdeltampiz < 30.);
    
    bool nuCut;
    if(NUCUT==1)
      nuCut = mm2 < MNUSQHIGH && mm2 > MNUSQLOW;
    else 
      nuCut = (emiss-pmiss) < EMPMHIGH && (emiss-pmiss) > EMPMLOW;
    
    // ANTONIO: ALL CUT STUDY, TO REMOVE!
    switch(SERVICE){
    case 0: AllCut = lPYesSig; break;
    case 1: AllCut = lPYesSig && nuCut; break;
    case 2: AllCut = lPYesSig && nuCut && ch < CHHIGH && ch > CHLOW; break;
    case 3: AllCut = lPYesSig && nuCut && ch < CHHIGH && ch > CHLOW && ksele == DEPL; break;
    case 4: AllCut = lPYesSig && nuCut && ch < CHHIGH && ch > CHLOW && ksele == DEPL && !(WdeltaCut) && !(WdeltaCutPiz) ; break;
    default: AllCut = lPYesSig && nuCut &&  ch < CHHIGH && ch > CHLOW &&  ksele == DEPL && !(WdeltaCut) && !(WdeltaCutPiz); break;
    // AllCut = lPYesSig && nuCut &&  ch < CHHIGH && ch > CHLOW &&  ksele == DEPL && !(WdeltaCut) && !(WdeltaCutPiz) && (pchi2xlbs >0.0001);
    
    }
    
    //calculate bitmask for thecomparison LepCut;
    //
    //|_____|___|____|______|_____|__________|____|
    //|ksele|nle|flav|intpur|Btype|acceptance|pcms|
    // 

    if ( ksele > 0) lepYesBitmask = 0x1; 
    lepYesBitmask <<= 1;

    switch(LEPTTYPE) {
    case 0: if( nel > 0 ) lepYesBitmask |= 0x1; break;
    case 1: if( nmu > 0 ) lepYesBitmask |= 0x1; break;
    case 2: if( nle > 0 ) lepYesBitmask |= 0x1; break;
    default: lepYesBitmask |= 0x0; break;
    }
    lepYesBitmask <<= 1;
    if ( isChargeCor ) lepYesBitmask |= 0x1; 
    lepYesBitmask <<= 1;
    if ( intpur > MININTPUR ) lepYesBitmask |= 0x1; 
    lepYesBitmask <<= 1;
    if ( selectBType ) lepYesBitmask |= 0x1; 
    lepYesBitmask <<= 1;
    if ( acceptance ) lepYesBitmask |= 0x1; 
    lepYesBitmask <<= 1;
    if ( isHighELept ) lepYesBitmask |= 0x1;

    /*
    if(!isHighELept){
      cout << "nel " << nel <<" nmu "<<nmu<<" nle "<<nle<<endl;
      cout << "isChargeCor "<< isChargeCor << endl;
      cout << "intpur " << (intpur > MININTPUR) <<endl;
      cout << "selectBType "<< selectBType <<endl;
      cout << "acceptance " << acceptance <<endl;
      cout << "isHighELept " << isHighELept <<endl;
      cout << "lepYes " <<lPYes<<" Bitmask "<<lepYesBitmask<<endl;
    }
    */

    //calculate bitmask for thecomparison AllCut
    //
    // |_____|____________|_________|___________|___|___|____|______|___|______|____|
    // |ksele|wdeltacutpiz|wdeltaCut|totalcharge|mm2|nle|flav|intpur|chg|accept|pcms|
    //

    if ( ksele > 0) AllCutBitmask |= 0x1; 
    AllCutBitmask <<= 1;
    if ( ! WdeltaCutPiz ) AllCutBitmask |= 0x1; 
    AllCutBitmask <<= 1;
    if ( !WdeltaCut ) AllCutBitmask |= 0x1; 
    AllCutBitmask <<= 1;
    if ( ch < CHHIGH && ch > CHLOW ) AllCutBitmask |= 0x1;
    AllCutBitmask <<= 1;
    if ( nuCut ) AllCutBitmask |= 0x1; 
    AllCutBitmask <<= 1;
    switch(LEPTTYPE) {
    case 0: if ( nel == 1 ) AllCutBitmask |= 0x1; break;
    case 1: if ( nmu == 1 ) AllCutBitmask |= 0x1; break;
    case 2: if ( nle == 1 ) AllCutBitmask |= 0x1; break;
    default: AllCutBitmask |= 0x0;
    } 
    AllCutBitmask <<= 1;
    if ( isChargeCor ) AllCutBitmask |= 1; 
    AllCutBitmask <<= 1;
    if ( intpur > MININTPUR ) AllCutBitmask |= 1; 
    AllCutBitmask <<= 1;
    if ( selectBType ) AllCutBitmask |= 1; 
    AllCutBitmask <<= 1;
    if ( acceptance ) AllCutBitmask |= 1; 
    AllCutBitmask <<= 1;
    if (isHighELept ) AllCutBitmask |= 1;     
    

    if(lPYes)  mySL++;
    if(lPYesSig) mylpYesSig++;
    if(lPYesSig && nuCut) mynucut++;
    if(lPYesSig && nuCut && ch< CHHIGH && ch > CHLOW) mychcut++;
    if(lPYesSig && nuCut && ch< CHHIGH && ch > CHLOW && ksele ==DEPL) mydepl++;
    if(AllCut) myWdeltaM++;

    
    if(this->GetVarfit()) 
      AllCut= AllCut && (q2fit>Q2CUT && mxhadfit>0. && mxhadfit<MXCUT);
    else 
      AllCut= AllCut && (q2>Q2CUT && mxhad>0. && mxhad<MXCUT);
    
    if (MU) AllCut = AllCut && mult==MU;

    if(TMath::Abs(brecocharge)==0) ischarged=false;
    if(TMath::Abs(brecocharge)==1) ischarged=true;


    if(AllCut) myAllcut++;

    // --- prepare fitted or measure kinematic variables ---

     mxhadTmp = mxhad;
     q2Tmp    = q2;
     pplusTmp = pplus;

    if(this->GetVarfit()) {
      mxhadTmp = mxhadfit;
      q2Tmp    = q2fit;
      pplusTmp = pplusfit;
    }

    // test for NaN
    if(TMath::IsNaN(mxhadTmp)) {
      std::cout << "MEZZEGA: NAN in MXHAD" << (this->GetVarfit()? "FIT" : "") << "!!" << std::endl;
      mxhadTmp = -999.;
    }   
    if(FITQ2 == 1 && TMath::IsNaN(q2Tmp)) {
      std::cout << "MEZZEGA: NAN in Q2" << (this->GetVarfit()? "FIT" : "") << "!!" << std::endl;
      q2Tmp = -999.;
    }
    if(FITQ2 == 2 && TMath::IsNaN(pplusTmp)) {
      std::cout << "MEZZEGA: NAN in PPLUS" << (this->GetVarfit()? "FIT" : "") << "!!" << std::endl;
      pplusTmp = -999.;
    }


    //---- tag the event with the lepton type
    if (mesTmp > 5.22 && lPYes) { 
    //    if (mesTmp > 5.22) { 
    
      Vmes->setVal(mesTmp);
      Vmx->setVal(mxhadTmp);
      Vq2->setVal(q2Tmp);
      Vpplus->setVal(pplusTmp);

      if (FITQ2 == 0) 
	if (UNFMX2) 
	  Vchop->setVal(mxhadTmp*mxhadTmp);
	else 
	  Vchop->setVal(mxhadTmp);
      
      else if (FITQ2 == 1) 
	Vchop->setVal(q2Tmp);
      else if (FITQ2 == 2) 
	Vchop->setVal(pplusTmp);
      
	
      //----------------------

      VlepYes->setVal(lPYes);  
      VflavB->setVal(flavB);
      VlepYaSe->setVal(AllCut);
      Vmultcat->setVal(mult);
      
      if(isdata == 1) 
	Vtrumtch->setVal(0.);
      else {
	  // int trmtch = (truemodeB == modeB) ? 2 : 1; //OLDTRUTHMATCHING

	  // NEW TRUTH-MATCHING RECONSTRUTCTED TRACKS:
	GetMultiplicity(modeB,mnchg,mnneu,mnks,mnpi0);
	chgreco=mnchg+2*mnks;
	neureco=mnneu+2*mnpi0;
	
	mnchg=mnneu=mnks=mnpi0=0;
	
	// NEW TRUTH-MATCHING GENERATED TRACKS:
	GetMultiplicity(truemodeB,mnchg,mnneu,mnks,mnpi0);
	chgtrue=mnchg+2*mnks;
	neutrue=mnneu+2*mnpi0;
	
	trmtch = (truemodeB!=-1 && TMath::Abs(chgreco-ch1B)==0 && TMath::Abs(neureco-neu1B)<3 && TMath::Abs(chgtrue-ch1B)==0 && TMath::Abs(neutrue-neu1B)<3) ? 2 : 1;
	//	  trmtch = (truemodeB!=-1 && TMath::Abs(chgreco-ch1B)==0 && TMath::Abs(neureco-neu1B)<3 && TMath::Abs(chgtrue-ch1B)==0 && TMath::Abs(neutrue-neu1B)<3) || (mesTmp > 5.27 && TMath::Abs(de)<0.019) ? 2 : 1;
	Vtrumtch->setVal(trmtch);
      }
    
      VlepVub->setVal(lepVub);  
      VlepVcb->setVal(lepVcb); 
      VlepVubSB->setVal(lepVubSB);
      Vch->setVal(ischarged);
      Vintpur->setVal(intpur);
      Vksele->setVal(ksele);
      VmodeB->setVal(modeB);
      VtruemodeB->setVal(truemodeB);
      Vch1B->setVal(ch1B);
      Vneu1B->setVal(neu1B);
      Vde->setVal(de);
      //CB allow selection of modes w/o pi0s
      haspi0 = Bmode(modeB);
      Vhaspi0->setVal(haspi0); 
      Vbrecoid->setVal(SetFlav(brecoid));
      Vbrecoidtrue->setVal(SetFlav(brecoidtrue));

      Vpcms->setVal(pcms);
      Vmm2->setVal(mm2);
      Vemiss->setVal(emiss);
      Vpmiss->setVal(pmiss);
      Vnchg->setVal(nchg);
      Vnneu->setVal(nneu);
      Vnkp->setVal(nkp);
      Vnks->setVal(nks);
      Vqtot->setVal(xcharge+brecocharge);
      Vwdeltam->setVal(wdeltam);
      Vwdeltampiz->setVal(wdeltampiz);
      VlepYesBitmask->setVal(lepYesBitmask);
      VAllCutBitmask->setVal(AllCutBitmask);
    
      if (icat && isdata==0) {
	if (vub == 1) {
	  //CONCEZIO	  pstarsample->add(mySet, wferminok);
	  //CONCEZIO	  Double_t weightVUB = wfermi*LUMI_GENERIC/LUMI_SIGNAL;
	  //CONCEZIO	  pstarsample->add(mySet, weightVUB);

	  // fill histograms
          if (TMath::Abs(Gvxbtyp) != 7) ((TH1D*)gDirectory->Get("pstar_res"))->Fill(mxhadgen,wferminok);
          if (TMath::Abs(Gvxbtyp) == 7) ((TH1D*)gDirectory->Get("pstar_nre"))->Fill(mxhadgen,wferminok);
          ((TH1D*)gDirectory->Get("pstar_all"))->Fill(mxhadgen,wferminok);
	} else 
	  //CB in principle we should reweight D** also on lepton cuts
	  //CB in practice this is a second order effect on pstar since it enters only in ratios
	  //	  double wmcvcb = (((vcb && !(TMath::Abs(Gvxbtyp)==1 || TMath::Abs(Gvxbtyp)==2)) || (vcb+vub) ==0) && icat!=0)? w*dssR: w;
	  //	  pstarsample->add(mySet, wmcvcb);
	  pstarsample->add(*mySet, w);
      }
      //CONCEZIO this should make the pstarsample equal to vcb from generic and vub from signal, properly reweighted
      if (icat==0 && isdata==0 && vub==1) {
	  Double_t weightVUB = wfermi*LUMI_GENERIC/LUMI_SIGNAL;
	  pstarsample->add(*mySet, weightVUB);
      }
    } // if (mesTmp>5.22 && lPYes)
    

    // --- prepare RooDataSets for fitting ---

    if(truthonly == 0 && mesTmp > 5.22 && lPYes) {
      //    if(truthonly == 0 && mesTmp > 5.22) {

//       Vmes->setVal(mesTmp);
//       Vmx->setVal(mxhadTmp);
//       Vq2->setVal(q2Tmp);

//       if (FITQ2 == 0) {
// 	if (UNFMX2) {
// 	  Vchop->setVal(mxhadTmp*mxhadTmp);
// 	} else {
// 	  Vchop->setVal(mxhadTmp);
// 	}
//       }	else if (FITQ2 == 1) {
// 	Vchop->setVal(q2Tmp);
//       } else if (FITQ2 == 2) {
// 	Vchop->setVal(pplusTmp);
//       }
	
//       //----------------------

//       VlepYes->setVal(lPYes);  
//       VlepVub->setVal(lepVub);  
//       VlepVcb->setVal(lepVcb);  
//       VflavB->setVal(flavB);
//       VlepYaSe->setVal(AllCut);
//       Vmultcat->setVal(mult);
      
//       Vch->setVal(ischarged);
//       Vintpur->setVal(intpur);
//       Vksele->setVal(ksele);
//       VmodeB->setVal(modeB);
//       VtruemodeB->setVal(truemodeB);
//       Vch1B->setVal(ch1B);
//       Vneu1B->setVal(neu1B);

//       //CB allow selection of modes w/o pi0s
//       haspi0 = Bmode(modeB);
//       Vhaspi0->setVal(haspi0); 

//       if(isdata ==1) {
// 	Vtrumtch->setVal(0.);
//       } else {
// 	int trmtch = (truemodeB == modeB) ? 2 : 1;
// 	Vtrumtch->setVal(trmtch);
//       }
      
//       //CB 
//       RooArgSet mySet = RooArgSet(*Vmes,*Vchop,*Vmx,*Vq2, *VlepYes,*VflavB,*VlepYaSe,*Vtrumtch);
//       mySet.add(*Vmultcat);
//       //CB add more variables if datasets are written
//       if(isWriteData) {
// 	mySet.add(*Vmultcat);
// 	mySet.add(*Vintpur); mySet.add(*Vch); mySet.add(*Vpplus); mySet.add(*Vksele);
// 	//CB
// 	mySet.add(*VmodeB);mySet.add(*VtruemodeB); mySet.add(*Vhaspi0);
// 	mySet.add(*Vch1B); mySet.add(*Vneu1B);
//       }

      int mxbn = getMxBin(mxhadTmp);
    
      if (isdata == 1) 
	datadata->add(*mySet, 1.);
      else {
	
	if (!Sun){
	  if (type == 2) datamcvub->add(*mySet,wtotal);
	  if (type == 3) {
	    //CB check that the D** reweighting bin by bin makes sense
	    if(FITDSS==2){
	      if(dssRew.size() > 0)
		w*=dssRew[mxbn];
	      else 
		w*=dssR;
	    }
	    datamcoth->add(*mySet,w);
	  }
	  
	} else {  //( Sun case)
	  // 1D or 2D unfolding 
	  ///FILL RESONANT SIGNAL ONLY IF RES OPTION IS SELECTED
	  if(type ==1 || type ==3) {
	    if(vcb && (TMath::Abs(Gvxbtyp)==1 || TMath::Abs(Gvxbtyp)==2) && icat!=0) datavcboth->add(*mySet,w); //vcb without D**
	    if(((vcb && !(TMath::Abs(Gvxbtyp)==1 || TMath::Abs(Gvxbtyp)==2)) || (vcb+vub) ==0) && icat!=0) {   //other + D**
	      //CB reweight the D** contribution
	      //Dss ratio determined from Mx 3p fit	      w*=0.102527/0.185907;
	      //CB reweight bin by bin 
	      if(dssRew.size() > 0)
		w*=dssRew[mxbn];
	      else 
		w*=dssR;
	      datavcboth->add(*mySet,w);
	    }
	  }
	  //	  if(type == 1 || type == 3) datavcboth->add(mySet,w);
	  if(area == 1) datavubin->add(*mySet,wtotal); 
	  if(area == 2) datavubout->add(*mySet,wtotal);
	  // if(nres == 0 && TMath::Abs(Gvxbtyp)==7 ){
	  // if(area == 1) datavubin->add(RooArgSet(*Vmes,*Vchop,*Vmx,*Vq2,*VlepYes,*VflavB, *VlepYaSe, *Vmultcat),wtotal);
	  // if(area == 2) datavubout->add(RooArgSet(*Vmes,*Vchop,*Vmx,*Vq2,*VlepYes,*VflavB, *VlepYaSe, *Vmultcat),wtotal);	
	  //}
	} // if (!Sun)
	
	if (type == 1){
	  //CB reweight the D** contribution
	  //rewrite this	  double wmcvcb = (((vcb && !(TMath::Abs(Gvxbtyp)==1 || TMath::Abs(Gvxbtyp)==2)) || (vcb+vub) ==0) && icat!=0)? w*dssR: w;
	  double wmcvcb = w;
	  if( ( ( (vcb && !(TMath::Abs(Gvxbtyp)==1 || TMath::Abs(Gvxbtyp)==2) ) || (vcb+vub) ==0 ) && icat!=0 ) )
	    if(dssRew.size() > 0)
	      wmcvcb*=dssRew[mxbn];
	    else 
	      wmcvcb*=dssR;
	  ((TH1F*)gDirectory->Get("totalweightvcb"))->Fill(wmcvcb);
	  datamcvcb->add(*mySet,wmcvcb);
	  //	  datamcvcb->add(mySet,w);
	}
      } // if (isdata == 1)
    
    } // if(truthonly==0 && mesTmp>5.22 && lPYes)

    // ========== prepare RooDataSets for unfolding ==========================

    if (truthonly == 0 && UNFBINNING && type==2) {
      Vallmes->setVal(mesTmp);
      if (UNFMX2) {
	Vchop->setVal(mxhadTmp*mxhadTmp);
	Vmxgenwoph->setVal(mxhadgenwoph*mxhadgenwoph);
      } else {
	Vchop->setVal(mxhadTmp);
	Vmxgenwoph->setVal(mxhadgenwoph);
      }
      VlepYaSe->setVal(AllCut);
      int trmtch=0;
      if(BRECOQUAL)
	trmtch = (brecoqual==1 || brecoqual==5 || brecoqual==7) ? 2 : 1;
      else
	trmtch = (truemodeB == modeB) ? 2 : 1;
      Vtrumtch->setVal(trmtch);
      Vmultcat->setVal(mult);
      Vmultcatgen->setVal(multgen);
      Vvxbtyp->setVal(Gvxbtyp);
      Velmom->setVal(pcmsgen);

      unfmcvub->add(RooArgSet(*Vallmes,*Vchop,*VlepYaSe,*Vmxgenwoph,*Vtrumtch, *Vmultcat, *Vmultcatgen,*Vvxbtyp,*Velmom),wtotal);
    } // if (truthonly == 0 && UNFBINNING && type==2)

    if (truthonly) {
      if (UNFMX2) {
	Vmxgenwoph->setVal(mxhadgenwoph*mxhadgenwoph);
      } else {
	Vmxgenwoph->setVal(mxhadgenwoph);
      }
      Vmultcatgen->setVal(multgen);
      Vvxbtyp->setVal(Gvxbtyp);
      Velmom->setVal(pcmsgen);

      unftmcvub->add(RooArgSet(*Vmxgenwoph,*Vmultcatgen,*Vvxbtyp,*Velmom),wtotal);
    } // if (truthonly)

  } // for (Int_t jentry=0; jentry<nentries; jentry++)
  
  delete mySet;
  return nentries;
}


void VirClass::theFit(int comb,int multc,int su, int ck)
{
  TStopwatch timer; // for timing purpose 

  if (comb == 0) {
    if (FITQ2 == 0) 
      std::cout << "Here starts the real fit code. Chisquare fit to mX dist." << std::endl;
    else if (FITQ2 == 1) 
      std::cout << "Here starts the real fit code. Chisquare fit to q2 dist." << std::endl;
    else if (FITQ2 == 2) 
      std::cout << "Here starts the real fit code. Chisquare fit to P+ dist." << std::endl;
  } else 
    std::cout << "Here starts the real fit code. Chisquare fit to mX/q2 dist." << std::endl;
  
      
  //mx q2 combined analysis [2D]
  mixingCorr(comb,multc,su); //Corrects mX or q2 histos for B0 mixing and produces mx(or q2)All.eps 


  fHistFile->cd();
  if (!(UNFBINNING)) {
    if (comb) {
      //true number of vub MC events after all cuts
      std::string name = ( su ? "vubincomb" : "vubcomb" );
      std::cout << std::endl << "Starting final fit using tag " << name << std::endl << std::endl;
      vubmcaftercuts = ((TH2D*)gDirectory->Get(name.c_str()))->Integral();

      //Performs the fit to the studied distribution
      fitWithErrors(comb,su,"comb");
    } else {
      //true number of vub MC events after all cuts
      std::string name = ( su ? "vubinchop" : "vubchop");
      std::cout << std::endl << "Starting final fit using tag " << name << std::endl << std::endl;
      vubmcaftercuts = ((TH1D*)gDirectory->Get(name.c_str()))->Integral();

      //Performs the fit to the studied distribution
      fitWithErrors(comb,su,"chop");
    }

    std::cout << std::endl << "Doing now the chi2 fit" << std::endl << std::endl;
    compChisq(comb, su);     //Computes the fit's chisquare
  }
  
  // Perform the background subtraction
  std::cout << "WWW:::" << (comb+1) << "Db:::" << comb << " " << su << std::endl;
  if (comb==0) 
    doBkgSub(multc,su);
  else 
    doBkgSub2D(su);
  

  if (!(UNFBINNING)) {
    //Time for background subtracted results plots!
    std::cout << "WWW:::" << (comb+1) << "Dl::: " << comb << " " << su << std::endl;
    if (comb==0)  
      makeBkgSubPlot(su);
    else
      makeBkgSubPlot2D(su);
    
    std::cout << "Bkg subtraction successful" << std::endl;

    //Time for Extarction of ratio of BR
    std::cout << "WWW:::" << (comb+1)<< "Df:::" << comb << " " << su << " " << ck << std::endl;
    if (comb == 0) { // 1d
      if(su)
	ruslExtr(comb, su, 1);  // there is only DFN for the 1D unfolding
      else 
	ruslExtr(comb, su, 0); 	// not unfolded always uses ck==0
    } else 	     // 2d
      ruslExtr(comb, su, ck); // not-unfolded will ignore flag ck
    

    //Time for Dumping information in human format
    resultDumping(comb, su);
  }

  // print out timer
  std::cout << "Timer for fitMes: "; timer.Print();

  return;
}


void VirClass::theFitOnly(int cmb,int multc,int su, int ck){


  char name[200];
  if (FITQ2 == 0) {
    cout<<"Here starts the real fit code. Chisquare fit to mX dist."<<endl;
  } else if (FITQ2 == 1) {
    cout<<"Here starts the real fit code. Chisquare fit to q2 dist."<<endl;
  } else if (FITQ2 == 2) {
    cout<<"Here starts the real fit code. Chisquare fit to P+ dist."<<endl;
  }

  //mx q2 combined analysis [2D]
  //  mixingCorr(cmb,multc,su); //Corrects mX or q2 histos for B0 mixing and produces mx(or q2)All.eps 
  fHistFile->cd();
  if(!(UNFBINNING)){
    if(cmb) {
      //true number of vub MC events after all cuts
      if(su){
	sprintf(name, "%s","vubincomb");
	cout<<"entrati:::"<<endl;
      } else{
	sprintf(name, "%s","vubcomb");
      }
      vubmcaftercuts = ((TH2D*)gDirectory->Get(name))->Integral();
      //Performs the fit to the studied distribution
      fitWithErrors(cmb,su,"comb");
    } else {
      //true number of vub MC events after all cuts
      if(su){
	sprintf(name, "%s","vubinchop");
	cout<<"entrati:::"<<endl;
      } else{
	sprintf(name, "%s","vubchop");
      }
      vubmcaftercuts = ((TH1D*)gDirectory->Get(name))->Integral();
      //Performs the fit to the studied distribution
      fitWithErrors(cmb,su,"chop");
    }

    cout<<"Left chi2 fit"<<endl;
    if(su){
      compChisq(cmb,1);     //Computes the fit's chisquare
    }else{
      compChisq(cmb,0);     //Computes the fit's chisquare
    }
  }
  
  //Perform the background subtraction
  if (!cmb) {
    cout<<"WWW:::1Db:::"<<cmb<<su<<endl;
    if(su){
      doBkgSub(multc,1);
    } else{
      doBkgSub(multc,0);
    }
  } else {
    if(su){
      doBkgSub2D(1);
    } else{
      doBkgSub2D(0);
    }
  }

  if(!(UNFBINNING)){
    //Time for background subtracted results plots!
    if(!cmb) { 
      cout<<"WWW:::1Dl:::"<<cmb<<su<<endl;
      if(su){
	cout<<"WWW:::1Dll:::"<<cmb<<su<<endl;
	makeBkgSubPlot(1);
      } else {
	makeBkgSubPlot(0);
      }
      cout << "Bkg subtraction successful" << endl;
    } else {
      if(su){
	makeBkgSubPlot2D(1);
      } else{
	makeBkgSubPlot2D(0);
      }
    }
    /* CB
   
    //Time for Extarction of ratio of BR
    if(!cmb) { //1d
      if(su){
	cout<<"WWW:::1Df:::"<<cmb<<su<<ck<<endl;    //entra qui con 010
	//// if(ck){
	  ruslExtr(0, 1, 1);    //DFN   
	////} else {
	////ruslExtr(0, 1, 0);  //BLL
	////}
	  } else {
	ruslExtr(0, 0, 0);
	}
    } else {	//2d
      cout<<"WWW:::2Df:::"<<cmb<<su<<ck<<endl;
	if(su){
	  if(ck){
	    ruslExtr(1, 1, 1);  //DFN
	  } else {
	    ruslExtr(1, 1, 0);  //BLL
	  }
	} else {
	  ruslExtr(1, 0, 0);   
	}
    }
    */
      //    if(su){
      //      if(ck){
      //	ruslExtr(1,1);
      //      }else{
      //	ruslExtr(1,0);
      //      }
      //   } else{
      //      ruslExtr(0,0);
      //    }     
      //Time for Dumping information in human format
      if(su){
	resultDumping(cmb,1);
      }else{
	resultDumping(cmb,0);
      }
  }
}

// ######################################################

#include "VirVubFitter/VirInit.icc"
#include "VirVubFitter/VirUtil2D.icc"
#include "VirVubFitter/VirUtil.icc"
