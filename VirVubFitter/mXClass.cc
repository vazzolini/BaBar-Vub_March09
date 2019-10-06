#include "VirVubFitter/mXClass.hh"

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

#include <TRandom2.h>

#include "RooFitCore/RooGlobalFunc.hh"
#include "RooFitCore/RooArgSet.hh"

#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooDataSet.hh"
#include "RooFitCore/RooPlot.hh"

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
mXClass::mXClass() 
  : compmod()
{}

// ----------------------------------------------------------------------
mXClass::mXClass(TTree* tree,TString filename, int Sys, int SysD, int unfB, double hiunfB, int mx2u,int mu, bool iscm2,bool varfit,const vector<float>& wFermivec, int rel, int comp) 
  : compmod(iscm2,varfit,wFermivec)
{

  Init(tree);
  
  choplowB = 0.; chophighB = 5.; nB = 11;
  Vchop= new RooRealVar("chop","m_X(GeV)",-10000.,10000.);

  if(unfB!=0){ // for the unfolding binning, KT
    choplowB = 0.; chophighB = hiunfB; nB = unfB+1;
    cout << "Equidistant mX binning: " << unfB << " bins, endpoint " << hiunfB << endl;
    Vchop = new RooRealVar("chop","mx(GeV)",-10000.,10000.);
  }

  MU = mu;

  UNFBINNING = unfB;
  UNFMX2 = mx2u;
  REL = rel;
  COMP = comp;

  // needs to be called after setting main flags
  initRest(filename);

  char ddecayfile[20];
  if(iscm2){
    if(REL==22){
      sprintf(ddecayfile,"ddecay.table.CM2.R22");    
      cout << "Reading rel 22 version of the D BFs" << endl;
    }else{
      if(REL==18){
	sprintf(ddecayfile,"ddecay.table.CM2.R18");    
	cout << "Reading rel 18 version of the D BFs" << endl;
      }else{
	sprintf(ddecayfile,"ddecay.table.CM2");
      }
    }
  }else{
    sprintf(ddecayfile,"ddecay.table.CM1");
  }

  int therandom = 0;
  dImode = 2;
  if(SysD > 200) { //inclusive reweighting of D BF
    therandom = SysD-200;
    SysD = 1;
    dImode = 1;
    if(UNFBINNING){
      Dvar = new recoilDSys("dIdecay.table",0,SysD);
      Dvar->recoilDSys3(therandom);
    } else {
      Dvar = new recoilDSys("dIdecay.table",therandom,SysD);
    }
  } else if(SysD > 0){//exclusive reweighting of D BF with smearing
    therandom = SysD;
    SysD = 2;
    if(UNFBINNING){
      //      Dvar = new recoilDSys(ddecayfile,therandom,SysD);
      Dvar = new recoilDSys(ddecayfile,0,SysD);
      Dvar->recoilDSys3(therandom);
    } else {
      Dvar = new recoilDSys(ddecayfile,therandom,SysD);
    }
  } else {// default reweighting of D BF (exclusive, no smearing)
    Dvar = new recoilDSys(ddecayfile,therandom,2);
  }

  //THIS IS THE GVXBTYP REWEIGHTING
  therandom = 0;
//   if(Sys > 0) {//smearing of b->clnu BFs
//     therandom = Sys;
//     if(UNFBINNING){
//       Bsem = new recoilDSys(0,REL);
//       Bsem->recoilDSys2(therandom);
//     } else {
//       Bsem = new recoilDSys(therandom,REL);
//     }
//   } else {//reweighting to default b->clnu BFs
//     Bsem = new recoilDSys(therandom,REL);
//   }

  //THIS IS THE GVCBTYP REWEIGHTING
  if(Sys > 0) {//smearing of b->clnu BFs requested...
    cout << "Sys>0 not suported anymore with the fine B->Xclnu BF reweighting" << endl;
    assert(false);
  } else {//reweighting to default b->clnu BFs
    Bsem = new recoilDSys(REL,true);
  }

  //  cout<<"therandom::  "<<therandom<<"   Dvar::   "<<Dvar<<"   Bsem::  "<<Bsem<<endl;

  Vmes      = new RooRealVar("mes","mes(GeV)",5.22,5.3);
  VlepYes   = new RooRealVar("lepYes","lepYes",0,1);
  VlepVub   = new RooRealVar("lepVub","lepVub",0,1);
  VlepVcb   = new RooRealVar("lepVcb","lepVcb",0,1);
  VflavB    = new RooRealVar("flavB","flavB",0,5);
  VlepYaSe  = new RooRealVar("lepYaSe","lepYaSe",0,1);
  Vwe       = new RooRealVar("weight","weight",0.,100.);
  Vtrumtch  = new RooRealVar("trumtch","trumtch",0,3);


  //Unfolding
  Vallmes     = new RooRealVar("allmes","allmes(GeV)",-100.,10);
  Vmxgenwoph  = new RooRealVar("mxgenwoph","mxgenwoph(GeV)",0.,25.);
  Vmultcat    = new RooRealVar("multcat","multcat",0,5);
  Vmultcatgen = new RooRealVar("multcatgen","multcatgen",0,5);
  Vvxbtyp     = new RooRealVar("vxbtyp","vxbtyp",-30,30);
  Velmom      = new RooRealVar("elmom","elmom",0.,5.);

  //Other useful stuff 
  Vksele = new RooRealVar("ksele","ksele",0,1);
  Vintpur = new RooRealVar("intpur","intpur",0,1);
  Vch = new RooRealVar("isBch","brecochage",0,1);
  
  //Breco study
  VmodeB = new RooRealVar("modeB","modeB",-16000,16000);
  VtruemodeB = new RooRealVar("truemodeB","truemodeB",-16000,16000);
  Vhaspi0 = new RooRealVar("haspi0","haspi0",0,1);

  //new truth-matching 
  Vch1B = new RooRealVar("ch1B","ch1B",0,15);
  Vneu1B = new RooRealVar("neu1B","neu1B",0,15);

  //DeltaE study
  Vde = new RooRealVar("de","deltaE (GeV)",-0.12,0.12);


  // and RooDataSets

  //Standard variables in the dataset
  RooArgSet mySet = RooArgSet(*Vmes,*Vchop,*Vwe,*VlepYes,*VflavB,*VlepYaSe,*Vtrumtch);
  mySet.add(*Vmultcat); mySet.add(*Vde);
  

  datadata      = new RooDataSet("DATA",   "DATA",   mySet, "weight");
  datamcvub   = new RooDataSet("Vub",    "Vub",    mySet, "weight");
  datamcoth   = new RooDataSet("Oth",    "Oth",    mySet, "weight");
  if(UNFBINNING){
    RooArgSet unfRAS(*Vallmes,*Vchop,*Vwe,*VlepYaSe,*Vmxgenwoph,*Vtrumtch,*Vmultcat,*Vmultcatgen,*Vvxbtyp);
    unfRAS.add(*Velmom);
    unfRAS.add(*VflavB);
    unfmcvub  = new RooDataSet("VubUnf", "VubUnf", unfRAS, "weight");
    RooArgSet unfTRAS(*Vwe,*Vmxgenwoph,*Vmultcatgen,*Vvxbtyp,*Velmom);
    unftmcvub = new RooDataSet("VubTru", "VubTru", unfTRAS, "weight");
  }
  datamcvcb     = new RooDataSet("Vcb",    "Vcb",    mySet, "weight");
  
  fDatasetRootFile=NULL;
}

// ----------------------------------------------------------------------
mXClass::~mXClass()
{
  compmod.~CMClass();
}

bool mXClass::IsCM2() const
{
  return compmod.cm2;
}
char* mXClass::GetEv() const
{
  return compmod.ev;
}
bool mXClass::GetVarfit() const
{
  return compmod.varfit;
}
void mXClass::SetCM2(bool choice)
{
  compmod.SetCM(choice);
}
void mXClass::SetVarfit(bool choice)
{
  compmod.SetVF(choice);
}
Int_t mXClass::Loop(int isdata, int icat, int isMC, int nres, int truthonly)
{
  //---- flags legenda:
  // isdata  0 = MC to get the shapes and the efficiencies, 1 = data 
  // icat    0 = vub component, 1 = vcb and other component, 
  // nevents = number of events
  // isMC    1 = fit on fake Data(MC), 0 = fit on data
  // truthonly 0 = normal, 1 = tree contains only some truth variables needed for the spectral unfolding

  fHistFile->cd();

  if (fChain == 0) return 0;

  //Variables for truth-matching.
  Int_t chgreco(0), neureco(0), chgtrue(0), neutrue(0), mnchg(0), mnneu(0), mnpi0(0), mnks(0), trmtch(0);
  
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
   
  Int_t nbytes = 0, nb = 0;

  // DEFINE VARIALBE OUTSIDE THE LOOP
  Int_t haspi0,flav,flavB,type,area,ys,ch,ksele,mult,multgen;
  bool selectBType, acceptance, isLp, isLpSig, isChargeCor, isHighELept, lPYes, lPYesSig, lepVub, lepVcb;
  bool WdeltaCut,nuCut,AllCut,ischarged;
  Double_t wtotal, w, wfermi, wssbar, wSP, mxhadTmp, mm2Tmp, emissTmp, pmissTmp, q2Tmp;
  //Double_t wfermiG, wnre, w wNew, wOld;

  // Add a data point, with its coordinates specified in the 'mySet' argset, to the dataset (the callee).
  // Any variables present in 'mySet' but not in the dataset (the callee) is silently ignored

  RooArgSet mySet = RooArgSet(*Vmes,*Vchop,*VflavB,*VlepYes,*VlepYaSe,*Vtrumtch,*Vwe);
  
  mySet.add(*VlepVcb); mySet.add(*VlepVub); mySet.add(*Vintpur);
  mySet.add(*Vch); mySet.add(*Vhaspi0); mySet.add(*Vch1B); mySet.add(*Vneu1B);
  mySet.add(*VmodeB); mySet.add(*VtruemodeB); mySet.add(*Vksele);
  mySet.add(*Vmultcat);  

  // HERE BEGINS THE LOOP

  int count1(0), count2(0), count3(0), count4(0), count5(0), count6(0), count7(0), count8(0), count9(0), count10(0), count11(0), count71(0), count72(0), count73(0), count74(0);

  double effevts(0.), rewevts(0.);

  //for looking at B(B->K+-lnuX)/B(B->lnuX)
  double slB0withK(0.), allslB0(0.), slBcwithK(0.), allslBc(0.), slB0withrecoK(0.), slBcwithrecoK(0.), slB0norecoK(0.), slBcnorecoK(0.);

  //for the K0 production rate correction
  TRandom2 rand(5248595);

  for (Int_t jentry=0; jentry<nentries; jentry++) {     
    
    count1++;

    Int_t ientry = LoadTree(jentry); //---- in case of a TChain, ientry is the entry number in the current file
    nb = fChain->GetEntry(jentry);  
    nbytes += nb;

    // select only neutral or charged Bs or both
    // std::cout << "BTYPE = " << BTYPE << std::endl;
    if (!(TMath::Abs(brecocharge) == BTYPE || BTYPE == 2 || truthonly==1)) continue;
    //do not take any vub events from the generic sample for default use
    if(vub && icat==1 && COMP==0) continue;

//     //the K0 production rate corrections
//     if(tnumKS>0){
//       double r = rand.Rndm();
//       if(tpupsKS[0]<0.4 && r>0.78) continue;
//       if(tpupsKS[0]>0.4 && tpupsKS[0]<1.4 && r>0.99) continue;
//       if(tpupsKS[0]>1.4 && r>0.91) continue;
//     }
//     if(tnumKL>0){
//       double r = rand.Rndm();
//       if(tpupsKL[0]<0.4 && r>0.78) continue;
//       if(tpupsKL[0]>0.4 && tpupsKL[0]<1.4 && r>0.99) continue;
//       if(tpupsKL[0]>1.4 && r>0.91) continue;
//     }

    // some definitions
    selectBType = (TMath::Abs(brecocharge) == BTYPE || BTYPE == 2);

    //Semileptonic events selection
    acceptance = (tlab<2.37) && (tlab>0.36) && (plab>0.5);
    isLp = acceptance && selectBType  && intpur>MININTPUR;  
    if (LEPTTYPE == 0) isLp = isLp && (nel > 0);
    if (LEPTTYPE == 1) isLp = isLp && (nmu > 0);
    if (LEPTTYPE == 2) isLp = isLp && (nle > 0);

    //Signal Events Selection    
    isLpSig = acceptance && selectBType && intpur>MININTPUR;
    if (LEPTTYPE == 0) isLpSig = isLpSig && (nel == 1);
    if (LEPTTYPE == 1) isLpSig = isLpSig && (nmu == 1);
    if (LEPTTYPE == 2) isLpSig = isLpSig && (nle == 1);

    //---- charge correlation
    flav =  lcharge + brecoflav;

    isChargeCor = !(TMath::Abs(brecocharge)!=0 && (flav)!=0);
    isHighELept = pcms > LEPTONPCUT;

    lPYes    = isHighELept && isLp && isChargeCor;
    lPYesSig = isHighELept && isLpSig && isChargeCor;
    lepVub = lPYes && vub==1;
    lepVcb = lPYes && vcb==1;

    //---- flavor category 
    //   3 = charged B,  4 = neutral B OS,  5 = neutral B SS
    
    flavB = 5;
    if (TMath::Abs(brecocharge)==1) flavB = 3;
    if (TMath::Abs(brecocharge)==0 && flav==0) flavB = 4;

    type = 0;

    //dataset component type
    //at this point: isdata is = 0,icat !=0  
    if (isMC) {
      if (vcb && icat!=0)            type = 1; //vcb
      if (vub && icat==0)            type = 2; //vub
      if ((vcb + vub)==0 && icat!=0) type = 3; //other 
    }

    //correct mes if that is what we want
    mesTmp = mes;
    if(DOMESMEANCORR && isdata){
      int runTmp=run;
      if(run>76440){
	runTmp=76440;
      }
      mesTmp = PidCorrectMesMean::correctMes(mes,runTmp);
    }

    //get the Gvxbtyp right...
    if(Gvxbtyp==-777) Gvxbtyp=-7;
    if(Gvxbtyp==777) Gvxbtyp=7;
 
    // --- Calculate ev-by-ev reweightings --- 

    //    wtotal = w = wfermi = wfermiG = wnre = wferminok = wssbar = wSP = 1.; //Weight initialization  
    w = wtotal = wfermi = wssbar = wSP = 1.; //Weight initialization  

    count2++;

    if (!isdata) {

      count3++;

      // Vcb and oth weights
      if (truthonly) w = 1.; // reset reweighting for truth only processing

      if (icat == 1) {
	w *= getBsysweight(Gvcbtyp,vub);                                // Bdec weighting
	w *= getFFDstarlnuWeight(Gvxbtyp);                              // B->D*lnu FF weighting
	w *= getDsysweight(GfDpi,GfDk,GfDks,GfDpiz,GfDlep,GfDgam,dImode,vub);  // Ddec weighting
	wSP = getGenericSPWeight(brecoidtrue);
	w *= wSP;
// 	if(abs(Gvxbtyp)!=1 && abs(Gvxbtyp)!=2){
// 	  if(DEPL)
// 	    w *= 0.682;
// 	  else
// 	    w *= 0.863; 
// 	}
      }

//       w *= getmm2Weight(mm2);

      //reweight decays with a KL from a D
//       if(GfDkl>0){
// 	w *= 1.237;
//       }

//      w *= getchMultWeight(nchg);
      //      w *= getneuMultWeight(nneu);
      //      w *= getchneuMultWeight(10*nchg+nneu);

//       if(nkp==0) w *= 1.021;
//       if(nkp==1) w *= 0.975328;
//       if(nkp==2) w *= 1.13379;

      count4++;

      // --- Start vub specific part ---

      if (vub && icat==0) {

	count5++;

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
	  if (Gvxbtyp == 11 ) {           // B0bar ->pi+ l- nubar + CC
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
	    
	    ys = newrHistmx(mxhadgen)+rHistq2(q2Gen)*8 + rHistel(ecmsgen)*8*8;
	    if(fprlRew==11){
	      ys = newbinrHistmx(mxhadgen)+rHistq2(q2Gen)*8 + rHistel(ecmsgen)*8*8;
	    }
	    
	    wfermi *= newMatrixW[ys];
	    
	    if(this->IsCM2()) {
	      wfermi *= compmod.wfermivec[7];//wfermi *= 0.46;
	    } else {
	      wfermi *= compmod.wfermivec[8];//wfermi *= 2.18;
	    }
	    
	  } else if(Gvxbtyp == -7) {      // B- -> Xu0 l- nubar + CC
	    ys = (newrHistmx(mxhadgen)+rHistq2(q2Gen)*8 + rHistel(ecmsgen)*8*8)+512;
	    if(fprlRew==11){
	      ys = (newbinrHistmx(mxhadgen)+rHistq2(q2Gen)*8 + rHistel(ecmsgen)*8*8)+512;
	    }
	    //
	    wfermi *= newMatrixW[ys];

	    if(this->IsCM2()) {
	      wfermi *= compmod.wfermivec[9];// wfermi *= 0.427;
	    } else {
	      wfermi *= compmod.wfermivec[10];//wfermi *= 1.92;
	    }

	  } else {
	    //reject all other vubish decays
	    continue;
	  }

	} else {     //end of if(fprlRew==10 || fprlRew==11){  reweighting for vub events
	  if(fprlRew){
	    cout << "Old fprlrew not supported anymore!" << endl;
	    break;
	  }
	} // if (fprlRew==10 || fprlRew==11)

	if (icat==0 && truthonly==0 && COMP==0) {
	  wSP = getSignalSPWeight(brecoidtrue);
	  wfermi *= wSP;
	}

	// fill histograms for weight checking
	if (truthonly==0 && icat==0) {
	  ((TH2D*)gDirectory->Get("signal_wfermi_2d"))->Fill(Gvxbtyp, wfermi, 1.);
	  ((TH2D*)gDirectory->Get("signal_wfermiMxhadgen_2d"))->Fill(Gvxbtyp, mxhadgen, wfermi);

	  ((TH1D*)gDirectory->Get("signal_wfermi"))->Fill(wfermi);
	  ((TH1D*)gDirectory->Get("signal_Mxhadgen"))->Fill(mxhadgen);
	  ((TH1D*)gDirectory->Get("signal_wfermiMxhadgen"))->Fill(mxhadgen, wfermi);
	  if (std::abs(Gvxbtyp)==7) {
	    ((TH1D*)gDirectory->Get("signal_nre_Mxhadgen"))->Fill(mxhadgen);
	    ((TH1D*)gDirectory->Get("signal_nre_wfermiMxhadgen"))->Fill(mxhadgen,wfermi);
	  } else {
	    ((TH1D*)gDirectory->Get("signal_res_Mxhadgen"))->Fill(mxhadgen);
	    ((TH1D*)gDirectory->Get("signal_res_wfermiMxhadgen"))->Fill(mxhadgen,wfermi);
	  }
	}

	if (TMath::Abs(Gvxbtyp)==7) {
	  ((TH1D*)gDirectory->Get("inclwei"))->Fill(wfermi, 1.);
	} else {
	  ((TH1D*)gDirectory->Get("exclwei"))->Fill(wfermi, 1.);
	}

	count6++;

      } // if (vub && icat==0) 

      // build final weight
      wtotal     = w * wssbar * wfermi;

      if (truthonly == 0) {
	if(icat==0){
	  if(vub&&(Gvxbtyp!=7&&Gvxbtyp!=-7)&&mesTmp>5.27&&intpur>.5){
	    ((TH1D*)gDirectory->Get("plotres"))->Fill(mxhadgen,wtotal);
	    ((TH1D*)gDirectory->Get("plotall"))->Fill(mxhadgen,wtotal);
	  }
	  if(vub&&!(Gvxbtyp!=7&&Gvxbtyp!=-7)&&mesTmp>5.27&&intpur>.5){
	    ((TH1D*)gDirectory->Get("plotnres"))->Fill(mxhadgen,wtotal);
	    ((TH1D*)gDirectory->Get("plotall"))->Fill(mxhadgen,wtotal);
	  }
	}
      } // if (thruthonly == 0)

    } // if (!isdata)
    
    ch = TMath::Abs(xcharge + brecocharge);  // total charge   
    if (truthonly) ch = -99;

    ksele = (nkp + nks) > 0;  // fit on the depleted sample?
    //ksele = nkp > 0;  // fit on the depleted sample?
       
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
    //WdeltaCut true is a good partially reco'd D*
    WdeltaCut = (wdeltam>PRMM2CUT && brecocharge == 0) || (wdeltampiz>-2 && wdeltampiz<30.);

    //charged slow pions only
    //    WdeltaCut = (wdeltam>PRMM2CUT && brecocharge == 0);

    //neutral slow pions only
    //    WdeltaCut = wdeltampiz>-2 && wdeltampiz<30.;

    //    WdeltaCut = false;

//     if(DEPL){
//       WdeltaCut = false;
//     }

     mxhadTmp = mxhad;
     q2Tmp = q2;
     mm2Tmp = mm2;
     emissTmp = emiss;
     pmissTmp = pmiss;

     if(NEUCORR){
       mxhadTmp = mxhadneucor;
       q2Tmp = q2neucor;
       mm2Tmp = mm2neucor;
       emissTmp = emissneucor;
       pmissTmp = pmissneucor;
     }

    //compute the variables from the track and neutral particles directly
    Double_t e(0.), p2(0.), px(0.), py(0.), pz(0.);
    Double_t enu(0.), pxnu(0.), pynu(0.), pznu(0.), p2nu(0.);
    Double_t excms(0.), pxcms(0.), pycms(0.), pzcms(0.);
    Double_t enucms(0.), pxnucms(0.), pynucms(0.), pznucms(0.);
    Double_t mx_list, q2_list, mm2_list, emiss_list, pmiss_list;
    Int_t numneu(0), numtrk(0);
    Int_t numtrk1(0), numtrk2(0), numtrk3(0), numtrk4(0);
    Int_t numtrk1_th(0), numtrk2_th(0), numtrk3_th(0), numtrk4_th(0);
    Double_t kaonp(0.), kaontheta(0.);

    for(int ip=0; ip<numchtrk; ip++){
      if(trktheta[ip]>0.){
	e += trke[ip];
	px += trkp[ip]*TMath::Sin(trktheta[ip])*TMath::Cos(trkphi[ip]);
	py += trkp[ip]*TMath::Sin(trktheta[ip])*TMath::Sin(trkphi[ip]);
	pz += trkp[ip]*TMath::Cos(trktheta[ip]);
	excms += trkecms[ip];
	pxcms += trkpcms[ip]*TMath::Sin(trkthetacms[ip])*TMath::Cos(trkphicms[ip]);
	pycms += trkpcms[ip]*TMath::Sin(trkthetacms[ip])*TMath::Sin(trkphicms[ip]);
	pzcms += trkpcms[ip]*TMath::Cos(trkthetacms[ip]);
	numtrk++;
      }
      if(trkp[ip]<0.2) numtrk1++;
      if(trkp[ip]>0.2 && trkp[ip]<0.7) numtrk2++;
      if(trkp[ip]>0.7 && trkp[ip]<1.3) numtrk3++;
      if(trkp[ip]>1.3) numtrk4++;
      if(trktheta[ip] < 0.7) numtrk1_th++;
      if(trktheta[ip] > 0.7 && trktheta[ip] < 1.2) numtrk2_th++;
      if(trktheta[ip] > 1.2 && trktheta[ip] < 1.9) numtrk3_th++;
      if(trktheta[ip] > 1.9) numtrk4_th++;
      if(trkK[ip]){
	kaonp = trkp[ip];
	kaontheta = trktheta[ip];
//  	w *= getpKWeight(kaonp);
      }
    }
    for(int ip=0; ip<numtotneutrk; ip++){
      if(isneumx[ip] && neue[ip]>NEUCUT && neutheta[ip]>0.0 && neutheta[ip]<5.){
	e += neue[ip];
	px += neup[ip]*TMath::Sin(neutheta[ip])*TMath::Cos(neuphi[ip]);
	py += neup[ip]*TMath::Sin(neutheta[ip])*TMath::Sin(neuphi[ip]);
	pz += neup[ip]*TMath::Cos(neutheta[ip]);
	excms += neuecms[ip];
	pxcms += neupcms[ip]*TMath::Sin(neuthetacms[ip])*TMath::Cos(neuphicms[ip]);
	pycms += neupcms[ip]*TMath::Sin(neuthetacms[ip])*TMath::Sin(neuphicms[ip]);
	pzcms += neupcms[ip]*TMath::Cos(neuthetacms[ip]);
	numneu++;
      }
    }
    p2 = px*px + py*py + pz*pz;

    if(e*e - p2){//this is consistent with what is done by HepLorentzVector
      mx_list = TMath::Sqrt(e*e - p2);
    }else{
      mx_list = TMath::Sqrt(p2 - e*e);
    }

    //compute the neutrino momentum
    enu = eUps;
    pxnu = pUps*TMath::Sin(thetaUps)*TMath::Cos(phiUps);
    pynu = pUps*TMath::Sin(thetaUps)*TMath::Sin(phiUps);
    pznu = pUps*TMath::Cos(thetaUps);
    enu -= eB;
    pxnu -= pB*TMath::Sin(thetaB)*TMath::Cos(phiB);
    pynu -= pB*TMath::Sin(thetaB)*TMath::Sin(phiB);
    pznu -= pB*TMath::Cos(thetaB);
    enu -=elab;
    pxnu -= plab*TMath::Sin(tlab)*TMath::Cos(flab);
    pynu -= plab*TMath::Sin(tlab)*TMath::Sin(flab);
    pznu -= plab*TMath::Cos(tlab);
    enu -= e;
    pxnu -= px;
    pynu -= py;
    pznu -= pz;

    enucms = esigBcms;
    enucms -=ecms;
    pxnucms -= pcms*TMath::Sin(tcms)*TMath::Cos(fcms);
    pynucms -= pcms*TMath::Sin(tcms)*TMath::Sin(fcms);
    pznucms -= pcms*TMath::Cos(tcms);
    enucms -= excms;
    pxnucms -= pxcms;
    pynucms -= pycms;
    pznucms -= pzcms;

    p2nu = pxnu*pxnu + pynu*pynu + pznu*pznu;

    mm2_list = enu*enu - p2nu;

    emiss_list = enucms;
    pmiss_list = TMath::Sqrt(pxnucms*pxnucms + pynucms*pynucms + pznucms*pznucms);

//     cout << "Var testing " << endl;
//     cout << "mx " << mxhadTmp << "  " << mx_list << "  " << mxhad << "  " << mxhadneucor<< endl;
//     cout << "mm2 " << mm2Tmp << "  " << mm2_list << "  " << mm2 << "  " << mm2neucor << endl;
//     cout << "nchg " << nchg << "  " << numtrk << endl;
//     cout << "nneu " << nneu << "  " << numneu << endl;
//     cout << "exhad " << exhad << "  " << e << endl;
//     cout << "pxhad " << pxhad << "  " << TMath::Sqrt(p2) << endl;
//     cout << "nkp " << nkp << endl;
//     cout << "nks " << nks << endl;
//     cout << "nkl " << KLnum << endl;
//     if(KLnum>0) cout << KLprodratecorr[0] << endl;
//     cout << "In the B frame " << endl;
//     cout << "emiss " << emiss << "  " << emiss_list << "  " << emiss_list-emiss << " (" << esigBcms << ")" << endl;
//     cout << "pmiss " << pmiss << "  " << pmiss_list << endl;

    if(RECOMPUTEVARS){
      mxhadTmp = mx_list;
      mm2Tmp = mm2_list;
      emissTmp = emiss_list;
      pmissTmp = pmiss_list;
    }

    bool nuCut;
    if(NUCUT==1){
      nuCut = mm2Tmp < MNUSQHIGH && mm2Tmp > MNUSQLOW;
    } else if(NUCUT==2){
      nuCut = (emissTmp-pmissTmp) < EMPMHIGH && (emissTmp-pmissTmp) > EMPMLOW;
    } else {
      nuCut = (emissTmp-pmissTmp) < EMPMHIGH && (emissTmp-pmissTmp) > EMPMLOW && mm2Tmp < MNUSQHIGH && mm2Tmp > MNUSQLOW;
    }

    AllCut = lPYesSig && nuCut &&  ch < CHHIGH && ch > CHLOW &&  ksele == DEPL && !(WdeltaCut);
    if(DEPL){
      AllCut = lPYesSig && nuCut &&  ch < CHHIGH && ch > CHLOW &&  (ksele == DEPL  || (WdeltaCut));
    }

    //The cases where we need extra cuts for the non-mx vars
    if(UNFFIT==4 || UNFFIT==5 || UNFFIT==12 || UNFFIT==13 || UNFFIT==14 || UNFFIT==15 || UNFFIT==16){//nkp, nks, nks_T, nks_VT, nks_TVT, KSdecaylenSig, KSmass, KScosp
      AllCut = lPYesSig && nuCut &&  ch < CHHIGH && ch > CHLOW;
      if(DEPL){
	AllCut = lPYesSig && nuCut &&  ch < CHHIGH && ch > CHLOW;
      }
    }
    if(UNFFIT==52 || UNFFIT==53){//kaonp, kaontheta -> require a charged K
      AllCut = lPYesSig && nuCut &&  ch < CHHIGH && ch > CHLOW && nkp;
      if(DEPL){
	AllCut = lPYesSig && nuCut &&  ch < CHHIGH && ch > CHLOW && nkp;
      }
    }

    if(UNFFIT==8){//qtot
      AllCut = lPYesSig && nuCut && ksele == DEPL && !(WdeltaCut);
      if(DEPL){
	AllCut = lPYesSig && nuCut && (ksele == DEPL  || (WdeltaCut));
      }
    }
    if(UNFFIT==9 || UNFFIT==10 || UNFFIT==18 || UNFFIT==19 || UNFFIT==20 || UNFFIT==21 || UNFFIT==22 || UNFFIT==23 || UNFFIT==24 || UNFFIT==25 || UNFFIT==26 || UNFFIT==27 || UNFFIT==28 || UNFFIT==29 || UNFFIT==30 || UNFFIT==31 || UNFFIT==32 || UNFFIT==33 || (UNFFIT>33 && UNFFIT<52)){//emiss-pmiss, mm2, f=mm2/q2 and this we want without neutrino cuts for this study
      AllCut = lPYesSig &&  ch < CHHIGH && ch > CHLOW &&  ksele == DEPL && !(WdeltaCut);
      if(DEPL){
	AllCut = lPYesSig && ch < CHHIGH && ch > CHLOW &&  (ksele == DEPL  || (WdeltaCut));
      }
    }
    if(UNFFIT==11 || UNFFIT==54){//wdeltam, wdeltampiz
      AllCut = lPYesSig && nuCut &&  ch < CHHIGH && ch > CHLOW &&  ksele == DEPL;
      if(DEPL){
	AllCut = lPYesSig && nuCut &&  ch < CHHIGH && ch > CHLOW && ksele == DEPL;
      }
    }

    if(this->GetVarfit()) 
      AllCut= AllCut && (q2fit>Q2CUT && mxhadfit>0. && mxhadfit<MXCUT);
    else 
      AllCut= AllCut && (q2Tmp>Q2CUT && mxhadTmp>0. && mxhadTmp<MXCUT);

    if (MU) AllCut = AllCut && mult==MU;

    if(TMath::Abs(brecocharge)==0) ischarged=false;
    if(TMath::Abs(brecocharge)==1) ischarged=true;

    if(vub) 
      count7++;

    if(TMath::Abs(Gvxbtyp)<900 && icat==1){
      if(Gvxbtyp>0){//neutral B
	if(GfDk){
	  slB0withK += w;
	  if(AllCut && nkp) slB0withrecoK += w;
	  if(AllCut && nkp==0) slB0norecoK += w;
	}
	allslB0 += w;
      }else{
	if(GfDk){
	  slBcwithK += w;
	  if(AllCut && nkp) slBcwithrecoK += w;
	  if(AllCut && nkp==0) slBcnorecoK += w;
	}
	allslBc += w;
      }
    }


    // --- prepare fitted or measure kinematic variables ---

    if(this->GetVarfit()) {
      mxhadTmp = mxhadfit;
    }

    // test for NaN
    if(TMath::IsNaN(mxhadTmp)) {
      std::cout << "MEZZEGA: NAN in MXHAD" << (this->GetVarfit()? "FIT" : "") << "!!" << std::endl;
      mxhadTmp = -999.;
    }   

    int multreco(0), chmultreco(0);

    // NEW TRUTH-MATCHING RECONSTRUCTED TRACKS:
    GetMultiplicity(modeB,mnchg,mnneu,mnks,mnpi0);
    chmultreco=mnchg+2*mnks;
    multreco = chmultreco + mnneu+2*mnpi0;
    mnchg=mnneu=mnks=mnpi0=0;

    //this is a look over the number of KS in the event for some of the data-MC comparison variables.
    int max=1;
    if(UNFFIT==15 || UNFFIT==16 || UNFFIT==17) max = nks;
    if(UNFFIT==19 || UNFFIT==51) max = numtotneutrk;
    if(UNFFIT==20 || UNFFIT==29 || UNFFIT==30 || UNFFIT==31 || UNFFIT==32 || UNFFIT==33) max = numchtrk;

    if(vub) 
      count71++;

    for(int nfill=1; nfill<=max; nfill++){
 
    if(vub) 
      count72++;

     //---- tag the event with the lepton type
      if (mesTmp>5.22 && lPYes) {

	if(vub) 
	  count73++;

	Vmes->setVal(mesTmp);

	if(!UNFMX2 && UNFFIT==0){
	  Vchop->setVal(mxhadTmp);
	}
	if (UNFMX2 || UNFFIT==1){ 
	  Vchop->setVal(mxhadTmp*mxhadTmp);
	}
	if(UNFFIT==2){
	  Vchop->setVal((double)numtrk);
	}
	if(UNFFIT==3){
	  Vchop->setVal((double)numneu);
	}
	if(UNFFIT==4){
	  Vchop->setVal((double)nkp);
	}
	if(UNFFIT==5){
	  Vchop->setVal((double)nks);
	}
	if(UNFFIT==6){
	  Vchop->setVal(exhad);
	}
	if(UNFFIT==7){
	  Vchop->setVal(plab);
	}
	if(UNFFIT==8){
	  Vchop->setVal((double)(xcharge + brecocharge)); //total charge
	}
	if(UNFFIT==9){
	  Vchop->setVal(emissTmp-pmissTmp);
	}
	if(UNFFIT==10){
	  Vchop->setVal(mm2Tmp);
	}
	if(UNFFIT==11){
	  Vchop->setVal(wdeltam);
	}
	if(UNFFIT==12){
	  Vchop->setVal((double)nks_T);
	}
	if(UNFFIT==13){
	  Vchop->setVal((double)nks_VT);
	}
	if(UNFFIT==14){
	  Vchop->setVal((double)nks_TVT);
	}
	if(UNFFIT==15){
	  Vchop->setVal(KSdecaylenSig[nfill-1]);
	}
	if(UNFFIT==16){
	  Vchop->setVal(KSmass[nfill-1]);
	}
 	if(UNFFIT==17){
	  Vchop->setVal(KScosp[nfill-1]);
	}
 	if(UNFFIT==18){
	  Vchop->setVal(mm2Tmp/q2);
	}
	if(UNFFIT==19){
	  if(isneumx[nfill-1] && neue[nfill-1]>NEUCUT && neutheta[nfill-1]>0. && neutheta[nfill-1]<5.){
	    Vchop->setVal(neue[nfill-1]);
	  }else{
	    continue;
	  }
	}
	if(UNFFIT==20){
	  Vchop->setVal(trkp[nfill-1]);
	}
	if(UNFFIT==21){
	  Vchop->setVal((double)numtrk1);
	}
	if(UNFFIT==22){
	  Vchop->setVal((double)numtrk2);
	}
	if(UNFFIT==23){
	  Vchop->setVal((double)numtrk3);
	}
	if(UNFFIT==24){
	  Vchop->setVal((double)numtrk4);
	}
	if(UNFFIT==25){
	  Vchop->setVal((double)numtrk1_th);
	}
	if(UNFFIT==26){
	  Vchop->setVal((double)numtrk2_th);
	}
	if(UNFFIT==27){
	  Vchop->setVal((double)numtrk3_th);
	}
	if(UNFFIT==28){
	  Vchop->setVal((double)numtrk4_th);
	}
	if(UNFFIT==29){
	  if(trktheta[nfill-1]<0.7){
	    Vchop->setVal(trkp[nfill-1]);
	  }else{
	    continue;
	  }
	}
	if(UNFFIT==30){
	  if(trktheta[nfill-1]>0.7 && trktheta[nfill-1]<1.2){
	    Vchop->setVal(trkp[nfill-1]);
	  }else{
	    continue;
	  }
	}
	if(UNFFIT==31){
	  if(trktheta[nfill-1]>1.2 && trktheta[nfill-1]<1.9){
	    Vchop->setVal(trkp[nfill-1]);
	  }else{
	    continue;
	  }
	}
	if(UNFFIT==32){
	  if(trktheta[nfill-1]>1.9){
	    Vchop->setVal(trkp[nfill-1]);
	  }else{
	    continue;
	  }
	}
	if(UNFFIT==33){
	  Vchop->setVal(trktheta[nfill-1]);
	}
	if(UNFFIT==34){
	  Vchop->setVal(multreco);
	}
	if(UNFFIT==35){
	  Vchop->setVal(chmultreco);
	}
	if(UNFFIT==36){
	  if(multreco<4){
	    Vchop->setVal(nchg);
	  }else{
	    continue;
	  }
	}
	if(UNFFIT==37){
	  if(multreco>3 && multreco<7){
	    Vchop->setVal(nchg);
	  }else{
	    continue;
	  }
	}
	if(UNFFIT==38){
	  if(multreco>6 && multreco<10){
	    Vchop->setVal(nchg);
	  }else{
	    continue;
	  }
	}
	if(UNFFIT==39){
	  if(multreco>9){
	    Vchop->setVal(nchg);
	  }else{
	    continue;
	  }
	}
	if(UNFFIT==40){
	  if(chmultreco<3){
	    Vchop->setVal(nchg);
	  }else{
	    continue;
	  }
	}
	if(UNFFIT==41){
	  if(chmultreco>2 && chmultreco<5){
	    Vchop->setVal(nchg);
	  }else{
	    continue;
	  }
	}
	if(UNFFIT==42){
	  if(chmultreco>4){
	    Vchop->setVal(nchg);
	  }else{
	    continue;
	  }
	}
	if(UNFFIT==43){
	  if(multreco<4){
	    Vchop->setVal(mm2_list);
	  }else{
	    continue;
	  }
	}
	if(UNFFIT==44){
	  if(multreco>3 && multreco<7){
	    Vchop->setVal(mm2_list);
	  }else{
	    continue;
	  }
	}
	if(UNFFIT==45){
	  if(multreco>6 && multreco<10){
	    Vchop->setVal(mm2_list);
	  }else{
	    continue;
	  }
	}
	if(UNFFIT==46){
	  if(multreco>9){
	    Vchop->setVal(mm2_list);
	  }else{
	    continue;
	  }
	}
	if(UNFFIT==47){
	  if(chmultreco<3){
	    Vchop->setVal(mm2_list);
	  }else{
	    continue;
	  }
	}
	if(UNFFIT==48){
	  if(chmultreco>2 && chmultreco<5){
	    Vchop->setVal(mm2_list);
	  }else{
	    continue;
	  }
	}
	if(UNFFIT==49){
	  if(chmultreco>4){
	    Vchop->setVal(mm2_list);
	  }else{
	    continue;
	  }
	}
	if(UNFFIT==50){
	  Vchop->setVal(10*nchg+nneu);
	}
	if(UNFFIT==51){
	  if(isneumx[nfill-1] && neue[nfill-1]>NEUCUT && neutheta[nfill-1]>0. && neutheta[nfill-1]<5.){
	    Vchop->setVal(neutheta[nfill-1]);
	  }else{
	    continue;
	  }
	}
	if(UNFFIT==52){
	  Vchop->setVal(kaonp);
	}
	if(UNFFIT==53){
	  Vchop->setVal(kaontheta);
	}
	if(UNFFIT==54){
	  Vchop->setVal(wdeltampiz);
	}
      
	//----------------------

	VlepYes->setVal(lPYes);  
	VflavB->setVal(flavB);
	VlepYaSe->setVal(AllCut);
	Vmultcat->setVal(mult);

	if(vub) 
	  count74++;


	if(isdata ==1) 
	  Vtrumtch->setVal(0.);
	else {
	  // NEW TRUTH-MATCHING RECONSTRUTCTED TRACKS:
	  GetMultiplicity(modeB,mnchg,mnneu,mnks,mnpi0);
	  chgreco=mnchg+2*mnks;
	  neureco=mnneu+2*mnpi0;
	  mnchg=mnneu=mnks=mnpi0=0;
	  // NEW TRUTH-MATCHING GENERATED TRACKS:
	  GetMultiplicity(truemodeB,mnchg,mnneu,mnks,mnpi0);
	  chgtrue=mnchg+2*mnks;
	  neutrue=mnneu+2*mnpi0;
	  trmtch = ((truemodeB!=-1 && chgreco-ch1B==0 && TMath::Abs(neureco-neu1B)<3 && chgtrue-ch1B==0 && TMath::Abs(neutrue-neu1B)<3) || ((mes>5.27 && TMath::Abs(de)<0.020))) ? 2 : 1;
	  Vtrumtch->setVal(trmtch);
	}
	
	VlepVub->setVal(lepVub);  
	VlepVcb->setVal(lepVcb); 
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

	if(vub)
	  count8++;
	
      } // if (mesTmp>5.22 && lPYes)
    
      // --- prepare RooDataSets for fitting ---

      if(truthonly==0 && mesTmp>5.22 && lPYes) {

	int mxbn = getMxBin(mxhadTmp);
      
	if (isdata == 1){
	  datadata->add(mySet, 1.);
	} else {

	  if(vub)
	    count9++;

	  if ((COMP == 1 && vub) || (COMP == 0 && type == 2)) {
	    if(COMP){//no reweighting for the vub events in the generic MC except for lumi reweighting
	      wtotal=wSP;
	    }
            datamcvub->add(mySet,wtotal);
	    effevts += wSP;
	    rewevts += wtotal;
	  }


// 	  if (type == 3 || (type == 1 && abs(Gvcbtyp) != 5 && abs(Gvxbtyp) != 6)) {
// 	    if(COMP){
// 	      w=wSP;
// 	    }
// 	    datamcvcb->add(mySet,w);
// 	  }
	  
// 	  if (type == 1 && (abs(Gvcbtyp) == 5 || abs(Gvxbtyp) == 6)){
// 	    datamcoth->add(mySet,w);
// 	    count10++;
// 	    if(AllCut){
// 	      effevts += wSP;
// 	      rewevts += w;
// 	    }
// 	  }

	  if (type == 3) {
	    if(COMP){
	      w=wSP;
	    }
	    datamcoth->add(mySet,w);
	  }
	  
	  if (type == 1){
	    datamcvcb->add(mySet,w);
	    count10++;
	    if(AllCut){
	      effevts += wSP;
	      rewevts += w;
	    }
	  }
	  
	}// if (isdata == 1)
    
      } // if(truthonly==0 && mesTmp>5.22 && lPYes)

      // ========== prepare RooDataSets for unfolding ==========================
      if(UNFFIT<2){
	if (truthonly == 0 && UNFBINNING && type==2) {
	  Vallmes->setVal(mesTmp);
	  if (UNFMX2 || UNFFIT==1) {
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
	  else{
	    // NEW TRUTH-MATCHING RECONSTRUCTED TRACKS:
	    GetMultiplicity(modeB,mnchg,mnneu,mnks,mnpi0);
	    chgreco=mnchg+2*mnks;
	    neureco=mnneu+2*mnpi0;
	    mnchg=mnneu=mnks=mnpi0=0;
	    
	    // NEW TRUTH-MATCHING GENERATED TRACKS:
	    GetMultiplicity(truemodeB,mnchg,mnneu,mnks,mnpi0);
	    chgtrue=mnchg+2*mnks;
	    neutrue=mnneu+2*mnpi0;
	    trmtch = ((truemodeB!=-1 && chgreco-ch1B==0 && TMath::Abs(neureco-neu1B)<3 && chgtrue-ch1B==0 && TMath::Abs(neutrue-neu1B)<3) || ((mes>5.27 && TMath::Abs(de)<0.020))) ? 2 : 1;
	  }
	  Vtrumtch->setVal(trmtch);
	  Vmultcat->setVal(mult);
	  Vmultcatgen->setVal(multgen);
	  Vvxbtyp->setVal(Gvxbtyp);
	  Velmom->setVal(pcmsgenwph);
	  
	  RooArgSet unfRAS(*Vallmes,*Vchop,*VlepYaSe,*Vmxgenwoph,*Vtrumtch,*Vmultcat,*Vmultcatgen,*Vvxbtyp);
	  unfRAS.add(*Velmom);
	  unfRAS.add(*VflavB);
	  
	  unfmcvub->add(unfRAS,wtotal);
	} // if (truthonly == 0 && UNFBINNING && type==2)
	
	if (truthonly && vub) {
	  if (UNFMX2 || UNFFIT==1) {
	    Vmxgenwoph->setVal(mxhadgenwoph*mxhadgenwoph);
	  } else {
	    Vmxgenwoph->setVal(mxhadgenwoph);
	  }
	  Vmultcatgen->setVal(multgen);
	  Vvxbtyp->setVal(Gvxbtyp);
	  Velmom->setVal(pcmsgenwph);
	  if(Gvxbtyp!=7){//extra magic k for truth MC due to size of truth MC samples
	    wtotal *= 1.035;
	  }
	  unftmcvub->add(RooArgSet(*Vmxgenwoph,*Vmultcatgen,*Vvxbtyp,*Velmom),wtotal);
	} // if (truthonly)
      } // if(UNFFIT<2)
    } // loop over KS (if asked for)
  } // for (Int_t jentry=0; jentry<nentries; jentry++)

  cout << "Corrections for the C* computations" << endl;
  cout << "Effective number of evts " << effevts << endl;
  cout << "Reweighted number of evts " << rewevts << endl;
  cout << "Correction factor " << rewevts/effevts << endl;

  cout << "Poor man's debugging counts" << endl;
  cout << "count1 " << count1 << endl;
  cout << "count2 " << count2 << endl;
  cout << "count3 " << count3 << endl;
  cout << "count4 " << count4 << endl;
  cout << "count5 " << count5 << endl;
  cout << "count6 " << count6 << endl;
  cout << "count7 " << count7 << endl;
  cout << "count71 " << count71 << endl;
  cout << "count72 " << count72 << endl;
  cout << "count73 " << count73 << endl;
  cout << "count74 " << count74 << endl;
  cout << "count8 " << count8 << endl;
  cout << "count9 " << count9 << endl;
  cout << "count10 " << count10 << endl;
  cout << endl;

  cout << "Looking at K+- production in sl B decays" << endl;
  cout << "sl B0               " << allslB0 << endl;
  cout << "sl B0 with K        " << slB0withK << endl;
  cout << "sl B0 with reco K   " << slB0withrecoK << endl;
  cout << "sl B0 no reco K     " << slB0norecoK << endl;
  cout << "sl B+-              " << allslBc << endl;
  cout << "sl B+- with K       " << slBcwithK << endl;
  cout << "sl B+- with reco K  " << slBcwithrecoK << endl;
  cout << "sl B+- no reco K    " << slBcnorecoK << endl;
  return nentries;
  
}

void mXClass::theFit(int multc)
{
  TStopwatch timer; // for timing purpose 

  std::cout << "Here starts the real fit code. Chisquare fit to mX dist." << std::endl;

  mixingCorr(multc); //Corrects mX or q2 histos for B0 mixing and produces mx(or q2)All.eps 

  fHistFile->cd();
  if (!(UNFBINNING)) {
    //true number of vub MC events after all cuts
    std::string name = "vubchop";
    std::cout << std::endl << "Starting final fit using tag " << name << std::endl << std::endl;
    vubmcaftercuts = ((TH1D*)gDirectory->Get(name.c_str()))->Integral();
    
    //Performs the fit to the studied distribution
    fitWithErrors("chop");
    
    std::cout << std::endl << "Doing now the chi2 fit" << std::endl << std::endl;
    compChisq();     //Computes the fit's chisquare
  }
  
  doBkgSub(multc);

  if (!(UNFBINNING)) {
    //Time for background subtracted results plots!
    makeBkgSubPlot();

    std::cout << "Bkg subtraction successful" << std::endl;

    //Time for Dumping information in human format
    resultDumping();
  }

  // print out timer
  std::cout << "Timer for fitMes: "; timer.Print();

  return;
}


// ######################################################

#include "VirVubFitter/mXInit.icc"
#include "VirVubFitter/mXUtil.icc"
