#include <fstream>

#include "VirVubFitter/exclfitNtp.hh"
#include "RecoilAnalysis/recoilDSys.hh"
#include "RecoilAnalysis/recoilBuSys.hh"

#include "TH2.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLegendEntry.h"
#include "TLegend.h"
#include <TRandom.h>
#include <TVector2.h>
#include <TLatex.h>
#include <TMinuit.h>

#include "TPostScript.h"

#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooDataSet.hh"
#include "RooFitCore/RooPlot.hh"

using namespace std;
using namespace RooFit;
// ----------------------------------------------------------------------
exclfitNtp::exclfitNtp(TTree *tree, int Sys, TString filename) {
  Init(tree);
  initRest(filename); 
  TRandom random(0);
  random.SetSeed(0);
  int therandom = 0;
  dImode = 0;
  if(Sys > 0) {
    cout<<"Passed here "<< Sys <<endl;
    
    dImode = Sys;
    therandom = int(random.Rndm() * 1000000);
    if(Sys == 2) {
      Dvar = new recoilDSys("ddecay.table",therandom,Sys);
    } else {
      Dvar = new recoilDSys("dIdecay.table",therandom,Sys);
    }
    Bsem = new recoilDSys(therandom);
    Busem = new recoilBuSys("budecay.dat",therandom);
  } else {
    cout<<"Default reweighting "<< Sys <<endl;
    Dvar = new recoilDSys("ddecay.table",therandom,2);
    Bsem = new recoilDSys(therandom);
    Busem = new recoilBuSys("budecay.dat",therandom);
  }

  varMes     = new RooRealVar("mes","mes(GeV)",5.2,5.3);
  varVar    = new RooRealVar("var","var",-10000.,100000.);
  varMm2     = new RooRealVar("mm2","mm2(GeV^{2})",-10000.,10000.);
  varLepYes  = new RooRealVar("lepYes","lepYes",0,1);
  varFlavB   = new RooRealVar("flavB","flavB",-5,5);
  varAllnovar = new RooRealVar("allnovar","allnovar",0,1);
  varAllnomm2 = new RooRealVar("allnomm2","allnomm2",0,1);
  varAllcuts = new RooRealVar("allcuts","allcuts",0,1);
  varWe      = new RooRealVar("weight","weight",0.,100.);
  datadata     = new RooDataSet("","",RooArgSet(*varMes,*varVar,*varMm2,*varWe,*varLepYes,*varFlavB,*varAllnomm2,*varAllnovar,*varAllcuts),"weight");
  datamcsig    = new RooDataSet("","",RooArgSet(*varMes,*varVar,*varMm2,*varWe,*varLepYes,*varFlavB,*varAllnomm2,*varAllnovar,*varAllcuts),"weight");
  datamcvcb    = new RooDataSet("","",RooArgSet(*varMes,*varVar,*varMm2,*varWe,*varLepYes,*varFlavB,*varAllnomm2,*varAllnovar,*varAllcuts),"weight");
  datamcvub    = new RooDataSet("","",RooArgSet(*varMes,*varVar,*varMm2,*varWe,*varLepYes,*varFlavB,*varAllnomm2,*varAllnovar,*varAllcuts),"weight");
  datamcoth    = new RooDataSet("","",RooArgSet(*varMes,*varVar,*varMm2,*varWe,*varLepYes,*varFlavB,*varAllnomm2,*varAllnovar,*varAllcuts),"weight");
  datamc       = new RooDataSet("","",RooArgSet(*varMes,*varWe,*varLepYes,*varFlavB),"weight");
  datamcvcbvub = new RooDataSet("","",RooArgSet(*varMes,*varWe,*varLepYes,*varFlavB),"weight");
  datamcsigforpstar = new RooDataSet("","",RooArgSet(*varMes,*varWe,*varLepYes,*varFlavB),"weight");
  vcbmeanweight = vubmeanweight = othmeanweight = 0;
  for (int j=0; j<9; j++) vcbmeanw[j] = vubmeanw[j] = othmeanw[j] = 0;
  countvcb = 0;
  countvub = 0;
  countoth = 0; 
  for (int j=0; j<9; j++) covcb[j] = covub[j] = cooth[j] = 0;
  TFile dal("dalitz.root");
  mydalitz = (*(TH2D*)gDirectory->Get("truedalitzsigforcut"));  
}

// ----------------------------------------------------------------------
exclfitNtp::~exclfitNtp() {
}

void exclfitNtp::Loop(int isdata, int icat, int nevents, int isMC, int nres)
  
  // flags legenda:
  // isdata  0 = MC to get the shapes and the efficiencies, 1 = fill the histos to fit
  // icat    0 = vub component, 1 = vcb and other component, 
  // nevents = number of events
  // isMC    0 = fit on the MC, 1 = fit on data
{
  //int acutcounter, cutcounter;
  //acutcounter=cutcounter=0;
  fHistFile->cd();      
 
  if (fChain == 0) return;
  
  Int_t nentries = Int_t(fChain->GetEntries()); 
  
  if( nentries > nevents) nentries = nevents;
  
  cout <<  nentries << " Entries"<<endl;
  
  Int_t nbytes = 0, nb = 0;
  
  TF1 *f1all = new TF1("gaussall","gaus",-.250,.250);
  f1all->SetParameters(1.,0.,0.050);
  TRandom rand(0);
  
  for (Int_t jentry=0; jentry<nentries;jentry++) {     
    
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    double w =1;

    // only B+ events for x-feed evaluation
    //    if(isMC!=0&&icat==0&&Gvxbtyp>0) {continue;}

    // redefinition of signal if q2 cut is applied
    if (Q2LOWCUT > 0. || Q2HIGHCUT < 100.)
      if(q2Gen < Q2LOWCUT || q2Gen > Q2HIGHCUT) Gvxbtyp = 0;
    
    //calculate reweightings
    if( !isdata ){
      w =  getBsysweight(Gvxbtyp,vub);                        //Bdec weighting
      w *= getDsysweight(GfDpi,GfDk,GfDks,GfDpiz,GfDlep,dImode,vub); //Ddec weighting
      w *= getBusysweight(Gvxbtyp);                           //bulnudec weighting
//       w *= getTrackingWeight();                               //trk weighting
//       w *= getNeutralWeight();                                //neu weighting
//       w *= getBrecoWeight(intpur);                            //breco weighting  
      if(vub && icat==0)  {
	if(DOTHEO && TMath::Abs(Gvxbtyp)==7 ){
	  if(nres == 0) continue;
	  w *= getTrueMxWeight(mxhadgen,Gvxbtyp); 
	  //	  cout<<"Weight :: "<<getTrueMxWeight(mxhadgen,Gvxbtyp)<<" "<<mxhadgen<<" "<<Gvxbtyp<<endl;
	}     
	if(vub&&(Gvxbtyp!=7&&Gvxbtyp!=-7)&&mes>5.27&&intpur>.5){
	  ((TH1D*)gDirectory->Get("plotres"))->Fill(mxhadgen,w);
	  ((TH1D*)gDirectory->Get("plotall"))->Fill(mxhadgen,w);
	}
	if(vub&&!(Gvxbtyp!=7&&Gvxbtyp!=-7)&&mes>5.27&&intpur>.5){
	  ((TH1D*)gDirectory->Get("plotnres"))->Fill(mxhadgen,w);
	  ((TH1D*)gDirectory->Get("plotall"))->Fill(mxhadgen,w);
	  ((TH1D*)gDirectory->Get("plotallnonres"))->Fill(mxhadgen,1);
	}
      }
    }
    
    // now that the pstarfactor is calculated into the fit no weighting at all
    w = 1;

    // FF reweighting for the signal Vub MC
    // if(isMC!=0&&icat==0) w = mcreweight;

    TLorentzVector p4B_lab(pxBrecoilgen, pyBrecoilgen, pzBrecoilgen, eBrecoilgen);
    TLorentzVector p4Lep_lab(pxLeptongen, pyLeptongen, pzLeptongen, eLeptongen);
    TLorentzVector p4Had_lab(pxMesongen, pyMesongen, pzMesongen, eMesongen);
    TLorentzVector p4VDaug_lab(0., 0., 0., 0.);

    // tag the event with the lepton type
    
    if (!((mxhadfit>0) || (mxhadfit<0) || (mxhadfit == 0))) {
      cout << "MEZZEGA: NAN in MXHADFIT!!" << endl;
      mxhadfit = -999.;
    }   
    
    bool goodrange = 0;
    if(RUN == 1 && run < 17500) goodrange=1;
    if(RUN == 2 && run > 17500) goodrange=1;
    if(RUN == 0) goodrange=1;
    
    // purity range
    if((intpur > MAXINTPUR)||(intpur < MININTPUR)||(pur<MINPUR)) {continue;} 
    // run 1 or run2 ?
    if(isdata && !goodrange) {continue;} 
    
    
    // set of analisys cuts...
    
    // B type
    bool rCh = (TMath::Abs(brecocharge) == BTYPE || BTYPE == 2);
    if(BTYPE>2) {
      rCh = (TMath::Abs(brecocharge) == BTYPE-3);
      if(isMC && (vub||Gvxbtyp == FITCATEGORY) && BTYPE-3 == 0) rCh = rCh * (Gvxbtyp>0);
      if(isMC && (vub||Gvxbtyp == FITCATEGORY) && BTYPE-3 == 1) rCh = rCh * (Gvxbtyp<0);
    }

    // acceptance cuts for leptons
    bool isLp = ((nlept500> 0) && (tlab<2.37) && (tlab>0.36) && (plab>0.5) && rCh);
    
    if(LEPTTYPE == 0)
      isLp = ((nelec500> 0) && (tlab<2.37) && (tlab>0.36) && (plab>0.5) && rCh);
    if(LEPTTYPE == 1) 
      isLp = ((nmu500 > 0) && (tlab<2.37) && (tlab>0.36) && (plab>0.5) && rCh);
    
    // charge correlation
    int flav =  lcharge + brecoflav;
    
    // lepton requests
    bool lPYes = (isLp && !(TMath::Abs(brecocharge)!=0 && (flav)!=0) );
    if(isele == 1) lPYes = lPYes && (pcms > ELECTRONPCUT);
    if(isele == 0) lPYes = lPYes && (pcms > MUONPCUT);

    //    if(FITCATEGORY == -19 && FITCATEGORY == Gvxbtyp) ma0 = ma0 + f1all->GetRandom();
    //    if(FITCATEGORY == 19 && FITCATEGORY == Gvxbtyp) ma0p = ma0p + f1all->GetRandom();
    
    // cut on the number of leptons in the recoil
    bool onelep = nle <= NLEP;
    if(LEPTTYPE == 0)
      onelep = nel <= NLEP;
    if(LEPTTYPE == 1) 
      onelep = nmu <= NLEP;

    // flavor category 
    //   3 = charged B,  4 = neutral B OS,  5 = neutral B SS     
    double flavB = 5;
    if(TMath::Abs(brecocharge)==1)flavB = 3;
    if(TMath::Abs(brecocharge)==0 && flav==0)flavB = 4;
    
    // total charge
    int ch = TMath::Abs(xcharge + brecocharge);  
    
    // wdeltam cut
    bool WdeltaCut = (wdeltam>PRMM2CUT && brecocharge == 0); //last cut based on wdelta
    
    // mm2 (per exclusive mode)
    double themm2=mm2;
    if(FITCATEGORY == -11) themm2 = mm2bestPi0;
    if(FITCATEGORY == 11) themm2 = mm2bestPi;
    if(FITCATEGORY == -12) themm2 = mm2bestEta;
    if(FITCATEGORY == -13) themm2 = mm2bestRho0;
    if(FITCATEGORY == 13) themm2 = mm2bestRho;
    if(FITCATEGORY == -14) themm2 = mm2bestOmega;
    if(FITCATEGORY == -15) themm2 = mm2bestEtap;
    if(FITCATEGORY == -19) themm2 = mm2a0;
    if(FITCATEGORY == 19) themm2 = mm2a0p;
    bool Mm2 = (themm2 < MNUSQHIGH && themm2 > MNUSQLOW);
    if(FITCATEGORY == -11 && nneu ==1)  MM2 = (themm2 < 2.);
    
    // number of charged tracks
    bool NCHG = nchg-1 >= NCHGLOWEXCL && nchg-1 <= NCHGHIGHEXCL;
    
    // momentum of the daughters
    bool MOM1(1), MOM2(1); 
    if(FITCATEGORY == -11){
      MOM1 = Estar1dauPi0[indexbestPi0]>MOM1MIN;
      MOM2 = Estar2dauPi0[indexbestPi0]>MOM2MIN;        
    }
    if(FITCATEGORY == -13){
      MOM1 = Estar1dauRho0[indexbestRho0]>MOM1MIN;
      MOM2 = Estar2dauRho0[indexbestRho0]>MOM2MIN;        
    }
    if(FITCATEGORY == 13){
      MOM1 = PimomdauRho[indexbestRho]>MOM1MIN;
      MOM2 = Pi0momdauRho[indexbestRho]>MOM2MIN;        
    }
    if(FITCATEGORY == -14){
      MOM1 = Pi1momdauOmega[indexbestOmega]>MOM1MIN;
      MOM2 = Pi2momdauOmega[indexbestOmega]>MOM2MIN;        
    }

    // kaon rejection
    bool isK = 0;
    if(KSELE) isK = nkp; // to be improved
    if(KSELE && FITCATEGORY == -14) isK = nkp + nks;
    // mass cut
    bool Mass;
    double mass = -999;
    if(FITCATEGORY == -11){
      mass = barembestPi0;
     }
    if(FITCATEGORY == 11){
      mass = 0.13957;
    }
    if(FITCATEGORY == -12){
      mass = barembestEta;
    }
    if(FITCATEGORY == -13){
      mass = barembestRho0;
    }
    if(FITCATEGORY == 13){
      mass = barembestRho;
    }
    if(FITCATEGORY == -14){
      mass = barembestOmega;
    }
    if(FITCATEGORY == -15){
      mass = barembestEtap;
    }
    if(FITCATEGORY == -19){
      mass = ma0;
    }
    if(FITCATEGORY == 19){
      mass = ma0p;
    }
    Mass = (mass < MXCUTHIGHEXCL && mass > MXCUTLOWEXCL);
    
    if(FITCATEGORY == -12){ 
      if(modeEta[indexbestEta]==1) Mass = (mass < MXCUTHIGHEXCL1 && mass > MXCUTLOWEXCL1); 
      if(modeEta[indexbestEta]==2) Mass = (mass < MXCUTHIGHEXCL2 && mass > MXCUTLOWEXCL2); 
      if(modeEta[indexbestEta]==3) Mass = (mass < MXCUTHIGHEXCL3 && mass > MXCUTLOWEXCL3); 
    } 

    if(FITCATEGORY == -15){ 
      if(modeEtap[indexbestEtap]==1) Mass = (mass < MXCUTHIGHEXCL1 && mass > MXCUTLOWEXCL1); 
      if(modeEtap[indexbestEtap]==2) Mass = (mass < MXCUTHIGHEXCL2 && mass > MXCUTLOWEXCL2); 
      if(modeEtap[indexbestEtap]==3) Mass = (mass < MXCUTHIGHEXCL3 && mass > MXCUTLOWEXCL3); 
      if(modeEtap[indexbestEtap]==4) Mass = (mass < MXCUTHIGHEXCL4 && mass > MXCUTLOWEXCL4); 
    } 

    // number of combinations
    bool Comb = 1;
    if(FITCATEGORY == -13){
      Comb =  nrecoRho0 >= NCOMBLOWEXCL && nrecoRho0 <= NCOMBHIGHEXCL;
    }
    if(FITCATEGORY == 13){
      Comb =  nrecoRho >= NCOMBLOWEXCL && nrecoRho <= NCOMBHIGHEXCL;
    }
    if(FITCATEGORY == -14){
      Comb =  nrecoOmega >= NCOMBLOWEXCL && nrecoOmega <= NCOMBHIGHEXCL;
    }
    
    // number of pi0s
    bool Piz;
    //Piz =  nrecoPi0 >= NPI0LOWEXCL && nrecoPi0 <= NPI0HIGHEXCL;
    Piz = 1;

    // Dalitz plot
    bool DALITZ = 1;
    //    if(FITCATEGORY == -14){
    //      double xbindal = int(dalitzpi1pi0ome/0.016);
    //      double ybindal = int(dalitzpi1pi2ome/0.016);
    //      if(xbindal<0) xbindal = 1;
    //      if(ybindal<0) ybindal = 1;      
    //      if(mydalitz.GetBinContent(xbindal,ybindal)<DALITZCUT) DALITZ = 0;
    //    }     

    // mm2 cut in order to remove Vub crossfeed
    bool MM2CROSS = 1;
    bool MM2PI0 = (mm2bestPi0 >= MNUSQPI0HIGH || mm2bestPi0 <= MNUSQPI0LOW);
    bool MM2ETA = (mm2bestEta >= MNUSQETAHIGH || mm2bestEta <= MNUSQETALOW);
    bool MM2RHO = (mm2bestRho >= MNUSQRHOHIGH || mm2bestRho <= MNUSQRHOLOW);
    bool MM2RHO0 = (mm2bestRho0 >= MNUSQRHO0HIGH || mm2bestRho0 <= MNUSQRHO0LOW);
    bool MM2OMEGA = (mm2bestOmega >= MNUSQOMEGAHIGH || mm2bestOmega <= MNUSQOMEGALOW);
    MM2CROSS = MM2PI0*MM2ETA*MM2RHO0*MM2RHO*MM2OMEGA;
    //MM2CROSS = MM2PI0*MM2ETA;

    // Exact charged tracks on the recoil
    if(DORIGHTNCHG){
      if(FITCATEGORY == -12){
	NCHG = 0;
	if(modeEta[indexbestEta] == 2) { NCHG = nchg-1 == 2; }
	else { NCHG = nchg-1 == 0; }
      }
      if(FITCATEGORY == -15){
	NCHG = 0;
	if(modeEtap[indexbestEtap] == 3) { NCHG = nchg-1 == 4; }
	else if(modeEtap[indexbestEtap]==5||modeEtap[indexbestEtap]==7) { NCHG = nchg-1 == 0; }
	else { NCHG = nchg-1 == 2; }
      }
      if(FITCATEGORY == 19){
	NCHG = 0;
	if(modea0p == 2) { NCHG = nchg-1 == 3; }
	else { NCHG = nchg-1 == 1; }
      }
      if(FITCATEGORY == -19){
	NCHG = 0;
	if(modea0 == 2) { NCHG = nchg-1 == 2; }
	else { NCHG = nchg-1 == 0; }
      }
    }		
    
    // cut on minimum momentum for etap -> rhogamma
    bool DAUGAMMOM = 1;    
    if(FITCATEGORY == -15) DAUGAMMOM = GammamomdauEtap[indexbestEtap]>DAUGAMMAMOM || GammamomdauEtap[indexbestEtap]<-10.;

    // cut on eta mass is daughter
    bool DAUETA = 1;
    if(FITCATEGORY == -15 &&  modeEtap[indexbestEtap]  != 1) DAUETA = TMath::Abs(EtamassdauEtap[indexbestEtap]-0.54775)<DAUETAMASS;
    if(FITCATEGORY == -15 &&  modeEtap[indexbestEtap]  == 2) DAUETA = TMath::Abs(EtamassdauEtap[indexbestEtap]-0.54775)<DAUETAMASS2; 
    if(FITCATEGORY == -15 &&  modeEtap[indexbestEtap]  == 3) DAUETA = TMath::Abs(EtamassdauEtap[indexbestEtap]-0.54775)<DAUETAMASS3; 
    if(FITCATEGORY == -15 &&  modeEtap[indexbestEtap]  == 4) DAUETA = TMath::Abs(EtamassdauEtap[indexbestEtap]-0.54775)<DAUETAMASS4; 


    if(FITCATEGORY == -19) DAUETA = TMath::Abs(a0massetadau-0.54775)<DAUETAMASS; 
    if(FITCATEGORY == 19) DAUETA = TMath::Abs(a0pmassetadau-0.54775)<DAUETAMASS;

    // cut on rho mass is daughter
    bool DAURHO = 1;
    if(FITCATEGORY == -15 && modeEtap[indexbestEtap] == 1) DAURHO = TMath::Abs(Rho0massdauEtap[indexbestEtap]-0.775)<DAURHOMASS;

    // additional cuts for pi l nu
    bool exclcuts=1;
    if(FITCATEGORY == 11){
     bool Pich = (nrecoPi > 0); 
      exclcuts= Pich && (Eneualt<MAXENEU ) && (TMath::Abs(baremassjpsiPi[indexbestPi]-3.1)>JPSIWIN);
      exclcuts = exclcuts && ((chbestPi*lcharge)<0);
    }

    if(FITCATEGORY == -12){
      exclcuts = exclcuts && (nrecoEta>0);
    }    

    if(FITCATEGORY == 13){ 
      exclcuts = exclcuts && (nrecoRho>0); 
    }     

    if(FITCATEGORY == -13){ 
      exclcuts = exclcuts && (nrecoRho0>0); 
    }     

    if(FITCATEGORY == -14){ 
      exclcuts = exclcuts && (nrecoOmega>0); 
    }     

    if(FITCATEGORY == -15){
          exclcuts = exclcuts && (nrecoEtap>0);
     }    


    // q2 cut
    double q2 = -999.;
    if(FITCATEGORY == 11){
      q2 = q2bestPi;
    }
    if(FITCATEGORY == -11){
      q2 = q2bestPi0;
    }
    if(FITCATEGORY == -12){
      q2 = q2bestEta;
    }
    if(FITCATEGORY == 13){ 
      q2 = q2bestRho; 
    } 
    if(FITCATEGORY == -13){ 
      q2 = q2bestRho0; 
    } 
    if(FITCATEGORY == -14){ 
      q2 = q2bestOmega; 
    } 
    if(FITCATEGORY == -15){
      q2 = q2bestEtap;
    }


    // PI L NU
    // Loop on all selected Pi candidates ONLY for events with more than ONE Pi
    // Select best Pi using mm2 only for Pi candidates passing other analysis cuts
    
    if(FITCATEGORY==11 && nrecoPi>1){
      
      int nPiCand = 0;
      int indexPiGoodCand[100];
      
      int chargePiCand = 0;
      double massJPsi = 0.;
      
      for(int k=0; k<100; k++){ indexPiGoodCand[k] = 999; }
      
      for(int i=0; i<nrecoPi; i++){
	
	bool AllPiCuts = (TMath::Abs(baremassjpsiPi[i]-3.1)>JPSIWIN) && (chPi[i]*lcharge)<0;
	AllPiCuts = AllPiCuts && q2Pi[i]>Q2LOWCUT && q2Pi[i]<Q2HIGHCUT;
	
	if(AllPiCuts){
	  indexPiGoodCand[nPiCand] = i;
	  nPiCand++;
	}
	
      }
      
      double mm2tmpbest = 999.;
      int indexbestPiCand = 999;
      int indextmp = 999;
      
      for(int j=0; j<nPiCand; j++){
	indextmp = indexPiGoodCand[j];
	
	if(TMath::Abs(mm2Pi[indextmp]) < mm2tmpbest){
	  mm2tmpbest = TMath::Abs(mm2Pi[indextmp]);
	  indexbestPiCand = indextmp;
	}
	
      }
      
      themm2 = mm2Pi[indexbestPiCand];
      if(indexbestPiCand == 999) themm2 = 999.;
      Mm2 = (themm2 < MNUSQHIGH && themm2 > MNUSQLOW);
      
      chargePiCand = chPi[indexbestPiCand];
      if(indexbestPiCand == 999) chargePiCand = 0;
      
      q2 = q2Pi[indexbestPiCand];
      if(indexbestPiCand == 999) q2 = -999.;
      
      massJPsi = baremassjpsiPi[indexbestPiCand];
      if(indexbestPiCand == 999) massJPsi = 999.;
      
      exclcuts= (Eneualt<MAXENEU ) && (TMath::Abs(massJPsi-3.1)>JPSIWIN);
      exclcuts = exclcuts && (chargePiCand*lcharge)<0;
      
    }
    
    
    // PI0 L NU
    // Loop on all selected Pi0 candidates ONLY for events with more than ONE Pi0
    // Select best Pi0 using mm2 only for Pi0 candidates passing other analysis cuts
    
    if(FITCATEGORY==-11 && nrecoPi0>1){
      
      int nPi0Cand = 0;
      int indexPi0GoodCand[100];
      
      for(int k=0; k<100; k++){	indexPi0GoodCand[k] = 999; }
      
      for(int i=0; i<nrecoPi0; i++){
	
	MOM1 = Estar1dauPi0[i] > MOM1MIN;
	
	Mass = (baremPi0[i] < MXCUTHIGHEXCL && baremPi0[i] > MXCUTLOWEXCL);
	
	bool AllPi0Cuts = MOM1 && Mass && q2Pi0[i]>Q2LOWCUT && q2Pi0[i]<Q2HIGHCUT;
	
	if(AllPi0Cuts){	  
	  indexPi0GoodCand[nPi0Cand] = i;
	  nPi0Cand++;
	}
	
      }
      
      double mm2tmpbest = 999.;
      int indexbestPi0Cand = 999;
      int indextmp = 999;
      
      for(int j=0; j<nPi0Cand; j++){
	indextmp = indexPi0GoodCand[j];
	
	if(TMath::Abs(mm2Pi0[indextmp]) < mm2tmpbest){
	  mm2tmpbest = TMath::Abs(mm2Pi0[indextmp]);
	  indexbestPi0Cand = indextmp;
	}
	
      }
      
      themm2 = mm2Pi0[indexbestPi0Cand];
      if(indexbestPi0Cand == 999) themm2 = 999.;
      Mm2 = (themm2 < MNUSQHIGH && themm2 > MNUSQLOW);
      
      MOM1 = Estar1dauPi0[indexbestPi0Cand] > MOM1MIN;
      
      mass = baremPi0[indexbestPi0Cand];
      if(indexbestPi0Cand == 999) mass = 999.;
      Mass = (mass < MXCUTHIGHEXCL && mass > MXCUTLOWEXCL);
      
      q2 = q2Pi0[indexbestPi0Cand];
      if(indexbestPi0Cand == 999) q2 = -999.;
      
    }
    
    
    // ETA L NU
    // Loop on all selected Eta candidates ONLY for events with more than ONE Eta
    // Select best Eta using mm2 only for Eta candidates passing other analysis cuts

    
    if(FITCATEGORY==-12 && nrecoEta>1){
      
      int nEtaCand = 0;
      int indexEtaGoodCand[100];
      
      for(int k=0; k<100; k++){	indexEtaGoodCand[k] = 999; }
      
      for(int i=0; i<nrecoEta; i++){
	
	if(modeEta[i]==1) { 
	  NCHG = nchg-1 == 0; 
	  Mass = (baremEta[i] < MXCUTHIGHEXCL1 && baremEta[i] > MXCUTLOWEXCL1); 
	}
	
	if(modeEta[i]==2) { 
	  NCHG = nchg-1 == 2; 
	  Mass = (baremEta[i] < MXCUTHIGHEXCL2 && baremEta[i] > MXCUTLOWEXCL2); 
	}
	
	if(modeEta[i]==3) {
	  NCHG = nchg-1 == 0; 
	  Mass = (baremEta[i] < MXCUTHIGHEXCL3 && baremEta[i] > MXCUTLOWEXCL3); 
	}
	
	bool AllEtaCuts = NCHG && Mass && q2Eta[i]>Q2LOWCUT && q2Eta[i]<Q2HIGHCUT;
	
	if(AllEtaCuts){	  
	  indexEtaGoodCand[nEtaCand] = i;
	  nEtaCand++;
	}
	
      }
      
      double mm2tmpbest = 999.;
      int indexbestEtaCand = 999;
      int indextmp = 999;
      
      for(int j=0; j<nEtaCand; j++){
	indextmp = indexEtaGoodCand[j];
	
	if(TMath::Abs(mm2Eta[indextmp]) < mm2tmpbest){
	  mm2tmpbest = TMath::Abs(mm2Eta[indextmp]);
	  indexbestEtaCand = indextmp;
	}
	
      }
      
      themm2 = mm2Eta[indexbestEtaCand];
      if(indexbestEtaCand == 999) themm2 = 999.;
      Mm2 = (themm2 < MNUSQHIGH && themm2 > MNUSQLOW);
      
      mass = baremEta[indexbestEtaCand];
      if(indexbestEtaCand == 999) mass = 999.;    
      
      if(modeEta[indexbestEtaCand]==1) {
	NCHG = nchg-1 == 0; 
	Mass = (mass < MXCUTHIGHEXCL1 && mass > MXCUTLOWEXCL1);
      }
      
      if(modeEta[indexbestEtaCand]==2) {
	NCHG = nchg-1 == 2; 
	Mass = (mass < MXCUTHIGHEXCL2 && mass > MXCUTLOWEXCL2);
      }
      
      if(modeEta[indexbestEtaCand]==3) {
	NCHG = nchg-1 == 0; 
	Mass = (mass < MXCUTHIGHEXCL3 && mass > MXCUTLOWEXCL3);
      }      
      
      q2 = q2Eta[indexbestEtaCand]; 
      if(indexbestEtaCand == 999) q2 = -999.; 
      
    }
   
    
    // ETA PRIME L NU
    // Loop on all selected Etap candidates ONLY for events with more than ONE Etap
    // Select best Etap using mm2 only for Etap candidates passing other analysis cuts
    
    if(FITCATEGORY==-15 && nrecoEtap>1){
      
      int nEtapCand = 0;
      int indexEtapGoodCand[100];
      DAURHO = 1;
      DAUETA = 1;      
      
      for(int k=0; k<100; k++){	indexEtapGoodCand[k] = 999; }
      
      for(int i=0; i<nrecoEtap; i++){
	
	DAUGAMMOM = GammamomdauEtap[i]>DAUGAMMAMOM || GammamomdauEtap[i]<-10.;
	
	if(modeEtap[i]==1) { 
	  NCHG = nchg-1 == 2; 
	  Mass = (baremEtap[i] < MXCUTHIGHEXCL1 && baremEtap[i] > MXCUTLOWEXCL1); 
	  DAURHO = TMath::Abs(Rho0massdauEtap[i]-0.775)<DAURHOMASS;
	}
	
	if(modeEtap[i]==2) { 
	  NCHG = nchg-1 == 2; 
	  Mass = (baremEtap[i] < MXCUTHIGHEXCL2 && baremEtap[i] > MXCUTLOWEXCL2); 
	  DAUETA = TMath::Abs(EtamassdauEtap[i]-0.54775)<DAUETAMASS2; 
	}
	
	if(modeEtap[i]==3) {
	  NCHG = nchg-1 == 4; 
	  Mass = (baremEtap[i] < MXCUTHIGHEXCL3 && baremEtap[i] > MXCUTLOWEXCL3);
	  DAUETA = TMath::Abs(EtamassdauEtap[i]-0.54775)<DAUETAMASS3; 
	}
	
	if(modeEtap[i]==4) {
	  NCHG = nchg-1 == 2; 
	  Mass = (baremEtap[i] < MXCUTHIGHEXCL4 && baremEtap[i] > MXCUTLOWEXCL4); 
	  DAUETA = TMath::Abs(EtamassdauEtap[i]-0.54775)<DAUETAMASS4; 
	}
	
	bool AllEtapCuts = NCHG && Mass && DAURHO && DAUETA && DAUGAMMAMOM && q2Etap[i]>Q2LOWCUT && q2Etap[i]<Q2HIGHCUT;
	
	if(AllEtapCuts){	  
	  indexEtapGoodCand[nEtapCand] = i;
	  nEtapCand++;
	}
	
      }
      
      double mm2tmpbest = 999.;
      int indexbestEtapCand = 999;
      int indextmp = 999;
      
      for(int j=0; j<nEtapCand; j++){
	indextmp = indexEtapGoodCand[j];
	
	if(TMath::Abs(mm2Etap[indextmp]) < mm2tmpbest){
	  mm2tmpbest = TMath::Abs(mm2Etap[indextmp]);
	  indexbestEtapCand = indextmp;
	}
	
      }
      
      themm2 = mm2Etap[indexbestEtapCand];
      if(indexbestEtapCand == 999) themm2 = 999.;
      Mm2 = (themm2 < MNUSQHIGH && themm2 > MNUSQLOW);
      
      mass = baremEtap[indexbestEtapCand];
      if(indexbestEtapCand == 999) mass = 999.;    
      
      DAUGAMMOM = GammamomdauEtap[indexbestEtapCand]>DAUGAMMAMOM || GammamomdauEtap[indexbestEtapCand]<-10.;
      
      DAURHO = 1;
      DAUETA = 1;
      
      if(modeEtap[indexbestEtapCand]==1) {
	NCHG = nchg-1 == 2; 
	Mass = (mass < MXCUTHIGHEXCL1 && mass > MXCUTLOWEXCL1);
	DAURHO = TMath::Abs(Rho0massdauEtap[indexbestEtapCand]-0.775)<DAURHOMASS;
      }
      
      if(modeEtap[indexbestEtapCand]==2) {
	NCHG = nchg-1 == 2; 
	Mass = (mass < MXCUTHIGHEXCL2 && mass > MXCUTLOWEXCL2);
	DAUETA = TMath::Abs(EtamassdauEtap[indexbestEtap]-0.54775)<DAUETAMASS2; 
      }
      
      if(modeEtap[indexbestEtapCand]==3) {
	NCHG = nchg-1 == 4; 
	Mass = (mass < MXCUTHIGHEXCL3 && mass > MXCUTLOWEXCL3);
	DAUETA = TMath::Abs(EtamassdauEtap[indexbestEtap]-0.54775)<DAUETAMASS3; 
      }      
      
      if(modeEtap[indexbestEtapCand]==4) {
	NCHG = nchg-1 == 2; 
	Mass = (mass < MXCUTHIGHEXCL4 && mass > MXCUTLOWEXCL4);
	DAUETA = TMath::Abs(EtamassdauEtap[indexbestEtap]-0.54775)<DAUETAMASS4; 
      }      
      
      q2 = q2Etap[indexbestEtapCand]; 
      if(indexbestEtapCand == 999) q2 = -999.; 
      
    }


    // RHO L NU
    // Loop on all selected Rho candidates ONLY for events with more than ONE Rho
    // Select best Rho using mm2 only for Rho candidates passing other analysis cuts

    
    if(FITCATEGORY==13 && nrecoRho>1){
      
      int nRhoCand = 0;
      int indexRhoGoodCand[100];
      
      for(int k=0; k<100; k++){	indexRhoGoodCand[k] = 999; }
      
      for(int i=0; i<nrecoRho; i++){
	
	MOM2 = 	Pi0momdauRho[i] > MOM2MIN;
	Mass = (baremRho[i] < MXCUTHIGHEXCL && baremRho[i] > MXCUTLOWEXCL);

	bool AllRhoCuts = NCHG && Mass && !isK && q2Eta[i]>Q2LOWCUT && q2Eta[i]<Q2HIGHCUT;
	
	if(AllRhoCuts){	  
	  indexRhoGoodCand[nRhoCand] = i;
	  nRhoCand++;
	}
	
      }
      
      double mm2tmpbest = 999.;
      int indexbestRhoCand = 999;
      int indextmp = 999;
      
      for(int j=0; j<nRhoCand; j++){
	indextmp = indexRhoGoodCand[j];
	
	if(TMath::Abs(mm2Rho[indextmp]) < mm2tmpbest){
	  mm2tmpbest = TMath::Abs(mm2Rho[indextmp]);
	  indexbestRhoCand = indextmp;
	}
	
      }
      
      themm2 = mm2Rho[indexbestRhoCand];
      if(indexbestRhoCand == 999) themm2 = 999.;
      Mm2 = (themm2 < MNUSQHIGH && themm2 > MNUSQLOW);

      MOM1 = Pi0momdauRho[indexbestRhoCand] > MOM1MIN;
      
      mass = baremRho[indexbestRhoCand];
      if(indexbestRhoCand == 999) mass = 999.;    
      Mass = (mass < MXCUTHIGHEXCL && mass > MXCUTLOWEXCL);      
	
      q2 = q2Rho[indexbestRhoCand]; 
      if(indexbestRhoCand == 999) q2 = -999.; 
      
    }
    
    double thevar = 0;
    
    //    cout << q2 << " " << Q2HIGHCUT << " " << Q2LOWCUT << endl;

    // ALL TOGETHER no var
    bool AllCutnovar = (lPYes && onelep && Mm2 &&  ch < CHHIGH && ch > CHLOW && !(WdeltaCut) && NCHG && MOM1 && MOM2 && Comb && Piz && !isK && DALITZ && MM2CROSS && DAUETA && DAURHO && DAUGAMMOM && exclcuts);  
    if(!strcmp(VAR, "mass")) { 
      AllCutnovar = AllCutnovar && q2>Q2LOWCUT && q2<Q2HIGHCUT;
      thevar = mass;
    } else    if(!strcmp(VAR, "q2")) { 
      AllCutnovar = AllCutnovar && Mass;
      thevar = q2;
    } else    if(!strcmp(VAR, "pcms")) { 
      thevar = pcms;
    } else {
      AllCutnovar = AllCutnovar && q2>Q2LOWCUT && q2<Q2HIGHCUT;
      thevar = mass;      
    }
  
    // ALL TOGETHER no mm2
    bool AllCutnomm2 = (lPYes && onelep && q2>Q2LOWCUT && q2<Q2HIGHCUT &&  ch < CHHIGH && ch > CHLOW && !(WdeltaCut) && NCHG && MOM1 && MOM2 && Mass && Comb && Piz && !isK && DALITZ && MM2CROSS && DAUETA && DAURHO && DAUGAMMOM && exclcuts);  

    // ALL TOGETHER
    bool AllCut = (lPYes && onelep && q2>Q2LOWCUT && q2<Q2HIGHCUT && Mm2 &&  ch < CHHIGH && ch > CHLOW && !(WdeltaCut) && NCHG && MOM1 && MOM2 && Mass && Comb && Piz && !isK && DALITZ && MM2CROSS && DAUETA && DAURHO && DAUGAMMOM && exclcuts);  

    if(mes>5.27 && AllCut && isMC && Gvxbtyp!=FITCATEGORY){

      if(icat == 0){

	if(Gvxbtyp>-99){
	  cout << "# VUB BKG: " << Gvxbtyp << " " << GfDpi << " " << GfDk << " " << GfDks << " " << GfDpiz << " " << GfDlep << " " << vub << " " << vcb << " " << mes << endl;
	} else {
          cout << "# OTHER BKG: " << Gvxbtyp << " " << GfDpi << " " << GfDk << " " << GfDks << " " << GfDpiz << " " << GfDlep << " " << vub << " " << vcb << " " << mes << endl; 
	}

      } else {

        if(Gvxbtyp>-99){ 
        cout << "# VCB BKG: " << Gvxbtyp << " " << GfDpi << " " << GfDk << " " << GfDks << " " << GfDpiz << " " << GfDlep << " " << vub << " " << vcb << " " << mes << endl;
	} else {
	  cout << "# OTHER BKG: " << Gvxbtyp << " " << GfDpi << " " << GfDk << " " << GfDks << " " << GfDpiz << " " << GfDlep << " " << vub << " " << vcb << " " << mes << endl; 
	}

      }
    }
    

    if(mes>5.2) {
      
      varMes->setVal(mes);
      varVar->setVal(thevar);
      varMm2->setVal(themm2);	
      varFlavB->setVal(flavB);	      
      varLepYes->setVal(lPYes);	      
      varAllnovar->setVal(AllCutnovar);	      
      varAllnomm2->setVal(AllCutnomm2);	      
      varAllcuts->setVal(AllCut);	      
      varWe->setVal(w);	      
      
      // define MC type
      int type = 0;
      if(isMC) {
	if(vcb && icat!=0) type = 1;
	if(vub && icat ==0) type = 2;
	if((vcb + vub) == 0 && icat!=0) { 
	  type = 3;
	}
      } 
      if(lPYes){
	if(isMC == 0)  {
	  datadata->add(RooArgSet(*varMes,*varVar,*varMm2,*varLepYes,*varFlavB, *varAllnovar, *varAllnomm2, *varAllcuts),1.);
	  if(lPYes) ((TH1D*)gDirectory->Get("nsldata"))->Fill(mes);
	  if(AllCut) ((TH1D*)gDirectory->Get("allcutsdata"))->Fill(mes);
	} else {
	  if(icat!=0){
	    datamc->add(RooArgSet(*varMes,*varLepYes,*varFlavB),w);
	    if(vcb) vcbcounter++;
	    if(vub + vcb != 0)  datamcvcbvub->add(RooArgSet(*varMes,*varLepYes,*varFlavB),w);
	    if(Gvxbtyp == FITCATEGORY) datamcsigforpstar->add(RooArgSet(*varMes,*varLepYes,*varFlavB),w);
	    if(vcb && (AllCutnovar || AllCutnomm2))  {
	      datamcvcb->add(RooArgSet(*varMes,*varVar,*varMm2,*varLepYes,*varFlavB, *varAllnovar, *varAllnomm2, *varAllcuts),w);
	      if(mes>5.2 && AllCut) cout << Gvxbtyp << " " << GfDpi << " " << GfDk << " " << GfDks << " " << GfDpiz << " " << GfDlep << " " << dImode << " " << vub << " " <<  mes << endl;
	    }
	    if(vub + vcb == 0 && (AllCutnovar || AllCutnomm2))  datamcoth->add(RooArgSet(*varMes,*varVar,*varMm2,*varLepYes,*varFlavB, *varAllnovar, *varAllnomm2, *varAllcuts),w);	  
	  }else{
	    if(Gvxbtyp == FITCATEGORY)  datamcsig->add(RooArgSet(*varMes,*varVar,*varMm2,*varLepYes,*varFlavB, *varAllnovar, *varAllnomm2, *varAllcuts),w);
	    if(vub && Gvxbtyp != FITCATEGORY)  {
	      datamcvub->add(RooArgSet(*varMes,*varVar,*varMm2,*varLepYes,*varFlavB, *varAllnovar, *varAllnomm2, *varAllcuts),w);
 	      if(mes>5.2 && AllCut) cout << Gvxbtyp << " " << GfDpi << " " << GfDk << " " << GfDks << " " << GfDpiz << " " << GfDlep << " " << dImode << " " << vub << " " << mes << " " << endl;
	    }
	  }
	}
      }      
    }
  }

}


void exclfitNtp::theFit(){    
  
  gROOT->SetStyle("Plain");
  fHistFile->cd();
  char cuts[100];
  
  // calculate Nsl on MC with and without non-semieleptonic events 
  
  TVector2 signal;
  
  sprintf(cuts,"lepYes");  
  signal = mixCorr(datamc,cuts,mesNslMC,1,"nslMC");
  double totmc = signal.X();
  double errtotmc = signal.Y();

  signal = mixCorr(datamcvcbvub,cuts,mesNslMC,1,"nslMCnobkg");
  double totmcnobkg = signal.X();
  double errtotmcnobkg = signal.Y();

  // calculate Nsl on data
  signal = mixCorr(datadata,cuts,mesNsl,1,"nsl");
  double tot = signal.X();
  double errtot = signal.Y();

  double ndata, nerrdata;
  double nsiglep, nerrsiglep;
  double nsiglepforpstar, nerrsiglepforpstar;
  double nvub, nerrvub;
  double nvublep, nerrvublep;
  double nvcb, nerrvcb;
  double noth, nerroth;

  // calculate nsl on signal
  sprintf(cuts,"lepYes");  
  signal = mixCorr(datamcsig,cuts,messigleptcuts,1,"nslsig");
  nsiglep = signal.X();
  nerrsiglep = signal.Y();

  // calculate number of signal event after all cuts
  sprintf(cuts,"allcuts");
  signal = mixCorr(datamcsig,cuts,messigcuts,1,"sigallcuts");
  nsig = signal.X();
  nerrsig = signal.Y();

  // calculate nsl on vub events      
  sprintf(cuts,"lepYes");
  signal = mixCorr(datamcvub,cuts,mesvubcuts,1,"nslvub");    // mes par to be fixed
  nvublep = signal.X();
  nerrvublep = signal.Y();

  // calculate number of vub event after all cuts      
  sprintf(cuts,"allcuts");
  signal = mixCorr(datamcvub,cuts,mesvubcuts,1,"vuballcuts");
  nvub = signal.X();
  nerrvub = signal.Y();

  // calculate number of vcb event after all cuts      
  sprintf(cuts,"allcuts");
  signal = mixCorr(datamcvcb,cuts,mesvcbcuts,1,"vcballcuts");
  nvcb = signal.X();
  nerrvcb = signal.Y();

// calculate number of other event after all cuts      
  sprintf(cuts,"allcuts");
  signal = mixCorr(datamcoth,cuts,mesothcuts,1,"othallcuts");
  noth = signal.X();
  nerroth = signal.Y();

  // calculate number of events on data after all cuts      
  sprintf(cuts,"allcuts");
  signal = mixCorr(datadata,cuts,mesdatacuts,1,"dataallcuts");
  ndata = signal.X();
  nerrdata = signal.Y();
  
  char normcuts[300], var[100];
  double nvubmm2, nerrvubmm2;
  double nvcbmm2, nerrvcbmm2;
  double nothmm2, nerrothmm2;
  double ndatamm2, nerrdatamm2;

   // calculate number of vub bkg events in 1<mm2<4 GeV 
  double mincut = 1;
  double maxcut = 4;
  sprintf(var,"mm2");
  sprintf(cuts,"allnomm2");
  sprintf(normcuts,"%s%s%s%s%f%s%s%s%f",cuts,"&&",var,">",mincut,"&&",var,"<",maxcut); 
  signal = mixCorr(datamcvub, normcuts, mesvubcuts,1,"vuballcutnorm");
  nvubmm2 = signal.X();
  nerrvubmm2 = signal.Y();

   // calculate number of vcb bkg events in 1<mm2<4 GeV 
  sprintf(var,"mm2");
  sprintf(cuts,"allnomm2");
  sprintf(normcuts,"%s%s%s%s%f%s%s%s%f",cuts,"&&",var,">",mincut,"&&",var,"<",maxcut);  
  signal = mixCorr(datamcvcb, normcuts, mesvcbcuts,1,"vcballcutnorm");
  nvcbmm2 = signal.X();
  nerrvcbmm2 = signal.Y();

  // calculate number of other bkg events in 1<mm2<4 GeV 
  sprintf(var,"mm2");
  sprintf(cuts,"allnomm2");
  sprintf(normcuts,"%s%s%s%s%f%s%s%s%f",cuts,"&&",var,">",mincut,"&&",var,"<",maxcut);  
  signal = mixCorr(datamcoth, normcuts, mesothcuts,1,"othallcutnorm");
  nothmm2 = signal.X();
  nerrothmm2 = signal.Y();

  //calculate number of data events in 1<mm2<4 GeV
  sprintf(var,"mm2");
  sprintf(cuts,"allnomm2");
  sprintf(normcuts,"%s%s%s%s%f%s%s%s%f",cuts,"&&",var,">",mincut,"&&",var,"<",maxcut);
  signal = mixCorr(datadata,normcuts,mesdatacuts,1,"dataallcutnorm");
  ndatamm2 = signal.X();
  nerrdatamm2 = signal.Y();
  
  //*************************************************

  // nsl sig  for p* calculation

  sprintf(cuts,"lepYes");  
  signal = mixCorr(datamcsigforpstar,cuts,messigleptcuts,1,"nslsigforpstar");
  nsiglepforpstar = signal.X();
  nerrsiglepforpstar = signal.Y();


  nsl = tot;
  nslmc = totmc;
  nslvub = nvublep;

  // pstar factor
  //  calcpstarfact = getPstarFactor(LEPTONPCUT);
  double MCSIGBR(0.00018),MCXLNUBR(0.1061);
  cout << "Setting MC branching fractions for SP" << SPTYPE << "...\n";
  if(SPTYPE==8){
    // SP8 values taken from:
    // http://www.slac.stanford.edu/BFROOT/www/Physics/Analysis/AWG/Semileptonic/SPChanges_Overview.txt
    // http://www.slac.stanford.edu/BFROOT/www/Physics/Analysis/AWG/InclusiveSL/common/main.html
    if(FITCATEGORY == 11) { MCSIGBR = 0.000133; MCXLNUBR = 0.1044;} 
    if(FITCATEGORY == -11) { MCSIGBR = 0.000072;  MCXLNUBR = 0.1129;} 
    if(FITCATEGORY == -12) { MCSIGBR = 0.000084; MCXLNUBR = 0.1129;} 
    if(FITCATEGORY == 13) { MCSIGBR = 0.000269; MCXLNUBR = 0.1044;} 
    if(FITCATEGORY == -13) { MCSIGBR = 0.000145; MCXLNUBR = 0.1129;} 
    if(FITCATEGORY == -14) { MCSIGBR = 0.000145; MCXLNUBR = 0.1129;}
    if(FITCATEGORY == -15) { MCSIGBR = 0.000084;  MCXLNUBR = 0.1129;} 
    if(vcbcounter < 1 && FITCATEGORY > 0) { MCXLNUBR = 0.000697;}
    if(vcbcounter < 1 && FITCATEGORY < 0) { MCXLNUBR = 0.000820;}  
  } else {
    // SP5-SP6 values
    if(FITCATEGORY == 11) { MCSIGBR = 0.00018; MCXLNUBR = 0.1061;}
    if(FITCATEGORY == -11) { MCSIGBR = 0.00009;  MCXLNUBR = 0.1061;}
    if(FITCATEGORY == -12) { MCSIGBR = 0.00003; MCXLNUBR = 0.1061;}
    if(FITCATEGORY == 13) { MCSIGBR = 0.00026; MCXLNUBR = 0.1061;} 
    if(FITCATEGORY == -13) { MCSIGBR = 0.00013; MCXLNUBR = 0.1061;}
    if(FITCATEGORY == -14) { MCSIGBR = 0.00013; MCXLNUBR = 0.1061;}
    if(FITCATEGORY == -15) { MCSIGBR = 0.00006;  MCXLNUBR = 0.1061;}
    if(vcbcounter < 1) { MCXLNUBR = 0.000735;}
  }
  calcpstarfact = (nsiglepforpstar/totmcnobkg) * (MCXLNUBR/(MCSIGBR*Q2CORR));
  double errcalcpstarfact = sqrt(((nerrsiglepforpstar*nerrsiglepforpstar)/(totmcnobkg*totmcnobkg)) + nsiglepforpstar*nsiglepforpstar*errtotmcnobkg*errtotmcnobkg/(pow(totmcnobkg,4)));
  errcalcpstarfact = errcalcpstarfact*(MCXLNUBR/(MCSIGBR*Q2CORR));

  // correction for Nsl
  fact = totmcnobkg/totmc; 

  // efficiency for signal
  epssig = nsig/nsiglep;
  errepssig = sqrt(TMath::Abs(epssig*(1 - epssig)/nsiglep)) ;

  // efficiency for vub
  epsvub = nvub/nvublep;
  errepsvub = sqrt(TMath::Abs(epsvub*(1 - epsvub)/nvublep)) ;


  // efficiency for vcb
  epsvcb = nvcb/(nslmc);
  errepsvcb = sqrt(TMath::Abs(epsvcb*(1 - epsvcb)/nslmc)) ;
  
  // efficiency for other
  epsoth = noth/(nslmc);
  errepsoth = sqrt(TMath::Abs(epsoth*(1 - epsoth)/nslmc)) ;

  //efficiency for mm2 cuts

  double epsvubmm2 = nvubmm2/nvublep;
  double errepsvubmm2 = sqrt(TMath::Abs(epsvubmm2*(1- epsvubmm2)/nvublep));

  double epsvcbmm2 = nvcbmm2/nslmc;
  double errepsvcbmm2 = sqrt(TMath::Abs(epsvcbmm2*(1- epsvcbmm2)/nslmc));

  double epsothmm2 = nothmm2/nslmc;
  double errepsothmm2 = sqrt(TMath::Abs(epsothmm2*(1-epsothmm2)/nslmc)); 

  //bkg events on data signal box
  double nvubdata = epsvub * nsl * fact * RATIOBR;
  double nvcbdata = epsvcb * nsl;
  double nothdata = epsoth * nsl;

  //bkg events sideband region

  double nvubdmm2 = epsvubmm2 * nsl * fact * RATIOBR; 
  double nvcbdmm2 = epsvcbmm2 * nsl;
  double nothdmm2 = epsothmm2 * nsl;

  //scale factor data-MC for normalization
  scalefact = ndatamm2/(nvubdmm2+nvcbdmm2+nothdmm2);

  // signal events on data

  nsigdata = ndata -((nvcbdata + nvubdata + nothdata)*scalefact);
  double nbkgdata= (nvcbdata + nvubdata + nothdata)*scalefact;

  double errsigdata2 =(nerrdata*nerrdata)+pow(((nvubdata+nvcbdata+nothdata)*nerrdatamm2/(nvubdmm2+nvcbdmm2+nothdmm2)),2)+pow(ndatamm2*nvubdata,2)/pow(nvubdmm2+nvcbdmm2+nothdmm2,4)*(pow(errepsvubmm2 * nsl * fact * RATIOBR,2) + pow(ndatamm2*nvubdata,2)/pow(nvubdmm2+nvcbdmm2+nothdmm2,4)*pow(errepsvcbmm2 * nsl,2) + pow(ndatamm2*nvubdata,2)/pow(nvubdmm2+nvcbdmm2+nothdmm2,4)*pow(errepsothmm2 * nsl,2));

  double errsigdata = sqrt(errsigdata2);

  // mc error on bkg
  double errMCsigdata = pow(errepsvub * nsl * fact * RATIOBR,2) + pow(errepsvcb * nsl,2) + pow(errepsoth * nsl,2) ;
  errMCsigdata = sqrt(errMCsigdata);

  // signal corrected
  double siglepdata = nsigdata / epssig;
  double errsiglepdata = errsigdata / epssig;
  double errMCsiglepdata = sqrt((errMCsigdata*errMCsigdata/(epssig*epssig)) + nsigdata*nsigdata*errepssig*errepssig/(pow(epssig,4)));

  // calculate BRBR
  double BRBRsig = siglepdata / (tot*fact*calcpstarfact);
  double errBRBRsig = errsiglepdata / (tot*fact*calcpstarfact);
  double errMCBRBRsig = sqrt((errMCsiglepdata*errMCsiglepdata)/(siglepdata*siglepdata)+(errcalcpstarfact*errcalcpstarfact)/(calcpstarfact*calcpstarfact))*(BRBRsig);

  cout << "BRBRsig = " << BRBRsig << " +/- " << errBRBRsig << "(stat) +/- " << errMCBRBRsig  << "(MCstat) " <<   endl;
  

  // draw final plots
   char nameeps[100],namehist[100];

  // for the var
  // signal
  sprintf(var,"var");
  sprintf(cuts,"allnovar");
  sprintf(nameeps,"%s%s",VAR,"signsuballcuts");   
  sprintf(namehist,"%s%s",VAR,"sigexcl");   
  makesubplot(datamcsig,var,cuts,messigcuts,BINSVAR,MINVAR,MAXVAR,namehist,1,nameeps);
  // vub
  sprintf(nameeps,"%s%s",VAR,"vubsuballcuts");   
  sprintf(namehist,"%s%s",VAR,"vubexcl");   
  makesubplot(datamcvub,var,cuts,mesvubcuts,BINSVAR,MINVAR,MAXVAR,namehist,1,nameeps);   
  // vcb							  
  sprintf(nameeps,"%s%s",VAR,"vcbsuballcuts");   		  
  sprintf(namehist,"%s%s",VAR,"vcbexcl");   
  makesubplot(datamcvcb,var,cuts,mesvcbcuts,BINSVAR,MINVAR,MAXVAR,namehist,1,nameeps);   
  // oth							  
  sprintf(nameeps,"%s%s",VAR,"othsuballcuts");   		  
  sprintf(namehist,"%s%s",VAR,"othexcl");   
  makesubplot(datamcoth,var,cuts,mesothcuts,BINSVAR,MINVAR,MAXVAR,namehist,1,nameeps);   
  // data
  sprintf(nameeps,"%s%s",VAR,"datasuballcuts");   
  sprintf(namehist,"%s%s",VAR,"dataexcl");   
  makesubplot(datadata,var,cuts,mesdatacuts,BINSVAR,MINVAR,MAXVAR,namehist,1,nameeps);   

  // let's draw it
  makefinalplots(VAR, BINSVAR);

  
  // for mm2
  // signal
  sprintf(var,"mm2");
  sprintf(cuts,"allnomm2");
  sprintf(nameeps,"%s%s",var,"signsuballcuts");   
  makesubplot(datamcsig,var,cuts,messigcuts,BINSMM2,MINMM2,MAXMM2,"mm2sigexcl",1,nameeps);
  // vub
  sprintf(nameeps,"%s%s",var,"vubsuballcuts");   
  makesubplot(datamcvub,var,cuts,mesvubcuts,BINSMM2,MINMM2,MAXMM2,"mm2vubexcl",1,nameeps);   
  // vcb
  sprintf(nameeps,"%s%s",var,"vcbsuballcuts");   
  makesubplot(datamcvcb,var,cuts,mesvcbcuts,BINSMM2,MINMM2,MAXMM2,"mm2vcbexcl",1,nameeps);   
  // oth
  sprintf(nameeps,"%s%s",var,"othsuballcuts");   
  makesubplot(datamcoth,var,cuts,mesothcuts,BINSMM2,MINMM2,MAXMM2,"mm2othexcl",1,nameeps);   
  // data
  sprintf(nameeps,"%s%s",var,"datasuballcuts");   
  makesubplot(datadata,var,cuts,mesdatacuts,BINSMM2,MINMM2,MAXMM2,"mm2dataexcl",1,nameeps);   

  // let's draw it
  makefinalplots(var, BINSMM2);

  
  // print results on file
  char name[100];
  sprintf(name, "%s%s%s",DIRNAME.Data(),PREFIXOUT.Data(),"result_exclusive.dat");        
  ofstream outfileexcl(name);

  outfileexcl << "#NUMBERS FOR THE SCANS" << endl;
  outfileexcl << endl;
  outfileexcl << endl;
  outfileexcl << " pstarfact           "  << calcpstarfact   <<  " "  << errcalcpstarfact << endl;  
  outfileexcl << " q2lowCut            "  << Q2LOWCUT   <<  endl;  
  outfileexcl << " q2highCut           "  << Q2HIGHCUT   <<  endl;  
  outfileexcl << " q2corr              "  << Q2CORR   <<  endl;  

  outfileexcl << " leptonPCut          "  << LEPTONPCUT   <<  endl;  
  outfileexcl << " electronPCut        "  << ELECTRONPCUT   <<  endl; 
  outfileexcl << " muonPCut            "  << MUONPCUT   <<  endl;  
  outfileexcl << " nlep                "  << NLEP   <<  endl;  
  outfileexcl << " prmm2cut            "  << PRMM2CUT   <<  endl;  
  outfileexcl << " mnuSqLow            "  << MNUSQLOW   <<  endl;  
  outfileexcl << " mnuSqHigh           "  << MNUSQHIGH   <<  endl;  
  outfileexcl << " chLow               "  << CHLOW   <<  endl;  
  outfileexcl << " chHigh              "  << CHHIGH   <<  endl;  
  outfileexcl << " Btype               "  << BTYPE   <<  endl;  
  outfileexcl << " lepttype            "  << LEPTTYPE   <<  endl;  
  outfileexcl << " minpur              "  << MINPUR    <<  endl; 
  outfileexcl << " minintpur           "  << MININTPUR    <<  endl;
  outfileexcl << " maxintpur           "  << MAXINTPUR    <<  endl;
  outfileexcl << " mxCutLowExcl        "  << MXCUTLOWEXCL << endl;
  outfileexcl << " mxCutHighExcl       "  << MXCUTHIGHEXCL << endl;
  outfileexcl << " mxCutLowExcl1        "  << MXCUTLOWEXCL1 << endl;
  outfileexcl << " mxCutHighExcl1       "  << MXCUTHIGHEXCL1 << endl;
  outfileexcl << " mxCutLowExcl2        "  << MXCUTLOWEXCL2 << endl;
  outfileexcl << " mxCutHighExcl2       "  << MXCUTHIGHEXCL2 << endl;
  outfileexcl << " mxCutLowExcl3        "  << MXCUTLOWEXCL3 << endl;
  outfileexcl << " mxCutHighExcl3       "  << MXCUTHIGHEXCL3 << endl;
  outfileexcl << " mxCutLowExcl4        "  << MXCUTLOWEXCL4 << endl;
  outfileexcl << " mxCutHighExcl4       "  << MXCUTHIGHEXCL4 << endl;
  outfileexcl << " chgLowexcl          "  << NCHGLOWEXCL << endl;
  outfileexcl << " chgHighexcl         "  << NCHGHIGHEXCL << endl;
  outfileexcl << " deltam              "  << DELTAM << endl;
  outfileexcl << " ncombLowExcl        "  << NCOMBLOWEXCL << endl;      
  outfileexcl << " ncombHighExcl       "  << NCOMBHIGHEXCL << endl;    
  outfileexcl << " npi0LowExcl         "  << NPI0LOWEXCL << endl;    
  outfileexcl << " npi0HighExcl        "  << NPI0HIGHEXCL << endl;    
  outfileexcl << " mom1min             "  << MOM1MIN << endl;    
  outfileexcl << " mom2min             "  << MOM2MIN << endl;    
  outfileexcl << " ksele               "  << KSELE << endl;    
  outfileexcl << " dalitzcut           "  << DALITZCUT << endl;    
  outfileexcl << " mnuSqpi0Low         "  << MNUSQPI0LOW   <<  endl;  
  outfileexcl << " mnuSqpi0High        "  << MNUSQPI0HIGH   <<  endl;  
  outfileexcl << " mnuSqetaLow         "  << MNUSQETALOW   <<  endl;  
  outfileexcl << " mnuSqetaHigh        "  << MNUSQETAHIGH   <<  endl;  
  outfileexcl << " mnuSqetaLow         "  << MNUSQETALOW   <<  endl;  
  outfileexcl << " mnuSqetaHigh        "  << MNUSQETAHIGH   <<  endl;  
  outfileexcl << " mnuSqrhoLow         "  << MNUSQRHOLOW   <<  endl;  
  outfileexcl << " mnuSqrhoHigh        "  << MNUSQRHOHIGH   <<  endl;  
  outfileexcl << " mnuSqrho0Low        "  << MNUSQRHO0LOW   <<  endl;  
  outfileexcl << " mnuSqrho0High       "  << MNUSQRHO0HIGH   <<  endl;  
  outfileexcl << " mnuSqomegaLow       "  << MNUSQOMEGALOW   <<  endl;  
  outfileexcl << " mnuSqomegaHigh      "  << MNUSQOMEGAHIGH   <<  endl;  
  outfileexcl << " maxeneu             "  << MAXENEU   <<  endl;  
  outfileexcl << " jpsiwindow          "  << JPSIWIN   <<  endl;  
  outfileexcl << " daugammamom         "  << DAUGAMMAMOM << endl;

  outfileexcl << " nsl                 "  << tot << " " << errtot << endl;
  outfileexcl << " nslmc               "  << totmc << " " << errtotmc << endl;
  outfileexcl << " nslmcbkg            "  << totmcnobkg << " " << errtotmcnobkg << endl; 
  outfileexcl << " fact                "  << fact << endl; 
  outfileexcl << " sig                 "  << nsig << "  " << nerrsig <<  endl;   
  outfileexcl << " siglep              "  << nsiglep << "  " << nerrsiglep <<  endl;   
  outfileexcl << " vub                 "  << nvub << "  " << nerrvub <<  endl;   
  outfileexcl << " vublep              "  << nvublep << "  " << nerrvublep <<  endl;    
  outfileexcl << " vcb                 "  << nvcb << "  " << nerrvcb <<  endl;   
  outfileexcl << " oth                 "  << noth << "  " << nerroth <<  endl;   
  outfileexcl << " data                "  << ndata << "  " << nerrdata <<  endl;   
  outfileexcl << " Epssig              "  << epssig << "  " << errepssig <<  endl;   
  outfileexcl << " Epsvub              "  << epsvub << "  " << errepsvub <<  endl;    
  outfileexcl << " Epsvcb              "  << epsvcb << "  " << errepsvcb <<  endl;    
  outfileexcl << " Epsoth              "  << epsoth << "  " << errepsoth <<  endl;    
  outfileexcl << " ratioBR             "  << RATIOBR <<  endl; 
  outfileexcl << " SPtype              "  << SPTYPE << endl;
  outfileexcl << " vubdata             "  << nvubdata << endl;      
  outfileexcl << " vcbdata             "  << nvcbdata << endl;      
  outfileexcl << " othdata             "  << nothdata << endl;     
  outfileexcl << " sigdata             "  << nsigdata << "  " << errsigdata << "  " << errMCsigdata << endl;    
  outfileexcl << " siglepdata          "  << siglepdata << "  " << errsiglepdata << "  " << errMCsiglepdata << endl;   
  outfileexcl << "####### mm2 region #########" << endl;
  outfileexcl << "nvubmm2              "  << nvubmm2   << "  " << nerrvubmm2 << endl;
  outfileexcl << "nvcbmm2              "  << nvcbmm2   << "  " << nerrvcbmm2 << endl;
  outfileexcl << "nothmm2              "  << nothmm2   << "  " << nerrothmm2 << endl;
  outfileexcl << "ndatamm2             "  << ndatamm2  << "  " << nerrdatamm2 << endl;
  outfileexcl << "epsvubmm2            "  << epsvubmm2 << "  " << errepsvubmm2 << endl;
  outfileexcl << "epsvcbmm2            "  << epsvcbmm2 << "  " << errepsvcbmm2 << endl;
  outfileexcl << "epsothmm2            "  << epsothmm2 << "  " << errepsothmm2 << endl;
  outfileexcl << "nvubdmm2             "  << nvubdmm2  << "  " << endl;
  outfileexcl << "nvcbdmm2             "  << nvcbdmm2  << "  " << endl;
  outfileexcl << "nothdmm2             "  << nothdmm2  << "  " << endl;
  outfileexcl << "scalefactor          "  << scalefact << "  " << endl;
  outfileexcl << "nbkgdata             "  << nbkgdata  << "  " << endl;
  outfileexcl << "####### mm2 region #########" << endl;    

  outfileexcl << " BRBRsig             "  << BRBRsig << "  " << errBRBRsig << "  " << errMCBRBRsig << endl;
   
  outfileexcl << endl;
  outfileexcl << endl;
  outfileexcl << "INPUT FILES" << endl;
  outfileexcl << endl;
  outfileexcl << endl;
  outfileexcl << " fileVubTotal        "  << FILEVUBTOTAL   <<  endl;  
  outfileexcl << " fileVcb             "  << FILEVCB   <<  endl;  
  outfileexcl << " fileData            "  << FILEDATA   <<  endl;  

  outfileexcl << endl;
  outfileexcl << endl;
  outfileexcl.close();


  sprintf(name, "%s%s%s",DIRNAME.Data(),PREFIXOUT.Data(),"both.eps");   
  openEpsFile(name);
  ((TH1D*)gDirectory->Get("plotall"))->SetMarkerStyle(3005);
  ((TH1D*)gDirectory->Get("plotall"))->SetLineColor(kRed);
  ((TH1D*)gDirectory->Get("plotall"))->SetFillColor(kRed);
  ((TH1D*)gDirectory->Get("plotall"))->Draw();
  ((TH1D*)gDirectory->Get("plotnres"))->SetMarkerStyle(3006);
  ((TH1D*)gDirectory->Get("plotnres"))->SetLineColor(kBlue);
  ((TH1D*)gDirectory->Get("plotnres"))->SetFillColor(kBlue);
  ((TH1D*)gDirectory->Get("plotnres"))->DrawCopy("same");
  closeEpsFile();
  c1->Clear();
  
  TH2 *h;
  h = new TH2D("inthyb","Mx Fraction ",40,0.,5.,500,0.,1.2);
  h = new TH2D("intnres","Mx Fraction",40,0.,5.,500,0.,1.2);
  
  double integralhyb=0;
  double integralnres=0;  
  double totintegralhyb=((TH1D*)gDirectory->Get("plotall"))->Integral();
  double totintegralnres=((TH1D*)gDirectory->Get("plotallnonres"))->Integral();  
  
  for(int i=1; i<41; i++){
    
    integralhyb+=((TH1D*)gDirectory->Get("plotall"))->GetBinContent(i);
    integralnres+=((TH1D*)gDirectory->Get("plotallnonres"))->GetBinContent(i);
    ((TH1D*)gDirectory->Get("inthyb"))->Fill(i*.125+0.0625,integralhyb/totintegralhyb);
    ((TH1D*)gDirectory->Get("intnres"))->Fill(i*.125+0.0625,integralnres/totintegralnres);
    
  }
  
  sprintf(name, "%s%s%s",DIRNAME.Data(),PREFIXOUT.Data(),"integrals.eps");   
  openEpsFile(name);
  ((TH1D*)gDirectory->Get("inthyb"))->SetMarkerStyle(20);
  ((TH1D*)gDirectory->Get("inthyb"))->SetMarkerSize(.8);
  ((TH1D*)gDirectory->Get("inthyb"))->SetMarkerColor(kBlue); 
  ((TH1D*)gDirectory->Get("inthyb"))->SetStats(0);
  ((TH1D*)gDirectory->Get("intnres"))->SetStats(0);
  ((TH1D*)gDirectory->Get("intnres"))->SetMarkerStyle(20);
  ((TH1D*)gDirectory->Get("intnres"))->SetMarkerSize(.8);
  ((TH1D*)gDirectory->Get("intnres"))->SetMarkerColor(kRed); 
  ((TH1D*)gDirectory->Get("intnres"))->Draw();
  ((TH1D*)gDirectory->Get("inthyb"))->DrawCopy("same");
  TLegend *leg9; 
  TLegendEntry *legge9; 
  leg9 = new TLegend(0.45,0.5,0.69,0.59);
  leg9->SetFillStyle(0); leg9->SetBorderSize(0.); leg9->SetTextSize(0.05); 
  leg9->SetFillColor(0); 
  legge9 = leg9->AddEntry(((TH1D*)gDirectory->Get("inthyb")), "Hybrid", "p"); 
  legge9 = leg9->AddEntry(((TH1D*)gDirectory->Get("intnres")), "Non-res", "p"); 
  leg9->Draw();
  closeEpsFile();
  
//    sprintf(name, "%s%s%s.tex",DIRNAME.Data(),PREFIXOUT.Data(),texPrefix.Data());        
//    cout << "**********************************************************************" << endl;
//    cout << "DUMPING TEXFILE 3 " << name << endl;
//    cout << "**********************************************************************" << endl;
//    ofstream texfile(name);
//    char texline[200];
//    sprintf(texline, "\\def\\%ssemil{%7.0f \\pm %6.0f} ", texPrefix.Data(), tot, errtot); 
//    texfile << texline << endl;
//    sprintf(texline, "\\def\\%sbgsl{%5.0f \\pm %4.0f} ", texPrefix.Data(), (1.-fact)*tot, (1.-fact)*errtot); 
//    texfile << texline << endl;
//    sprintf(texline, "\\def\\%sepssl{%5.2f \\pm %4.2f} ", texPrefix.Data(), calcpstarfact, errcalcpstarfact); 
//    texfile << texline << endl;
//    sprintf(texline, "\\def\\%snu{%5.0f \\pm %4.0f} ", texPrefix.Data(), S, errS); 
//    texfile << texline << endl;
//    sprintf(texline, "\\def\\%sbgc{%5.0f \\pm %4.0f} ", texPrefix.Data(), vcbaftercutsbin1, errfitvcbaftercutsbin1); 
//    texfile << texline << endl;
//    sprintf(texline, "\\def\\%sbgo{%5.0f \\pm %4.0f} ", texPrefix.Data(), otheraftercutsbin1, errfitotheraftercutsbin1); 
//    texfile << texline << endl;
//    sprintf(texline, "\\def\\%sepsu{%4.3f} ", texPrefix.Data(), epsu); 
//    texfile << texline << endl;
//    sprintf(texline, "\\def\\%sepsmx{%4.3f} ", texPrefix.Data(), epsmx); 
//    texfile << texline << endl;
//    sprintf(texline, "\\def\\%sbrbr{%6.4f \\pm %5.4f} ", texPrefix.Data(), BRBR, errBRBR); 
//    texfile << texline << endl;
//    sprintf(texline, "\\def\\%sbrbrerrmc{%5.4f} ", texPrefix.Data(), errBRBRMCstat); 
//    texfile << texline << endl;
//    texfile.close(); 

} 

#include "VirVubFitter/exclfitUtil.icc"
#include "VirVubFitter/exclfitInit.icc"
