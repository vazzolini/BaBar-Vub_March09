#include <fstream.h>

#include "VubAnalysis/fitNtp.hh"
#include "VubAnalysis/recoilDSys.hh"

#include "TH2.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLegendEntry.h"
#include "TLegend.h"
#include <TRandom.h>
#include <TVector2.h>
#include <TLatex.h>
#include <TMinuit.h> 

TH1D* vubHist, *vcbHist, * othHist, *dataHist; // global histograms

// ----------------------------------------------------------------------
fitNtp::fitNtp(TTree *tree, int Sys, TString filename, int prl) {
  Init(tree);
  initRest(filename);
  
  TRandom random(0);
  random.SetSeed(0); 
  int therandom = 0;
  dImode = 2;
  if(Sys > 0) {
    cout<<"Passed here "<< Sys <<endl;

    dImode = Sys;
    therandom = random.Rndm() * 1000000;
    if(Sys == 2) {
      Dvar = new recoilDSys("ddecay.table",therandom,Sys);
    } else {
      Dvar = new recoilDSys("dIdecay.table",therandom,Sys);
    }
    Bsem = new recoilDSys(therandom);
  } else {
    cout<<"Default reweighting "<< Sys <<endl;
    Dvar = new recoilDSys("ddecay.table",therandom,2);
    Bsem = new recoilDSys(therandom);
  }
  fprlRew = prl;
}

// ----------------------------------------------------------------------
fitNtp::~fitNtp() {
}

void fitNtp::Loop(int isdata, int icat, int nevents, int isMC, int nres)

  // flags legenda:
  // isdata  0 = MC to get the shapes and the efficiencies, 1 = fill the histos to fit
  // icat    0 = vub component, 1 = vcb and other component, 
  // nevents = number of events
  // isMC    0 = fit on the MC, 1 = fit on data
{

   fHistFile->cd();
      
   char name[200], preine[100], preich[100], prefix[100];
   char name2[100]; 
   char flavcat[100],lgroup[100],truthgroup[100];
   char le[100];
  // multiplicity categoriy definition 
   double chargecut[2];
   chargecut[0] = .5;
   chargecut[1] = 2.5;
   double neucut = .5;
   int ich;
   int ine;
   int iother(0);
   bool islept;
   double w, wfermi;
   //There you can define the various weights.
   double cutflav = 0;    

   if (fChain == 0) return;
  
   Int_t nentries = Int_t(fChain->GetEntries());

   if( nentries > nevents) nentries = nevents;

   cout <<  nentries << " Entries"<<endl;
   
   Int_t nbytes = 0, nb = 0;

  
   f1all = new TF1("gaussall","gaus",-5*SMEARALLSIGMA,5*SMEARALLSIGMA);
   f1bkg = new TF1("gaussbck","gaus",-5*SMEARBKGSIGMA,5*SMEARBKGSIGMA);
   TRandom rand(0);
   for (Int_t jentry=0; jentry<nentries;jentry++) {     
     Int_t ientry = LoadTree(jentry); //in case of a TChain, ientry is the entry number in the current file
     nb = fChain->GetEntry(jentry);   nbytes += nb;

     double toy(isToyMC());
     if(toy>0 && rand.Rndm()>toy)continue;
     // tag the event for sys. study
     w =1;
     wfermi = 1;
     double wnre = 1;
     //Scan and differential BR for different variables
     if((Q2LOCUT != 0) || (Q2HICUT != 30)) {
       if(q2fit<Q2LOCUT || q2fit>Q2HICUT) continue;
     }
     //Scan and differential BR for different variables
     if((CSILOCUT != 0) || (CSIHICUT != 1)) {
       if(csiCiuc<CSILOCUT || csiCiuc>CSIHICUT) continue;
     }
     if((XLOCUT != 0) || (XHICUT != 1.5)) {
       if(xCiuc<XLOCUT || xCiuc>XHICUT) continue;
     }
     if((WLOCUT != 0) || (WHICUT != 2)) {
       if(wCiuc<WLOCUT || wCiuc>WHICUT) continue;
     }
     if((EWPWLOCUT != 3) || (EWPWHICUT != 5.5)) {
       if(EwPwfit<EWPWLOCUT || EwPwfit>EWPWHICUT) continue;
     }

     //calculate reweightings
     if( !isdata ){
       w =  getBsysweight(Gvxbtyp,vub);                        //Bdec weighting
       w *= getDsysweight(GfDpi,GfDk,GfDks,GfDpiz,GfDlep,dImode,vub); //Ddec weighting
       w *= getTrackingWeight();                               //trk weighting
       w *= getNeutralWeight();                                //neu weighting
       w *= getBrecoWeight(intpur);                            //breco weighting  

       if(DOTHEO && !icat && vub)  {
	 //	 cout<<"Calling weight fun:: "<<mxhadfit<<endl;
	 if(TMath::Abs(Gvxbtyp)==7 ){
	   // just select resonant vub
	   if( nres == 0) continue;
 	   if(FERMIAPP) {
 	     wfermi *= FermiWeight(fkplus,DELTAMB,DELTAA);
	     wnre = wfermi;
 	   }

	   if(fprlRew) {
	     wfermi *= getTrueMxWeight(mxhadgen,Gvxbtyp); 
	     if((DELTAMB)>0.1) {
	       wfermi *= 1.14063009693489213123;
	     }
	     if((DELTAMB)<-0.1) {
	       wfermi *= 0.94723753679861141874;
	     } 
	   } else {
	     int ys = rHistmx(mxhadgen)+ rHistq2(q2Gen)*14 + rHistel(ecmsgen)*14*8;
	     wfermi *= 1.7755*MatrixW[ys];
 	     if(MULFAC != 0) {
	       wfermi *= 1/MULFAC;
	     }
	   }
	 }
       }
     }
     //Non reso 2353432; Reso 2242354
     //Pesi gen:: sta gen = Pesi ric :: sta ric
     //Pesi gen = Pesi Ric * sta gen / sta ric


     wfermi *= getTrackingWeight();
     wfermi *= getNeutralWeight();
     wfermi *= getBrecoWeight(intpur);  
     // tag the event with the multiplicity category
     
     if(nchg-nle<chargecut[0]) {ich = 1;}
     else if(nchg-nle<chargecut[1]){ich = 2;}
     else{ich = 3;};
      
     double neucat = nneu;
     if(MULTIFIT==2) neucat = nneu-nneu80_160;
     if(MULTIFIT==3) neucat = nneu-nneu80_160-nneu160_320;
     if(neucat<neucut) {ine = 1;}
     else{ine = 2;};      
      
     // tag the event with the lepton type
     
     if (!((mxhadfit>0) || (mxhadfit<0) || (mxhadfit == 0))) {
       cout << "MEZZEGA: NAN in MXHADFIT!!" << endl;
       mxhadfit = -999.;
     }   

     //     mxhadfit = mxhad;
     double prevalue = mxhadfit;

     sprintf(preich, "%s%d", "ch", ich);
     sprintf(preine, "%s%d", "ne", ine);
     sprintf(prefix, "%s%s", preich , preine);          

     bool goodrange = 0;
     if(RUN == 1 && run < 17500) goodrange=1;
     if(RUN == 2 && run > 17500) goodrange=1;
     if(RUN == 0) goodrange=1;

     if(TMath::Abs(brecocharge) == BTYPE || BTYPE == 2) {
       
       if(isdata && ((intpur > MAXINTPUR)||(intpur < MININTPUR))) {continue;} 
       if(isdata && nnpi0!=CUTNNPI0 && CUTNNPI0>=0) {continue;} 
       if(isdata && !goodrange) {continue;} 
       
       
       // tag the event category vub, vcb, other 
       // fill histograms on fitted sample as a crosscheck
       
       iother = 0;
       
       if(vcb && icat!=0) {
	 sprintf(lgroup, "vcb");	  
       }
       if(vub && icat==0) {
	 sprintf(lgroup, "vub");	  
       }	
       if(vcb + vub == 0 && icat!=0) {
	 sprintf(lgroup, "oth");	  
	 iother = 1;
       }	
       
       if (isdata != 1){
	 if(icat==0){ 
	   if(vub&&(Gvxbtyp!=7&&Gvxbtyp!=-7)&&mes>5.27&&intpur>.5){
	     ((TH1D*)gDirectory->Get("plotres"))->Fill(mxhadgen,wfermi);
	     ((TH1D*)gDirectory->Get("plotall"))->Fill(mxhadgen,wfermi);
	   }
	   if(vub&&!(Gvxbtyp!=7&&Gvxbtyp!=-7)&&mes>5.27&&intpur>.5){
	     ((TH1D*)gDirectory->Get("plotnres"))->Fill(mxhadgen,wfermi);
	     ((TH1D*)gDirectory->Get("plotall"))->Fill(mxhadgen,wfermi);
	     ((TH1D*)gDirectory->Get("plotallnonres"))->Fill(mxhadgen,wnre);
	   }
	 }
       }
       
       // extract the original number of events in each category
       
       if(isdata){
	 if(FITMC){
	   if(vcb && icat!=0) {
	     ((TH1D*)gDirectory->Get("vcbmcnocutdata"))->Fill(mes); 
	     ((TH1D*)gDirectory->Get("allmcnocutdata"))->Fill(mes);
	   }
	   if(vub && icat==0) {
	     ((TH1D*)gDirectory->Get("vubmcnocutdata"))->Fill(mes);
	     ((TH1D*)gDirectory->Get("allmcnocutdata"))->Fill(mes);
	   }	   
	   if(vcb + vub == 0 && icat!=0) {
	     ((TH1D*)gDirectory->Get("othmcnocutdata"))->Fill(mes);
	     ((TH1D*)gDirectory->Get("allmcnocutdata"))->Fill(mes);
	   }
	 }else{
	   ((TH1D*)gDirectory->Get("allmcnocutdata"))->Fill(mes);
	   if(vcb) ((TH1D*)gDirectory->Get("vcbmcnocutdata"))->Fill(mes); 
	   if(vub) ((TH1D*)gDirectory->Get("vubmcnocutdata"))->Fill(mes); 
	   if(vub+vcb == 0) ((TH1D*)gDirectory->Get("othmcnocutdata"))->Fill(mes); 
	 }
       }else{
	 // number of true vub and vcb on mc to get the pstarfact (eps_u^l/eps_sl^l)
	 if(icat==0){ 
	   if(vub) {		
	     ((TH1D*)gDirectory->Get("vubmcnocut"))->Fill(mes,wfermi); 		
	   }
	 }else{
	   if(vcb){
	     ((TH1D*)gDirectory->Get("vcbmcnocut"))->Fill(mes,w); 	
	   }
	   if(vcb+vub==0){
	     ((TH1D*)gDirectory->Get("othmcnocut"))->Fill(mes,w); 	
	   }
	 }
       }
       
       //id = Mx category == Mx bin
       int id;
       
       //number of leptons cut 
       islept = ((nle > 0) && (tlab<2.37) && (tlab>0.36) && (plab>0.5));  
       if(LEPTTYPE == 0) islept = ((nel > 0) && (tlab<2.37) && (tlab>0.36) && (plab>0.5));
       if(LEPTTYPE == 1) islept = ((nmu > 0) && (tlab<2.37) && (tlab>0.36) && (plab>0.5));
       
       if(pcms > LEPTONPCUT && islept && !(TMath::Abs(brecocharge)!=0 && (lcharge + brecoflav)!=0)){
	 int flav =  lcharge + brecoflav; // charge correlation
	 // flavor category (3 = charged B, 4 = neutral B OS, 5 = neutral B SS)  
	 cutflav = 5;
	 if(TMath::Abs(brecocharge)) cutflav = 3;
	 if(TMath::Abs(brecocharge)==0 && flav==0) cutflav = 4;	  
	 
	 if(cutflav==3) sprintf(flavcat, "%s", "bch");
	 if(cutflav==4) sprintf(flavcat, "%s", "bos");
	 if(cutflav==5) sprintf(flavcat, "%s", "bss");
	 
	 id = hist(mxhadfit);        

	 //fitted sample
	 if(isdata == 1) {
	   //events after lepton cut (on data) to estimate b->clnu contribution 
	   sprintf(name, "%s%s","leptondata",flavcat);
	   if(isMC != 1) {
	     // datalike fit (fitMC = 0) 
	     ((TH1D*)gDirectory->Get(name))->Fill(mes);
	   }else{
	     // MC fit (fitMC = 1) 
	     if(icat==0){ 
	       if(vub) {
		 ((TH1D*)gDirectory->Get(name))->Fill(mes);
	       }
	     }else{
	       if(vcb || iother){
		 ((TH1D*)gDirectory->Get(name))->Fill(mes);
	       }
	     } 
	   }
	   
	   //model
	   
	 }else{
	   //events after lepton cut (on mc) to estimate b->clnu contribution, vub efficiencies,...
	   if(icat==0){ 
	     if(vub) {		
	       sprintf(name, "%s%s%s",lgroup,"leptonmc",flavcat);
	       ((TH1D*)gDirectory->Get(name))->Fill(mes,wfermi);
	       sprintf(name, "%s%s%s%s",lgroup,"leptonmc","eff",flavcat);
	       ((TH1D*)gDirectory->Get(name))->Fill(mes,wfermi);
	       sprintf(name, "%s%s%s%s%s",prefix,lgroup,"leptonmc","eff",flavcat);
	       ((TH1D*)gDirectory->Get(name))->Fill(mes,wfermi);
	       if(mxhadfit<MXCUT){
		 sprintf(name, "%s%s%s%s",lgroup,"leptonmc","mxcut",flavcat);
		 ((TH1D*)gDirectory->Get(name))->Fill(mes,wfermi);
		 sprintf(name, "%s%s%s%s%s",prefix,lgroup,"leptonmc","mxcut",flavcat);
		 ((TH1D*)gDirectory->Get(name))->Fill(mes,wfermi);
	       }
	     }
	   }else{
	     if(vcb || iother){
	       //Weight for vcb is adedd
	       sprintf(name, "%s%s%s","all","leptonmc",flavcat);      
	       ((TH1D*)gDirectory->Get(name))->Fill(mes,w);
	       sprintf(name, "%s%s%s",lgroup,"leptonmc",flavcat);
	       ((TH1D*)gDirectory->Get(name))->Fill(mes,w);
 	       sprintf(name, "%s%s%s",prefix,lgroup,"leptonmc");
	       ((TH1D*)gDirectory->Get(name))->Fill(mes,w);
	     }
	   } 	  	
	 }
	 
       }
       
       // variables for the final selection
       
       //number of leptons cut 
       islept = ((nle == 1) && (tlab<2.37) && (tlab>0.36) && (plab>0.5));
       if(LEPTTYPE == 0) islept = (((nel == 1) && (nle == 1)) && (tlab<2.37) && (tlab>0.36) && (plab>0.5));
       if(LEPTTYPE == 1) islept = (((nmu == 1) && (nle == 1)) && (tlab<2.37) && (tlab>0.36) && (plab>0.5));
       
       int flav =  lcharge + brecoflav; // charge correlation
       bool ksele = (nkp + nks) > 0;    // fit on the depleted sample?
       int ch = TMath::Abs(xcharge + brecocharge);  // total charge
       double cutflav = 0;              // flavor category (3 = charged B, 4 = neutral B OS, 5 = neutral B SS)  
       int ibch = 1;                     
       
       if(TMath::Abs(brecocharge) && flav != 0) ibch=0;
       //       cout<<" event:: "<<jentry<<"2ONLY :: "<<nres<<" "<<Gvxbtyp<<endl;
       if(pcms> LEPTONPCUT && islept && ibch && mm2 < MNUSQHIGH && mm2 > MNUSQLOW &&  ch < CHHIGH && ch > CHLOW &&  ksele == DEPL && mxhadfit>0. && !(TMath::Abs(brecocharge)!=0 && (lcharge + brecoflav)!=0) && q2fit>=Q2CUT ){  

	 if(wdeltam>PRMM2CUT && brecocharge == 0) continue;

	 if(ISSMEARALL && isdata) mxhadfit = smeargauss(mxhadfit,SMEARALLMEANVALUE,SMEARALLSIGMA);    
	 if(ISSMEARBKG && !isdata) mxhadfit = smeargauss(mxhadfit,SMEARALLMEANVALUE,SMEARALLSIGMA); 
	 
	 if((ISSMEARALL && isdata) || (ISSMEARBKG && !isdata)) {
	   ((TH1D*)gDirectory->Get("h8888"))->Fill(mxhadfit - prevalue);
	   if(mxhadfit - prevalue == 0) cout << mxhadfit << endl;
	   if (mxhadfit<0) mxhadfit=0;
	 }	 
	 
	 // assign flavor category
	 cutflav = 5;
	 if(TMath::Abs(brecocharge)) cutflav = 3;
	 if(TMath::Abs(brecocharge)==0 && flav==0) cutflav = 4;	  
         
	 if(cutflav==3) sprintf(flavcat, "%s", "bch");
	 if(cutflav==4) sprintf(flavcat, "%s", "bos");
	 if(cutflav==5) sprintf(flavcat, "%s", "bss");
	 
	 id = hist(mxhadfit);    
	 if  (isdata != 1){

	   //Crosscheck on reweightings
	   //Tobedeleted
	   if(icat==0){ 
	     if(vub&&(Gvxbtyp!=7&&Gvxbtyp!=-7)&&mes>5.27&&intpur>.5){
	       ((TH1D*)gDirectory->Get("ACplotres"))->Fill(mxhadgen,wfermi);
	       ((TH1D*)gDirectory->Get("ACplotall"))->Fill(mxhadgen,wfermi);
	     }
	     if(vub&&!(Gvxbtyp!=7&&Gvxbtyp!=-7)&&mes>5.27&&intpur>.5){
	       ((TH1D*)gDirectory->Get("ACplotnres"))->Fill(mxhadgen,wfermi);
	       ((TH1D*)gDirectory->Get("ACplotall"))->Fill(mxhadgen,wfermi);
	       ((TH1D*)gDirectory->Get("ACplotallnonres"))->Fill(mxhadgen,wnre);
	     }
	   }

	   if(icat==0){ 
	     if(vub) {		
	       sprintf(name, "%s%s",lgroup,"allcutsmc");
	       if(mxhadfit<2.5) ((TH1D*)gDirectory->Get(name))->Fill(mes,wfermi); 
	       sprintf(name, "%s%s%s%s%d",lgroup,"allcutsmc",flavcat,"bin",id);
	       ((TH1D*)gDirectory->Get(name))->Fill(mes,wfermi); 
	       if(mxhadfit<MXCUT) {
		 if(DOVARSTU) fillMesStu(pcms,q2fit,EwPwfit,csiCiuc,xCiuc,wCiuc,lgroup,flavcat,mes,wfermi,"mc");
	       }
	       sprintf(name, "%s%s%s%s%s%d",prefix,lgroup,"allcutsmc",flavcat,"bin",id);
	       ((TH1D*)gDirectory->Get(name))->Fill(mes,wfermi);  
	       sprintf(name, "%s%s%s%s",lgroup,"allcutsmc","eff",flavcat);
	       ((TH1D*)gDirectory->Get(name))->Fill(mes,wfermi);
	       if (hist(mxhadfit)<5) {
		 sprintf(name, "%s%s%s%s",lgroup,"allcutsmc","1bin",flavcat);
		 ((TH1D*)gDirectory->Get(name))->Fill(mes,wfermi);
		 sprintf(name, "%s%s%s%s%s",prefix,lgroup,"allcutsmc","1bin",flavcat);
		 ((TH1D*)gDirectory->Get(name))->Fill(mes,wfermi);
	       }
	     }
	   }else{
	     if(vcb || iother){
	       sprintf(name, "%s%s",lgroup,"allcutsmc");
	       if(mxhadfit<2.5) ((TH1D*)gDirectory->Get(name))->Fill(mes,w); 
	       sprintf(name, "%s%s%s%s%d",lgroup,"allcutsmc",flavcat,"bin",id);
	       ((TH1D*)gDirectory->Get(name))->Fill(mes,w);  
	       if(mxhadfit<MXCUT) {
		 if(DOVARSTU) fillMesStu(pcms,q2fit,EwPwfit,csiCiuc,xCiuc,wCiuc,lgroup,flavcat,mes,w,"mc");
	       }
	       sprintf(name, "%s%s%s%s%s%d",prefix,lgroup,"allcutsmc",flavcat,"bin",id);
	       ((TH1D*)gDirectory->Get(name))->Fill(mes,w);  
	       if (hist(mxhadfit)<5) {
		 sprintf(name, "%s%s%s%s",lgroup,"allcutsmc","1bin",flavcat);
		 ((TH1D*)gDirectory->Get(name))->Fill(mes,w);
		 sprintf(name, "%s%s%s%s%s",prefix,lgroup,"allcutsmc","1bin",flavcat);
		 ((TH1D*)gDirectory->Get(name))->Fill(mes,w);
	       }
	     }
	   } 	  
	 }else{
	   if(isMC != 1) {
	     ((TH1D*)gDirectory->Get("allcutsdata2"))->Fill(mes); 	      
	     if(mxhadfit<2.5) ((TH1D*)gDirectory->Get("allcutsdata"))->Fill(mes); 
	     sprintf(name, "%s%s%s%d","allcutsdata",flavcat,"bin",id);
	     ((TH1D*)gDirectory->Get(name))->Fill(mes);  

	     if(mxhadfit<MXCUT) {
	       //	       cout<<"Values3:: "<<pcms<< " "<<q2fit<< " "<<EwPwfit<< " "<<csiCiuc<< " " << wCiuc<< " "<<xCiuc<<endl;
		 if(DOVARSTU) fillMesStu(pcms,q2fit,EwPwfit,csiCiuc,xCiuc,wCiuc,lgroup,flavcat,mes,1,"data");
	     }
	     sprintf(name, "%s%s%s%s%d",prefix,"allcutsdata",flavcat,"bin",id);
	     ((TH1D*)gDirectory->Get(name))->Fill(mes);  
	     if (hist(mxhadfit)<5) {
	       sprintf(name, "%s%s%s","allcutsdata","1bin",flavcat);
	       ((TH1D*)gDirectory->Get(name))->Fill(mes);
	       sprintf(name, "%s%s%s%s",prefix,"allcutsdata","1bin",flavcat);
	       ((TH1D*)gDirectory->Get(name))->Fill(mes);		
	       sprintf(truthgroup,"vcb");
	       if(vcb + vub == 0) sprintf(truthgroup,"oth");
	       if(vub) sprintf(truthgroup,"vub");
	       sprintf(name, "%s%s%s%s",truthgroup,"allcutsdata","1bintruth",flavcat);
	       ((TH1D*)gDirectory->Get(name))->Fill(mes);
	       sprintf(name, "%s%s%s%s%s",prefix,truthgroup,"allcutsdata","1bintruth",flavcat);
	       ((TH1D*)gDirectory->Get(name))->Fill(mes);
	     }

	   } else if (isMC == 1) {	      
	     if(icat==0){ 
	       if(vub) {		
		 ((TH1D*)gDirectory->Get("allcutsdata2"))->Fill(mes);  
		 if(mxhadfit<2.5) ((TH1D*)gDirectory->Get("allcutsdata"))->Fill(mes); 
		 sprintf(name, "%s%s%s%d","allcutsdata",flavcat,"bin",id);
		 ((TH1D*)gDirectory->Get(name))->Fill(mes);  
		 if(mxhadfit<MXCUT) {
		   //		   cout<<"Values4:: "<<pcms<< " "<<q2fit<< " "<<EwPwfit<< " "<<csiCiuc<< " " << wCiuc<< " "<<xCiuc<<endl;
		   if(DOVARSTU) fillMesStu(pcms,q2fit,EwPwfit,csiCiuc,xCiuc,wCiuc,lgroup,flavcat,mes,1,"data");
		 }
		 sprintf(name, "%s%s%s%s%d",prefix,"allcutsdata",flavcat,"bin",id);
		 ((TH1D*)gDirectory->Get(name))->Fill(mes);  
		 if (hist(mxhadfit)<5) {
		   sprintf(name, "%s%s%s","allcutsdata","1bin",flavcat);
		   ((TH1D*)gDirectory->Get(name))->Fill(mes);
		   sprintf(name, "%s%s%s%s",prefix,"allcutsdata","1bin",flavcat);
		   ((TH1D*)gDirectory->Get(name))->Fill(mes);		
		   sprintf(truthgroup,"vcb");
		   if(vcb + vub == 0) sprintf(truthgroup,"oth");
		   if(vub) sprintf(truthgroup,"vub");
		   sprintf(name, "%s%s%s%s",truthgroup,"allcutsdata","1bintruth",flavcat);
		   ((TH1D*)gDirectory->Get(name))->Fill(mes);
		   sprintf(name, "%s%s%s%s%s",prefix,truthgroup,"allcutsdata","1bintruth",flavcat);
		   ((TH1D*)gDirectory->Get(name))->Fill(mes);
		 }
	       }	        
	     }else{
	       if(vcb || iother){
		 ((TH1D*)gDirectory->Get("allcutsdata2"))->Fill(mes);  
		 if(mxhadfit<2.5) ((TH1D*)gDirectory->Get("allcutsdata"))->Fill(mes); 
		 sprintf(name, "%s%s%s%d","allcutsdata",flavcat,"bin",id);
		 ((TH1D*)gDirectory->Get(name))->Fill(mes);  
		 if(mxhadfit<MXCUT) {
		   if(DOVARSTU) fillMesStu(pcms,q2fit,EwPwfit,csiCiuc,xCiuc,wCiuc,lgroup,flavcat,mes,1,"data");
		 } 
		 sprintf(name, "%s%s%s%s%d",prefix,"allcutsdata",flavcat,"bin",id);
		 ((TH1D*)gDirectory->Get(name))->Fill(mes);  
		 if (hist(mxhadfit)<5) {
		   sprintf(name, "%s%s%s","allcutsdata","1bin",flavcat);
		   ((TH1D*)gDirectory->Get(name))->Fill(mes);
		   sprintf(name, "%s%s%s%s",prefix,"allcutsdata","1bin",flavcat);
		   ((TH1D*)gDirectory->Get(name))->Fill(mes);		
		   sprintf(truthgroup,"vcb");
		   if(vcb + vub == 0) sprintf(truthgroup,"oth");
		   if(vub) sprintf(truthgroup,"vub");
		   sprintf(name, "%s%s%s%s",truthgroup,"allcutsdata","1bintruth",flavcat);
		   ((TH1D*)gDirectory->Get(name))->Fill(mes);
		   sprintf(name, "%s%s%s%s%s",prefix,truthgroup,"allcutsdata","1bintruth",flavcat);
		   ((TH1D*)gDirectory->Get(name))->Fill(mes);
		 }
	       }		
	     }
	   }
	 }
       }
     }
   }  
}




Double_t ftotal(Double_t *x, Double_t *par) {
  
  // histo model to fit the Mx distribution 
  // par[0] = fitted integral of Vcb component
  // par[1] = fitted integral of other component 
  // par[2] = integral of MC Vcb histogram used as model
  // par[3] = integral of MC Vub histogram used as model 
  // par[4] = charged multipl. category 
  // par[5] = neutral multipl. category 
  // par[6] = switch to use multiplicity fit
  // par[7] = ratio of BR to determine Vub tail at Mx > Mx
  // par[8] = option to fit 3 components: vub, vcb, other
  // par[9] = fitted integral of Vub component

  Double_t xx = x[0];
  char name2[100], preine2[100], preich2[100], prefix2[100];
  int ich2 = par[4];  // charged multipl. categories
  int ine2 = par[5];  // neutral multipl. categories
  sprintf(preich2, "%s%d", "ch", ich2);
  sprintf(preine2, "%s%d", "ne", ine2);
  if (par[6]) { sprintf(prefix2, "%s%s",preich2 , preine2);} // use multiplicities? (0 = no, 1 = yes)
  else{sprintf(prefix2, "");}               
  sprintf(name2, "%s%s",prefix2,"mxdata");     
  if(par[8]) sprintf(name2, "%s%s",prefix2,"mxonebdata");

  Int_t bin = ((TH1D*)gDirectory->Get(name2))->GetXaxis()->FindBin(xx);
  double thereturn;
  if(par[8]) {
    sprintf(name2, "%s%s",prefix2,"mxonebvub");   
    Double_t vub = par[9] * ((TH1D*)gDirectory->Get(name2))->GetBinContent(bin);
    sprintf(name2, "%s%s",prefix2,"mxonebvcb");   
    Double_t vcb = par[0] * ((TH1D*)gDirectory->Get(name2))->GetBinContent(bin);
    sprintf(name2, "%s%s",prefix2,"mxoneboth");    
    Double_t oth = par[1] * ((TH1D*)gDirectory->Get(name2))->GetBinContent(bin);
    thereturn = vub + vcb + oth;    
  }else{
    sprintf(name2, "%s%s",prefix2,"mxvub");   
    Double_t vub = par[0] * ((TH1D*)gDirectory->Get(name2))->GetBinContent(bin);
    sprintf(name2, "%s%s",prefix2,"mxvcb");   
    Double_t vcb = par[0] * ((TH1D*)gDirectory->Get(name2))->GetBinContent(bin);
    sprintf(name2, "%s%s",prefix2,"mxoth");    
    Double_t oth = par[1] * ((TH1D*)gDirectory->Get(name2))->GetBinContent(bin);
    thereturn = vub * (par[7]/.104) * (par[2] / par[3]) + vcb + oth;
  }      
  return  thereturn;
}
 
TH1D* fitNtp::makeMesHist(const char * name, bool doToy){    
  // if TOY is requested each Mes fit is randomized in a poisson way in each bin
  TH1D* origHist=(TH1D*)gDirectory->Get(name);
  TH1D* newHist=new TH1D(*((TH1D*)gDirectory->Get(name)));
  if(doToy) {
    TRandom random(0);
    
    for(int k=1;k<40;k++){
      int newBinContent=random.Poisson(origHist->GetBinContent(k));
      newHist->SetBinContent(k,newBinContent);
      newHist->SetBinError(k,sqrt(newBinContent));
    }
  }
  return newHist;
}


void chi2Hist(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
    int ioth=TMath::Abs(par[3])>2 ? 1 : 2 ;
    Double_t chi2(0.);
    int nBin = vubHist->GetNbinsX();
    for (int j=0;j<nBin;j++){
      double x=par[0]*vubHist->GetBinContent(j)+par[1]*vcbHist->GetBinContent(j)+par[ioth]*othHist->GetBinContent(j);
      if(x>0 && dataHist->GetBinContent(j)>0){
	double s2=pow(dataHist->GetBinError(j),2);
	if(par[3]<0.1) s2+= // add MC stat
	 (pow(par[0]*vubHist->GetBinError(j),2)+
	  pow(par[1]*vcbHist->GetBinError(j),2)+
	  pow(par[1]*othHist->GetBinError(j),2));

	chi2+=pow((x-dataHist->GetBinContent(j)),2)/s2;
      }
    }
    f=chi2;
}

TH1D* fitNtp::fitWithErrors(int ich,int ine,char * typ){

  TMinuit aMinuit(4);
  int ierflg;
  // initialize minuit
  aMinuit.SetFCN(chi2Hist);
  aMinuit.mnparm(0, "vubcomp",  1.0, 0.001, 0., 10., ierflg); 
  aMinuit.mnparm(1, "vcbcomp",  1.0, 0.001, 0., 10., ierflg); 
  aMinuit.mnparm(2, "othcomp",  1.0, 0.001, 0., 10., ierflg); 
  aMinuit.mnparm(3, "flag",  0., 0., 0., 100., ierflg);
  if(FITOPT==1){
    aMinuit.mnparm(3, "flag",  -2., 0., -100., 100., ierflg);
    aMinuit.FixParameter(2); 
  }
  if(FITOPT==3){
    //VALUE TAKEN FROM DEFAULT FIT (Mx<1.75, no Q2 CUT)
    aMinuit.mnparm(1, "vcbcomp", 0.318434, 0.001, 0., 10., ierflg); 
    aMinuit.FixParameter(1); 
    aMinuit.mnparm(2, "othcomp", 0.0747699, 0.001, 0., 10., ierflg); 
    aMinuit.FixParameter(2); 
  }
  if(FITOPT==4){
    //VALUE TAKEN FROM DEFAULT FIT (Mx<1.55, no Q2 CUT)
    aMinuit.mnparm(1, "vcbcomp", 0.291005, 0.001, 0., 10., ierflg); 
    aMinuit.FixParameter(1); 
    aMinuit.mnparm(2, "othcomp", 0.0632736, 0.001, 0., 10., ierflg); 
    aMinuit.FixParameter(2); 
  }
  if(FITOPT==5){
    cout<<VCBCOMP<<" "<<OTHCOMP<<endl;
    aMinuit.mnparm(1, "vcbcomp", VCBCOMP, 0.001, 0., 10., ierflg); 
    aMinuit.FixParameter(1); 
    aMinuit.mnparm(2, "othcomp", OTHCOMP, 0.001, 0., 10., ierflg); 
    aMinuit.FixParameter(2); 
  }
  aMinuit.FixParameter(3); 
  // get histograms
  char name[200],name2[100],  preine2[100], preich2[100] , prefix2[100];
  if(ich>=0 && ine>=0){
    sprintf(preich2, "%s%d", "ch", ich);
    sprintf(preine2, "%s%d", "ne", ine);
    sprintf(prefix2, "%s%s",preich2 , preine2);} // use multiplicities? (0 = no, 1 = yes)
  else{sprintf(prefix2, "");}               

  sprintf(name, "%s%s%s",prefix2,typ,"data");
  
  if(FITTOTSHAPE) {
    sprintf(name2, "%s%s%s",prefix2,typ,"onebvub");   
    vubHist = (TH1D*)gDirectory->Get(name2);
    sprintf(name2, "%s%s%s",prefix2,typ,"onebvcb");   
    vcbHist = (TH1D*)gDirectory->Get(name2);
    sprintf(name2, "%s%s%s",prefix2,typ,"oneboth");    
    othHist = (TH1D*)gDirectory->Get(name2);
    sprintf(name2, "%s%s%s",prefix2,typ,"onebdata");    
    dataHist = (TH1D*)gDirectory->Get(name2);
    
  }else{
    sprintf(name2, "%s%s%s",prefix2,typ,"vub");   
    vubHist = (TH1D*)gDirectory->Get(name2);
    sprintf(name2, "%s%s%s",prefix2,typ,"vcb");   
    vcbHist = (TH1D*)gDirectory->Get(name2);
    sprintf(name2, "%s%s%s",prefix2,typ,"oth");    
    othHist = (TH1D*)gDirectory->Get(name2);
    sprintf(name2, "%s%s%s",prefix2,typ,"data");    
    dataHist = (TH1D*)gDirectory->Get(name2);
    vubHist->Scale(BRRATIOVALUETAIL/0.104*vcbmc/vcbmc);	
    
  } 
  int ioth=FITOPT==1 ? 1 : 2 ;
  // fit
  Double_t arglis[10];
  arglis[0] = 10000; // maxcalls
  arglis[1] = 0.001;  // tolerance
  aMinuit.mnexcm("MINI", arglis, 2, ierflg);
  aMinuit.mnexcm("MINOS", arglis, 2, ierflg);
  if(ierflg!=0){
    cout <<" fit with MC truth: minos error "<<ierflg<<endl;
  }
  aMinuit.GetParameter(0,vubcomp,errvubcomp);
  aMinuit.GetParameter(1,vcbcomp,errvcbcomp);  
  aMinuit.GetParameter(ioth,othcomp,errothcomp);
// no mc truth 
  aMinuit.mnparm(3, "flag",  1., 0., 0., 100., ierflg);
  aMinuit.FixParameter(3); 
  if(FITOPT==1){
    aMinuit.mnparm(3, "flag",  2., 0., 0., 100., ierflg);
    aMinuit.FixParameter(2); 
  }
  if(FITOPT==3){
    //VALUE TAKEN FROM DEFAULT FIT (Mx<1.75, no Q2 CUT)
    aMinuit.mnparm(1, "vcbcomp", 0.318434, 0.001, 0., 10., ierflg); 
    aMinuit.FixParameter(1); 
    aMinuit.mnparm(2, "othcomp", 0.0747699, 0.001, 0., 10., ierflg); 
    aMinuit.FixParameter(2); 
  }
  if(FITOPT==4){
    //VALUE TAKEN FROM DEFAULT FIT (Mx<1.55, no Q2 CUT)
    aMinuit.mnparm(1, "vcbcomp", 0.291005, 0.001, 0., 10., ierflg); 
    aMinuit.FixParameter(1); 
    aMinuit.mnparm(2, "othcomp", 0.0632736, 0.001, 0., 10., ierflg); 
    aMinuit.FixParameter(2); 
  }
  if(FITOPT==5){
    cout<<VCBCOMP<<" "<<OTHCOMP<<endl;
    aMinuit.mnparm(1, "vcbcomp", VCBCOMP, 0.001, 0., 10., ierflg); 
    aMinuit.FixParameter(1); 
    aMinuit.mnparm(2, "othcomp", OTHCOMP, 0.001, 0., 10., ierflg); 
    aMinuit.FixParameter(2); 
  }
  aMinuit.mnexcm("MINI", arglis, 2, ierflg);
  if(ierflg!=0){
    cout <<" fit without MC truth: minos error "<<ierflg<<endl;
  }
  aMinuit.mnexcm("MINOS", arglis, 2, ierflg);
  aMinuit.GetParameter(0,vubcompNOMC,errvubcompNOMC);
  aMinuit.GetParameter(1,vcbcompNOMC,errvcbcompNOMC);
  aMinuit.GetParameter(ioth,othcompNOMC,errothcompNOMC);
  if(vubcompNOMC>vubcomp || vubcompNOMC>vubcomp || vubcompNOMC>vubcomp ){
    cout <<" removing MC truth improves error!!! "<<" vub "<<vubcompNOMC<<" "<<vubcomp
	 <<" vcb "<<vcbcompNOMC<<" "<<vcbcomp<<" oth "<<othcompNOMC<<" "<<othcomp<<endl;
    //    assert(0);
  }
  int nBin=vubHist->GetNbinsX();

  TH1D* result = new TH1D(*vubHist);
  for (int j=0;j<nBin;j++){
 	double x=vubcomp*vubHist->GetBinContent(j)+vcbcomp*vcbHist->GetBinContent(j)+othcomp*othHist->GetBinContent(j);
 	double s=sqrt(
		      pow(vubcomp*vubHist->GetBinError(j),2)+
		      pow(vcbcomp*vcbHist->GetBinContent(j),2)+
		      pow(othcomp*othHist->GetBinContent(j),2)
		      );
	result->SetBinContent(j,x); result->SetBinError(j,s);
  }

  TCanvas *cla1 = new TCanvas(); 
  TGraph *g1 = scanParameter(0, 2, aMinuit, chi2Hist); 
  g1->SetTitle("scan b2u"); 
  g1->SetMarkerSize(0.5); 
  g1->Draw("ap"); 
  cla1->SaveAs(Form("%s%s%s",DIRNAME.Data(),PREFIXOUT.Data(),"vub_scan.eps"));

  cla1->Clear(); 
  TGraph *g2 = scanParameter(1, 2, aMinuit, chi2Hist); 
  g2->SetTitle("scan b2c"); 
  g2->SetMarkerSize(0.5); 
  g2->Draw("ap"); 
  cla1->SaveAs(Form("%s%s%s",DIRNAME.Data(),PREFIXOUT.Data(),"vcb_scan.eps"));

  cla1->Clear(); 
  TGraph *g3 = scanParameter(2, 2, aMinuit, chi2Hist); 
  g3->SetTitle("scan other"); 
  g3->SetMarkerSize(0.5); 
  g3->Draw("ap"); 
  cla1->SaveAs(Form("%s%s%s",DIRNAME.Data(),PREFIXOUT.Data(),"oth_scan.eps"));
  cout<<"********** :: "<<vcbcomp<<" "<<othcomp<<endl;
  delete g1; 
  delete g2; 
  delete g3; 
  delete cla1; 
  
  return result;   
}


    
// ----------------------------------------------------------------------
TGraph* fitNtp::scanParameter(int parnum, int nsig, TMinuit &a, void (*func)(int &, double *, double &, double *, int)) {

#define MAXPAR 100
#define NSTEP  400

  double par[MAXPAR], parE[MAXPAR];

  int npar = a.GetNumPars(); 

  for (int ipar = 0; ipar < npar; ++ipar) {
    a.GetParameter(ipar, par[ipar], parE[ipar]); 
  }

  int iflag; 
  double x[NSTEP], y[NSTEP], dummy[NSTEP]; 
  char ctitle[100]; 

  double min  = par[parnum] - nsig*parE[parnum]; 
  double max  = par[parnum] + nsig*parE[parnum]; 
  double step = (max - min)/NSTEP; 
  double val(0.); 

  for (int ix = 0; ix < NSTEP; ++ix) {
    par[parnum] = min + ix*step; 
    (*func)(iflag, dummy, val, par, iflag); 
    x[ix] = par[parnum]; 
    y[ix] = val; 
  }
    
  TGraph *tg = new TGraph(NSTEP, x, y); 
  tg->Draw("ap"); 

  return tg; 
}



TH1D* fitNtp::fitWithoutErrors(int ich ,int ine){
  double highvalue = 5.;
  Double_t lowvalue = MXCUT;
  if(FITTOTSHAPE) lowvalue=0;

  TF1 *ftot = new TF1("ftot", ftotal, lowvalue, highvalue, 10);
  int mult= (ich>=0 && ine>=0)? 1 : 0;
  ftot->SetParameters(1,1,vcbmc,vubmc,ich,ine,mult,BRRATIOVALUETAIL,FITTOTSHAPE,1);
  ftot->FixParameter(2, vcbmc);
  ftot->FixParameter(3, vubmc);
  ftot->FixParameter(4, ich);
  ftot->FixParameter(5, ine);
  ftot->FixParameter(6, mult);         
  ftot->FixParameter(7, BRRATIOVALUETAIL);   
  ftot->FixParameter(8, FITTOTSHAPE);
  //Need to fix vcb component
  ftot->SetParLimits(0, 0.0, 3.0);
  ftot->SetParLimits(1, 0.0, 3.0);
  ftot->SetParLimits(9, 0.0, 3.0);
  if(FITOPT==5){
    cout<<VCBCOMP<<" "<<OTHCOMP<<endl;
    ftot->FixParameter(0, VCBCOMP);
    ftot->FixParameter(1, OTHCOMP);
  }

  if(!(FITTOTSHAPE)) ftot->FixParameter(9, 0);
  char name[200],name2[100],  preine2[100], preich2[100] , prefix2[100];
  if(ich>=0 && ine>=0){
    sprintf(preich2, "%s%d", "ch", ich);
    sprintf(preine2, "%s%d", "ne", ine);
    sprintf(prefix2, "%s%s",preich2 , preine2);} // use multiplicities? (0 = no, 1 = yes)
  else{sprintf(prefix2, "");}               
  sprintf(name, "%s%s",prefix2,"mxdata");
	         
  if(FITTOTSHAPE) sprintf(name, "%s%s",prefix2,"mxonebdata");	         
  if(BLINDING){
    ((TH1D*)gDirectory->Get(name))->Fit("ftot","bcemqn"); 
  }else{
    ((TH1D*)gDirectory->Get(name))->Fit("ftot","bcemn"); 
  }
  
  //component and their errors from the fit
  
  vubcomp = ftot->GetParameter(9);
  errvubcomp = ftot->GetParError(9);
  vcbcomp = ftot->GetParameter(0);
  errvcbcomp = ftot->GetParError(0);
  othcomp = ftot->GetParameter(1);
  errothcomp = ftot->GetParError(1);
  
  
  if(FITTOTSHAPE) {
    sprintf(name2, "%s%s",prefix2,"mxonebvub");   
    vubHist = (TH1D*)gDirectory->Get(name2);
    sprintf(name2, "%s%s",prefix2,"mxonebvcb");   
    vcbHist = (TH1D*)gDirectory->Get(name2);
    sprintf(name2, "%s%s",prefix2,"mxoneboth");    
    othHist = (TH1D*)gDirectory->Get(name2);

  }else{
    sprintf(name2, "%s%s",prefix2,"mxvub");   
    vubHist = (TH1D*)gDirectory->Get(name2);
    sprintf(name2, "%s%s",prefix2,"mxvcb");   
    vcbHist = (TH1D*)gDirectory->Get(name2);
    sprintf(name2, "%s%s",prefix2,"mxoth");    
    othHist = (TH1D*)gDirectory->Get(name2);
    vubHist->Scale(ftot->GetParameter(7)/0.104*ftot->GetParameter(2)/ftot->GetParameter(3)/vubcomp);	
    
  } 
  int nBin=vubHist->GetNbinsX();
  TH1D* result = new TH1D(*vubHist);
  for (int j=0;j<nBin;j++){
 	double x=vubcomp*vubHist->GetBinContent(j)+vcbcomp*vcbHist->GetBinContent(j)+othcomp*othHist->GetBinContent(j);
 	double s=sqrt(
		      pow(vubcomp*vubHist->GetBinError(j),2)+
		      pow(vcbcomp*vcbHist->GetBinError(j),2)+
		      pow(othcomp*othHist->GetBinError(j),2)
		      );
	result->SetBinContent(j,x); result->SetBinError(j,s);
  }
  delete ftot;
  return result;   
}

void fitNtp::theFit(){    


   fHistFile->cd();
   char name[200], name2[100], preine[100], preich[100], prefix[100];
   char namebch[100], nameb0os[100], nameb0ss[100];
   char nameps[100];  
   char prefixps[100];
   double vubafter[3][2] = {{0,0} , {0,0} , {0,0}};
   double errvubafter[3][2] = {{0,0} , {0,0} , {0,0}};
   double errvubafterMCstat[3][2] = {{0,0} , {0,0} , {0,0}};
   double errchisq;
   double tempbin = 0;
   double temperr = 0;
   double tempbinvub = 0;
   double tempbinvcb = 0;
   double tempbinoth = 0;
   double temperrvcb = 0;
   double temperroth = 0;
   double tempbinerrvub = 0;
   double tempbinerrvcb = 0;
   double tempbinerroth = 0;
   double vubaftercutsbin1 = 0;
   double errvubaftercutsbin1 = 0;
   double MCerrvubaftercutsbin1 = 0;
   areavcbaftercutsbin1 = 0;
   vcbaftercutsbin1 = 0;
   errvcbaftercutsbin1 = 0;
   errfitvcbaftercutsbin1 = 0;
   otheraftercutsbin1 = 0;
   errotheraftercutsbin1 = 0;
   errfitotheraftercutsbin1 = 0;
   double tempbinchb;
   double tempbinb0os;
   double tempbinb0ss;
   double temperrchb;
   double temperrb0os;
   double temperrb0ss;
   double chid = 0.181;
   double fact;
   double mxcut = MXCUT;
   Int_t ich;
   Int_t ine;
   double B;
   double Bvcb;
   double Bother;
   double errB;
   double errBvcb;
   double errBother;
   double vubmcselected2;
   double vcbmcselected2;
   double othermcselected2;
   double vubmccut;
   double vubmc2cut;
   double vubmc2aftercutsbin1;
   double vubmcorig;
   double vcbmcorig;
   double othermcorig;
   double foomean;
   double foosigma;
   double fooalpha;
   double foon;
   double min, max;
   blindfactor = getblindfact();
   
   for (int is=1;is<11;is++) {
     if(is==1)  {sprintf(name, "vcbleptonmc");}
     if(is==2)  {sprintf(name, "vubleptonmc");}
     if(is==3)  {sprintf(name, "allleptonmc");}
     if(is==4)  {sprintf(name, "leptondata");}
     if(is==5)  {sprintf(name, "vubleptonmcmxcut");}
     if(is==6)  {sprintf(name, "vcballcutsdata1bintruth");}
     if(is==7)  {sprintf(name, "vuballcutsdata1bintruth");}
     if(is==8)  {sprintf(name, "othallcutsdata1bintruth");}
     if(is==9)  {sprintf(name, "vubleptonmceff");}
     if(is==10)  {sprintf(name, "vuballcutsmceff");}
     sprintf(namebch, "%s%s",name,"bch");  
     sprintf(nameb0os, "%s%s",name,"bos");  
     sprintf(nameb0ss, "%s%s",name,"bss");       
     
  // Correction to have the right Bch/B0 ratio (as data)
     double correctionratio = 1;
     if(BTYPE == 2){
       if(is==1 || is==3) {
	 correctionratio = correctionratiovcb;
       }
       if(is == 2 || is ==5 || is == 9 || is == 10) {
	 correctionratio = correctionratiovub;
       }
     }

     for(int k=1;k<40;k++){
       tempbinchb = ((TH1D*)gDirectory->Get(namebch))->GetBinContent(k) * correctionratio;
       tempbinb0os = ((TH1D*)gDirectory->Get(nameb0os))->GetBinContent(k);
       tempbinb0ss = ((TH1D*)gDirectory->Get(nameb0ss))->GetBinContent(k);
       temperrchb = ((TH1D*)gDirectory->Get(namebch))->GetBinError(k) * correctionratio;
       temperrb0os = ((TH1D*)gDirectory->Get(nameb0os))->GetBinError(k);
       temperrb0ss = ((TH1D*)gDirectory->Get(nameb0ss))->GetBinError(k);
       tempbin = tempbinchb + ((1-chid)/(1-2*chid)) * tempbinb0os - (chid/(1-2*chid)) * tempbinb0ss;
       if(MIXCORR==0){
	 tempbin = tempbinchb + tempbinb0os;
	 temperr = sqrt(temperrchb * temperrchb + temperrb0os * temperrb0os);
       }else{
	 if(tempbin<0) tempbin=0;
	 temperr = sqrt(temperrchb*temperrchb + ((1-chid)/(1-2*chid)) * ((1-chid)/(1-2*chid)) * temperrb0os * temperrb0os + (chid/(1-2*chid)) * (chid/(1-2*chid)) * temperrb0ss * temperrb0ss);
       }
       ((TH1D*)gDirectory->Get(name))->SetBinContent(k, tempbin);
       ((TH1D*)gDirectory->Get(name))->SetBinError(k, temperr); 
     }	     
   }
   
   if(MULTIFIT){
     
     for (ich=1;ich<4;ich++){
       for (ine=1;ine<3;ine++){      
	 sprintf(preich, "%s%d", "ch", ich);
	 sprintf(preine, "%s%d", "ne", ine);
	 sprintf(prefix, "%s%s",preich , preine);     
	 for (int is=1;is<7;is++) {
	   if(is==1)	 {sprintf(name, "%s%s",prefix, "vubleptonmcmxcut");}
	   if(is==2)	 {sprintf(name, "%s%s",prefix, "vubleptonmceff");}
	   if(is==3)	 {sprintf(name, "%s%s",prefix, "vuballcutsmc1bin");}
	   if(is==4)	 {sprintf(name, "%s%s",prefix, "vcballcutsdata1bintruth");}
	   if(is==5)	 {sprintf(name, "%s%s",prefix, "vuballcutsdata1bintruth");}
	   if(is==6)	 {sprintf(name, "%s%s",prefix, "othallcutsdata1bintruth");}
	   sprintf(namebch, "%s%s",name,"bch");  
	   sprintf(nameb0os, "%s%s",name,"bos");  
	   sprintf(nameb0ss, "%s%s",name,"bss");       
	   cout << name << endl;
	   for(int k=1;k<40;k++){
	     tempbinchb = ((TH1D*)gDirectory->Get(namebch))->GetBinContent(k);
	     tempbinb0os = ((TH1D*)gDirectory->Get(nameb0os))->GetBinContent(k);
	     tempbinb0ss = ((TH1D*)gDirectory->Get(nameb0ss))->GetBinContent(k);
	     temperrchb = ((TH1D*)gDirectory->Get(namebch))->GetBinError(k);
	     temperrb0os = ((TH1D*)gDirectory->Get(nameb0os))->GetBinError(k);
	     temperrb0ss = ((TH1D*)gDirectory->Get(nameb0ss))->GetBinError(k);
	     tempbin = tempbinchb + ((1-chid)/(1-2*chid)) * tempbinb0os - (chid/(1-2*chid)) * tempbinb0ss;
	     if(MIXCORR==0){
	       tempbin = tempbinchb + tempbinb0os;
	       temperr = sqrt(temperrchb * temperrchb + temperrb0os * temperrb0os);
	     }else{
	       if(tempbin<0) tempbin=0;
	       temperr = sqrt(temperrchb*temperrchb + ((1-chid)/(1-2*chid)) * ((1-chid)/(1-2*chid)) * temperrb0os * temperrb0os + (chid/(1-2*chid)) * (chid/(1-2*chid)) * temperrb0ss * temperrb0ss);
	     }
	     ((TH1D*)gDirectory->Get(name))->SetBinContent(k, tempbin);
	     ((TH1D*)gDirectory->Get(name))->SetBinError(k, temperr); 
	   }	
	 }
       }
     }
   }   
   
   c1 = new TCanvas("c1"," ",200,10,1200,1000); 
   c1->Clear();
   c1->Divide(4, 4);

   sprintf(prefixps, "%s", "mes_misc");

   sprintf(nameps, "%s%s%s%s",DIRNAME.Data(),PREFIXOUT.Data(),prefixps,".eps");  

   openEpsFile(nameps);

   c1->cd(1);

   //true number of vcb MC  
   TVector2 signal = sighisto(((TH1D*)gDirectory->Get("vcbmcnocutdata")),foomean,foosigma,fooalpha,foon,1,-1111111.,-1111111., -1111111.,5.,-1111111.);
   vcbmcorig = signal.X();

   c1->cd(2);

   //true number of vub MC     
   signal = sighisto(((TH1D*)gDirectory->Get("vubmcnocutdata")),foomean,foosigma,fooalpha,foon,1,-1111111.,-1111111.,-1111111., 5.,-1111111.);
   vubmcorig = signal.X();

   c1->cd(3);

   //true number of other MC     
   signal = sighisto(((TH1D*)gDirectory->Get("othmcnocutdata")),foomean,foosigma,fooalpha,foon,1,-1111111.,-1111111.,-1111111., 5.,-1111111.);
   othermcorig = signal.X();

   c1->cd(4);

   //true number of vcb MC events after lepton cut    
   signal = sighisto(((TH1D*)gDirectory->Get("vcbleptonmc")),foomean,foosigma,fooalpha,foon,1,mesvcbMC[0],mesvcbMC[1],mesvcbMC[2], 5.,-1111111.);
   vcbmc = signal.X();
   errvcbmc = signal.Y();

   c1->cd(5);

   //true number of vub MC events after lepton cut    
   signal = sighisto(((TH1D*)gDirectory->Get("vubleptonmc")),foomean,foosigma,fooalpha,foon,1,mesvubMC[0],mesvubMC[1],mesvubMC[2], 5.,-1111111.);
   vubmc = signal.X();
   errvubmc = signal.Y();

   c1->cd(6);

   //true number of vcb MC events after all cuts (crooscheck to find biases)    
   signal = sighisto(((TH1D*)gDirectory->Get("vcballcutsdata1bintruth")),foomean,foosigma,fooalpha,foon,1,-1111111.,-1111111.,-1111111., 5.,-1111111.);
   vcbmcselected = signal.X();

   c1->cd(7);

   //true number of vub MC events after all cuts (crooscheck to find biases)    
   signal = sighisto(((TH1D*)gDirectory->Get("vuballcutsdata1bintruth")),foomean,foosigma,fooalpha,foon,1,-1111111.,-1111111.,-1111111.,5.,-1111111.);
   vubmcselected = signal.X();

   c1->cd(8);

   //true number of other MC events after all cuts (crooscheck to find biases)    
   signal = sighisto(((TH1D*)gDirectory->Get("othallcutsdata1bintruth")),foomean,foosigma,fooalpha,foon,1,-1111111.,-1111111.,-1111111., 5.,-1111111.);
   othermcselected = signal.X();

   c1->cd(9);

   //true number of vub MC events after lepton cut for Mx<MXCUT  
   signal = sighisto(((TH1D*)gDirectory->Get("vubleptonmcmxcut")),foomean,foosigma,fooalpha,foon,1,mesvubMCmx[0],mesvubMCmx[1],mesvubMCmx[2], 5.,-1111111.);
   vubmccut = signal.X();

   c1->cd(10);

   //true number of vcb+other MC events after lepton cut (no vub event, correction after)    
   signal = sighisto(((TH1D*)gDirectory->Get("allleptonmc")),foomean,foosigma,fooalpha,foon,1,mesNslMC[0],mesNslMC[1],mesNslMC[2], 5.,-1111111.);
   totmc = signal.X();

   c1->cd(11);

   //total number (data) of events after lepton cut    
   signal = sighisto(((TH1D*)gDirectory->Get("leptondata")),foomean,foosigma,fooalpha,foon,1,mesNsl[0],mesNsl[1],mesNsl[2],mesNsl[3],-1111111);
   tot = signal.X();
   errtot = signal.Y();

   c1->cd(12);

   //true number of vub MC events after lepton cut for efficiency calculation (only vub-generic)    
   signal = sighisto(((TH1D*)gDirectory->Get("vubleptonmceff")),foomean,foosigma,fooalpha,foon,1,mesvubMClepteff[0],mesvubMClepteff[1],mesvubMClepteff[2], 5.,-1111111.);
   vubmcleptforeff = signal.X();

   c1->cd(13);

   //true number of vub MC events after all cuts for efficiency calculation (only vub-generic)    
   signal = sighisto(((TH1D*)gDirectory->Get("vuballcutsmceff")),foomean,foosigma,fooalpha,foon,1,mesvubMCalleff[0],mesvubMCalleff[1],mesvubMCalleff[2], 5.,-1111111.);
   vubmcallforeff = signal.X();

   c1->cd(14);

   //true number of vub MC events no cut
   signal = sighisto(((TH1D*)gDirectory->Get("vubmcnocut")),foomean,foosigma,fooalpha,foon,1,mesvubMCall[0],mesvubMCall[1],mesvubMCall[2], 5.,-1111111.);
   vubmcnocut = signal.X();

   c1->cd(15);

   //true number of vcb MC events no cut
   signal = sighisto(((TH1D*)gDirectory->Get("vcbmcnocut")),foomean,foosigma,fooalpha,foon,1,mesvcbMCall[0],mesvcbMCall[1],mesvcbMCall[2], 5.,-1111111.);
   vcbmcnocut = signal.X();

   c1->cd(16);

   //true number of oth MC events no cut
   signal = sighisto(((TH1D*)gDirectory->Get("othmcnocut")),foomean,foosigma,fooalpha,foon,1,-1111111.,-1111111.,-1111111., 5.,-1111111.);
   othmcnocut = signal.X();
   closeEpsFile();

   c1 = new TCanvas("c1"," ",200,10,1200,800); 
   c1->Clear();
   c1->Divide(1);
   sprintf(nameps, "%s%s%s",DIRNAME.Data(),PREFIXOUT.Data(),"nsl.eps");  
   openEpsFile(nameps);
   //total number (data) of events after lepton cut    
   ((TH1D*)gDirectory->Get("leptondata"))->SetTitle("NSL: lepton cuts");
   signal = sighisto(((TH1D*)gDirectory->Get("leptondata")),foomean,foosigma,fooalpha,foon,1,-1111111.,-1111111.,-1111111., -1111111.,-1111111);
   closeEpsFile();

   calcpstarfact  = (vubmc/vubmcnocut)/(vcbmc/vcbmcnocut);
   errcalcpstarfact = TMath::Sqrt((errvubmc/vubmc)*(errvubmc/vubmc) +  (errvcbmc/vcbmc)*(errvcbmc/vcbmc)); // approximately
   calcpstarfacttemp = calcpstarfact; 
   if(PSTARFACT>0) calcpstarfact = getPstarFactor(LEPTONPCUT);

   if(MULTIFIT){
     
     for (ich=1;ich<4;ich++){
       for (ine=1;ine<3;ine++){    
	 
	 sprintf(preich, "%s%d", "ch", ich);
	 sprintf(preine, "%s%d", "ne", ine);
	 sprintf(prefix, "%s%s",preich , preine);     
	 
  	 //true number of vub MC events after lepton cut    
 	 sprintf(name, "%s%s",prefix,"vubleptonmceff");	         
 	 signal = sighisto(((TH1D*)gDirectory->Get(name)),foomean,foosigma,fooalpha,foon,1,-1111111.,-1111111.,-1111111.,-1111111.,-1111111.);
 	 vubmc2 = signal.X();

 	 //true number of vub MC events after lepton cut for Mx<MXCUT  
 	 sprintf(name, "%s%s",prefix,"vubleptonmcmxcut");	         
 	 signal = sighisto(((TH1D*)gDirectory->Get(name)),foomean,foosigma,fooalpha,foon,1,-1111111.,-1111111.,-1111111.,-1111111.,-1111111.);
 	 vubmc2cut = signal.X();

 	 //true number of vub MC events after all cuts for Mx<MXCUT  
 	 sprintf(name, "%s%s",prefix,"vuballcutsmc1bin");
 	 signal = sighisto(((TH1D*)gDirectory->Get(name)),foomean,foosigma,fooalpha,foon,1,-1111111.,-1111111.,-1111111.,-1111111.,-1111111.);
 	 vubmc2aftercutsbin1 = signal.X();

  	 //true number of vcb MC events after all cuts (crooscheck to find biases)    
  	 sprintf(name, "%s%s",prefix,"vcballcutsdata1bintruth");	         
  	 signal = sighisto(((TH1D*)gDirectory->Get(name)),foomean,foosigma,fooalpha,foon,1,-1111111.,-1111111.,-1111111.,-1111111.,-1111111.);
  	 vcbmcselected2 = signal.X();
     
  	 //true number of vub MC events after all cuts (crooscheck to find biases)    
  	 sprintf(name, "%s%s",prefix,"vuballcutsdata1bintruth");	         
  	 signal = sighisto(((TH1D*)gDirectory->Get(name)),foomean,foosigma,fooalpha,foon,1,-1111111.,-1111111.,-1111111.,-1111111.,-1111111.);
  	 vubmcselected2 = signal.X();

  	 //true number of other MC events after all cuts (crooscheck to find biases)    
  	 sprintf(name, "%s%s",prefix,"othallcutsdata1bintruth");	         
  	 signal = sighisto(((TH1D*)gDirectory->Get(name)),foomean,foosigma,fooalpha,foon,1,-1111111.,-1111111.,-1111111.,-1111111.,-1111111.);
  	 othermcselected2 = signal.X();

 	 //true number of vub MC events after all cuts
 	 sprintf(name, "%s%s",prefix,"mxvub");	   
 	 vubmcaftertemp = ((TH1D*)gDirectory->Get(name))->Integral();


	 TH1D* htotal;
	 if(FITTOTSHAPE==2){
	   htotal = fitWithErrors(ich,ine,"mx");
	 } else {
	   htotal = fitWithoutErrors(ich,ine);
	 }
	 
	 // chisq calculation 
	 chisq = 0; 
	 double tempchisq = 0;
         sprintf(preich, "%s%d", "ch", ich);
         sprintf(preine, "%s%d", "ne", ine);
         sprintf(prefix, "%s%s",preich , preine);     

	 NDOF = 0;
	 
	 tempbinvub = 0;
	 for(int i=1;i<11;i++){
	   tempchisq = 0;
	   sprintf(name, "%s%s",prefix,"mxonebdata");	         
	   if(((TH1D*)gDirectory->Get(name))->GetBinContent(i)) {
	     tempbin = ((TH1D*)gDirectory->Get(name))->GetBinContent(i);
	     temperr = ((TH1D*)gDirectory->Get(name))->GetBinError(i);
	     sprintf(name, "%s%s",prefix,"mxonebvub");	         
	     tempbinvub = ((TH1D*)gDirectory->Get(name))->GetBinContent(i) * vubcomp;
	     tempbinerrvub = ((TH1D*)gDirectory->Get(name))->GetBinError(i) * vubcomp;
	     sprintf(name, "%s%s",prefix,"mxonebvcb");	         
	     tempbinvcb = ((TH1D*)gDirectory->Get(name))->GetBinContent(i) * vcbcomp;
	     tempbinerrvcb = ((TH1D*)gDirectory->Get(name))->GetBinError(i) * vcbcomp;
	     sprintf(name, "%s%s",prefix,"mxoneboth");	     
	     tempbinoth = ((TH1D*)gDirectory->Get(name))->GetBinContent(i) * othcomp;
	     tempbinerroth = ((TH1D*)gDirectory->Get(name))->GetBinError(i) * othcomp;
	     temperr = sqrt(temperr*temperr + tempbinerrvub*tempbinerrvub + tempbinerrvcb * tempbinerrvcb + tempbinerroth*tempbinerroth);
	     if(temperr<1) temperr = 1; 
	     sprintf(name, "%s%s",prefix,"mxonebdata");	
	     tempchisq = (tempbinvub + tempbinvcb + tempbinoth - ((TH1D*)gDirectory->Get(name))->GetBinContent(i))/temperr;     
	     NDOF++;
	   }   
	   chisq = chisq + tempchisq*tempchisq;
	 }
	 
	 chisq = chisq / (NDOF-3);
	 
	 cout << "Chi Square of the Fit :" << chisq << endl;
	 cout << "NDOF " << NDOF-3 << endl;

	 
	 // background subtraction (plot with 13 bins)
	 tempbin = 0;
	 temperr = 0;
	 tempbinvcb = 0;
	 tempbinoth = 0;
	 temperrvcb = 0;
	 temperroth = 0;
   
	 for(int i=1;i<14;i++){
	   sprintf(name, "%s%s",prefix,"mxdata");	         
	   tempbin = ((TH1D*)gDirectory->Get(name))->GetBinContent(i);
	   temperr = ((TH1D*)gDirectory->Get(name))->GetBinError(i);
	   sprintf(name, "%s%s",prefix,"mxvcb");	         
	   tempbinvcb = ((TH1D*)gDirectory->Get(name))->GetBinContent(i) * vcbcomp;
	   tempbinerrvcb = ((TH1D*)gDirectory->Get(name))->GetBinError(i) * vcbcomp;
	   temperrvcb = tempbinvcb * errvcbcomp / vcbcomp;
	   //temperrvcb = sqrt(temperrvcb * temperrvcb + tempbinerrvcb * tempbinerrvcb);     //THE MC STAT IS NOT INCLUDED 
	   sprintf(name, "%s%s",prefix,"mxoth");	     
	   tempbinoth = ((TH1D*)gDirectory->Get(name))->GetBinContent(i) * othcomp;
	   tempbinerroth = ((TH1D*)gDirectory->Get(name))->GetBinError(i) * othcomp;
	   temperroth = tempbinoth * errothcomp/othcomp;
	   //temperroth = sqrt(temperroth * temperroth + tempbinerroth * tempbinerroth);     //THE MC STAT IS NOT INCLUDED 
	   tempbin = tempbin - tempbinvcb - tempbinoth;
	   temperr = sqrt(temperr*temperr + tempbinerrvcb*tempbinerrvcb + tempbinerroth*tempbinerroth);
	   sprintf(name, "%s%s",prefix,"mxsubdata");
	   ((TH1D*)gDirectory->Get(name))->SetBinContent(i, tempbin);
	   ((TH1D*)gDirectory->Get(name))->SetBinError(i, temperr); 
	 }	
	 
	 int y;
	 sprintf(name, "%s%s",prefix,"mxscalevcb");
	 sprintf(name2, "%s%s",prefix,"mxvcb");   
	 for(y=1;y<15;y++){
	   ((TH1D*)gDirectory->Get(name))->SetBinContent(y,((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*vcbcomp);
	 }
   
	 sprintf(name, "%s%s",prefix,"mxscaleoth");
	 sprintf(name2, "%s%s",prefix,"mxoth");
	 for(y=1;y<15;y++){
	   ((TH1D*)gDirectory->Get(name))->SetBinContent(y,((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*othcomp);
	 }
	 
	 sprintf(name, "%s%s",prefix,"mxscalevub");
	 sprintf(name2, "%s%s",prefix,"mxvub");
	 for(y=1;y<15;y++){
	   ((TH1D*)gDirectory->Get(name))->SetBinContent(y,((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*vubcomp);
	 }
	 
	 sprintf(name, "%s%s",prefix,"mxallbkg");
	 sprintf(name2, "%s%s",prefix,"mxvcb");
	 for(y=1;y<11;y++){
	   ((TH1D*)gDirectory->Get(name))->SetBinContent(y,((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*vcbcomp);
	 }
	 sprintf(name2, "%s%s",prefix,"mxoth");
	 double it;
	 for(y=1;y<15;y++){
	   it = ((TH1D*)gDirectory->Get(name))->GetBinContent(y);
	   ((TH1D*)gDirectory->Get(name))->SetBinContent(y,it + ((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*othcomp);
	 }
	 sprintf(name, "%s%s",prefix,"mxallmc");
	 sprintf(name2, "%s%s",prefix,"mxvcb");
	 for(y=1;y<11;y++){
	   ((TH1D*)gDirectory->Get(name))->SetBinContent(y,((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*vcbcomp);
	 }
	 sprintf(name2, "%s%s",prefix,"mxoth");
	 for(y=1;y<15;y++){
	   it = ((TH1D*)gDirectory->Get(name))->GetBinContent(y);
	   ((TH1D*)gDirectory->Get(name))->SetBinContent(y,it + ((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*othcomp);
	 }
	 sprintf(name2, "%s%s",prefix,"mxvub");
	 for(y=1;y<15;y++){
	   it = ((TH1D*)gDirectory->Get(name))->GetBinContent(y);
	   ((TH1D*)gDirectory->Get(name))->SetBinContent(y,it + ((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*vubcomp);
	 }
	 
	 
	 // background subtraction (plot with 10 bins)
	 
	 for(int ik=1;ik<11;ik++){
	   sprintf(name, "%s%s",prefix,"mxonebdata");	         
	   tempbin = ((TH1D*)gDirectory->Get(name))->GetBinContent(ik);
	   temperr = ((TH1D*)gDirectory->Get(name))->GetBinError(ik);
	   sprintf(name, "%s%s",prefix,"mxonebvcb");	         
	   tempbinvcb = ((TH1D*)gDirectory->Get(name))->GetBinContent(ik) * vcbcomp;
	   tempbinerrvcb = ((TH1D*)gDirectory->Get(name))->GetBinError(ik) * vcbcomp;
	   temperrvcb = tempbinvcb * errvcbcomp / vcbcomp;
	   //temperrvcb = sqrt(temperrvcb * temperrvcb + tempbinerrvcb * tempbinerrvcb);     //THE MC STAT IS NOT INCLUDED 
	   sprintf(name, "%s%s",prefix,"mxoneboth");	         
	   tempbinoth = ((TH1D*)gDirectory->Get(name))->GetBinContent(ik) * othcomp;
	   tempbinerroth = ((TH1D*)gDirectory->Get(name))->GetBinError(ik) * othcomp;
	   temperroth = tempbinoth * errothcomp / othcomp;     
	   //temperroth = sqrt(temperroth * temperroth + tempbinerroth * tempbinerroth);     //THE MC STAT IS NOT INCLUDED 
	   tempbin = tempbin - tempbinvcb - tempbinoth;
	   temperr = sqrt(temperr*temperr + tempbinerrvcb*tempbinerrvcb + tempbinerroth*tempbinerroth);	
	   sprintf(name, "%s%s",prefix,"mxonebsubdata");	 
	   ((TH1D*)gDirectory->Get(name))->SetBinContent(ik, tempbin);
	   ((TH1D*)gDirectory->Get(name))->SetBinError(ik, temperr); 
	 }	
	 
	 sprintf(name, "%s%s",prefix,"mxonebvcb");	         
	 vcbaftercutsbin1 = ((TH1D*)gDirectory->Get(name))->GetBinContent(1) * vcbcomp;
	 errvcbaftercutsbin1 = ((TH1D*)gDirectory->Get(name))->GetBinError(1) * vcbcomp;
	 errfitvcbaftercutsbin1 = ((TH1D*)gDirectory->Get(name))->GetBinContent(1) * errvcbcomp;
	 sprintf(name, "%s%s",prefix,"mxoneboth");	         
	 otheraftercutsbin1 = ((TH1D*)gDirectory->Get(name))->GetBinContent(1) * othcomp;
	 errotheraftercutsbin1 = ((TH1D*)gDirectory->Get(name))->GetBinError(1) * othcomp;
	 errfitotheraftercutsbin1 = ((TH1D*)gDirectory->Get(name))->GetBinContent(1) * errothcomp;
	 
	 sprintf(name, "%s%s",prefix,"mxonebscalevcb");
	 sprintf(name2, "%s%s",prefix,"mxonebvcb");   
	 for(y=1;y<11;y++){
	   ((TH1D*)gDirectory->Get(name))->SetBinContent(y,((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*vubcomp);
	 }
	 
	 sprintf(name, "%s%s",prefix,"mxonebscaleoth");
	 sprintf(name2, "%s%s",prefix,"mxoneboth");
	 for(y=1;y<11;y++){
	   ((TH1D*)gDirectory->Get(name))->SetBinContent(y,((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*othcomp);
	 }
	 
	 sprintf(name, "%s%s",prefix,"mxonebscalevub");
	 sprintf(name2, "%s%s",prefix,"mxonebvub");
	 for(y=1;y<11;y++){
	   ((TH1D*)gDirectory->Get(name))->SetBinContent(y,((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*vcbcomp);
	 }
	 
	 sprintf(name, "%s%s",prefix,"mxoneballbkg");
	 sprintf(name2, "%s%s",prefix,"mxonebvcb");
	 for(y=1;y<11;y++){
	   ((TH1D*)gDirectory->Get(name))->SetBinContent(y,((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*vcbcomp);
	 }
	 sprintf(name2, "%s%s",prefix,"mxoneboth");
	 for(y=1;y<11;y++){
	   it = ((TH1D*)gDirectory->Get(name))->GetBinContent(y);
	   ((TH1D*)gDirectory->Get(name))->SetBinContent(y,it + ((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*othcomp);
	 }
	 
	 sprintf(name, "%s%s",prefix,"mxoneballmc");
	 sprintf(name2, "%s%s",prefix,"mxonebvcb");
	 for(y=1;y<11;y++){
	   ((TH1D*)gDirectory->Get(name))->SetBinContent(y,((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*vcbcomp);
	 }
	 sprintf(name2, "%s%s",prefix,"mxoneboth");
	 for(y=1;y<11;y++){
	   it = ((TH1D*)gDirectory->Get(name))->GetBinContent(y);
	   ((TH1D*)gDirectory->Get(name))->SetBinContent(y,it + ((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*othcomp);
	 }
	 sprintf(name2, "%s%s",prefix,"mxonebvub");
	 for(y=1;y<11;y++){
	   it = ((TH1D*)gDirectory->Get(name))->GetBinContent(y);
	   ((TH1D*)gDirectory->Get(name))->SetBinContent(y,it + ((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*vubcomp);
	 }
	 

	 
	 c1 = new TCanvas("c1"," ",200,10,1300,520); 
	 
	 sprintf(name, "%s%s%s%s",DIRNAME.Data(),PREFIXOUT.Data(),prefix,"fitresults_nocat.eps");
	 
	 openEpsFile(name);
	 
	 
	 c1->Clear();
	 
	 if(BLINDING){c1->Divide(1, 1);}else{c1->Divide(3, 1);}
	 
	 c1->cd(1);
	 
	 sprintf(name, "%s%s",prefix,"mxdata");
	 if(FITTOTSHAPE) sprintf(name, "%s%s",prefix,"mxonebdata");
	 ((TH1D*)gDirectory->Get(name))->SetMarkerSize(1.);	
	 ((TH1D*)gDirectory->Get(name))->SetMarkerColor(kBlack);	
	 ((TH1D*)gDirectory->Get(name))->SetStats(0);
	 double themax = ((TH1D*)gDirectory->Get(name))->GetMaximum() * 1.2; 
	 //	 if(BLINDING) sprintf(name,"%s%s",prefix,"mxblinddata" );
	 ((TH1D*)gDirectory->Get(name))->SetMaximum(themax); 
	 ((TH1D*)gDirectory->Get(name))->SetStats(0);
	 ((TH1D*)gDirectory->Get(name))->SetXTitle("Mx(GeV)");
	 ((TH1D*)gDirectory->Get(name))->Draw();
	 sprintf(name, "%s%s",prefix,"mxdata");
	 if(FITTOTSHAPE) sprintf(name, "%s%s",prefix,"mxonebdata");
	 ((TH1D*)gDirectory->Get(name))->SetMarkerStyle(20);
	 ((TH1D*)gDirectory->Get(name))->SetStats(0);
	 ((TH1D*)gDirectory->Get(name))->DrawCopy("same");
	 sprintf(name, "%s%s",prefix,"mxallmc");
	 if(FITTOTSHAPE) sprintf(name, "%s%s",prefix,"mxoneballmc");
	 //	 if(BLINDING)  ((TH1D*)gDirectory->Get(name))->SetAxisRange(1.7,5.);
	 //	 ((TH1D*)gDirectory->Get(name))->SetFillStyle(3004);
	 ((TH1D*)gDirectory->Get(name))->SetLineColor(38);
	 ((TH1D*)gDirectory->Get(name))->SetFillColor(38);
	 ((TH1D*)gDirectory->Get(name))->SetStats(0);
	 ((TH1D*)gDirectory->Get(name))->DrawCopy("same");  
	 ((TH1D*)gDirectory->Get(name))->SetAxisRange(0.,5.);
	 //	 if(!BLINDING) {
	   sprintf(name, "%s%s",prefix,"mxallbkg");
	   if(FITTOTSHAPE) sprintf(name, "%s%s",prefix,"mxoneballbkg");
	   //	   ((TH1D*)gDirectory->Get(name))->SetFillStyle(3005);
	   ((TH1D*)gDirectory->Get(name))->SetLineColor(kYellow);
	   ((TH1D*)gDirectory->Get(name))->SetFillColor(kYellow);
	   ((TH1D*)gDirectory->Get(name))->SetStats(0);
	   ((TH1D*)gDirectory->Get(name))->DrawCopy("same");
	   //	 }
	 sprintf(name, "%s%s",prefix,"mxscaleoth");
	 if(FITTOTSHAPE) sprintf(name, "%s%s",prefix,"mxonebscaleoth");
	 ((TH1D*)gDirectory->Get(name))->SetLineColor(13);
	 ((TH1D*)gDirectory->Get(name))->SetFillColor(13);
	 ((TH1D*)gDirectory->Get(name))->SetStats(0);
	 ((TH1D*)gDirectory->Get(name))->DrawCopy("same");
	 sprintf(name, "%s%s",prefix,"mxdata");
	 if(FITTOTSHAPE) sprintf(name, "%s%s",prefix,"mxonebdata");
	 ((TH1D*)gDirectory->Get(name))->DrawCopy("same");
	 TLegendEntry *leggemu; 
	 TLegend *legmu;
	 legmu = new TLegend(0.6,0.7,0.88,0.89);
	 legmu->SetFillStyle(0); legmu->SetBorderSize(0.); legmu->SetTextSize(0.05); 
	 legmu->SetFillColor(0); 
	 sprintf(name, "%s%s",prefix,"mxallmc");
	 leggemu = legmu->AddEntry(((TH1D*)gDirectory->Get(name)), "b->ulnu", "f"); 
	   sprintf(name, "%s%s",prefix,"mxallbkg");
	 leggemu = legmu->AddEntry(((TH1D*)gDirectory->Get(name)), "b->clnu", "f"); 
	 sprintf(name, "%s%s",prefix,"mxscaleoth");
	 leggemu = legmu->AddEntry(((TH1D*)gDirectory->Get(name)), "other", "f"); 
	 sprintf(name, "%s%s",prefix,"mxonebdata");
	 leggemu = legmu->AddEntry(((TH1D*)gDirectory->Get(name)), "data", "p"); 
	 legmu->Draw();
	 
//	 if(!BLINDING){
	   
	   c1->cd(2);
	   
	   sprintf(name, "%s%s",prefix,"mxdata");
	   themax = ((TH1D*)gDirectory->Get(name))->GetMaximum() * 1.2; 
	   //	   if(BLINDING) sprintf(name,"%s%s",prefix,"mxblinddata" );
	   ((TH1D*)gDirectory->Get(name))->SetMaximum(themax); 
	   ((TH1D*)gDirectory->Get(name))->SetMarkerStyle(20);
	   ((TH1D*)gDirectory->Get(name))->SetStats(0);
	   ((TH1D*)gDirectory->Get(name))->SetXTitle("Mx(GeV)");
	   ((TH1D*)gDirectory->Get(name))->Draw();
	   sprintf(name, "%s%s",prefix,"mxallmc");
	   //	   if(BLINDING)  ((TH1D*)gDirectory->Get(name))->SetAxisRange(1.7,5.);
	   //	   ((TH1D*)gDirectory->Get(name))->SetFillStyle(3004);
	   ((TH1D*)gDirectory->Get(name))->SetLineColor(38);
	   ((TH1D*)gDirectory->Get(name))->SetFillColor(38);
	   ((TH1D*)gDirectory->Get(name))->SetStats(0);
	   ((TH1D*)gDirectory->Get(name))->DrawCopy("same");  
	   ((TH1D*)gDirectory->Get(name))->SetAxisRange(0.,5.);
	   //	   if(!BLINDING) {
	     sprintf(name, "%s%s",prefix,"mxallbkg");
	     //	     ((TH1D*)gDirectory->Get(name))->SetFillStyle(3005);
	     ((TH1D*)gDirectory->Get(name))->SetLineColor(kYellow);
	     ((TH1D*)gDirectory->Get(name))->SetFillColor(kYellow);
	     ((TH1D*)gDirectory->Get(name))->SetStats(0);
	     ((TH1D*)gDirectory->Get(name))->DrawCopy("same");
	     //	   }
	   sprintf(name, "%s%s",prefix,"mxscaleoth");
	   ((TH1D*)gDirectory->Get(name))->SetLineColor(13);
	   ((TH1D*)gDirectory->Get(name))->SetFillColor(13);
	   ((TH1D*)gDirectory->Get(name))->SetStats(0);
	   ((TH1D*)gDirectory->Get(name))->DrawCopy("same");
	   sprintf(name, "%s%s",prefix,"mxdata");
	   ((TH1D*)gDirectory->Get(name))->DrawCopy("same");
	   legmu->Draw();
	   
	   c1->cd(3);
	   
	   sprintf(name, "%s%s",prefix,"mxsubdata");
	   max = 1.2*((TH1D*)gDirectory->Get(name))->GetMaximum();
	   min = 1.2*((TH1D*)gDirectory->Get(name))->GetMinimum();
	   ((TH1D*)gDirectory->Get(name))->SetMaximum(max);
	   ((TH1D*)gDirectory->Get(name))->SetMinimum(min);
	   ((TH1D*)gDirectory->Get(name))->SetMarkerStyle(20);
	   ((TH1D*)gDirectory->Get(name))->SetStats(0);
	   ((TH1D*)gDirectory->Get(name))->SetXTitle("Mx(GeV)");
	   sprintf(name, "%s%s",prefix,"mxscalevub");
	   //	   ((TH1D*)gDirectory->Get(name))->SetFillStyle(3004);
	   ((TH1D*)gDirectory->Get(name))->SetLineColor(38);
	   ((TH1D*)gDirectory->Get(name))->SetFillColor(38);
	   ((TH1D*)gDirectory->Get(name))->SetStats(0);
	   ((TH1D*)gDirectory->Get(name))->DrawCopy("same");  
	   sprintf(name, "%s%s",prefix,"mxsubdata");
	   ((TH1D*)gDirectory->Get(name))->DrawCopy("same");
	   TLegendEntry *legge2mul; 
	   TLegend *leg2mul;
	   leg2mul = new TLegend(0.5,0.7,0.88,0.89);
	   leg2mul->SetFillStyle(0); leg2mul->SetBorderSize(0.); leg2mul->SetTextSize(0.05); 
	   leg2mul->SetFillColor(0); 
	   sprintf(name, "%s%s",prefix,"mxscalevub");
	   legge2mul = leg2mul->AddEntry(((TH1D*)gDirectory->Get(name)), "scaled MC", "f"); 
	   sprintf(name, "%s%s",prefix,"mxsubdata");
	   legge2mul = leg2mul->AddEntry(((TH1D*)gDirectory->Get("mxonebdata")), "data subtr.", "p"); 
	   leg2mul->Draw();
	   TLine line(0.,0.,5.,0.);
	   line.SetLineColor(kRed); line.SetLineWidth(2); line.SetLineStyle(2); line.Draw();
	   ((TH1D*)gDirectory->Get(name))->Draw();
	   
	   //	 }
	 
	 closeEpsFile();
	 
	 delete c1;       
	 
	 // calculate ratio of BR
	 
 	 sprintf(name, "%s%s",prefix,"mxonebvub");
 	 vubmcaftercutsbin1 = ((TH1D*)gDirectory->Get(name))->GetBinContent(1);
 	 areavubmcaftercutsbin1 = ((TH1D*)gDirectory->Get(name))->Integral();
 	 errvubmcaftercutsbin1 = ((TH1D*)gDirectory->Get(name))->GetBinError(1);
	 
 	 sprintf(name, "%s%s",prefix,"mxonebsubdata");       
 	 vubaftercutsbin1 = ((TH1D*)gDirectory->Get(name))->GetBinContent(1);
 	 errvubaftercutsbin1 = ((TH1D*)gDirectory->Get(name))->GetBinError(1);
 	 if(FITTOTSHAPE){
 	   sprintf(name, "%s%s",prefix,"mxonebdata");
	   dataFirstBin= ((TH1D*)gDirectory->Get(name))->GetBinContent(1);
	   dataErrFirstBin= ((TH1D*)gDirectory->Get(name))->GetBinError(1);
 	   sprintf(name, "%s%s",prefix,"mxonebvcb");
	   areavcbaftercutsbin1 = ((TH1D*)gDirectory->Get(name))->Integral();
	   vcbFirstBin= ((TH1D*)gDirectory->Get(name))->GetBinContent(1);
	   vcbErrFirstBin= ((TH1D*)gDirectory->Get(name))->GetBinError(1);
 	   sprintf(name, "%s%s",prefix,"mxoneboth");
	   areaothaftercutsbin1 = ((TH1D*)gDirectory->Get(name))->Integral();
	   othFirstBin= ((TH1D*)gDirectory->Get(name))->GetBinContent(1);
	   othErrFirstBin= ((TH1D*)gDirectory->Get(name))->GetBinError(1);


 	   vubaftercutsbin1 = dataFirstBin-vcbFirstBin*vcbcomp-othFirstBin*othcomp;
	   MCerrvubaftercutsbin1 = sqrt(pow(vcbErrFirstBin*vcbcomp,2)+pow(othErrFirstBin*othcomp,2));
	   if(FITTOTSHAPE==2) {
	     errvubaftercutsbin1 = sqrt(pow(dataErrFirstBin,2)+pow(vcbFirstBin*errvcbcompNOMC,2)+pow(othFirstBin*errothcompNOMC,2));
	     MCerrvubaftercutsbin1 = sqrt(
					  pow(MCerrvubaftercutsbin1,2)
					  +pow(vcbFirstBin*errvcbcomp,2)+pow(othFirstBin*errothcomp,2)
					  -pow(vcbFirstBin*errvcbcompNOMC,2)-pow(othFirstBin*errothcompNOMC,2));
	   }else{
	     errvubaftercutsbin1 = sqrt(pow(dataErrFirstBin,2)+pow(vcbFirstBin*errvcbcomp,2)+pow(othFirstBin*errothcomp,2));
	   }
	 }
	 
	 if(BLINDING){
	   vubaftercutsbin1 = vubaftercutsbin1 * blindfactor;
	   errvubaftercutsbin1 = errvubaftercutsbin1 * blindfactor;
	 }
	 
 	 S = vubaftercutsbin1;
 	 errS = errvubaftercutsbin1;
 	 Bvcb = vcbaftercutsbin1;
 	 temperrvcb = sqrt(errvcbaftercutsbin1 * errvcbaftercutsbin1 + errfitvcbaftercutsbin1 * errfitvcbaftercutsbin1);
 	 errBvcb = temperrvcb;
 	 Bother = otheraftercutsbin1;
 	 temperroth = sqrt(errotheraftercutsbin1 * errotheraftercutsbin1 + errfitotheraftercutsbin1 * errfitotheraftercutsbin1);
 	 errBother = temperroth;
 	 B = Bvcb + Bother;
 	 errB = sqrt(errBvcb * errBvcb + errBother * errBother);
 
 	 double epsucat = vubmc2aftercutsbin1/vubmc2cut;
 	 double errepsucat = sqrt(epsucat*(1 - epsucat)/vubmc2cut) ;
 	 double epsmxcat = vubmc2cut/vubmc2;
 	 double errepsmxcat = sqrt(epsmxcat*(1 - epsmxcat)/vubmc2) ;   
 	 double epstotcat = epsucat * epsmxcat;
 	 double errepstotcat = sqrt(epsucat*epsucat*errepsmxcat*errepsmxcat + epsmxcat*epsmxcat*errepsucat*errepsucat);
	 
 	 fact = vcbmc * (1 + BRRATIOVALUETAIL/0.104) / (totmc + vcbmc * BRRATIOVALUETAIL/0.104); //factor to take into account the vub contribute into all events after lept cut
 	 //BRRATIOVALUETAIL for BR(b->ulnu) is assumed
	 
 	 //now pstarfact calculated in the fit	 
 	 BRBR = S * (1/(tot*fact*calcpstarfact)) / (epsucat * epsmxcat);
 	 errBRBR = errS * (1/(tot*fact*calcpstarfact)) / (epsucat * epsmxcat);	
	 errBRBRMCstat=BRBR*
	   sqrt(pow(MCerrvubaftercutsbin1/S,2)+pow(errepstotcat/epstotcat,2));
 	 errBRBRtheo = (2.1002 -1.0995 * mxcut +  0.29050 * mxcut*mxcut);
 	 errBRBRtheo = (errBRBRtheo-1.) * 2. * BRBR ;
 	 errtotalBRBR = sqrt(errBRBR * errBRBR + errBRBRtheo * errBRBRtheo);
	 
	 //number to be combined
	 vubafter[ich-1][ine-1] = S / epsucat;
	 errvubafter[ich-1][ine-1] = errS / epsucat;
	 errvubafterMCstat[ich-1][ine-1] = sqrt(errvcbaftercutsbin1*errvcbaftercutsbin1 + errotheraftercutsbin1+errotheraftercutsbin1 + (errepsucat/epsucat) * (errepsucat/epsucat) * S * S);
	 
 	 // print results on the screen
 	 cout << endl;
 	 cout << endl;
 	 cout << endl;
 	 cout << "BRBR = " << BRBR << " +- " << errBRBR << "(stat) +- " << errBRBRMCstat << "(MC stat)" << endl;
 	 cout << endl;
 	 cout << endl;
 	 cout << endl;
       
 	 // print results on file
 	 sprintf(name, "%s%s%s%s",DIRNAME.Data(),PREFIXOUT.Data(),prefix,"result_nocat.dat");        
 	 ofstream outfile(name);

 	 outfile << "ALL NUMBERS" << endl;
 	 outfile << endl;
 	 outfile << "MX FIT" << endl;
 	 outfile << endl;
	 outfile << "Vub comp = " << vubcomp << " +- " << errvubcomp << endl;
 	 outfile << "Vcb comp = " << vcbcomp << " +- " << errvcbcomp << endl;
 	 outfile << "oth comp = " << othcomp << " +- " << errothcomp << endl;
 	 outfile << endl;
 	 outfile << "Vub tot area = " << ((TH1D*)gDirectory->Get("mxonebvub"))->Integral() << endl;
 	 outfile << "Vcb tot area = " << ((TH1D*)gDirectory->Get("mxonebvcb"))->Integral() << endl;
 	 outfile << "oth tot area = " << ((TH1D*)gDirectory->Get("mxoneboth"))->Integral() << endl;
 	 outfile << endl;
 	 outfile << "Vub 1' bin = " << ((TH1D*)gDirectory->Get("mxonebvub"))->GetBinContent(1) << " +- " << ((TH1D*)gDirectory->Get("mxonebvub"))->GetBinError(1) << endl;
 	 outfile << "Vcb 1' bin = " << ((TH1D*)gDirectory->Get("mxonebvcb"))->GetBinContent(1) << " +- " << ((TH1D*)gDirectory->Get("mxonebvcb"))->GetBinError(1) << endl;
 	 outfile << "oth 1' bin = " << ((TH1D*)gDirectory->Get("mxoneboth"))->GetBinContent(1) << " +- " << ((TH1D*)gDirectory->Get("mxoneboth"))->GetBinError(1) << endl;
 	 outfile << endl;
 	 outfile << "WEIGHTED NUMBERS" << endl;
 	 outfile << endl;
 	 if(!BLINDING)      outfile << "Vub fitted 1' bin = " << S << " +- " <<  ((TH1D*)gDirectory->Get("mxonebvub"))->GetBinError(1)*vubcomp << "(stat MC) +- " << errS << "(err fit)" << endl;
 	 outfile << "Vcb fitted 1' bin = " << vcbaftercutsbin1 << " +- " << errvcbaftercutsbin1 << "(stat MC) +- " << errfitvcbaftercutsbin1 << "(err fit)" << endl;
 	 outfile << "oth fitted 1' bin = " << otheraftercutsbin1 << " +- " << errotheraftercutsbin1 << "(stat MC) +- " << errfitotheraftercutsbin1 << "(err fit)"  << endl;
 	 outfile << endl;
 	 outfile << "Vub 1' true = " << vubmcselected2 << endl;
 	 outfile << "Vcb 1' true = " << vcbmcselected2 << endl;
 	 outfile << "oth 1' true = " << othermcselected2 << endl;
 	 outfile << endl;
 	 outfile << "EFFICIENCY Vub" << endl;
 	 outfile << endl;
 	 outfile << "Vub total MC (lepton cut) = " << vubmc2 << endl;
 	 outfile << "Vub total MC (lepton cut + Mx cut) = " << vubmc2cut << endl;
 	 outfile << "Vub MC (all cuts) = " << vubmcaftertemp << endl;
 	 outfile << "Vub MC (all cuts + Mx cut) = " << vubmcaftercutsbin1 << endl;
 	 outfile << "Eps_u =  " << epsucat << " +- " << errepsucat << endl;
 	 outfile << "Eps_Mx = " << epsmxcat << " +- " << errepsmxcat << endl;
 	 outfile << "Eps_tot = " << epstotcat << " +- " << errepstotcat << endl;
 	 outfile << endl; 
 	 outfile << endl;
 	 outfile << "NSL NUMBERS" << endl;
 	 outfile << endl;
 	 outfile << "Nsl = " << tot << endl;
 	 outfile << endl;
 	 outfile << "Nsl MC = " << totmc << endl;
 	 outfile << "Nsl - BG = " << vcbmc << endl;
 	 outfile << "(Nsl - BG)/Nsl (vub corrected) = " << fact << endl;
 	 outfile << endl;
 	 outfile << "Pstar fact = " << PSTARFACT << endl;
 	 outfile << endl;
 	 outfile << endl; 
 	 outfile << "Eff on Vcb (Nvcb/(Nsl-BGsl)) = " << ((TH1D*)gDirectory->Get("mxonebvcb"))->Integral() * vcbcomp / (tot * fact) << " +- " <<  ((TH1D*)gDirectory->Get("mxonebvcb"))->Integral() * errvcbcomp / (tot * fact) << endl;  
 	 outfile << "Eff on other (other/(Nsl-BGsl)) = " << ((TH1D*)gDirectory->Get("mxoneboth"))->Integral() * othcomp / (tot * fact) << " +- " << ((TH1D*)gDirectory->Get("mxoneboth"))->Integral() * errothcomp / (tot * fact)  << endl;  
 	 outfile << endl; 
 	 outfile << endl;
 	 outfile << "Vub gene MC (no cuts) = " << vubmcnocut << endl;
 	 outfile << "Vcb gene MC (no cuts) = " << vcbmcnocut << endl;
 	 outfile << "Oth gene MC (no cuts) = " << othmcnocut << endl;
 	 outfile << "Recalculated Pstar fact = " << calcpstarfacttemp << endl;
 	 outfile << "Used Pstar fact = " << calcpstarfact << endl;
 	 outfile << endl;
 	 outfile << endl; 
 	 outfile << "Nsig = " << S/epstotcat  << " +- " <<  errS/epstotcat   << endl;
       
 	 outfile << endl;
 	 outfile << endl; 
 	 outfile << "BRBR = " << BRBR << " +- " << errBRBR << "(stat) +- " << errBRBRMCstat << "(MC stat)" << endl;
 	 outfile << endl;
 	 outfile << "Chi Square of the Fit = " << chisq << endl;
 	 outfile << "NDOF = " << NDOF-3 << endl;
 	 outfile << endl;
 	 outfile << endl;
 	 double Vub = 0.00445 * sqrt((BRBR * .104 * 1.55) / (0.002 * 1.622));
 	 outfile << "Vub(*10-3) = " << Vub*1000 << " +- " << Vub*1000*(errBRBR/(2*BRBR)) << "(stat) +- " << Vub*1000*errBRBRMCstat/(2*BRBR) << "(MC stat)" << endl;
 	 outfile << endl;
 	 outfile << endl;
 	 outfile << endl;
 	 outfile << endl;         
 	 outfile.close();
	 delete htotal;

 	 sprintf(name, "%s%s%s%s.tex",DIRNAME.Data(),PREFIXOUT.Data(),prefix,texPrefix.Data());        
	 cout << "**********************************************************************" << endl;
	 cout << "DUMPING TEXFILE 1 " << name << endl;
	 cout << "**********************************************************************" << endl;
 	 ofstream texfile(name);
	 char texline[200];
	 sprintf(texline, "\\def\\%ssemil{%7.0f \\pm %6.0f} ", texPrefix.Data(), tot, errtot); 
	 texfile << texline << endl;
	 sprintf(texline, "\\def\\%sslsub{%5.0f \\pm %4.0f} ", texPrefix.Data(), tot - (1.-fact)*tot, errtot); 
	 texfile << texline << endl;
	 sprintf(texline, "\\def\\%sbgsl{%5.0f \\pm %4.0f} ", texPrefix.Data(), (1.-fact)*tot, (1.-fact)*errtot); 
	 texfile << texline << endl;
	 sprintf(texline, "\\def\\%sepssl{%5.2f \\pm %4.2f} ", texPrefix.Data(), calcpstarfact, errcalcpstarfact); 
	 texfile << texline << endl;
	 sprintf(texline, "\\def\\%snu{%5.0f \\pm %4.0f} ", texPrefix.Data(), S, errS); 
	 texfile << texline << endl;
	 sprintf(texline, "\\def\\%sbgc{%5.0f \\pm %4.0f} ", texPrefix.Data(), vcbaftercutsbin1, errfitvcbaftercutsbin1); 
	 texfile << texline << endl;
	 sprintf(texline, "\\def\\%sbgo{%5.0f \\pm %4.0f} ", texPrefix.Data(), otheraftercutsbin1, errfitotheraftercutsbin1); 
	 texfile << texline << endl;
	 sprintf(texline, "\\def\\%sepsu{%4.3f} ", texPrefix.Data(), epsucat); 
	 texfile << texline << endl;
	 sprintf(texline, "\\def\\%sepsmx{%4.3f} ", texPrefix.Data(), epsmxcat); 
	 texfile << texline << endl;
	 sprintf(texline, "\\def\\%sbrbr{%6.4f \\pm %5.4f} ", texPrefix.Data(), BRBR, errBRBR); 
	 texfile << texline << endl;
	 sprintf(texline, "\\def\\%sbrbrerrmc{%5.4f} ", texPrefix.Data(), errBRBRMCstat); 
	 texfile << texline << endl;
	 texfile.close(); 
	 
       }
     }
     
  
     // put together all the multiplicity results
   
     double vubaftertot = 0;
     double errvubaftertot = 0;
     double errvubafterMCstattot = 0;
     
     for (ich=1;ich<4;ich++){
       for (ine=1;ine<3;ine++){    
	 if(((vubafter[ich-1][ine-1]>0) || (vubafter[ich-1][ine-1]<0) || ( vubafter[ich-1][ine-1]== 0))&&TMath::Abs(vubafter[ich-1][ine-1])<10000000){
	   vubaftertot =  vubafter[ich-1][ine-1]+ vubaftertot;
	   errvubaftertot = sqrt(errvubafter[ich-1][ine-1]*errvubafter[ich-1][ine-1] + errvubaftertot*errvubaftertot);
	   errvubafterMCstattot = sqrt(errvubafterMCstat[ich-1][ine-1]*errvubafterMCstat[ich-1][ine-1] + errvubaftertot*errvubaftertot);
	 }

       }
     }

     S = vubaftertot;
     errS = errvubaftertot;

     double epsmxtot = vubmccut/vubmcleptforeff;
     double errepsmxtot = sqrt(epsmxtot*(1 - epsmxtot)/vubmcleptforeff) ;   

     fact = vcbmc * (1 + BRRATIOVALUETAIL/0.104) / (totmc + vcbmc * BRRATIOVALUETAIL/0.104); //factor to take into account the vub contribute into all events after lept cut
     //BRRATIOVALUETAIL for BR(b->ulnu) is assumed
 
     //now pstarfact calculated in the fit	 
     BRBR = S * (1/(tot*fact*calcpstarfact)) / (epsmxtot);
     errBRBR = errS * (1/(tot*fact*calcpstarfact)) / (epsmxtot);
     errBRBRMCstat = BRBR*errvubafterMCstattot/S;
     errBRBRtheo = (2.1002 -1.0995 * mxcut +  0.29050 * mxcut*mxcut);
     errBRBRtheo = (errBRBRtheo-1.) * 2. * BRBR ;
     errtotalBRBR = sqrt(errBRBR * errBRBR + errBRBRtheo * errBRBRtheo);
    
     // print results on the screen
     cout << endl;
     cout << endl;
     cout << "COMBINED RESULT:" << endl;
     cout << endl;
     cout << "BRBR = " << BRBR << " +- " << errBRBR << "(stat) +- " << errBRBRMCstat << "(MC stat)" << endl;
     cout << endl;
     cout << endl;
     cout << endl;
     
     sprintf(name, "%s%s%s",DIRNAME.Data(),PREFIXOUT.Data(),"resultotal.dat");        
     ofstream outfilefin(name);
     
     // print results on file
  
     outfilefin << "ALL NUMBERS" << endl;
     outfilefin << endl;
     outfilefin << "MX FIT" << endl;
     outfilefin << endl;
     outfilefin << "Vub comp tot= " << S << " +- " << errS <<  " +- " << errvubafterMCstattot <<endl;
     outfilefin << endl;
     outfilefin << "EFFICIENCY Vub" << endl;
     outfilefin << endl;
     outfilefin << "Vub total MC (lepton cut) = " << vubmcleptforeff << endl;
     outfilefin << "Vub total MC (lepton cut + Mx cut) = " << vubmccut << endl;
     outfilefin << "Eps_Mx = " << epsmxtot << " +- " << errepsmxtot << endl;
     outfilefin << endl; 
     outfilefin << endl;
     outfilefin << "NSL NUMBERS" << endl;
     outfilefin << endl;
     outfilefin << "Nsl = " << tot << endl;
     outfilefin << endl;
     outfilefin << "Nsl MC = " << totmc << endl;
     outfilefin << "Nsl - BG = " << vcbmc << endl;
     outfilefin << "(Nsl - BG)/Nsl (vub corrected) = " << fact << endl;
     outfilefin << endl;
     outfilefin << "Pstar fact = " << PSTARFACT << endl;
     outfilefin << endl;
     outfilefin << endl; 
     outfilefin << "Vub gene MC (no cuts) = " << vubmcnocut << endl;
     outfilefin << "Vcb gene MC (no cuts) = " << vcbmcnocut << endl;
     outfilefin << "Oth gene MC (no cuts) = " << othmcnocut << endl;
     outfilefin << "Recalculated Pstar fact = " << calcpstarfacttemp << endl;
     outfilefin << endl;
     outfilefin << endl; 
     outfilefin << "Nsig = " << S/epsmxtot  << " +- " <<  errS/epsmxtot   << endl;
     
     outfilefin << endl;
     outfilefin << endl; 
     outfilefin << "BRBR = " << BRBR << " +- " << errBRBR << "(stat) +- " << errBRBRMCstat << "(MC stat)" << endl;
     outfilefin << endl;
     outfilefin << endl;
     double Vub = 0.00445 * sqrt((BRBR * .104 * 1.55) / (0.002 * 1.622));
     outfilefin << "Vub(*10-3) = " << Vub*1000 << " +- " << Vub*1000*(errBRBR/(2*BRBR)) << "(stat) +- " << Vub*1000*errBRBRMCstat/(2*BRBR) << "(MC stat)" << endl;
     outfilefin << endl;
     outfilefin << endl;
     outfilefin << endl;
     outfilefin << endl;         
     outfilefin.close();

     sprintf(name, "%s%s%s.tex",DIRNAME.Data(),PREFIXOUT.Data(),texPrefix.Data());        
     cout << "**********************************************************************" << endl;
     cout << "DUMPING TEXFILE 2 " << name << endl;
     cout << "**********************************************************************" << endl;
     ofstream texfile(name);
     char texline[200];
     sprintf(texline, "\\def\\%ssemil{%7.0f \\pm %6.0f} ", texPrefix.Data(), tot, errtot); 
     texfile << texline << endl;
     sprintf(texline, "\\def\\%sslsub{%5.0f \\pm %4.0f} ", texPrefix.Data(), tot - (1.-fact)*tot, errtot); 
     texfile << texline << endl;
     sprintf(texline, "\\def\\%sbgsl{%5.0f \\pm %4.0f} ", texPrefix.Data(), (1.-fact)*tot, (1.-fact)*errtot); 
     texfile << texline << endl;
     sprintf(texline, "\\def\\%sepssl{%5.2f \\pm %4.2f} ", texPrefix.Data(), calcpstarfact, errcalcpstarfact); 
     texfile << texline << endl;
     sprintf(texline, "\\def\\%snu{%5.0f \\pm %4.0f} ", texPrefix.Data(), S, errS); 
     texfile << texline << endl;
     sprintf(texline, "\\def\\%sbgc{%5.0f \\pm %4.0f} ", texPrefix.Data(), vcbaftercutsbin1, errfitvcbaftercutsbin1); 
     texfile << texline << endl;
     sprintf(texline, "\\def\\%sbgo{%5.0f \\pm %4.0f} ", texPrefix.Data(), otheraftercutsbin1, errfitotheraftercutsbin1); 
     texfile << texline << endl;
     sprintf(texline, "\\def\\%sepsu{%4.3f} ", texPrefix.Data(), 0.); 
     texfile << texline << endl;
     sprintf(texline, "\\def\\%sepsmx{%4.3f} ", texPrefix.Data(), epsmxtot);
     texfile << texline << endl;
     sprintf(texline, "\\def\\%sbrbr{%6.4f \\pm %5.4f} ", texPrefix.Data(), BRBR, errBRBR); 
     texfile << texline << endl;
     sprintf(texline, "\\def\\%sbrbrerrmc{%5.4f} ", texPrefix.Data(), errBRBRMCstat); 
     texfile << texline << endl;
     texfile.close(); 

   }
     

   // FIT WITH NO CATEGORIES to Q2 distribution
   if(FITQ2) {     theq2fit();   }

   // FIT WITH NO CATEGORIES to mx distribution
   
   //true number of vub MC events after all cuts
   sprintf(name, "%s","mxvub");	   
   vubmcaftercuts = ((TH1D*)gDirectory->Get(name))->Integral();
   
   TH1D* htotal;
   if(FITTOTSHAPE==2){
     htotal = fitWithErrors(-1,-1,"mx");
   } else {
     htotal = fitWithoutErrors(-1,-1);
   }

   cout<<"************ :: "<<vubcomp<<" "<<vcbcomp<<" "<<othcomp<<endl;

   // chisq calculation 
   chisq = 0; 
   double tempchisq = 0;
   
   cout << "CHISQ CALCULATION" << endl;

   NDOF = 0;

   for(int i=1;i<11;i++){
     tempchisq = 0;
     sprintf(name, "%s","mxonebdata");	         
     if(((TH1D*)gDirectory->Get(name))->GetBinContent(i)) {
       tempbin = ((TH1D*)gDirectory->Get(name))->GetBinContent(i);
       temperr = ((TH1D*)gDirectory->Get(name))->GetBinError(i);
       sprintf(name, "%s","mxonebvub");	         
       tempbinvub = ((TH1D*)gDirectory->Get(name))->GetBinContent(i) * vubcomp;
       tempbinerrvub = ((TH1D*)gDirectory->Get(name))->GetBinError(i) * vubcomp;
       sprintf(name, "%s","mxonebvcb");	         
       tempbinvcb = ((TH1D*)gDirectory->Get(name))->GetBinContent(i) * vcbcomp;
       tempbinerrvcb = ((TH1D*)gDirectory->Get(name))->GetBinError(i) * vcbcomp;
       sprintf(name, "%s","mxoneboth");	     
       tempbinoth = ((TH1D*)gDirectory->Get(name))->GetBinContent(i) * othcomp;
       tempbinerroth = ((TH1D*)gDirectory->Get(name))->GetBinError(i) * othcomp;
       temperr = sqrt(temperr*temperr + tempbinerrvub*tempbinerrvub + tempbinerrvcb * tempbinerrvcb + tempbinerroth*tempbinerroth);
       if(temperr<1) temperr = 1; 
       sprintf(name, "%s","mxonebdata");	
       tempchisq = (tempbinvub + tempbinvcb + tempbinoth - ((TH1D*)gDirectory->Get(name))->GetBinContent(i))/temperr;     
       NDOF++;
     }   
     chisq = chisq + tempchisq*tempchisq;
   }
   
   chisq = chisq / (NDOF-3);
   
   cout << " Chi Square of the Fit :" << chisq << endl;
   cout << "NDOF " << NDOF-3 << endl;

   // background subtraction (plot with 13 bins)
   tempbin = 0;
   temperr = 0;
   tempbinvcb = 0;
   tempbinoth = 0;
   temperrvcb = 0;
   temperroth = 0;
   
   for(int i=1;i<14;i++){
     sprintf(name, "%s","mxdata");	         
     tempbin = ((TH1D*)gDirectory->Get(name))->GetBinContent(i);
     temperr = ((TH1D*)gDirectory->Get(name))->GetBinError(i);
     sprintf(name, "%s","mxvcb");	         
     tempbinvcb = ((TH1D*)gDirectory->Get(name))->GetBinContent(i) * vcbcomp;
     tempbinerrvcb = ((TH1D*)gDirectory->Get(name))->GetBinError(i) * vcbcomp;
     temperrvcb = tempbinvcb * errvcbcomp / vcbcomp;
     //temperrvcb = sqrt(temperrvcb * temperrvcb + tempbinerrvcb * tempbinerrvcb);     //THE MC STAT IS NOT INCLUDED 
     sprintf(name, "%s","mxoth");	     
     tempbinoth = ((TH1D*)gDirectory->Get(name))->GetBinContent(i) * othcomp;
     tempbinerroth = ((TH1D*)gDirectory->Get(name))->GetBinError(i) * othcomp;
     temperroth = tempbinoth * errothcomp / othcomp;
     //temperroth = sqrt(temperroth * temperroth + tempbinerroth * tempbinerroth);     //THE MC STAT IS NOT INCLUDED 
     tempbin = tempbin - tempbinvcb - tempbinoth;
     temperr = sqrt(temperr*temperr + tempbinerrvcb*tempbinerrvcb + tempbinerroth*temperroth); 
     sprintf(name, "%s","mxsubdata");
     ((TH1D*)gDirectory->Get(name))->SetBinContent(i, tempbin);
     ((TH1D*)gDirectory->Get(name))->SetBinError(i, temperr); 
   }	
      
   int y;
   sprintf(name, "%s","mxscalevcb");
   sprintf(name2, "%s","mxvcb");   
   for(y=1;y<15;y++){
     ((TH1D*)gDirectory->Get(name))->SetBinContent(y,((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*vcbcomp);
   }
   
   sprintf(name, "%s","mxscaleoth");
   sprintf(name2, "%s","mxoth");
   for(y=1;y<15;y++){
     ((TH1D*)gDirectory->Get(name))->SetBinContent(y,((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*othcomp);
   }
   
   sprintf(name, "%s","mxscalevub");
   sprintf(name2, "%s","mxvub");
   for(y=1;y<15;y++){
     ((TH1D*)gDirectory->Get(name))->SetBinContent(y,((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*vubcomp);
   }
   
   sprintf(name, "%s","mxallbkg");
   sprintf(name2, "%s","mxvcb");
   for(y=1;y<15;y++){
     ((TH1D*)gDirectory->Get(name))->SetBinContent(y,((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*vcbcomp);
   }
   sprintf(name2, "%s","mxoth");
   double it;
   for(y=1;y<15;y++){
     it = ((TH1D*)gDirectory->Get(name))->GetBinContent(y);
     ((TH1D*)gDirectory->Get(name))->SetBinContent(y,it + ((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*othcomp);
   }
   sprintf(name, "%s","mxallmc");
   sprintf(name2, "%s","mxvcb");
   for(y=1;y<15;y++){
     ((TH1D*)gDirectory->Get(name))->SetBinContent(y,((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*vcbcomp);
   }
   sprintf(name2, "%s","mxoth");
   for(y=1;y<15;y++){
     it = ((TH1D*)gDirectory->Get(name))->GetBinContent(y);
     ((TH1D*)gDirectory->Get(name))->SetBinContent(y,it + ((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*othcomp);
   }
   sprintf(name2, "%s","mxvub");
   for(y=1;y<15;y++){
     it = ((TH1D*)gDirectory->Get(name))->GetBinContent(y);
     ((TH1D*)gDirectory->Get(name))->SetBinContent(y,it + ((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*vubcomp);
   }
   
   
   // background subtraction (plot with 10 bins)

   for(int ik=1;ik<11;ik++){
     sprintf(name, "%s","mxonebdata");	         
     tempbin = ((TH1D*)gDirectory->Get(name))->GetBinContent(ik);
     temperr = ((TH1D*)gDirectory->Get(name))->GetBinError(ik);
     sprintf(name, "%s","mxonebvcb");	         
     tempbinvcb = ((TH1D*)gDirectory->Get(name))->GetBinContent(ik) * vcbcomp;
     tempbinerrvcb = ((TH1D*)gDirectory->Get(name))->GetBinError(ik) * vcbcomp;
     temperrvcb = tempbinvcb * errvcbcomp / vcbcomp;
     //temperrvcb = sqrt(temperrvcb * temperrvcb + tempbinerrvcb * tempbinerrvcb);     //THE MC STAT IS NOT INCLUDED 
     sprintf(name, "%s","mxoneboth");	         
     tempbinoth = ((TH1D*)gDirectory->Get(name))->GetBinContent(ik) * othcomp;
     tempbinerroth = ((TH1D*)gDirectory->Get(name))->GetBinError(ik) * othcomp;
     temperroth = tempbinoth * errothcomp / othcomp;     
     //temperroth = sqrt(temperroth * temperroth + tempbinerroth * tempbinerroth);     //THE MC STAT IS NOT INCLUDED 
     tempbin = tempbin - tempbinvcb - tempbinoth;
     temperr = sqrt(temperr*temperr + tempbinerrvcb*tempbinerrvcb + tempbinerroth*tempbinerroth);
     sprintf(name, "%s","mxonebsubdata");	 
     ((TH1D*)gDirectory->Get(name))->SetBinContent(ik, tempbin);
     ((TH1D*)gDirectory->Get(name))->SetBinError(ik, temperr); 
   }	
   
   sprintf(name, "%s","mxonebvcb");	         
   vcbaftercutsbin1 = ((TH1D*)gDirectory->Get(name))->GetBinContent(1) * vcbcomp;
   errvcbaftercutsbin1 = ((TH1D*)gDirectory->Get(name))->GetBinError(1) * vcbcomp;
   errfitvcbaftercutsbin1 = ((TH1D*)gDirectory->Get(name))->GetBinContent(1) * errvcbcomp;
   sprintf(name, "%s","mxoneboth");	         
   otheraftercutsbin1 = ((TH1D*)gDirectory->Get(name))->GetBinContent(1) * othcomp;
   errotheraftercutsbin1 = ((TH1D*)gDirectory->Get(name))->GetBinError(1) * othcomp;
   errfitotheraftercutsbin1 = ((TH1D*)gDirectory->Get(name))->GetBinContent(1) * errothcomp;

   sprintf(name, "%s","mxonebscalevcb");
   sprintf(name2, "%s","mxonebvcb");   
   for(y=1;y<11;y++){
     ((TH1D*)gDirectory->Get(name))->SetBinContent(y,((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*vcbcomp);
   }
   
   sprintf(name, "%s","mxonebscaleoth");
   sprintf(name2, "%s","mxoneboth");
   for(y=1;y<11;y++){
     ((TH1D*)gDirectory->Get(name))->SetBinContent(y,((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*othcomp);
   }
   
   sprintf(name, "%s","mxonebscalevub");
   sprintf(name2, "%s","mxonebvub");
   for(y=1;y<11;y++){
     ((TH1D*)gDirectory->Get(name))->SetBinContent(y,((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*vubcomp);
   }
   
   sprintf(name, "%s","mxoneballbkg");
   sprintf(name2, "%s","mxonebvcb");
   for(y=1;y<11;y++){
     ((TH1D*)gDirectory->Get(name))->SetBinContent(y,((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*vcbcomp);
   }
   sprintf(name2, "%s","mxoneboth");
   for(y=1;y<11;y++){
     it = ((TH1D*)gDirectory->Get(name))->GetBinContent(y);
     ((TH1D*)gDirectory->Get(name))->SetBinContent(y,it + ((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*othcomp);
   }

   sprintf(name, "%s","mxoneballmc");
   sprintf(name2, "%s","mxonebvcb");
   for(y=1;y<11;y++){
     ((TH1D*)gDirectory->Get(name))->SetBinContent(y,((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*vcbcomp);
   }
   sprintf(name2, "%s","mxoneboth");
   for(y=1;y<11;y++){
     it = ((TH1D*)gDirectory->Get(name))->GetBinContent(y);
     ((TH1D*)gDirectory->Get(name))->SetBinContent(y,it + ((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*othcomp);
   }
   sprintf(name2, "%s","mxonebvub");
   for(y=1;y<11;y++){
     it = ((TH1D*)gDirectory->Get(name))->GetBinContent(y);
     ((TH1D*)gDirectory->Get(name))->SetBinContent(y,it + ((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*vubcomp);
   }

   if(DOVARSTU) {
     char varNm[200];
     for (int iVarm = 0;  iVarm < 6; iVarm++) {       
       if(iVarm == 0) {
	 sprintf(varNm,"mqb");
       } else if(iVarm == 1) {
	 sprintf(varNm,"q2");
       } else if(iVarm == 2) {
	 sprintf(varNm,"csi");
       } else if(iVarm == 3) {
	 sprintf(varNm,"w");
       } else if(iVarm == 4) {
	 sprintf(varNm,"x");
       } else if(iVarm == 5) {
	 sprintf(varNm,"lep");
       } 
       cout<< varNm << " "<< iVarm<<endl;
       varStudy(varNm);
     }
   }
   c1 = new TCanvas("c1"," ",200,10,1300,520); 
   
   sprintf(name, "%s%s%s",DIRNAME.Data(),PREFIXOUT.Data(),"fitresults_nocat.eps");
   
   openEpsFile(name);
   
   
   c1->Clear();

   c1->Divide(3, 1);
   c1->cd(1);
   
   sprintf(name, "%s","mxdata");
   if(FITTOTSHAPE) sprintf(name, "%s","mxonebdata");
   ((TH1D*)gDirectory->Get(name))->SetMarkerSize(1.);	
   ((TH1D*)gDirectory->Get(name))->SetMarkerColor(kBlack);	
   ((TH1D*)gDirectory->Get(name))->SetStats(0);
   ((TH1D*)gDirectory->Get(name))->SetXTitle("Mx(GeV)");
   double themax = ((TH1D*)gDirectory->Get(name))->GetMaximum() * 1.2; 
   ((TH1D*)gDirectory->Get(name))->SetMaximum(themax); 
   ((TH1D*)gDirectory->Get(name))->SetStats(0);
   ((TH1D*)gDirectory->Get(name))->Draw();
   sprintf(name, "%s","mxdata");
   if(FITTOTSHAPE) sprintf(name, "%s","mxonebdata");
   ((TH1D*)gDirectory->Get(name))->SetMarkerStyle(20);
   ((TH1D*)gDirectory->Get(name))->SetStats(0);
   ((TH1D*)gDirectory->Get(name))->DrawCopy("same");
   sprintf(name, "%s","mxallmc");
   if(FITTOTSHAPE) sprintf(name, "%s","mxoneballmc");
   //   ((TH1D*)gDirectory->Get(name))->SetFillStyle(3004);
   ((TH1D*)gDirectory->Get(name))->SetLineColor(38);
   ((TH1D*)gDirectory->Get(name))->SetFillColor(38);
   ((TH1D*)gDirectory->Get(name))->SetStats(0);
   ((TH1D*)gDirectory->Get(name))->DrawCopy("same");  
   ((TH1D*)gDirectory->Get(name))->SetAxisRange(0.,5.);
   sprintf(name, "%s","mxallbkg");
   if(FITTOTSHAPE) sprintf(name, "%s","mxoneballbkg");
   //   ((TH1D*)gDirectory->Get(name))->SetFillStyle(3005);
   ((TH1D*)gDirectory->Get(name))->SetLineColor(kYellow);
   ((TH1D*)gDirectory->Get(name))->SetFillColor(kYellow);
   ((TH1D*)gDirectory->Get(name))->SetStats(0);
   ((TH1D*)gDirectory->Get(name))->DrawCopy("same");
   sprintf(name, "%s","mxscaleoth");
   if(FITTOTSHAPE) sprintf(name, "%s","mxonebscaleoth");
   ((TH1D*)gDirectory->Get(name))->SetLineColor(13);
   ((TH1D*)gDirectory->Get(name))->SetFillColor(13);
   ((TH1D*)gDirectory->Get(name))->SetStats(0);
   ((TH1D*)gDirectory->Get(name))->DrawCopy("same");
   sprintf(name, "%s","mxdata");
   if(FITTOTSHAPE) sprintf(name, "%s","mxonebdata");
   ((TH1D*)gDirectory->Get(name))->DrawCopy("same");
   TLegendEntry *legge; 
   TLegend *leg;
   leg = new TLegend(0.6,0.7,0.88,0.89);
   leg->SetFillStyle(0); leg->SetBorderSize(0.); leg->SetTextSize(0.05); 
   leg->SetFillColor(0); 
   legge = leg->AddEntry(((TH1D*)gDirectory->Get("mxoneballmc")), "b->ulnu", "f"); 
   legge = leg->AddEntry(((TH1D*)gDirectory->Get("mxallbkg")), "b->clnu", "f"); 
   legge = leg->AddEntry(((TH1D*)gDirectory->Get("mxscaleoth")), "other", "f"); 
   legge = leg->AddEntry(((TH1D*)gDirectory->Get("mxonebdata")), "data", "p"); 
   leg->Draw();
   

   c1->cd(2);
   
   sprintf(name, "%s","mxdata");
   themax = ((TH1D*)gDirectory->Get(name))->GetMaximum() * 1.2; 
   ((TH1D*)gDirectory->Get(name))->SetMaximum(themax); 
   ((TH1D*)gDirectory->Get(name))->SetMarkerStyle(20);
   ((TH1D*)gDirectory->Get(name))->SetStats(0);
   ((TH1D*)gDirectory->Get(name))->SetXTitle("Mx(GeV)");
   ((TH1D*)gDirectory->Get(name))->Draw();
   sprintf(name, "%s","mxallmc");
   //   ((TH1D*)gDirectory->Get(name))->SetFillStyle(3004);
   ((TH1D*)gDirectory->Get(name))->SetLineColor(38);
   ((TH1D*)gDirectory->Get(name))->SetFillColor(38);
   ((TH1D*)gDirectory->Get(name))->SetStats(0);
   ((TH1D*)gDirectory->Get(name))->DrawCopy("same");  
   ((TH1D*)gDirectory->Get(name))->SetAxisRange(0.,5.);
   sprintf(name, "%s","mxallbkg");
   //   ((TH1D*)gDirectory->Get(name))->SetFillStyle(3005);
   ((TH1D*)gDirectory->Get(name))->SetLineColor(kYellow);
   ((TH1D*)gDirectory->Get(name))->SetFillColor(kYellow);
   ((TH1D*)gDirectory->Get(name))->SetStats(0);
   ((TH1D*)gDirectory->Get(name))->DrawCopy("same");
   sprintf(name, "%s","mxscaleoth");
   ((TH1D*)gDirectory->Get(name))->SetLineColor(13);
   ((TH1D*)gDirectory->Get(name))->SetFillColor(13);
   ((TH1D*)gDirectory->Get(name))->SetStats(0);
   ((TH1D*)gDirectory->Get(name))->DrawCopy("same");
   sprintf(name, "%s","mxdata");
   ((TH1D*)gDirectory->Get(name))->DrawCopy("same");
   leg->Draw();
   
   c1->cd(3);
   
   sprintf(name, "%s","mxsubdata");
   max = 1.2*((TH1D*)gDirectory->Get(name))->GetMaximum();
   min = 1.2*((TH1D*)gDirectory->Get(name))->GetMinimum();
   ((TH1D*)gDirectory->Get(name))->SetMaximum(max);
   ((TH1D*)gDirectory->Get(name))->SetMinimum(min);
   ((TH1D*)gDirectory->Get(name))->SetMarkerStyle(20);
   ((TH1D*)gDirectory->Get(name))->SetStats(0);
   ((TH1D*)gDirectory->Get(name))->SetXTitle("Mx(GeV)");  
   ((TH1D*)gDirectory->Get(name))->Draw();
   sprintf(name, "%s","mxscalevub");
   //   ((TH1D*)gDirectory->Get(name))->SetFillStyle(3004);
   ((TH1D*)gDirectory->Get(name))->SetLineColor(38);
   ((TH1D*)gDirectory->Get(name))->SetFillColor(38);
   ((TH1D*)gDirectory->Get(name))->SetStats(0);
   ((TH1D*)gDirectory->Get(name))->DrawCopy("same");  
   sprintf(name, "%s","mxsubdata");
   ((TH1D*)gDirectory->Get(name))->DrawCopy("same");
   TLegendEntry *legge2; 
   TLegend *leg2;
   leg2 = new TLegend(0.5,0.7,0.88,0.89);
   leg2->SetFillStyle(0); leg2->SetBorderSize(0.); leg2->SetTextSize(0.05); 
   leg2->SetFillColor(0); 
   legge2 = leg2->AddEntry(((TH1D*)gDirectory->Get("mxscalevub")), "scaled MC", "f"); 
   legge2 = leg2->AddEntry(((TH1D*)gDirectory->Get("mxonebdata")), "data subtr.", "p"); 
   leg2->Draw();
   TLine line(0.,0.,5.,0.);
   line.SetLineColor(kRed); line.SetLineWidth(2); line.SetLineStyle(2); line.Draw();
   
   closeEpsFile();
   
   delete c1;

   c1 = new TCanvas("c1"," ",200,10,800,520); 
   
   sprintf(name, "%s%s%s",DIRNAME.Data(),PREFIXOUT.Data(),"mxspectrum.eps");
   
   openEpsFile(name);
 
   c1->Clear(); 
   c1->Divide(1);
   c1->cd(1);

   sprintf(name, "%s","mxsubdata");
   ((TH1D*)gDirectory->Get(name))->Draw();
   sprintf(name, "%s","mxscalevub");
   ((TH1D*)gDirectory->Get(name))->DrawCopy("same");  
   sprintf(name, "%s","mxsubdata");
   ((TH1D*)gDirectory->Get(name))->DrawCopy("same");
   line.SetLineColor(kRed); line.SetLineWidth(2); line.SetLineStyle(2); line.Draw();   
   leg2->Draw();

   closeEpsFile();


   c1 = new TCanvas("c1"," ",472,0,800,900);
   c1->Clear();
      
   sprintf(name, "%s%s%s",DIRNAME.Data(),PREFIXOUT.Data(),"res.eps");   
   openEpsFile(name);
   ((TH1D*)gDirectory->Get("plotres"))->SetMarkerStyle(3005);
   ((TH1D*)gDirectory->Get("plotres"))->SetLineColor(kRed);
   ((TH1D*)gDirectory->Get("plotres"))->SetFillColor(kRed);
   ((TH1D*)gDirectory->Get("plotres"))->Draw();
   closeEpsFile();
   c1->Clear();

   sprintf(name, "%s%s%s",DIRNAME.Data(),PREFIXOUT.Data(),"nores.eps");   
   openEpsFile(name);
   ((TH1D*)gDirectory->Get("plotnres"))->SetMarkerStyle(3006);
   ((TH1D*)gDirectory->Get("plotnres"))->SetLineColor(kBlue);
   ((TH1D*)gDirectory->Get("plotnres"))->SetFillColor(kBlue);
   ((TH1D*)gDirectory->Get("plotnres"))->Draw();
   closeEpsFile();
   c1->Clear();


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

   //Crosschecktobedeleted
   sprintf(name, "%s%s%s",DIRNAME.Data(),PREFIXOUT.Data(),"ACboth.eps");   
   openEpsFile(name);
   ((TH1D*)gDirectory->Get("ACplotall"))->SetMarkerStyle(3005);
   ((TH1D*)gDirectory->Get("ACplotall"))->SetLineColor(kRed);
   ((TH1D*)gDirectory->Get("ACplotall"))->SetFillColor(kRed);
   ((TH1D*)gDirectory->Get("ACplotall"))->Draw();
   ((TH1D*)gDirectory->Get("ACplotnres"))->SetMarkerStyle(3006);
   ((TH1D*)gDirectory->Get("ACplotnres"))->SetLineColor(kBlue);
   ((TH1D*)gDirectory->Get("ACplotnres"))->SetFillColor(kBlue);
   ((TH1D*)gDirectory->Get("ACplotnres"))->DrawCopy("same");
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

   delete htotal;
   delete c1;       
   
   // calculate ratio of BR

   sprintf(name, "%s","mxonebvub");
   vubmcaftercutsbin1 = ((TH1D*)gDirectory->Get(name))->GetBinContent(1);
   errvubmcaftercutsbin1 = ((TH1D*)gDirectory->Get(name))->GetBinError(1);
   areavubmcaftercutsbin1 = ((TH1D*)gDirectory->Get(name))->Integral();
   
   sprintf(name, "%s","mxonebsubdata");       
   vubaftercutsbin1 = ((TH1D*)gDirectory->Get(name))->GetBinContent(1);
   errvubaftercutsbin1 = ((TH1D*)gDirectory->Get(name))->GetBinError(1);

   if(FITTOTSHAPE){
     sprintf(name, "%s","mxonebdata");
     dataFirstBin= ((TH1D*)gDirectory->Get(name))->GetBinContent(1);
     dataErrFirstBin= ((TH1D*)gDirectory->Get(name))->GetBinError(1);
     sprintf(name, "%s","mxonebvcb");
     areavcbaftercutsbin1 = ((TH1D*)gDirectory->Get(name))->Integral();
     vcbFirstBin= ((TH1D*)gDirectory->Get(name))->GetBinContent(1);
     vcbErrFirstBin= ((TH1D*)gDirectory->Get(name))->GetBinError(1);
     sprintf(name, "%s","mxoneboth");
     areaothaftercutsbin1 = ((TH1D*)gDirectory->Get(name))->Integral();
     othFirstBin= ((TH1D*)gDirectory->Get(name))->GetBinContent(1);
     othErrFirstBin= ((TH1D*)gDirectory->Get(name))->GetBinError(1);
     
     
     vubaftercutsbin1 = dataFirstBin-vcbFirstBin*vcbcomp-othFirstBin*othcomp;
     MCerrvubaftercutsbin1 =sqrt( pow(vcbErrFirstBin*vcbcomp,2)+pow(othErrFirstBin*othcomp,2));
     if(FITTOTSHAPE==2) {
       errvubaftercutsbin1 = sqrt(pow(dataErrFirstBin,2)+pow(vcbFirstBin*errvcbcompNOMC,2)+pow(othFirstBin*errothcompNOMC,2));
       MCerrvubaftercutsbin1 = sqrt(
				    pow(MCerrvubaftercutsbin1,2)
				    +pow(vcbFirstBin*errvcbcomp,2)+pow(othFirstBin*errothcomp,2)
				    -pow(vcbFirstBin*errvcbcompNOMC,2)-pow(othFirstBin*errothcompNOMC,2));
     }else{
       errvubaftercutsbin1 = sqrt(pow(dataErrFirstBin,2)+pow(vcbFirstBin*errvcbcomp,2)+pow(othFirstBin*errothcomp,2));
     }
   }
   
   if(BLINDING){
     vubaftercutsbin1 *=  blindfactor;
     errvubaftercutsbin1 *=  blindfactor;
     MCerrvubaftercutsbin1  *=  blindfactor;
   }else{
     blindfactor = 1;
   }

   S = vubaftercutsbin1;
   errS = errvubaftercutsbin1;
   Bvcb = vcbaftercutsbin1;
   temperrvcb = sqrt(errvcbaftercutsbin1 * errvcbaftercutsbin1 + errfitvcbaftercutsbin1 * errfitvcbaftercutsbin1);
   errBvcb = temperrvcb;
   Bother = otheraftercutsbin1;
   temperroth = sqrt(errotheraftercutsbin1 * errotheraftercutsbin1 + errfitotheraftercutsbin1 * errfitotheraftercutsbin1);
   errBother = temperroth;
   B = Bvcb + Bother;
   errB = sqrt(errBvcb * errBvcb + errBother * errBother);
   //Various eff. calculations
   epsu = vubmcallforeff/vubmcleptforeff;
   errepsu = sqrt(epsu*(1 - epsu)/vubmcleptforeff) ;
   // this is the old code, if we will use a different MC for vub efficiecies we will uncomment
//    double epsmx = vubmcaftercutsbin1/vubmcaftercuts;    
//    double errepsmx = sqrt(epsmx*(1 - epsmx)/vubmcaftercuts) ;   
   epsmx = vubmcaftercutsbin1/vubmcallforeff;
   cout<<"Normalization/eff:: "<<vubmcaftercutsbin1<<" "<<vubmcallforeff<<" "<<vubmcleptforeff<<endl;
   errepsmx = sqrt(epsmx*(1 - epsmx)/vubmcallforeff) ;   
   epstot = epsu * epsmx;
   errepstot = sqrt(epsu*epsu*errepsmx*errepsmx + epsmx*epsmx*errepsu*errepsu);
   
   fact = vcbmc * (1 + BRRATIOVALUETAIL/0.104) / (totmc + vcbmc * BRRATIOVALUETAIL/0.104); //factor to take into account the vub contribute into all events after lept cut
   //BRRATIOVALUETAIL for BR(b->ulnu) is assumed
   
   
   //now pstarfact calculated in the fit	 
   BRBR = S * (1/(tot*fact*calcpstarfact)) / (epsu * epsmx);
   errBRBR = errS * (1/(tot*fact*calcpstarfact)) / (epsu * epsmx);
   errBRBRMCstat=BRBR*
     sqrt(pow(MCerrvubaftercutsbin1/S,2)+pow(errepstot/epstot,2));
   errBRBRtheo = (2.1002 -1.0995 * mxcut +  0.29050 * mxcut*mxcut);
   errBRBRtheo = (errBRBRtheo-1.) * 2. * BRBR ;
   errtotalBRBR = sqrt(errBRBR * errBRBR + errBRBRtheo * errBRBRtheo);
   
   // print results on the screen
   cout << endl;
   cout << endl;
   cout << endl;
   cout << "BRBR = " << BRBR << " +- " << errBRBR << "(stat) +- " << errBRBRMCstat << "(MC stat)" << endl;
   cout << endl;
   cout << endl;
   cout << endl;
   
   // print results on file
   sprintf(name, "%s%s%s",DIRNAME.Data(),PREFIXOUT.Data(),"result_nocat.dat");        
   ofstream outfile(name);
   
   outfile << "#NUMBERS FOR THE SCANS" << endl;
   outfile << endl;
   outfile << endl;
   outfile << " pstarfact           "  << PSTARFACT   <<  endl;  
   outfile << " mxCut               "  << MXCUT   <<  endl;  
   outfile << " q2Cut               "  << Q2CUT   <<  endl;  
   outfile << " csiloCut            "  << CSILOCUT   <<  endl;  
   outfile << " csihiCut            "  << CSIHICUT   <<  endl;  
   outfile << " xloCut              "  << XLOCUT   <<  endl;  
   outfile << " xhiCut              "  << XHICUT   <<  endl;  
   outfile << " wloCut              "  << WLOCUT   <<  endl;  
   outfile << " whiCut              "  << WHICUT   <<  endl;  
   outfile << " q2loCut             "  << Q2LOCUT   <<  endl;  
   outfile << " q2hiCut             "  << Q2HICUT   <<  endl;  
   outfile << " ewpwloCut           "  << EWPWLOCUT   <<  endl;  
   outfile << " ewpwhiCut           "  << EWPWHICUT   <<  endl;  
   outfile << " leptonPCut          "  << LEPTONPCUT   <<  endl;  
   outfile << " prmm2cut            "  << PRMM2CUT   <<  endl;  
   outfile << " mnuSqLow            "  << MNUSQLOW   <<  endl;  
   outfile << " mnuSqHigh           "  << MNUSQHIGH   <<  endl;  
   outfile << " chLow               "  << CHLOW   <<  endl;  
   outfile << " chHigh              "  << CHHIGH   <<  endl;  
   outfile << " depl                "  << DEPL   <<  endl;  
   outfile << " Btype               "  << BTYPE   <<  endl;  
   outfile << " lepttype            "  << LEPTTYPE   <<  endl;  
   outfile << " minintpur           "  << MININTPUR    <<  endl;
   outfile << " maxintpur           "  << MAXINTPUR    <<  endl;
   outfile << " nnpi0               "  << CUTNNPI0  <<  endl;
   outfile << " vubcomp             "  << vubcomp << "  " << errvubcomp << endl;
   outfile << " vcbcomp             "  << vcbcomp << "  " << errvcbcomp << endl;
   outfile << " othcomp             "  << othcomp << "  " << errothcomp << endl;
   outfile << " data1bin            "  << dataFirstBin << "  " << dataErrFirstBin <<  "  " << endl;
   outfile << " vub1bin             "  << S/blindfactor << "  " << errS/blindfactor <<  "  " << ((TH1D*)gDirectory->Get("mxonebvub"))->GetBinError(1)*vubcomp/blindfactor << endl;
   outfile << " vcb1bin             "  << vcbaftercutsbin1 << "  " <<  errfitvcbaftercutsbin1 << "  " << errvcbaftercutsbin1 << endl;
   outfile << " oth1bin             "  << otheraftercutsbin1 << "  " << errfitotheraftercutsbin1 <<  "  " << errotheraftercutsbin1 << endl;
   outfile << " epsu                " << epsu << "  " << errepsu << endl;
   outfile << " epsmx               " << epsmx << "  " << errepsmx << endl;
   outfile << " epstot              " << epstot << "  " << errepstot << endl;
   outfile << " nsl                 " << tot << endl;
   outfile << " nslmc               " << totmc << endl;
   outfile << " fact                " << fact << endl;
   outfile << " pstarfact           " << calcpstarfact << endl;
   outfile << " BRBR                " << BRBR << "  " << errBRBR << "  " << errBRBRMCstat << endl;
   outfile << " chisq               " << chisq << endl;
   


   outfile << endl;
   outfile << endl;
   outfile << "INPUT FILES" << endl;
   outfile << endl;
   outfile << endl;
   outfile << " fileVubTotal        "  << FILEVUBTOTAL   <<  endl;  
   outfile << " fileVcb             "  << FILEVCB   <<  endl;  
   outfile << " fileData            "  << FILEDATA   <<  endl;  
   
   outfile << endl;
   outfile << endl;
   outfile << "ALL NUMBERS" << endl;
   outfile << endl;
   outfile << "MX FIT" << endl;
   outfile << endl;
   outfile << "Vub comp = " << vubcomp << " +- " << errvubcomp << endl;
   outfile << "Vcb comp = " << vcbcomp << " +- " << errvcbcomp << endl;
   outfile << "oth comp = " << othcomp << " +- " << errothcomp << endl;
   outfile << endl;
   outfile << "Vub tot area = " << ((TH1D*)gDirectory->Get("mxonebvub"))->Integral() << endl;
   outfile << "Vcb tot area = " << ((TH1D*)gDirectory->Get("mxonebvcb"))->Integral() << endl;
   outfile << "oth tot area = " << ((TH1D*)gDirectory->Get("mxoneboth"))->Integral() << endl;
   outfile << endl;
   outfile << "data 1' bin = " << ((TH1D*)gDirectory->Get("mxonebdata"))->GetBinContent(1) << " +- " << ((TH1D*)gDirectory->Get("mxonebdata"))->GetBinError(1) << endl;
   outfile << "Vub 1' bin = " << ((TH1D*)gDirectory->Get("mxonebvub"))->GetBinContent(1) << " +- " << ((TH1D*)gDirectory->Get("mxonebvub"))->GetBinError(1) << endl;
   outfile << "Vcb 1' bin = " << ((TH1D*)gDirectory->Get("mxonebvcb"))->GetBinContent(1) << " +- " << ((TH1D*)gDirectory->Get("mxonebvcb"))->GetBinError(1) << endl;
   outfile << "oth 1' bin = " << ((TH1D*)gDirectory->Get("mxoneboth"))->GetBinContent(1) << " +- " << ((TH1D*)gDirectory->Get("mxoneboth"))->GetBinError(1) << endl;
   outfile << endl;
   outfile << "WEIGHTED NUMBERS" << endl;
   outfile << endl;
   outfile << "Vub fitted 1' bin = " << S/blindfactor << " +- " <<  ((TH1D*)gDirectory->Get("mxonebvub"))->GetBinError(1)*vubcomp/blindfactor << "(stat MC) +- " << errS/blindfactor << "(err fit)" << endl;
   outfile << "Vcb fitted 1' bin = " << vcbaftercutsbin1 << " +- " << errvcbaftercutsbin1 << "(stat MC) +- " << errfitvcbaftercutsbin1 << "(err fit)" << endl;
   outfile << "oth fitted 1' bin = " << otheraftercutsbin1 << " +- " << errotheraftercutsbin1 << "(stat MC) +- " << errfitotheraftercutsbin1 << "(err fit)"  << endl;
   outfile << endl;
   outfile << "Vub 1' true = " << vubmcselected << endl;
   outfile << "Vcb 1' true = " << vcbmcselected << endl;
   outfile << "oth 1' true = " << othermcselected << endl;
   outfile << endl;
   outfile << "EFFICIENCY Vub" << endl;
   outfile << endl;
   outfile << "Vub total MC (lepton cut) = " << vubmc << endl;
   outfile << "Vub MC (all cuts) = " << vubmcaftercuts << endl;
   outfile << "Vub MC (all cuts + Mx cut) = " << vubmcaftercutsbin1 << endl;
   outfile << "Vub gene total MC (lepton cut) = " << vubmcleptforeff << endl;
   outfile << "Vub gene MC (all cuts) = " << vubmcallforeff << endl;
   outfile << "Eps_u =  " << epsu << " +- " << errepsu << endl;
   outfile << "Eps_Mx = " << epsmx << " +- " << errepsmx << endl;
   outfile << "Eps_tot = " << epstot << " +- " << errepstot << endl;
   outfile << endl; 
   outfile << endl;
   outfile << "NSL NUMBERS" << endl;
   outfile << endl;
   outfile << "Nsl = " << tot << endl;
   outfile << endl;
   outfile << "Nsl MC = " << totmc << endl;
   outfile << "Nsl - BG = " << vcbmc << endl;
   outfile << "(Nsl - BG)/Nsl (vub corrected) = " << fact << endl;
   outfile << endl;
   outfile << "Pstar fact = " << PSTARFACT << endl;
   outfile << endl;
   outfile << endl; 
   outfile << "Eff on Vcb (Nvcb/(Nsl-BGsl)) = " << ((TH1D*)gDirectory->Get("mxonebvcb"))->Integral() * vcbcomp / (tot * fact) << " +- " <<  ((TH1D*)gDirectory->Get("mxonebvcb"))->Integral() * errvcbcomp / (tot * fact) << endl;  
   outfile << "Eff on other (other/(Nsl-BGsl)) = " << ((TH1D*)gDirectory->Get("mxoneboth"))->Integral() * othcomp / (tot * fact) << " +- " << ((TH1D*)gDirectory->Get("mxoneboth"))->Integral() * errothcomp / (tot * fact)  << endl;  
   outfile << endl; 
   outfile << endl;
   outfile << "Vub gene MC (no cuts) = " << vubmcnocut << endl;
   outfile << "Vcb gene MC (no cuts) = " << vcbmcnocut << endl;
   outfile << "Oth gene MC (no cuts) = " << othmcnocut << endl;
   outfile << "Recalculated Pstar fact = " << calcpstarfacttemp << endl;
   outfile << endl;
   outfile << endl; 
   outfile << "Nsig = " << S/epstot  << " +- " <<  errS/epstot   << endl;

   outfile << endl;
   outfile << endl; 
   outfile << "BRBR = " << BRBR << " +- " << errBRBR << "(stat) +- " << errBRBRMCstat << "(MC stat)" << endl;
   outfile << endl;
   outfile << "Chi Square of the Fit = " << chisq << endl;
   outfile << "NDOF = " << NDOF-3 << endl;
   outfile << endl;
   outfile << endl;
   double Vub = 0.00445 * sqrt((BRBR * .104 * 1.55) / (0.002 * 1.622));
   outfile << "Vub(*10-3) = " << Vub*1000 << " +- " << Vub*1000*(errBRBR/(2*BRBR)) << "(stat) +- " << Vub*1000*errBRBRMCstat/(2*BRBR) << "(MC stat)" << endl;
   outfile << endl;
   outfile << endl;
   outfile << endl;
   outfile << endl;         
   outfile.close();

   sprintf(name, "%s%s%s.tex",DIRNAME.Data(),PREFIXOUT.Data(),texPrefix.Data());        
   cout << "**********************************************************************" << endl;
   cout << "DUMPING TEXFILE 3 " << name << endl;
   cout << "**********************************************************************" << endl;
   ofstream texfile(name);
   char texline[200];
   sprintf(texline, "\\def\\%ssemil{%7.0f \\pm %6.0f} ", texPrefix.Data(), tot, errtot); 
   texfile << texline << endl;
   sprintf(texline, "\\def\\%sbgsl{%5.0f \\pm %4.0f} ", texPrefix.Data(), (1.-fact)*tot, (1.-fact)*errtot); 
   texfile << texline << endl;
   sprintf(texline, "\\def\\%sslsub{%5.0f \\pm %4.0f} ", texPrefix.Data(), tot - (1.-fact)*tot, errtot); 
   texfile << texline << endl;
   sprintf(texline, "\\def\\%sepssl{%5.2f \\pm %4.2f} ", texPrefix.Data(), calcpstarfact, errcalcpstarfact); 
   texfile << texline << endl;
   sprintf(texline, "\\def\\%snu{%5.0f \\pm %4.0f} ", texPrefix.Data(), S, errS); 
   texfile << texline << endl;
   sprintf(texline, "\\def\\%sbgc{%5.0f \\pm %4.0f} ", texPrefix.Data(), vcbaftercutsbin1, errfitvcbaftercutsbin1); 
   texfile << texline << endl;
   sprintf(texline, "\\def\\%sbgo{%5.0f \\pm %4.0f} ", texPrefix.Data(), otheraftercutsbin1, errfitotheraftercutsbin1); 
   texfile << texline << endl;
   sprintf(texline, "\\def\\%sepsu{%4.3f} ", texPrefix.Data(), epsu); 
   texfile << texline << endl;
   sprintf(texline, "\\def\\%sepsmx{%4.3f} ", texPrefix.Data(), epsmx); 
   texfile << texline << endl;
   sprintf(texline, "\\def\\%sbrbr{%6.4f \\pm %5.4f} ", texPrefix.Data(), BRBR, errBRBR); 
   texfile << texline << endl;
   sprintf(texline, "\\def\\%sbrbrerrmc{%5.4f} ", texPrefix.Data(), errBRBRMCstat); 
   texfile << texline << endl;
   texfile.close(); 

} 

// ----------------------------------------------------------------------
Double_t fitNtp::smeargauss(double invalue, double mean, double sigma){    

  // smears using a gaussian distribution

  double returnvalue=0.;
  f1all->SetParameters(1.,0.,sigma);
  returnvalue = invalue + mean + f1all->GetRandom(); 
  //  ((TH1D*)gDirectory->Get("h8888"))->Fill(returnvalue - invalue);
  return returnvalue;

}

void fitNtp::theq2fit(){    

  Double_t vubmcaftercuts; char name[200], namet[200];
  double tempbin = 0;
  double temperr = 0;
  double tempbinvub = 0;
  double tempbinvcb = 0;
  double tempbinoth = 0;
  double temperrvcb = 0;
  double temperroth = 0;
  double tempbinerrvub = 0;
  double tempbinerrvcb = 0;
  double tempbinerroth = 0;

  //true number of vub MC events after all cuts
  sprintf(name, "%s","q2vub");	   
  vubmcaftercuts = ((TH1D*)gDirectory->Get(name))->Integral();
  
  TH1D* htotal;
  if(FITTOTSHAPE==2){
    htotal = fitWithErrors(-1,-1,"q2");
  } else {
    htotal = fitWithoutErrors();
  }
   // chisq calculation 
  chisq = 0; 
  double tempchisq = 0;
  
  NDOF = 0;
  cout<<"************ Q2 fit *********:: "<<vubcomp<<" "<<vcbcomp<<" "<<othcomp<<endl;

  cout << "CHISQ CALCULATION for Q2 dist" << endl;
  for(int i=1;i<11;i++){
    tempchisq = 0;
    sprintf(name, "%s%s","q2","onebdata");	         
    if(((TH1D*)gDirectory->Get(name))->GetBinContent(i)) {
      tempbin = ((TH1D*)gDirectory->Get(name))->GetBinContent(i);
      temperr = ((TH1D*)gDirectory->Get(name))->GetBinError(i);
      sprintf(name, "%s%s","q2","onebvub");	         
      tempbinvub = ((TH1D*)gDirectory->Get(name))->GetBinContent(i) * vubcomp;
      tempbinerrvub = ((TH1D*)gDirectory->Get(name))->GetBinError(i) * vubcomp;
      sprintf(name, "%s%s","q2","onebvcb");	         
      tempbinvcb = ((TH1D*)gDirectory->Get(name))->GetBinContent(i) * vcbcomp;
      tempbinerrvcb = ((TH1D*)gDirectory->Get(name))->GetBinError(i) * vcbcomp;
      sprintf(name, "%s%s","q2","oneboth");	     
      tempbinoth = ((TH1D*)gDirectory->Get(name))->GetBinContent(i) * othcomp;
      tempbinerroth = ((TH1D*)gDirectory->Get(name))->GetBinError(i) * othcomp;
      temperr = sqrt(temperr*temperr + tempbinerrvub*tempbinerrvub + tempbinerrvcb * tempbinerrvcb + tempbinerroth*tempbinerroth);
      if(temperr<1) temperr = 1; 
      sprintf(name, "%s%s","q2","onebdata");	
      tempchisq = (tempbinvub + tempbinvcb + tempbinoth - ((TH1D*)gDirectory->Get(name))->GetBinContent(i))/temperr;     
      NDOF++;
    }   
    chisq = chisq + tempchisq*tempchisq;
  }
  
  chisq = chisq / (NDOF-3);
  
  cout << " Chi Square of the Q2 Fit :" << chisq << endl;
  cout << "NDOF " << NDOF-3 << endl;

//   ///TOBECHECKED
//    sprintf(name, "%s","mxonebvcb");	         
//    vcbaftercutsbin1 = ((TH1D*)gDirectory->Get(name))->GetBinContent(1) * vcbcomp;
//    errvcbaftercutsbin1 = ((TH1D*)gDirectory->Get(name))->GetBinError(1) * vcbcomp;
//    errfitvcbaftercutsbin1 = ((TH1D*)gDirectory->Get(name))->GetBinContent(1) * errvcbcomp;
//    sprintf(name, "%s","mxoneboth");	         
//    otheraftercutsbin1 = ((TH1D*)gDirectory->Get(name))->GetBinContent(1) * othcomp;
//    errotheraftercutsbin1 = ((TH1D*)gDirectory->Get(name))->GetBinError(1) * othcomp;
//    errfitotheraftercutsbin1 = ((TH1D*)gDirectory->Get(name))->GetBinContent(1) * errothcomp;


  // background subtraction q2 spectrum

   for(int ik=1;ik<11;ik++){
     sprintf(name, "%s","q2data");	         
     tempbin = ((TH1D*)gDirectory->Get(name))->GetBinContent(ik);
     temperr = ((TH1D*)gDirectory->Get(name))->GetBinError(ik);
     sprintf(name, "%s","q2vcb");	         
     tempbinvcb = ((TH1D*)gDirectory->Get(name))->GetBinContent(ik) * vcbcomp;
     tempbinerrvcb = ((TH1D*)gDirectory->Get(name))->GetBinError(ik) * vcbcomp;
     temperrvcb = tempbinvcb * errvcbcomp / vcbcomp;
     //temperrvcb = sqrt(temperrvcb * temperrvcb + tempbinerrvcb * tempbinerrvcb);     //THE MC STAT IS NOT INCLUDED 
     sprintf(name, "%s","q2oth");	         
     tempbinoth = ((TH1D*)gDirectory->Get(name))->GetBinContent(ik) * othcomp;
     tempbinerroth = ((TH1D*)gDirectory->Get(name))->GetBinError(ik) * othcomp;
     temperroth = tempbinoth * errothcomp / othcomp;     
     //temperroth = sqrt(temperroth * temperroth + tempbinerroth * tempbinerroth);     //THE MC STAT IS NOT INCLUDED 
     tempbin = tempbin - tempbinvcb - tempbinoth;
     temperr = sqrt(temperr*temperr + tempbinerrvcb*tempbinerrvcb + tempbinerroth*tempbinerroth);
     sprintf(name, "%s","q2subdata");	 
     ((TH1D*)gDirectory->Get(name))->SetBinContent(ik, tempbin);
     ((TH1D*)gDirectory->Get(name))->SetBinError(ik, temperr); 
   }	
   int y; double it; char name2[200];
   sprintf(name, "%s","q2scalevcb");
   sprintf(name2, "%s","q2vcb");   
   for(y=1;y<11;y++){
     ((TH1D*)gDirectory->Get(name))->SetBinContent(y,((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*vcbcomp);
   }
   
   sprintf(name, "%s","q2scaleoth");
   sprintf(name2, "%s","q2oth");
   for(y=1;y<11;y++){
     ((TH1D*)gDirectory->Get(name))->SetBinContent(y,((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*othcomp);
   }
   
   sprintf(name, "%s","q2scalevub");
   sprintf(name2, "%s","q2vub");
   for(y=1;y<11;y++){
     ((TH1D*)gDirectory->Get(name))->SetBinContent(y,((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*vubcomp);
   }
   
   sprintf(name, "%s","q2allbkg");
   sprintf(name2, "%s","q2vcb");
   for(y=1;y<11;y++){
     ((TH1D*)gDirectory->Get(name))->SetBinContent(y,((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*vcbcomp);
   }
   sprintf(name2, "%s","q2oth");
   for(y=1;y<11;y++){
     it = ((TH1D*)gDirectory->Get(name))->GetBinContent(y);
     ((TH1D*)gDirectory->Get(name))->SetBinContent(y,it + ((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*othcomp);
   }

   sprintf(name, "%s","q2allmc");
   sprintf(name2, "%s","q2vcb");
   for(y=1;y<11;y++){
     ((TH1D*)gDirectory->Get(name))->SetBinContent(y,((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*vcbcomp);
   }
   sprintf(name2, "%s","q2oth");
   for(y=1;y<11;y++){
     it = ((TH1D*)gDirectory->Get(name))->GetBinContent(y);
     ((TH1D*)gDirectory->Get(name))->SetBinContent(y,it + ((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*othcomp);
   }
   sprintf(name2, "%s","q2vub");
   for(y=1;y<11;y++){
     it = ((TH1D*)gDirectory->Get(name))->GetBinContent(y);
     ((TH1D*)gDirectory->Get(name))->SetBinContent(y,it + ((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*vubcomp);
   }

//    NDOF = 0;
//    newchisq = 0; 
//    newtempchisq = 0;
//    for(int i=1;i<11;i++){
//      tempchisq = 0;
//      sprintf(name, "%s","q2data");	         
//      if(((TH1D*)gDirectory->Get(name))->GetBinContent(i)) {
//        tempbin = ((TH1D*)gDirectory->Get(name))->GetBinContent(i);
//        temperr = ((TH1D*)gDirectory->Get(name))->GetBinError(i);
//        sprintf(name, "%s","q2vub");	         
//        tempbinvub = ((TH1D*)gDirectory->Get(name))->GetBinContent(i) * vubcomp;
//        tempbinerrvub = ((TH1D*)gDirectory->Get(name))->GetBinError(i) * vubcomp;
//        sprintf(name, "%s","q2vcb");	         
//        tempbinvcb = ((TH1D*)gDirectory->Get(name))->GetBinContent(i) * vcbcomp;
//        tempbinerrvcb = ((TH1D*)gDirectory->Get(name))->GetBinError(i) * vcbcomp;
//        sprintf(name, "%s","q2oth");	     
//        tempbinoth = ((TH1D*)gDirectory->Get(name))->GetBinContent(i) * othcomp;
//        tempbinerroth = ((TH1D*)gDirectory->Get(name))->GetBinError(i) * othcomp;
//        temperr = sqrt(temperr*temperr + tempbinerrvub*tempbinerrvub + tempbinerrvcb*tempbinerrvcb + tempbinerroth*tempbinerroth);
//        if(temperr<1) temperr = 1; 
//        sprintf(name, "%s","q2data");	
//        newtempchisq = (tempbinvub + tempbinvcb + tempbinoth - ((TH1D*)gDirectory->Get(name))->GetBinContent(i))/temperr;     
//        NDOF++;
//      }   
//      newchisq = newchisq + newtempchisq*newtempchisq;
//    }
   
//    newchisq = newchisq / (NDOF-3);
   
//    cout << "NEW (lepton) Chi Square of the Fit :" << newchisq << endl;
//    cout << "NDOF " << NDOF-1 << endl;

//    //Just cut and paste some histos

//    // draw q2 spectrum

//    c1 = new TCanvas("c1"," ",200,10,1000,520); 
   
//    sprintf(name, "%s%s%s",DIRNAME.Data(),PREFIXOUT.Data(),"q2spectrum.eps");
   
//    openEpsFile(name);
   
   
//    c1->Clear();

//    c1->Divide(2, 1);
//    c1->cd(1);
   
//    sprintf(name, "%s","q2data");
//     ((TH1D*)gDirectory->Get(name))->SetMarkerSize(1.);	
//    ((TH1D*)gDirectory->Get(name))->SetMarkerColor(kBlack);	
//    ((TH1D*)gDirectory->Get(name))->SetStats(0);
//    ((TH1D*)gDirectory->Get(name))->SetXTitle("Q^2(GeV)");
//    themax = ((TH1D*)gDirectory->Get(name))->GetMaximum() * 1.2; 
//    ((TH1D*)gDirectory->Get(name))->SetMaximum(themax); 
//    ((TH1D*)gDirectory->Get(name))->SetStats(0);
//    ((TH1D*)gDirectory->Get(name))->Draw();
//    sprintf(name, "%s","q2data");
//    ((TH1D*)gDirectory->Get(name))->SetMarkerStyle(20);
//    ((TH1D*)gDirectory->Get(name))->SetStats(0);
//    ((TH1D*)gDirectory->Get(name))->DrawCopy("same");
//    sprintf(name, "%s","q2allmc");
//    ((TH1D*)gDirectory->Get(name))->SetLineColor(38);
//    ((TH1D*)gDirectory->Get(name))->SetFillColor(38);
//    ((TH1D*)gDirectory->Get(name))->SetStats(0);
//    ((TH1D*)gDirectory->Get(name))->DrawCopy("same");  
//    ((TH1D*)gDirectory->Get(name))->SetAxisRange(0.,5.);
//    sprintf(name, "%s","q2allbkg");
//    ((TH1D*)gDirectory->Get(name))->SetLineColor(kYellow);
//    ((TH1D*)gDirectory->Get(name))->SetFillColor(kYellow);
//    ((TH1D*)gDirectory->Get(name))->SetStats(0);
//    ((TH1D*)gDirectory->Get(name))->DrawCopy("same");
//    sprintf(name, "%s","q2scaleoth");
//    ((TH1D*)gDirectory->Get(name))->SetLineColor(13);
//    ((TH1D*)gDirectory->Get(name))->SetFillColor(13);
//    ((TH1D*)gDirectory->Get(name))->SetStats(0);
//    ((TH1D*)gDirectory->Get(name))->DrawCopy("same");
//    sprintf(name, "%s","q2data");
//    ((TH1D*)gDirectory->Get(name))->DrawCopy("same");
//    TLegendEntry *legge3q2; 
//    TLegend *leg3q2;
//    leg3q2 = new TLegend(0.1,0.7,0.88,0.89);
//    leg3q2->SetFillStyle(0); leg3q2->SetBorderSize(0.); leg3q2->SetTextSize(0.05); 
//    leg3q2->SetFillColor(0); 
//    legge3q2 = leg3q2->AddEntry(((TH1D*)gDirectory->Get("q2allmc")), "b->ulnu", "f"); 
//    legge3q2 = leg3q2->AddEntry(((TH1D*)gDirectory->Get("q2allbkg")), "b->clnu", "f"); 
//    legge3q2 = leg3q2->AddEntry(((TH1D*)gDirectory->Get("q2scaleoth")), "other", "f"); 
//    legge3q2 = leg3q2->AddEntry(((TH1D*)gDirectory->Get("q2data")), "data", "p"); 
//    leg3q2->Draw();

//    c1->cd(2);
   
//    sprintf(name, "%s","q2subdata");
//    max = 1.2*((TH1D*)gDirectory->Get(name))->GetMaximum();
//    min = 1.2*((TH1D*)gDirectory->Get(name))->GetMinimum();
//    ((TH1D*)gDirectory->Get(name))->SetMaximum(max);
//    ((TH1D*)gDirectory->Get(name))->SetMinimum(0.);
//    ((TH1D*)gDirectory->Get(name))->SetMarkerStyle(20);
//    ((TH1D*)gDirectory->Get(name))->SetStats(0);
//    ((TH1D*)gDirectory->Get(name))->SetXTitle("Q^2(GeV)");
//    ((TH1D*)gDirectory->Get(name))->Draw();
//    sprintf(name, "%s","q2scalevub");
//    ((TH1D*)gDirectory->Get(name))->SetLineColor(38);
//    ((TH1D*)gDirectory->Get(name))->SetFillColor(38);
//    ((TH1D*)gDirectory->Get(name))->SetStats(0);
//    ((TH1D*)gDirectory->Get(name))->DrawCopy("same");  
//    sprintf(name, "%s","q2subdata");
//    ((TH1D*)gDirectory->Get(name))->DrawCopy("same");  
//    TLegendEntry *legge4q2; 
//    TLegend *leg4q2;
//    leg4q2 = new TLegend(0.1,0.7,0.88,0.89);
//    leg4q2->SetFillStyle(0); leg4q2->SetBorderSize(0.); leg4q2->SetTextSize(0.05); 
//    leg4q2->SetFillColor(0); 
//    legge4q2 = leg4q2->AddEntry(((TH1D*)gDirectory->Get("q2scalevub")), "scaled MC", "f"); 
//    legge4q2 = leg4q2->AddEntry(((TH1D*)gDirectory->Get("q2subdata")), "data subt.", "p"); 
//    leg4q2->Draw();

   
//    closeEpsFile();

//    //Dumping some of the information
//   // calculate ratio of BR

//    sprintf(name, "%s","mxonebvub");
//    vubmcaftercutsbin1 = ((TH1D*)gDirectory->Get(name))->GetBinContent(1);
   
//    sprintf(name, "%s","mxonebsubdata");       
//    vubaftercutsbin1 = ((TH1D*)gDirectory->Get(name))->GetBinContent(1);
//    errvubaftercutsbin1 = ((TH1D*)gDirectory->Get(name))->GetBinError(1);

//    if(FITTOTSHAPE){
//      sprintf(name, "%s","mxonebdata");
//      dataFirstBin= ((TH1D*)gDirectory->Get(name))->GetBinContent(1);
//      dataErrFirstBin= ((TH1D*)gDirectory->Get(name))->GetBinError(1);
//      sprintf(name, "%s","mxonebvcb");
//      double vcbFirstBin= ((TH1D*)gDirectory->Get(name))->GetBinContent(1);
//      double vcbErrFirstBin= ((TH1D*)gDirectory->Get(name))->GetBinError(1);
//      sprintf(name, "%s","mxoneboth");
//      double othFirstBin= ((TH1D*)gDirectory->Get(name))->GetBinContent(1);
//      double othErrFirstBin= ((TH1D*)gDirectory->Get(name))->GetBinError(1);
     
     
//      vubaftercutsbin1 = dataFirstBin-vcbFirstBin*vcbcomp-othFirstBin*othcomp;
//      MCerrvubaftercutsbin1 =sqrt( pow(vcbErrFirstBin*vcbcomp,2)+pow(othErrFirstBin*othcomp,2));
//      if(FITTOTSHAPE==2) {
//        errvubaftercutsbin1 = sqrt(pow(dataErrFirstBin,2)+pow(vcbFirstBin*errvcbcompNOMC,2)+pow(othFirstBin*errothcompNOMC,2));
//        MCerrvubaftercutsbin1 = sqrt(
// 				    pow(MCerrvubaftercutsbin1,2)
// 				    +pow(vcbFirstBin*errvcbcomp,2)+pow(othFirstBin*errothcomp,2)
// 				    -pow(vcbFirstBin*errvcbcompNOMC,2)-pow(othFirstBin*errothcompNOMC,2));
//      }else{
//        errvubaftercutsbin1 = sqrt(pow(dataErrFirstBin,2)+pow(vcbFirstBin*errvcbcomp,2)+pow(othFirstBin*errothcomp,2));
//      }
//    }
   
//    if(BLINDING){
//      vubaftercutsbin1 *=  blindfactor;
//      errvubaftercutsbin1 *=  blindfactor;
//      MCerrvubaftercutsbin1  *=  blindfactor;
//    }else{
//      blindfactor = 1;
//    }

//    S = vubaftercutsbin1;
//    errS = errvubaftercutsbin1;
//    Bvcb = vcbaftercutsbin1;
//    temperrvcb = sqrt(errvcbaftercutsbin1 * errvcbaftercutsbin1 + errfitvcbaftercutsbin1 * errfitvcbaftercutsbin1);
//    errBvcb = temperrvcb;
//    Bother = otheraftercutsbin1;
//    temperroth = sqrt(errotheraftercutsbin1 * errotheraftercutsbin1 + errfitotheraftercutsbin1 * errfitotheraftercutsbin1);
//    errBother = temperroth;
//    B = Bvcb + Bother;
//    errB = sqrt(errBvcb * errBvcb + errBother * errBother);
   
//    double epsu = vubmcallforeff/vubmcleptforeff;
//    double errepsu = sqrt(epsu*(1 - epsu)/vubmcleptforeff) ;
//    // this is the old code, if we will use a different MC for vub efficiecies we will uncomment
// //    double epsmx = vubmcaftercutsbin1/vubmcaftercuts;    
// //    double errepsmx = sqrt(epsmx*(1 - epsmx)/vubmcaftercuts) ;   
//    double epsmx = vubmcaftercutsbin1/vubmcallforeff;
//    double errepsmx = sqrt(epsmx*(1 - epsmx)/vubmcallforeff) ;   
//    double epstot = epsu * epsmx;
//    double errepstot = sqrt(epsu*epsu*errepsmx*errepsmx + epsmx*epsmx*errepsu*errepsu);
   
//    fact = vcbmc * (1 + BRRATIOVALUETAIL/0.104) / (totmc + vcbmc * BRRATIOVALUETAIL/0.104); //factor to take into account the vub contribute into all events after lept cut
//    //BRRATIOVALUETAIL for BR(b->ulnu) is assumed
   
   
//    //now pstarfact calculated in the fit	 
//    BRBR = S * (1/(tot*fact*calcpstarfact)) / (epsu * epsmx);
//    errBRBR = errS * (1/(tot*fact*calcpstarfact)) / (epsu * epsmx);
//    errBRBRMCstat=BRBR*
//      sqrt(pow(MCerrvubaftercutsbin1/S,2)+pow(errepstot/epstot,2));
//    errBRBRtheo = (2.1002 -1.0995 * mxcut +  0.29050 * mxcut*mxcut);
//    errBRBRtheo = (errBRBRtheo-1.) * 2. * BRBR ;
//    errtotalBRBR = sqrt(errBRBR * errBRBR + errBRBRtheo * errBRBRtheo);
   
//    // print results on the screen
//    cout << endl;
//    cout << endl;
//    cout << endl;
//    cout << "BRBR = " << BRBR << " +- " << errBRBR << "(stat) +- " << errBRBRMCstat << "(MC stat)" << endl;
//    cout << endl;
//    cout << endl;
//    cout << endl;
   // print results on file
   sprintf(name, "%s%s%s",DIRNAME.Data(),PREFIXOUT.Data(),"result_nocat.dat");        
   sprintf(namet, "%s%s%s.tex",DIRNAME.Data(),PREFIXOUT.Data(),texPrefix.Data());
   DumpOFile(name,namet);
   
}
 


#include "VubAnalysis/fitInit.icc"
#include "VubAnalysis/fitUtil.icc"
