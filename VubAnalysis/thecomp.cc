#include <fstream.h>

#include "thecomp.hh"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "../RecoilAnalysis/mesData.hh" 
#include <TVector2.h>
thecomp::thecomp(char *var, char *dir, double mi, double ma, int b)
{
  char name[100];
  sprintf(name, "%s%s%s%s",dir,"/",var,".root");
  fHistFile = new TFile(name, "RECREATE");
   themin = mi;
   themax = ma;
   thebins = b;
   thevar = var ;
   TFile dal("dalitz.root");
   mydalitz = (*(TH2D*)gDirectory->Get("truedalitzsigforcut"));
   //    Dvar = new recoilDSys("ddecay.table",0,0);
   //    Bsem = new recoilDSys(0);
   
}
 
void thecomp::Loop(int nevents, int cat)
{

  fHistFile->cd();
  
  char name[100];
  char le[100];
  int sel=0;  
  //id = Mx category == Mx bin
  int id(0);
  int idflav;
  int group;        
  
  if (fChain == 0) return;
  
  Int_t nentries = Int_t(fChain->GetEntries());
  
  if( nentries > nevents) nentries = nevents;
  
  cout <<  nentries << endl;

  Int_t nbytes = 0, nb = 0;
  for (Int_t jentry=0; jentry<nentries;jentry++) {
    LoadTree(jentry); //in case of a TChain, ientry is the entry number in the current file
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%10000 == 0) cout << "Event " << jentry << endl;

     // set of analisys cuts...
    
    // B type
    bool rCh = (TMath::Abs(brecocharge) == BTYPE || BTYPE == 2);
    
    // acceptance cuts for leptons
    bool isLp = ((nle > 0) && (tlab<2.37) && (tlab>0.36) && (plab>0.5) && rCh);  
    if(LEPTTYPE == 0)
      isLp = ((nel > 0) && (tlab<2.37) && (tlab>0.36) && (plab>0.5) && rCh);
    if(LEPTTYPE == 1) 
      isLp = ((nmu > 0) && (tlab<2.37) && (tlab>0.36) && (plab>0.5) && rCh);
    
    // charge correlation
    int flav =  lcharge + brecoflav;
    
    
    // cut on the number of leptons in the recoil
    bool onelep = nle <= NLEP;
    if(LEPTTYPE == 0)
      onelep = nel <= NLEP;
    if(LEPTTYPE == 1) 
      onelep = nmu <= NLEP;

    
    // total charge
    int ch = TMath::Abs(xcharge + brecocharge);  
    
    // wdeltam cut
    bool WdeltaCut = (wdeltam>PRMM2CUT && brecocharge == 0); //last cut based on wdelta
    
    // mm2 (per exclusive mode)
    double themm2=mm2;
    if(FITCATEGORY == -11 && nneu !=1) themm2 = mm2pi0;
    if(FITCATEGORY == -11 && nneu ==1) themm2 = mm2gamma;
    if(FITCATEGORY == 11) themm2 = mm2pi;
    if(FITCATEGORY == -12) themm2 = mm2eta;
    if(FITCATEGORY == -13) themm2 = mm2rho0;
    if(FITCATEGORY == 13) themm2 = mm2rho;
    if(FITCATEGORY == -14) themm2 = mm2omega;
    if(FITCATEGORY == -15) themm2 = mm2etap;
    if(FITCATEGORY == -19) themm2 = mm2a0;
    if(FITCATEGORY == 19) themm2 = mm2a0p;
    bool Mm2 = (themm2 < MNUSQHIGH && themm2 > MNUSQLOW);
    if(FITCATEGORY == -11 && nneu ==1)  MM2 = (themm2 < 2.);
    
    // number of charged tracks
    bool NCHG = nchg-1 >= NCHGLOWEXCL && nchg-1 <= NCHGHIGHEXCL;
    
    // momentum of the daughters
    bool MOM1(1), MOM2(1); 
    if(FITCATEGORY == -11){
      MOM1 = mom1phpi0>MOM1MIN;
      MOM2 = mom2phpi0>MOM2MIN;        
    }
    if(FITCATEGORY == -13){
      MOM1 = mom1pirho0>MOM1MIN;
      MOM2 = mom2pirho0>MOM2MIN;        
    }
    if(FITCATEGORY == 13){
      MOM1 = mompirho>MOM1MIN;
      MOM2 = mompi0rho>MOM2MIN;        
    }
    if(FITCATEGORY == -14){
      MOM1 = mom1piome>MOM1MIN;
      MOM2 = mom2piome>MOM2MIN;        
    }

    // kaon rejection
    bool isK = 0;
    if(KSELE) isK = nkp + nks; // to be improved
    
    // mass cut
    bool Mass;
    double mass = -999;
    if(FITCATEGORY == -11){
      mass = mpi0;
     }
    if(FITCATEGORY == 11){
      mass = 0.13957;
    }
    if(FITCATEGORY == -12){
      mass = meta;
    }
    if(FITCATEGORY == -13){
      mass = mrho0;
    }
    if(FITCATEGORY == 13){
      mass = mrho;
    }
    if(FITCATEGORY == -14){
      mass = momega;
    }
    if(FITCATEGORY == -15){
      mass = metap;
    }
    if(FITCATEGORY == -19){
      mass = ma0;
    }
    if(FITCATEGORY == 19){
      mass = ma0p;
    }
    Mass = (mass < MXCUTHIGHEXCL && mass > MXCUTLOWEXCL);
    
    // number of combinations
    bool Comb = 1;
    if(FITCATEGORY == -13){
      Comb =  nrecorho0 >= NCOMBLOWEXCL && nrecorho0 <= NCOMBHIGHEXCL;
    }
    if(FITCATEGORY == 13){
      Comb =  nrecorho >= NCOMBLOWEXCL && nrecorho <= NCOMBHIGHEXCL;
    }
    if(FITCATEGORY == -14){
      Comb =  nrecoomega >= NCOMBLOWEXCL && nrecoomega <= NCOMBHIGHEXCL;
    }
    
    // number of pi0s
    bool Piz;
    Piz =  nrecopi0 >= NPI0LOWEXCL && nrecopi0 <= NPI0HIGHEXCL;
    
    // Dalitz plot
    bool DALITZ = 1;
    if(FITCATEGORY == -14){
      double xbindal = int(dalitzpi1pi0ome/0.016);
      double ybindal = int(dalitzpi1pi2ome/0.016);
      if(xbindal<0) xbindal = 1;
      if(ybindal<0) ybindal = 1;      
      if(mydalitz.GetBinContent(xbindal,ybindal)<DALITZCUT) DALITZ = 0;
    }     

    // mm2 cut in order to remove Vub crossfeed
    bool MM2CROSS = 1;
    if(FITCATEGORY == 11){
      bool MM2RHO = (mm2rho > MNUSQRHOHIGH || mm2rho < MNUSQRHOLOW);
      MM2CROSS = MM2RHO;
    }  
    if(FITCATEGORY == -12){
      bool MM2PI0 = (mm2pi0 > MNUSQPI0HIGH || mm2pi0 < MNUSQPI0LOW);
      MM2CROSS = MM2PI0;
    }  
    if(FITCATEGORY == -15){
      bool MM2RHO0 = (mm2rho0 > MNUSQRHO0HIGH || mm2rho0 < MNUSQRHO0LOW);
      bool MM2OMEGA = (mm2omega > MNUSQOMEGAHIGH || mm2omega < MNUSQOMEGALOW);
      MM2CROSS = MM2RHO0 * MM2OMEGA;
    }    
    if(FITCATEGORY == -19){
      bool MM2PI0 = (mm2pi0 > MNUSQPI0HIGH || mm2pi0 < MNUSQPI0LOW);
      bool MM2ETA = (mm2eta > MNUSQETAHIGH || mm2eta < MNUSQETALOW);
      bool MM2OMEGA = (mm2omega > MNUSQOMEGAHIGH || mm2omega < MNUSQOMEGALOW);      
      MM2CROSS = MM2PI0 * MM2ETA * MM2OMEGA;
    }
    if(FITCATEGORY == 19){
      bool MM2RHO = (mm2rho > MNUSQRHOHIGH || mm2rho < MNUSQRHOLOW);
      MM2CROSS = MM2RHO;
    }
    // additional cuts for pi l nu
    bool exclcuts=1;
    if(FITCATEGORY == 11){
      exclcuts=(  eneu<MAXENEU )&& (TMath::Abs(mtrkpi-3.1)>JPSIWIN);
      exclcuts=exclcuts&&((lepmap&3)!=3); 
    }
    bool PCMS=(pcms>LEPTONPCUT);
    bool CH=(ch < CHHIGH && ch > CHLOW);
    double myvar = -1;
    if (strcmp(thevar,"mes") == 0) {myvar = mes; sprintf(xname,"m_{es} [GeV]");}
    if (strcmp(thevar,"mm2") == 0) {myvar = mm2; MM2=1;  sprintf(xname,"m_{miss}^{2} [GeV^{2}]");}
    if (strcmp(thevar,"prmm2") == 0) {myvar = wdeltam; WdeltaCut=0;  sprintf(xname,"m_{missPR}^{2} [GeV^{2}]");}
    if (strcmp(thevar,"nchg") == 0) { myvar = nchg-1; NCHG = 1; CH=0;  sprintf(xname,"N_{chg}");}
    if (strcmp(thevar,"nneu") == 0) { myvar = nneu; sprintf(xname,"N_{neu}");}
    if (strcmp(thevar,"nrecopi0") == 0) { myvar = nrecopi0; Piz = 1; sprintf(xname,"N_{#pi0}");}
    if (strcmp(thevar,"pcms") == 0) {myvar = pcms;PCMS=1;  sprintf(xname,"P_{cms} [GeV]");}
    if (strcmp(thevar,"qtot") == 0) {myvar = ch; CH = 0;NCHG = 1; sprintf(xname,"#Sigma_{i} chg_{i}");}
    if (strcmp(thevar,"mm2pi") == 0) {myvar = mm2pi; Mm2 =1;sprintf(xname,"mm2_{#pi} [GeV^{2}]");}
    if (strcmp(thevar,"mpi0") == 0) {myvar = mpi0; Mass =1;sprintf(xname,"m_{#pi0} [GeV]");}
    if (strcmp(thevar,"resompi0") == 0) {myvar = mpi0-mxhadgen; Mass =1;sprintf(xname,"m_{#pi0}-m^{true}_{#pi0} [GeV]");}
    if (strcmp(thevar,"mm2pi0") == 0) {myvar = mm2pi0; Mm2 =1; MM2CROSS=1; sprintf(xname,"mm2_{#pi0} [GeV^{2}]");}
    if (strcmp(thevar,"mm2gamma") == 0) {myvar = mm2gamma; Mm2 =1;sprintf(xname,"mm2_{#gamma} [GeV^{2}]");}
    if (strcmp(thevar,"mom1phpi0") == 0) {myvar = mom1phpi0; MOM1 =1;sprintf(xname,"mom1#gamma_{#pi0} [GeV]");}
    if (strcmp(thevar,"mom2phpi0") == 0) {myvar = mom2phpi0; MOM2 =1;sprintf(xname,"mom2#gamma_{#pi0} [GeV]");}
    if (strcmp(thevar,"truemom1phpi0") == 0) {myvar = truemom1phpi0; MOM1 =1;sprintf(xname,"truemom1#gamma_{#pi0} [GeV]");}
    if (strcmp(thevar,"truemom2phpi0") == 0) {myvar = truemom2phpi0; MOM2 =1;sprintf(xname,"truemom2#gamma_{#pi0} [GeV]");}
    if (strcmp(thevar,"nrecorho0") == 0) { myvar = nrecorho0; Comb = 1; sprintf(xname,"N_{#rho0}");}
    if (strcmp(thevar,"mrho0") == 0) {myvar = mrho0; Mass =1;sprintf(xname,"m_{#rho0} [GeV]");}
    if (strcmp(thevar,"resomrho0") == 0) {myvar = mrho0-mxhadgen; Mass =1;sprintf(xname,"m_{#rho0}-m^{true}_{#rho0} [GeV]");}
    if (strcmp(thevar,"mm2rho0") == 0) {myvar = mm2rho0; Mm2 =1; MM2CROSS=1; sprintf(xname,"mm2_{#rho0} [GeV^{2}]");}
    if (strcmp(thevar,"mom1pirho0") == 0) {myvar = mom1pirho0; MOM1 =1;sprintf(xname,"mom1#pi_{#rho0} [GeV]");}
    if (strcmp(thevar,"mom2pirho0") == 0) {myvar = mom2pirho0; MOM2 =1;sprintf(xname,"mom2#pi_{#rho0} [GeV]");}
    if (strcmp(thevar,"truemom1pirho0") == 0) {myvar = truemom1pirho0; MOM1 =1;sprintf(xname,"truemom1#pi_{#rho0} [GeV]");}
    if (strcmp(thevar,"truemom2pirho0") == 0) {myvar = truemom2pirho0; MOM2 =1;sprintf(xname,"truemom2#pi_{#rho0} [GeV]");}
    if (strcmp(thevar,"nrecorho") == 0) { myvar = nrecorho; Comb = 1; sprintf(xname,"N_{#rho}");}
    if (strcmp(thevar,"mrho") == 0) {myvar = mrho; Mass =1;sprintf(xname,"m_{#rho} [GeV]");}
    if (strcmp(thevar,"resomrho") == 0) {myvar = mrho-mxhadgen; Mass =1;sprintf(xname,"m_{#rho}-m^{true}_{#rho} [GeV]");}
    if (strcmp(thevar,"mm2rho") == 0) {myvar = mm2rho; Mm2 =1; MM2CROSS=1; sprintf(xname,"mm2_{#rho} [GeV^{2}]");}
    if (strcmp(thevar,"mompirho") == 0) {myvar = mompirho; MOM1 =1;sprintf(xname,"mom#pi_{#rho} [GeV]");}
    if (strcmp(thevar,"mompi0rho") == 0) {myvar = mompi0rho; MOM2 =1;sprintf(xname,"mom#pi0_{#rho} [GeV]");}
    if (strcmp(thevar,"truemompirho") == 0) {myvar = truemompirho; MOM1 =1;sprintf(xname,"truemom#pi_{#rho} [GeV]");}
    if (strcmp(thevar,"truemompi0rho") == 0) {myvar = truemompi0rho; MOM2 =1;sprintf(xname,"truemom#pi0_{#rho} [GeV]");}
    if (strcmp(thevar,"nrecoomega") == 0) { myvar = nrecoomega; Comb = 1; sprintf(xname,"N_{#omega}");}
    if (strcmp(thevar,"momega") == 0) {myvar = momega; Mass =1;sprintf(xname,"m_{#omega} [GeV]");}
    if (strcmp(thevar,"resomomega") == 0) {myvar = momega-mxhadgen; Mass =1;sprintf(xname,"m_{#omega}-m^{true}_{#omega} [GeV]");}
    if (strcmp(thevar,"mm2omega") == 0) {myvar = mm2omega; Mm2 =1; MM2CROSS=1; sprintf(xname,"mm2_{#omega} [GeV^{2}]");}
    if (strcmp(thevar,"mom1piome") == 0) {myvar = mom1piome; MOM1 =1;sprintf(xname,"mom1#pi_{#omega} [GeV]");}
    if (strcmp(thevar,"mom2piome") == 0) {myvar = mom2piome; MOM2 =1;sprintf(xname,"mom2#pi_{#omega} [GeV]");}
    if (strcmp(thevar,"truemom1piome") == 0) {myvar = truemom1piome; MOM1 =1;sprintf(xname,"truemom1#pi_{#omega} [GeV]");}
    if (strcmp(thevar,"truemom2piome") == 0) {myvar = truemom2piome; MOM2 =1;sprintf(xname,"truemom2#pi_{#omega} [GeV]");}
    if (strcmp(thevar,"costhome") == 0) {myvar = costhome;sprintf(xname,"cos(#theta) ");}
    if (strcmp(thevar,"meta") == 0) {myvar = meta; Mass =1;sprintf(xname,"m_{#eta} [GeV]");}
    if (strcmp(thevar,"mm2eta") == 0) {myvar = mm2eta; Mm2 =1; MM2CROSS=1; sprintf(xname,"mm2_{#eta} [GeV^{2}]");}
    if (strcmp(thevar,"modeeta") == 0) {myvar = modeeta;sprintf(xname,"mode_{#eta}");}
    if (strcmp(thevar,"resometa") == 0) {myvar = meta-mxhadgen; Mass =1;sprintf(xname,"m_{#eta}-m^{true}_{#eta} [GeV]");}
    if (strcmp(thevar,"metap") == 0) {myvar = metap; Mass =1;sprintf(xname,"m_{#eta'} [GeV]");}
    if (strcmp(thevar,"mm2etap") == 0) {myvar = mm2etap; Mm2 =1;sprintf(xname,"mm2_{#eta'} [GeV^{2}]");}
    if (strcmp(thevar,"modeetap") == 0) {myvar = modeetap;sprintf(xname,"mode_{#eta'}");}
    if (strcmp(thevar,"resometap") == 0) {myvar = metap-mxhadgen; Mass =1;sprintf(xname,"m_{#eta'}-m^{true}_{#eta'} [GeV]");}
    if (strcmp(thevar,"momrho0ph") == 0) {myvar = momrho0ph;MOM1=MOM2 = 1; sprintf(xname,"mom#gamma_{#eta'}[GeV]");}
    if (strcmp(thevar,"etapmassetadau") == 0) {myvar = etapmassetadau; sprintf(xname,"m_{#eta}#eta' dau[GeV]");}
    if (strcmp(thevar,"ma0p") == 0) {myvar = ma0p; Mass =1;sprintf(xname,"m_{a_{0}^{#pm}}} [GeV]");}
    if (strcmp(thevar,"mm2a0p") == 0) {myvar = mm2a0p; Mm2 =1;sprintf(xname,"mm2_{a_{0}^{#pm}} [GeV^{2}]");}
    if (strcmp(thevar,"resoma0p") == 0) {myvar = ma0p-mxhadgen; Mass =1;sprintf(xname,"m_{a_{0}^{#pm}-m^{true}_{a_{0}^{#pm'}}  [GeV]");}
    if (strcmp(thevar,"a0massetadau") == 0) {myvar = a0massetadau;sprintf(xname,"m_{#eta}a_{0} dau[GeV]");}
    if (strcmp(thevar,"modea0") == 0) {myvar = modea0;sprintf(xname,"mode_{a_{0}^{0}}");}
    if (strcmp(thevar,"ma0") == 0) {myvar = ma0; Mass =1;sprintf(xname,"m_{a_{0}'} [GeV]");}
    if (strcmp(thevar,"mm2a0") == 0) {myvar = mm2a0; Mm2 =1;sprintf(xname,"mm2_{a_{0}'} [GeV^{2}]");}
    if (strcmp(thevar,"resoma0") == 0) {myvar = ma0-mxhadgen; Mass =1;sprintf(xname,"m_{a_{0}'}-m^{true}_{a_{0}'} [GeV]");}
    if (strcmp(thevar,"a0pmassetadau") == 0) {myvar = a0pmassetadau;sprintf(xname,"m_{#eta}a_{0}^{+-} dau[GeV]");}
    if (strcmp(thevar,"modea0p") == 0) {myvar = modea0p;sprintf(xname,"mode_{a_{0}^{+-}}");}
    
    // lepton requests
    bool lPYes = (PCMS&& isLp && !(TMath::Abs(brecocharge)!=0 && (flav)!=0) );
    // ALL TOGETHER no mass
    bool AllCutnomass = (lPYes && onelep && q2>Q2CUT && Mm2 &&  CH && !(WdeltaCut) && 
                         NCHG && MOM1 && MOM2 && Comb && Piz && !isK && DALITZ && MM2CROSS && exclcuts);  

    // ALL TOGETHER
    bool AllCut = (lPYes && onelep && q2>Q2CUT && Mm2 && CH && !(WdeltaCut) && 
                   NCHG && MOM1 && MOM2 && Mass && Comb && Piz && !isK && DALITZ && MM2CROSS && exclcuts);  
   
      
      // assign flavor category
//      int cutflav = 5;
//      if(TMath::Abs(brecocharge)) cutflav = 3;
//      if(TMath::Abs(brecocharge)==0 && flav==0) cutflav = 4;	  
      int cutflav=0;
      if( AllCut){  	
      
	sel++;
      id = hist(myvar);
      idflav = hist(myvar) + cutflav * 10000;    
      group = 200 + cat * 1000 + 100000;
      sprintf(name, "%s%d",le,group+idflav);
    //  if((!KSELE) && MM2 && CH && NLE1 && NCHG && !WdeltaCut && DELTA && MX && MOM1 && MOM2 && COMB && PIZ && DALITZ && MM2PI0 && MM2ETA && MM2RHO && MM2RHO0 && MM2OMEGA && DAUETA && DAUGAMMOM&&exclcuts){  	
        sprintf(le, "h");
	sprintf(name, "%s%d",le,group+idflav);
	if (id!=-1) ((TH1D*)gDirectory->Get(name))->Fill(mes); 
      }
    
  //  if(PCMS && NLE && !flav && BCH && IPUR){  
 //     sprintf(le, "h");
  //    group = 90206 + CATEGORY * 1000 + 100000;
   //   sprintf(name, "%s%d",le,group);
  //    if (id!=-1) ((TH1D*)gDirectory->Get(name))->Fill(mes);  
 //     if((!KSELE) && MM2 && CH && NLE1 && NCHG && !WdeltaCut && DELTA && MX && MOM1 && MOM2 && COMB && PIZ && DALITZ && MM2PI0 && MM2ETA && MM2RHO && MM2RHO0 && MM2OMEGA && DAUETA && DAUGAMMOM){  
//	group = 90206 + CATEGORY * 1000;
//	sprintf(name, "%s%d",le,group);
//	if (id!=-1) ((TH1D*)gDirectory->Get(name))->Fill(mes);  
//      }
//    }
    
  }
   

cout<<"selected "<<sel<<endl;
}
// ----------------------------------------------------------------------
double thecomp::getBsysweight(int decType,int thevub) {
  double theweight;
  theweight = 1.;  
  theweight = Bsem->weight(decType); 
  if(thevub) theweight = 1.;
  return theweight;
}
// ----------------------------------------------------------------------
double thecomp::getDsysweight(int decDpi,int decDk,int decDks,int decDpiz,int decDlep,int decImode,int thevub) {
  double theweight;
  theweight = 1.;  
  int bla(1);
  theweight = Dvar->weight(decDpi,decDk,decDks,decDpiz,decDlep,decImode,bla);    
  if(thevub) theweight = 1.;
  return theweight;
}

double thecomp::FermiWeight(double kp, double deltamb, double deltaa){

  Double_t BMASS, BQMASS, A0;
  BMASS   = 5.2792;
  BQMASS  = 4.800;
  A0      = 1.29;
  double mb1 = BQMASS + deltamb;
  
  double a1 = A0 + deltaa ;
  double w81;

  Int_t nIntegral = 10000;
  const double kmin = -5.;
  const double kmax =  1.;
  double Nold = 0;
  double Nnew1(0);

  // -- Compute reweighting ratios
  for (int i = 0; i < nIntegral; ++i) {
    double kplus = (i+0.5)/((double)nIntegral)*(kmax-kmin)+kmin;
    Nold += (kmax-kmin)/((double)nIntegral)*fermi(kplus, BQMASS, A0);
    Nnew1 += (kmax-kmin)/((double)nIntegral)*fermi(kplus, mb1, a1);
  }
  double kpold = kp; 
  double kpnew1 = 5.279 - mb1 - (5.279 - BQMASS - kp);
  
  double fold = fermi(kpold, BQMASS, A0)/Nold; 
  double fnew1 = fermi(kpnew1, mb1, a1)/Nnew1;
  if(fold<=0)cout<<"?? "<<fold<<" "<<kpold<<" "<<BQMASS<<" "<<A0<<endl;
  w81 = fnew1/fold; 
  return w81;
}

  // -- Shape function
double thecomp::fermi(double kp, double m, double a) {
  double BMASS   = 5.2792;
  double x = kp/(BMASS - m);
  if ((kp>-m) && (x <= 1.)) {
    return TMath::Power((1-x), a) * TMath::Exp((1+a)*x); 
  } 
  return 0.;
}


void thecomp::Bookhist()
{
  fHistFile->cd();
  TH1 *h;
  char name[100], title[100], number[100];
  int lo; 
 
  sprintf(name, "h1100000");  sprintf(title, "%s%s", thevar, " data");  h = new TH1D(name, title, thebins, themin, themax  );  h->Sumw2();
  sprintf(name, "h1200000");  sprintf(title, "%s%s", thevar, " MC");  h = new TH1D(name, title, thebins, themin, themax );  h->Sumw2();

  for (lo=1;lo<thebins+1;lo++) {
    
    sprintf(number, "%d",  200+lo);     

    sprintf(name,"%s%s" , "h101",number);  sprintf(title, "mes data");  h = new TH1D(name, title, 40, 5.2, 5.3);         h->Sumw2();
    sprintf(name,"%s%s" , "h131",number);  sprintf(title, "mes data bch");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h141",number);  sprintf(title, "mes data b0 os");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h151",number);  sprintf(title, "mes data b0 ss");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h102",number);  sprintf(title, "mes MC");  h = new TH1D(name, title, 40, 5.2, 5.3);         h->Sumw2();
    sprintf(name,"%s%s" , "h132",number);  sprintf(title, "mes MC bch");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h142",number);  sprintf(title, "mes MC b0 os");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h152",number);  sprintf(title, "mes MC b0 ss");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();

  }
   sprintf(name, "h8888");  sprintf(title, "test gauss");  h = new TH1D(name, title, 100, -.05, .05);    h->Sumw2(); 
}

int thecomp::hist(double mx){

  // categories
  int bin = int((mx-themin)/((themax-themin)/thebins)+1);
  if (mx>themax || mx==themax) bin = -1;
  if (mx<themin || mx==themin) bin = -1;
  return bin;
}

void thecomp::Fitmes(int cat, int cut){

  fHistFile->cd();
  // fit to mes distribution and fill of the Mx plots...

  //  gROOT->cd();

  int group = 200 + cat * 1000;
   
  char name[100];
  char namebch[100];
  char nameb0os[100];
  char nameb0ss[100];
  char le[100];
  double sigs[1000];
  double errsigs[1000];
  double tempbin;
  double tempbinchb;
  double tempbinb0os;
  double tempbinb0ss;
  double temperr;
  double temperrchb;
  double temperrb0os;
  double temperrb0ss;
  double chid = 0.174;
  int is;
  int thegroup = group;      
  double foomean;
  double foosigma;
  double fooalpha;
  double foon;
  double usedmean;
  double usedsigma;
  double usedalpha;
  double usedn;
  int i;
  int addcut = 1; 
//  if(cut) addcut = 1;
  int const theb = thebins + 1;
  sprintf(le, "h");

  //extracting the signal fit parameters 
//  sprintf(name, "%s%d",le,group+90006+addcut*100000);  
//  cout << name << endl;
//  cout << ((TH1D*)gDirectory->Get(name))->Integral() << endl;
//  double dummy1,dummy2;
//  sighisto(dummy1,dummy2,(TH1D*)gDirectory->Get(name),usedmean,usedsigma,usedalpha,usedn,1,-11111111.,-1111111.,-1111111.,-1111111.,-1111111.);
//  cout << "mes result for: " << name << " is, MEAN " << usedmean << " SIGMA " << usedsigma << " ALPHA " << usedalpha << " N " << usedn << endl;
  
  group = thegroup;
  for (i = 1;  i < theb; i++) { 
    is = i;
    sprintf(name, "%s%d",le,group+is+addcut*100000);  
    
//    sprintf(namebch, "%s%d",le,group+is+30000+addcut*100000);  
//    sprintf(nameb0os, "%s%d",le,group+is+40000+addcut*100000);  
//    sprintf(nameb0ss, "%s%d",le,group+is+50000+addcut*100000);  
    
    // mixing correction      
    
//    for(int k=1;k<40;k++){
//      tempbinchb = ((TH1D*)gDirectory->Get(namebch))->GetBinContent(k);
//      tempbinb0os = ((TH1D*)gDirectory->Get(nameb0os))->GetBinContent(k);
//      tempbinb0ss = ((TH1D*)gDirectory->Get(nameb0ss))->GetBinContent(k);
//      temperrchb = ((TH1D*)gDirectory->Get(namebch))->GetBinError(k);
//      temperrb0os = ((TH1D*)gDirectory->Get(nameb0os))->GetBinError(k);
//      temperrb0ss = ((TH1D*)gDirectory->Get(nameb0ss))->GetBinError(k);
//      tempbin = tempbinchb + ((1-chid)/(1-2*chid)) * tempbinb0os - (chid/(1-2*chid)) * tempbinb0ss;
//      temperr = sqrt(temperrchb*temperrchb + ((1-chid)/(1-2*chid)) * ((1-chid)/(1-2*chid)) * tempbinb0os + (chid/(1-2*chid)) * (chid/(1-2*chid))* tempbinb0ss);
//      ((TH1D*)gDirectory->Get(name))->SetBinContent(k, tempbin);
//      ((TH1D*)gDirectory->Get(name))->SetBinError(k, temperr); 
//    }	
    
    // put the result of the fit in each Mx bin
    sighisto(sigs[is-1],errsigs[is-1],(TH1D*)gDirectory->Get(name),foomean,foosigma,fooalpha,foon,0,usedmean,usedsigma,usedalpha,usedn,-1111111);
    if(!strcmp(thevar,"mes")){      
      sigs[is-1] =((TH1D*)gDirectory->Get(name))->Integral();
      errsigs[is-1] = sqrt(((TH1D*)gDirectory->Get(name))->Integral());     
    }
    int title = cat * 100000 + addcut * 1000000;
    sprintf(name, "%s%d",le,title);  
    ((TH1D*)gDirectory->Get(name))->SetBinContent(is, sigs[is-1]);
    ((TH1D*)gDirectory->Get(name))->SetBinError(is, errsigs[is-1]);
  }
}

void 
thecomp::sighisto(double&signal, double& signalErr, TH1D *histo, double &resmean, double &ressigma, double &resalpha, double &resn, int fixpar, double mean, double sigma, double alpha,  double n, double argus){

  recoilAnalysis b;
  mesData themes;
  double thesigma = -1111111;
  double themean = -1111111;
  double theargus = -1111111;
  double thealpha = -1111111;
  double then = -1111111;
  if (fixpar){
    themean = mean;
    thesigma = sigma;
    thealpha = alpha;
    then = n;
    mesData tempthemes(*(b.vubMes(histo, resmean, ressigma, resalpha, resn, 1, 1, themean, thesigma, thealpha, then, theargus)));
    themes = tempthemes;
  }else{
    mesData tempthemes(*(b.vubMes(histo, resmean, ressigma, resalpha, resn, 1,1, themean, thesigma, thealpha, then, theargus)));
    themes = tempthemes;
  }
  signal=themes.theSig(); signalErr=themes.theErrSig(); 

}

void thecomp::overlap( int norm, TString dir){
   
   char mycut[100];
   int addcut = 1; 

   c6 = new TCanvas("c6", "c6", 300.,   0., 400,400);  
   // -- top left
   fPads[1]= new TPad("pad1", "", 0.01, 0.25, 0.99, 0.99);   fPads[1]->Draw(); 
   fPads[2]= new TPad("pad2", "", 0.01, 0.01, 0.99, 0.24);   fPads[2]->Draw();  
   // -- bottom left
   char name[100], name2[100], preine[100], preich[100], prefix[100], shift[100], smear[100], line[100];
   int h = 0 ; 

   sprintf(name,"%s%d" ,"h", 100000 +  1000000);
   TH1D y1;
   ((TH1D*)gDirectory->Get(name))->Copy(y1);
   sprintf(name,"%s%d" ,"h", 200000 +  1000000);
   TH1D y2;
   ((TH1D*)gDirectory->Get(name))->Copy(y2);

   fPads[1]->cd();
   y1.SetMarkerStyle(20);
   double inter1;
   double max = 1.4*y1.GetMaximum();
   y1.SetLabelSize(0.07, "Y");
   y1.SetMaximum(max);
   y1.SetMinimum(0.);
   inter1 = y1.Integral();
   cout<<"Integral 1 "<<inter1<<endl;
   y1.SetNormFactor(inter1);
   y1.SetLineColor(kRed);
   y1.Draw();
   cout<<y2.GetBinContent(1)<<" xxx "<<y1.GetBinContent(2)<<endl;
   cout<<((TH1D*)gDirectory->Get("h1200000"))->GetBinContent(1)<<endl;
   
 
   double inter2;
   inter2 = y2.Integral();
   cout<<"Integral 2 "<<inter2<<endl;  
   intdataen = inter1;
   intMCen = inter2;
   y2.SetFillColor(2);
   y2.SetFillStyle(3004);
   if (norm) {
      if(intMCen>0)y2.Scale(intdataen/intMCen);
   }else{
      if(inter2>0)y2.Scale(inter1/inter2);
   }   
   y2.Draw("samehist");
   
   double g = chisq(&y1,&y2);
   TLatex *tl=new TLatex;
   sprintf(line, "#chi^{2} = %5.4f", g); tl->SetTextSizePixels(25); tl->DrawLatex(0.22, 0.75, line);
  
   fPads[2]->cd();
   //shrinkPad(0.001, 0.2); 
   //shrinkPad(0.4, 0.2, 0.1, 0.001);

   gPad->SetGridx(1);  gPad->SetGridy(1);
   TH1D* hratio=new TH1D("hratio", "hratio", thebins, themin, themax  ); 
   //TH1D *hratio = new TH1D(y1); hratio->SetName("hratio"); 
   //hratio->Reset();
   hratio->Divide( &y1,&y2);
   hratio->SetMinimum(0.5); hratio->SetMaximum(1.5);
   hratio->SetMarkerStyle(24);
   hratio->SetNdivisions(504, "Y");
   hratio->SetLabelSize(0.22, "X");  hratio->SetLabelSize(0.17, "Y");
   hratio->SetStats(0);
   hratio->SetTitle();
   hratio->Draw();

//    cout << "the chi square is " << g << endl;
//    SHIFTNEUT = SHIFTNEUT *1000;
//    SIGMANEUT = SIGMANEUT *1000;
//    sprintf(shift, "%s%d","-", SHIFTNEUT);
//    sprintf(smear, "%s%d","-", SIGMANEUT);   
  
   sprintf(mycut, "allcuts");   
   if (norm) {   
      sprintf(name, "%s%s%s%s", dir.Data(), "/comparisonnorm",thevar,".eps");
   }else{
      sprintf(name, "%s%s%s%s", dir.Data(), "/comparison",thevar,".eps");
   }
   
   c6->SaveAs(name);  
//   sprintf(name, "%s%s%s%s%s%s%s", "comparison",thevar,"smear",smear,"shift",shift,".dat");
}

double thecomp::chisq(TH1 *h1, TH1 *h2){
  double chisq = 0; 
  double tempchisq = 0;
  double inter;
  double temperr1;
  double temperr2;
  inter = h1->Integral();
  double inter2;
  inter2 = h2->Integral();
  int nbins = h1->GetNbinsX();
  for(int i=1;i<nbins;i++) {
    tempchisq = 0;
    temperr1 = 1;
    temperr2 = 1;
    if (h1->GetBinError(i)>.9) temperr1 = h1->GetBinError(i);
    if (h2->GetBinError(i)>.9) temperr2 = h2->GetBinError(i);
    if (((h1->GetBinError(i))*(h1->GetBinError(i))+(h2->GetBinError(i))*(h2->GetBinError(i)))>0) tempchisq = ((h1->GetBinContent(i)) - (h2->GetBinContent(i))*inter/inter2) / sqrt((temperr1*temperr1)+(temperr2*temperr2)*(inter/inter2)*(inter/inter2));  
    chisq = chisq + tempchisq*tempchisq;
  }
  chisq = chisq / nbins;
  return chisq;
}


// ----------------------------------------------------------------------
Double_t thecomp::smeargauss(double invalue, double mean, double sigma){    

  // smears using a gaussian distribution

  double returnvalue=0.;
  double xsmear = gRandom->Gaus(SHIFTNEUT, SIGMANEUT); 
  returnvalue = invalue + xsmear;
  //cout << sigma << mean << " " << invalue << " " << returnvalue-mean-invalue << endl;
  ((TH1D*)gDirectory->Get("h8888"))->Fill(returnvalue - invalue);
  return returnvalue;
}

void thecomp::Init(TTree *tree)
{
//   Set branch addresses
   if (tree == 0) return;
   fChain    = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("mes",&mes);
   fChain->SetBranchAddress("de",&de);
   fChain->SetBranchAddress("pur",&pur);
   fChain->SetBranchAddress("intpur",&intpur);
   fChain->SetBranchAddress("mode",&mode);
   fChain->SetBranchAddress("brecoflav",&brecoflav);
   fChain->SetBranchAddress("brecocharge",&brecocharge);
   fChain->SetBranchAddress("mxhad",&mxhad);
   fChain->SetBranchAddress("mxhadgen",&mxhadgen);
   fChain->SetBranchAddress("q2",&q2);
   fChain->SetBranchAddress("q2fit",&q2fit);
   fChain->SetBranchAddress("q2gen",&q2gen);
   fChain->SetBranchAddress("pcmsgen",&pcmsgen);
   fChain->SetBranchAddress("vub",&vub);
   fChain->SetBranchAddress("vcb",&vcb);
   fChain->SetBranchAddress("other",&other);
   fChain->SetBranchAddress("xcharge",&xcharge);
   fChain->SetBranchAddress("mxhad",&mxhad);
   fChain->SetBranchAddress("mxhadfit",&mxhadfit);
   fChain->SetBranchAddress("mxhadchg",&mxhadchg);
   fChain->SetBranchAddress("lcharge",&lcharge);
   fChain->SetBranchAddress("plab",&plab);
   fChain->SetBranchAddress("tlab",&tlab);
   fChain->SetBranchAddress("flab",&flab);
   fChain->SetBranchAddress("pcms",&pcms);
   fChain->SetBranchAddress("tcms",&tcms);
   fChain->SetBranchAddress("fcms",&fcms);
   fChain->SetBranchAddress("ecms",&ecms);
   fChain->SetBranchAddress("nle",&nle);
   fChain->SetBranchAddress("nel",&nel);
   fChain->SetBranchAddress("nmu",&nmu);
   fChain->SetBranchAddress("nchg",&nchg);
   fChain->SetBranchAddress("nneu",&nneu);
   fChain->SetBranchAddress("nkp",&nkp);
   fChain->SetBranchAddress("nks",&nks);
   fChain->SetBranchAddress("mm2",&mm2);
   fChain->SetBranchAddress("wdeltam",&wdeltam);
   fChain->SetBranchAddress("Gvxbtyp",&Gvxbtyp);
   fChain->SetBranchAddress("nrecopi0",&nrecopi0);
   fChain->SetBranchAddress("ENeu",&eneu);
   fChain->SetBranchAddress("mtrkpi",&mtrkpi);
   fChain->SetBranchAddress("pnupi",&pnupi);
   fChain->SetBranchAddress("lepmap",&lepmap);
   fChain->SetBranchAddress("mm2pi",&mm2pi);
   fChain->SetBranchAddress("mpi0",&mpi0);
   fChain->SetBranchAddress("mm2pi0",&mm2pi0);
   fChain->SetBranchAddress("mm2gamma",&mm2gamma);
   fChain->SetBranchAddress("truemom1phpi0",&truemom1phpi0);
   fChain->SetBranchAddress("truemom2phpi0",&truemom2phpi0);
   fChain->SetBranchAddress("truemomlab1phpi0",&truemomlab1phpi0);
   fChain->SetBranchAddress("truemomlab2phpi0",&truemomlab2phpi0);
   fChain->SetBranchAddress("trueth1phpi0",&trueth1phpi0);
   fChain->SetBranchAddress("trueth2phpi0",&trueth2phpi0);
   fChain->SetBranchAddress("mom1phpi0",&mom1phpi0);
   fChain->SetBranchAddress("mom2phpi0",&mom2phpi0);
   fChain->SetBranchAddress("truemrho",&truemrho);
   fChain->SetBranchAddress("nrecorho",&nrecorho);
   fChain->SetBranchAddress("mrho",&mrho);
   fChain->SetBranchAddress("mm2rho",&mm2rho);
   fChain->SetBranchAddress("truemompirho",&truemompirho);
   fChain->SetBranchAddress("truemompi0rho",&truemompi0rho);
   fChain->SetBranchAddress("mompirho",&mompirho);
   fChain->SetBranchAddress("mompi0rho",&mompi0rho);
   fChain->SetBranchAddress("truemrho0",&truemrho0);
   fChain->SetBranchAddress("nrecorho0",&nrecorho0);
   fChain->SetBranchAddress("mrho0",&mrho0);
   fChain->SetBranchAddress("mm2rho0",&mm2rho0);
   fChain->SetBranchAddress("truemom1pirho0",&truemom1pirho0);
   fChain->SetBranchAddress("truemom2pirho0",&truemom2pirho0);
   fChain->SetBranchAddress("mom1pirho0",&mom1pirho0);
   fChain->SetBranchAddress("mom2pirho0",&mom2pirho0);
   fChain->SetBranchAddress("truemomega",&truemomega);
   fChain->SetBranchAddress("nrecoomega",&nrecoomega);
   fChain->SetBranchAddress("momega",&momega);
   fChain->SetBranchAddress("mm2omega",&mm2omega);
   fChain->SetBranchAddress("truemom1piome",&truemom1piome);
   fChain->SetBranchAddress("truemom2piome",&truemom2piome);
   fChain->SetBranchAddress("truemompi0ome",&truemompi0ome);
   fChain->SetBranchAddress("truedalitzpi1pi2ome",&truedalitzpi1pi2ome);
   fChain->SetBranchAddress("truedalitzpi1pi0ome",&truedalitzpi1pi0ome);
   fChain->SetBranchAddress("truecosthome",&truecosthome);
   fChain->SetBranchAddress("mom1piome",&mom1piome);
   fChain->SetBranchAddress("mom2piome",&mom2piome);
   fChain->SetBranchAddress("mompi0ome",&mompi0ome);
   fChain->SetBranchAddress("dalitzpi1pi2ome",&dalitzpi1pi2ome);
   fChain->SetBranchAddress("dalitzpi1pi0ome",&dalitzpi1pi0ome);
   fChain->SetBranchAddress("costhome",&costhome);
   fChain->SetBranchAddress("meta",&meta);
   fChain->SetBranchAddress("mm2eta",&mm2eta);
   fChain->SetBranchAddress("modeeta",&modeeta);
   fChain->SetBranchAddress("metap",&metap);
   fChain->SetBranchAddress("mm2etap",&mm2etap);
   fChain->SetBranchAddress("modeetap",&modeetap);
   fChain->SetBranchAddress("etapmassetadau",&etapmassetadau);
   fChain->SetBranchAddress("momrho0ph",&momrho0ph);
   fChain->SetBranchAddress("ma0",&ma0);
   fChain->SetBranchAddress("mm2a0",&mm2a0);
   fChain->SetBranchAddress("a0massetadau",&a0massetadau);
   fChain->SetBranchAddress("modea0",&modea0);
   fChain->SetBranchAddress("ma0p",&ma0p);
   fChain->SetBranchAddress("mm2a0p",&mm2a0p);
   fChain->SetBranchAddress("a0pmassetadau",&a0pmassetadau);
   fChain->SetBranchAddress("modea0p",&modea0p);
   Notify();
}

Bool_t thecomp::Notify()
{
//   called when loading a new file
//   get branch pointers
   b_mes = fChain->GetBranch("mes");
   b_de = fChain->GetBranch("de");
   b_pur = fChain->GetBranch("pur");
   b_intpur = fChain->GetBranch("intpur");
   b_mode = fChain->GetBranch("mode");
   b_brecoflav = fChain->GetBranch("brecoflav");
   b_brecocharge = fChain->GetBranch("brecocharge");
   b_mxhadgen = fChain->GetBranch("mxhadgen");
   b_q2 = fChain->GetBranch("q2");
   b_q2fit = fChain->GetBranch("q2fit");
   b_q2gen = fChain->GetBranch("q2gen");
   b_pcmsgen = fChain->GetBranch("pcmsgen");
   b_vub = fChain->GetBranch("vub");
   b_vcb = fChain->GetBranch("vcb");
   b_other = fChain->GetBranch("other");
   b_xcharge = fChain->GetBranch("xcharge");
   b_mxhadfit = fChain->GetBranch("mxhadfit");
   b_mxhadchg = fChain->GetBranch("mxhadchg");
   b_lcharge = fChain->GetBranch("lcharge");
   b_plab = fChain->GetBranch("plab");
   b_tlab = fChain->GetBranch("tlab");
   b_flab = fChain->GetBranch("flab");
   b_pcms = fChain->GetBranch("pcms");
   b_tcms = fChain->GetBranch("tcms");
   b_fcms = fChain->GetBranch("fcms");
   b_ecms = fChain->GetBranch("ecms");
   b_nle = fChain->GetBranch("nle");
   b_nel = fChain->GetBranch("nel");
   b_nmu = fChain->GetBranch("nmu");
   b_nchg = fChain->GetBranch("nchg");
   b_nneu = fChain->GetBranch("nneu");
   b_nkp = fChain->GetBranch("nkp");
   b_nks = fChain->GetBranch("nks");
   b_mm2 = fChain->GetBranch("mm2");
   b_wdeltam = fChain->GetBranch("wdeltam");
   b_vxbtyp = fChain->GetBranch("Gvxbtyp");
   b_nrecopi0 = fChain->GetBranch("nrecopi0");
   b_mm2pi = fChain->GetBranch("mm2pi");
   b_eneu = fChain->GetBranch("ENeu");
   b_mtrkpi = fChain->GetBranch("mtrkpi");
   b_pnupi = fChain->GetBranch("pnupi");
   b_lepmap = fChain->GetBranch("lepmap");
   b_mpi0 = fChain->GetBranch("mpi0");
   b_mm2pi0 = fChain->GetBranch("mm2pi0");
   b_mm2gamma = fChain->GetBranch("mm2gamma");
   b_truemom1phpi0 = fChain->GetBranch("truemom1phpi0");
   b_truemom2phpi0 = fChain->GetBranch("truemom2phpi0");
   b_truemomlab1phpi0 = fChain->GetBranch("truemomlab1phpi0");
   b_truemomlab2phpi0 = fChain->GetBranch("truemomlab2phpi0");
   b_trueth1phpi0 = fChain->GetBranch("trueth1phpi0");
   b_trueth2phpi0 = fChain->GetBranch("trueth2phpi0");
   b_mom1phpi0 = fChain->GetBranch("mom1phpi0");
   b_mom2phpi0 = fChain->GetBranch("mom2phpi0");
   b_truemrho = fChain->GetBranch("truemrho");
   b_nrecorho = fChain->GetBranch("nrecorho");
   b_mrho = fChain->GetBranch("mrho");
   b_mm2rho = fChain->GetBranch("mm2rho");
   b_truemompirho = fChain->GetBranch("truemompirho");
   b_truemompi0rho = fChain->GetBranch("truemompi0rho");
   b_mompirho = fChain->GetBranch("mompirho");
   b_mompi0rho = fChain->GetBranch("mompi0rho");
   b_truemrho0 = fChain->GetBranch("truemrho0");
   b_nrecorho0 = fChain->GetBranch("nrecorho0");
   b_mrho0 = fChain->GetBranch("mrho0");
   b_mm2rho0 = fChain->GetBranch("mm2rho0");
   b_truemom1pirho0 = fChain->GetBranch("truemom1pirho0");
   b_truemom2pirho0 = fChain->GetBranch("truemom2pirho0");
   b_mom1pirho0 = fChain->GetBranch("mom1pirho0");
   b_mom2pirho0 = fChain->GetBranch("mom2pirho0");
   b_truemomega = fChain->GetBranch("truemomega");
   b_nrecoomega = fChain->GetBranch("nrecoomega");
   b_momega = fChain->GetBranch("momega");
   b_mm2omega = fChain->GetBranch("mm2omega");
   b_truemom1piome = fChain->GetBranch("truemom1piome");
   b_truemom2piome = fChain->GetBranch("truemom2piome");
   b_truemompi0ome = fChain->GetBranch("truemompi0ome");
   b_truedalitzpi1pi2ome = fChain->GetBranch("truedalitzpi1pi2ome");
   b_truedalitzpi1pi0ome = fChain->GetBranch("truedalitzpi1pi0ome");
   b_truecosthome = fChain->GetBranch("truecosthome");
   b_mom1piome = fChain->GetBranch("mom1piome");
   b_mom2piome = fChain->GetBranch("mom2piome");
   b_mompi0ome = fChain->GetBranch("mompi0ome");
   b_dalitzpi1pi2ome = fChain->GetBranch("dalitzpi1pi2ome");
   b_dalitzpi1pi0ome = fChain->GetBranch("dalitzpi1pi0ome");
   b_costhome = fChain->GetBranch("costhome");
   b_meta = fChain->GetBranch("meta");
   b_mm2eta = fChain->GetBranch("mm2eta");
   b_modeeta = fChain->GetBranch("modeeta");
   b_metap = fChain->GetBranch("metap");
   b_mm2etap = fChain->GetBranch("mm2etap");
   b_modeetap = fChain->GetBranch("modeetap");
   b_etapmassetadau = fChain->GetBranch("etapmassetadau");
   b_momrho0ph = fChain->GetBranch("momrho0ph");
   b_ma0 = fChain->GetBranch("ma0");
   b_mm2a0 = fChain->GetBranch("mm2a0");
   b_a0massetadau = fChain->GetBranch("a0massetadau");
   b_modea0 = fChain->GetBranch("modea0");
   b_ma0p = fChain->GetBranch("ma0p");
   b_mm2a0p = fChain->GetBranch("mm2a0p");
   b_a0pmassetadau = fChain->GetBranch("a0pmassetadau");
   b_modea0p = fChain->GetBranch("modea0p");
   return kTRUE;
}
// ----------------------------------------------------------------------
void thecomp::readCuts(TString filename) {
  char  buffer[200];
  sprintf(buffer, "%s", filename.Data());
  ifstream is(buffer);
  char CutName[100];
  float CutValue;
  int ok(0);
  TOTALSTAT = 270. ;
  TOTALSTATMODEL = 270.;
  BRRATIOGENVALUE = 0.0017 ;
  BRRATIOVALUETAIL = 0.0017 ;
  PSTARFACT = 1.17;
  FITCATEGORY = 11;
  Q2CUT = 0;
  MXCUT = 1.6;
  USECB = 1;
  FIXMEANVALUE = 1;
  FIXSIGMA = 0;
  FIXARGUS1 = 0;
  FIXARGUS2 = 0;
  FIXCB1 = 0;
  FIXCB2 = 0;
  LEPTONPCUT = 1.;
  NLEP = 1;
  DELTAM = 0.1;
  PRMM2CUT = 1000000.;
  MNUSQLOW = -1000.;
  MNUSQHIGH = 0.5;
  CHLOW = -.5;
  CHHIGH = .5;
  DEPL = 1.;
  SSBAR = 0;
  BTYPE = 2;
  LEPTTYPE = 2;
  FITTOTSHAPE = 0;
  MIXCORR = 2;
  FITMC = 1;
  FITOPT = 0;
  TOYMC = 0.;
  MULTIFIT = 0;
  BLINDING = 1;
  RANDOMSEED = 990717;
  BLINDSIZE = 0.51;
  ISSMEARALL = 0;              
  SMEARALLMEANVALUE = 0;
  SMEARALLSIGMA   = 1;      
  ISSMEARBKG      = 0;        
  SMEARBKGMEANVALUE = 0;
  SMEARBKGSIGMA    = 1; 
  DOBRECOWEIGHT = 0;
  DOBDECWEIGHT = 0;
  DODDECWEIGHT = 0;
  DOTRKWEIGHT = 0;
  DONEUWEIGHT = 0;
  MAXINTPUR = 1000.;
  MININTPUR = 0.;
  RUN = 0;
  CUTNNPI0 = -1000;
  DOFERMI = 0;
  DOTHEO = 0;
  FERMIAPP = 0;
  DELTAMB = 0;
  DELTAA = 0;


  MXCUTLOWEXCL=-1000.;
  MXCUTHIGHEXCL=1000.;
  NCHGLOWEXCL=0;
  NCHGHIGHEXCL=1000;
  MOM1MIN=-1000.;
  MOM2MIN=-1000.;
  NCOMBLOWEXCL=0;
  NCOMBHIGHEXCL=1000;
  NPI0LOWEXCL=0;
  NPI0HIGHEXCL=1000;
  DALITZCUT=1.;
  DAUETAMASS=1000.;
  DAUGAMMAMOM = -1000.;
  DORIGHTNCHG = 0;

  MNUSQPI0LOW =  0.;
  MNUSQPI0HIGH = 0.;
  MNUSQETALOW =  0.;
  MNUSQETAHIGH = 0.;
  MNUSQRHOLOW =  0.;
  MNUSQRHOHIGH = 0.;
  MNUSQRHO0LOW =  0.;
  MNUSQRHO0HIGH = 0.;
  MNUSQOMEGALOW =  0.;
  MNUSQOMEGAHIGH = 0.;
  MAXENEU = 1000.;
  JPSIWIN = 1000.;

  while (is.getline(buffer, 200, '\n')) {
    ok = 0;
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);
    // -- 
    if (!strcmp(CutName, "totalStat")) { TOTALSTAT  = CutValue; ok = 1;}
    if (!strcmp(CutName, "totalStatModel")) { TOTALSTATMODEL = CutValue; ok = 1;}
    if (!strcmp(CutName, "BRRatioGenValue")) { BRRATIOGENVALUE  = CutValue; ok = 1;}
    if (!strcmp(CutName, "BRRatioValueTail")) { BRRATIOVALUETAIL  = CutValue; ok = 1;}
    if (!strcmp(CutName, "pstarfact")) { PSTARFACT  = CutValue; ok = 1;}
    if (!strcmp(CutName, "fitcategory")) { FITCATEGORY= int(CutValue); ok = 1;}
    if (!strcmp(CutName, "mxCut")) { MXCUT  = CutValue; ok = 1;}
    if (!strcmp(CutName, "q2Cut")) { Q2CUT = CutValue; ok = 1;}
    if (!strcmp(CutName, "useCB")) { USECB= int(CutValue); ok = 1;}
    if (!strcmp(CutName, "fixMeanValue")) { FIXMEANVALUE= int(CutValue); ok = 1;}
    if (!strcmp(CutName, "fixSigma")) { FIXSIGMA= int(CutValue); ok = 1;}
    if (!strcmp(CutName, "fixArgus1")) { FIXARGUS1= int(CutValue); ok = 1;}
    if (!strcmp(CutName, "fixArgus2")) { FIXARGUS2 = int(CutValue); ok = 1;}
    if (!strcmp(CutName, "fixCB1")) { FIXCB1 = int(CutValue); ok = 1;}
    if (!strcmp(CutName, "fixCB2")) { FIXCB2 = int(CutValue); ok = 1;}
    if (!strcmp(CutName, "leptonPCut")) { LEPTONPCUT = CutValue; ok = 1;}
    if (!strcmp(CutName, "nlep")) { NLEP= int(CutValue); ok = 1;}
    if (!strcmp(CutName, "deltam")) { DELTAM= CutValue; ok = 1;}
    if (!strcmp(CutName, "prmm2Cut")) { PRMM2CUT = CutValue; ok = 1;}
    if (!strcmp(CutName, "mnuSqLow")) { MNUSQLOW= CutValue; ok = 1;}
    if (!strcmp(CutName, "mnuSqHigh")) { MNUSQHIGH= CutValue; ok = 1;}
    if (!strcmp(CutName, "chLow")) { CHLOW= CutValue; ok = 1;}
    if (!strcmp(CutName, "chHigh")) { CHHIGH= CutValue; ok = 1;}
    if (!strcmp(CutName, "depl")) { DEPL = CutValue; ok = 1;}
    if (!strcmp(CutName, "ssbar")) { SSBAR = CutValue; ok = 1;}
    if (!strcmp(CutName, "Btype")) { BTYPE = CutValue; ok = 1;}
    if (!strcmp(CutName, "lepttype")) { LEPTTYPE = CutValue; ok = 1;}
    if (!strcmp(CutName, "fittotshape")) { FITTOTSHAPE= int(CutValue); ok = 1;}
    if (!strcmp(CutName, "mixcorr")) { MIXCORR = int(CutValue); ok = 1;}
    if (!strcmp(CutName, "fitMC")) { FITMC= int(CutValue); ok = 1;}
    if (!strcmp(CutName, "fitOption")) { FITOPT= int(CutValue); ok = 1;}
    if (!strcmp(CutName, "toyMC")) { TOYMC= CutValue; ok = 1;}
    if (!strcmp(CutName, "multifit")) { MULTIFIT = int(CutValue); ok = 1;}
    if (!strcmp(CutName, "blinding")) { BLINDING = int(CutValue); ok = 1;}
    if (!strcmp(CutName, "randomseed")) { RANDOMSEED = int(CutValue); ok = 1;}
    if (!strcmp(CutName, "blindsize")) { BLINDSIZE = CutValue; ok = 1; }
    if (!strcmp(CutName, "issmearAll")) {ISSMEARALL = int(CutValue); ok = 1;}       
    if (!strcmp(CutName, "smearAllMeanValue")) {SMEARALLMEANVALUE = CutValue; ok = 1;} 
    if (!strcmp(CutName, "smearAllSigma")) {SMEARALLSIGMA = CutValue; ok = 1;}     
    if (!strcmp(CutName, "issmearBkg")) {ISSMEARBKG = int(CutValue); ok = 1;}       
    if (!strcmp(CutName, "smearBkgMeanValue")) {SMEARBKGMEANVALUE = CutValue; ok = 1;} 
    if (!strcmp(CutName, "smearBkgSigma")) {SMEARBKGSIGMA = CutValue; ok = 1;}    
    if (!strcmp(CutName, "dobrecoreweight")) { DOBRECOWEIGHT = int(CutValue); ok = 1;}
    if (!strcmp(CutName, "doBdecreweight")) { DOBDECWEIGHT = int(CutValue); ok = 1;}
    if (!strcmp(CutName, "doDdecreweight")) { DODDECWEIGHT = int(CutValue); ok = 1;}
    if (!strcmp(CutName, "dotrkreweight")) { DOTRKWEIGHT = int(CutValue); ok = 1;}
    if (!strcmp(CutName, "doneureweight")) { DONEUWEIGHT = int(CutValue); ok = 1;}
    if (!strcmp(CutName, "maxintpur")) {MAXINTPUR = CutValue; ok = 1;}
    if (!strcmp(CutName, "minintpur")) {MININTPUR = CutValue; ok = 1;}
    if (!strcmp(CutName, "run")) {RUN = CutValue; ok = 1;}
    if (!strcmp(CutName, "nnpi0")) {CUTNNPI0 = int(CutValue); ok = 1;}
    if (!strcmp(CutName, "fermiapp")) {FERMIAPP = CutValue; ok = 1;}
    if (!strcmp(CutName, "deltamb")) {DELTAMB = CutValue; ok = 1;}
    if (!strcmp(CutName, "deltaa")) {DELTAA = CutValue; ok = 1;}

    if (!strcmp(CutName, "mxCutLowExcl")) { MXCUTLOWEXCL  = CutValue; ok = 1;}
    if (!strcmp(CutName, "mxCutHighExcl")) { MXCUTHIGHEXCL  = CutValue; ok = 1;}
    if (!strcmp(CutName, "nchgLowExcl")) { NCHGLOWEXCL = CutValue; ok = 1;}
    if (!strcmp(CutName, "nchgHighExcl")) { NCHGHIGHEXCL = CutValue; ok = 1;}
    if (!strcmp(CutName, "ncombLowExcl")) { NCOMBLOWEXCL = CutValue; ok = 1;}
    if (!strcmp(CutName, "ncombHighExcl")) { NCOMBHIGHEXCL = CutValue; ok = 1;}
    if (!strcmp(CutName, "npi0LowExcl")) { NPI0LOWEXCL = CutValue; ok = 1;}
    if (!strcmp(CutName, "npi0HighExcl")) { NPI0HIGHEXCL = CutValue; ok = 1;}
    if (!strcmp(CutName, "mom1min")) { MOM1MIN = CutValue; ok = 1;}
    if (!strcmp(CutName, "mom2min")) { MOM2MIN = CutValue; ok = 1;}
    if (!strcmp(CutName, "dalitzcut")) { DALITZCUT = CutValue; ok = 1;}
    if (!strcmp(CutName, "dauetamass")) { DAUETAMASS = CutValue; ok = 1;}    
    if (!strcmp(CutName, "daugammamom")) { DAUGAMMAMOM = CutValue; ok = 1;}    
    if (!strcmp(CutName, "dorightnchg")) { DORIGHTNCHG = CutValue; ok = 1;}    

    if (!strcmp(CutName, "mnuSqpi0Low")) { MNUSQPI0LOW= CutValue; ok = 1;}
    if (!strcmp(CutName, "mnuSqpi0High")) { MNUSQPI0HIGH= CutValue; ok = 1;}
    if (!strcmp(CutName, "mnuSqetaLow")) { MNUSQETALOW= CutValue; ok = 1;}
    if (!strcmp(CutName, "mnuSqetaHigh")) { MNUSQETAHIGH= CutValue; ok = 1;}
    if (!strcmp(CutName, "mnuSqrhoLow")) { MNUSQRHOLOW= CutValue; ok = 1;}
    if (!strcmp(CutName, "mnuSqrhoHigh")) { MNUSQRHOHIGH= CutValue; ok = 1;}
    if (!strcmp(CutName, "mnuSqrho0Low")) { MNUSQRHO0LOW= CutValue; ok = 1;}
    if (!strcmp(CutName, "mnuSqrho0High")) { MNUSQRHO0HIGH= CutValue; ok = 1;}
    if (!strcmp(CutName, "mnuSqomegaLow")) { MNUSQOMEGALOW= CutValue; ok = 1;}
    if (!strcmp(CutName, "mnuSqomegaHigh")) { MNUSQOMEGAHIGH= CutValue; ok = 1;}
    if (!strcmp(CutName, "maxeneu")) {MAXENEU  = CutValue; ok = 1;}
    if (!strcmp(CutName, "jpsiwindow")) {JPSIWIN  = CutValue; ok = 1;}

    if (ok == 0)  cout << "==> exclfitNtp::readCuts() Error: Don't know about variable " << CutName << endl;
  }
  dumpCuts();

}

// ----------------------------------------------------------------------
void thecomp::dumpCuts() {
  cout << "====================================" << endl;
  cout << " fitcategory         :"  << FITCATEGORY   <<  endl;  
  cout << " q2Cut               :"  << Q2CUT   <<  endl;  
  cout << " leptonPCut          :"  << LEPTONPCUT   <<  endl;  
  cout << " deltam              :"  << DELTAM  <<  endl;  
  cout << " prmm2Cut            :"  << PRMM2CUT   <<  endl;  
  cout << " mnuSqLow            :"  << MNUSQLOW   <<  endl;  
  cout << " mnuSqHigh           :"  << MNUSQHIGH   <<  endl;  
  cout << " chLow               :"  << CHLOW   <<  endl;  
  cout << " chHigh              :"  << CHHIGH   <<  endl;  
  cout << " Btype               :"  << BTYPE   <<  endl;  
  cout << " lepttype            :"  << LEPTTYPE   <<  endl;  
  cout << " maxintpur           :"  << MAXINTPUR    <<  endl;
  cout << " minintpur           :"  << MININTPUR    <<  endl;
  cout << " run                 :"  << RUN    <<  endl;
  cout << " mxCutLowExcl        :"  << MXCUTLOWEXCL   <<  endl;  
  cout << " mxCutHighExcl       :"  << MXCUTHIGHEXCL   <<  endl;  
  cout << " nchgLowExcl         :"  << NCHGLOWEXCL   <<  endl;  
  cout << " nchgHighExcl        :"  << NCHGHIGHEXCL   <<  endl;  
  cout << " ncombLowExcl        :"  << NCOMBLOWEXCL   <<  endl;  
  cout << " ncombHighExcl       :"  << NCOMBHIGHEXCL   <<  endl;  
  cout << " npi0LowExcl         :"  << NPI0LOWEXCL   <<  endl;  
  cout << " npi0HighExcl        :"  << NPI0HIGHEXCL   <<  endl;  
  cout << " mom1min             :"  << MOM1MIN   <<  endl;  
  cout << " mom2min             :"  << MOM2MIN   <<  endl;  
  cout << " dalitzcut           :"  << DALITZCUT   <<  endl;  
  cout << " dauetamass          :"  << DAUETAMASS   <<  endl;  
  cout << " daugammamom         :"  << DAUGAMMAMOM   <<  endl;  
  cout << " dorightnchg         :"  << DORIGHTNCHG   <<  endl;  
  cout << " mnuSqpi0Low         :"  << MNUSQPI0LOW   <<  endl;  
  cout << " mnuSqpi0High        :"  << MNUSQPI0HIGH   <<  endl;  
  cout << " mnuSqetaLow         :"  << MNUSQETALOW   <<  endl;  
  cout << " mnuSqetaHigh        :"  << MNUSQETAHIGH   <<  endl;  
  cout << " mnuSqrhoLow         :"  << MNUSQRHOLOW   <<  endl;  
  cout << " mnuSqrhoHigh        :"  << MNUSQRHOHIGH   <<  endl;  
  cout << " mnuSqrho0Low        :"  << MNUSQRHO0LOW   <<  endl;  
  cout << " mnuSqrho0High       :"  << MNUSQRHO0HIGH   <<  endl;  
  cout << " mnuSqomegaLow       :"  << MNUSQOMEGALOW   <<  endl;  
  cout << " mnuSqomegaHigh      :"  << MNUSQOMEGAHIGH   <<  endl;  
  cout << " maxeneu             :"  << MAXENEU   <<  endl;  
  cout << " jpsiwindow          :"  << JPSIWIN   <<  endl;  
  cout << "====================================" << endl;
  cout << "" << endl;
  cout << "" << endl;
}


thecomp::~thecomp()
{
  fHistFile->cd();
  fHistFile->Write();
  fHistFile->Close();
  delete fHistFile;
  //delete Dvar;
  //delete Bsem;
   if (!fChain) return;
   //   delete fChain->GetCurrentFile();
}

Int_t thecomp::GetEntry(Int_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Int_t thecomp::LoadTree(Int_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Int_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->IsA() != TChain::Class()) return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void thecomp::Show(Int_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t thecomp::Cut(Int_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

