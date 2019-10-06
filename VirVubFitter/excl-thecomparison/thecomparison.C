#include "thecomparison.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include <TVector2.h>
#include "util.hh"
#include "util.cc"
#include "TH1D.h"

ClassImp(thecomparison)
  
thecomparison::thecomparison()
{
  themin=0;
  themax=0;
  thebins=0;
  thevar="";
  Dvar=NULL;
  Bsem=NULL;
}

thecomparison::thecomparison(char *var, double mi, double ma, int b)
{
   char name[100];
   themin = mi;
   themax = ma;
   thebins = b;
   thevar = TString(var) ;
   Dvar = new recoilDSys("ddecay.table.CM2",0,2);
   Bsem = new recoilDSys(0);
}
 
void thecomparison::Loop(int nevents, int cat, double shift, double smear )
{

  char name[100];
  char le[100];
  //id = Mx category == Mx bin         
  int id(0);
  int idflav;
  //  double thecat =0;
  
  int group;  
  
  SHIFTNEUT = shift;
  SIGMANEUT = smear;   
      
  if (fChain == 0) return;
  
  Int_t nentries = Int_t(fChain->GetEntries());
  
  if( nentries > nevents) nentries = nevents;
  
  cout <<"Nentries = "<< nentries << endl;
  
  Int_t nbytes = 0, nb = 0;
  for (Int_t jentry=0; jentry<nentries;jentry++) {
        //if(cat==5)cout <<" entry # "<<jentry;
    Int_t ientry = LoadTree(jentry); //in case of a TChain, ientry is the entry number in the current file
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    // vub - vcb selection
    //    if(bsel<0)
    //  if( (!vub && cat==4) || (!vcb && cat==5) )continue;
    
    int flav =  lcharge + brecoflav; // charge correlation
    bool ksele = nkp + nks;          // fit on the depleted sample?
    
    int ch = xcharge + brecocharge;  // total charge
    double cutflav = 0;              // flavor category (3 = charged B, 4 = neutral B OS, 5 = neutral B SS)  
    // int ibch = 1;                     
    
    // cuts
    
    bool isLp = ((nlept500> 0) && (tlab<2.37) && (tlab>0.36) && (plab>0.5));       
    bool PCMS = pcms>LEPTONPCUT;
    if(isele == 1) PCMS = pcms > ELECTRONPCUT;   
    if(isele == 0) PCMS = pcms > MUONPCUT; 
    
    bool NLE = nlept500 <=1;
    bool CH = (ch >CHLOW  && ch < CHHIGH);
    bool KSELE = ksele;
    
    int FLAV = !(TMath::Abs(brecocharge)!=0 && (lcharge + brecoflav)!=0);
    bool BCH = TMath::Abs(brecocharge) == BTYPE;     
    bool IPUR = intpur>MININTPUR && intpur<MAXINTPUR; 
    bool NCHG = nchg-1 >= NCHGLOWEXCL && nchg-1 <= NCHGHIGHEXCL;       
    bool OTHER(1);
    if(FITCATEGORY == -12) OTHER = nrecoEta>0;
    if(FITCATEGORY == -15) OTHER = nrecoEtap>0;    
    
    bool MOM1(1), MOM2(1);
    bool COMB = 1;
    bool PIZ = 1;
    bool DALITZ = 1;
    bool WdeltaCut = (wdeltam>PRMM2CUT && brecocharge == 0); //last cut based on wdelta                  
    bool DAUETA = 1;
    bool DAURHO = 1;
    bool DAUGAMMOM = 1; 

    double themm2=-999;
    bool MM2 = (themm2 < MNUSQHIGH && themm2 > MNUSQLOW);  
    if(FITCATEGORY == -12) {themm2 = mm2bestEta; 
    MM2 = (mm2bestEta < MNUSQHIGH && mm2bestEta > MNUSQLOW);}  
    if(FITCATEGORY == -15) {themm2 = mm2bestEtap;
    MM2 = (mm2bestEtap < MNUSQHIGH && mm2bestEtap > MNUSQLOW);} 

    double mass = -999;  
    bool MX = (mass < MXCUTHIGHEXCL && mass > MXCUTLOWEXCL);
    
    if(FITCATEGORY == -12){
      mass = barembestEta;   
      //      MX = (barembestEta < MXCUTHIGHEXCL && barembestEta > MXCUTLOWEXCL);
      if(modeEta[indexbestEta]==1) MX = (barembestEta < MXCUTHIGHEXCL1 && barembestEta > MXCUTLOWEXCL1);
      if(modeEta[indexbestEta]==2) MX = (barembestEta < MXCUTHIGHEXCL2 && barembestEta > MXCUTLOWEXCL2);
      if(modeEta[indexbestEta]==3) MX = (barembestEta < MXCUTHIGHEXCL3 && barembestEta > MXCUTLOWEXCL3);
      
      if(DORIGHTNCHG){
        NCHG = 0;
        if(modeEta[indexbestEta] == 2) { NCHG = nchg-1 == 2; }   
        else { NCHG = nchg-1 == 0; }
      }
    }         
    if(FITCATEGORY == -15){
      mass = barembestEtap;   
      //      MX = (barembestEtap < MXCUTHIGHEXCL && barembestEtap > MXCUTLOWEXCL);
      
      if(modeEtap[indexbestEtap]==1) MX = (barembestEtap < MXCUTHIGHEXCL1 && barembestEtap > MXCUTLOWEXCL1);
      if(modeEtap[indexbestEtap]==2) MX = (barembestEtap < MXCUTHIGHEXCL2 && barembestEtap > MXCUTLOWEXCL2);
      if(modeEtap[indexbestEtap]==3) MX = (barembestEtap < MXCUTHIGHEXCL3 && barembestEtap > MXCUTLOWEXCL3);
      if(modeEtap[indexbestEtap]==4) MX = (barembestEtap < MXCUTHIGHEXCL4 && barembestEtap > MXCUTLOWEXCL4);
      DAUGAMMOM = GammamomdauEtap[indexbestEtap]>DAUGAMMAMOM || GammamomdauEtap[indexbestEtap]<-10.;
      if(modeEtap[indexbestEtap] == 1) {
        DAURHO = TMath::Abs(Rho0massdauEtap[indexbestEtap]-0.775) < DAURHOMASS;
      }else{
        DAUETA = TMath::Abs(EtamassdauEtap[indexbestEtap]-0.54775)<DAUETAMASS;
	if(modeEtap[indexbestEtap]==2) DAUETA= TMath::Abs(EtamassdauEtap[indexbestEtap]-0.54775)<DAUETAMASS2;
	if(modeEtap[indexbestEtap]==3) DAUETA= TMath::Abs(EtamassdauEtap[indexbestEtap]-0.54775)<DAUETAMASS3; 
	if(modeEtap[indexbestEtap]==4) DAUETA= TMath::Abs(EtamassdauEtap[indexbestEtap]-0.54775)<DAUETAMASS4;
      }
      if(DORIGHTNCHG){
        NCHG = 0;
        if(modeEtap[indexbestEtap] == 3) { NCHG = nchg-1 == 4; }
        else { NCHG = nchg-1 == 2; }
      }
    }    
    
    // mm2 cut in order to remove Vub crossfeed   
    bool MM2CROSS = 1;     
    bool MM2PI0 = (mm2bestPi0 >= MNUSQPI0HIGH || mm2bestPi0 <= MNUSQPI0LOW);
    bool MM2ETA = (mm2bestEta >= MNUSQETAHIGH || mm2bestEta <= MNUSQETALOW);   
    MM2CROSS = MM2PI0*MM2ETA; 

    if (BTYPE == 2) BCH = 1;        
    
    //    bool CAT = thecat;
    
    double myvar = -1;
    //eta 
    if (thevar == "barembestEta") {myvar = barembestEta; MX=1; MM2=1;}
    if (thevar == "modeEta") {myvar = modeEta[indexbestEta];MM2=1;}
    if (thevar == "mm2bestEta"){myvar = mm2bestEta; MM2=1; }
    if (thevar == "pcmsEta") {myvar = pcms; PCMS=1;MM2=1; }
    if (thevar == "nchgEta") {myvar = nchg; NCHG=1; CH=1;MM2=1;}
//    if (thevar == "Gvxbtypeta") {myvar = Gvxbtyp;}
    //etap
    if (thevar == "pcmsEtap") {myvar = pcms; PCMS=1; MM2=1;}
    if (thevar == "nchgEtap") {myvar = nchg; NCHG=1; CH=1;MM2=1;}
    if (thevar == "barembestEtap") {myvar = barembestEtap; MX=1; MM2=1;}
    if (thevar == "modeEtap") {myvar = modeEtap[indexbestEtap];MM2=1;}
    if (thevar == "mm2bestEtap") {myvar = mm2bestEtap; MM2=1;}
    if (thevar == "EtamassdauEtap") {myvar = EtamassdauEtap[indexbestEtap];DAUETA=1;MM2=1;}
    if (thevar == "GammamomdauEtap") {myvar = GammamomdauEtap[indexbestEtap];DAUGAMMOM=1;MM2=1;}
    if (thevar == "Rho0massdauEtap") {myvar = Rho0massdauEtap[indexbestEtap];DAURHO=1;MM2=1;}

    //if (thevar == "kmin") {if(nkp<=0)continue;myvar = kminmom;}
    //if (thevar == "kmax") {if(nkp<=0)continue;myvar = kmaxmom;}
    if (!((myvar>0) || (myvar<0) || (myvar == 0))) {
      myvar = -999.;
    } 
    double totweight = 1;
    
    if(PCMS && NLE && FLAV && BCH && IPUR && OTHER && isLp) {
      sprintf(le, "h");
      // assign flavor category
      cutflav = 5;
      if(TMath::Abs(brecocharge)) cutflav = 3;
      if(TMath::Abs(brecocharge)==0 && flav==0) cutflav = 4;	  
      id = hist(myvar);
      idflav =(int)(hist(myvar) + cutflav * 10000);    
      group = 200 + cat * 1000 + 100000;
      sprintf(name, "%s%d",le,group+idflav);
      if (id!=-1)((TH1D*)gDirectory->Get(name))->Fill(mes,totweight);
      group = 90206 + cat * 1000 + 100000;
      sprintf(name, "%s%d",le,group);
      if (id!=-1)((TH1D*)gDirectory->Get(name))->Fill(mes,totweight);   
     
      if((!KSELE) && MM2 && CH && NCHG && (!WdeltaCut) && MX && MOM1 && MOM2 && COMB && PIZ && DALITZ && MM2CROSS && 
DAUETA && DAURHO && DAUGAMMOM){  
      	group = 200 + cat * 1000;
	sprintf(name, "%s%d",le,group+idflav);
        if (id!=-1)((TH1D*)gDirectory->Get(name))->Fill(mes,totweight);  	
	group = 90206 + cat * 1000;
	sprintf(name, "%s%d",le,group);
	if (id!=-1)((TH1D*)gDirectory->Get(name))->Fill(mes,totweight);  
      }
    }
  }
}
// ----------------------------------------------------------------------
double thecomparison::getBsysweight(int decType,int thevub) {
  double theweight;
  theweight = 1.;  
  theweight = Bsem->weight(decType); 
  if(thevub) theweight = 1.;
  return theweight;
}
// ----------------------------------------------------------------------
double thecomparison::getDsysweight(int decDpi,int decDk,int decDks,int decDpiz,int decDlep,int decImode,int thevub) {
  double theweight;
  theweight = 1.;  
  int bla(1);
  theweight = Dvar->weight(decDpi,decDk,decDks,decDpiz,decDlep,decImode,bla);    
  if(thevub) theweight = 1.;
  return theweight;
}

double thecomparison::FermiWeight(double kp, double deltamb, double deltaa){

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
double thecomparison::fermi(double kp, double m, double a) {
  double BMASS   = 5.2792;
  double x = kp/(BMASS - m);
  if ((kp>-m) && (x <= 1.)) {
    return TMath::Power((1-x), a) * TMath::Exp((1+a)*x); 
  } 
  return 0.;
}


void thecomparison::Bookhist()
{
  //  gROOT->cd();
  TH1 *h;
  char name[100], title[100],  number[100];
  //Int_t nbins = 13;
  //Int_t nbins2 = 10;
  //char  preine[100], preich[100], prefix[100]
  //int ine; int ich; 
  int lo; 
  
  sprintf(name, "h400000");  sprintf(title, "%s%s", thevar.Data(), " data events after all cuts");  h = new TH1D(name, title, thebins, themin, themax );  h->Sumw2();
  sprintf(name, "h500000");  sprintf(title, "%s%s", thevar.Data(), " MC events  after all cuts");  h = new TH1D(name, title, thebins, themin, themax  );  h->Sumw2();
  sprintf(name, "h94206");  sprintf(title, "mes data after all cuts");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  sprintf(name, "h95206");  sprintf(title, "mes MC after all cuts");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  sprintf(name, "h1400000");  sprintf(title, "%s%s", thevar.Data(), " data events after lepton cuts");  h = new TH1D(name, title, thebins, themin, themax );  h->Sumw2();
  sprintf(name, "h1500000");  sprintf(title, "%s%s", thevar.Data(), " MC events after lepton cuts");  h = new TH1D(name, title, thebins, themin, themax  );  h->Sumw2();
  sprintf(name, "h194206");  sprintf(title, "mes data after lepton cuts");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  sprintf(name, "h195206");  sprintf(title, "mes MC after lepton cuts");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  
  for (lo=1;lo<thebins+1;lo++) {
    
    sprintf(number, "%d",  200+lo);     
    sprintf(name,"%s%s" , "h4",number);  sprintf(title, "mes data");  h = new TH1D(name, title, 40, 5.2, 5.3);         h->Sumw2();
    sprintf(name,"%s%s" , "h34",number);  sprintf(title, "mes data bch");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h44",number);  sprintf(title, "mes data b0 os");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h54",number);  sprintf(title, "mes data b0 ss");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h5",number);  sprintf(title, "mes MC");  h = new TH1D(name, title, 40, 5.2, 5.3);         h->Sumw2();
    sprintf(name,"%s%s" , "h35",number);  sprintf(title, "mes MC bch");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h45",number);  sprintf(title, "mes MC b0 os");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h55",number);  sprintf(title, "mes MC b0 ss");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();

    sprintf(name,"%s%s" , "h14",number);  sprintf(title, "mes data");  h = new TH1D(name, title, 40, 5.2, 5.3);         h->Sumw2();
    sprintf(name,"%s%s" , "h134",number);  sprintf(title, "mes data bch");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h144",number);  sprintf(title, "mes data b0 os");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h154",number);  sprintf(title, "mes data b0 ss");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h15",number);  sprintf(title, "mes MC");  h = new TH1D(name, title, 40, 5.2, 5.3);         h->Sumw2();
    sprintf(name,"%s%s" , "h135",number);  sprintf(title, "mes MC bch");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h145",number);  sprintf(title, "mes MC b0 os");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h155",number);  sprintf(title, "mes MC b0 ss");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  }
   sprintf(name, "h8888");  sprintf(title, "test gauss");  h = new TH1D(name, title, 100, -.05, .05);    h->Sumw2(); 
}

int thecomparison::hist(double mx){

  // categories
   int bin = int((mx-themin)/((themax-themin)/thebins)+1);
   if (mx==themax) bin = thebins;
   if (mx==themin) bin = 1;
   if (mx<themin || mx>themax) bin = -1;
   return bin;
}

void thecomparison::Fitmes(int cat, int cut){

  // fit to mes distribution and fill of the Mx plots...

  //  gROOT->cd();

  int group = 200 + cat * 1000;
   
  //  char nameps[100], preine[100];
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
  //int ia;
  int is;
  //int ich;
  //int ine;
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
  int addcut = 0; 
  if(cut) addcut = 1;
  int const theb = thebins + 1;
  for (int y=0; y<2; y++) {
    sprintf(le, "h");
//    if (y == 1) sprintf(le, "d");
    
    //extracting the signal fit parameters 
    sprintf(name, "%s%d",le,group+90006+addcut*100000);  
    cout << ((TH1D*)gDirectory->Get(name))->Integral() << endl;
    double dummy1,dummy2;
    sighisto(dummy1,dummy2,(TH1D*)gDirectory->Get(name),usedmean,usedsigma,usedalpha,usedn,1,-11111111.,-1111111.,-1111111.,-1111111.,-1111111.);
    cout << "mes result for: " << name << " is, MEAN " << usedmean << " SIGMA " << usedsigma << " ALPHA " << usedalpha << " N " << usedn << endl;
  
    group = thegroup;
    for (i = 1;  i < theb; i++) { 
      is = i;
      sprintf(name, "%s%d",le,group+is+addcut*10000);  
      
      sprintf(namebch, "%s%d",le,group+is+30000+addcut*100000);  
      sprintf(nameb0os, "%s%d",le,group+is+40000+addcut*100000);  
      sprintf(nameb0ss, "%s%d",le,group+is+50000+addcut*100000);  
    
    // mixing correction      
    
      for(int k=1;k<40;k++){
	tempbinchb = ((TH1D*)gDirectory->Get(namebch))->GetBinContent(k);
	tempbinb0os = ((TH1D*)gDirectory->Get(nameb0os))->GetBinContent(k);
	tempbinb0ss = ((TH1D*)gDirectory->Get(nameb0ss))->GetBinContent(k);
	temperrchb = ((TH1D*)gDirectory->Get(namebch))->GetBinError(k);
	temperrb0os = ((TH1D*)gDirectory->Get(nameb0os))->GetBinError(k);
	temperrb0ss = ((TH1D*)gDirectory->Get(nameb0ss))->GetBinError(k);
	tempbin = tempbinchb + ((1-chid)/(1-2*chid)) * tempbinb0os - (chid/(1-2*chid)) * tempbinb0ss;
	temperr = sqrt(temperrchb*temperrchb + ((1-chid)/(1-2*chid)) * ((1-chid)/(1-2*chid)) * tempbinb0os + (chid/(1-2*chid)) * (chid/(1-2*chid))* tempbinb0ss);
	((TH1D*)gDirectory->Get(name))->SetBinContent(k, tempbin);
	((TH1D*)gDirectory->Get(name))->SetBinError(k, temperr); 
      }	
    
    // put the result of the fit in each Mx bin
    
      sighisto(sigs[is-1],errsigs[is-1],(TH1D*)gDirectory->Get(name),foomean,foosigma,fooalpha,foon,1,usedmean,usedsigma,usedalpha,usedn,-1111111);
      int title = cat * 100000 + addcut * 1000000;
      sprintf(name, "%s%d",le,title);  
      ((TH1D*)gDirectory->Get(name))->SetBinContent(is, sigs[is-1]);
      ((TH1D*)gDirectory->Get(name))->SetBinError(is, errsigs[is-1]);
    }
  }  
}

void thecomparison::sighisto(double&signal, double& signalErr, TH1D *histo, double &resmean, double &ressigma, double &resalpha, double &resn, int fixpar, double mean, double sigma, double alpha,  double n, double argus){

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

void thecomparison::overlap(int cut, int norm, TString dir){
   
   char mycut[100];
   int addcut = 0; 
   if(cut) addcut = 1;
   
   //   c6 = new TCanvas("c6", "c6", 300, 0, 400,800);  
   c6 = new TCanvas("c6", "c6");  
      // -- top left
      fPads[1]= new TPad("pad1", "", 0.00, 0.30, 0.99, 0.99);   fPads[1]->Draw(); 
      fPads[2]= new TPad("pad2", "", 0.00, 0.00, 0.99, 0.29);   fPads[2]->Draw();  
   // -- bottom left
//   fPads[3]= new TPad("pad3", "", 0.00, 0.15, 0.99, 0.49);   fPads[3]->Draw(); 
//   fPads[4]= new TPad("pad4", "", 0.00, 0.00, 0.99, 0.14);   fPads[4]->Draw(); 
   char name[100], line[100];
     //char name2[100], preine[100], preich[100], prefix[100], shift[100], smear[100], 
     //int h = 0 ; 

   sprintf(name,"%s%d" ,"h", 400000 + addcut * 1000000);
   TH1D y1;
   ((TH1D*)gDirectory->Get(name))->Copy(y1);
   sprintf(name,"%s%d" ,"h", 500000 + addcut * 1000000);
   TH1D y2;
   ((TH1D*)gDirectory->Get(name))->Copy(y2);
 //   sprintf(name,"%s%d" ,"d", 400000 + addcut * 1000000);
//    TH1D y3;
//    ((TH1D*)gDirectory->Get(name))->Copy(y3);
//    sprintf(name,"%s%d" ,"d", 500000 + addcut * 1000000);
//    TH1D y4;
//    ((TH1D*)gDirectory->Get(name))->Copy(y4);

   fPads[1]->cd(); shrinkPad(0.001, 0.2); 

   setFilledHist(&y1, kBlack, kBlue, 3004);
   double inter1;
   double max1 = 1.4*y1.GetMaximum();
   y1.SetLabelSize(0.07, "Y");
   y1.SetMaximum(max1);
   y1.SetMinimum(0.);
   y1.SetMarkerSize(0.5);
   inter1 = y1.Integral();
// y1.SetNormFactor(inter1);
   y1.SetLineColor(kRed);
   y1.Draw();
 
   double inter2;
   inter2 = y2.Integral();
   y2.SetMarkerSize(0.5);
   if(cut) {
      intdataen = inter1;
      intMCen = inter2;
   }
   if (norm) {
      if(intMCen>0)y2.Scale(intdataen/intMCen);
   }else{
      if(inter2>0)y2.Scale(inter1/inter2);
   }   
   setFilledHist(&y2 , kBlack, kRed, 3005);
   y2.DrawCopy("histsame");

   double g1 = chisq(&y1,&y2);
   sprintf(line, "#chi^{2} = %5.4f", g1);

   TLatex tl;
   tl.SetNDC(kTRUE);
   tl.SetTextSizePixels(50); tl.DrawLatex(0.22, 0.75, line);
  
    fPads[2]->cd(); shrinkPad(0.001, 0.2); 
    shrinkPad(0.4, 0.2, 0.1, 0.001);

   gPad->SetGridx(1);  gPad->SetGridy(1);
   TH1D *hratio = new TH1D(y1); hratio->SetName("hratio"); hratio->Reset();
   hratio->Divide(&y1, &y2);
   hratio->SetMinimum(0.5); hratio->SetMaximum(1.5);
   hratio->SetMarkerStyle(24);
   hratio->SetNdivisions(504, "Y");
   hratio->SetLabelSize(0.22, "X");  hratio->SetLabelSize(0.17, "Y");
   hratio->SetStats(0);
   hratio->SetTitle("");
   hratio->Draw();

//    fPads[3]->cd(); shrinkPad(0.001, 0.2); 

//    setFilledHist(&y3, kBlack, kBlue, 3004);
//    double inter3;
//    double max2 = 1.4*y3.GetMaximum();
//    y3.SetLabelSize(0.07, "Y");
//    y3.SetMaximum(max2);
//    y3.SetMinimum(0.);
//    y3.SetMarkerSize(0.5);
//    inter3 = y3.Integral();
// // y3.SetNormFactor(inter1);
//    y3.SetLineColor(kRed);
//    y3.Draw();
 
//    double inter4;
//    inter4 = y4.Integral();
//    y4.SetMarkerSize(0.5);
//    if(cut) {
//       intdatadepl = inter3;
//       intMCdepl = inter4;
//    }
//    if (norm) {
//       y4.Scale(intdatadepl/intMCdepl);
//    }else{
//       y4.Scale(inter3/inter4);
//    }   

//    setFilledHist(&y4 , kBlack, kRed, 3005);
//    y4.DrawCopy("histsame");

//    double g2 = chisq(&y3,&y4);
//    sprintf(line, "#chi^{2} = %5.4f", g2); 
//    tl.SetTextSizePixels(50); tl.DrawLatex(0.22, 0.75, line);
  
//    fPads[4]->cd(); shrinkPad(0.001, 0.2); 
//    shrinkPad(0.4, 0.2, 0.1, 0.001);

//    gPad->SetGridx(1);  gPad->SetGridy(1);
//    TH1D *hratio2 = new TH1D(y4); hratio2->SetName("hratio2"); hratio2->Reset();
//    hratio2->Divide(&y3, &y4);
//    hratio2->SetMinimum(0.5); hratio2->SetMaximum(1.5);
//    hratio2->SetMarkerStyle(24);
//    hratio2->SetNdivisions(504, "Y");
//    hratio2->SetLabelSize(0.22, "X");  hratio2->SetLabelSize(0.17, "Y");
//    hratio2->SetStats(0);
//    hratio2->SetTitle("");
//    hratio2->Draw();

//    cout << "the chi square is " << g << endl;
//    SHIFTNEUT = SHIFTNEUT *1000;
//    SIGMANEUT = SIGMANEUT *1000;
//    sprintf(shift, "%s%d","-", SHIFTNEUT);
//    sprintf(smear, "%s%d","-", SIGMANEUT);   
  
   sprintf(mycut, "allcuts");   
   if (cut) sprintf(mycut, "leptoncuts");   
   if (norm) {   
      sprintf(name, "%s%s%d%s%s%s", dir.Data(), "/comparisonnorm",multipl,thevar.Data(),mycut,".eps");
   }else{
      sprintf(name, "%s%s%d%s%s%s", dir.Data(), "/comparison",multipl,thevar.Data(),mycut,".eps");
   }
   c6->SaveAs(name);  
//   sprintf(name, "%s%s%s%s%s%s%s", "comparison",thevar,"smear",smear,"shift",shift,".dat");
//    ofstream outfile(name);   
//    outfile << "the chi square for " << thevar << " with smear " << smear << " and shift " << shift << " is " << endl;
//    outfile << "  CHISQ "  << g << endl; 
}

double thecomparison::chisq(TH1 *h1, TH1 *h2){
  double chisq = 0; 
  double tempchisq = 0;
  double inter;
  double temperr1;
  double temperr2;
  inter = h1->Integral();
  double inter2;
  inter2 = h2->Integral();
  int nbins = h1->GetNbinsX();
  for(int i=1;i<nbins+1;i++) {
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

// void thecomparison::effplots(TString dir){

//    char name[100];
//    sprintf(name,"%s%d" ,"h", 400000);
//    TH1D *y1(((TH1D*)gDirectory->Get(name)));
//    sprintf(name,"%s%d" ,"h", 1400000);
//    TH1D *y2(((TH1D*)gDirectory->Get(name)));
//    sprintf(name,"%s%d" ,"d", 400000);
//    TH1D *y3(((TH1D*)gDirectory->Get(name)));
//    sprintf(name,"%s%d" ,"d", 1400000);
//    TH1D *y4(((TH1D*)gDirectory->Get(name)));
//    sprintf(name,"%s%d" ,"h", 500000);
//    TH1D *y5(((TH1D*)gDirectory->Get(name)));
//    sprintf(name,"%s%d" ,"h", 1500000);
//    TH1D *y6(((TH1D*)gDirectory->Get(name)));
//    sprintf(name,"%s%d" ,"d", 500000);
//    TH1D *y7(((TH1D*)gDirectory->Get(name)));
//    sprintf(name,"%s%d" ,"d", 1500000);
//    TH1D *y8(((TH1D*)gDirectory->Get(name)));

//    TH1D *hratio1 = new TH1D(*y1); hratio1->SetName("hratio1"); hratio1->Reset();
//    TH1D *hratio2 = new TH1D(*y5); hratio2->SetName("hratio2"); hratio2->Reset();
//    TH1D *hratio3 = new TH1D(*y3); hratio3->SetName("hratio3"); hratio3->Reset();
//    TH1D *hratio4 = new TH1D(*y7); hratio4->SetName("hratio4"); hratio4->Reset();
//    TH1D *hratioratio1 = new TH1D(*hratio1); hratioratio1->SetName("hratioratio1"); hratioratio1->Reset();
//    TH1D *hratioratio2 = new TH1D(*hratio3); hratioratio2->SetName("hratioratio2"); hratioratio2->Reset();

  
//    hratio1->Divide(y1, y2, 1, 1, "B");
//    hratio2->Divide(y5, y6, 1, 1, "B");
//    hratio3->Divide(y3, y4, 1, 1, "B");
//    hratio4->Divide(y7, y8, 1, 1, "B");
//    hratioratio1->Divide(hratio1, hratio2);
//    hratioratio2->Divide(hratio3, hratio4);

//    c7 = new TCanvas("c7", "c7", 300,   0, 400,800);  
//    // -- top left
//    fPads[1]= new TPad("pad1", "", 0.00, 0.65, 0.99, 0.99);   fPads[1]->Draw(); 
//    fPads[2]= new TPad("pad2", "", 0.00, 0.50, 0.99, 0.64);   fPads[2]->Draw();  
//    // -- bottom left
//    fPads[3]= new TPad("pad3", "", 0.00, 0.15, 0.99, 0.49);   fPads[3]->Draw(); 
//    fPads[4]= new TPad("pad4", "", 0.00, 0.00, 0.99, 0.14);   fPads[4]->Draw(); 

//    fPads[1]->cd(); shrinkPad(0.001, 0.2); 

// //   double max = 1.4*hratio1->GetMaximum();
//    hratio1->SetLabelSize(0.07, "Y");
// //   hratio1->SetMaximum(max);
// //cout << max << endl;
//    hratio1->SetMaximum(1.);
//    hratio1->SetMinimum(-0.1);
//    setFilledHist(hratio1 , kBlack, kBlue, 3005);
//    hratio1->Draw();

//    setFilledHist(hratio2 , kBlack, kRed, 3005);
//    hratio2->DrawCopy("histsame");
   
//    fPads[2]->cd(); shrinkPad(0.001, 0.2); 
//    shrinkPad(0.4, 0.2, 0.1, 0.001);

//    gPad->SetGridx(1);  gPad->SetGridy(1);

//    hratioratio1->SetMinimum(0.5); hratioratio1->SetMaximum(1.5);
//    hratioratio1->SetMarkerStyle(24);
//    hratioratio1->SetNdivisions(504, "Y");
//    hratioratio1->SetLabelSize(0.22, "X");  hratioratio1->SetLabelSize(0.17, "Y");
//    hratioratio1->SetStats(0);
//    hratioratio1->SetTitle("");
//    hratioratio1->Draw();

//    fPads[3]->cd(); shrinkPad(0.001, 0.2); 

// //   double max = 1.4*hratio1->GetMaximum();
//    hratio3->SetLabelSize(0.07, "Y");
// //   hratio1->SetMaximum(max);
// //cout << max << endl;
//    hratio3->SetMaximum(1.);
//    hratio3->SetMinimum(-0.1);
//    setFilledHist(hratio3 , kBlack, kBlue, 3005);
//    hratio3->Draw();

//    setFilledHist(hratio4 , kBlack, kRed, 3005);
//    hratio4->DrawCopy("histsame");
   
//    fPads[4]->cd(); shrinkPad(0.001, 0.2); 
//    shrinkPad(0.4, 0.2, 0.1, 0.001);

//    gPad->SetGridx(1);  gPad->SetGridy(1);

//    hratioratio2->SetMinimum(0.5); hratioratio2->SetMaximum(1.5);
//    hratioratio2->SetMarkerStyle(24);
//    hratioratio2->SetNdivisions(504, "Y");
//    hratioratio2->SetLabelSize(0.22, "X");  hratioratio2->SetLabelSize(0.17, "Y");
//    hratioratio2->SetStats(0);
//    hratioratio2->SetTitle("");
//    hratioratio2->Draw();

//    sprintf(name, "%s%s%d%s%s", dir.Data(), "/comparisoneff",multipl,thevar.Data(),".ps");
//    c7->SaveAs(name);  
// }

// ----------------------------------------------------------------------
Double_t thecomparison::smeargauss(double invalue, double mean, double sigma){    

  // smears using a gaussian distribution

  double returnvalue=0.;
  double xsmear = gRandom->Gaus(SHIFTNEUT, SIGMANEUT); 
  returnvalue = invalue + xsmear;
  //cout << sigma << mean << " " << invalue << " " << returnvalue-mean-invalue << endl;
  ((TH1D*)gDirectory->Get("h8888"))->Fill(returnvalue - invalue);
  return returnvalue;
}

void thecomparison::Init(TTree *tree)
{
  //   Set branch addresses
   if (tree == 0) return;
   fChain    = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run",&run);
   fChain->SetBranchAddress("lower",&lower);
   fChain->SetBranchAddress("upper",&upper);
   fChain->SetBranchAddress("nchg",&nchg);
   fChain->SetBranchAddress("nneu",&nneu);
   fChain->SetBranchAddress("nB",&nB);
   fChain->SetBranchAddress("IdB",&IdB);  
   fChain->SetBranchAddress("brecoflav",&brecoflav);
   fChain->SetBranchAddress("brecocharge",&brecocharge);
   fChain->SetBranchAddress("xcharge",&xcharge);
   fChain->SetBranchAddress("mode",&mode);
   fChain->SetBranchAddress("mes",&mes);
   fChain->SetBranchAddress("de",&de);
   fChain->SetBranchAddress("pur",&pur);
   fChain->SetBranchAddress("intpur",&intpur);
   //   fChain->SetBranchAddress("modeB",&modeB);
  fChain->SetBranchAddress("nle",&nle);
  fChain->SetBranchAddress("nel",&nel);
  fChain->SetBranchAddress("nmu",&nmu);
  fChain->SetBranchAddress("nkp",&nkp);
  fChain->SetBranchAddress("nks",&nks);
  fChain->SetBranchAddress("nlept500",&nlept500); 
  fChain->SetBranchAddress("plab",&plab);
  fChain->SetBranchAddress("tlab",&tlab);
  fChain->SetBranchAddress("flab",&flab);
  fChain->SetBranchAddress("elab",&elab);
  fChain->SetBranchAddress("pcms",&pcms);
  fChain->SetBranchAddress("ecms",&ecms);
  fChain->SetBranchAddress("tcms",&tcms);
  fChain->SetBranchAddress("fcms",&fcms);
  fChain->SetBranchAddress("lcharge",&lcharge);
  fChain->SetBranchAddress("isele",&isele);  
  fChain->SetBranchAddress("wdeltam",&wdeltam);
  fChain->SetBranchAddress("Eneualt",&Eneualt);
  //pi l nu
  fChain->SetBranchAddress("indexbestPi",&indexbestPi);
  fChain->SetBranchAddress("chbestPi",&chbestPi);
  fChain->SetBranchAddress("barembestPi",&barembestPi);
  fChain->SetBranchAddress("mm2bestPi",&mm2bestPi);
  fChain->SetBranchAddress("q2bestPi",&q2bestPi); 
  fChain->SetBranchAddress("nrecoPi",&nrecoPi);
  //pi0 l nu
  fChain->SetBranchAddress("indexbestPi0",&indexbestPi0);
  fChain->SetBranchAddress("chbestPi0",&chbestPi0);
  fChain->SetBranchAddress("barembestPi0",&barembestPi0);
  fChain->SetBranchAddress("mm2bestPi0",&mm2bestPi0);
  fChain->SetBranchAddress("q2bestPi0",&q2bestPi0); 
  fChain->SetBranchAddress("nrecoPi0",&nrecoPi0);
  fChain->SetBranchAddress("Estar1dauPi0",&Estar1dauPi0);
  fChain->SetBranchAddress("Estar2dauPi0",&Estar2dauPi0);
  //eta l nu
  fChain->SetBranchAddress("indexbestEta",&indexbestEta);
  fChain->SetBranchAddress("barembestEta",&barembestEta);
  fChain->SetBranchAddress("mm2bestEta",&mm2bestEta);
  fChain->SetBranchAddress("nrecoEta",&nrecoEta);
  fChain->SetBranchAddress("modeEta",&modeEta);
  fChain->SetBranchAddress("Estar1dauEta",&Estar1dauEta);
  fChain->SetBranchAddress("Estar2dauEta",&Estar2dauEta);
  fChain->SetBranchAddress("Estar3dauEta",&Estar3dauEta);
  fChain->SetBranchAddress("Elab1dauEta",&Elab1dauEta);
  fChain->SetBranchAddress("Elab2dauEta",&Elab2dauEta);
  fChain->SetBranchAddress("Elab3dauEta",&Elab3dauEta);
  //etap l nu
  fChain->SetBranchAddress("indexbestEtap",&indexbestEtap);
  fChain->SetBranchAddress("barembestEtap",&barembestEtap);
  fChain->SetBranchAddress("mm2bestEtap",&mm2bestEtap);
  fChain->SetBranchAddress("nrecoEtap",&nrecoEtap);
  fChain->SetBranchAddress("modeEtap",&modeEtap);
  fChain->SetBranchAddress("EtamassdauEtap",&EtamassdauEtap);
  fChain->SetBranchAddress("GammamomdauEtap",&GammamomdauEtap);
  fChain->SetBranchAddress("Rho0massdauEtap",&Rho0massdauEtap);
  fChain->SetBranchAddress("Estar1dauEtap",&Estar1dauEtap);
  fChain->SetBranchAddress("Estar2dauEtap",&Estar2dauEtap);
  fChain->SetBranchAddress("Estar3dauEtap",&Estar3dauEtap);
  fChain->SetBranchAddress("Elab1dauEtap",&Elab1dauEtap);
  fChain->SetBranchAddress("Elab2dauEtap",&Elab2dauEtap);
  fChain->SetBranchAddress("Elab3dauEtap",&Elab3dauEtap);
  //MC truth
//   fChain->SetBranchAddress("isassocB",&isassocB);
//   fChain->SetBranchAddress("ass_deltapB",&ass_deltapB);
//   fChain->SetBranchAddress("vub",&vub);
//   fChain->SetBranchAddress("vcb",&vcb);
//   fChain->SetBranchAddress("other",&other);
//   fChain->SetBranchAddress("nvubexcl",&nvubexcl);
//   fChain->SetBranchAddress("mxhadgen",&mxhadgen);
//   fChain->SetBranchAddress("pcmsgen",&pcmsgen);
//   fChain->SetBranchAddress("ecmsgen",&ecmsgen);
//   fChain->SetBranchAddress("pxhadgen",&pxhadgen);
//   fChain->SetBranchAddress("exhadgen",&exhadgen);
//   fChain->SetBranchAddress("q2Gen",&q2Gen);
//   fChain->SetBranchAddress("Gvxbtyp",&Gvxbtyp);
//   fChain->SetBranchAddress("GSem",&GSem);
//   fChain->SetBranchAddress("GfDpi",&GfDpi);
//   fChain->SetBranchAddress("GfDpiz",&GfDpiz);
//   fChain->SetBranchAddress("GfDk",&GfDk);
//   fChain->SetBranchAddress("GfDks",&GfDks);
//   fChain->SetBranchAddress("GfDkl",&GfDkl);
//   fChain->SetBranchAddress("GfDlep",&GfDlep);
//   fChain->SetBranchAddress("GfDgam",&GfDgam);
//   fChain->SetBranchAddress("GfDnu",&GfDnu);
//   fChain->SetBranchAddress("GfD0Ds",&GfD0Ds);
//   fChain->SetBranchAddress("GfDDs",&GfDDs);
//   fChain->SetBranchAddress("GfDkspiopio",&GfDkspiopio);

   Notify();
}

Bool_t thecomparison::Notify()
{
//   called when loading a new file
//   get branch pointers
   b_run = fChain->GetBranch("run");
   b_lower = fChain->GetBranch("lower");
   b_upper = fChain->GetBranch("upper");
   b_nchg = fChain->GetBranch("nchg");
   b_nneu = fChain->GetBranch("nneu");
   b_nB = fChain->GetBranch("nB");
   b_IdB = fChain->GetBranch("IdB");
   b_brecoflav = fChain->GetBranch("brecoflav");
   b_brecocharge = fChain->GetBranch("brecocharge");
   b_xcharge = fChain->GetBranch("xcharge");
   b_mode = fChain->GetBranch("mode");
   b_mes = fChain->GetBranch("mes");
   b_de = fChain->GetBranch("de");
   b_pur = fChain->GetBranch("pur");
   b_intpur = fChain->GetBranch("intpur");
   b_nle = fChain->GetBranch("nle");
   b_nel = fChain->GetBranch("nel");
   b_nmu = fChain->GetBranch("nmu");
   b_nkp = fChain->GetBranch("nkp");
   b_nks = fChain->GetBranch("nks");
   b_nlept500 = fChain->GetBranch("nlept500");
   b_plab = fChain->GetBranch("plab");
   b_tlab = fChain->GetBranch("tlab");
   b_flab = fChain->GetBranch("flab");
   b_elab = fChain->GetBranch("elab");
   b_pcms = fChain->GetBranch("pcms");
   b_ecms = fChain->GetBranch("ecms");
   b_tcms = fChain->GetBranch("tcms");
   b_fcms = fChain->GetBranch("fcms");
   b_lcharge = fChain->GetBranch("lcharge");
   b_isele = fChain->GetBranch("isele");
   b_wdeltam = fChain->GetBranch("wdealtam");
   b_Eneualt = fChain->GetBranch("Eneualt");

   //pi l nu
   b_indexbestPi = fChain->GetBranch("indexbestPi");
   b_chbestPi = fChain->GetBranch("chbestPi");
   b_barembestPi = fChain->GetBranch("barembestPi");
   b_mm2bestPi = fChain->GetBranch("mm2bestPi");
   b_nrecoPi = fChain->GetBranch("nrecoPi");

   //pi0 l nu
   b_indexbestPi0 = fChain->GetBranch("indexbestPi0");
   b_chbestPi0 = fChain->GetBranch("chbestPi0");
   b_barembestPi0 = fChain->GetBranch("barembestPi0");
   b_mm2bestPi0 = fChain->GetBranch("mm2bestPi0");
   b_q2bestPi0 = fChain->GetBranch("q2bestPi0");
   b_nrecoPi0 = fChain->GetBranch("nrecoPi0");
   b_ndauPi0 = fChain->GetBranch("ndauPi0");
   b_Estar1dauPi0 = fChain->GetBranch("Estar1dauPi0");
   b_Estar2dauPi0 = fChain->GetBranch("Estar2dauPi0");

   //eta l nu 
   b_indexbestEta = fChain->GetBranch("indexbestEta");
   b_barembestEta = fChain->GetBranch("barembestEta");
   b_mm2bestEta = fChain->GetBranch("mm2bestEta");
   b_nrecoEta = fChain->GetBranch("nrecoEta");
   b_ndauEta = fChain->GetBranch("ndauEta");
   b_modeEta = fChain->GetBranch("modeEta");
   b_Estar1dauEta = fChain->GetBranch("Estar1dauEta");
   b_Estar2dauEta = fChain->GetBranch("Estar2dauEta");
   b_Estar3dauEta = fChain->GetBranch("Estar3dauEta");
   b_Elab1dauEta = fChain->GetBranch("Elab1dauEta");
   b_Elab2dauEta = fChain->GetBranch("Elab2dauEta");
   b_Elab3dauEta = fChain->GetBranch("Elab3dauEta");

   //etap l nu
   b_indexbestEtap = fChain->GetBranch("indexbestEtap");
   b_barembestEtap = fChain->GetBranch("barembestEtap");
   b_mm2bestEtap = fChain->GetBranch("mm2bestEtap");
   b_nrecoEtap = fChain->GetBranch("nrecoEtap");
   b_ndauEtap = fChain->GetBranch("ndauEtap");
   b_modeEtap = fChain->GetBranch("modeEtap");
   b_EtamassdauEtap = fChain->GetBranch("EtamassdauEtap");
   b_Rho0massdauEtap = fChain->GetBranch("Rho0massdauEtap");
   b_GammamomdauEtap = fChain->GetBranch("GammamomdauEtap");
   b_Estar1dauEtap = fChain->GetBranch("Estar1dauEtap");
   b_Estar2dauEtap = fChain->GetBranch("Estar2dauEtap");
   b_Estar3dauEtap = fChain->GetBranch("Estar3dauEtap");

   //MC truth
//    b_isassocB = fChain->GetBranch("isassocB");
//    b_ass_deltapB = fChain->GetBranch("ass_deltapB");
//    b_vub = fChain->GetBranch("vub");
//    b_vcb = fChain->GetBranch("vcb");
//    b_other = fChain->GetBranch("other");
//    b_nvubexcl = fChain->GetBranch("nvubexcl");
//    b_mxhadgen = fChain->GetBranch("mxhadgen");
//    b_pcmsgen = fChain->GetBranch("pcmsgen");
//    b_ecmsgen = fChain->GetBranch("ecmsgen");
//    b_pxhadgen = fChain->GetBranch("pxhadgen");
//    b_exhadgen = fChain->GetBranch("exhadgen");
//    b_q2Gen = fChain->GetBranch("q2Gen");
//    b_Gvxbtyp = fChain->GetBranch("Gvxbtyp");
//    b_GSem = fChain->GetBranch("GSem");
//    b_GfDpi = fChain->GetBranch("GfDpi");
//    b_GfDpiz = fChain->GetBranch("GfDpiz");
//    b_GfDk = fChain->GetBranch("GfDk");
//    b_GfDks = fChain->GetBranch("GfDks");
//    b_GfDkl = fChain->GetBranch("GfDkl");
//    b_GfDlep = fChain->GetBranch("GfDlep");
//    b_GfDgam = fChain->GetBranch("GfDgam");
//    b_GfDnu = fChain->GetBranch("GfDnu");
//    b_GfD0Ds = fChain->GetBranch("GfD0Ds");
//    b_GfDDs = fChain->GetBranch("GfDDs");
//    b_GfDkspiopio = fChain->GetBranch("GfDkspiopio");
   return kTRUE;
}
// ----------------------------------------------------------------------
void thecomparison::readCuts(TString filename) {
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
  ELECTRONPCUT = 1.;
  MUONPCUT = 1.;
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
  MXCUTLOWEXCL1=-1000.;
  MXCUTHIGHEXCL1=1000.;
  MXCUTLOWEXCL2=-1000.;
  MXCUTHIGHEXCL2=1000.;
  MXCUTLOWEXCL3=-1000.;    
  MXCUTHIGHEXCL3=1000.;
  MXCUTLOWEXCL4=-1000.;
  MXCUTHIGHEXCL4=1000.;
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
  DAUETAMASS2=1000.;
  DAUETAMASS3=1000.;
  DAUETAMASS4=1000.;
  DAURHOMASS=1000.;
  DAUGAMMAMOM = -1000.;
  DORIGHTNCHG = 0;
  MNUSQPI0LOW =  0.;          
  MNUSQPI0HIGH = 0.;
  MNUSQETALOW = 0.;
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
    if (!strcmp(CutName, "electronPCut")) { ELECTRONPCUT = CutValue; ok = 1;}
    if (!strcmp(CutName, "muonPCut")) { MUONPCUT = CutValue; ok = 1;}    
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
    if (!strcmp(CutName, "mxCutHighExcl")) { MXCUTHIGHEXCL  = CutValue; ok = 1;}    if (!strcmp(CutName, "mxCutLowExcl1")) { MXCUTLOWEXCL1  = CutValue; ok = 1;}    if (!strcmp(CutName, "mxCutHighExcl1")) { MXCUTHIGHEXCL1  = CutValue; ok = 1;}
    if (!strcmp(CutName, "mxCutLowExcl2")) { MXCUTLOWEXCL2  = CutValue; ok = 1;}    if (!strcmp(CutName, "mxCutHighExcl2")) { MXCUTHIGHEXCL2  = CutValue; ok = 1;}
    if (!strcmp(CutName, "mxCutLowExcl3")) { MXCUTLOWEXCL3  = CutValue; ok = 1;}    if (!strcmp(CutName, "mxCutHighExcl3")) { MXCUTHIGHEXCL3  = CutValue; ok = 1;}
    if (!strcmp(CutName, "mxCutLowExcl4")) { MXCUTLOWEXCL4  = CutValue; ok = 1;}    if (!strcmp(CutName, "mxCutHighExcl4")) { MXCUTHIGHEXCL4  = CutValue; ok = 1;}
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
    if (!strcmp(CutName, "dauetamass2")) { DAUETAMASS2 = CutValue; ok = 1;}
    if (!strcmp(CutName, "dauetamass3")) { DAUETAMASS3 = CutValue; ok = 1;}
    if (!strcmp(CutName, "dauetamass4")) { DAUETAMASS4 = CutValue; ok = 1;}
    if (!strcmp(CutName, "daurhomass")) { DAURHOMASS = CutValue; ok = 1;}
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
 
    if (ok == 0)  cout << "==> readCuts() Error: Don't know about variable " << CutName << endl;
  }
  dumpCuts();
 
}                                                                                                 
// ----------------------------------------------------------------------
void thecomparison::dumpCuts() {
  cout << "====================================" << endl;
  cout << " fitcategory         :"  << FITCATEGORY   <<  endl;
  cout << " q2Cut               :"  << Q2CUT   <<  endl;
  cout << " leptonPCut          :"  << LEPTONPCUT   <<  endl;
  cout << " electronPCut        :"  << ELECTRONPCUT   <<  endl;
  cout << " muonPCut            :"  << MUONPCUT   <<  endl;
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
  cout << " mxCutLowExcl1        :"  << MXCUTLOWEXCL1   <<  endl;
  cout << " mxCutHighExcl1       :"  << MXCUTHIGHEXCL1   <<  endl;
  cout << " mxCutLowExcl2        :"  << MXCUTLOWEXCL2   <<  endl;
  cout << " mxCutHighExcl2       :"  << MXCUTHIGHEXCL2   <<  endl;
  cout << " mxCutLowExcl3        :"  << MXCUTLOWEXCL3   <<  endl;
  cout << " mxCutHighExcl3       :"  << MXCUTHIGHEXCL3   <<  endl;
  cout << " mxCutLowExcl4        :"  << MXCUTLOWEXCL4   <<  endl;
  cout << " mxCutHighExcl4       :"  << MXCUTHIGHEXCL4   <<  endl;
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
  cout << " dauetamass2          :"  << DAUETAMASS2   <<  endl;
  cout << " dauetamass3          :"  << DAUETAMASS3   <<  endl;   
  cout << " dauetamass4          :"  << DAUETAMASS4   <<  endl;
  cout << " daurhomass          :"  << DAURHOMASS   <<  endl;
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



// ----------------------------------------------------------------------
TChain *  thecomparison::getchain(char *thechain, char* treename) {

  TChain *chain = new TChain(treename);
  cout << "Chaining ... " << thechain << endl;
  char pName[2000];
  char buffer[200];
  sprintf(buffer, "%s", thechain);
  ifstream is(buffer);
  cout << "files " << buffer <<  endl;
  while(is.getline(buffer, 200, '\n')){
    //    if (buffer[0] == '#') continue;
    sscanf(buffer, "%s", pName);
    cout << "   Add: " << buffer << endl;
    chain->Add(pName);
  }
  is.close();
  return chain;

}
thecomparison::~thecomparison()
{
  delete Dvar;
  delete Bsem;
  if (!fChain) return;
   //   delete fChain->GetCurrentFile();
}

Int_t thecomparison::GetEntry(Int_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Int_t thecomparison::LoadTree(Int_t entry)
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

void thecomparison::Show(Int_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t thecomparison::Cut(Int_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
