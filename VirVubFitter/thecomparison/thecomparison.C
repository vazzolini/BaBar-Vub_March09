#include "thecomparison.h"

#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include <TVector2.h>
#include "util.hh"
#include "TH1D.h"

#include <iostream>
#include <fstream>

#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooAddPdf.hh"
#include "RooFitCore/RooExtendPdf.hh"
#include "RooFitCore/RooPlot.hh"

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
   themin = mi;
   themax = ma;
   thebins = b;
   thevar = TString(var) ;
   Dvar = new recoilDSys("ddecay.table.CM2",0,2);
   Bsem = new recoilDSys(0);
}
 
void thecomparison::Loop(int nevents, int cat, double shift, double smear, int bsel, int multcat, int seed,int sys)
{
  char name[100];
  char le[100];
  int id;
  int idflav;
  int ich, ine;
  //double cutflav = 0;    
  double thecat =0;
  //int number;
  int group;        
  
  multipl = multcat;

  SHIFTNEUT = shift;
  SIGMANEUT = smear;
  // restore B0/B+ flag
  int isbch =  bsel>=0? bsel: -(1+bsel);
      
  if (fChain == 0) return;
  
  Int_t nentries = Int_t(fChain->GetEntries());
  
  if( nentries > nevents) nentries = nevents;
  
  cout <<"Nentries = "<< nentries << endl;
  
  Int_t nbytes = 0, nb = 0;
  for (Int_t jentry=0; jentry<nentries;jentry++) {
        //if(cat==5)cout <<" entry # "<<jentry;
    Int_t ientry = LoadTree(jentry); //in case of a TChain, ientry is the entry number in the current file
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    //    if(cat==5)cout <<" loaded "<<endl;
    // vub - vcb selection
    if(bsel<0)
      if( (!vub && cat==4) || (!vcb && cat==5) )continue;

    //    if(mxhadfit > 1.55) continue;
    //    if(!vcb && cat==5) continue;

    // if (Cut(ientry) < 0) continue;
    int flav =  lcharge + brecoflav; // charge correlation
    bool ksele = nkp + nks;          // fit on the depleted sample?
    //bool ksele = nkp>1;          // fit on the depleted sample?
    int ch = xcharge + brecocharge;  // total charge
    
    double cutflav = 0;              // flavor category (3 = charged B, 4 = neutral B OS, 5 = neutral B SS)  
    // int ibch = 1;                     
    
    if(nchg-nle<1) {ich = 1;}
    else if(nchg-nle<3){ich = 2;}
    else{ich = 3;};
    
    if(nneu80_160<1) {ine = 1;}
    else{ine = 2;};      

    if (ich==1 && ine ==1) thecat = 1;
    if (ich==2 && ine ==1) thecat = 2;
    if (ich==3 && ine ==1) thecat = 3;
    if (ich==1 && ine ==2) thecat = 4;
    if (ich==2 && ine ==2) thecat = 5;
    if (ich==3 && ine ==2) thecat = 6;

    // cuts
    
    bool PCMS =  pcms>1.;
    // bool NPI0 = nnpi0<2;

    bool RUN = 1;
       //bool RUN = eneu<0.15*exhad;
    bool NLE = nle > 0;
    //bool NEL = nel;
    //bool NMU = nmu;
    bool NLE1 = (nle == 1 &&pcms>1.) ;
    bool MM2 =  mm2 < .5;
    bool CH = ch == 0;
    bool KSELE = ksele >0;
    int FLAV = !(TMath::Abs(brecocharge)!=0 && (lcharge + brecoflav)!=0);
    bool BCH = TMath::Abs(brecocharge) == isbch;     
    bool IPUR = intpur>0;
    bool CAT = multcat == thecat;

    bool SEED = 1;
    if (seed == 2) SEED = (modeB>13000 && modeB<14000);
    if (seed == 3) SEED = (modeB>14000 && modeB<16000);
    if (seed == 4) SEED = (modeB>12000 && modeB<13000);
    if (seed == 5) SEED = (modeB>11000 && modeB<12000);    

    if (isbch == 2) BCH = 1;
    if (multcat == 7) CAT = 1;

    double myvar = -1;
    if (thevar == "mxhadfit") {myvar = mxhadfit;}
    /*if (thevar == "csiCiuc") {myvar = csiCiuc;}
    if (thevar == "xCiuc") {myvar = xCiuc;}
    if (thevar == "wCiuc") {myvar = wCiuc;}
    if (thevar == "EwPwfit") {myvar = EwPwfit;}
    if (thevar == "GcsiCiuc") {myvar = GcsiCiuc;}
    if (thevar == "GxCiuc") {myvar = GxCiuc;}
    if (thevar == "GwCiuc") {myvar = GwCiuc;}
    if (thevar == "EwPwG") {myvar = EwPwG;}*/
    if (thevar == "q2fit") {myvar = q2fit;}
    if (thevar == "q2") {myvar = q2;}
    if (thevar == "q2Gen") {myvar = q2Gen;}
    //if (thevar == "rescsiCiuc") {myvar = csiCiuc-GcsiCiuc;}
    //if (thevar == "resxCiuc")   {myvar = xCiuc-GxCiuc;}
    //if (thevar == "reswCiuc")   {myvar = wCiuc-GwCiuc;}
    //if (thevar == "resEwPwfit") {myvar = EwPwfit-EwPwG;}
    if (thevar == "resq2fit")   {myvar = q2fit-q2Gen;}
    if (thevar == "mxhad") {myvar = mxhad;}
    if (thevar == "mm2") {myvar = mm2; MM2=1; }
    if (thevar == "pplus") myvar = pplus; 
    if (thevar == "wdeltam") myvar = wdeltam; 
    if (thevar == "ksele") myvar = ksele;
    //if (thevar == "efneu") {myvar =exhad>0 ? eneu/exhad:-1;}
    //if (thevar == "eneu") {myvar =eneu;}
    //if (thevar == "epi0") {myvar =epiz;}
    //if (thevar == "esneu") {myvar =nneu>0 ? eneu/nneu:-1;}
    //if (thevar == "etrk") {myvar = exhad-eneu;}
    //if (thevar == "kmin") {if(nkp<=0)continue;myvar = kminmom;}
    //if (thevar == "kmax") {if(nkp<=0)continue;myvar = kmaxmom;}
    if (thevar == "nneu") myvar = nneu;
    if (thevar == "pnu") myvar = pnu;
    /*if (thevar == "dxy")  {myvar = sqrt(dx*dx+dy*dy);if(myvar>1.)continue;}
    if (thevar == "sxy")  {
      myvar = dx*dx+dy*dy;
      if(myvar>1. || myvar<= 0)continue;
      myvar=myvar/sqrt(s2dxx*dx*dx+s2dyy*dy*dy+2*s2dxy*dx*dy);
    }
    if (thevar == "d3d") {myvar = sqrt(dx*dx+dy*dy+dz*dz);if(myvar>1.)continue;}
    if (thevar == "s3d")  {
      myvar = dx*dx+dy*dy+dz*dz;
      if(myvar>1. || myvar<= 0)continue;
      myvar=myvar/sqrt(s2dxx*dx*dx+s2dyy*dy*dy+2*s2dxy*dx*dy+s2dzz*dz*dz+2*s2dxz*dx*dz+2*s2dyz*dy*dz);
      }*/
    if (thevar == "nneu80_160") myvar = nneu80_160;
    if (thevar == "nneu160_320") myvar = nneu160_320;
    if (thevar == "nneu320") myvar = nneu - nneu80_160;
    if (thevar == "nchg") myvar = nchg;
    if (thevar == "pcms") {myvar = pcms; PCMS = 1;}
    if (thevar == "nle") {myvar = nle; NLE = 1; NLE1 = 1;}
    if (thevar == "qtot") {myvar = ch; CH = 1;}
    //    if (thevar == "npi0") {myvar = npi0; }
    if (thevar == "nkp") {myvar = nkp; KSELE = nks > 0; }
    if (thevar == "nks") {myvar = nks; KSELE = nkp > 0; }
    if (thevar == "intpur") {myvar = intpur; }
    if (thevar == "pur") {myvar = pur; }
    if (thevar == "modeB") {myvar = modeB; }
    if (!((myvar>0) || (myvar<0) || (myvar == 0))) {
      myvar = -999.;
    } 
    totweight = 1;
    
    if(cat!=4 && sys==1){
      int dImode=0;
      //      cout<<" "<<Gvxbtyp<<" "<<vub<<" "<<GfDpi<<" "<<GfDk<<" "<<endl;
      totweight *= getBsysweight(Gvxbtyp,vub);//Bdec weighting
      totweight *= getDsysweight(GfDpi,GfDk,GfDks,GfDpiz,GfDlep,dImode,vub); //Ddec weighting

    } 
    // COMMENTED OUT ON 24-MAY-2006
    //else if( cat!=4 && sys==2&&vub&&fkplus>0){
    //  totweight *= (1.5*FermiWeight(fkplus,-0.15,0.));// 1.5 to account for the possible difference in yields
      
    // }

    //totweight = 1; //ATTENZIONE AGGIUNTO PER IL COMP DATI DATI e MC MC

    if(PCMS && NLE && FLAV && BCH && IPUR && CAT && SEED && RUN) {
      sprintf(le, "h");
      if(KSELE) sprintf(le, "d");
      if(SIGMANEUT>0.0) myvar = smeargauss(myvar, SHIFTNEUT, SIGMANEUT); 
      // assign flavor category
      cutflav = 5;
      if(TMath::Abs(brecocharge)) cutflav = 3;
      if(TMath::Abs(brecocharge)==0 && flav==0) cutflav = 4;	  
      id = hist(myvar);
      idflav =(int)(hist(myvar) + cutflav * 10000);    
      group = 200 + cat * 1000 + 100000;
      sprintf(name, "%s%d",le,group+idflav);
      ((TH1D*)gDirectory->Get(name))->Fill(mes,totweight);  
      if(MM2 && CH && NLE1){  
	group = 200 + cat * 1000;
	sprintf(name, "%s%d",le,group+idflav);
	((TH1D*)gDirectory->Get(name))->Fill(mes,totweight);  
      }
      //       group = 90206 + cat * 1000;
      //        sprintf(name, "%s%d",le,group);
      //        ((TH1D*)gDirectory->Get(name))->Fill(mes,totweight);  
    }
    
    if(PCMS && NLE && !flav && BCH && IPUR && CAT && SEED &&RUN){  
      
      sprintf(le, "h");
      if(KSELE) sprintf(le, "d");
      group = 90206 + cat * 1000 + 100000;
      sprintf(name, "%s%d",le,group);
      ((TH1D*)gDirectory->Get(name))->Fill(mes,totweight);  
      if(MM2 && CH && NLE1){  
	group = 90206 + cat * 1000;
	sprintf(name, "%s%d",le,group);
	((TH1D*)gDirectory->Get(name))->Fill(mes,totweight);  
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
 
  sprintf(name, "h400000");  sprintf(title, "%s%s", thevar.Data(), " data events after all cuts: enriched");  h = new TH1D(name, title, thebins, themin, themax );  h->Sumw2();
  sprintf(name, "h500000");  sprintf(title, "%s%s", thevar.Data(), " MC events  after all cuts: enriched");  h = new TH1D(name, title, thebins, themin, themax  );  h->Sumw2();
  sprintf(name, "h94206");  sprintf(title, "mes data after all cuts: enriched");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  sprintf(name, "h95206");  sprintf(title, "mes MC after all cuts: enriched");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  sprintf(name, "h1400000");  sprintf(title, "%s%s", thevar.Data(), " data events after lepton cuts: enriched");  h = new TH1D(name, title, thebins, themin, themax );  h->Sumw2();
  sprintf(name, "h1500000");  sprintf(title, "%s%s", thevar.Data(), " MC events after lepton cuts: enriched");  h = new TH1D(name, title, thebins, themin, themax  );  h->Sumw2();
  sprintf(name, "h194206");  sprintf(title, "mes data after lepton cuts: enriched");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  sprintf(name, "h195206");  sprintf(title, "mes MC after lepton cuts: enriched");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();

  sprintf(name, "d400000");  sprintf(title, "%s%s", thevar.Data(), " data events after all cuts: depleted");  h = new TH1D(name, title, thebins, themin, themax );  h->Sumw2();
  sprintf(name, "d500000");  sprintf(title, "%s%s", thevar.Data(), " MC events  after all cuts: depleted");  h = new TH1D(name, title, thebins, themin, themax  );  h->Sumw2();
  sprintf(name, "d94206");  sprintf(title, "mes data after all cuts: depleted");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  sprintf(name, "d95206");  sprintf(title, "mes MC after all cuts: depleted");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  sprintf(name, "d1400000");  sprintf(title, "%s%s", thevar.Data(), " data events after lepton cuts: depleted");  h = new TH1D(name, title, thebins, themin, themax );  h->Sumw2();
  sprintf(name, "d1500000");  sprintf(title, "%s%s", thevar.Data(), " MC events after lepton cuts: depleted");  h = new TH1D(name, title, thebins, themin, themax  );  h->Sumw2();
  sprintf(name, "d194206");  sprintf(title, "mes data after lepton cuts: depleted");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  sprintf(name, "d195206");  sprintf(title, "mes MC after lepton cuts: depleted");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();

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

    sprintf(name,"%s%s" , "d4",number);  sprintf(title, "mes data");  h = new TH1D(name, title, 40, 5.2, 5.3);         h->Sumw2();
    sprintf(name,"%s%s" , "d34",number);  sprintf(title, "mes data bch");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "d44",number);  sprintf(title, "mes data b0 os");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "d54",number);  sprintf(title, "mes data b0 ss");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "d5",number);  sprintf(title, "mes MC");  h = new TH1D(name, title, 40, 5.2, 5.3);         h->Sumw2();
    sprintf(name,"%s%s" , "d35",number);  sprintf(title, "mes MC bch");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "d45",number);  sprintf(title, "mes MC b0 os");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "d55",number);  sprintf(title, "mes MC b0 ss");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();

    sprintf(name,"%s%s" , "d14",number);  sprintf(title, "mes data");  h = new TH1D(name, title, 40, 5.2, 5.3);         h->Sumw2();
    sprintf(name,"%s%s" , "d134",number);  sprintf(title, "mes data bch");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "d144",number);  sprintf(title, "mes data b0 os");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "d154",number);  sprintf(title, "mes data b0 ss");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "d15",number);  sprintf(title, "mes MC");  h = new TH1D(name, title, 40, 5.2, 5.3);         h->Sumw2();
    sprintf(name,"%s%s" , "d135",number);  sprintf(title, "mes MC bch");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "d145",number);  sprintf(title, "mes MC b0 os");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "d155",number);  sprintf(title, "mes MC b0 ss");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  }
   sprintf(name, "h8888");  sprintf(title, "test gauss");  h = new TH1D(name, title, 100, -.05, .05);    h->Sumw2(); 
}

int thecomparison::hist(double mx){

  // categories
  int bin = int((mx-themin)/((themax-themin)/thebins)+1);
  if (mx>themax || mx==themax) bin = thebins;
if (mx<themin || mx==themin) bin = 1;
  return bin;
}

void thecomparison::Fitmes(int cat, int cut, bool istmodel){

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
  
  vector<double> results(14);
  vector<double> inputPar(14);
  vector<double> fooresults(14);

  int i;
  int addcut = 0; 
  if(cut) addcut = 1;
  
  bool isData = cat == 5 ? 1 : 0;

  int const theb = thebins + 1;
  for (int y=0; y<2; y++) {
    sprintf(le, "h");
    if (y == 1) sprintf(le, "d");
    
    //extracting the signal fit parameters 
    sprintf(name, "%s%d",le,group+90006+addcut*100000);  
    cout << ((TH1D*)gDirectory->Get(name))->Integral() << endl;
  
    // this fit is to get good starting values, do a fit with not all parameters fixed for tmodel
    double dummy1,dummy2;
    if(istmodel){
      
      // put in starting values, take them from ../mesparsetting_thorsten.dat
      inputPar[thecomparison::iThoSigR] =0.936655;
      inputPar[thecomparison::iSigma_r1] =0.00158;
      inputPar[thecomparison::iThoSigXc] =5.27975;
      inputPar[thecomparison::iSigma_r2] =0.00226044;
      inputPar[thecomparison::iSigma_l] =0.00192;
      inputPar[thecomparison::iThoSigN] =1.24962;
      inputPar[thecomparison::iThoSigAlpha] =4.19962;
					  
      inputPar[thecomparison::iMean] = 5.27992;
      inputPar[thecomparison::iSigma] = 0.00510;
      inputPar[thecomparison::iN] = 20.;
      inputPar[thecomparison::iAlpha]= 0.46452;
      
      if(isData)
	inputPar[thecomparison::iEndpoint] =5.28938;
      else
	inputPar[thecomparison::iEndpoint] =5.2899;

      inputPar[thecomparison::iArgus] =-60.;
      inputPar[thecomparison::iCutOff]=5.28927;

      sighisto_newmodel(dummy1,dummy2,(TH1D*)gDirectory->Get(name),results,0,inputPar,isData);
      cout << "mes result for: " << name << " is, xc " << results[thecomparison::iThoSigXc] << " SIGMA_L " << results[thecomparison::iSigma_l] << " ALPHA "
	   << results[thecomparison::iThoSigAlpha] << " N " << results[thecomparison::iThoSigN]
	   << "R "<< results[thecomparison::iThoSigR]<<" SIGMA_R1 "<<results[thecomparison::iSigma_r1]<<" SIGMA_R2 "<<results[thecomparison::iSigma_r2]<<endl;
    } else {
      sighisto(dummy1,dummy2,(TH1D*)gDirectory->Get(name),usedmean,usedsigma,usedalpha,usedn,1,-11111111.,-1111111.,-1111111.,-1111111.,-1111111.);
      cout << "mes result for: " << name << " is, MEAN " << usedmean << " SIGMA " << usedsigma << " ALPHA " << usedalpha << " N " << usedn << endl;
    }
    
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
    
      // put the result of the fit in each Mx bin (fix all parameters, except for the argus shape parameter)
    
      if(istmodel)
	sighisto_newmodel(sigs[is-1],errsigs[is-1],(TH1D*)gDirectory->Get(name),fooresults,1,results,isData);
      else
	sighisto(sigs[is-1],errsigs[is-1],(TH1D*)gDirectory->Get(name),foomean,foosigma,fooalpha,foon,1,usedmean,usedsigma,usedalpha,usedn,-1111111);
      
      int title = cat * 100000 + addcut * 1000000;
      sprintf(name, "%s%d",le,title);  
      ((TH1D*)gDirectory->Get(name))->SetBinContent(is, sigs[is-1]);
      ((TH1D*)gDirectory->Get(name))->SetBinError(is, errsigs[is-1]);
    }
  }
  
}

mesData* thecomparison::vubMesUnb(RooDataHist *data, RooRealVar *x, vector<double>& results, int fixpar, const vector<double>& inputPar, bool isData)
{
  RooRealVar* mes(x);
  
  //define mesint
  mes->setRange("mesint",5.27,5.29);
  cout << "dataset title = " << data->GetTitle() << std::endl;

  //Signal and Background events
  RooRealVar* nsig = new RooRealVar("S","number of sig events", 100., 0., 800000.);
  RooRealVar* nbpk = new RooRealVar("P","number of bpk events", 1., 0., 800000.);
  RooRealVar* nbkg = new RooRealVar("B","number of bkg events", 1., 0., 800000.);

  RooAbsPdf* model = 0;

  RooAbsPdf* pArgus = createArgus(*mes);

  /* --- Build modified Crystal Ball peaking background PDF --- */
  RooAbsPdf* pCcb = createCCB(*mes);

  /* --- Build modified Gauz signal PDF --- */
  RooAbsPdf* pSignal = createThorsten(*mes);

  // build extended pdfs
  RooExtendPdf* ae = new RooExtendPdf("ae", "argus extended", *pArgus, *nbkg, "mesint");
  RooExtendPdf* ccbe = new RooExtendPdf("ccbe", "ccb extended", *pCcb, *nbpk, "mesint");
  RooExtendPdf* se = new RooExtendPdf("se","signal extended", *pSignal, *nsig, "mesint");
    
  // build final model
  model = new RooAddPdf("model","a+ccb+s", RooArgList(*ae,*ccbe,*se));
  

  Double_t histoentries=data->numEntries(kTRUE); //To Get the right number of entries
  cout<<"======================> "<<histoentries<<endl;
  // get starting values for nsig, nbkg, nbpk
  if(isData){
    //scaling with Entries in Dataset taken from web page 11 May
    nsig->setVal(histoentries*90138/176470);
    nbpk->setVal(histoentries*30650/176470);
    nbkg->setVal(histoentries*56698/176470);
  } else if(!isData){
    //scaling with Entries in Dataset taken from web page 
    nsig->setVal(histoentries*558553/1563108);
    nbpk->setVal(histoentries*191111/1563108);
    nbkg->setVal(histoentries*813444/1563108);
  }


  if (fixpar == 0) { // preset variables, but don't fix all
       
    if (inputPar[thecomparison::iThoSigR] > 0.)     { Thor->setVal(inputPar[thecomparison::iThoSigR]);         Thor->setConstant(); }
    if (inputPar[thecomparison::iSigma_r1] > 0.)    { Thosigma_r1->setVal(inputPar[thecomparison::iSigma_r1]); }
    if (inputPar[thecomparison::iThoSigXc] > 0.)    { Thoxc->setVal(inputPar[thecomparison::iThoSigXc]); }          
    if (inputPar[thecomparison::iSigma_r2] > 0.)    { Thosigma_r2->setVal(inputPar[thecomparison::iSigma_r2]); Thosigma_r2->setConstant(); }
    if (inputPar[thecomparison::iSigma_l]  > 0.)    { Thosigma_l->setVal(inputPar[thecomparison::iSigma_l]); }
    if (inputPar[thecomparison::iThoSigN] > 0.)     { Thon->setVal(inputPar[thecomparison::iThoSigN]);         Thon->setConstant(); }
    if (inputPar[thecomparison::iThoSigAlpha] > 0.) { Thoalpha->setVal(inputPar[thecomparison::iThoSigAlpha]); Thoalpha->setConstant(); }
    
    if (inputPar[thecomparison::iMean]  > 0.)       { Rm->setVal(inputPar[thecomparison::iMean]);  Rm->setConstant(); }
    if (inputPar[thecomparison::iSigma] > 0.)       { Rs->setVal(inputPar[thecomparison::iSigma]); Rs->setConstant(); }
    if (inputPar[thecomparison::iN] > 0.)           { Rn->setVal(inputPar[thecomparison::iN]);     Rn->setConstant(); }
    if (inputPar[thecomparison::iAlpha] > 0.)       { Ra->setVal(inputPar[thecomparison::iAlpha]); Ra->setConstant(); }
    if (isData) { // free on data
      if (inputPar[thecomparison::iEndpoint] > 0.)  { Rendpoint->setVal(inputPar[thecomparison::iEndpoint]); }
    } else {
      if (inputPar[thecomparison::iEndpoint] > 0.)  { Rendpoint->setVal(5.2891); Rendpoint->setConstant(); }
    }
  
    if (inputPar[thecomparison::iArgus] > -55.)     { pArgPar->setVal(inputPar[thecomparison::iArgus]); }
    if (isData) { // free on data
      if (inputPar[thecomparison::iCutOff] > 0.)    { pCutOff->setVal(inputPar[thecomparison::iCutOff]); }
    } else {      // fixed on MC at 5.2891
      if (inputPar[thecomparison::iCutOff] > 0.)    { pCutOff->setVal(5.2891); pCutOff->setConstant(); pCutOff->setConstant(); }
    }
  
  } else if(fixpar == 1) { // set and fix all parameters (if reasonable value is given)
    
    if (inputPar[thecomparison::iThoSigR] > 0.)     { Thor->setVal(inputPar[thecomparison::iThoSigR]);         Thor->setConstant(); }
    if (inputPar[thecomparison::iSigma_r1] > 0.)    { Thosigma_r1->setVal(inputPar[thecomparison::iSigma_r1]); Thosigma_r1->setConstant(); }
    if (inputPar[thecomparison::iThoSigXc] > 0.)    { Thoxc->setVal(inputPar[thecomparison::iThoSigXc]);       Thoxc->setConstant(); }
    if (inputPar[thecomparison::iSigma_r2] > 0.)    { Thosigma_r2->setVal(inputPar[thecomparison::iSigma_r2]); Thosigma_r2->setConstant(); }
    if (inputPar[thecomparison::iSigma_l]  > 0.)    { Thosigma_l->setVal(inputPar[thecomparison::iSigma_l]);   Thosigma_l->setConstant(); }
    if (inputPar[thecomparison::iThoSigN] > 0.)     { Thon->setVal(inputPar[thecomparison::iThoSigN]);         Thon->setConstant(); }
    if (inputPar[thecomparison::iThoSigAlpha] > 0.) { Thoalpha->setVal(inputPar[thecomparison::iThoSigAlpha]); Thoalpha->setConstant(); }
     
    if (inputPar[thecomparison::iMean]  > 0.)       { Rm->setVal(inputPar[thecomparison::iMean]);  Rm->setConstant(); }
    if (inputPar[thecomparison::iSigma] > 0.)       { Rs->setVal(inputPar[thecomparison::iSigma]); Rs->setConstant(); }
    if (inputPar[thecomparison::iN] > 0.)           { Rn->setVal(inputPar[thecomparison::iN]);     Rn->setConstant(); }
    if (inputPar[thecomparison::iAlpha] > 0.)       { Ra->setVal(inputPar[thecomparison::iAlpha]); Ra->setConstant(); }
    if (inputPar[thecomparison::iEndpoint] > 0.)    { Rendpoint->setVal(inputPar[thecomparison::iEndpoint]); Rendpoint->setConstant(); }
  
    // argus will always be free!!!
    if (inputPar[thecomparison::iArgus] > -55.)     { pArgPar->setVal(inputPar[thecomparison::iArgus]); }
    if (inputPar[thecomparison::iCutOff] > 0.)      { pCutOff->setVal(inputPar[thecomparison::iCutOff]); pCutOff->setConstant(); }
  }

  std::cout<< "Parameters before fitting::" 
	   << " " << Thor->getVal() << " " << Thosigma_r1->getVal() << " " << Thoxc->getVal() << " " << Thosigma_r2->getVal() 
	   << " " << Thosigma_l->getVal() << " " << Thon->getVal() << " " << Thoalpha->getVal()
	   << " " << Rm->getVal() << " " << Rs->getVal() << " " << Ra->getVal() << " " << Rn->getVal() << " " << Rendpoint->getVal()
	   << " " << pArgPar->getVal() << " " << pCutOff->getVal()
	   << " " << nsig->getVal() << " " << nbpk->getVal() << " " << nbkg->getVal() << std::endl;
  

  RooFitResult* r = model->fitTo(*data, "rmhe");

  
  // get results and fill mesData object
  double p0 = nsig->getVal(); double Dp0 = nsig->getError();
  double p1 = nbkg->getVal(); double Dp1 = nbkg->getError();  
  p1 += nbpk->getVal(); 
  Dp1 = sqrt(Dp1*Dp1 + nbpk->getError()*nbpk->getError());
  if (data->sumEntries() == 0) p0 = Dp0 = p1 = Dp1 = 0.;
  
  mesData *pD = new mesData(data->GetName(),p0,Dp0,p1,Dp1,p0/(p0+p1),dBinomial(p0, p0+p1));
 
  // get back the best fit

  results.resize(14);
  if (results.size() != inputPar.size()) {
    std::cout << "Error: inconsistency in input/output parameters for mes fit! Exiting!" << std::cout;
    exit(EXIT_FAILURE);
  }

  results[thecomparison::iThoSigR]     = getVal(model, "ThoSigR");
  results[thecomparison::iSigma_r1]    = getVal(model, "sigma_r1");
  results[thecomparison::iThoSigXc]    = getVal(model, "ThoSigXc");
  results[thecomparison::iSigma_r2]    = getVal(model, "sigma_r2");
  results[thecomparison::iSigma_l]     = getVal(model, "sigma_l");
  results[thecomparison::iThoSigN]     = getVal(model, "ThoSigN");
  results[thecomparison::iThoSigAlpha] = getVal(model, "ThoSigAlpha");
     
  results[thecomparison::iMean]        = getVal(model, "ccbmean");
  results[thecomparison::iSigma]       = getVal(model, "ccbsigma");
  results[thecomparison::iN]           = getVal(model, "ccbalpha");
  results[thecomparison::iAlpha]       = getVal(model, "ccbn");
  results[thecomparison::iEndpoint]    = getVal(model, "ccbcutoff");

  results[thecomparison::iArgus]       = getVal(model, "ar");
  results[thecomparison::iCutOff]      = getVal(model, "cutoff");
  // clean up

  //  delete bdatared; bdatared = 0;

  // clean up model and associated pdfs/RooRealVars

  // get list of parameters and delete them
  RooArgSet* params = model->getParameters(RooArgSet());

   TIterator* iter = params->createIterator();
   RooAbsArg* var = 0;
   while ((var=(RooAbsArg*)iter->Next())) {
     if (std::string(var->GetName()) == "mes") continue;
     delete var;
   }
   delete params;;

   // get list of pdfs and delete them
   RooArgSet* comps = model->getComponents();

   iter = comps->createIterator();
   RooAbsArg* pdf = 0;
   while ((pdf=(RooAbsArg*)iter->Next())) {
     delete pdf;
   }
   delete iter; delete comps;

   model = 0;

  // done

   return pD;
}

void thecomparison::sighisto_newmodel(double& signal, double& signalErr,TH1D* histo, vector<double>& results, int fixpar, const vector<double>& inputPar, bool isData)
{
  mesData themes;
  
  // this should be mes
  RooRealVar *x=new RooRealVar("mes","mes", 5.2, 5.3);
  RooDataHist* rdh = new RooDataHist("binned","binned",RooArgList(*x),histo);  
 
  mesData tempthemes(*vubMesUnb(rdh,x,results,fixpar,inputPar,isData));
  themes=tempthemes;
     
  signal=themes.theSig(); signalErr=themes.theErrSig();
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
   
   c6 = new TCanvas("c6", "c6", 300, 0, 400,800);  
   // -- top left
   fPads[1]= new TPad("pad1", "", 0.00, 0.65, 0.99, 0.99);   fPads[1]->Draw(); 
   fPads[2]= new TPad("pad2", "", 0.00, 0.50, 0.99, 0.64);   fPads[2]->Draw();  
   // -- bottom left
   fPads[3]= new TPad("pad3", "", 0.00, 0.15, 0.99, 0.49);   fPads[3]->Draw(); 
   fPads[4]= new TPad("pad4", "", 0.00, 0.00, 0.99, 0.14);   fPads[4]->Draw(); 
   char name[100], line[100];
     //char name2[100], preine[100], preich[100], prefix[100], shift[100], smear[100], 
     //int h = 0 ; 

   sprintf(name,"%s%d" ,"h", 400000 + addcut * 1000000);
   TH1D y1;
   ((TH1D*)gDirectory->Get(name))->Copy(y1);
   sprintf(name,"%s%d" ,"h", 500000 + addcut * 1000000);
   TH1D y2;
   ((TH1D*)gDirectory->Get(name))->Copy(y2);
   sprintf(name,"%s%d" ,"d", 400000 + addcut * 1000000);
   TH1D y3;
   ((TH1D*)gDirectory->Get(name))->Copy(y3);
   sprintf(name,"%s%d" ,"d", 500000 + addcut * 1000000);
   TH1D y4;
   ((TH1D*)gDirectory->Get(name))->Copy(y4);

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

   double g1 = chisquared(&y1,&y2);
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

   fPads[3]->cd(); shrinkPad(0.001, 0.2); 

   setFilledHist(&y3, kBlack, kBlue, 3004);
   double inter3;
   double max2 = 1.4*y3.GetMaximum();
   y3.SetLabelSize(0.07, "Y");
   y3.SetMaximum(max2);
   y3.SetMinimum(0.);
   y3.SetMarkerSize(0.5);
   inter3 = y3.Integral();
// y3.SetNormFactor(inter1);
   y3.SetLineColor(kRed);
   y3.Draw();
 
   double inter4;
   inter4 = y4.Integral();
   y4.SetMarkerSize(0.5);
   if(cut) {
      intdatadepl = inter3;
      intMCdepl = inter4;
   }
   if (norm) {
      y4.Scale(intdatadepl/intMCdepl);
   }else{
      y4.Scale(inter3/inter4);
   }   

   setFilledHist(&y4 , kBlack, kRed, 3005);
   y4.DrawCopy("histsame");

   double g2 = chisquared(&y3,&y4);
   sprintf(line, "#chi^{2} = %5.4f", g2); 
   tl.SetTextSizePixels(50); tl.DrawLatex(0.22, 0.75, line);
  
   fPads[4]->cd(); shrinkPad(0.001, 0.2); 
   shrinkPad(0.4, 0.2, 0.1, 0.001);

   gPad->SetGridx(1);  gPad->SetGridy(1);
   TH1D *hratio2 = new TH1D(y4); hratio2->SetName("hratio2"); hratio2->Reset();
   hratio2->Divide(&y3, &y4);
   hratio2->SetMinimum(0.5); hratio2->SetMaximum(1.5);
   hratio2->SetMarkerStyle(24);
   hratio2->SetNdivisions(504, "Y");
   hratio2->SetLabelSize(0.22, "X");  hratio2->SetLabelSize(0.17, "Y");
   hratio2->SetStats(0);
   hratio2->SetTitle("");
   hratio2->Draw();

//    cout << "the chi square is " << g << endl;
//    SHIFTNEUT = SHIFTNEUT *1000;
//    SIGMANEUT = SIGMANEUT *1000;
//    sprintf(shift, "%s%d","-", SHIFTNEUT);
//    sprintf(smear, "%s%d","-", SIGMANEUT);   
  
   sprintf(mycut, "allcuts");   
   if (cut) sprintf(mycut, "leptoncuts");   
   if (norm) {   
      sprintf(name, "%s%s%d%s%s%s", dir.Data(), "/comparisonnorm",multipl,thevar.Data(),mycut,".ps");
   }else{
      sprintf(name, "%s%s%d%s%s%s", dir.Data(), "/comparison",multipl,thevar.Data(),mycut,".ps");
   }
   c6->SaveAs(name);  
//   sprintf(name, "%s%s%s%s%s%s%s", "comparison",thevar,"smear",smear,"shift",shift,".dat");
//    ofstream outfile(name);   
//    outfile << "the chi square for " << thevar << " with smear " << smear << " and shift " << shift << " is " << endl;
//    outfile << "  CHISQ "  << g << endl; 
}

double thecomparison::chisquared(TH1 *h1, TH1 *h2){
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

void thecomparison::effplots(TString dir){

   char name[100];
   sprintf(name,"%s%d" ,"h", 400000);
   TH1D *y1(((TH1D*)gDirectory->Get(name)));
   sprintf(name,"%s%d" ,"h", 1400000);
   TH1D *y2(((TH1D*)gDirectory->Get(name)));
   sprintf(name,"%s%d" ,"d", 400000);
   TH1D *y3(((TH1D*)gDirectory->Get(name)));
   sprintf(name,"%s%d" ,"d", 1400000);
   TH1D *y4(((TH1D*)gDirectory->Get(name)));
   sprintf(name,"%s%d" ,"h", 500000);
   TH1D *y5(((TH1D*)gDirectory->Get(name)));
   sprintf(name,"%s%d" ,"h", 1500000);
   TH1D *y6(((TH1D*)gDirectory->Get(name)));
   sprintf(name,"%s%d" ,"d", 500000);
   TH1D *y7(((TH1D*)gDirectory->Get(name)));
   sprintf(name,"%s%d" ,"d", 1500000);
   TH1D *y8(((TH1D*)gDirectory->Get(name)));

   TH1D *hratio1 = new TH1D(*y1); hratio1->SetName("hratio1"); hratio1->Reset();
   TH1D *hratio2 = new TH1D(*y5); hratio2->SetName("hratio2"); hratio2->Reset();
   TH1D *hratio3 = new TH1D(*y3); hratio3->SetName("hratio3"); hratio3->Reset();
   TH1D *hratio4 = new TH1D(*y7); hratio4->SetName("hratio4"); hratio4->Reset();
   TH1D *hratioratio1 = new TH1D(*hratio1); hratioratio1->SetName("hratioratio1"); hratioratio1->Reset();
   TH1D *hratioratio2 = new TH1D(*hratio3); hratioratio2->SetName("hratioratio2"); hratioratio2->Reset();

  
   hratio1->Divide(y1, y2, 1, 1, "B");
   hratio2->Divide(y5, y6, 1, 1, "B");
   hratio3->Divide(y3, y4, 1, 1, "B");
   hratio4->Divide(y7, y8, 1, 1, "B");
   hratioratio1->Divide(hratio1, hratio2);
   hratioratio2->Divide(hratio3, hratio4);

   c7 = new TCanvas("c7", "c7", 300,   0, 400,800);  
   // -- top left
   fPads[1]= new TPad("pad1", "", 0.00, 0.65, 0.99, 0.99);   fPads[1]->Draw(); 
   fPads[2]= new TPad("pad2", "", 0.00, 0.50, 0.99, 0.64);   fPads[2]->Draw();  
   // -- bottom left
   fPads[3]= new TPad("pad3", "", 0.00, 0.15, 0.99, 0.49);   fPads[3]->Draw(); 
   fPads[4]= new TPad("pad4", "", 0.00, 0.00, 0.99, 0.14);   fPads[4]->Draw(); 

   fPads[1]->cd(); shrinkPad(0.001, 0.2); 

//   double max = 1.4*hratio1->GetMaximum();
   hratio1->SetLabelSize(0.07, "Y");
//   hratio1->SetMaximum(max);
//cout << max << endl;
   hratio1->SetMaximum(1.);
   hratio1->SetMinimum(-0.1);
   setFilledHist(hratio1 , kBlack, kBlue, 3005);
   hratio1->Draw();

   setFilledHist(hratio2 , kBlack, kRed, 3005);
   hratio2->DrawCopy("histsame");
   
   fPads[2]->cd(); shrinkPad(0.001, 0.2); 
   shrinkPad(0.4, 0.2, 0.1, 0.001);

   gPad->SetGridx(1);  gPad->SetGridy(1);

   hratioratio1->SetMinimum(0.5); hratioratio1->SetMaximum(1.5);
   hratioratio1->SetMarkerStyle(24);
   hratioratio1->SetNdivisions(504, "Y");
   hratioratio1->SetLabelSize(0.22, "X");  hratioratio1->SetLabelSize(0.17, "Y");
   hratioratio1->SetStats(0);
   hratioratio1->SetTitle("");
   hratioratio1->Draw();

   fPads[3]->cd(); shrinkPad(0.001, 0.2); 

//   double max = 1.4*hratio1->GetMaximum();
   hratio3->SetLabelSize(0.07, "Y");
//   hratio1->SetMaximum(max);
//cout << max << endl;
   hratio3->SetMaximum(1.);
   hratio3->SetMinimum(-0.1);
   setFilledHist(hratio3 , kBlack, kBlue, 3005);
   hratio3->Draw();

   setFilledHist(hratio4 , kBlack, kRed, 3005);
   hratio4->DrawCopy("histsame");
   
   fPads[4]->cd(); shrinkPad(0.001, 0.2); 
   shrinkPad(0.4, 0.2, 0.1, 0.001);

   gPad->SetGridx(1);  gPad->SetGridy(1);

   hratioratio2->SetMinimum(0.5); hratioratio2->SetMaximum(1.5);
   hratioratio2->SetMarkerStyle(24);
   hratioratio2->SetNdivisions(504, "Y");
   hratioratio2->SetLabelSize(0.22, "X");  hratioratio2->SetLabelSize(0.17, "Y");
   hratioratio2->SetStats(0);
   hratioratio2->SetTitle("");
   hratioratio2->Draw();

   sprintf(name, "%s%s%d%s%s", dir.Data(), "/comparisoneff",multipl,thevar.Data(),".ps");
   c7->SaveAs(name);  
}

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

void thecomparison::Init(TTree *tree)
{
 // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses of the tree
   // will be set. It is normaly not necessary to make changes to the
   // generated code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running with PROOF.

   // Set branch addresses
   if (tree == 0) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run",&run);
   fChain->SetBranchAddress("upper",&upper);
   fChain->SetBranchAddress("lower",&lower);
   fChain->SetBranchAddress("Gvxbtyp",&Gvxbtyp);
   fChain->SetBranchAddress("GSem",&GSem);
   fChain->SetBranchAddress("GfDpi",&GfDpi);
   fChain->SetBranchAddress("GfDpiz",&GfDpiz);
   fChain->SetBranchAddress("GfDk",&GfDk);
   fChain->SetBranchAddress("GfDks",&GfDks);
   fChain->SetBranchAddress("GfDkl",&GfDkl);
   fChain->SetBranchAddress("GfDlep",&GfDlep);
   fChain->SetBranchAddress("GfDgam",&GfDgam);
   fChain->SetBranchAddress("GfDnu",&GfDnu);
   fChain->SetBranchAddress("GfD0Ds",&GfD0Ds);
   fChain->SetBranchAddress("GfDDs",&GfDDs);
   fChain->SetBranchAddress("GfDkspiopio",&GfDkspiopio);
   fChain->SetBranchAddress("GfK",&GfK);
   fChain->SetBranchAddress("isassocB",&isassocB);
   fChain->SetBranchAddress("isassocB_GHIT",&isassocB_GHIT);
   fChain->SetBranchAddress("ass_deltapB",&ass_deltapB);
   fChain->SetBranchAddress("isGoodMatch",&isGoodMatch);
   fChain->SetBranchAddress("ch1B",&ch1B);
   fChain->SetBranchAddress("ch2B",&ch2B);
   fChain->SetBranchAddress("chunm",&chunm);
   fChain->SetBranchAddress("neu1B",&neu1B);
   fChain->SetBranchAddress("neu2B",&neu2B);
   fChain->SetBranchAddress("neuunm",&neuunm);
   fChain->SetBranchAddress("brecoqual",&brecoqual);
   fChain->SetBranchAddress("brqual",&brqual);
   fChain->SetBranchAddress("brecoqualangle",&brecoqualangle);
   fChain->SetBranchAddress("chgdaugen",&chgdaugen);
   fChain->SetBranchAddress("neudaugen",&neudaugen);
   fChain->SetBranchAddress("nchg",&nchg);
   fChain->SetBranchAddress("nneu",&nneu);
   fChain->SetBranchAddress("xcharge",&xcharge);
   fChain->SetBranchAddress("nB",&nB);
   fChain->SetBranchAddress("brecoid",&brecoid);
   fChain->SetBranchAddress("brecoidtrue",&brecoidtrue);
   fChain->SetBranchAddress("brecoflav",&brecoflav);
   fChain->SetBranchAddress("brecocharge",&brecocharge);
   fChain->SetBranchAddress("modeB",&modeB);
   fChain->SetBranchAddress("truemodeB",&truemodeB);
   fChain->SetBranchAddress("isdoubleD",&isdoubleD);
   fChain->SetBranchAddress("mes",&mes);
   fChain->SetBranchAddress("mesendpoint",&mesendpoint);
   fChain->SetBranchAddress("de",&de);
   fChain->SetBranchAddress("pB",&pB);
   fChain->SetBranchAddress("eB",&eB);
   fChain->SetBranchAddress("eUps",&eUps);
   fChain->SetBranchAddress("pUps",&pUps);
   fChain->SetBranchAddress("thetaUps",&thetaUps);
   fChain->SetBranchAddress("phiUps",&phiUps);
   fChain->SetBranchAddress("thetaB",&thetaB);
   fChain->SetBranchAddress("phiB",&phiB);
   fChain->SetBranchAddress("pBtrue",&pBtrue);
   fChain->SetBranchAddress("eBtrue",&eBtrue);
   fChain->SetBranchAddress("tBtrue",&tBtrue);
   fChain->SetBranchAddress("fBtrue",&fBtrue);
   fChain->SetBranchAddress("pur",&pur);
   fChain->SetBranchAddress("intpur",&intpur);
   fChain->SetBranchAddress("nle_nopcut",&nle_nopcut);
   fChain->SetBranchAddress("nle",&nle);
   fChain->SetBranchAddress("nel",&nel);
   fChain->SetBranchAddress("nmu",&nmu);
   fChain->SetBranchAddress("nkp",&nkp);
   fChain->SetBranchAddress("nks",&nks);
   fChain->SetBranchAddress("nlept500",&nlept500);
   fChain->SetBranchAddress("nelec500",&nelec500);
   fChain->SetBranchAddress("nmu500",&nmu500);
   fChain->SetBranchAddress("nlept1000",&nlept1000);
   fChain->SetBranchAddress("nelec1000",&nelec1000);
   fChain->SetBranchAddress("nmu1000",&nmu1000);
   fChain->SetBranchAddress("deltam",&deltam);
   fChain->SetBranchAddress("MM1pr",&MM1pr);
   fChain->SetBranchAddress("MM2pr",&MM2pr);
   fChain->SetBranchAddress("MM3pr",&MM3pr);
   fChain->SetBranchAddress("OA1",&OA1);
   fChain->SetBranchAddress("OA2",&OA2);
   fChain->SetBranchAddress("OA3",&OA3);
   fChain->SetBranchAddress("PiMin1",&PiMin1);
   fChain->SetBranchAddress("PiMin2",&PiMin2);
   fChain->SetBranchAddress("PiMin3",&PiMin3);
   fChain->SetBranchAddress("plab",&plab);
   fChain->SetBranchAddress("elab",&elab);
   fChain->SetBranchAddress("tlab",&tlab);
   fChain->SetBranchAddress("flab",&flab);
   fChain->SetBranchAddress("plabgen",&plabgen);
   fChain->SetBranchAddress("elabgen",&elabgen);
   fChain->SetBranchAddress("tlabgen",&tlabgen);
   fChain->SetBranchAddress("flabgen",&flabgen);
   fChain->SetBranchAddress("lchargegen",&lchargegen);
   fChain->SetBranchAddress("pcms",&pcms);
   fChain->SetBranchAddress("ecms",&ecms);
   fChain->SetBranchAddress("tcms",&tcms);
   fChain->SetBranchAddress("fcms",&fcms);
   fChain->SetBranchAddress("lcharge",&lcharge);
   fChain->SetBranchAddress("leptidgen",&leptidgen);
   fChain->SetBranchAddress("leptorg",&leptorg);
   fChain->SetBranchAddress("isele",&isele);
   fChain->SetBranchAddress("vub",&vub);
   fChain->SetBranchAddress("vcb",&vcb);
   fChain->SetBranchAddress("other",&other);
   fChain->SetBranchAddress("nvubexcl",&nvubexcl);
   fChain->SetBranchAddress("nvubnres",&nvubnres);
   fChain->SetBranchAddress("ntau",&ntau);
   fChain->SetBranchAddress("mxhadgen",&mxhadgen);
   fChain->SetBranchAddress("mxhadgenwoph",&mxhadgenwoph);
   fChain->SetBranchAddress("xchargegen",&xchargegen);
   fChain->SetBranchAddress("pcmsgen",&pcmsgen);
   fChain->SetBranchAddress("ecmsgen",&ecmsgen);
   fChain->SetBranchAddress("tcmsgen",&tcmsgen);
   fChain->SetBranchAddress("fcmsgen",&fcmsgen);
   fChain->SetBranchAddress("pxhadgen",&pxhadgen);
   fChain->SetBranchAddress("exhadgen",&exhadgen);
   fChain->SetBranchAddress("exhadgencms",&exhadgencms);
   fChain->SetBranchAddress("txhadgen",&txhadgen);
   fChain->SetBranchAddress("fxhadgen",&fxhadgen);
   fChain->SetBranchAddress("q2Gen",&q2Gen);
   fChain->SetBranchAddress("ctvgen",&ctvgen);
   fChain->SetBranchAddress("ctlgen",&ctlgen);
   fChain->SetBranchAddress("chigen",&chigen);
   fChain->SetBranchAddress("enugencms",&enugencms);
   fChain->SetBranchAddress("pnugencms",&pnugencms);
   fChain->SetBranchAddress("mxhad",&mxhad);
   fChain->SetBranchAddress("emiss",&emiss);
   fChain->SetBranchAddress("pmiss",&pmiss);
   fChain->SetBranchAddress("mm2",&mm2);
   fChain->SetBranchAddress("q2",&q2);
   fChain->SetBranchAddress("q2new",&q2new);
   fChain->SetBranchAddress("exhadcms",&exhadcms);
   fChain->SetBranchAddress("pxhad",&pxhad);
   fChain->SetBranchAddress("exhad",&exhad);
   fChain->SetBranchAddress("txhad",&txhad);
   fChain->SetBranchAddress("fxhad",&fxhad);
   fChain->SetBranchAddress("pnu",&pnu);
   fChain->SetBranchAddress("tnu",&tnu);
   fChain->SetBranchAddress("fnu",&fnu);
   fChain->SetBranchAddress("pplus",&pplus);
   fChain->SetBranchAddress("pminus",&pminus);
   fChain->SetBranchAddress("pplusgen",&pplusgen);
   fChain->SetBranchAddress("pminusgen",&pminusgen);
   fChain->SetBranchAddress("pplusfit",&pplusfit);
   fChain->SetBranchAddress("pminusfit",&pminusfit);
   fChain->SetBranchAddress("wdeltam",&wdeltam);
   fChain->SetBranchAddress("mxhadfit",&mxhadfit);
   fChain->SetBranchAddress("mm2fit",&mm2fit);
   fChain->SetBranchAddress("q2fit",&q2fit);
   fChain->SetBranchAddress("chisq",&chisq);
   fChain->SetBranchAddress("globchisq",&globchisq);
   fChain->SetBranchAddress("probchisq",&probchisq);
   fChain->SetBranchAddress("ndof",&ndof);
   fChain->SetBranchAddress("fitstatus",&fitstatus);
   fChain->SetBranchAddress("nnpi0",&nnpi0);
   fChain->SetBranchAddress("nneu80_160",&nneu80_160);
   fChain->SetBranchAddress("nneu160_320",&nneu160_320);
   fChain->SetBranchAddress("pstarfitlept",&pstarfitlept);
   fChain->SetBranchAddress("pfitX",&pfitX);
   fChain->SetBranchAddress("thetafitX",&thetafitX);
   fChain->SetBranchAddress("phifitX",&phifitX);
   fChain->SetBranchAddress("pfitlept",&pfitlept);
   fChain->SetBranchAddress("thetafitlept",&thetafitlept);
   fChain->SetBranchAddress("phifitlept",&phifitlept);
   fChain->SetBranchAddress("pfitB",&pfitB);
   fChain->SetBranchAddress("thetafitB",&thetafitB);
   fChain->SetBranchAddress("phifitB",&phifitB);
   fChain->SetBranchAddress("totweight",&totweight);
   fChain->SetBranchAddress("totweightNutMult",&totweightNutMult);
   fChain->SetBranchAddress("totweightTrkMult",&totweightTrkMult);
   fChain->SetBranchAddress("kplus",&kplus);
   fChain->SetBranchAddress("nBrems",&nBrems);
   fChain->SetBranchAddress("eBrems",eBrems);
   fChain->SetBranchAddress("mxhadfit0",&mxhadfit0);
   fChain->SetBranchAddress("q2fit0",&q2fit0);
   fChain->SetBranchAddress("pplusfit0",&pplusfit0);
   fChain->SetBranchAddress("chisqfit0",&chisqfit0);
   fChain->SetBranchAddress("fitstatusfit0",&fitstatusfit0);
   fChain->SetBranchAddress("chisqfit",&chisqfit);
   fChain->SetBranchAddress("fitstatusfit",&fitstatusfit);
   fChain->SetBranchAddress("mxhadfit1",&mxhadfit1);
   fChain->SetBranchAddress("q2fit1",&q2fit1);
   fChain->SetBranchAddress("pplusfit1",&pplusfit1);
   fChain->SetBranchAddress("chisqfit1",&chisqfit1);
   fChain->SetBranchAddress("fitstatusfit1",&fitstatusfit1);
   fChain->SetBranchAddress("mxhadfit2",&mxhadfit2);
   fChain->SetBranchAddress("q2fit2",&q2fit2);
   fChain->SetBranchAddress("pplusfit2",&pplusfit2);
   fChain->SetBranchAddress("chisqfit2",&chisqfit2);
   fChain->SetBranchAddress("fitstatusfit2",&fitstatusfit2);
   fChain->SetBranchAddress("mxhadfit3",&mxhadfit3);
   fChain->SetBranchAddress("q2fit3",&q2fit3);
   fChain->SetBranchAddress("pplusfit3",&pplusfit3);
   fChain->SetBranchAddress("chisqfit3",&chisqfit3);
   fChain->SetBranchAddress("fitstatusfit3",&fitstatusfit3);
   fChain->SetBranchAddress("mxhadfit4",&mxhadfit4);
   fChain->SetBranchAddress("q2fit4",&q2fit4);
   fChain->SetBranchAddress("pplusfit4",&pplusfit4);
   fChain->SetBranchAddress("chisqfit4",&chisqfit4);
   fChain->SetBranchAddress("fitstatusfit4",&fitstatusfit4);
   fChain->SetBranchAddress("mxhadfit5",&mxhadfit5);
   fChain->SetBranchAddress("q2fit5",&q2fit5);
   fChain->SetBranchAddress("pplusfit5",&pplusfit5);
   fChain->SetBranchAddress("chisqfit5",&chisqfit5);
   fChain->SetBranchAddress("fitstatusfit5",&fitstatusfit5);
   fChain->SetBranchAddress("mxhadfit6",&mxhadfit6);
   fChain->SetBranchAddress("q2fit6",&q2fit6);
   fChain->SetBranchAddress("pplusfit6",&pplusfit6);
   fChain->SetBranchAddress("chisqfit6",&chisqfit6);
   fChain->SetBranchAddress("fitstatusfit6",&fitstatusfit6);
   fChain->SetBranchAddress("mxhadfit7",&mxhadfit7);
   fChain->SetBranchAddress("q2fit7",&q2fit7);
   fChain->SetBranchAddress("pplusfit7",&pplusfit7);
   fChain->SetBranchAddress("chisqfit7",&chisqfit7);
   fChain->SetBranchAddress("fitstatusfit7",&fitstatusfit7);
   fChain->SetBranchAddress("mxhadfit8",&mxhadfit8);
   fChain->SetBranchAddress("q2fit8",&q2fit8);
   fChain->SetBranchAddress("pplusfit8",&pplusfit8);
   fChain->SetBranchAddress("chisqfit8",&chisqfit8);
   fChain->SetBranchAddress("fitstatusfit8",&fitstatusfit8);
   Notify();
}

Bool_t thecomparison::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. Typically here the branch pointers
   // will be retrieved. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed.

   // Get branch pointers
   b_run = fChain->GetBranch("run");
   b_upper = fChain->GetBranch("upper");
   b_lower = fChain->GetBranch("lower");
   b_Gvxbtyp = fChain->GetBranch("Gvxbtyp");
   b_GSem = fChain->GetBranch("GSem");
   b_GfDpi = fChain->GetBranch("GfDpi");
   b_GfDpiz = fChain->GetBranch("GfDpiz");
   b_GfDk = fChain->GetBranch("GfDk");
   b_GfDks = fChain->GetBranch("GfDks");
   b_GfDkl = fChain->GetBranch("GfDkl");
   b_GfDlep = fChain->GetBranch("GfDlep");
   b_GfDgam = fChain->GetBranch("GfDgam");
   b_GfDnu = fChain->GetBranch("GfDnu");
   b_GfD0Ds = fChain->GetBranch("GfD0Ds");
   b_GfDDs = fChain->GetBranch("GfDDs");
   b_GfDkspiopio = fChain->GetBranch("GfDkspiopio");
   b_GfK = fChain->GetBranch("GfK");
   b_isassocB = fChain->GetBranch("isassocB");
   b_isassocB_GHIT = fChain->GetBranch("isassocB_GHIT");
   b_ass_deltapB = fChain->GetBranch("ass_deltapB");
   b_isGoodMatch = fChain->GetBranch("isGoodMatch");
   b_ch1B = fChain->GetBranch("ch1B");
   b_ch2B = fChain->GetBranch("ch2B");
   b_chunm = fChain->GetBranch("chunm");
   b_neu1B = fChain->GetBranch("neu1B");
   b_neu2B = fChain->GetBranch("neu2B");
   b_neuunm = fChain->GetBranch("neuunm");
   b_brecoqual = fChain->GetBranch("brecoqual");
   b_brqual = fChain->GetBranch("brqual");
   b_brecoqualangle = fChain->GetBranch("brecoqualangle");
   b_chgdaugen = fChain->GetBranch("chgdaugen");
   b_neudaugen = fChain->GetBranch("neudaugen");
   b_nchg = fChain->GetBranch("nchg");
   b_nneu = fChain->GetBranch("nneu");
   b_xcharge = fChain->GetBranch("xcharge");
   b_nB = fChain->GetBranch("nB");
   b_brecoid = fChain->GetBranch("brecoid");
   b_brecoidtrue = fChain->GetBranch("brecoidtrue");
   b_brecoflav = fChain->GetBranch("brecoflav");
   b_brecocharge = fChain->GetBranch("brecocharge");
   b_modeB = fChain->GetBranch("modeB");
   b_truemodeB = fChain->GetBranch("truemodeB");
   b_isdoubleD = fChain->GetBranch("isdoubleD");
   b_mes = fChain->GetBranch("mes");
   b_mesendpoint = fChain->GetBranch("mesendpoint");
   b_de = fChain->GetBranch("de");
   b_pB = fChain->GetBranch("pB");
   b_eB = fChain->GetBranch("eB");
   b_eUps = fChain->GetBranch("eUps");
   b_pUps = fChain->GetBranch("pUps");
   b_thetaUps = fChain->GetBranch("thetaUps");
   b_phiUps = fChain->GetBranch("phiUps");
   b_thetaB = fChain->GetBranch("thetaB");
   b_phiB = fChain->GetBranch("phiB");
   b_pBtrue = fChain->GetBranch("pBtrue");
   b_eBtrue = fChain->GetBranch("eBtrue");
   b_tBtrue = fChain->GetBranch("tBtrue");
   b_fBtrue = fChain->GetBranch("fBtrue");
   b_pur = fChain->GetBranch("pur");
   b_intpur = fChain->GetBranch("intpur");
   b_nle_nopcut = fChain->GetBranch("nle_nopcut");
   b_nle = fChain->GetBranch("nle");
   b_nel = fChain->GetBranch("nel");
   b_nmu = fChain->GetBranch("nmu");
   b_nkp = fChain->GetBranch("nkp");
   b_nks = fChain->GetBranch("nks");
   b_nlept500 = fChain->GetBranch("nlept500");
   b_nelec500 = fChain->GetBranch("nelec500");
   b_nmu500 = fChain->GetBranch("nmu500");
   b_nlept1000 = fChain->GetBranch("nlept1000");
   b_nelec1000 = fChain->GetBranch("nelec1000");
   b_nmu1000 = fChain->GetBranch("nmu1000");
   b_deltam = fChain->GetBranch("deltam");
   b_MM1pr = fChain->GetBranch("MM1pr");
   b_MM2pr = fChain->GetBranch("MM2pr");
   b_MM3pr = fChain->GetBranch("MM3pr");
   b_OA1 = fChain->GetBranch("OA1");
   b_OA2 = fChain->GetBranch("OA2");
   b_OA3 = fChain->GetBranch("OA3");
   b_PiMin1 = fChain->GetBranch("PiMin1");
   b_PiMin2 = fChain->GetBranch("PiMin2");
   b_PiMin3 = fChain->GetBranch("PiMin3");
   b_plab = fChain->GetBranch("plab");
   b_elab = fChain->GetBranch("elab");
   b_tlab = fChain->GetBranch("tlab");
   b_flab = fChain->GetBranch("flab");
   b_plabgen = fChain->GetBranch("plabgen");
   b_elabgen = fChain->GetBranch("elabgen");
   b_tlabgen = fChain->GetBranch("tlabgen");
   b_flabgen = fChain->GetBranch("flabgen");
   b_lchargegen = fChain->GetBranch("lchargegen");
   b_pcms = fChain->GetBranch("pcms");
   b_ecms = fChain->GetBranch("ecms");
   b_tcms = fChain->GetBranch("tcms");
   b_fcms = fChain->GetBranch("fcms");
   b_lcharge = fChain->GetBranch("lcharge");
   b_leptidgen = fChain->GetBranch("leptidgen");
   b_leptorg = fChain->GetBranch("leptorg");
   b_isele = fChain->GetBranch("isele");
   b_vub = fChain->GetBranch("vub");
   b_vcb = fChain->GetBranch("vcb");
   b_other = fChain->GetBranch("other");
   b_nvubexcl = fChain->GetBranch("nvubexcl");
   b_nvubnres = fChain->GetBranch("nvubnres");
   b_ntau = fChain->GetBranch("ntau");
   b_mxhadgen = fChain->GetBranch("mxhadgen");
   b_mxhadgenwoph = fChain->GetBranch("mxhadgenwoph");
   b_xchargegen = fChain->GetBranch("xchargegen");
   b_pcmsgen = fChain->GetBranch("pcmsgen");
   b_ecmsgen = fChain->GetBranch("ecmsgen");
   b_tcmsgen = fChain->GetBranch("tcmsgen");
   b_fcmsgen = fChain->GetBranch("fcmsgen");
   b_pxhadgen = fChain->GetBranch("pxhadgen");
   b_exhadgen = fChain->GetBranch("exhadgen");
   b_exhadgencms = fChain->GetBranch("exhadgencms");
   b_txhadgen = fChain->GetBranch("txhadgen");
   b_fxhadgen = fChain->GetBranch("fxhadgen");
   b_q2Gen = fChain->GetBranch("q2Gen");
   b_ctvgen = fChain->GetBranch("ctvgen");
   b_ctlgen = fChain->GetBranch("ctlgen");
   b_chigen = fChain->GetBranch("chigen");
   b_enugencms = fChain->GetBranch("enugencms");
   b_pnugencms = fChain->GetBranch("pnugencms");
   b_mxhad = fChain->GetBranch("mxhad");
   b_emiss = fChain->GetBranch("emiss");
   b_pmiss = fChain->GetBranch("pmiss");
   b_mm2 = fChain->GetBranch("mm2");
   b_q2 = fChain->GetBranch("q2");
   b_q2new = fChain->GetBranch("q2new");
   b_exhadcms = fChain->GetBranch("exhadcms");
   b_pxhad = fChain->GetBranch("pxhad");
   b_exhad = fChain->GetBranch("exhad");
   b_txhad = fChain->GetBranch("txhad");
   b_fxhad = fChain->GetBranch("fxhad");
   b_pnu = fChain->GetBranch("pnu");
   b_tnu = fChain->GetBranch("tnu");
   b_fnu = fChain->GetBranch("fnu");
   b_pplus = fChain->GetBranch("pplus");
   b_pminus = fChain->GetBranch("pminus");
   b_pplusgen = fChain->GetBranch("pplusgen");
   b_pminusgen = fChain->GetBranch("pminusgen");
   b_pplusfit = fChain->GetBranch("pplusfit");
   b_pminusfit = fChain->GetBranch("pminusfit");
   b_wdeltam = fChain->GetBranch("wdeltam");
   b_mxhadfit = fChain->GetBranch("mxhadfit");
   b_mm2fit = fChain->GetBranch("mm2fit");
   b_q2fit = fChain->GetBranch("q2fit");
   b_chisq = fChain->GetBranch("chisq");
   b_globchisq = fChain->GetBranch("globchisq");
   b_probchisq = fChain->GetBranch("probchisq");
   b_ndof = fChain->GetBranch("ndof");
   b_fitstatus = fChain->GetBranch("fitstatus");
   b_nnpi0 = fChain->GetBranch("nnpi0");
   b_nneu80_160 = fChain->GetBranch("nneu80_160");
   b_nneu160_320 = fChain->GetBranch("nneu160_320");
   b_pstarfitlept = fChain->GetBranch("pstarfitlept");
   b_pfitX = fChain->GetBranch("pfitX");
   b_thetafitX = fChain->GetBranch("thetafitX");
   b_phifitX = fChain->GetBranch("phifitX");
   b_pfitlept = fChain->GetBranch("pfitlept");
   b_thetafitlept = fChain->GetBranch("thetafitlept");
   b_phifitlept = fChain->GetBranch("phifitlept");
   b_pfitB = fChain->GetBranch("pfitB");
   b_thetafitB = fChain->GetBranch("thetafitB");
   b_phifitB = fChain->GetBranch("phifitB");
   b_totweight = fChain->GetBranch("totweight");
   b_totweightNutMult = fChain->GetBranch("totweightNutMult");
   b_totweightTrkMult = fChain->GetBranch("totweightTrkMult");
   b_kplus = fChain->GetBranch("kplus");
   b_nBrems = fChain->GetBranch("nBrems");
   b_eBrems = fChain->GetBranch("eBrems");
   b_mxhadfit0 = fChain->GetBranch("mxhadfit0");
   b_q2fit0 = fChain->GetBranch("q2fit0");
   b_pplusfit0 = fChain->GetBranch("pplusfit0");
   b_chisqfit0 = fChain->GetBranch("chisqfit0");
   b_fitstatusfit0 = fChain->GetBranch("fitstatusfit0");
   b_chisqfit = fChain->GetBranch("chisqfit");
   b_fitstatusfit = fChain->GetBranch("fitstatusfit");
   b_mxhadfit1 = fChain->GetBranch("mxhadfit1");
   b_q2fit1 = fChain->GetBranch("q2fit1");
   b_pplusfit1 = fChain->GetBranch("pplusfit1");
   b_chisqfit1 = fChain->GetBranch("chisqfit1");
   b_fitstatusfit1 = fChain->GetBranch("fitstatusfit1");
   b_mxhadfit2 = fChain->GetBranch("mxhadfit2");
   b_q2fit2 = fChain->GetBranch("q2fit2");
   b_pplusfit2 = fChain->GetBranch("pplusfit2");
   b_chisqfit2 = fChain->GetBranch("chisqfit2");
   b_fitstatusfit2 = fChain->GetBranch("fitstatusfit2");
   b_mxhadfit3 = fChain->GetBranch("mxhadfit3");
   b_q2fit3 = fChain->GetBranch("q2fit3");
   b_pplusfit3 = fChain->GetBranch("pplusfit3");
   b_chisqfit3 = fChain->GetBranch("chisqfit3");
   b_fitstatusfit3 = fChain->GetBranch("fitstatusfit3");
   b_mxhadfit4 = fChain->GetBranch("mxhadfit4");
   b_q2fit4 = fChain->GetBranch("q2fit4");
   b_pplusfit4 = fChain->GetBranch("pplusfit4");
   b_chisqfit4 = fChain->GetBranch("chisqfit4");
   b_fitstatusfit4 = fChain->GetBranch("fitstatusfit4");
   b_mxhadfit5 = fChain->GetBranch("mxhadfit5");
   b_q2fit5 = fChain->GetBranch("q2fit5");
   b_pplusfit5 = fChain->GetBranch("pplusfit5");
   b_chisqfit5 = fChain->GetBranch("chisqfit5");
   b_fitstatusfit5 = fChain->GetBranch("fitstatusfit5");
   b_mxhadfit6 = fChain->GetBranch("mxhadfit6");
   b_q2fit6 = fChain->GetBranch("q2fit6");
   b_pplusfit6 = fChain->GetBranch("pplusfit6");
   b_chisqfit6 = fChain->GetBranch("chisqfit6");
   b_fitstatusfit6 = fChain->GetBranch("fitstatusfit6");
   b_mxhadfit7 = fChain->GetBranch("mxhadfit7");
   b_q2fit7 = fChain->GetBranch("q2fit7");
   b_pplusfit7 = fChain->GetBranch("pplusfit7");
   b_chisqfit7 = fChain->GetBranch("chisqfit7");
   b_fitstatusfit7 = fChain->GetBranch("fitstatusfit7");
   b_mxhadfit8 = fChain->GetBranch("mxhadfit8");
   b_q2fit8 = fChain->GetBranch("q2fit8");
   b_pplusfit8 = fChain->GetBranch("pplusfit8");
   b_chisqfit8 = fChain->GetBranch("chisqfit8");
   b_fitstatusfit8 = fChain->GetBranch("fitstatusfit8");

   return kTRUE;
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
RooAbsPdf* thecomparison::createThorsten(RooRealVar& mes)
{
  // modified gauz parameters
  Thor =        new RooRealVar("ThoSigR","thosig r",0.9,0.2,1);
  Thosigma_r1 = new RooRealVar("sigma_r1","sigma_r1",0.00235,0.00130,0.00190);
  Thoxc =       new RooRealVar("ThoSigXc","Thosig xc",5.27966,5.2790,5.2810);
  Thosigma_r2 = new RooRealVar("sigma_r2","sigma_r2",0.00304,0.0015,0.0035);
  Thosigma_l  = new RooRealVar("sigma_l","sigma_l",0.0017,0.00160,0.00280);
  Thon =        new RooRealVar("ThoSigN","thosig n",1.2,0.5,2.5);
  Thoalpha =    new RooRealVar("ThoSigAlpha","thosig alpha",4.2,3.4,4.5);

  RooThorstenSig* pFunction = new RooThorstenSig("ThoSig","Thorstens Triple Function", mes,
						 *Thor, *Thosigma_r1, *Thoxc, *Thosigma_r2, *Thosigma_l, *Thon, *Thoalpha);

  return pFunction;
}
RooAbsPdf* thecomparison::createCCB(RooRealVar& mes)
{
  // return function
  RooAbsPdf* pFunction = 0;

  // modified crystal ball parameters (same for all functions)
  Rm = new RooRealVar("ccbmean","mean of gaussian 1", 5.2795, 5.2789, 5.2806);
  Rs = new RooRealVar("ccbsigma","width of gaussians", 0.003, 0.001, 0.007);
  Ra = new RooRealVar("ccbalpha","alpha parameter", 0.529, 0.1, 0.7);
  Rn = new RooRealVar("ccbn","n parameter", 5., 0.1, 20.);
  Rendpoint = new RooRealVar("ccbcutoff", "ccb cutoff", 5.2891, 5.288, 5.2944);

  pFunction = new RooCCB("ccb", "Modified Crystal Ball (ccb)", mes, *Rm, *Rs, *Ra, *Rn, *Rendpoint);

  return pFunction;
}

RooAbsPdf* thecomparison::createArgus(RooRealVar& mes)
{
  // return function
  RooAbsPdf* pFunction = 0;

  // argus parameters (same for all functions)
  pArgPar = new RooRealVar("ar", "argus shape parameter", -60., -100., -10.);
  pCutOff = new RooRealVar("cutoff", "argus cutoff", 5.2891, 5.28, 5.292);
  pFunction = new RooArgusBG("a", "Argus PDF", mes, *pCutOff, *pArgPar);

  return pFunction;
}

//! helper access method used by vubMesUnb
Double_t getVal(RooAbsPdf* pdf, const char* name)
{
  RooRealVar* var = getPointer(pdf, name);
  return (var != 0 ? var->getVal() : 0.);
}
RooRealVar* getPointer(RooAbsPdf* pdf, const char* name)
{
  RooArgSet* list = pdf->getParameters(RooArgSet());
  RooRealVar* var = dynamic_cast<RooRealVar*>(list->find(name));
  delete list;

  return var;
}

// HERE ARE LISTED THE INIT() AND MODIFY METHODS FOR NTUPLES BEFORE SUMMER 06 PRODUCTION

// void thecomparison::Init(TTree *tree)
// {
//   //   Set branch addresses
//    if (tree == 0) return;
//    fChain    = tree;
//    fCurrent = -1;
//    fChain->SetMakeClass(1);

//    fChain->SetBranchAddress("run",&run);
//    fChain->SetBranchAddress("lower",&lower);
//    fChain->SetBranchAddress("upper",&upper);
//    /*   fChain->SetBranchAddress("bmass",&bmass);
//    fChain->SetBranchAddress("bmassfit",&bmassfit);
//    fChain->SetBranchAddress("sbox",&sbox);*/
//    fChain->SetBranchAddress("mes",&mes);
//    fChain->SetBranchAddress("de",&de);
//    fChain->SetBranchAddress("pur",&pur);
//    fChain->SetBranchAddress("intpur",&intpur);
//    fChain->SetBranchAddress("modeB",&modeB);
//    fChain->SetBranchAddress("nnpi0",&nnpi0);
//    fChain->SetBranchAddress("GSem",&GSem);  
//    fChain->SetBranchAddress("GfDpi",&GfDpi); 
//    fChain->SetBranchAddress("GfDpiz",&GfDpiz);
//    fChain->SetBranchAddress("GfDk",&GfDk);  
//    fChain->SetBranchAddress("GfDks",&GfDks); 
//    fChain->SetBranchAddress("GfDlep",&GfDlep);
//    fChain->SetBranchAddress("GfDgam",&GfDgam);
//    fChain->SetBranchAddress("brecoflav",&brecoflav);
//    fChain->SetBranchAddress("brecocharge",&brecocharge);
//    //   fChain->SetBranchAddress("brecomc",&brecomc);
//    fChain->SetBranchAddress("mxhadgen",&mxhadgen);
//    fChain->SetBranchAddress("pcmsgen",&pcmsgen);
//    fChain->SetBranchAddress("tcmsgen",&tcmsgen);
//    fChain->SetBranchAddress("fcmsgen",&fcmsgen);
//    fChain->SetBranchAddress("ecmsgen",&ecmsgen);
//    fChain->SetBranchAddress("pxhadgen",&pxhadgen);
//    fChain->SetBranchAddress("txhadgen",&txhadgen);
//    fChain->SetBranchAddress("fxhadgen",&fxhadgen);
//    fChain->SetBranchAddress("exhadgen",&exhadgen);
//    fChain->SetBranchAddress("kplus",&fkplus);
//    /*   fChain->SetBranchAddress("GoodEvent",&GoodEvent);
//    fChain->SetBranchAddress("isDupli",&isDupli);
//    fChain->SetBranchAddress("ValMap",&ValMap);*/
//    fChain->SetBranchAddress("vub",&vub);
//    fChain->SetBranchAddress("vcb",&vcb);
//    fChain->SetBranchAddress("Gvxbtyp",&Gvxbtyp);
//    fChain->SetBranchAddress("other",&other);
//    //   fChain->SetBranchAddress("bgcat",&bgcat);
//    fChain->SetBranchAddress("xcharge",&xcharge);
//    fChain->SetBranchAddress("pxhad",&pxhad);
//    fChain->SetBranchAddress("txhad",&txhad);
//    fChain->SetBranchAddress("fxhad",&fxhad);
//    fChain->SetBranchAddress("exhad",&exhad);
//    fChain->SetBranchAddress("mxhad",&mxhad);
//    fChain->SetBranchAddress("mxhadfit",&mxhadfit);
//    /*
//      fChain->SetBranchAddress("gmax",&gmax);
//      fChain->SetBranchAddress("EwPwfit",&EwPwfit);
//    fChain->SetBranchAddress("csiCiuc",&csiCiuc);
//    fChain->SetBranchAddress("xCiuc",&xCiuc);
//    fChain->SetBranchAddress("wCiuc",&wCiuc);
//    fChain->SetBranchAddress("GcsiCiuc",&GcsiCiuc);
//    fChain->SetBranchAddress("GxCiuc",&GxCiuc);
//    fChain->SetBranchAddress("GwCiuc",&GwCiuc);
//    fChain->SetBranchAddress("EwPwG",&EwPwG);*/
//    fChain->SetBranchAddress("q2fit",&q2fit);
//    fChain->SetBranchAddress("q2",&q2);
//    fChain->SetBranchAddress("q2Gen",&q2Gen);

//    fChain->SetBranchAddress("lcharge",&lcharge);
//    fChain->SetBranchAddress("plab",&plab);
//    fChain->SetBranchAddress("tlab",&tlab);
//    fChain->SetBranchAddress("flab",&flab);
//    fChain->SetBranchAddress("pcms",&pcms);
//    fChain->SetBranchAddress("tcms",&tcms);
//    fChain->SetBranchAddress("fcms",&fcms);
//    fChain->SetBranchAddress("ecms",&ecms);
//    fChain->SetBranchAddress("nle",&nle);
//    fChain->SetBranchAddress("nel",&nel);
//    fChain->SetBranchAddress("nmu",&nmu);
//    fChain->SetBranchAddress("nchg",&nchg);
//    fChain->SetBranchAddress("nneu",&nneu);
//    fChain->SetBranchAddress("nneu80_160",&nneu80_160);
//    fChain->SetBranchAddress("nneu160_320",&nneu160_320);
//    //   fChain->SetBranchAddress("nneufromB",&nneufromB);
//    // fChain->SetBranchAddress("nneufromB80_160",&nneufromB80_160);
//    // fChain->SetBranchAddress("nneufromB160_320",&nneufromB160_320);
//    fChain->SetBranchAddress("nkp",&nkp);
//    fChain->SetBranchAddress("nks",&nks);
//    //   fChain->SetBranchAddress("npi0",&npi0);
//    fChain->SetBranchAddress("pnu",&pnu);
//    /*   fChain->SetBranchAddress("dx",&dx);
//    fChain->SetBranchAddress("dy",&dy);
//    fChain->SetBranchAddress("dz",&dz);
//    fChain->SetBranchAddress("s2dxx",&s2dxx);
//    fChain->SetBranchAddress("s2dyy",&s2dyy);
//    fChain->SetBranchAddress("s2dzz",&s2dzz);
//    fChain->SetBranchAddress("s2dxy",&s2dxy);
//    fChain->SetBranchAddress("s2dyz",&s2dyz);
//    fChain->SetBranchAddress("s2dxz",&s2dxz);
//    fChain->SetBranchAddress("EPiz",&epiz);
//    fChain->SetBranchAddress("MinKMom",&kminmom);
//    fChain->SetBranchAddress("MaxKMom",&kmaxmom);
//    fChain->SetBranchAddress("mm2nc",&mm2nc);
// */
//    fChain->SetBranchAddress("tnu",&tnu);
//    fChain->SetBranchAddress("fnu",&fnu);
//    //   fChain->SetBranchAddress("ENeu",&eneu);
   
//    fChain->SetBranchAddress("mm2",&mm2);

//    fChain->SetBranchAddress("mm2fit",&mm2fit);
//    /*   fChain->SetBranchAddress("allksm0",allksm0);
//    fChain->SetBranchAddress("allksp",allksp);
//    fChain->SetBranchAddress("allksmc",allksmc);
//    fChain->SetBranchAddress("allchkp",allchkp);
//    fChain->SetBranchAddress("allchkmc",allchkmc);
//    fChain->SetBranchAddress("m0ks",&m0ks);
//    fChain->SetBranchAddress("pks",&pks);
//    fChain->SetBranchAddress("pksmc",&pksmc);
//    fChain->SetBranchAddress("ntkl",&ntkl);
//    fChain->SetBranchAddress("tklp",&tklp);
//    fChain->SetBranchAddress("tklth",&tklth);
//    fChain->SetBranchAddress("tklph",&tklph);
//    fChain->SetBranchAddress("tklisol",&tklisol);
//    fChain->SetBranchAddress("ntks",&ntks);
//    fChain->SetBranchAddress("tksp",&tksp);
//    fChain->SetBranchAddress("tksth",&tksth);
//    fChain->SetBranchAddress("tksph",&tksph);
//    fChain->SetBranchAddress("tksdec",&tksdec);
//    fChain->SetBranchAddress("ntchk",&ntchk);
//    fChain->SetBranchAddress("tchkp",&tchkp);
//    fChain->SetBranchAddress("tchkth",&tchkth);
//    fChain->SetBranchAddress("tchkph",&tchkph);
//    fChain->SetBranchAddress("nklres",&nklres);
//    fChain->SetBranchAddress("klresth",&klresth);
//    fChain->SetBranchAddress("klresph",&klresph);
//    fChain->SetBranchAddress("klid",&klid);
//    fChain->SetBranchAddress("klcone",&klcone);
//    fChain->SetBranchAddress("emckl",&emckl);
//    fChain->SetBranchAddress("emckl0",&emckl0);
//    fChain->SetBranchAddress("emckl22",&emckl22);
//    fChain->SetBranchAddress("mxks",&mxks);
//    fChain->SetBranchAddress("mm2ks",&mm2ks);
//    fChain->SetBranchAddress("mxksfit",&mxksfit);
//    fChain->SetBranchAddress("mm2misk",&mm2misk);
//    fChain->SetBranchAddress("mxmisk",&mxmisk);
//    fChain->SetBranchAddress("mxmiskfit",&mxmiskfit);
//    fChain->SetBranchAddress("mm2chk",&mm2mchk);
//    fChain->SetBranchAddress("mxchk",&mxchk);
//    fChain->SetBranchAddress("mxchkfit",&mxchkfit);*/
//    fChain->SetBranchAddress("totweight",&totweight);
//    fChain->SetBranchAddress("totweightNutMult",&totweightNutMult);
//    fChain->SetBranchAddress("totweightTrkMult",&totweightTrkMult);
//    Notify();
// }

// Bool_t thecomparison::Notify()
// {
// //   called when loading a new file
// //   get branch pointers
//    b_run = fChain->GetBranch("run");
//    b_lower = fChain->GetBranch("lower");
//    b_upper = fChain->GetBranch("upper");
//    /*   b_bmass = fChain->GetBranch("bmass");
//    b_bmassfit = fChain->GetBranch("bmassfit");
//    b_sbox = fChain->GetBranch("sbox");*/
//    b_mes = fChain->GetBranch("mes");
//    b_de = fChain->GetBranch("de");
//    b_pur = fChain->GetBranch("pur");
//    b_intpur = fChain->GetBranch("intpur");
//    b_modeB = fChain->GetBranch("modeB");
//    b_nnpi0 = fChain->GetBranch("nnpi0");
//    b_brecoflav = fChain->GetBranch("brecoflav");
//    b_brecocharge = fChain->GetBranch("brecocharge");
//    //   b_brecomc = fChain->GetBranch("brecomc");
//    b_mxhadgen = fChain->GetBranch("mxhadgen");
//    b_pcmsgen = fChain->GetBranch("pcmsgen");
//    b_tcmsgen = fChain->GetBranch("tcmsgen");
//    b_fcmsgen = fChain->GetBranch("fcmsgen");
//    b_ecmsgen = fChain->GetBranch("ecmsgen");
//    b_pxhadgen = fChain->GetBranch("pxhadgen");
//    b_txhadgen = fChain->GetBranch("txhadgen");
//    b_fxhadgen = fChain->GetBranch("fxhadgen");
//    b_exhadgen = fChain->GetBranch("exhadgen");
//    b_kplus = fChain->GetBranch("kplus");
//    /*   b_GoodEvent = fChain->GetBranch("GoodEvent");
//    b_isDupli = fChain->GetBranch("isDupli");
//    b_ValMap = fChain->GetBranch("ValMap");*/
//    b_vub = fChain->GetBranch("vub"); 
//    b_vcb = fChain->GetBranch("vcb");
//    b_Gvxbtyp = fChain->GetBranch("Gvxbtyp");
//    b_other = fChain->GetBranch("other");
//    //   b_bgcat = fChain->GetBranch("bgcat");
//    b_xcharge = fChain->GetBranch("xcharge");
//    b_pxhad = fChain->GetBranch("pxhad");
//    b_txhad = fChain->GetBranch("txhad");
//    b_fxhad = fChain->GetBranch("fxhad");
//    b_exhad = fChain->GetBranch("exhad");
//    b_mxhad = fChain->GetBranch("mxhad");
//    //   b_gmax = fChain->GetBranch("gmax");
//    b_mxhadfit = fChain->GetBranch("mxhadfit");
//    b_q2fit = fChain->GetBranch("q2fit");
//    b_q2 = fChain->GetBranch("q2");
//    //   b_EwPwfit = fChain->GetBranch("EwPwfit");
//    //b_csiCiuc = fChain->GetBranch("csiCiuc");
//    //b_xCiuc = fChain->GetBranch("xCiuc");
//    //b_wCiuc = fChain->GetBranch("wCiuc");
//    b_q2Gen = fChain->GetBranch("q2Gen");
//    //b_EwPwG = fChain->GetBranch("EwPwG");
//    //b_GcsiCiuc = fChain->GetBranch("GcsiCiuc");
//    //b_GxCiuc = fChain->GetBranch("GxCiuc");
//    //b_GwCiuc = fChain->GetBranch("GwCiuc");
//    b_lcharge = fChain->GetBranch("lcharge");
//    b_plab = fChain->GetBranch("plab");
//    b_tlab = fChain->GetBranch("tlab");
//    b_flab = fChain->GetBranch("flab");
//    b_pcms = fChain->GetBranch("pcms");
//    b_tcms = fChain->GetBranch("tcms");
//    b_fcms = fChain->GetBranch("fcms");
//    b_ecms = fChain->GetBranch("ecms");
//    b_nle = fChain->GetBranch("nle");
//    b_nel = fChain->GetBranch("nel");
//    b_nmu = fChain->GetBranch("nmu");
//    b_nchg = fChain->GetBranch("nchg");
//    b_nneu = fChain->GetBranch("nneu");
//    b_nneu80_160 = fChain->GetBranch("nneu80_160");
//    b_nneu160_320 = fChain->GetBranch("nneu160_320");
//    //   b_nneufromB = fChain->GetBranch("nneufromB");
//    // b_nneufromB80_160 = fChain->GetBranch("nneufromB80_160");
//    // b_nneufromB160_320 = fChain->GetBranch("nneufromB160_320");
//    b_nkp = fChain->GetBranch("nkp");
//    b_nks = fChain->GetBranch("nks");
//    //   b_npi0 = fChain->GetBranch("npi0");
//    b_pnu = fChain->GetBranch("pnu");
//    /*   b_dx = fChain->GetBranch("dx");
//    b_dy = fChain->GetBranch("dy");
//    b_dz = fChain->GetBranch("dz");
//    b_s2dxx = fChain->GetBranch("s2dxx");
//    b_s2dyy = fChain->GetBranch("s2dyy");
//    b_s2dzz = fChain->GetBranch("s2dzz");
//    b_s2dxy = fChain->GetBranch("s2dxy");
//    b_s2dyz = fChain->GetBranch("s2dyz");
//    b_s2dxz = fChain->GetBranch("s2dxz");
//    b_epiz = fChain->GetBranch("EPiz");
//    b_kminmom = fChain->GetBranch("minKMinmom");
//    b_kmaxmom = fChain->GetBranch("minKMom");
//    b_mm2nc = fChain->GetBranch("mm2nc");
// */
//    b_tnu = fChain->GetBranch("tnu");
//    b_fnu = fChain->GetBranch("fnu");
//    //   b_eneu = fChain->GetBranch("ENeu");

//    b_mm2 = fChain->GetBranch("mm2");

//    b_mm2fit = fChain->GetBranch("mm2fit");
//    /*b_allksm0 = fChain->GetBranch("allksm0");
//       b_allksp = fChain->GetBranch("allksp");
//    b_allksmc = fChain->GetBranch("allksmc");
//    b_allchkp = fChain->GetBranch("allchkp");
//    b_allchkmc = fChain->GetBranch("allchkmc");
//    b_m0ks = fChain->GetBranch("m0ks");
//    b_pks = fChain->GetBranch("pks");
//    b_pksmc = fChain->GetBranch("pksmc");
//    b_ntkl = fChain->GetBranch("ntkl");
//    b_tklp = fChain->GetBranch("tklp");
//    b_tklth = fChain->GetBranch("tklth");
//    b_tklph = fChain->GetBranch("tklph");
//    b_tklisol = fChain->GetBranch("tklisol");
//    b_ntks = fChain->GetBranch("ntks");
//    b_tksp = fChain->GetBranch("tksp");
//    b_tksth = fChain->GetBranch("tksth");
//    b_tksph = fChain->GetBranch("tksph");
//    b_tksdec = fChain->GetBranch("tksdec");
//    b_ntchk = fChain->GetBranch("ntchk");
//    b_tchkp = fChain->GetBranch("tchkp");
//    b_tchkth = fChain->GetBranch("tchkth");
//    b_tchkph = fChain->GetBranch("tchkph");
//    b_nklres = fChain->GetBranch("nklres");
//    b_klresth = fChain->GetBranch("klresth");
//    b_klresph = fChain->GetBranch("klresph");
//    b_klid = fChain->GetBranch("klid");
//    b_klcone = fChain->GetBranch("klcone");
//    b_emckl = fChain->GetBranch("emckl");
//    b_emckl0 = fChain->GetBranch("emckl0");
//    b_emckl22 = fChain->GetBranch("emckl22");
//    b_mxks = fChain->GetBranch("mxks");
//    b_mm2ks = fChain->GetBranch("mm2ks");
//    b_mxksfit = fChain->GetBranch("mxksfit");
//    b_mm2misk = fChain->GetBranch("mm2misk");
//    b_mxmisk = fChain->GetBranch("mxmisk");
//    b_mxmiskfit = fChain->GetBranch("mxmiskfit");
//    b_mm2chk = fChain->GetBranch("mm2chk");
//    b_mxchk = fChain->GetBranch("mxchk");
//    b_mxchkfit = fChain->GetBranch("mxchkfit");*/
   
//    b_totweight = fChain->GetBranch("totweight");
//    b_totweightNutMult = fChain->GetBranch("totweightNutMult");
//    b_totweightTrkMult = fChain->GetBranch("totweightTrkMult");
//    return kTRUE;
// }
