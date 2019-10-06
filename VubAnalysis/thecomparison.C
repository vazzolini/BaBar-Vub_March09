#define thecomparison_cxx
#include "thecomparison.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "../RecoilAnalysis/mesData.hh" 
#include <TVector2.h>
thecomparison::thecomparison(char *var, double mi, double ma, int b)
{
   themin = mi;
   themax = ma;
   thebins = b;
   thevar = var ;
   Dvar = new recoilDSys("ddecay.table",0,2);
   Bsem = new recoilDSys(0);
   
}
 
void thecomparison::Loop(int nevents, int cat, double shift, double smear, int bsel, int multcat, int seed,int sys)
{
//   In a Root session, you can do:
//      Root > .L thecomparison.C
//      Root > thecomparison t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
  //      Root > t.Show(16);     // Read and show values of entry 16
  //      Root > t.Loop();       // Loop on all entries
  //
  
  //     This is the loop skeleton
  //       To read only selected branches, Insert statements like:
  // METHOD1:
  //    fChain->SetBranchStatus("*",0);  // disable all branches
  //    fChain->SetBranchStatus("branchname",1);  // activate branchname
  // METHOD2: replace line
  //    fChain->GetEntry(i);  // read all branches
  //by  b_branchname->GetEntry(i); //read only this branch

  //  gROOT->cd();
  
  char name[100], preine[100], preich[100], prefix[100];
  char name2[100];
  char le[100];
  
  //id = Mx category == Mx bin
  int id;
  int idflav;
  
  int ich, ine;

  double cutflav = 0;    
  double thecat;
  int number;
  int group;        
  
  multipl = multcat;

  SHIFTNEUT = shift;
  SIGMANEUT = smear;
  // restore B0/B+ flag
  int isbch =  bsel>=0? bsel: -(1+bsel);
      
  if (fChain == 0) return;
  
  Int_t nentries = Int_t(fChain->GetEntries());
  
  if( nentries > nevents) nentries = nevents;
  
  cout <<  nentries << endl;
  
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
    int ibch = 1;                     
    
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
    bool NPI0 = nnpi0<2;

    bool RUN = 1;
       //bool RUN = eneu<0.15*exhad;
    bool NLE = nle > 0;
    bool NEL = nel;
    bool NMU = nmu;
    bool NLE1 = (nle == 1 &&pcms>1.) ;
    bool MM2 =  mm2 < .5;
    bool CH = ch == 0;
    bool KSELE = ksele >0;
    int FLAV = !(TMath::Abs(brecocharge)!=0 && (lcharge + brecoflav)!=0);
    bool BCH = TMath::Abs(brecocharge) == isbch;     
    bool IPUR = intpur>0;
    bool CAT = multcat == thecat;

    bool SEED = 1;
    if (seed == 2) SEED = (mode>13000 && mode<14000);
    if (seed == 3) SEED = (mode>14000 && mode<16000);
    if (seed == 4) SEED = (mode>12000 && mode<13000);
    if (seed == 5) SEED = (mode>11000 && mode<12000);    

    if (isbch == 2) BCH = 1;
    if (multcat == 7) CAT = 1;

    double myvar = -1;
    if (thevar == "mxhadfit") {myvar = mxhadfit;}
    if (thevar == "csiCiuc") {myvar = csiCiuc;}
    if (thevar == "xCiuc") {myvar = xCiuc;}
    if (thevar == "wCiuc") {myvar = wCiuc;}
    if (thevar == "EwPwfit") {myvar = EwPwfit;}
    if (thevar == "q2fit") {myvar = q2fit;}

    if (thevar == "GcsiCiuc") {myvar = GcsiCiuc;}
    if (thevar == "GxCiuc") {myvar = GxCiuc;}
    if (thevar == "GwCiuc") {myvar = GwCiuc;}
    if (thevar == "EwPwG") {myvar = EwPwG;}
    if (thevar == "q2Gen") {myvar = q2Gen;}

    if (thevar == "rescsiCiuc") {myvar = csiCiuc-GcsiCiuc;}
    if (thevar == "resxCiuc")   {myvar = xCiuc-GxCiuc;}
    if (thevar == "reswCiuc")   {myvar = wCiuc-GwCiuc;}
    if (thevar == "resEwPwfit") {myvar = EwPwfit-EwPwG;}
    if (thevar == "resq2fit")   {myvar = q2fit-q2Gen;}
    if (thevar == "mxhad") {myvar = mxhad;}
    if (thevar == "mm2") {myvar = mm2; MM2=1; }
    if (thevar == "efneu") {myvar =exhad>0 ? eneu/exhad:-1;}
    if (thevar == "eneu") {myvar =eneu;}
    if (thevar == "epi0") {myvar =epiz;}
    if (thevar == "esneu") {myvar =nneu>0 ? eneu/nneu:-1;}
    if (thevar == "etrk") {myvar = exhad-eneu;}
    if (thevar == "kmin") {if(nkp<=0)continue;myvar = kminmom;}
    if (thevar == "kmax") {if(nkp<=0)continue;myvar = kmaxmom;}
    if (thevar == "nneu") myvar = nneu;
    if (thevar == "pnu") myvar = pnu;
    if (thevar == "dxy")  {myvar = sqrt(dx*dx+dy*dy);if(myvar>1.)continue;}
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
    }
    if (thevar == "nneu80_160") myvar = nneu80_160;
    if (thevar == "nneu160_320") myvar = nneu160_320;
    if (thevar == "nneu320") myvar = nneu - nneu80_160;
    if (thevar == "nchg") myvar = nchg;
    if (thevar == "pcms") {myvar = pcms; PCMS = 1;}
    if (thevar == "nle") {myvar = nle; NLE = 1; NLE1 = 1;}
    if (thevar == "qtot") {myvar = ch; CH = 1;}
    if (thevar == "npi0") {myvar = npi0; }
    if (thevar == "nkp") {myvar = nkp; KSELE = nks > 0; }
    if (thevar == "nks") {myvar = nks; KSELE = nkp > 0; }
    if (thevar == "intpur") {myvar = intpur; }
    if (thevar == "pur") {myvar = pur; }
    if (thevar == "mode") {myvar = mode; }
    if (!((myvar>0) || (myvar<0) || (myvar == 0))) {
      myvar = -999.;
    } 
    totweight = 1;

    if(cat!=4 && sys==1){
      int dImode;
      //      cout<<" "<<vxbtyp<<" "<<vub<<" "<<GfDpi<<" "<<GfDk<<" "<<endl;
      totweight *= getBsysweight(vxbtyp,vub);//Bdec weighting
      totweight *= getDsysweight(GfDpi,GfDk,GfDks,GfDpiz,GfDlep,dImode,vub); //Ddec weighting

    } else if( cat!=4 && sys==2&&vub&&fkplus>0){
      totweight *= (1.5*FermiWeight(fkplus,-0.15,0.));// 1.5 to account for the possible difference in yields
      
    }
    if(PCMS && NLE && FLAV && BCH && IPUR && CAT && SEED && RUN) {
      
      sprintf(le, "h");
      if(KSELE) sprintf(le, "d");
      if(SIGMANEUT>0.0) myvar = smeargauss(myvar, SHIFTNEUT, SIGMANEUT); 
      
      // assign flavor category
      cutflav = 5;
      if(TMath::Abs(brecocharge)) cutflav = 3;
      if(TMath::Abs(brecocharge)==0 && flav==0) cutflav = 4;	  
      
      id = hist(myvar);
      idflav = hist(myvar) + cutflav * 10000;    
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
  char name[100], title[100], preine[100], preich[100], prefix[100], number[100];
  Int_t nbins = 13;
  Int_t nbins2 = 10;

  int ine; int ich; int lo; 
 
  sprintf(name, "h400000");  sprintf(title, "%s%s", thevar, " data events after all cuts: enriched");  h = new TH1D(name, title, thebins, themin, themax );  h->Sumw2();
  sprintf(name, "h500000");  sprintf(title, "%s%s", thevar, " MC events  after all cuts: enriched");  h = new TH1D(name, title, thebins, themin, themax  );  h->Sumw2();
  sprintf(name, "h94206");  sprintf(title, "mes data after all cuts: enriched");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  sprintf(name, "h95206");  sprintf(title, "mes MC after all cuts: enriched");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  sprintf(name, "h1400000");  sprintf(title, "%s%s", thevar, " data events after lepton cuts: enriched");  h = new TH1D(name, title, thebins, themin, themax );  h->Sumw2();
  sprintf(name, "h1500000");  sprintf(title, "%s%s", thevar, " MC events after lepton cuts: enriched");  h = new TH1D(name, title, thebins, themin, themax  );  h->Sumw2();
  sprintf(name, "h194206");  sprintf(title, "mes data after lepton cuts: enriched");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  sprintf(name, "h195206");  sprintf(title, "mes MC after lepton cuts: enriched");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();

  sprintf(name, "d400000");  sprintf(title, "%s%s", thevar, " data events after all cuts: depleted");  h = new TH1D(name, title, thebins, themin, themax );  h->Sumw2();
  sprintf(name, "d500000");  sprintf(title, "%s%s", thevar, " MC events  after all cuts: depleted");  h = new TH1D(name, title, thebins, themin, themax  );  h->Sumw2();
  sprintf(name, "d94206");  sprintf(title, "mes data after all cuts: depleted");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  sprintf(name, "d95206");  sprintf(title, "mes MC after all cuts: depleted");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  sprintf(name, "d1400000");  sprintf(title, "%s%s", thevar, " data events after lepton cuts: depleted");  h = new TH1D(name, title, thebins, themin, themax );  h->Sumw2();
  sprintf(name, "d1500000");  sprintf(title, "%s%s", thevar, " MC events after lepton cuts: depleted");  h = new TH1D(name, title, thebins, themin, themax  );  h->Sumw2();
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

void thecomparison::Fitmes(int cat, int cut){

  // fit to mes distribution and fill of the Mx plots...

  //  gROOT->cd();

  int group = 200 + cat * 1000;
   
  char nameps[100], preine[100];
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
  int ia;
  int is;
  int ich;
  int ine;
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
  if (y == 1) sprintf(le, "d");

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

void 
thecomparison::sighisto(double&signal, double& signalErr, TH1D *histo, double &resmean, double &ressigma, double &resalpha, double &resn, int fixpar, double mean, double sigma, double alpha,  double n, double argus){

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

   c6 = new TCanvas("c6", "c6", 300.,   0., 400,800);  
   // -- top left
   fPads[1]= new TPad("pad1", "", 0.00, 0.65, 0.99, 0.99);   fPads[1]->Draw(); 
   fPads[2]= new TPad("pad2", "", 0.00, 0.50, 0.99, 0.64);   fPads[2]->Draw();  
   // -- bottom left
   fPads[3]= new TPad("pad3", "", 0.00, 0.15, 0.99, 0.49);   fPads[3]->Draw(); 
   fPads[4]= new TPad("pad4", "", 0.00, 0.00, 0.99, 0.14);   fPads[4]->Draw(); 
   char name[100], name2[100], preine[100], preich[100], prefix[100], shift[100], smear[100], line[100];
   int h = 0 ; 

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
   double max = 1.4*y1.GetMaximum();
   y1.SetLabelSize(0.07, "Y");
   y1.SetMaximum(max);
   y1.SetMinimum(0.);
   inter1 = y1.Integral();
// y1.SetNormFactor(inter1);
   y1.SetLineColor(kRed);
   y1.Draw();
 
   double inter2;
   inter2 = y2.Integral();
    
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

   double g = chisq(&y1,&y2);
   sprintf(line, "#chi^{2} = %5.4f", g); tl->SetTextSizePixels(50); tl->DrawLatex(0.22, 0.75, line);
  
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
   hratio->SetTitle();
   hratio->Draw();

   fPads[3]->cd(); shrinkPad(0.001, 0.2); 

   setFilledHist(&y3, kBlack, kBlue, 3004);
   double inter3;
   double max = 1.4*y3->GetMaximum();
   y3.SetLabelSize(0.07, "Y");
   y3.SetMaximum(max);
   y3.SetMinimum(0.);
   inter3 = y3.Integral();
// y3.SetNormFactor(inter1);
   y3.SetLineColor(kRed);
   y3.Draw();
 
   double inter4;
   inter4 = y4.Integral();
    
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

   double g = chisq(&y3,&y4);
   sprintf(line, "#chi^{2} = %5.4f", g); tl->SetTextSizePixels(50); tl->DrawLatex(0.22, 0.75, line);
  
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
   hratio2->SetTitle();
   hratio2->Draw();

//    cout << "the chi square is " << g << endl;
//    SHIFTNEUT = SHIFTNEUT *1000;
//    SIGMANEUT = SIGMANEUT *1000;
//    sprintf(shift, "%s%d","-", SHIFTNEUT);
//    sprintf(smear, "%s%d","-", SIGMANEUT);   
  
   sprintf(mycut, "allcuts");   
   if (cut) sprintf(mycut, "leptoncuts");   
   if (norm) {   
      sprintf(name, "%s%s%d%s%s%s", dir.Data(), "/comparisonnorm",multipl,thevar,mycut,".ps");
   }else{
      sprintf(name, "%s%s%d%s%s%s", dir.Data(), "/comparison",multipl,thevar,mycut,".ps");
   }
   c6.SaveAs(name);  
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

   TH1D *hratio1 = new TH1D(*y1); hratio->SetName("hratio1"); hratio->Reset();
   TH1D *hratio2 = new TH1D(*y5); hratio->SetName("hratio2"); hratio->Reset();
   TH1D *hratio3 = new TH1D(*y3); hratio->SetName("hratio3"); hratio->Reset();
   TH1D *hratio4 = new TH1D(*y7); hratio->SetName("hratio4"); hratio->Reset();
   TH1D *hratioratio1 = new TH1D(*hratio1); hratio->SetName("hratioratio1"); hratio->Reset();
   TH1D *hratioratio2 = new TH1D(*hratio3); hratio->SetName("hratioratio2"); hratio->Reset();

  
   hratio1->Divide(y1, y2, 1, 1, "B");
   hratio2->Divide(y5, y6, 1, 1, "B");
   hratio3->Divide(y3, y4, 1, 1, "B");
   hratio4->Divide(y7, y8, 1, 1, "B");
   hratioratio1->Divide(hratio1, hratio2);
   hratioratio2->Divide(hratio3, hratio4);

   c7 = new TCanvas("c7", "c7", 300.,   0., 400,800);  
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
   hratio1->SetMinimum(-0.1.);
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
   hratioratio1->SetTitle();
   hratioratio1->Draw();

   fPads[3]->cd(); shrinkPad(0.001, 0.2); 

//   double max = 1.4*hratio1->GetMaximum();
   hratio3->SetLabelSize(0.07, "Y");
//   hratio1->SetMaximum(max);
//cout << max << endl;
   hratio3->SetMaximum(1.);
   hratio3->SetMinimum(-0.1.);
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
   hratioratio2->SetTitle();
   hratioratio2->Draw();

   sprintf(name, "%s%s%d%s%s", dir.Data(), "/comparisoneff",multipl,thevar,".ps");
   c7.SaveAs(name);  


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
   fChain->SetBranchAddress("bmass",&bmass);
   fChain->SetBranchAddress("bmassfit",&bmassfit);
   fChain->SetBranchAddress("sbox",&sbox);
   fChain->SetBranchAddress("mes",&mes);
   fChain->SetBranchAddress("de",&de);
   fChain->SetBranchAddress("pur",&pur);
   fChain->SetBranchAddress("intpur",&intpur);
   fChain->SetBranchAddress("mode",&mode);
   fChain->SetBranchAddress("nnpi0",&nnpi0);
   fChain->SetBranchAddress("GSem",&GSem);  
   fChain->SetBranchAddress("GfDpi",&GfDpi); 
   fChain->SetBranchAddress("GfDpiz",&GfDpiz);
   fChain->SetBranchAddress("GfDk",&GfDk);  
   fChain->SetBranchAddress("GfDks",&GfDks); 
   fChain->SetBranchAddress("GfDlep",&GfDlep);
   fChain->SetBranchAddress("GfDgam",&GfDgam);
   fChain->SetBranchAddress("brecoflav",&brecoflav);
   fChain->SetBranchAddress("brecocharge",&brecocharge);
   fChain->SetBranchAddress("brecomc",&brecomc);
   fChain->SetBranchAddress("mxhadgen",&mxhadgen);
   fChain->SetBranchAddress("pcmsgen",&pcmsgen);
   fChain->SetBranchAddress("tcmsgen",&tcmsgen);
   fChain->SetBranchAddress("fcmsgen",&fcmsgen);
   fChain->SetBranchAddress("ecmsgen",&ecmsgen);
   fChain->SetBranchAddress("pxhadgen",&pxhadgen);
   fChain->SetBranchAddress("txhadgen",&txhadgen);
   fChain->SetBranchAddress("fxhadgen",&fxhadgen);
   fChain->SetBranchAddress("exhadgen",&exhadgen);
   fChain->SetBranchAddress("kplus",&fkplus);
   fChain->SetBranchAddress("GoodEvent",&GoodEvent);
   fChain->SetBranchAddress("isDupli",&isDupli);
   fChain->SetBranchAddress("ValMap",&ValMap);
   fChain->SetBranchAddress("vub",&vub);
   fChain->SetBranchAddress("vcb",&vcb);
   fChain->SetBranchAddress("vxbtyp",&vxbtyp);
   fChain->SetBranchAddress("other",&other);
   fChain->SetBranchAddress("bgcat",&bgcat);
   fChain->SetBranchAddress("xcharge",&xcharge);
   fChain->SetBranchAddress("pxhad",&pxhad);
   fChain->SetBranchAddress("txhad",&txhad);
   fChain->SetBranchAddress("fxhad",&fxhad);
   fChain->SetBranchAddress("exhad",&exhad);
   fChain->SetBranchAddress("mxhad",&mxhad);
   fChain->SetBranchAddress("gmax",&gmax);
   fChain->SetBranchAddress("mxhadfit",&mxhadfit);
   fChain->SetBranchAddress("EwPwfit",&EwPwfit);
   fChain->SetBranchAddress("q2fit",&q2fit);
   fChain->SetBranchAddress("csiCiuc",&csiCiuc);
   fChain->SetBranchAddress("xCiuc",&xCiuc);
   fChain->SetBranchAddress("wCiuc",&wCiuc);
   fChain->SetBranchAddress("EwPwG",&EwPwG);
   fChain->SetBranchAddress("q2Gen",&q2Gen);
   fChain->SetBranchAddress("GcsiCiuc",&GcsiCiuc);
   fChain->SetBranchAddress("GxCiuc",&GxCiuc);
   fChain->SetBranchAddress("GwCiuc",&GwCiuc);
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
   fChain->SetBranchAddress("nneu80_160",&nneu80_160);
   fChain->SetBranchAddress("nneu160_320",&nneu160_320);
   fChain->SetBranchAddress("nneufromB",&nneufromB);
   fChain->SetBranchAddress("nneufromB80_160",&nneufromB80_160);
   fChain->SetBranchAddress("nneufromB160_320",&nneufromB160_320);
   fChain->SetBranchAddress("nkp",&nkp);
   fChain->SetBranchAddress("nks",&nks);
   fChain->SetBranchAddress("npi0",&npi0);
   fChain->SetBranchAddress("pnu",&pnu);
   fChain->SetBranchAddress("dx",&dx);
   fChain->SetBranchAddress("dy",&dy);
   fChain->SetBranchAddress("dz",&dz);
   fChain->SetBranchAddress("s2dxx",&s2dxx);
   fChain->SetBranchAddress("s2dyy",&s2dyy);
   fChain->SetBranchAddress("s2dzz",&s2dzz);
   fChain->SetBranchAddress("s2dxy",&s2dxy);
   fChain->SetBranchAddress("s2dyz",&s2dyz);
   fChain->SetBranchAddress("s2dxz",&s2dxz);
   fChain->SetBranchAddress("tnu",&tnu);
   fChain->SetBranchAddress("fnu",&fnu);
   fChain->SetBranchAddress("ENeu",&eneu);
   fChain->SetBranchAddress("EPiz",&epiz);
   fChain->SetBranchAddress("MinKMom",&kminmom);
   fChain->SetBranchAddress("MaxKMom",&kmaxmom);
   fChain->SetBranchAddress("mm2",&mm2);
   fChain->SetBranchAddress("mm2nc",&mm2nc);
   fChain->SetBranchAddress("mm2fit",&mm2fit);
   fChain->SetBranchAddress("allksm0",allksm0);
   fChain->SetBranchAddress("allksp",allksp);
   fChain->SetBranchAddress("allksmc",allksmc);
   fChain->SetBranchAddress("allchkp",allchkp);
   fChain->SetBranchAddress("allchkmc",allchkmc);
   fChain->SetBranchAddress("m0ks",&m0ks);
   fChain->SetBranchAddress("pks",&pks);
   fChain->SetBranchAddress("pksmc",&pksmc);
   fChain->SetBranchAddress("ntkl",&ntkl);
   fChain->SetBranchAddress("tklp",&tklp);
   fChain->SetBranchAddress("tklth",&tklth);
   fChain->SetBranchAddress("tklph",&tklph);
   fChain->SetBranchAddress("tklisol",&tklisol);
   fChain->SetBranchAddress("ntks",&ntks);
   fChain->SetBranchAddress("tksp",&tksp);
   fChain->SetBranchAddress("tksth",&tksth);
   fChain->SetBranchAddress("tksph",&tksph);
   fChain->SetBranchAddress("tksdec",&tksdec);
   fChain->SetBranchAddress("ntchk",&ntchk);
   fChain->SetBranchAddress("tchkp",&tchkp);
   fChain->SetBranchAddress("tchkth",&tchkth);
   fChain->SetBranchAddress("tchkph",&tchkph);
   fChain->SetBranchAddress("nklres",&nklres);
   fChain->SetBranchAddress("klresth",&klresth);
   fChain->SetBranchAddress("klresph",&klresph);
   fChain->SetBranchAddress("klid",&klid);
   fChain->SetBranchAddress("klcone",&klcone);
   fChain->SetBranchAddress("emckl",&emckl);
   fChain->SetBranchAddress("emckl0",&emckl0);
   fChain->SetBranchAddress("emckl22",&emckl22);
   fChain->SetBranchAddress("mxks",&mxks);
   fChain->SetBranchAddress("mm2ks",&mm2ks);
   fChain->SetBranchAddress("mxksfit",&mxksfit);
   fChain->SetBranchAddress("mm2misk",&mm2misk);
   fChain->SetBranchAddress("mxmisk",&mxmisk);
   fChain->SetBranchAddress("mxmiskfit",&mxmiskfit);
   fChain->SetBranchAddress("mm2chk",&mm2mchk);
   fChain->SetBranchAddress("mxchk",&mxchk);
   fChain->SetBranchAddress("mxchkfit",&mxchkfit);
   fChain->SetBranchAddress("totweight",&totweight);
   fChain->SetBranchAddress("totweightNutMult",&totweightNutMult);
   fChain->SetBranchAddress("totweightTrkMult",&totweightTrkMult);
   Notify();
}

Bool_t thecomparison::Notify()
{
//   called when loading a new file
//   get branch pointers
   b_run = fChain->GetBranch("run");
   b_lower = fChain->GetBranch("lower");
   b_upper = fChain->GetBranch("upper");
   b_bmass = fChain->GetBranch("bmass");
   b_bmassfit = fChain->GetBranch("bmassfit");
   b_sbox = fChain->GetBranch("sbox");
   b_mes = fChain->GetBranch("mes");
   b_de = fChain->GetBranch("de");
   b_pur = fChain->GetBranch("pur");
   b_intpur = fChain->GetBranch("intpur");
   b_mode = fChain->GetBranch("mode");
   b_nnpi0 = fChain->GetBranch("nnpi0");
   b_brecoflav = fChain->GetBranch("brecoflav");
   b_brecocharge = fChain->GetBranch("brecocharge");
   b_brecomc = fChain->GetBranch("brecomc");
   b_mxhadgen = fChain->GetBranch("mxhadgen");
   b_pcmsgen = fChain->GetBranch("pcmsgen");
   b_tcmsgen = fChain->GetBranch("tcmsgen");
   b_fcmsgen = fChain->GetBranch("fcmsgen");
   b_ecmsgen = fChain->GetBranch("ecmsgen");
   b_pxhadgen = fChain->GetBranch("pxhadgen");
   b_txhadgen = fChain->GetBranch("txhadgen");
   b_fxhadgen = fChain->GetBranch("fxhadgen");
   b_exhadgen = fChain->GetBranch("exhadgen");
   b_kplus = fChain->GetBranch("kplus");
   b_GoodEvent = fChain->GetBranch("GoodEvent");
   b_isDupli = fChain->GetBranch("isDupli");
   b_ValMap = fChain->GetBranch("ValMap");
   b_vub = fChain->GetBranch("vub");
   b_vcb = fChain->GetBranch("vcb");
   b_vxbtyp = fChain->GetBranch("vxbtyp");
   b_other = fChain->GetBranch("other");
   b_bgcat = fChain->GetBranch("bgcat");
   b_xcharge = fChain->GetBranch("xcharge");
   b_pxhad = fChain->GetBranch("pxhad");
   b_txhad = fChain->GetBranch("txhad");
   b_fxhad = fChain->GetBranch("fxhad");
   b_exhad = fChain->GetBranch("exhad");
   b_mxhad = fChain->GetBranch("mxhad");
   b_gmax = fChain->GetBranch("gmax");
   b_mxhadfit = fChain->GetBranch("mxhadfit");
   b_q2fit = fChain->GetBranch("q2fit");
   b_EwPwfit = fChain->GetBranch("EwPwfit");
   b_csiCiuc = fChain->GetBranch("csiCiuc");
   b_xCiuc = fChain->GetBranch("xCiuc");
   b_wCiuc = fChain->GetBranch("wCiuc");
   b_q2Gen = fChain->GetBranch("q2Gen");
   b_EwPwG = fChain->GetBranch("EwPwG");
   b_GcsiCiuc = fChain->GetBranch("GcsiCiuc");
   b_GxCiuc = fChain->GetBranch("GxCiuc");
   b_GwCiuc = fChain->GetBranch("GwCiuc");
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
   b_nneu80_160 = fChain->GetBranch("nneu80_160");
   b_nneu160_320 = fChain->GetBranch("nneu160_320");
   b_nneufromB = fChain->GetBranch("nneufromB");
   b_nneufromB80_160 = fChain->GetBranch("nneufromB80_160");
   b_nneufromB160_320 = fChain->GetBranch("nneufromB160_320");
   b_nkp = fChain->GetBranch("nkp");
   b_nks = fChain->GetBranch("nks");
   b_npi0 = fChain->GetBranch("npi0");
   b_pnu = fChain->GetBranch("pnu");
   b_dx = fChain->GetBranch("dx");
   b_dy = fChain->GetBranch("dy");
   b_dz = fChain->GetBranch("dz");
   b_s2dxx = fChain->GetBranch("s2dxx");
   b_s2dyy = fChain->GetBranch("s2dyy");
   b_s2dzz = fChain->GetBranch("s2dzz");
   b_s2dxy = fChain->GetBranch("s2dxy");
   b_s2dyz = fChain->GetBranch("s2dyz");
   b_s2dxz = fChain->GetBranch("s2dxz");
   b_tnu = fChain->GetBranch("tnu");
   b_fnu = fChain->GetBranch("fnu");
   b_eneu = fChain->GetBranch("ENeu");
   b_epiz = fChain->GetBranch("EPiz");
   b_kminmom = fChain->GetBranch("minKMinmom");
   b_kmaxmom = fChain->GetBranch("minKMom");
   b_mm2 = fChain->GetBranch("mm2");
   b_mm2nc = fChain->GetBranch("mm2nc");
   b_mm2fit = fChain->GetBranch("mm2fit");
   b_allksm0 = fChain->GetBranch("allksm0");
   b_allksp = fChain->GetBranch("allksp");
   b_allksmc = fChain->GetBranch("allksmc");
   b_allchkp = fChain->GetBranch("allchkp");
   b_allchkmc = fChain->GetBranch("allchkmc");
   b_m0ks = fChain->GetBranch("m0ks");
   b_pks = fChain->GetBranch("pks");
   b_pksmc = fChain->GetBranch("pksmc");
   b_ntkl = fChain->GetBranch("ntkl");
   b_tklp = fChain->GetBranch("tklp");
   b_tklth = fChain->GetBranch("tklth");
   b_tklph = fChain->GetBranch("tklph");
   b_tklisol = fChain->GetBranch("tklisol");
   b_ntks = fChain->GetBranch("ntks");
   b_tksp = fChain->GetBranch("tksp");
   b_tksth = fChain->GetBranch("tksth");
   b_tksph = fChain->GetBranch("tksph");
   b_tksdec = fChain->GetBranch("tksdec");
   b_ntchk = fChain->GetBranch("ntchk");
   b_tchkp = fChain->GetBranch("tchkp");
   b_tchkth = fChain->GetBranch("tchkth");
   b_tchkph = fChain->GetBranch("tchkph");
   b_nklres = fChain->GetBranch("nklres");
   b_klresth = fChain->GetBranch("klresth");
   b_klresph = fChain->GetBranch("klresph");
   b_klid = fChain->GetBranch("klid");
   b_klcone = fChain->GetBranch("klcone");
   b_emckl = fChain->GetBranch("emckl");
   b_emckl0 = fChain->GetBranch("emckl0");
   b_emckl22 = fChain->GetBranch("emckl22");
   b_mxks = fChain->GetBranch("mxks");
   b_mm2ks = fChain->GetBranch("mm2ks");
   b_mxksfit = fChain->GetBranch("mxksfit");
   b_mm2misk = fChain->GetBranch("mm2misk");
   b_mxmisk = fChain->GetBranch("mxmisk");
   b_mxmiskfit = fChain->GetBranch("mxmiskfit");
   b_mm2chk = fChain->GetBranch("mm2chk");
   b_mxchk = fChain->GetBranch("mxchk");
   b_mxchkfit = fChain->GetBranch("totweight");
   b_mxchkfit = fChain->GetBranch("totweightNutMult");
   b_mxchkfit = fChain->GetBranch("totweightTrkMult");
   return kTRUE;
}
