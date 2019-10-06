#include "CompTool.hh"
#include "TH2.h"
#include "TH1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"
#include "../RecoilAnalysis/recoilAnalysis.hh" 
#include "../RecoilAnalysis/mesData.hh" 
#include <TVector2.h>
#include <TRandom.h>
#include <fstream.h>
#include <iostream.h>
#include <strstream.h>
#include <iomanip.h>
#include "recoilDSys.hh"

void CompTool::Loop(int nevents, int cat, double shift, double smear, int bsel, int multcat, int seed)
{
  gROOT->cd();
  
  char name[100];  char le[100];
  
  //id = Mx category == Mx bin
  int id, idflav, ich, ine;
  double thecat;  int group;        
  
  multipl = multcat;  SHIFTNEUT = shift;  SIGMANEUT = smear;
  // restore B0/B+ flag
  int isbch =  bsel>=0? bsel: -(1+bsel);
      
  if (fChain == 0) return;
  Int_t nentries = Int_t(fChain->GetEntries());
  if( nentries > nevents) nentries = nevents;
  cout <<  nentries << endl;  Int_t nbytes = 0, nb = 0;  double w; 

  for (Int_t jentry=0; jentry<nentries;jentry++) {
    //    Int_t ientry = LoadTree(jentry); //in case of a TChain, ientry is the entry number in the current file
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // vub - vcb selection
    if(bsel<0)
      if( (!vub && cat==4) || (!vcb && cat==5) )continue;
    
    // if (Cut(ientry) < 0) continue;
    int flav =  lcharge + brecoflav; // charge correlation
    bool ksele = nkp + nks;          // fit on the depleted sample?
    //bool ksele = nkp>1;          // fit on the depleted sample?
    int ch = xcharge + brecocharge;  // total charge
    double cutflav = 0;              // flavor category (3 = charged B, 4 = neutral B OS, 5 = neutral B SS)  
    
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
    
    bool PCMS =  pcms>CUTPCMS;
    //    bool NPI0 = nnpi0<2;

    bool NLE = nle > 0;
    //    bool NEL = nel;
    //    bool NMU = nmu;
    bool NLE1 = nle == 1;
    bool MM2 =  mm2 < .5;  bool CH = ch == 0;    bool KSELE = ksele >0;
    int FLAV = !(TMath::Abs(brecocharge)!=0 && (lcharge + brecoflav)!=0);
    bool BCH = TMath::Abs(brecocharge) == isbch;     
    bool IPUR = intpur>INTPURITY;    bool CAT = multcat == thecat;

    bool SEED = 1;
    if (seed == 2) SEED = (mode>13000 && mode<14000);
    if (seed == 3) SEED = (mode>14000 && mode<16000);
    if (seed == 4) SEED = (mode>12000 && mode<13000);
    if (seed == 5) SEED = (mode>11000 && mode<12000);    
    char MYVar[100]; char V1[100]; char  V2[100]; char  V3[100];
    char V4[100];  char  V5[100];  char  V6[100]; char  V7[100]; 
    char V8[100];  char  V9[100];  char  V16[100];
    char V10[100]; char  V11[100]; char  V12[100]; 
    char V13[100]; char  V14[100]; char  V15[100];
    sprintf(MYVar,"%s",thevar);
    sprintf(V1,"mxhad");   sprintf(V2,"mm2");
    sprintf(V3,"efneu");   sprintf(V4,"eneu");
    sprintf(V5,"epi0");    sprintf(V6,"esneu");
    sprintf(V7,"etrk");    sprintf(V8,"kmin");
    sprintf(V9,"kmax");    sprintf(V10,"nneu");
    sprintf(V11,"pnu");    sprintf(V12,"dxy");
    sprintf(V13,"sxy");    sprintf(V14,"d3d");
    sprintf(V15,"s3d");    sprintf(V16,"qtot");
    w = 1;  //calculate reweightings
    w = getBsysweight(Gvxbtyp,vub);                        //Bdec weighting
    w *= getDsysweight(GfDpi,GfDk,GfDks,GfDpiz,GfDlep,dImode,vub); //Ddec weighting
    if (isbch == 2) BCH = 1;
    if (multcat == 7) CAT = 1;
    double myvar = mxhadfit;
    if (!strcmp(MYVar,V1)) {myvar = mxhad;}
    if (!strcmp(MYVar,V2)) {myvar = mm2; MM2 = 1;}
    if (!strcmp(MYVar,V3)) {myvar =exhad>0 ? eneu/exhad:-1;}
    if (!strcmp(MYVar,V4)) {myvar =eneu;}
    if (!strcmp(MYVar,V5)) {myvar =epiz;}
    if (!strcmp(MYVar,V6)) {myvar =nneu>0 ? eneu/nneu:-1;}
    if (!strcmp(MYVar,V7)) {myvar = exhad-eneu;}
    if (!strcmp(MYVar,V8)) {if(nkp<=0)continue;myvar = kminmom;}
    if (!strcmp(MYVar,V9)) {if(nkp<=0)continue;myvar = kmaxmom;}
    if (!strcmp(MYVar,V10)) myvar = nneu;
    if (!strcmp(MYVar,V11)) myvar = pnu;
    if (!strcmp(MYVar,V16)) {myvar = xcharge + brecocharge; CH = 1;}
    if (!strcmp(MYVar,V12))  {myvar = sqrt(dx*dx+dy*dy);if(myvar>1.)continue;}
    if (!strcmp(MYVar,V13))  {
      myvar = dx*dx+dy*dy;
      if(myvar>1. || myvar<= 0)continue;
      myvar=myvar/sqrt(s2dxx*dx*dx+s2dyy*dy*dy+2*s2dxy*dx*dy);
    }
    if (!strcmp(MYVar,V14)) {myvar = sqrt(dx*dx+dy*dy+dz*dz);if(myvar>1.)continue;}
    if (!strcmp(MYVar,V15))  {
      myvar = dx*dx+dy*dy+dz*dz;
      if(myvar>1. || myvar<= 0)continue;
      myvar=myvar/sqrt(s2dxx*dx*dx+s2dyy*dy*dy+2*s2dxy*dx*dy+s2dzz*dz*dz+2*s2dxz*dx*dz+2*s2dyz*dy*dz);
    }
    if (!((myvar>0) || (myvar<0) || (myvar == 0))) {
      myvar = -999.;
    }   
    //if(brecocharge==0) totweight = 1.165;
    totweight *= w;
    if(PCMS && NLE && FLAV && BCH && IPUR && CAT && SEED) {
      
      sprintf(le, "h");
      if(KSELE) sprintf(le, "d");
      //      if(SIGMANEUT>0.0) myvar = smeargauss(myvar, SHIFTNEUT, SIGMANEUT); 
      
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
    }
    
    //    if(PCMS && NLE && !flav && BCH && IPUR && CAT && SEED){  
    if(PCMS && NLE && FLAV && BCH && IPUR && CAT && SEED){  
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
  cout<< "Leaving Loop"<<endl;
}

void CompTool::Bookhist()
{
  gROOT->cd();
  TH1 *h;  char name[100], title[100], number[100]; int lo; 
  for(int iflav=3; iflav<7; iflav++) {
    sprintf(name,"%s%d", "h",400000+iflav*100);  sprintf(title, "%s%s", thevar, " data events after all cuts: enriched");  h = new TH1D(name, title, thebins, themin, themax );  h->Sumw2();
    sprintf(name,"%s%d", "h",500000+iflav*100);  sprintf(title, "%s%s", thevar, " MC events  after all cuts: enriched");  h = new TH1D(name, title, thebins, themin, themax  );  h->Sumw2();
    sprintf(name,"%s%d", "h",1400000+iflav*100);  sprintf(title, "%s%s", thevar, " data events after lepton cuts: enriched");  h = new TH1D(name, title, thebins, themin, themax );  h->Sumw2();
    sprintf(name,"%s%d", "h",1500000+iflav*100);  sprintf(title, "%s%s", thevar, " MC events after lepton cuts: enriched");  h = new TH1D(name, title, thebins, themin, themax  );  h->Sumw2();
    sprintf(name,"%s%d", "d",400000+iflav*100);  sprintf(title, "%s%s", thevar, " data events after all cuts: depleted");  h = new TH1D(name, title, thebins, themin, themax );  h->Sumw2();
    sprintf(name,"%s%d", "d",500000+iflav*100);  sprintf(title, "%s%s", thevar, " MC events  after all cuts: depleted");  h = new TH1D(name, title, thebins, themin, themax  );  h->Sumw2();
    sprintf(name,"%s%d", "d",1400000+iflav*100);  sprintf(title, "%s%s", thevar, " data events after lepton cuts: depleted");  h = new TH1D(name, title, thebins, themin, themax );  h->Sumw2();
    sprintf(name,"%s%d", "d",1500000+iflav*100);  sprintf(title, "%s%s", thevar, " MC events after lepton cuts: depleted");  h = new TH1D(name, title, thebins, themin, themax  );  h->Sumw2();
  }
  sprintf(name, "h94206");  sprintf(title, "mes data after all cuts: enriched");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  sprintf(name, "h95206");  sprintf(title, "mes MC after all cuts: enriched");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  sprintf(name, "h194206");  sprintf(title, "mes data after lepton cuts: enriched");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  sprintf(name, "h195206");  sprintf(title, "mes MC after lepton cuts: enriched");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  sprintf(name, "d94206");  sprintf(title, "mes data after all cuts: depleted");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  sprintf(name, "d95206");  sprintf(title, "mes MC after all cuts: depleted");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
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

void CompTool::Fitmes(int cat, int cut, TString flag){

  // fit to mes distribution and fill of the Mx plots...

  gROOT->cd();
  cout<<"Entered Fitmes function"<<endl;
  int group = 200 + cat * 1000;
   
  char name[100], namebch[100], namedi[100], nameb0os[100], nameb0ss[100], le[100];
  double sigs[1000], errsigs[1000], tempbin, tempbinchb, tempbinb0os;
  double tempbinb0ss, temperr, temperrchb, temperrb0os, temperrb0ss, chid = 0.174;
  int is, thegroup = group;      
  double foomean, foosigma, fooalpha, foon, usedmean, usedsigma, usedalpha, usedn;
  int i;  char dumname[100];  char dummydummy[300];

  int addcut = 0; 
  if(cut) addcut = 1;
  int const theb = thebins + 1;

  TPad * duPads[theb];
  c8 = new TCanvas("c8", "c8", 300.,0., 1500,1500);

  for (int y=0; y<2; y++) {
    sprintf(le, "h");
    if (y == 1) sprintf(le, "d");
    
    //extracting the signal fit parameters 
    sprintf(name, "%s%d",le,group+90006+addcut*100000);  
    cout << ((TH1D*)gDirectory->Get(name))->Integral() << endl;
    double dummy1,dummy2;
    sighisto(dummy1,dummy2,(TH1D*)gDirectory->Get(name),usedmean,usedsigma,usedalpha,usedn,1,-11111111.,-1111111.,-1111111.,-1111111.,-1111111.);
    cout << "mes result for: " << name << " is, MEAN " << usedmean << " SIGMA " << usedsigma << " ALPHA " << usedalpha << " N " << usedn << endl;
    c8->Clear(); c8->cd();
    for (int ith = 1;  ith < theb; ith++) { 
      sprintf(dumname,"%s%d","padDum",ith);
      if(ith<6) {
	duPads[ith]= new TPad(dumname, "", 0.00+0.2*(ith-1), 0.8, 0.2+0.2*(ith-1), 0.99);   duPads[ith]->Draw();  shrinkPad(0.001, 0.2);
      } else if (ith>5 && ith<11) {
	duPads[ith]= new TPad(dumname, "", 0.00+0.2*(ith-6), 0.6, 0.2+0.2*(ith-6), 0.79);   duPads[ith]->Draw();  shrinkPad(0.001, 0.2);
      } else if (ith>10 && ith<16) {
	duPads[ith]= new TPad(dumname, "", 0.00+0.2*(ith-11), 0.4, 0.2+0.2*(ith-11), 0.59);   duPads[ith]->Draw();  shrinkPad(0.001, 0.2);
      } else if (ith>15 && ith<21) {
	duPads[ith]= new TPad(dumname, "", 0.00+0.2*(ith-16), 0.2, 0.2+0.2*(ith-16), 0.39);   duPads[ith]->Draw();  shrinkPad(0.001, 0.2);
      } else if (ith>20 && ith<26) {
	duPads[ith]= new TPad(dumname, "", 0.00+0.2*(ith-21), 0.0, 0.2+0.2*(ith-21), 0.19);   duPads[ith]->Draw();  shrinkPad(0.001, 0.2);
      }
    }
    c8->Update(); 
    
    group = thegroup;
    for(int iflav =3; iflav<6; iflav++) {
      //      cout<<iflav<<endl;
      for (i = 1;  i < theb; i++) { 
	c8->cd();
	duPads[i]->cd(); gPad->Update();
	is = i;
	sprintf(namedi, "%s%d",le,group+is+addcut*100000+iflav*10000);  
	//	cout<<namedi<<endl;
	sighisto(sigs[is-1],errsigs[is-1],(TH1D*)gDirectory->Get(namedi),foomean,foosigma,fooalpha,foon,1,usedmean,usedsigma,usedalpha,usedn,-1111111);
	((TH1D*)gDirectory->Get(namedi))->Draw();
	int title = cat * 100000 + addcut * 1000000 + iflav*100;
	sprintf(name, "%s%d",le,title);  
	((TH1D*)gDirectory->Get(name))->SetBinContent(is, sigs[is-1]);
	((TH1D*)gDirectory->Get(name))->SetBinError(is, errsigs[is-1]);
      }
      //      sprintf(dummydummy,"%s%s%s%d%s",name,"_",flag.Data(),iflav,".eps");    c8->Print(dummydummy);  
    }
    
    sprintf(name,     "%s%d",le, cat * 100000 + addcut * 1000000 + 600);  
    sprintf(namebch,  "%s%d",le, cat * 100000 + addcut * 1000000 + 300);  
    sprintf(nameb0os, "%s%d",le, cat * 100000 + addcut * 1000000 + 400);  
    sprintf(nameb0ss, "%s%d",le, cat * 100000 + addcut * 1000000 + 500);  
    for (i = 1;  i < theb; i++) { 
      tempbinchb = ((TH1D*)gDirectory->Get(namebch))->GetBinContent(i);
      tempbinb0os = ((TH1D*)gDirectory->Get(nameb0os))->GetBinContent(i);
      tempbinb0ss = ((TH1D*)gDirectory->Get(nameb0ss))->GetBinContent(i);
      temperrchb = ((TH1D*)gDirectory->Get(namebch))->GetBinError(i);
      temperrb0os = ((TH1D*)gDirectory->Get(nameb0os))->GetBinError(i);
      temperrb0ss = ((TH1D*)gDirectory->Get(nameb0ss))->GetBinError(i);
      tempbin = tempbinchb + ((1-chid)/(1-2*chid)) * tempbinb0os - (chid/(1-2*chid)) * tempbinb0ss;
      temperr = sqrt(temperrchb*temperrchb + ((1-chid)/(1-2*chid)) * ((1-chid)/(1-2*chid)) * tempbinb0os + (chid/(1-2*chid)) * (chid/(1-2*chid))* tempbinb0ss);
      ((TH1D*)gDirectory->Get(name))->SetBinContent(i, tempbin);
      ((TH1D*)gDirectory->Get(name))->SetBinError(i, temperr); 
    }    
    cout<<"name "<<name<<endl;
  }
}

void 
CompTool::sighisto(double&signal, double& signalErr, TH1D *histo, double &resmean, double &ressigma, double &resalpha, double &resn, int fixpar, double mean, double sigma, double alpha,  double n, double argus){

  mesData *themes;
  double thesigma = -1111111;  double themean = -1111111;  double theargus = -1111111;
  double thealpha = -1111111;  double then = -1111111;
  if (fixpar){
    themean = mean;    thesigma = sigma;    thealpha = alpha;    then = n;
    themes = b.vubMes(histo, resmean, ressigma, resalpha, resn, 1, 1, themean, thesigma, thealpha, then, theargus);
  } else {
    themes = b.vubMes(histo, resmean, ressigma, resalpha, resn, 1, 1, themean, thesigma, thealpha, then, theargus);
  }
  signal=themes->theSig(); signalErr=themes->theErrSig(); 

}

void CompTool::overlap(int cut, int norm, TString dir, TString flag){
  
  char mycut[100];  char name[200], line[100];
  int addcut = 0; 
  if(cut) addcut = 1;
  tl.SetNDC(kTRUE);
  c6 = new TCanvas("c6", "c6", 300.,   0., 400,800);  
  
  sprintf(name,"%s%d" ,"h", 400600 + addcut * 1000000);
  TH1D y1;  ((TH1D*)gDirectory->Get(name))->Copy(y1);
  sprintf(name,"%s%d" ,"h", 500600 + addcut * 1000000);
  TH1D y2;  ((TH1D*)gDirectory->Get(name))->Copy(y2);
  sprintf(name,"%s%d" ,"d", 400600 + addcut * 1000000);
  TH1D y3;  ((TH1D*)gDirectory->Get(name))->Copy(y3);
  sprintf(name,"%s%d" ,"d", 500600 + addcut * 1000000);
  TH1D y4;  ((TH1D*)gDirectory->Get(name))->Copy(y4);
  c6->Clear(); c6->cd();
  // -- top left
  fPads[1]= new TPad("pad1", "", 0.00, 0.65, 0.99, 0.99);   fPads[1]->Draw(); 
  fPads[2]= new TPad("pad2", "", 0.00, 0.50, 0.99, 0.64);   fPads[2]->Draw();  
  // -- bottom left
  fPads[3]= new TPad("pad3", "", 0.00, 0.15, 0.99, 0.49);   fPads[3]->Draw(); 
  fPads[4]= new TPad("pad4", "", 0.00, 0.00, 0.99, 0.14);   fPads[4]->Draw(); 
  fPads[1]->cd(); shrinkPad(0.001, 0.2); 

  b.setFilledHist(&y1, kBlack, kBlue, 3004);
  double inter1;  double max = 1.4*y1.GetMaximum();
  y1.SetLabelSize(0.07, "Y");  y1.SetMaximum(max);  y1.SetMinimum(0.);
  inter1 = y1.Integral();  // y1.SetNormFactor(inter1);
  y1.SetLineColor(kRed);  y1.Draw();
 
  double inter2;  inter2 = y2.Integral();
  
  if(cut) {
    intdataen = inter1;
    intMCen = inter2;
  }
  if (norm) {
    if(intMCen>0)y2.Scale(intdataen/intMCen);
  }else{
    if(inter2>0)y2.Scale(inter1/inter2);
  }   
  b.setFilledHist(&y2 , kBlack, kRed, 3005);
  y2.Draw("histsame");

  sprintf(mycut, "allcuts");   
  if (cut) sprintf(mycut, "leptoncuts");   

  flag.Append(thevar);  flag.Append(mycut);
  if(norm) flag.Append("norm");

  sprintf(name, "%s%s%s", dir.Data(),"/",flag.Data());
  double g = chisqerr(&y1,&y2,name);
  sprintf(line, "#chi^{2} = %5.4f", g); tl.SetTextSizePixels(50); tl.DrawLatex(0.22, 0.75, line);
  
  fPads[2]->cd(); shrinkPad(0.001, 0.2); 
  shrinkPad(0.4, 0.2, 0.1, 0.001);
  
  gPad->SetGridx(1);  gPad->SetGridy(1);
  TH1D *hratio = new TH1D(y1); hratio->SetName("hratio"); hratio->Reset();
  hratio->Divide(&y1, &y2);
  hratio->SetMinimum(0.5); hratio->SetMaximum(1.5);
  hratio->SetMarkerStyle(24);  hratio->SetNdivisions(504, "Y");
  hratio->SetLabelSize(0.22, "X");  hratio->SetLabelSize(0.17, "Y");
  hratio->SetStats(0);  hratio->SetTitle();  hratio->Draw();
  
  fPads[3]->cd(); shrinkPad(0.001, 0.2); 
  
  b.setFilledHist(&y3, kBlack, kBlue, 3004);
  double inter3;
  max = 1.4*y3.GetMaximum();
  y3.SetLabelSize(0.07, "Y");  y3.SetMaximum(max);   y3.SetMinimum(0.); 
  inter3 = y3.Integral();  y3.SetLineColor(kRed);  y3.Draw();
  
  double inter4;  inter4 = y4.Integral();
  
  if(cut) {
    intdatadepl = inter3;
    intMCdepl = inter4;
  }
  if (norm) {
    y4.Scale(intdatadepl/intMCdepl);
  }else{
    y4.Scale(inter3/inter4);
  }   
  
  b.setFilledHist(&y4 , kBlack, kRed, 3005);
  y4.Draw("histsame");

  sprintf(name, "%s%s%s", dir.Data(),"/",flag.Data());
  g = chisqerr(&y3,&y4,name);
  sprintf(line, "#chi^{2} = %5.4f", g); tl.SetTextSizePixels(50); tl.DrawLatex(0.22, 0.75, line);
  
  fPads[4]->cd(); shrinkPad(0.001, 0.2); 
  shrinkPad(0.4, 0.2, 0.1, 0.001);
  
  gPad->SetGridx(1);  gPad->SetGridy(1);
  TH1D *hratio2 = new TH1D(y4); hratio2->SetName("hratio2"); hratio2->Reset();
  hratio2->Divide(&y3, &y4);
  hratio2->SetMinimum(0.5); hratio2->SetMaximum(1.5);
  hratio2->SetMarkerStyle(24);  hratio2->SetNdivisions(504, "Y");
  hratio2->SetLabelSize(0.22, "X");  hratio2->SetLabelSize(0.17, "Y");
  hratio2->SetStats(0);  hratio2->SetTitle();  hratio2->Draw();
  
  if (norm) {   
    sprintf(name, "%s%s%s%d%s", dir.Data(),"/",flag.Data(),multipl,".ps");
  } else {
    sprintf(name, "%s%s%s%d%s", dir.Data(),"/",flag.Data(),multipl,".ps");
  }
  c6->Print(name);  
  //  delete hratio2; delete hratio;  delete c6;  delete fPads;
}

double CompTool::chisq(TH1 *h1, TH1 *h2, TString outfile) {
  int nbins = h1->GetNbinsX();
  if (nbins != h2->GetNbinsX()) {
    cout << "chi2Test: Number of bins not the same" << endl;
    return -99.;
  }
  char name[200];
  sprintf(name,"%s%s%s%s","chisq_",outfile.Data(),h1->GetName(),h2->GetName());
  ofstream Name(name);
  int constrain =0;
  double df = nbins - 1 - constrain; 
  double chsq(0.), a1(0.), a2(0.);
  for (int i = 1; i <= nbins; ++i) {
    a1 = h1->GetBinContent(i);
    a2 = h2->GetBinContent(i);
    if ((TMath::Abs(a1) < 1.e-6) && (TMath::Abs(a2) < 1.e-6)) {
      //      cout << "Zero entries in bin " << i << endl;
      df -= 1.;
    } else if ((a1 < 0.) || (a2 < 0.)) {
      //      cout << "Below zero entries in bin " << i << endl;
      df -= 1.;
    } else {
      Name << "Adding " << ((a1-a2)*(a1-a2))/(a1+a2)  << " from " << a1 << "  " << a2 << " for xmin =  " << h1->GetBinLowEdge(i) << endl;
      chsq += ((a1-a2)*(a1-a2))/(a1+a2);
    }
  }
  Name.close();
  double gamma = 1. - TMath::Gamma(0.5*df, 0.5*chsq);
  //  chi2 = chsq;
  //  ndof = df;
  //  return gamma;
  return chsq/df;
}


// ----------------------------------------------------------------------
double CompTool::chisqerr(TH1 *h1, TH1 *h2, char * outfile) {
  int nbins = h1->GetNbinsX();
  if (nbins != h2->GetNbinsX()) {
    cout << "chi2Test: Number of bins not the same" << endl;
    return -99.;
  }
  char name[200];
  sprintf(name,"%s%s%s%s",outfile,"chisq_",h1->GetName(),h2->GetName());
  ofstream Name(name);
  int constrain =0;
  double df = nbins - 1 - constrain; 
  double chsq(0.), a1(0.), a2(0.), e1(0.), e2(0.);
  for (int i = 1; i <= nbins; ++i) {
    a1 = h1->GetBinContent(i);
    e1 = h1->GetBinError(i) * h1->GetBinError(i);
    a2 = h2->GetBinContent(i);
    e2 = h2->GetBinError(i) *  h2->GetBinError(i);
    if ((TMath::Abs(a1) < 1.e-8) && (TMath::Abs(a2) < 1.e-8)) {
      //      cout << "Zero entries in bin " << i << endl;
      df -= 1.;
    } else if ((a1 < 0.) || (a2 < 0.)) {
      //      cout << "Below zero entries in bin " << i << endl;
      df -= 1.;
    } else {
      Name << "Adding " << ((a1-a2)*(a1-a2))/(e1+e2)  << " from " << a1 << "+/-" << TMath::Sqrt(e1)
      	   << "  " << a2 << "+/-" << TMath::Sqrt(e2)
	   << " for xmin =  " << h1->GetBinLowEdge(i) << endl;
      chsq += ((a1 - a2) * (a1 - a2)) / (e1 + e2);
    }
  }
  //  double gamma = 1. - TMath::Gamma(0.5*df, 0.5*chsq);
  //  chi2 = chsq;
  //  ndof = df;
  Name <<"Chisquare: "<< chsq/df <<endl;
  Name.close();
  return chsq/df;
}


void CompTool::sameP(TString dir, TString flag){
  
  //  gROOT->cd();
  char name[200];
  sprintf(name,"%s%d","h", 400600);  TH1D yD1; ((TH1D*)gDirectory->Get(name))->Copy(yD1);
  sprintf(name,"%s%d","h", 1400600); TH1D yD2; ((TH1D*)gDirectory->Get(name))->Copy(yD2);
  sprintf(name,"%s%d","d", 400600);  TH1D yD3; ((TH1D*)gDirectory->Get(name))->Copy(yD3);
  sprintf(name,"%s%d","d", 1400600); TH1D yD4; ((TH1D*)gDirectory->Get(name))->Copy(yD4);
  sprintf(name,"%s%d","h", 500600);  TH1D yM1; ((TH1D*)gDirectory->Get(name))->Copy(yM1);
  sprintf(name,"%s%d","h", 1500600); TH1D yM2; ((TH1D*)gDirectory->Get(name))->Copy(yM2);
  sprintf(name,"%s%d","d", 500600);  TH1D yM3; ((TH1D*)gDirectory->Get(name))->Copy(yM3);
  sprintf(name,"%s%d","d", 1500600); TH1D yM4; ((TH1D*)gDirectory->Get(name))->Copy(yM4);


  sprintf(name, "%s%s%s%s%s%s", dir.Data(),"/",flag.Data(),"Env",thevar,".root"); 
  fHistFile = new TFile(name, "RECREATE"); 
  cout << "Opened " << fHistFile->GetName() << endl;

  fHistFile->cd();
  TH1D *hr1 = new TH1D(yD1); hr1->SetName("EnvDaHAc"); 
  TH1D *hr2 = new TH1D(yD2); hr2->SetName("EnvDaHBc"); 
  TH1D *hr3 = new TH1D(yD3); hr3->SetName("EnvDaDAc"); 
  TH1D *hr4 = new TH1D(yD4); hr4->SetName("EnvDaDBc"); 
  TH1D *hr5 = new TH1D(yM1); hr5->SetName("EnvMcHAc"); 
  TH1D *hr6 = new TH1D(yM2); hr6->SetName("EnvMcHBc"); 
  TH1D *hr7 = new TH1D(yM3); hr7->SetName("EnvMcDAc"); 
  TH1D *hr8 = new TH1D(yM4); hr8->SetName("EnvMcDBc"); 
  
  cout << "Writing " << fHistFile->GetName() << endl;
  fHistFile->Write();  fHistFile->Close();
  delete fHistFile;

}


void CompTool::effplots(TString dir, TString flag){

  //   recoilAnalysis b;
   char name[100], line[100];
   sprintf(name,"%s%d" ,"h", 400600);   TH1D *y1(((TH1D*)gDirectory->Get(name)));
   sprintf(name,"%s%d" ,"h", 1400600);   TH1D *y2(((TH1D*)gDirectory->Get(name)));
   sprintf(name,"%s%d" ,"d", 400600);   TH1D *y3(((TH1D*)gDirectory->Get(name)));
   sprintf(name,"%s%d" ,"d", 1400600);   TH1D *y4(((TH1D*)gDirectory->Get(name)));
   sprintf(name,"%s%d" ,"h", 500600);   TH1D *y5(((TH1D*)gDirectory->Get(name)));
   sprintf(name,"%s%d" ,"h", 1500600);   TH1D *y6(((TH1D*)gDirectory->Get(name)));
   sprintf(name,"%s%d" ,"d", 500600);   TH1D *y7(((TH1D*)gDirectory->Get(name)));
   sprintf(name,"%s%d" ,"d", 1500600);   TH1D *y8(((TH1D*)gDirectory->Get(name)));

   TH1D *hratio1 = new TH1D(*y1); hratio1->SetName("hratio1"); hratio1->Reset();
   TH1D *hratio2 = new TH1D(*y5); hratio2->SetName("hratio2"); hratio2->Reset();
   TH1D *hratio3 = new TH1D(*y3); hratio3->SetName("hratio3"); hratio3->Reset();
   TH1D *hratio4 = new TH1D(*y7); hratio4->SetName("hratio4"); hratio4->Reset();
   TH1D *hratioratio1 = new TH1D(*hratio1); hratioratio1->SetName("hratioratio1"); hratioratio1->Reset();
   TH1D *hratioratio2 = new TH1D(*hratio3); hratioratio2->SetName("hratioratio2"); hratioratio2->Reset();

   tl.SetNDC(kTRUE);
   hratio1->Divide(y1, y2, 1, 1, "B");   hratio2->Divide(y5, y6, 1, 1, "B");
   hratio3->Divide(y3, y4, 1, 1, "B");   hratio4->Divide(y7, y8, 1, 1, "B");
   hratioratio1->Divide(hratio1, hratio2);
   hratioratio2->Divide(hratio3, hratio4);

   c7 = new TCanvas("c7", "c7", 300.,   0., 400,800);  
   // -- top left
   TPad *fuPads[5];

   fuPads[1]= new TPad("pad1", "", 0.00, 0.65, 0.99, 0.99);   fuPads[1]->Draw(); 
   fuPads[2]= new TPad("pad2", "", 0.00, 0.50, 0.99, 0.64);   fuPads[2]->Draw();  
   // -- bottom left
   fuPads[3]= new TPad("pad3", "", 0.00, 0.15, 0.99, 0.49);   fuPads[3]->Draw(); 
   fuPads[4]= new TPad("pad4", "", 0.00, 0.00, 0.99, 0.14);   fuPads[4]->Draw(); 

   fuPads[1]->cd(); shrinkPad(0.001, 0.2); 

//   double max = 1.4*hratio1->GetMaximum();
   hratio1->SetLabelSize(0.07, "Y");   hratio1->SetMaximum(1.);
   hratio1->SetMinimum(-0.1);   b.setFilledHist(hratio1 , kBlack, kBlue, 3005);
   hratio1->Draw();

   b.setFilledHist(hratio2 , kBlack, kRed, 3005);
   hratio2->Draw("histsame");
   sprintf(name, "%s%s%s", dir.Data(),"/",flag.Data());
   double g = chisqerr(hratio1,hratio2,name);
   sprintf(line, "#chi^{2} = %5.4f", g); tl.SetTextSizePixels(50); tl.DrawLatex(0.22, 0.75, line);
   
   fuPads[2]->cd(); shrinkPad(0.001, 0.2); 
   shrinkPad(0.4, 0.2, 0.1, 0.001);

   gPad->SetGridx(1);  gPad->SetGridy(1);

   hratioratio1->SetMinimum(0.5); hratioratio1->SetMaximum(1.5);
   hratioratio1->SetMarkerStyle(24);   hratioratio1->SetNdivisions(504, "Y");
   hratioratio1->SetLabelSize(0.22, "X");  hratioratio1->SetLabelSize(0.17, "Y");
   hratioratio1->SetStats(0);   hratioratio1->SetTitle();   hratioratio1->Draw();

   fuPads[3]->cd(); shrinkPad(0.001, 0.2); 

   hratio3->SetLabelSize(0.07, "Y");   hratio3->SetMaximum(1.);
   hratio3->SetMinimum(-0.1);   b.setFilledHist(hratio3 , kBlack, kBlue, 3005);
   hratio3->Draw();

   b.setFilledHist(hratio4 , kBlack, kRed, 3005);
   hratio4->Draw("histsame");
   sprintf(name, "%s%s%s", dir.Data(),"/",flag.Data());
   g = chisqerr(hratio3,hratio4,name);
   sprintf(line, "#chi^{2} = %5.4f", g); tl.SetTextSizePixels(50); tl.DrawLatex(0.22, 0.75, line);
   
   
   fuPads[4]->cd(); shrinkPad(0.001, 0.2); 
   shrinkPad(0.4, 0.2, 0.1, 0.001);

   gPad->SetGridx(1);  gPad->SetGridy(1);

   hratioratio2->SetMinimum(0.5); hratioratio2->SetMaximum(1.5);
   hratioratio2->SetMarkerStyle(24);   hratioratio2->SetNdivisions(504, "Y");
   hratioratio2->SetLabelSize(0.22, "X");  hratioratio2->SetLabelSize(0.17, "Y");
   hratioratio2->SetStats(0);   hratioratio2->SetTitle();   hratioratio2->Draw();

   sprintf(name, "%s%s%s%s%d%s%s", dir.Data(),"/",flag.Data(),"Ceff",multipl,thevar,".ps");
   c7->Print(name);  
}

// ----------------------------------------------------------------------
Double_t CompTool::smeargauss(double invalue, double mean, double sigma){    

  // smears using a gaussian distribution

  double returnvalue=0.;
  double xsmear = gRandom->Gaus(SHIFTNEUT, SIGMANEUT); 
  returnvalue = invalue + xsmear;
  //cout << sigma << mean << " " << invalue << " " << returnvalue-mean-invalue << endl;
  ((TH1D*)gDirectory->Get("h8888"))->Fill(returnvalue - invalue);
  return returnvalue;
}

CompTool::CompTool(const char *var, double mi, double ma, int b, int Sys)
{
   themin = mi;
   themax = ma;
   thebins = b;
   thevar = var ;
   initRest();
   TRandom random(0);
   random.SetSeed(0);
   int therandom = 0;
   dImode = 0;
   if(Sys >= 0) {
     dImode = Sys;
     therandom = random.Rndm() * 1000000;
     if(Sys == 0) {
       Dvar = new recoilDSys("ddecay.table",therandom,Sys);
     } else {
       Dvar = new recoilDSys("dIdecay.table",therandom,Sys);
     }
     Bsem = new recoilDSys(therandom);
   } else {
     Dvar = new recoilDSys("ddecay.table",therandom,0);
     Bsem = new recoilDSys(therandom);
   }
}

CompTool::~CompTool()
{
  //   if (!fChain) return;
  //   delete fChain->GetCurrentFile();
}

Int_t CompTool::GetEntry(Int_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Int_t CompTool::LoadTree(Int_t entry)
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

void CompTool::Show(Int_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t CompTool::Cut(Int_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

int CompTool::hist(double mx){

  // categories
  int bin = int((mx-themin)/((themax-themin)/thebins)+1);
  if (mx>themax || mx==themax) bin = thebins;
if (mx<themin || mx==themin) bin = 1;
  return bin;
}

void CompTool::Init(TTree *tree)
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
   fChain->SetBranchAddress("Gvxbtyp",&Gvxbtyp);
   fChain->SetBranchAddress("GSem",&GSem);  
   fChain->SetBranchAddress("GfDpi",&GfDpi); 
   fChain->SetBranchAddress("GfDpiz",&GfDpiz);
   fChain->SetBranchAddress("GfDk",&GfDk);  
   fChain->SetBranchAddress("GfDks",&GfDks); 
   fChain->SetBranchAddress("GfDlep",&GfDlep);
   fChain->SetBranchAddress("GfD0Ds",&GfD0Ds);
   fChain->SetBranchAddress("GfDDs",&GfDDs);
   fChain->SetBranchAddress("GfDgam",&GfDgam);
   fChain->SetBranchAddress("intpur",&intpur);
   fChain->SetBranchAddress("mode",&mode);
   fChain->SetBranchAddress("nnpi0",&nnpi0);
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

Bool_t CompTool::Notify()
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
   b_Gvxbtyp=    fChain->GetBranch("Gvxbtyp");
   b_GSem=       fChain->GetBranch("GSem");  
   b_GfDpi=      fChain->GetBranch("GfDpi");  
   b_GfDpiz=     fChain->GetBranch("GfDpiz");
   b_GfDk=       fChain->GetBranch("GfDk");  
   b_GfDks=      fChain->GetBranch("GfDks"); 
   b_GfDlep=     fChain->GetBranch("GfDlep");
   b_GfD0Ds=     fChain->GetBranch("GfD0Ds");
   b_GfDDs=      fChain->GetBranch("GfDDs");
   b_GfDgam=     fChain->GetBranch("GfDgam");
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


// ----------------------------------------------------------------------
void CompTool::dumpCuts() {
  cout << "====================================" << endl;
  cout << "Cut file " << fCutFile << endl; 
  cout << "------------------------------------" << endl;
  cout << "deSignal:          " << DESIGNALLO << " ... " << DESIGNALHI << endl;
  cout << "pcms:              " << CUTPCMS   << endl;
  cout << "intPurity:         " << INTPURITY << endl;
  cout << "Mx max cut:        " << MXMAXCUT  << endl;
  cout << "Mm2 cut:           " << MM2CUT    << endl;
  cout << "D from D:          " << DFROMD    << endl;
  cout << "D0B:               " << DOBDECWEIGHT << endl;
  cout << "D0D:               " << DODDECWEIGHT << endl;
  cout << "====================================" << endl;
}


// ----------------------------------------------------------------------
void CompTool::readCuts(TString filename, int dump) {
  char  buffer[200];
  fCutFile = filename;
  sprintf(buffer, "%s", filename.Data());
  ifstream is(buffer);
  char CutName[100];
  float CutValue;   int ok(0);

  while (is.getline(buffer, 200, '\n')) {
    ok = 0;
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);
    // -- breco 
    if (!strcmp(CutName, "deSignalLo")) {DESIGNALLO = CutValue; ok = 1;}
    if (!strcmp(CutName, "deSignalHi")) {DESIGNALHI = CutValue; ok = 1;}
    if (!strcmp(CutName, "pcms"))       {CUTPCMS    = CutValue; ok = 1;}
    if (!strcmp(CutName, "intPurity"))  {INTPURITY  = CutValue; ok = 1;}
    if (!strcmp(CutName, "mxcut"))      {MXMAXCUT   = CutValue; ok = 1;}
    if (!strcmp(CutName, "mm2cut"))     {MM2CUT   = CutValue; ok = 1;}
    if (!strcmp(CutName, "dfromd"))     {DFROMD     = CutValue; ok = 1;}
    if (!strcmp(CutName, "d0B"))        {DOBDECWEIGHT  = CutValue; ok = 1;}
    if (!strcmp(CutName, "d0D"))        {DODDECWEIGHT  = CutValue; ok = 1;}
    if (ok == 0)  cout << "==> recoilNtp::readCuts() Error: Don't know about variable " << CutName << endl;
  }

  if (dump == 1) dumpCuts();

  readintpur();
}

// ----------------------------------------------------------------------
void CompTool::readintpur(){    

   char buffer[200];
   float bmode, dmode, sig, bkg, pur, sb;
   ifstream is("tables/tablepurity.dat");
   int mode;
   while (is.getline(buffer, 200, '\n')) {
     if (buffer[0] == '#') {continue;}
     sscanf(buffer, "%f %f %f %f %f %f", &bmode, &dmode, &sig, &bkg, &pur, &sb);
     mode = (dmode+100) * 100 + bmode-10000;
     brecosig[mode] = sig;	
     brecobkg[mode] = bkg;
     brecointpur[mode] = pur; 
     
   }
	
}
void CompTool::initRest() {
  // -- CUTS
  DESIGNALLO = -0.1; 
  DESIGNALHI =  0.1; 
  CUTPCMS    = 1.;
  INTPURITY  = 0.;
  MXMAXCUT  = 5.;
  MM2CUT  = .5;
  DFROMD  = 1;
  DOBDECWEIGHT = 0;
  DOBDECWEIGHT = 0;
}

// ----------------------------------------------------------------------
double CompTool::getBsysweight(int decType,int thevub) {
  double theweight;
  theweight = 1.;  
  if(DOBDECWEIGHT) theweight = Bsem->weight(decType); 
  if(thevub) theweight = 1.;
  return theweight;
}
// ----------------------------------------------------------------------
double CompTool::getDsysweight(int decDpi,int decDk,int decDks,int decDpiz,int decDlep,int decImode,int thevub) {
  double theweight;
  theweight = 1.;  
  int bla(1);
  if(DODDECWEIGHT){
   theweight = Dvar->weight(decDpi,decDk,decDks,decDpiz,decDlep,decImode,bla);    
  }
  if(thevub) theweight = 1.;
  return theweight;
}
