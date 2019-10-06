#include "TheorTool.hh"
#include "recoilDSys.hh"
#include "TCut.h"
#include <TMath.h>
#include <TTree.h>
#include <TRandom.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLegendEntry.h>
#include <TH1.h>
#include <TH2.h>
#include <TPaveStats.h>
#include <fstream.h>
#include <iostream.h>
#include <iomanip.h>
#include "util.hh"


TheorTool::TheorTool(const char *var, double mi, int ma, const char *b, double fit, int seed)
{
  BMASS   = 5.2792;  BQMASS  = 4.800;  A0      = 1.29;
  themin = mi; //Reweight Excl Decays
  themax = ma; //Change A0 +/- one sigma
  thebins = 21;  thevar = var;
  theseed = seed;
  const double scalfDAMC = 81.9/900;
  double histbA[21] = {0.0, 1.2, 1.264, 1.341, 1.414, 1.483, 1.549, 1.612, 1.673, 1.732, 1.788, 1.843, 1.897, 1.949, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0};
  for(int ist=0; ist<21; ist++) {
    histbins[ist] = histbA[ist];
  } 
  TRandom random(theseed);
  rndm = random;
  //  rnd = randomized(BQMASS,0.090);
  // rnd = BQMASS;
  rnd = BQMASS;
  //Now the BRexcl becomes something that we can 'change'
  //Error is needed
  const double BRxulnu = 20.34; //Fit outp
  //const double BRxulnu = 20.34; //Fit outp
  const double BRxulnuErr = 5.61;//Within stat error.
  const double BRxlnu = 0.1087; //Ref
  const double BRexclMeasB0 = 0.735;
  //  const double BRexclMeasB0 = 0.565;
  const double BRexclErrMeasB0 = 0.17;//Tobe changed
  const double BRexclMeasBch = 0.73;
  //const double BRexclMeasBch = 0.55;
  const double BRexclErrMeasBch = 0.18;//Tobe changed
  BRexclB0 = BRexclMeasB0;
  BRexclBch = BRexclMeasBch;
  fittedBR = BRxulnu*BRxlnu;
  if(themin) {
    BRexclB0 = randomized(BRexclMeasB0,BRexclErrMeasB0);
    BRexclBch = randomized(BRexclMeasBch,BRexclErrMeasBch);
  }
  if(fit) fittedBR = randomized(BRxulnu*BRxlnu,BRxulnuErr*BRxlnu); 
  //  cout<<"Excl and Xulnu BR:: "<<BRexclB0<< " "<<BRexclBch<< " "<<fittedBR<<" "<<rnd<<endl;
  char buffer[200]; char CutName[100];
  float CutValue;
  sprintf(buffer,"%s%s","sysWd/plot_",b);
  //  cout<<"buffer:: "<<buffer<<endl;
  ifstream is(buffer);
  for (int y=0;y<42;y++){
    is.getline(buffer, 200, '\n');
    sscanf(buffer, "%s %f", CutName, &CutValue);
    TrueMxWeight[y] = CutValue;
  }
  initRest();
  char fermname[200];
  sprintf(fermname,"%s%d%s","fermiparam_",themax,".dat");
  ofstream outFile(fermname,ios::app);
  outFile<<"deltamb  "<<rnd - BQMASS<<endl;
  if(themax == 1) {
    outFile<<"deltaa  "<<0<<endl;
  } else if (themax == 2) {
    outFile<<"deltaa  "<<2.31<<endl;
  } else if (themax == 3) {
    outFile<<"deltaa  "<<-0.91<<endl;
  }

}


TheorTool::~TheorTool()
{
  //   if (!fChain) return;
  //   delete fChain->GetCurrentFile();
}

void TheorTool::Loop(int nevents, int cat, int isbch, int doMbaRew)
{
  gROOT->cd();
  //id = Mx category == Mx bin
  int id, id2, idflav, dcat, group;        
  char name[100]; char le[100];
    if (fChain == 0) return;
    Int_t nentries = Int_t(fChain->GetEntries());
  if( nentries > nevents) nentries = nevents;
  Int_t nbytes = 0, nb = 0;
  
  for (Int_t jentry=0; jentry<nentries;jentry++) {

    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // cuts
    if(!vub) continue;
    int ch = xcharge + brecocharge;  // total charge

    char MyVar[100];  char tmpA1[100];
    char tmpB1[100];  char tmpNo[100]; 
    char statr[100];    char stattmp[100];
    double myvar = mxhadfit;
    //    double myvar2 = mm2;

    sprintf(MyVar,"%s",thevar);
    sprintf(tmpA1,"mxhadgen");    sprintf(tmpB1,"mxhadfit");
    sprintf(tmpNo,"NoComp");
    if (!strcmp(MyVar,tmpA1)) {myvar = mxhadgen;}
    if (!strcmp(MyVar,tmpB1)) {myvar = mxhadfit;}
    //    cout<<myvar<<" "<<mxhadfit<<" "<<mxhadgen<<endl;
    double w = 1;    //calculate reweightings
    
    //Mb and a reweighting
    //Function of deltaM on b quark and kplus
    //A is taken and varied +/- one sigma

    if (doMbaRew && cat == 4) {
      if(fkplus>-1) {
	//	cout<<"The max :: "<<themax<<endl;
	w *= reweightMba(fkplus,mxhadgen,mxhadfit,themax);
      }
    }
    //    if(cat == 6 || cat == 7){
    if(intpur<0.2) continue;
    //    } 
    double weight; char mesf[100];char mesftmp[100];
    double weightNomba;
    //Adding Fermi reweighting
    weight = totweight * w;
    weightNomba = totweight;
    intStudy(vub,myvar,Gvxbtyp,fBchgen,cat);
    //    cout<<"Weight: "<<w<<endl;
    if(cat == 5 || cat == 7){
      //      cout<<"Value of Gvxbtyp:: "<<Gvxbtyp<<endl;
      if(TMath::Abs(Gvxbtyp)<25) {
	sprintf(mesf,"%s","hy");
	if(Gvxbtyp!=7&&Gvxbtyp!=-7) {
	  sprintf(stattmp,"%s","re");
	} else {
	  sprintf(stattmp,"%s","nr");
	}
      }
    } else {
      sprintf(mesf,"%s","pn");
      sprintf(stattmp,"%s","da");
    }

    sprintf(statr,"%s%s","nu",stattmp);
    if(Gvxbtyp<0)     sprintf(statr,"%s%s","ch",stattmp);

    if (!((myvar>0) || (myvar<0) || (myvar == 0))) {
      myvar = -999.;
      cout<<"Something weird"<<endl;
    }   
    
    if(cat == 4) sprintf(le,"%s","DA");
    if(cat == 5) sprintf(le,"%s","MC");
    if(cat == 6) sprintf(le,"%s","GDA");
    if(cat == 7) sprintf(le,"%s","GMC");

    if(cat < 6) {
      id = hist(myvar);
      idflav = 200 + id;    
      sprintf(name, "%s%s%d%s","sc",le,idflav,statr);
      ((TH1D*)gDirectory->Get(name))->Fill(mes,weight);  
      sprintf(name, "%s%s%d%s","scNMB",le,idflav,statr);
      ((TH1D*)gDirectory->Get(name))->Fill(mes,weightNomba);  
    }

    sprintf(name, "%s%s%s","Allsc",le,statr);
    ((TH1D*)gDirectory->Get(name))->Fill(mes,weight);  
    sprintf(name, "%s%s%s","AllscNMB",le,statr);
    ((TH1D*)gDirectory->Get(name))->Fill(mes,weightNomba);  

    sprintf(mesftmp, "%s%s",le,mesf);
    ((TH1D*)gDirectory->Get(mesftmp))->Fill(mes,weightNomba);  

  }
}

void TheorTool::Bookhist()
{
  gROOT->cd();
  TH1 *h;  char name[100], title[200], number[100]; TH2 *h2;
  char le[100];  char ler[100]; char fch[100];
  Float_t ranges[21] = {0.0, 1.2, 1.264, 1.341, 1.414, 1.483, 1.549, 1.612, 1.673, 1.732, 1.788, 1.843, 1.897, 1.949, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0};
  //  sprintf(name,"%s%d%s%d", "B",300000+iflav*100,"ct",lcat);  sprintf(title, "%s%s", thevar, " data events after all cuts");  h2 = new TH2D(name, title, thebins, themin, themax, thebins2, themin2, themax2 );  h2->Sumw2();

  //  for(int ist=0; ist<16; ist++) {
  //  ranges[ist] = histbins[ist];
  //} 
  //Loop for B0, B+
  for (int lchg = 0; lchg <2; lchg++) {
    if(lchg == 0) sprintf(fch,"%s","nu");
    if(lchg == 1) sprintf(fch,"%s","ch");
    //Loop for resonant, non reso, data
    for (int catr = 0; catr <3; catr++) {
      if(catr == 0) sprintf(ler,"%s%s",fch,"re");
      if(catr == 1) sprintf(ler,"%s%s",fch,"nr");
      if(catr == 2) sprintf(ler,"%s%s",fch,"da");
      
      //Loop for data and MC
      for (int cat = 0; cat <4; cat++) {
	if(cat == 0) sprintf(le,"%s","DA");
	if(cat == 1) sprintf(le,"%s","MC");
	if(cat == 2) sprintf(le,"%s","GDA");
	if(cat == 3) sprintf(le,"%s","GMC");
	
	
	//Mes plots for all the candidates
	sprintf(name,"%s%s%s","Allsc",le,ler);  sprintf(title, "mes (lepton cuts");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
	sprintf(name,"%s%s%s","AllscNMB",le,ler);  sprintf(title, "mes (lepton cuts");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
	sprintf(name,"%s%s%s","Allac",le,ler);  sprintf(title, "mes (all cuts)");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();

	//Mes plots per bin
	if(cat <2) {
	  for(int ic=0; ic<thebins+1; ic++) {
	    sprintf(name,"%s%s%d%s","sc",le,200+ic,ler);  sprintf(title, "mes (lepton cuts");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
	    sprintf(name,"%s%s%d%s","scNMB",le,200+ic,ler);  sprintf(title, "mes (lepton cuts");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
	    //      cout<<"Hist : "<<name<<endl;
	    sprintf(name,"%s%s%d%s","ac",le,200+ic,ler);  sprintf(title, "mes (all cuts)");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
	  }
	}

	//    cout<<thebins<<" bins"<<endl;
	sprintf(name,"%s%s%s","Mxsc",le,ler);  sprintf(title, "Mx (lepton cuts)");  h = new TH1D(name, title, 20, ranges);  h->Sumw2();
	sprintf(name,"%s%s%s","MxscNMB",le,ler);  sprintf(title, "Mx (lepton cuts)");  h = new TH1D(name, title, 20, ranges);  h->Sumw2();
	sprintf(name,"%s%s%s","Mxac",le,ler);  sprintf(title, "Mx (all cuts)");  h = new TH1D(name, title, 20, ranges);  h->Sumw2();
	//  sprintf(name,"%s%d","B93206ct",lcat);  sprintf(title, "mes data after all cuts: depleted");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
      }    
    }
  }
  sprintf(name,"%s","MxNRSub");  sprintf(title, "Mx (all cuts)");  h = new TH1D(name, title, 20, ranges);  h->Sumw2();
  sprintf(name,"%s","GMxNRSub");  sprintf(title, "Mx (all cuts)");  h = new TH1D(name, title, 20, ranges);  h->Sumw2();

  //Mb and a reweighting
  h = new TH1D("nw1", "weights", 200, 0., 5.);
  h = new TH1D("nw2", "weights", 200, 0., 5.);
  h = new TH1D("nw3", "weights", 200, 0., 5.);
  h = new TH1D("nw4", "weights", 200, 0., 5.);
  h = new TH1D("nk0", "f(kPlus)", 100, -2., 1.);
  h = new TH1D("nk1", "f1(kPlus)", 100, -2., 1.);
  h = new TH1D("nk2", "f2(kPlus)", 100, -2., 1.);
  h = new TH1D("nk3", "f3(kPlus)", 100, -2., 1.);
  h = new TH1D("nk4", "f4(kPlus)", 100, -2., 1.);
  
  h = new TH1D("nh0", "Mxgen", 80, 0., 4.);
  h = new TH1D("nh1", "Mxgen", 80, 0., 4.);
  h = new TH1D("nh2", "Mxgen", 80, 0., 4.);
  h = new TH1D("nh3", "Mxgen", 80, 0., 4.);
  h = new TH1D("nh4", "Mxgen", 80, 0., 4.);
  h = new TH1D("ni0", "Mxgen (integrated)", 80, 0., 4.);
  h = new TH1D("ni1", "Mxgen (integrated)", 80, 0., 4.);
  h = new TH1D("ni2", "Mxgen (integrated)", 80, 0., 4.);
  h = new TH1D("ni3", "Mxgen (integrated)", 80, 0., 4.);
  h = new TH1D("ni4", "Mxgen (integrated)", 80, 0., 4.);
  h = new TH1D("nm1", "Mxfit (integrated)", 80, 0., 4.);
  h = new TH1D("nm2", "Mxfit (integrated)", 80, 0., 4.);
  h = new TH1D("nm3", "Mxfit (integrated)", 80, 0., 4.);
  h = new TH1D("nm4", "Mxfit (integrated)", 80, 0., 4.);
  h = new TH1D("ne0", "Rel. Error due to cut on Mxgen", 80, 0., 4.);
  h = new TH1D("ne02", "Rel. Error due to cut on Mxgen", 80, 0., 4.);
  h = new TH1D("ne1", "Rel. Error due to cut on Mxgen", 80, 0., 4.);
  h = new TH1D("ne2", "Rel. Error due to cut on Mxgen", 80, 0., 4.);
  h = new TH1D("ne3", "Rel. Error due to cut on Mxgen", 80, 0., 4.);
  h = new TH1D("ne4", "Rel. Error due to cut on Mxgen", 80, 0., 4.);
  
  h = new TH1D("nH0", "Mxfit", 80, 0., 4.);
  h = new TH1D("nH1", "Mxfit", 80, 0., 4.);
  h = new TH1D("nH2", "Mxfit", 80, 0., 4.);
  h = new TH1D("nH3", "Mxfit", 80, 0., 4.);
  h = new TH1D("nH4", "Mxfit", 80, 0., 4.);
  h = new TH1D("nI0", "Mxfit (integrated)", 80, 0., 4.);
  h = new TH1D("nI1", "Mxfit (integrated)", 80, 0., 4.);
  h = new TH1D("nI2", "Mxfit (integrated)", 80, 0., 4.);
  h = new TH1D("nI3", "Mxfit (integrated)", 80, 0., 4.);
  h = new TH1D("nI4", "Mxfit (integrated)", 80, 0., 4.);
  h = new TH1D("nM1", "Mxfit (integrated)", 80, 0., 4.);
  h = new TH1D("nM2", "Mxfit (integrated)", 80, 0., 4.);
  h = new TH1D("nM3", "Mxfit (integrated)", 80, 0., 4.);
  h = new TH1D("nM4", "Mxfit (integrated)", 80, 0., 4.);
  h = new TH1D("nE0", "Rel. Error due to cut on Mxfit", 80, 0., 4.);
  h = new TH1D("nE02", "Rel. Error due to cut on Mxfit", 80, 0., 4.);
  h = new TH1D("nE1", "Rel. Error due to cut on Mxfit", 80, 0., 4.);
  h = new TH1D("nE2", "Rel. Error due to cut on Mxfit", 80, 0., 4.);
  h = new TH1D("nE3", "Rel. Error due to cut on Mxfit", 80, 0., 4.);
  h = new TH1D("nE4", "Rel. Error due to cut on Mxfit", 80, 0., 4.);
  //Integrals study
  h = new TH1D("plotall","Mx Hibryd",100,0.,5.);
  h = new TH1D("plotnres","Mx Hybrid",100,0.,5.);
  h = new TH1D("plotres","Mx Resonant",100,0.,5.);
  
  h = new TH1D("plotallb0","Mx Hibryd",100,0.,5.);
  h = new TH1D("plotnresb0","Mx Hybrid",100,0.,5.);
  h = new TH1D("plotresb0","Mx Resonant",100,0.,5.);
  
  h = new TH1D("plotallchb","Mx Hibryd",100,0.,5.);
  h = new TH1D("plotnreschb","Mx Hybrid",100,0.,5.);
  h = new TH1D("plotreschb","Mx Resonant",100,0.,5.);

  h = new TH1D("plotnonres","Mx Non-resonant",100,0.,5.);
  h = new TH1D("plotnonresb0","Mx Non-resonant",100,0.,5.);
  h = new TH1D("plotnonreschb","Mx Non-resonant",100,0.,5.);

  h2 = new TH2D("inthyb","Mx Fraction ",100,0.,5.,500,0.,1.2);
  h2 = new TH2D("intres","Mx Fraction",100,0.,5.,500,0.,1.2);
  h2 = new TH2D("intnres","Mx Fraction",100,0.,5.,500,0.,1.2);

  h2 = new TH2D("inthybchb","Mx Fraction ",100,0.,5.,500,0.,1.2);
  h2 = new TH2D("intreschb","Mx Fraction",100,0.,5.,500,0.,1.2);
  h2 = new TH2D("intnreschb","Mx Fraction",100,0.,5.,500,0.,1.2);

  h2 = new TH2D("inthybb0","Mx Fraction ",100,0.,5.,500,0.,1.2);
  h2 = new TH2D("intresb0","Mx Fraction",100,0.,5.,500,0.,1.2);
  h2 = new TH2D("intnresb0","Mx Fraction",100,0.,5.,500,0.,1.2);

  //Mes fits
  
  sprintf(name,"%s","DApn",le,ler);  sprintf(title, "mes no cuts");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  sprintf(name,"%s","DAhy",le,ler);  sprintf(title, "mes no cuts");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  sprintf(name,"%s","MCpn",le,ler);  sprintf(title, "mes no cuts");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  sprintf(name,"%s","MChy",le,ler);  sprintf(title, "mes no cuts");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  sprintf(name,"%s","GDApn",le,ler);  sprintf(title, "mes no cuts");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  sprintf(name,"%s","GDAhy",le,ler);  sprintf(title, "mes no cuts");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  sprintf(name,"%s","GMCpn",le,ler);  sprintf(title, "mes no cuts");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  sprintf(name,"%s","GMChy",le,ler);  sprintf(title, "mes no cuts");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();

}

void TheorTool::Fitmes(int cat, int cut, TString flag){

  // fit to mes distribution and fill of the Mx plots...
  gROOT->cd();
   
  char name[100], namedi[100], le[100], cutn[100];
  double sigs[1000], errsigs[1000];
  /*
  double tempbin, tempbinchb, tempbinb0os;
  double tempbinb0ss, temperr, temperrchb, temperrb0os;
  double temperrb0ss, chid = 0.174; 
  */
  double foomean, foosigma, fooalpha, foon, usedmean, usedsigma;
  double usedalpha, usedn;
  //  char nameou[100], namebch[100], nameb0os[100], nameb0ss[100];
  char ler[100]; int Ncat =1;
  //Separing Data and MC
  if(cat == 4) {
    sprintf(le,"%s","DA");
  } else if(cat == 5) {
    sprintf(le,"%s","MC");
    Ncat =2;
  } else if(cat == 6) {
    sprintf(le,"%s","GDA");
  } else if(cat == 7) {
    sprintf(le,"%s","GMC");
    Ncat =2;
  }
  char fch[100];
  c4 = new TCanvas("c4", "c4",600,800);
  for (int lchg = 0; lchg <2; lchg++) {
    if(lchg == 0) sprintf(fch,"%s","nu");
    if(lchg == 1) sprintf(fch,"%s","ch");

    for (int catr = 0; catr <Ncat; catr++) {
      if(catr == 0 && cat == 4) {
	sprintf(ler,"%s%s",fch,"da");
      } else if(catr== 0 && cat ==5) {
	sprintf(ler,"%s%s",fch,"re");
      } else if(catr== 0 && cat ==6) {
	sprintf(ler,"%s%s",fch,"da");
      } else if(catr== 0 && cat ==7) {
	sprintf(ler,"%s%s",fch,"re");
      }
      if(catr == 1) sprintf(ler,"%s%s",fch,"nr");
      
      //      cout<<"In mes fit:: "<<ler<<endl;
      int addcut = 0;  int group = 200 + cat * 1000;
      sprintf(cutn,"%s","sc");
      if(cut) sprintf(cutn,"%s","scNMB");
      int const theb = thebins+1;
      //extracting the signal fit parameters 
      
      sprintf(name, "%s%s%s%s","All",cutn,le,ler);  
      double dummy1,dummy2;
      sighisto(dummy1,dummy2,(TH1D*)gDirectory->Get(name),usedmean,usedsigma,usedalpha,usedn,1,-11111111.,-1111111.,-1111111.,-1111111.,-1111111.);
      cout <<name<<" :: "<<dummy1<< " "<< dummy2 << endl;
      sprintf(name, "%s%s%s%s","Mx",cutn,le,ler);  
      c4->Clear(); c4->Divide(5,5); 
      if(cat < 6) {
	for (int i = 1;  i < theb; i++) { 
	  c4->cd(i);
	  sprintf(namedi,"%s%s%d%s",cutn,le,200+i,ler);  
	  sighisto(sigs[i-1],errsigs[i-1],(TH1D*)gDirectory->Get(namedi),foomean,foosigma,fooalpha,foon,1,usedmean,usedsigma,usedalpha,usedn,-1111111);
	  //    cout<<"name "<<name<<" "<<sigs[i-1]<<" "<<errsigs[i-1]<<endl;
	  ((TH1D*)gDirectory->Get(namedi))->Draw();
	  ((TH1D*)gDirectory->Get(name))->SetBinContent(i, sigs[i-1]);
	  ((TH1D*)gDirectory->Get(name))->SetBinError(i, errsigs[i-1]);
	}
      }
      char tmptes[100]; sprintf(tmptes,"%s%s",name,"_mesfit.eps");
      c4->Print(tmptes);
    }
  }

  c4 = new TCanvas("c4", "c4",600,800);   c4->Clear(); c4->cd();
  char nameMf[100];  char nameMfeps[100];
  double dum1,dum2;
  if(cat == 4 || cat == 6) {
    sprintf(nameMf,"%s%s",le,"pn");  
    sighisto(dum1,dum2,(TH1D*)gDirectory->Get(nameMf),usedmean,usedsigma,usedalpha,usedn,1,-11111111.,-1111111.,-1111111.,-1111111.,-1111111.);
    cout <<"Mes fits 'pure n':: "<<nameMf<<" :: "<<dum1<< " "<< dum2 << endl;
  } else {
    sprintf(nameMf,"%s%s",le,"hy");  
    sighisto(dum1,dum2,(TH1D*)gDirectory->Get(nameMf),usedmean,usedsigma,usedalpha,usedn,1,-11111111.,-1111111.,-1111111.,-1111111.,-1111111.);
    cout <<"Mes fits 'hybrid':: "<<nameMf<<" :: "<<dum1<< " "<< dum2 << endl;
  }
  ((TH1D*)gDirectory->Get(nameMf))->Draw();
  sprintf(nameMfeps,"%s%s",nameMf,".eps");
  c4->Print(nameMfeps);
}

void TheorTool::intStudy(int Avub, double Amxhadgen, int AGvxbtyp, int AfBchgen, int Acat) {

  double w =1;  
  if(CUTPCMS && Avub && Acat==4 )  {
    w *= getTrueMxWeight(Amxhadgen,AGvxbtyp);                            //True Mx reweighting only for NonReso MC
  }

  if(Acat == 5) {
    if(w != 1) cout<<"Excl weight :: "<<w<<endl;
    if(Avub)  ((TH1D*)gDirectory->Get("plotall"))->Fill(Amxhadgen);
    if(Avub&&(AGvxbtyp!=7&&AGvxbtyp!=-7)) ((TH1D*)gDirectory->Get("plotres"))->Fill(Amxhadgen,w);
    
    if(Avub && AfBchgen==1) ((TH1D*)gDirectory->Get("plotallchb"))->Fill(Amxhadgen);
    if(Avub&&(AGvxbtyp!=7&&AGvxbtyp!=-7) && AfBchgen==1) ((TH1D*)gDirectory->Get("plotreschb"))->Fill(Amxhadgen,w);
    
    if(Avub && AfBchgen==2)((TH1D*)gDirectory->Get("plotallb0"))->Fill(Amxhadgen);
    if(Avub&&(AGvxbtyp!=7&&AGvxbtyp!=-7) && AfBchgen==2) ((TH1D*)gDirectory->Get("plotresb0"))->Fill(Amxhadgen,w);

  } else if (Acat == 4) {
    if(Avub) ((TH1D*)gDirectory->Get("plotnonres"))->Fill(mxhadgen);
    if(Avub && AfBchgen==1) ((TH1D*)gDirectory->Get("plotnonreschb"))->Fill(mxhadgen);
    if(Avub && AfBchgen==2) ((TH1D*)gDirectory->Get("plotnonresb0"))->Fill(mxhadgen);
    //Non resonant IS nonresonant pure rescaled
    if(Avub)  ((TH1D*)gDirectory->Get("plotnres"))->Fill(Amxhadgen,w);
    if(Avub && AfBchgen==1)((TH1D*)gDirectory->Get("plotnreschb"))->Fill(Amxhadgen,w);
    if(Avub && AfBchgen==2) ((TH1D*)gDirectory->Get("plotnresb0"))->Fill(Amxhadgen,w);
  }
}

void TheorTool::PintStudy() {
  TH1D *hplotall  = ((TH1D*)gDirectory->Get("plotall"));
  TH1D *hplotnres = ((TH1D*)gDirectory->Get("plotnres"));
  TH1D *hplotres  = ((TH1D*)gDirectory->Get("plotres"));

  TH1D *hplotallchb  = ((TH1D*)gDirectory->Get("plotallchb"));
  TH1D *hplotnreschb = ((TH1D*)gDirectory->Get("plotnreschb"));
  TH1D *hplotreschb  = ((TH1D*)gDirectory->Get("plotreschb"));

  TH1D *hplotallb0  = ((TH1D*)gDirectory->Get("plotallb0"));
  TH1D *hplotnresb0 = ((TH1D*)gDirectory->Get("plotnresb0"));
  TH1D *hplotresb0  = ((TH1D*)gDirectory->Get("plotresb0"));

  TH1D *hplotnonres = ((TH1D*)gDirectory->Get("plotnonres"));
  TH1D *hplotnonreschb = ((TH1D*)gDirectory->Get("plotnonreschb"));
  TH1D *hplotnonresb0 = ((TH1D*)gDirectory->Get("plotnonresb0"));

  TH2D *hinthyb =((TH2D*)gDirectory->Get("inthyb"));
  TH2D *hintres =((TH2D*)gDirectory->Get("intres"));
  TH2D *hintnres =((TH2D*)gDirectory->Get("intnres"));

  TH2D *hinthybchb =((TH2D*)gDirectory->Get("inthybchb"));
  TH2D *hintreschb =((TH2D*)gDirectory->Get("intreschb"));
  TH2D *hintnreschb =((TH2D*)gDirectory->Get("intnreschb"));

  TH2D *hinthybb0 =((TH2D*)gDirectory->Get("inthybb0"));
  TH2D *hintresb0 =((TH2D*)gDirectory->Get("intresb0"));
  TH2D *hintnresb0 =((TH2D*)gDirectory->Get("intnresb0"));

  Float_t ranges[21] = {0.0,  1.2, 1.264, 1.341, 1.414, 1.483, 1.549, 1.612, 1.673, 1.732, 1.788, 1.843, 1.897, 1.949, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0};
  TH1D *hhy = new TH1D("Sum", "Sum" ,100,0.,5.);//Mymodel
  tl.SetNDC(kTRUE);   tl.SetTextSizePixels(50.1);
  c4 = new TCanvas("c4", "c4",600,800);   c4->Clear(); c4->cd();
  setHist(hplotall,1,0,1,1,3005);   
  setHist(hplotnonres,2,0,1,1,3004);   
  setHist(hplotnres,3,0,1,1,3004);   
  setHist(hplotres,4,0,1,1,3005);
  setHist(hhy,5,0,1,1,3005);

  cout<<"Reso, Non reso, Fraction:: "<<hplotres->Integral()<<" "<< hplotnres->Integral()<<" " <<hplotres->Integral()/hplotnres->Integral()<< endl;
  hhy->Add(hplotres, hplotnres, 1., 1.); 
  hhy->Draw();
  hplotnonres->Draw("same");
  TLegend *leg; 
  TLegendEntry *legge; 
  leg = new TLegend(0.45,0.5,0.69,0.59);
  setLegend(leg,0,0.,0.05,0,hplotnonres,hhy,"Pure Non reso","My model","f");
  leg->Draw();
  //  c1->SaveAs("mxplothibselex_mxgen.eps");
  c4->SaveAs("mxplotmymodel.eps");

  double integralhyb,integralres,integralnres;  
  integralhyb = integralres = integralnres =0;  
  double totintegralhyb=hplotall->Integral();
  double totintegralres=hplotres->Integral();  
  double totintegralNres=hplotnres->Integral();  
  double totintegralnres=hplotnonres->Integral();  
  
  double integralhybb0,integralresb0,integralnresb0;  
  integralhybb0 = integralresb0 = integralnresb0=0;  
  double totintegralhybb0=hplotallb0->Integral();
  double totintegralresb0=hplotresb0->Integral();  
  double totintegralNresb0=hplotnresb0->Integral();  
  double totintegralnresb0=hplotnonresb0->Integral();  
  
  double integralhybchb,integralreschb,integralnreschb;  
  integralhybchb = integralreschb = integralnreschb=0;  
  double totintegralhybchb=hplotallchb->Integral();
  double totintegralreschb=hplotreschb->Integral();  
  double totintegralNreschb=hplotnreschb->Integral();  
  double totintegralnreschb=hplotnonreschb->Integral();  

  double bin;

  for(int iloop=0; iloop<100; iloop++){
    bin = 0.05 * iloop+0.025;
    //    integralhyb+=hplotall->GetBinContent(iloop);
    integralhyb+=(hplotres->GetBinContent(iloop) + hplotnres->GetBinContent(iloop))/(totintegralres + totintegralNres);
    integralres+=hplotres->GetBinContent(iloop);
    integralnres+=hplotnonres->GetBinContent(iloop);
    hinthyb->Fill(bin,integralhyb);
    hintres->Fill(bin,integralres/totintegralres);
    hintnres->Fill(bin,integralnres/totintegralnres);
    //    integralhybb0+=hplotallb0->GetBinContent(iloop);
    integralhybb0+=(hplotresb0->GetBinContent(iloop) + hplotnresb0->GetBinContent(iloop))/(totintegralresb0 + totintegralNresb0);
    integralresb0+=hplotresb0->GetBinContent(iloop);
    integralnresb0+=hplotnonresb0->GetBinContent(iloop);
    hinthybb0->Fill(bin,integralhybb0);
    hintresb0->Fill(bin,integralresb0/totintegralresb0);
    hintnresb0->Fill(bin,integralnresb0/totintegralnresb0);
    //    integralhybchb  += hplotallchb->GetBinContent(iloop);
    integralhybchb  += (hplotreschb->GetBinContent(iloop) + hplotnreschb->GetBinContent(iloop))/(totintegralreschb + totintegralNreschb);
    integralreschb  += hplotreschb->GetBinContent(iloop);
    integralnreschb += hplotnonreschb->GetBinContent(iloop);
    hinthybchb->Fill(bin,integralhybchb);
    hintreschb->Fill(bin,integralreschb/totintegralreschb);
    hintnreschb->Fill(bin,integralnreschb/totintegralnreschb);
  }
  c4->Clear();  c4->cd();
  hinthyb->SetMarkerColor(kBlue); 
  hinthyb->SetMarkerSize(5); 
  hinthyb->SetStats(0);
  hintnres->SetStats(0);
  hintnres->SetMarkerSize(5); 
  hintnres->SetMarkerColor(kRed); 
  hintnres->Draw();
  hinthyb->DrawCopy("same");
  TLegend *leg2; 
  TLegendEntry *legge2; 
  leg2 = new TLegend(0.45,0.5,0.69,0.59);
  leg2->SetFillStyle(0); leg2->SetBorderSize(0.); leg2->SetTextSize(0.05); 
  leg2->SetFillColor(0); 
  legge2 = leg2->AddEntry(hinthyb, "Hybrid", "p"); 
  legge2 = leg2->AddEntry(hintnres, "Non-res", "p"); 
  leg2->Draw();
  //  c1->SaveAs("integralsselex_mxgen.eps");
  c4->SaveAs("integralsselex_mxgen.eps");

  c4->Clear();  c4->cd();
  hinthybchb->SetMarkerColor(kBlue); 
  hinthybchb->SetMarkerSize(5); 
  hinthybchb->SetStats(0);

  hinthybb0->SetMarkerColor(kRed); 
  hinthybb0->SetMarkerSize(5); 
  hinthybb0->SetStats(0);

  hinthybb0->Draw();
  hinthybchb->DrawCopy("same");

  leg2 = new TLegend(0.45,0.5,0.69,0.59);
  leg2->SetFillStyle(0); leg2->SetBorderSize(0.); leg2->SetTextSize(0.05); 
  leg2->SetFillColor(0); 
  legge2 = leg2->AddEntry(hinthybb0, "B0", "p"); 
  legge2 = leg2->AddEntry(hinthybchb, "B", "p"); 
  leg2->Draw();
  //  c1->SaveAs("integralsselex_mxgen.eps");
  c4->SaveAs("integralshy_supselex_mxgen.eps");
  c4->Clear();  c4->cd();
  hintnreschb->SetStats(0);
  hintnreschb->SetMarkerSize(5); 
  hintnreschb->SetMarkerColor(kBlue); 

  hintnresb0->SetStats(0);
  hintnresb0->SetMarkerSize(5); 
  hintnresb0->SetMarkerColor(kRed); 

  hintnresb0->Draw();
  hintnreschb->DrawCopy("same");

  leg2 = new TLegend(0.45,0.5,0.69,0.59);
  leg2->SetFillStyle(0); leg2->SetBorderSize(0.); leg2->SetTextSize(0.05); 
  leg2->SetFillColor(0); 
  legge2 = leg2->AddEntry(hintnresb0, "B0", "p"); 
  legge2 = leg2->AddEntry(hintnreschb, "B", "p"); 
  leg2->Draw();
  //  c1->SaveAs("integralsselex_mxgen.eps");
  c4->SaveAs("integralsnr_supselex_mxgen.eps");

  c4->Clear();  c4->cd();
  hintreschb->SetStats(0);
  hintreschb->SetMarkerSize(5); 
  hintreschb->SetMarkerColor(kBlue); 

  hintresb0->SetStats(0);
  hintresb0->SetMarkerSize(5); 
  hintresb0->SetMarkerColor(kRed); 

  hintresb0->Draw();
  hintreschb->DrawCopy("same");

  leg2 = new TLegend(0.45,0.5,0.69,0.59);
  leg2->SetFillStyle(0); leg2->SetBorderSize(0.); leg2->SetTextSize(0.05); 
  leg2->SetFillColor(0); 
  legge2 = leg2->AddEntry(hintreschb, "B", "p"); 
  legge2 = leg2->AddEntry(hintresb0, "B0", "p"); 
  leg2->Draw();
  //  c1->SaveAs("integralsselex_mxgen.eps");
  c4->SaveAs("integralsre_supselex_mxgen.eps");

}

void TheorTool::overlap(int cut, TString dir, TString flag, int scal){

  gROOT->cd();
  char name[200];   int addcut = 0;   char cutn[200]; char flagD[200]; 
  char rename[200], nrname[200], hyname[200];
  char moname[200],moname2[200], hymoname[200]; 
  char Grename[200], Gnrname[200], Ghyname[200];
  char Gmoname[200], Gmoname2[200], Ghymoname[200]; 
  sprintf(cutn,"%s","sc");

  tl.SetNDC(kTRUE);   tl.SetTextSizePixels(50.1);
  char line[200];  char fchar[100];
  gROOT->cd();gStyle->SetOptStat(0);
  Float_t ranges[21] = {0.0,  1.2, 1.264, 1.341, 1.414, 1.483, 1.549, 1.612, 1.673, 1.732, 1.788, 1.843, 1.897, 1.949, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0};
  //Loop for resonant, non reso, data

  sprintf(fchar,"%s","nu");    
  if(scal)   sprintf(fchar,"%s","ch");    
  sprintf(moname,"%s%s%s%s%s","Mx",cutn,"DA",fchar,"da");
  sprintf(moname2,"%s%s%s%s%s","Mx",cutn,"NMBDA",fchar,"da");
  sprintf(rename,"%s%s%s%s%s","Mx",cutn,"MC",fchar,"re");
  sprintf(nrname,"%s%s%s%s%s","Mx",cutn,"MC",fchar,"nr");
  sprintf(hyname,"%s%s%s","Mx",cutn,"MChy");
  sprintf(hymoname,"%s%s%s","Mx",cutn,"MChyMO");
  sprintf(Ghyname,"%s%s%s","Mx",cutn,"GMChy");
  sprintf(Ghymoname,"%s%s%s","Mx",cutn,"GMChyMO");

  TH1D *hmo = (TH1D*)gDirectory->Get(moname);//My Model
  TH1D *hmo2;
  if(cut) hmo2 = (TH1D*)gDirectory->Get(moname2);//My Model NoMbaRew
  TH1D *hre = (TH1D*)gDirectory->Get(rename);//Resonant
  TH1D *hnr = (TH1D*)gDirectory->Get(nrname);//Nonresonant
  TH1D *hnrsub = (TH1D*)gDirectory->Get("MxNRSub");//Nonreso-Sub
  TH1D *hhy = new TH1D(hyname, "",20,ranges); //Hybrid (mine/def)
  TH1D *hhymo = new TH1D(hymoname, "",20,ranges); //Hybrid (mine/def)

  double BRexcltmp;
  BRexcltmp = BRexclB0;
  if(scal)    BRexcltmp = BRexclBch;
  cout<<BRexcltmp<<" vales for :: mode "<<scal<<" "<<endl;

  //Entry based quantities. (bre)
  double RefVal = hmo->Integral();  //Reference value (our model)
  double ResoVal = hre->Integral(); //Reso entries 
  double NonResoVal = hnr->Integral(); //NonReso entries 
  /*
    [gen]                            [bre]                            
    DAnuda :: 2134.44 69.7026	     AllscDAnuda :: 7220.15 119.523   
    MCnure :: 1922.31 62.4104	     AllscMCnure :: 3313.1 76.968     
    MCnunr :: 1406.72 48.7989	     AllscMCnunr :: 2685.87 73.7469   
    DAchda :: 3089.67 80.1119	     AllscDAchda :: 7703.09 122.244   
    MCchre :: 2527.13 72.35	     AllscMCchre :: 3544.89 82.7385   
    MCchnr :: 2127.27 67.2583	     AllscMCchnr :: 3032.1 77.4689    

                (14923 / 12574)        (5224  /  7983) 
    Genfactor:: (Nonres / Nres)_cock / (Nonres / Nres)_gen
    Genfactor:: 1.8136
  */

  double Genfactor = 1.74;
  //Genfactor =1;

  cout<<"******************** "<< Genfactor <<" ******************"<<endl;
  sprintf(flagD, "%s%s%s",flag.Data(),"dum",fchar);
  sprintf(name, "%s%s%s%s", dir.Data(),"/",flagD,"model.dat");
  
  ofstream outFile(name,ios::app);
  ofstream outFileW("sysWeights.dat",ios::app);
  //bre
  //  outFile<<"Reso Scale factor:: "<< ResoRefScalFac<<endl;
  setHist(hre,1,0,1,1,3005);   
  setHist(hnr,2,0,1,1,3004);   
  setHist(hmo,4,0,1,1,3004);   
  if(cut)  setHist(hmo2,2,0,1,1,3005);   
  setHist(hhy,3,0,1,1,3005);   
  setHist(hhymo,6,0,1,1,3005);   
  setHist(hnrsub,5,0,1,1,3005);   


  //My model
  c6 = new TCanvas("c6", "c6",1200,900); c6->Clear();
  c6->cd();
  //  sprintf(line,"%s","re/nrMC ");
  TLegend *leg2; 
  TLegendEntry *legge2; 
  hmo->Draw("h");   
  if(cut)  hmo2->Draw("sameh");   
  leg2 = new TLegend(0.55,0.7,0.99,0.89);
  leg2->SetFillStyle(0); leg2->SetBorderSize(0.); leg2->SetTextSize(0.05); 
  leg2->SetFillColor(0); 
  legge2 = leg2->AddEntry(hmo, "Non reso model", "f"); 
  if(cut)  legge2 = leg2->AddEntry(hmo2, "Non reso model (No mba rew)", "f"); 
  leg2->Draw(); gPad->Update(); 
  //Definition of FLAG
  sprintf(flagD, "%s%s",flag.Data(),"NOSC");

  sprintf(name, "%s%s%s%s%s", dir.Data(),"/",flagD,fchar,"model.ps");
  c6->Print(name);  

  //Comparison model/reso before and after rescaling
  c6 = new TCanvas("c6", "c6",1200,900); c6->Clear(); c6->cd();
  hmo->Draw("h");  hre->Draw("sameh");
  leg2 = new TLegend(0.55,0.7,0.99,0.89);
  setLegend(leg2,0,0.,0.05,0,hmo,hre,"Non reso model","Resonant","f");
  leg2->Draw(); gPad->Update();
  sprintf(name, "%s%s%s%s%s", dir.Data(),"/",flagD,fchar,"model_vs_resoNO.ps");
  c6->Print(name);  

  //Comparison model/reso before and after rescaling
  c6 = new TCanvas("c6", "c6",1200,900); c6->Clear(); c6->cd();
  hmo->Draw("h");  hnr->Draw("sameh");
  leg2 = new TLegend(0.55,0.7,0.99,0.89);
  setLegend(leg2,0,0.,0.05,0,hmo,hnr,"Non reso model","Non Resonant","f");
  leg2->Draw(); gPad->Update();
  sprintf(name, "%s%s%s%s%s", dir.Data(),"/",flagD,fchar,"model_vs_nonresoNO.ps");
  c6->Print(name);  

  //My Hybrid:
  hhy->Add(hre, hmo, 1., 1.); //Resonant already scaled
  c6 = new TCanvas("c6", "c6",1200,900); c6->Clear(); c6->cd();
  hhy->Draw("h");   hmo->Draw("sameh");   
  leg2 = new TLegend(0.55,0.7,0.99,0.89);
  setLegend(leg2,0,0.,0.05,0,hmo,hhy,"Model","Hybrid","f");
  leg2->Draw(); gPad->Update(); 

  sprintf(name, "%s%s%s%s%s", dir.Data(),"/",flagD,fchar,"model_vs_mymod.ps");
  c6->Print(name);  

  //Difference btw Model and Pure non-resonant rescaled
  double nnr[21], nhy[21], nre[21], nmo[21], nmonw[21];
  double nnrexp;
  outFileW<<"bin0"<<"    "<<1<<endl;
  double FinalWe[21]; double tmpval; double nnrMe;
  double ResoRefScalFac;
  for (int ibin = 1; ibin<21; ibin++) {
    nmo[ibin] = hmo->GetBinContent(ibin);
    nre[ibin] = hre->GetBinContent(ibin);
    
    ResoRefScalFac = (RefVal/ResoVal)*(BRexcltmp/fittedBR);
    //    cout<<"Checks:: "<<nmo[ibin]<< " "<<nre[ibin]<< " "<<ResoRefScalFac<<endl;
    nnrexp = nmo[ibin] - nre[ibin] * ResoRefScalFac;
    nnrMe =  nmo[ibin] * ResoRefScalFac;

    hnrsub->SetBinContent(ibin,nnrexp);
    if(nnrMe==0) {
      tmpval = 1;
      if(nnrexp != 0){ 
	cout<<"WARNING: nnrexp/nnrMe --> infinity "<<ibin<<endl;
      }
    } else {
      if(nnrexp/nnrMe < 0) {
	tmpval = 1;
      } else if (nnrexp/nnrMe > 1000) {
	tmpval = 1;
      } else {
	tmpval = nnrexp/nnrMe;
      }
    }
    tmpval = tmpval * Genfactor;
    outFileW<<"bin"<<ibin<<"    "<<tmpval<<endl;
    FinalWe[ibin] = tmpval;
  }
  
  double NonResoRescVal = hnrsub->Integral(); //Non Reso Sub entries 
  double ResoRescVal = hre->Integral(); //Reso Resc entries 
  outFile<<"Refer. Reso. Non Reso Sum"<<RefVal<<" "<<ResoRescVal<<" "<<NonResoRescVal<<" "<<NonResoRescVal+ResoRescVal<<endl;
  outFile<<"Reso/Non reso fractions"<<ResoRescVal/RefVal<<" "<<NonResoRescVal/RefVal<<endl;

  c6->Clear(); c6->cd();
  hmo->Draw("h");  
  hnrsub->Draw("sameh");
  leg2 = new TLegend(0.55,0.7,0.99,0.89);
  setLegend(leg2,0,0.,0.05,0,hnrsub,hmo,"Non reso Sub","Non reso","f");
  leg2->Draw(); gPad->Update(); 
  sprintf(name, "%s%s%s%s%s", dir.Data(),"/",flagD,fchar,"model_vs_nonres.ps");
  c6->Print(name);  
}


int TheorTool::hist(double mx){

  // Mx categories
  if(mx<=histbins[0]) return 0;
  if(mx>histbins[0] && mx<=histbins[1]) return 1;
  if(mx>histbins[1] && mx<=histbins[2]) return 2;
  if(mx>histbins[2] && mx<=histbins[3]) return 3;
  if(mx>histbins[3] && mx<=histbins[4]) return 4;
  if(mx>histbins[4] && mx<=histbins[5]) return 5;
  if(mx>histbins[5] && mx<=histbins[6]) return 6;
  if(mx>histbins[6] && mx<=histbins[7]) return 7;
  if(mx>histbins[7] && mx<=histbins[8]) return 8;
  if(mx>histbins[8] && mx<=histbins[9]) return 9;
  if(mx>histbins[9] && mx<=histbins[10]) return 10;
  if(mx>histbins[10] && mx<=histbins[11]) return 11;
  if(mx>histbins[11] && mx<=histbins[12]) return 12;
  if(mx>histbins[12] && mx<=histbins[13]) return 13;
  if(mx>histbins[13] && mx<=histbins[14]) return 14;
  if(mx>histbins[14] && mx<=histbins[15]) return 15;
  if(mx>histbins[15] && mx<=histbins[16]) return 16;
  if(mx>histbins[16] && mx<=histbins[17]) return 17;
  if(mx>histbins[17] && mx<=histbins[18]) return 18;
  if(mx>histbins[18] && mx<=histbins[19]) return 19;
  if(mx>histbins[19] && mx<=histbins[20]) return 20;
  if(mx>histbins[20]) {
    //Mx >5 ??
    return 21;
  }
}
// ----------------------------------------------------
// Mb and a rewighting functions
// ----------------------------------------------------

void TheorTool::initMba(double dmb)
{
  //A0 is taken in one sigma range
  a1 = A0 - 0.91;
  a2 = A0 + 2.31;
  //Function of the Quark mass
  mb1 = BQMASS - dmb;
  mb2 = BQMASS + dmb;
  Nold = Nnew1 = Nnew2 = Nnew3 = Nnew4 = Nsys = 0;

  // -- Compute reweighting ratios
  Int_t nIntegral = 10000;
  const double kmin = -5.;
  const double kmax =  1.;
  for (int i = 0; i < nIntegral; ++i) {
    double kplus = (i+0.5)/((double)nIntegral)*(kmax-kmin)+kmin;
    Nold += (kmax-kmin)/((double)nIntegral)*fermi(kplus, BQMASS, A0);
    Nsys += (kmax-kmin)/((double)nIntegral)*fermi(kplus, rnd, A0);
    Nnew1 += (kmax-kmin)/((double)nIntegral)*fermi(kplus, mb1, A0);
    Nnew2 += (kmax-kmin)/((double)nIntegral)*fermi(kplus, mb2, A0);
    //    Nnew3 += (kmax-kmin)/((double)nIntegral)*fermi(kplus, rnd, a1);
    //    Nnew4 += (kmax-kmin)/((double)nIntegral)*fermi(kplus, rnd, a2);
    Nnew3 += (kmax-kmin)/((double)nIntegral)*fermi(kplus, rnd, a1);
    Nnew4 += (kmax-kmin)/((double)nIntegral)*fermi(kplus, rnd, a2);
  }
  //  cout<<"YAM "<<Nold<<Nnew1<<Nnew2<<Nnew3<<Nnew4<<endl;
}

// ----------------------------------------------------
// Mb and a rewighting functions
// ----------------------------------------------------

double TheorTool::reweightMba(double fKp, double mxHG, double mxHF, int sysWh)
{

  //Function of the Quark mass
  gROOT->cd();

  double w81, w82, w83, w84, w85;
  double wsys =1;
  double kpold = fKp; 
  double kpsys  = 5.279 - rnd - (5.279 - BQMASS - fKp);
  double kpnew1 = 5.279 - mb1 - (5.279 - BQMASS - fKp);
  double kpnew2 = 5.279 - mb2 - (5.279 - BQMASS - fKp);
  double kpnew3 = fKp;
  double kpnew4 = fKp;
  
  double fold = fermi(kpold, BQMASS, A0)/Nold; 
  double fsys = fermi(kpsys, rnd, A0)/Nsys;
  double fnew1 = fermi(kpnew1, mb1, A0)/Nnew1;
  double fnew2 = fermi(kpnew2, mb2, A0)/Nnew2;
  double fnew3 = fermi(kpsys, rnd, a1)/Nnew3;
  double fnew4 = fermi(kpsys, rnd, a2)/Nnew4;

  w81 = fnew1/fold; 
  w82 = fnew2/fold; 
  w83 = fnew3/fold; 
  w84 = fnew4/fold;
  w85 = fsys/fold;

  ((TH1D*)gDirectory->Get("nw1"))->Fill(w81);  
  ((TH1D*)gDirectory->Get("nw2"))->Fill(w82);         
  ((TH1D*)gDirectory->Get("nw3"))->Fill(w83);         
  ((TH1D*)gDirectory->Get("nw4"))->Fill(w84);
  ((TH1D*)gDirectory->Get("nk0"))->Fill(kpold, 1.);   
  ((TH1D*)gDirectory->Get("nk1"))->Fill(kpnew1, w81);   
  ((TH1D*)gDirectory->Get("nk2"))->Fill(kpnew2, w82);  
  ((TH1D*)gDirectory->Get("nk3"))->Fill(kpnew3, w83);  
  ((TH1D*)gDirectory->Get("nk4"))->Fill(kpnew4, w84);
  ((TH1D*)gDirectory->Get("nh0"))->Fill(mxHG, 1.); 
  ((TH1D*)gDirectory->Get("nh1"))->Fill(mxHG, w81);
  ((TH1D*)gDirectory->Get("nh2"))->Fill(mxHG, w82);
  ((TH1D*)gDirectory->Get("nh3"))->Fill(mxHG, w83);
  ((TH1D*)gDirectory->Get("nh4"))->Fill(mxHG, w84);

  if (mxHF > 0.01) {
    ((TH1D*)gDirectory->Get("nH0"))->Fill(mxHF, 1.); 
    ((TH1D*)gDirectory->Get("nH1"))->Fill(mxHF, w81);
    ((TH1D*)gDirectory->Get("nH2"))->Fill(mxHF, w82); 
    ((TH1D*)gDirectory->Get("nH3"))->Fill(mxHF, w83); 
    ((TH1D*)gDirectory->Get("nH4"))->Fill(mxHF, w84);
  }
  wsys = w85;
  if(sysWh == 2) {
    wsys = w84;
  } else if (sysWh == 3) {
    wsys = w83;
  }
  return wsys;
}

void TheorTool::computeWMba(double dmb) {
  gROOT->cd();
  char line[200];
	     			         
  TH1D *h0  = (TH1D*)gDirectory->Get("nh0");     TH1D *H0  = (TH1D*)gDirectory->Get("nH0");  
  TH1D *h1  = (TH1D*)gDirectory->Get("nh1");     TH1D *H1  = (TH1D*)gDirectory->Get("nH1");  
  TH1D *h2  = (TH1D*)gDirectory->Get("nh2");     TH1D *H2  = (TH1D*)gDirectory->Get("nH2");  
  TH1D *h3  = (TH1D*)gDirectory->Get("nh3");     TH1D *H3  = (TH1D*)gDirectory->Get("nH3");  
  TH1D *h4  = (TH1D*)gDirectory->Get("nh4");     TH1D *H4  = (TH1D*)gDirectory->Get("nH4");  
  TH1D *i0  = (TH1D*)gDirectory->Get("ni0");     TH1D *I0  = (TH1D*)gDirectory->Get("nI0");  
  TH1D *i1  = (TH1D*)gDirectory->Get("ni1");     TH1D *I1  = (TH1D*)gDirectory->Get("nI1");  
  TH1D *i2  = (TH1D*)gDirectory->Get("ni2");     TH1D *I2  = (TH1D*)gDirectory->Get("nI2");  
  TH1D *i3  = (TH1D*)gDirectory->Get("ni3");     TH1D *I3  = (TH1D*)gDirectory->Get("nI3");  
  TH1D *i4  = (TH1D*)gDirectory->Get("ni4");     TH1D *I4  = (TH1D*)gDirectory->Get("nI4");  
  TH1D *m1  = (TH1D*)gDirectory->Get("nm1");     TH1D *M1  = (TH1D*)gDirectory->Get("nM1");  
  TH1D *m2  = (TH1D*)gDirectory->Get("nm2");     TH1D *M2  = (TH1D*)gDirectory->Get("nM2");  
  TH1D *m3  = (TH1D*)gDirectory->Get("nm3");     TH1D *M3  = (TH1D*)gDirectory->Get("nM3");  
  TH1D *m4  = (TH1D*)gDirectory->Get("nm4");     TH1D *M4  = (TH1D*)gDirectory->Get("nM4");  
  TH1D *e0  = (TH1D*)gDirectory->Get("ne0");     TH1D *E0  = (TH1D*)gDirectory->Get("nE0");  
  TH1D *e02 = (TH1D*)gDirectory->Get("ne02");    TH1D *E02 = (TH1D*)gDirectory->Get("nE02"); 
  TH1D *e1  = (TH1D*)gDirectory->Get("ne1");     TH1D *E1  = (TH1D*)gDirectory->Get("nE1");  
  TH1D *e2  = (TH1D*)gDirectory->Get("ne2");     TH1D *E2  = (TH1D*)gDirectory->Get("nE2");  
  TH1D *e3  = (TH1D*)gDirectory->Get("ne3");     TH1D *E3  = (TH1D*)gDirectory->Get("nE3");  
  TH1D *e4  = (TH1D*)gDirectory->Get("ne4");     TH1D *E4  = (TH1D*)gDirectory->Get("nE4");  


  h0->Scale(1./h0->GetSumOfWeights());
  h1->Scale(1./h1->GetSumOfWeights());
  h2->Scale(1./h2->GetSumOfWeights());
  h3->Scale(1./h3->GetSumOfWeights());
  h4->Scale(1./h4->GetSumOfWeights());
  H0->Scale(1./H0->GetSumOfWeights());
  H1->Scale(1./H1->GetSumOfWeights());
  H2->Scale(1./H2->GetSumOfWeights());
  H3->Scale(1./H3->GetSumOfWeights());
  H4->Scale(1./H4->GetSumOfWeights());

  for (int i=1; i <= h0->GetNbinsX(); ++i) {
    i0->SetBinContent(i, h0->Integral(0, i));
    i1->SetBinContent(i, h1->Integral(0, i));
    i2->SetBinContent(i, h2->Integral(0, i));
    m1->SetBinContent(i, TMath::Max(i1->GetBinContent(i), i2->GetBinContent(i))); 
    m2->SetBinContent(i, TMath::Min(i1->GetBinContent(i), i2->GetBinContent(i))); 
 
    i3->SetBinContent(i, h3->Integral(0, i));
    i4->SetBinContent(i, h4->Integral(0, i));
    m3->SetBinContent(i, TMath::Min(i3->GetBinContent(i), i4->GetBinContent(i))); 
    m4->SetBinContent(i, TMath::Max(i3->GetBinContent(i), i4->GetBinContent(i))); 
    cout<<"Values :: "<<h0->Integral(0, i)<<" "<<h1->Integral(0, i)<<" "<<h2->Integral(0, i)<<" "<<h3->Integral(0, i)<<" "<<h4->Integral(0, i)<<" "<<endl;


    if (i0->GetBinContent(i)) e1->SetBinContent(i, TMath::Abs((i2->GetBinContent(i) - i0->GetBinContent(i)))/i0->GetBinContent(i));
    if (i0->GetBinContent(i)) e2->SetBinContent(i, TMath::Abs((i1->GetBinContent(i) - i0->GetBinContent(i)))/i0->GetBinContent(i));
    if (i0->GetBinContent(i)) e3->SetBinContent(i, TMath::Abs((i3->GetBinContent(i) - i0->GetBinContent(i)))/i0->GetBinContent(i));
    if (i0->GetBinContent(i)) e4->SetBinContent(i, TMath::Abs((i4->GetBinContent(i) - i0->GetBinContent(i)))/i0->GetBinContent(i));
    e0->SetBinContent(i, 0.5*(e1->GetBinContent(i) + e2->GetBinContent(i)));
    e0->SetBinContent(i, TMath::Max(e1->GetBinContent(i),e2->GetBinContent(i)));
    e02->SetBinContent(i, 0.5*(e3->GetBinContent(i) + e4->GetBinContent(i)));
    e02->SetBinContent(i, TMath::Max(e3->GetBinContent(i),e4->GetBinContent(i)));

    I0->SetBinContent(i, H0->Integral(0, i));
    I1->SetBinContent(i, H1->Integral(0, i));
    I2->SetBinContent(i, H2->Integral(0, i));
    M1->SetBinContent(i, TMath::Max(I1->GetBinContent(i), I2->GetBinContent(i))); 
    M2->SetBinContent(i, TMath::Min(I1->GetBinContent(i), I2->GetBinContent(i))); 

    I3->SetBinContent(i, H3->Integral(0, i));
    I4->SetBinContent(i, H4->Integral(0, i));
    M3->SetBinContent(i, TMath::Min(I3->GetBinContent(i), I4->GetBinContent(i))); 
    M4->SetBinContent(i, TMath::Max(I3->GetBinContent(i), I4->GetBinContent(i))); 

    if (I0->GetBinContent(i)) E1->SetBinContent(i, TMath::Abs((I1->GetBinContent(i) - I0->GetBinContent(i)))/I0->GetBinContent(i));
    if (I0->GetBinContent(i)) E2->SetBinContent(i, TMath::Abs((I2->GetBinContent(i) - I0->GetBinContent(i)))/I0->GetBinContent(i));
    if (I0->GetBinContent(i)) E3->SetBinContent(i, TMath::Abs((I3->GetBinContent(i) - I0->GetBinContent(i)))/I0->GetBinContent(i));
    if (I0->GetBinContent(i)) E4->SetBinContent(i, TMath::Abs((I4->GetBinContent(i) - I0->GetBinContent(i)))/I0->GetBinContent(i));
    E0->SetBinContent(i, 0.5*(E1->GetBinContent(i) + E2->GetBinContent(i)));
    E0->SetBinContent(i, TMath::Max(E1->GetBinContent(i),E2->GetBinContent(i)));
    E02->SetBinContent(i, 0.5*(E3->GetBinContent(i) + E4->GetBinContent(i)));
    E02->SetBinContent(i, TMath::Max(E3->GetBinContent(i),E4->GetBinContent(i)));
  }

  gStyle->SetOptStat(0);  gStyle->SetOptTitle(0);
  c0 = new TCanvas("c1"," ",200,10,800,1100); 
  c0->Clear(); c0->Divide(2,3); 

  // -- Mxfit
  c0->cd(1); shrinkPad(0.12, 0.10);
  setTitles(H0, "M_{X}^{fit} [GeV]", " ", 0.06);
  H0->SetMaximum(1.4*H0->GetMaximum());
  H0->SetLineWidth(3);
  H0->SetTitleOffset(.5, "X");
  H0->Draw();
  setFilledHist(H1, kBlue, kBlue, 3005);
  setFilledHist(H2, kRed, kRed, 3004);
  H1->Draw("same");
  H2->Draw("same");
  TLegend *leg;
  TLegendEntry *legge;
  leg = new TLegend(0.25,0.7,0.8,0.89);
  leg->SetFillStyle(0); leg->SetBorderSize(0.); leg->SetTextSize(0.05);  leg->SetFillColor(0); 
  legge = leg->AddEntry(H0, Form("m_{} = 4.800, a_{} = 1.290"), "l"); legge->SetTextColor(kBlack);
  sprintf(line, "m_{b} = %4.3f, a = %4.3f", mb1, A0);
  legge = leg->AddEntry(H1, line, "l"); legge->SetTextColor(kBlue);
  sprintf(line, "m_{b} = %4.3f, a = %4.3f", mb2, A0);
  legge = leg->AddEntry(H2, line, "l"); legge->SetTextColor(kRed);
  leg->Draw();

  // -- Mxgen
  c0->cd(2); shrinkPad(0.12, 0.10);
  setTitles(h0, "M_{X}^{gen} [GeV]", " ", 0.06);
  h0->SetMaximum(1.4*h0->GetMaximum());
  h0->SetLineWidth(3);
  h0->Draw();
  setFilledHist(h1, kBlue, kBlue, 3005);
  setFilledHist(h2, kRed, kRed, 3004);
  h1->Draw("same");
  h2->Draw("same");

  
  // -- Mxfit
  c0->cd(3); shrinkPad(0.12, 0.10);
  setTitles(H0, "M_{X}^{fit} [GeV]", " ", 0.06);
  H0->SetLineWidth(3);
  H0->Draw();
  setFilledHist(H3, kBlue, kBlue, 3005);
  setFilledHist(H4, kRed, kRed, 3004);
  H3->Draw("same");
  H4->Draw("same");

  leg = new TLegend(0.25,0.7,0.8,0.89);
  leg->SetFillStyle(0); leg->SetBorderSize(0.); leg->SetTextSize(0.05);  leg->SetFillColor(0); 
  legge = leg->AddEntry(H0, Form("m_{} = 4.800,a_{} = 1.290"), "l"); legge->SetTextColor(kBlack);
  sprintf(line, "m_{b} = %4.3f, a = %4.3f", BQMASS, a1);
  legge = leg->AddEntry(H3, line, "l"); legge->SetTextColor(kBlue);
  sprintf(line, "m_{b} = %4.3f, a = %4.3f", BQMASS, a2);
  legge = leg->AddEntry(H4, line, "l"); legge->SetTextColor(kRed);
  leg->Draw();

  // -- Mxgen
  c0->cd(4); shrinkPad(0.12, 0.10);
  setTitles(h0, "M_{X}^{gen} [GeV]", " ", 0.06);
  h0->SetLineWidth(3);
  h0->Draw();
  setFilledHist(h3, kBlue, kBlue, 3005);
  setFilledHist(h4, kRed, kRed, 3004);
  h3->Draw("same");
  h4->Draw("same");

  
  // -- integral of Mxfit
  c0->cd(5); shrinkPad(0.12, 0.10); gPad->SetGridx(1); gPad->SetGridy(1);
  I0->SetLineWidth(3);
  setTitles(I2, "M_{X}^{fit} [GeV]", "f(u)", 0.06);

  setFilledHist(M1, kRed, kYellow, 1000); 
  setFilledHist(M4, kBlue, kBlue, 1000); 
  setFilledHist(M3, kBlue, kYellow, 1000); 
  setFilledHist(M2, kRed, 10, 1000); 
  setTitles(M1, "M_{X}^{fit} [GeV]", " ", 0.06);

  M1->Smooth(200); M1->Draw("c");
  M4->Smooth(200); M4->Draw("csame");
  M3->Smooth(200); M3->Draw("csame");
  M2->Smooth(200); M2->Draw("cSAME");
  I0->Smooth(200); I0->Draw("csame");
  M1->Draw("axissame");


  leg = new TLegend(0.45,0.2,0.8,0.5);
  leg->SetFillStyle(0); leg->SetBorderSize(0.); leg->SetTextSize(0.05);  leg->SetFillColor(0); 
  sprintf(line, "m_{b} = %4.3f ^{+%4.3f}_{-%4.3f}", BQMASS, mb2-BQMASS, BQMASS-mb1);
  legge = leg->AddEntry(M1, line, "f"); 
  sprintf(line, "a = %4.3f^{+%4.3f}_{-%4.3f}", A0, a2-A0, A0-a1);
  legge = leg->AddEntry(M4, line, "f"); 
  leg->Draw();


  // -- integral of Mxgen
  c0->cd(6); shrinkPad(0.12, 0.10); gPad->SetGridx(1); gPad->SetGridy(1);
  i0->SetLineWidth(3);
  setTitles(i2, "M_{X}^{gen} [GeV]", "f(u)", 0.06);

  setFilledHist(m1, kRed, kYellow, 1000); 
  setFilledHist(m4, kBlue, kBlue, 1000); 
  setFilledHist(m3, kBlue, kYellow, 1000); 
  setFilledHist(m2, kRed, 10, 1000); 
  setTitles(m1, "M_{X}^{gen} [GeV]", " ", 0.06);

  m1->Smooth(200); m1->Draw("c");
  m4->Smooth(200); m4->Draw("csame");
  m3->Smooth(200); m3->Draw("csame");
  m2->Smooth(200); m2->Draw("cSAME");
  i0->Smooth(200); i0->Draw("csame");
  m1->Draw("axissame");

  sprintf(line, "sfDependency-%4.2f-allcuts.eps", dmb); 
  c0->Print(line);

  c1 = new TCanvas("c2"," ",200,10,1200,900); 
  c1->Clear();
  c1->Divide(2,2);

  c1->cd(1);  shrinkPad(0.15); gPad->SetGridx(1); gPad->SetGridy(1);
  gStyle->SetOptFit(0);
  setTitles(E0, "M_{X}^{fit} [GeV]", "rel. Error", 0.06);
  E0->SetMinimum(0.); E0->SetMaximum(0.5); E0->SetAxisRange(0.5, 4.0);
  E0->Draw("");
  E0->Fit("pol4", "r", "", 0.5, 3.0);
  tl.DrawLatex(0.6, 0.62, Form("M_{X}^{fit} < 1.55 GeV:"));
  tl.DrawLatex(0.6, 0.53, Form("#Delta R/ R(mb) = %4.3f", E0->GetBinContent(E0->FindBin(1.55))));


  c1->cd(2); shrinkPad(0.15); gPad->SetGridx(1); gPad->SetGridy(1);
  //  e1->SetLineColor(kBlue);  e1->Draw("");
  //  e2->SetLineColor(kRed);   e2->Draw("same");
  setTitles(e0, "M_{X}^{gen} [GeV]", "rel. Error", 0.06);
  e0->SetMinimum(0.); e0->SetMaximum(0.5); e0->SetAxisRange(0.5, 4.0);
  e0->Draw("");
  e0->Fit("pol4", "r", "", 0.5, 3.0);
  tl.DrawLatex(0.6, 0.62, Form("M_{X}^{gen} < 1.55 GeV:"));
  tl.DrawLatex(0.6, 0.53, Form("#Delta R/ R(mb) = %4.3f", e0->GetBinContent(e0->FindBin(1.55))));

  c1->cd(3);  shrinkPad(0.15); gPad->SetGridx(1); gPad->SetGridy(1);
  gStyle->SetOptFit(0);
  setTitles(E02, "M_{X}^{fit} [GeV]", "rel. Error", 0.06);
  E02->SetMinimum(0.); E02->SetMaximum(0.5); E02->SetAxisRange(0.5, 4.0);
  E02->Draw("");
  E02->Fit("pol4", "r", "", 0.5, 3.0);
  tl.DrawLatex(0.6, 0.62, Form("M_{X}^{fit} < 1.55 GeV:"));
  tl.DrawLatex(0.6, 0.53, Form("#Delta R/ R(a) = %4.3f", E02->GetBinContent(E02->FindBin(1.55))));


  c1->cd(4); shrinkPad(0.15); gPad->SetGridx(1); gPad->SetGridy(1);
  //  e1->SetLineColor(kBlue);  e1->Draw("");
  //  e2->SetLineColor(kRed);   e2->Draw("same");
  setTitles(e02, "M_{X}^{gen} [GeV]", "rel. Error", 0.06);
  e02->SetMinimum(0.); e02->SetMaximum(0.5); e02->SetAxisRange(0.5, 4.0);
  e02->Draw("");
  e02->Fit("pol4", "r", "", 0.5, 3.0);
  tl.DrawLatex(0.6, 0.62, Form("M_{X}^{gen} < 1.55 GeV:"));
  tl.DrawLatex(0.6, 0.53, Form("#Delta R/ R (a) = %4.3f", e02->GetBinContent(e02->FindBin(1.55))));

  sprintf(line, "effscan-allcuts.txt");
  
  ofstream OUT(line);
  E0->Smooth(50);   e0->Smooth(50);
  for (int ij = 20; ij < 40; ++ij) {
    /*
    cout << "mx(fit) <  " << E0->GetBinLowEdge(ij) << ": " 
	 << Form("  %4.3f", I0->GetBinContent(ij)) << "  " 
	 << Form(" +/- %4.3f", E0->GetBinContent(ij)) << "  " 

	 << "mx(gen) <  " << e0->GetBinLowEdge(ij) << ": " 
	 << Form("  %4.3f", i0->GetBinContent(ij))
	 << Form(" +/- %4.3f", e0->GetBinContent(ij))
	 << endl;
    */
    OUT  << "mx(fit) <  " << E0->GetBinLowEdge(ij) << ": " 
	 << Form("  %4.3f", I0->GetBinContent(ij)) << "  " 
	 << Form(" +/- %4.3f", E0->GetBinContent(ij)) << "  " 

	 << "mx(gen) <  " << e0->GetBinLowEdge(ij) << ": " 
	 << Form("  %4.3f", i0->GetBinContent(ij))
	 << Form(" +/- %4.3f", e0->GetBinContent(ij))
	 << endl;

  }

  OUT.close();

  sprintf(line, "sfError-%4.2f-allcuts.eps", dmb); 

  c1->Print(line);
}

  // -- Shape function
double TheorTool::fermi(double kp, double m, double a) {
  double BMASS   = 5.2792;
  double x = kp/(BMASS - m);
  if ((kp>-m) && (x <= 1.)) {
    return TMath::Power((1-x), a) * TMath::Exp((1+a)*x); 
  } 
  return 0.;
}


// ----------------------------------------------------
// Utility functions
// ----------------------------------------------------

void TheorTool::sighisto(double&signal, double& signalErr, TH1D *histo, double &resmean, double &ressigma, double &resalpha, double &resn, int fixpar, double mean, double sigma, double alpha,  double n, double argus){

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

// ----------------------------------------------------------------------
void TheorTool::setHist(TH1 *h, Int_t color, Int_t symbol, Double_t size, Double_t width, Int_t style) {
  h->SetLineWidth(width);
  if(symbol != 0) {
    h->SetMarkerColor(color); h->SetMarkerStyle(symbol); 
    h->SetMarkerSize(size); 
  } else if (style != 0) {
    h->SetFillStyle(style);
    h->SetFillColor(color);
  }
  h->SetLineColor(color); 
}
// ----------------------------------------------------------------------

float TheorTool::randomized(float mean,float sig){
  if(theseed == 0 ) return mean;// no smearing required
  if(sig == 0 ) return mean;// just the mean
  float result=rndm.Gaus(mean,sig);
  if(result < 0) result = randomized(mean,sig); // if negative weight throw the dice again
  return result;
}

// ----------------------------------------------------------------------
double TheorTool::getTrueMxWeight(double thetrumx, int index) {
  //  cout<<"Entered function: "<<thetrumx<<endl;
  int thebin = hist(thetrumx);
  if(index<0) thebin = thebin + 21;
  return TrueMxWeight[thebin];
}

// ----------------------------------------------------------------------
void TheorTool::setLegend(TLegend *legend, Int_t style, Double_t bsize, Double_t tsize, Double_t fcolor, TH1 *fHi, TH1 *sHi, char * fTit, char * sTit, char * type) {
  TLegendEntry *leggent; 
  legend->SetFillStyle(style); legend->SetBorderSize(bsize);
  legend->SetTextSize(tsize);  legend->SetFillColor(fcolor); 
  leggent = legend->AddEntry(fHi,fTit,type); 
  leggent = legend->AddEntry(sHi,sTit,type); 
}
// ----------------------------------------------------------------------

// ----------------------------------------------------------------------
void TheorTool::dumpCuts() {
  cout << "====================================" << endl;
  cout << "Cut file " << fCutFile << endl; 
  cout << "------------------------------------" << endl;
  cout << "deSignal:          " << DESIGNALLO << " ... " << DESIGNALHI << endl;
  cout << "pcms:              " << CUTPCMS   << endl;
  cout << "intPurity:         " << INTPURITY << endl;
  cout << "====================================" << endl;
}

void TheorTool::initRest() {
  // -- CUTS
  DESIGNALLO = -0.1; 
  DESIGNALHI =  0.1; 
  CUTPCMS    = 1;
  INTPURITY  = 0.;
  MXMAXCUT  = 5.;
  MM2CUT  = .5;
  DFROMD  = 1;
}

// ----------------------------------------------------------------------
void TheorTool::readCuts(TString filename, int dump) {
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
    if (ok == 0)  cout << "==> recoilNtp::readCuts() Error: Don't know about variable " << CutName << endl;
  }

  if (dump == 1) dumpCuts();

}

void TheorTool::Init(TTree *tree)
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
   /*
     fChain->SetBranchAddress("GfDkmiss",&GfDkmiss);  
     fChain->SetBranchAddress("GfDkl",&GfDkl); 
     fChain->SetBranchAddress("GfDkspipi",&GfDkspipi); 
     fChain->SetBranchAddress("GfDkspiopio",&GfDkspiopio); 
     fChain->SetBranchAddress("bgcat",&bgcat);
   */
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
   fChain->SetBranchAddress("vub",&vub);
   fChain->SetBranchAddress("vcb",&vcb);
   fChain->SetBranchAddress("vxbtyp",&vxbtyp);
   fChain->SetBranchAddress("other",&other);
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
   fChain->SetBranchAddress("fBchgen",&fBchgen); 
   fChain->SetBranchAddress("nkp",&nkp);
   fChain->SetBranchAddress("nks",&nks);
   fChain->SetBranchAddress("npi0",&npi0);
   fChain->SetBranchAddress("pnu",&pnu);
   /*
     fChain->SetBranchAddress("dx",&dx);
     fChain->SetBranchAddress("dy",&dy);
     fChain->SetBranchAddress("dz",&dz);
     fChain->SetBranchAddress("s2dxx",&s2dxx);
     fChain->SetBranchAddress("s2dyy",&s2dyy);
     fChain->SetBranchAddress("s2dzz",&s2dzz);
     fChain->SetBranchAddress("s2dxy",&s2dxy);
     fChain->SetBranchAddress("s2dyz",&s2dyz);
     fChain->SetBranchAddress("s2dxz",&s2dxz);
   */
   fChain->SetBranchAddress("tnu",&tnu);
   fChain->SetBranchAddress("fnu",&fnu);
   fChain->SetBranchAddress("ENeu",&eneu);
   fChain->SetBranchAddress("EPiz",&epiz);
   fChain->SetBranchAddress("MinKMom",&kminmom);
   fChain->SetBranchAddress("MaxKMom",&kmaxmom);
   fChain->SetBranchAddress("mm2",&mm2);
   fChain->SetBranchAddress("mm2nc",&mm2nc);
   fChain->SetBranchAddress("mm2fit",&mm2fit);
   /*
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
   */
   fChain->SetBranchAddress("totweight",&totweight);
   fChain->SetBranchAddress("totweightNutMult",&totweightNutMult);
   fChain->SetBranchAddress("totweightTrkMult",&totweightTrkMult);
   Notify();
}

Bool_t TheorTool::Notify()
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
   b_fBchgen=      fChain->GetBranch("fBchgen");  
   b_GfDpi=      fChain->GetBranch("GfDpi");  
   b_GfDpiz=     fChain->GetBranch("GfDpiz");
   b_GfDk=       fChain->GetBranch("GfDk");  
   b_GfDkmiss=      fChain->GetBranch("GfDkmiss"); 
   b_GfDks=      fChain->GetBranch("GfDks"); 
   b_GfDkl=      fChain->GetBranch("GfDkl"); 
   b_GfDkspipi=      fChain->GetBranch("GfDkspipi"); 
   b_GfDkspiopio=      fChain->GetBranch("GfDkspiopio"); 
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

