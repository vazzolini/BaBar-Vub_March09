#include "BkgbrekTool.hh"
#include "recoilDSys.hh"
#include "TCut.h"
#include <TMath.h>
#include <TTree.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLegendEntry.h>
#include <TH1.h>
#include <TH2.h>
#include <TPaveStats.h>
#include <fstream.h>
#include <iostream.h>
#include <strstream.h>
#include <iomanip.h>


BkgbrekTool::BkgbrekTool(const char *var, double mi, double ma, int b, const char *var2)
{
  for(int ik = 0; ik<31; ik++){
    for(int ip = 0; ip<6; ip++){
      MA[ik][ip] = 0;
    }
  }
  //Pi         K            KS       Pi0      
  MA[0][0]=1; MA[0][2]=1; MA[0][3]=0; MA[0][1]=0; MA[0][4]=0; MA[0][5]=0;
  MA[1][0]=0; MA[1][2]=0; MA[1][3]=1; MA[1][1]=1; MA[1][4]=0; MA[1][5]=0;
  MA[2][0]=2; MA[2][2]=0; MA[2][3]=1; MA[2][1]=0; MA[2][4]=0; MA[2][5]=0;
  MA[3][0]=1; MA[3][2]=1; MA[3][3]=0; MA[3][1]=1; MA[3][4]=0; MA[3][5]=0;
  MA[4][0]=0; MA[4][2]=0; MA[4][3]=1; MA[4][1]=2; MA[4][4]=0; MA[4][5]=0;
  MA[5][0]=3; MA[5][2]=1; MA[5][3]=0; MA[5][1]=0; MA[5][4]=0; MA[5][5]=0;
  MA[6][0]=2; MA[6][2]=0; MA[6][3]=1; MA[6][1]=1; MA[6][4]=0; MA[6][5]=0;
  MA[7][0]=1; MA[7][2]=1; MA[7][3]=0; MA[7][1]=2; MA[7][4]=0; MA[7][5]=0;
  MA[8][0]=3; MA[8][2]=1; MA[8][3]=0; MA[8][1]=1; MA[8][4]=0; MA[8][5]=0;
  MA[9][0]=2; MA[9][2]=0; MA[9][3]=0; MA[9][1]=0; MA[9][4]=0; MA[9][5]=0;
  MA[10][0]=0;MA[10][2]=0;MA[10][3]=0;MA[10][1]=2;MA[10][4]=0;MA[10][5]=0;
  MA[11][0]=2;MA[11][2]=0;MA[11][3]=0;MA[11][1]=1;MA[11][4]=0;MA[11][5]=0;
  MA[12][0]=4;MA[12][2]=0;MA[12][3]=0;MA[12][1]=0;MA[12][4]=0;MA[12][5]=0;
  MA[13][0]=4;MA[13][2]=0;MA[13][3]=0;MA[13][1]=1;MA[13][4]=0;MA[13][5]=0;
  MA[14][0]=6;MA[14][2]=0;MA[14][3]=0;MA[14][1]=0;MA[14][4]=0;MA[14][5]=0;
                    		 	    
  MA[15][0]=1;MA[15][2]=0;MA[15][3]=1;MA[15][1]=0;MA[15][4]=0;MA[15][5]=0;
  MA[16][0]=2;MA[16][2]=1;MA[16][3]=0;MA[16][1]=0;MA[16][4]=0;MA[16][5]=0;
  MA[17][0]=1;MA[17][2]=0;MA[17][3]=1;MA[17][1]=1;MA[17][4]=0;MA[17][5]=0;
  MA[18][0]=2;MA[18][2]=1;MA[18][3]=0;MA[18][1]=1;MA[18][4]=0;MA[18][5]=0;
  MA[19][0]=3;MA[19][2]=0;MA[19][3]=1;MA[19][1]=0;MA[19][4]=0;MA[19][5]=0;
  MA[20][0]=4;MA[20][2]=1;MA[20][3]=0;MA[20][1]=0;MA[20][4]=0;MA[20][5]=0;
  MA[21][0]=2;MA[21][2]=1;MA[21][3]=0;MA[21][1]=2;MA[21][4]=0;MA[21][5]=0;
  MA[22][0]=3;MA[22][2]=0;MA[22][3]=1;MA[22][1]=1;MA[22][4]=0;MA[22][5]=0;
  MA[23][0]=5;MA[23][2]=0;MA[23][3]=1;MA[23][1]=0;MA[23][4]=0;MA[23][5]=0;
  MA[24][0]=4;MA[24][2]=1;MA[24][3]=0;MA[24][1]=1;MA[24][4]=0;MA[24][5]=0;
  MA[25][0]=0;MA[25][2]=1;MA[25][3]=2;MA[25][1]=0;MA[25][4]=0;MA[25][5]=0;
  MA[26][0]=1;MA[26][2]=0;MA[26][3]=0;MA[26][1]=1;MA[26][4]=0;MA[26][5]=0;
  MA[27][0]=3;MA[27][2]=0;MA[27][3]=0;MA[27][1]=0;MA[27][4]=0;MA[27][5]=0;
  MA[28][0]=3;MA[28][2]=0;MA[28][3]=0;MA[28][1]=1;MA[28][4]=0;MA[28][5]=0;
  MA[29][0]=5;MA[29][2]=0;MA[29][3]=0;MA[29][1]=0;MA[29][4]=0;MA[29][5]=0;
  MA[30][0]=5;MA[30][2]=0;MA[30][3]=0;MA[30][1]=1;MA[30][4]=0;MA[30][5]=0;
  themin = mi;
  themax = ma;
  thebins = b;
  thevar = var;
  thevar2 = var2;
  for(int id=0; id<3;id++) {
    for(int idd=0; idd<2;idd++) {
      fNall[id][idd] = 0;
    }
  }
  //  cout<<thevar<<endl;
  initRest();
}
BkgbrekTool::BkgbrekTool(const char *var, double mi, double ma, int b, const char *var2, double mi2, double ma2, int b2)
{
  for(int ik = 0; ik<31; ik++){
    for(int ip = 0; ip<6; ip++){
      MA[ik][ip] = 0;
    }
  }
  //Pi         K            KS       Pi0      
  MA[0][0]=1; MA[0][2]=1; MA[0][3]=0; MA[0][1]=0; MA[0][4]=0; MA[0][5]=0;
  MA[1][0]=0; MA[1][2]=0; MA[1][3]=1; MA[1][1]=1; MA[1][4]=0; MA[1][5]=0;
  MA[2][0]=2; MA[2][2]=0; MA[2][3]=1; MA[2][1]=0; MA[2][4]=0; MA[2][5]=0;
  MA[3][0]=1; MA[3][2]=1; MA[3][3]=0; MA[3][1]=1; MA[3][4]=0; MA[3][5]=0;
  MA[4][0]=0; MA[4][2]=0; MA[4][3]=1; MA[4][1]=2; MA[4][4]=0; MA[4][5]=0;
  MA[5][0]=3; MA[5][2]=1; MA[5][3]=0; MA[5][1]=0; MA[5][4]=0; MA[5][5]=0;
  MA[6][0]=2; MA[6][2]=0; MA[6][3]=1; MA[6][1]=1; MA[6][4]=0; MA[6][5]=0;
  MA[7][0]=1; MA[7][2]=1; MA[7][3]=0; MA[7][1]=2; MA[7][4]=0; MA[7][5]=0;
  MA[8][0]=3; MA[8][2]=1; MA[8][3]=0; MA[8][1]=1; MA[8][4]=0; MA[8][5]=0;
  MA[9][0]=2; MA[9][2]=0; MA[9][3]=0; MA[9][1]=0; MA[9][4]=0; MA[9][5]=0;
  MA[10][0]=0;MA[10][2]=0;MA[10][3]=0;MA[10][1]=2;MA[10][4]=0;MA[10][5]=0;
  MA[11][0]=2;MA[11][2]=0;MA[11][3]=0;MA[11][1]=1;MA[11][4]=0;MA[11][5]=0;
  MA[12][0]=4;MA[12][2]=0;MA[12][3]=0;MA[12][1]=0;MA[12][4]=0;MA[12][5]=0;
  MA[13][0]=4;MA[13][2]=0;MA[13][3]=0;MA[13][1]=1;MA[13][4]=0;MA[13][5]=0;
  MA[14][0]=6;MA[14][2]=0;MA[14][3]=0;MA[14][1]=0;MA[14][4]=0;MA[14][5]=0;
                    		 	    
  MA[15][0]=1;MA[15][2]=0;MA[15][3]=1;MA[15][1]=0;MA[15][4]=0;MA[15][5]=0;
  MA[16][0]=2;MA[16][2]=1;MA[16][3]=0;MA[16][1]=0;MA[16][4]=0;MA[16][5]=0;
  MA[17][0]=1;MA[17][2]=0;MA[17][3]=1;MA[17][1]=1;MA[17][4]=0;MA[17][5]=0;
  MA[18][0]=2;MA[18][2]=1;MA[18][3]=0;MA[18][1]=1;MA[18][4]=0;MA[18][5]=0;
  MA[19][0]=3;MA[19][2]=0;MA[19][3]=1;MA[19][1]=0;MA[19][4]=0;MA[19][5]=0;
  MA[20][0]=4;MA[20][2]=1;MA[20][3]=0;MA[20][1]=0;MA[20][4]=0;MA[20][5]=0;
  MA[21][0]=2;MA[21][2]=1;MA[21][3]=0;MA[21][1]=2;MA[21][4]=0;MA[21][5]=0;
  MA[22][0]=3;MA[22][2]=0;MA[22][3]=1;MA[22][1]=1;MA[22][4]=0;MA[22][5]=0;
  MA[23][0]=5;MA[23][2]=0;MA[23][3]=1;MA[23][1]=0;MA[23][4]=0;MA[23][5]=0;
  MA[24][0]=4;MA[24][2]=1;MA[24][3]=0;MA[24][1]=1;MA[24][4]=0;MA[24][5]=0;
  MA[25][0]=0;MA[25][2]=1;MA[25][3]=2;MA[25][1]=0;MA[25][4]=0;MA[25][5]=0;
  MA[26][0]=1;MA[26][2]=0;MA[26][3]=0;MA[26][1]=1;MA[26][4]=0;MA[26][5]=0;
  MA[27][0]=3;MA[27][2]=0;MA[27][3]=0;MA[27][1]=0;MA[27][4]=0;MA[27][5]=0;
  MA[28][0]=3;MA[28][2]=0;MA[28][3]=0;MA[28][1]=1;MA[28][4]=0;MA[28][5]=0;
  MA[29][0]=5;MA[29][2]=0;MA[29][3]=0;MA[29][1]=0;MA[29][4]=0;MA[29][5]=0;
  MA[30][0]=5;MA[30][2]=0;MA[30][3]=0;MA[30][1]=1;MA[30][4]=0;MA[30][5]=0;
  themin = mi;
  themax = ma;
  thebins = b;
  thevar = var ;
  themin2 = mi2;
  themax2 = ma2;
  thebins2 = b2;
  thevar2 = var2 ;
  for(int id=0; id<3;id++) {
    for(int idd=0; idd<2;idd++) {
      fNall[id][idd] = 0;
    }
  }
  //  cout<<thevar<<endl;
  initRest();
}

BkgbrekTool::~BkgbrekTool()
{
  //   if (!fChain) return;
  //   delete fChain->GetCurrentFile();
}

void BkgbrekTool::Loop(int nevents, int cat, int bsel, int semilep)
{

  gROOT->cd();
  //id = Mx category == Mx bin
  int id, id2, idflav, ich, ine, dcat, group;        
  char name[100];  char le[100];
  
  // restore B0/B+ flag
  int isbch =  bsel>=0? bsel: -(1+bsel);
  if (fChain == 0) return;
  Int_t nentries = Int_t(fChain->GetEntries());
  if( nentries > nevents) nentries = nevents;
  Int_t nbytes = 0, nb = 0;
  
  for (Int_t jentry=0; jentry<nentries;jentry++) {
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if(bsel<0)
      if( (!vub && cat==4) || (!vcb && cat==5) )continue;
    
    int flav =  lcharge + brecoflav; // charge correlation
    //    bool ksele = nkp + nks;          // fit on the depleted sample?
    int ch = xcharge + brecocharge;  // total charge
    int cutflav = 0;              // flavor category (3 = charged B, 4 = neutral B OS, 5 = neutral B SS)  
    
    if(nchg-nle<1) {ich = 1;}
    else if(nchg-nle<3){ich = 2;}
    else{ich = 3;};
    
    if(nneu80_160<1) {ine = 1;}
    else{ine = 2;};      
    if(wdeltam>-3 && brecocharge == 0) continue;
    // cuts
    bool PCMS =  pcms>CUTPCMS;
    //    bool NPI0 = nnpi0<2;
    bool NLE = nle > 0; bool NLE1 = nle == 1; bool IPUR = intpur>INTPURITY;
    bool MM2 =  mm2 < MM2CUT; bool CH = ch == 0;
    int FLAV = !(TMath::Abs(brecocharge)!=0 && (lcharge + brecoflav)!=0);
    bool BCH = TMath::Abs(brecocharge) == isbch;     
    bool MXCUT = mxhadfit < MXMAXCUT;
    if (isbch == 2) BCH = 1;
    bool Dfromd = ((GfD0Ds == DFROMD) || (GfDDs == DFROMD));
    if(DFROMD == 0) Dfromd =1;
    char MyVar[100]; char MyVar2[100]; char tmpA1[100]; char tmpA2[100];
    char tmpB1[100]; char tmpB2[100]; char tmpNo[100]; 

    double myvar = mxhadfit;
    double myvar2 = mm2;
    sprintf(MyVar,"%s",thevar);
    sprintf(MyVar2,"%s",thevar2);
    sprintf(tmpA1,"mxhadfit");    sprintf(tmpB1,"mxhadfit");
    sprintf(tmpA2,"mm2"); sprintf(tmpB2,"mm2"); sprintf(tmpNo,"NoComp");
    //    if (!strcmp(MyVar,tmpA1)) {myvar = mxhadfit;}
    if (!strcmp(MyVar,tmpA2)) {myvar = mm2;}
    if (!strcmp(MyVar2,tmpB1)) {myvar2 = mxhadfit;}
    //    if (!strcmp(MyVar2,tmpB2)) {myvar2 = mm2;}
    
    if (!((myvar>0) || (myvar<0) || (myvar == 0))) {
      myvar = -999.;
    }   
    if(!Dfromd)  continue;    
    double w;
    w = 1;    //calculate reweightings
    //    w =  getBsysweight(Gvxbtyp,vub);                        //Bdec weighting
    //    w *= getDsysweight(GfDpi,GfDk,GfDks,GfDpiz,GfDlep,dImode,vub); //Ddec weighting

    if(semilep) {
      dcat =  retDScat(Gvxbtyp);
      sprintf(le,"S");
      if(dcat <0) dcat = 9;
    } else {
      dcat =  retDcat(GfDk,GfDks,GfDpi,GfDpiz);
      sprintf(le,"B");
      if(dcat <0) {
	dcat = 33;
	if((GfD0Ds == 1) || (GfD0Ds == 2)) dcat = 31;
	if((GfDDs == 1) || (GfDDs == 2)) dcat = 32;
      }
    }
    
    //    cout<<dcat<<endl;
    if(PCMS && NLE && FLAV && BCH && IPUR && MXCUT) {
      cutflav = 5;
      if(TMath::Abs(brecocharge)) cutflav = 3;
      if(TMath::Abs(brecocharge)==0 && flav==0) cutflav = 4;	  
      if(strcmp(MyVar2,tmpNo)) {
	//	cout<<" New stuff"<<endl;
	//Mx vs Mnu
	id = hist(myvar);
	id2 = hist2(myvar2);
	//	cout<<id<<" "<<id2<<" "<<cutflav<<endl;
	if(MM2 && CH && NLE1){  
	  group = 0;
	  sprintf(name, "%s%d%s%d",le,cat * 100000 + cutflav*100,"ct",dcat);  
	  //	  ((TH2D*)gDirectory->Get(name))->Fill(id,id2);
	  ((TH2D*)gDirectory->Get(name))->Fill(myvar,myvar2);
	}
      } else {
	//	cout<<" Old stuff"<<endl;
	id = hist(myvar);
	idflav = 200 + id + cat * 1000+ cutflav * 10000;    
	group = 100000;
	sprintf(name, "%s%d%s%d",le,group+idflav,"ct",dcat);
	((TH1D*)gDirectory->Get(name))->Fill(mes,totweight);  
	if(MM2 && CH && NLE1){  
	  group = 0;
	  sprintf(name, "%s%d%s%d",le,group+idflav,"ct",dcat);
	  ((TH1D*)gDirectory->Get(name))->Fill(mes,totweight);  
	}
      }
    }
    if(PCMS && NLE && FLAV && BCH && IPUR && MXCUT){  
      //    if(PCMS && NLE && !flav && BCH && IPUR) {???
      //Mx vs Mnu
      group = 90206 + cat *1000 + 100000;
      sprintf(name, "%s%d%s%d",le,group,"ct",dcat);
      ((TH1D*)gDirectory->Get(name))->Fill(mes,totweight);  
      if(MM2 && CH && NLE1){  
	group = 90206 + cat *1000;
	sprintf(name, "%s%d%s%d",le,group,"ct",dcat);
	((TH1D*)gDirectory->Get(name))->Fill(mes,totweight);  
      }
    }
  }
}

void BkgbrekTool::LoopT(int nevents, int cat, int bsel)

{
  gROOT->cd();
  //id = Mx category == Mx bin
  int id, idflav, dcat, group;        
  char name[100];  char le[100];char named[100]; char namesl[100]; 
  // restore B0/B+ flag
  int isbch =  bsel>=0? bsel: -(1+bsel);
  if (fChain == 0) return;
  Int_t nentries = Int_t(fChain->GetEntries());
  if( nentries > nevents) nentries = nevents;
  Int_t nbytes = 0, nb = 0;
  TH2 *hsc = new TH2D("Other","npi,pi0",6,0.,6.,6,0.,6.);  hsc->Sumw2();
      
  for (Int_t jentry=0; jentry<nentries;jentry++) {
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    //    if(bsel<0)
    //      if( (!vub && cat==4) || (!vcb && cat==5) )continue;
    //Requesting depleted region
    if(nks + nkp > 0) continue;
    if(!vcb) continue;
    int flav =  lcharge + brecoflav; // charge correlation
    //    bool ksele = nkp + nks;          // fit on the depleted sample?
    int ch = xcharge + brecocharge;  // total charge
    int cutflav = 0;              // flavor category (3 = charged B, 4 = neutral B OS, 5 = neutral B SS)  
    
    // cuts
    bool PCMS =  pcms>CUTPCMS;
    //    bool NPI0 = nnpi0<2;
    bool NLE = nle > 0; bool NLE1 = nle == 1; bool IPUR = intpur>INTPURITY;
    bool MM2 =  mm2 < MM2CUT; bool CH = ch == 0;
    int FLAV = !(TMath::Abs(brecocharge)!=0 && (lcharge + brecoflav)!=0);
    bool BCH = TMath::Abs(brecocharge) == isbch;     
    char MYVar[100]; char V1[100]; char V2[100];
    bool MXCUT = mxhadfit < MXMAXCUT;
    sprintf(V1,"mm2");
    sprintf(V2,"reso");
    sprintf(MYVar,"%s",thevar);

    if (isbch == 2) BCH = 1;
    double myvar = mxhadfit;
    
    if (!strcmp(MYVar,V1)) myvar = mm2;     
    if (!strcmp(MYVar,V2)) myvar = mxhadfit-mxhadgen;     

    if (!((myvar>0) || (myvar<0) || (myvar == 0))) {
      myvar = -999.;
    }   

    sprintf(le,"B"); dcat =0;
    if((GfD0Ds == 1) || (GfD0Ds == 2)) { dcat =1;} 
    if((GfDDs == 1) || (GfDDs == 2)) { dcat =2;}
    //    Remove to restore old ana
    if(GfDkl > 0)  {
      dcat =5;
    } else if(GfDnu > 0) {
      dcat =6;
    } else if(GfDkspiopio > 0) {
      dcat =4;
    } else if(GfDkspipi > 0) {
      dcat =3;
    } else if(GfDkmiss > 0) {
      dcat =7;
    } else if(GfDk >= GfDkmiss+1)  {
      dcat =8;
      //    } else if(GfDlep > 0)  {
      //      dcat =6;
    } else {
      if(GfDpi >0 || GfDpiz >0) {
	dcat = 9;
      } 
      /*
	else {
	cout<<"Weird Event "<<GfDkspiopio<<" "<<GfDkspipi<<" "<<GfDkmiss<<" "<<GfDk<<" "<<GfDlep<<" "<<GfDpi<<" "<<GfDpiz<<endl;
	}
      */
    }
    int scat =10;
    int avxbtyp = TMath::Abs(Gvxbtyp);

    if(avxbtyp >=0 && avxbtyp < 7) {
      if(avxbtyp ==2 || avxbtyp ==1) {
	scat = scat + avxbtyp;
      } else { 
	scat = 13;
      }
    }    
    /* else {
       cout<<avxbtyp<<" :: Vxb Category"<<endl;
       }*/
    //    cout<<scat<<" :: Semileptonic Category"<<endl;
    if(PCMS && NLE && FLAV && BCH && IPUR && MXCUT) {
      cutflav = 5;
      if(TMath::Abs(brecocharge)) cutflav = 3;
      if(TMath::Abs(brecocharge)==0 && flav==0) cutflav = 4;	  
      id = hist(myvar);
      idflav = 200 + id + cat * 1000+ cutflav * 10000;    
      group = 100000;
      //Comment
      sprintf(name, "%s%d%s%d",le,group+idflav,"ct",dcat);
      ((TH1D*)gDirectory->Get(name))->Fill(mes,totweight);  
      //Semileptonic Breakdown
      sprintf(namesl, "%s%d%s%d",le,group+idflav,"ct",scat);
      ((TH1D*)gDirectory->Get(namesl))->Fill(mes,totweight);  
      //Uncomment
      sprintf(named, "%s%d%s%d",le,group+idflav,"ct",0);
      ((TH1D*)gDirectory->Get(named))->Fill(mes,totweight);  

      if(MM2 && CH && NLE1){  
	group = 0;
 	sprintf(named, "%s%d%s%d",le,group+idflav,"ct",0);
	((TH1D*)gDirectory->Get(named))->Fill(mes,totweight);  
	//COMMENT
	sprintf(name, "%s%d%s%d",le,group+idflav,"ct",dcat);
	((TH1D*)gDirectory->Get(name))->Fill(mes,totweight);  
	//Semileptonic Breakdown
	sprintf(namesl, "%s%d%s%d",le,group+idflav,"ct",scat);
	((TH1D*)gDirectory->Get(namesl))->Fill(mes,totweight);  
	//Uncomment

	if(dcat == 9) {
	  hsc->Fill(GfDpi,GfDpiz,totweight);
	}
      }
    }
    if(PCMS && NLE && FLAV && BCH && IPUR && MXCUT){  
      group = 90206 + cat *1000 + 100000;
      sprintf(named, "%s%d%s%d",le,group,"ct",0);
      ((TH1D*)gDirectory->Get(named))->Fill(mes,totweight);  
      //Comment
      sprintf(name, "%s%d%s%d",le,group,"ct",dcat);
      ((TH1D*)gDirectory->Get(name))->Fill(mes,totweight);  
      //Semileptonic Breakdown
      sprintf(namesl, "%s%d%s%d",le,group,"ct",scat);
      ((TH1D*)gDirectory->Get(namesl))->Fill(mes,totweight);  
      //Uncomment
      
      if(MM2 && CH && NLE1){  
	group = 90206 + cat *1000;
	//Comment
	if(dcat<2 || dcat >9) {
	  cout<<"* pippopippopippopippopippopippo *"<<endl;
	  cout<<"* pippo!* "<<dcat<<" pippo!* *"<<endl;
	  cout<<"* pippopippopippopippopippopippo *"<<endl;
	} 
 	sprintf(name, "%s%d%s%d",le,group,"ct",dcat);
 	((TH1D*)gDirectory->Get(name))->Fill(mes,totweight);  
	//Semileptonic Breakdown
	sprintf(namesl, "%s%d%s%d",le,group,"ct",scat);
	((TH1D*)gDirectory->Get(namesl))->Fill(mes,totweight);  
	//Comment finish
	sprintf(named, "%s%d%s%d",le,group,"ct",0);
 	((TH1D*)gDirectory->Get(named))->Fill(mes,totweight);  
      }
    }
  }
  c6 = new TCanvas("c6", "c6",600,800);   
  c6->Clear(); 
  hsc->SetMarkerSize(5);
  hsc->Draw();
  char names [100];
  sprintf(names,"%s%d%s","sc_",cat,"plot.eps");
  c6->Print(names);
}

void BkgbrekTool::Bookhist(int mxmns)
{
  gROOT->cd();
  TH1 *h;  char name[100], title[200], number[100]; TH2 *h2;

  if(mxmns) {
    //First attempt to have just D decays
    //Now I have all the mes plots!
    for (int lcat=0;lcat<34;lcat++) {
      for(int iflav=3; iflav<7; iflav++) {
	sprintf(name,"%s%d%s%d", "B",300000+iflav*100,"ct",lcat);  sprintf(title, "%s%s", thevar, " data events after all cuts");  h2 = new TH2D(name, title, thebins, themin, themax, thebins2, themin2, themax2 );  h2->Sumw2();
	sprintf(name,"%s%d%s%d", "B",1300000+iflav*100,"ct",lcat);  sprintf(title, "%s%s", thevar, " MC events  after all cuts");  h2 = new TH2D(name, title, thebins, themin, themax, thebins2, themin2, themax2  );  h2->Sumw2();
	sprintf(name,"%s%d%s%d", "B",400000+iflav*100,"ct",lcat);  sprintf(title, "%s%s", thevar, " data events after all cuts");  h2 = new TH2D(name, title, thebins, themin, themax, thebins2, themin2, themax2 );  h2->Sumw2();
	sprintf(name,"%s%d%s%d", "B",500000+iflav*100,"ct",lcat);  sprintf(title, "%s%s", thevar, " MC events  after all cuts");  h2 = new TH2D(name, title, thebins, themin, themax, thebins2, themin2, themax2  );  h2->Sumw2();
	sprintf(name,"%s%d%s%d", "B",1400000+iflav*100,"ct",lcat);  sprintf(title, "%s%s", thevar, " data events after lepton cuts");  h2 = new TH2D(name, title, thebins, themin, themax, thebins2, themin2, themax2 );  h2->Sumw2();
	sprintf(name,"%s%d%s%d", "B",1500000+iflav*100,"ct",lcat);  sprintf(title, "%s%s", thevar, " MC events after lepton cuts");  h2 = new TH2D(name, title, thebins, themin, themax, thebins2, themin2, themax2  );  h2->Sumw2();
      }
      //Plots for Mes for all bins divided per D decay
      sprintf(name,"%s%d","B93206ct",lcat);  sprintf(title, "mes data after all cuts: depleted");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
      sprintf(name,"%s%d","B94206ct",lcat);  sprintf(title, "mes data after all cuts: depleted");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
      sprintf(name,"%s%d","B95206ct",lcat);  sprintf(title, "mes data after all cuts: depleted");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
      sprintf(name,"%s%d","B194206ct",lcat);  sprintf(title, "mes data after all cuts: depleted");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
      sprintf(name,"%s%d","B195206ct",lcat);  sprintf(title, "mes data after all cuts: depleted");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
      sprintf(name,"%s%d","B193206ct",lcat);  sprintf(title, "mes data after all cuts: depleted");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    }
  } else {
    //First attempt to have just D decays
    //Now I have all the mes plots!
    for (int lcat=0;lcat<34;lcat++) {
      for(int iflav=3; iflav<7; iflav++) {
	sprintf(name,"%s%d%s%d", "B",300000+iflav*100,"ct",lcat);  sprintf(title, "%s%s", thevar, " data events after all cuts");  h = new TH1D(name, title, thebins, themin, themax );  h->Sumw2();
	sprintf(name,"%s%d%s%d", "B",1300000+iflav*100,"ct",lcat);  sprintf(title, "%s%s", thevar, " MC events  after all cuts");  h = new TH1D(name, title, thebins, themin, themax  );  h->Sumw2();
	sprintf(name,"%s%d%s%d", "B",400000+iflav*100,"ct",lcat);  sprintf(title, "%s%s", thevar, " data events after all cuts");  h = new TH1D(name, title, thebins, themin, themax );  h->Sumw2();
	sprintf(name,"%s%d%s%d", "B",500000+iflav*100,"ct",lcat);  sprintf(title, "%s%s", thevar, " MC events  after all cuts");  h = new TH1D(name, title, thebins, themin, themax  );  h->Sumw2();
	sprintf(name,"%s%d%s%d", "B",1400000+iflav*100,"ct",lcat);  sprintf(title, "%s%s", thevar, " data events after lepton cuts");  h = new TH1D(name, title, thebins, themin, themax );  h->Sumw2();
	sprintf(name,"%s%d%s%d", "B",1500000+iflav*100,"ct",lcat);  sprintf(title, "%s%s", thevar, " MC events after lepton cuts");  h = new TH1D(name, title, thebins, themin, themax  );  h->Sumw2();
      }
      //Plots for Mes for all bins divided per D decay
      sprintf(name,"%s%d","B93206ct",lcat);  sprintf(title, "mes data after all cuts: depleted");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
      sprintf(name,"%s%d","B94206ct",lcat);  sprintf(title, "mes data after all cuts: depleted");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
      sprintf(name,"%s%d","B95206ct",lcat);  sprintf(title, "mes data after all cuts: depleted");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
      sprintf(name,"%s%d","B194206ct",lcat);  sprintf(title, "mes data after all cuts: depleted");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
      sprintf(name,"%s%d","B195206ct",lcat);  sprintf(title, "mes data after all cuts: depleted");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
      sprintf(name,"%s%d","B193206ct",lcat);  sprintf(title, "mes data after all cuts: depleted");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
      //Plots for Mes for each bin divided per D decay
      for (int lo=1;lo<thebins+1;lo++) {
	
	sprintf(number, "%d",  200+lo);   //Bin flag
	sprintf(name,"%s%s%s%d" , "B4",number,"ct",lcat);  sprintf(title, "mes data");  h = new TH1D(name, title, 40, 5.2, 5.3);   h->Sumw2();
	sprintf(name,"%s%s%s%d" , "B44",number,"ct",lcat);  sprintf(title, "mes data");  h = new TH1D(name, title, 40, 5.2, 5.3);   h->Sumw2();
	sprintf(name,"%s%s%s%d" , "B54",number,"ct",lcat);  sprintf(title, "mes data");  h = new TH1D(name, title, 40, 5.2, 5.3);   h->Sumw2();
	sprintf(name,"%s%s%s%d" , "B34",number,"ct",lcat);  sprintf(title, "mes data");  h = new TH1D(name, title, 40, 5.2, 5.3);   h->Sumw2();
	sprintf(name,"%s%s%s%d" , "B14",number,"ct",lcat);  sprintf(title, "mes data");  h = new TH1D(name, title, 40, 5.2, 5.3);   h->Sumw2();
	sprintf(name,"%s%s%s%d" , "B154",number,"ct",lcat);  sprintf(title, "mes data");  h = new TH1D(name, title, 40, 5.2, 5.3);   h->Sumw2();
	sprintf(name,"%s%s%s%d" , "B144",number,"ct",lcat);  sprintf(title, "mes data");  h = new TH1D(name, title, 40, 5.2, 5.3);   h->Sumw2();
	sprintf(name,"%s%s%s%d" , "B134",number,"ct",lcat);  sprintf(title, "mes data");  h = new TH1D(name, title, 40, 5.2, 5.3);   h->Sumw2();
	sprintf(name,"%s%s%s%d" , "B5",number,"ct",lcat);  sprintf(title, "mes data");  h = new TH1D(name, title, 40, 5.2, 5.3);   h->Sumw2();
	sprintf(name,"%s%s%s%d" , "B45",number,"ct",lcat);  sprintf(title, "mes data");  h = new TH1D(name, title, 40, 5.2, 5.3);   h->Sumw2();
	sprintf(name,"%s%s%s%d" , "B55",number,"ct",lcat);  sprintf(title, "mes data");  h = new TH1D(name, title, 40, 5.2, 5.3);   h->Sumw2();
	sprintf(name,"%s%s%s%d" , "B35",number,"ct",lcat);  sprintf(title, "mes data");  h = new TH1D(name, title, 40, 5.2, 5.3);   h->Sumw2();
	sprintf(name,"%s%s%s%d" , "B15",number,"ct",lcat);  sprintf(title, "mes data");  h = new TH1D(name, title, 40, 5.2, 5.3);   h->Sumw2();
	sprintf(name,"%s%s%s%d" , "B155",number,"ct",lcat);  sprintf(title, "mes data");  h = new TH1D(name, title, 40, 5.2, 5.3);   h->Sumw2();
	sprintf(name,"%s%s%s%d" , "B145",number,"ct",lcat);  sprintf(title, "mes data");  h = new TH1D(name, title, 40, 5.2, 5.3);   h->Sumw2();
	sprintf(name,"%s%s%s%d" , "B135",number,"ct",lcat);  sprintf(title, "mes data");  h = new TH1D(name, title, 40, 5.2, 5.3);   h->Sumw2();
	//      cout<< name <<" "<<endl;
      }
    }
    
    for (int lcat=0;lcat<10;lcat++) {
      for (int iflav=3;iflav<7;iflav++) {
	sprintf(name,"%s%d%s%d", "S",300000+iflav*100,"ct",lcat);  sprintf(title, "%s%s", thevar, " data events after all cuts");  h = new TH1D(name, title, thebins, themin, themax );  h->Sumw2();
	sprintf(name,"%s%d%s%d", "S",1300000+iflav*100,"ct",lcat);  sprintf(title, "%s%s", thevar, " MC events  after all cuts");  h = new TH1D(name, title, thebins, themin, themax  );  h->Sumw2();
	sprintf(name,"%s%d%s%d", "S",400000+iflav*100,"ct",lcat);  sprintf(title, "%s%s", thevar, " data events after all cuts");  h = new TH1D(name, title, thebins, themin, themax );  h->Sumw2();
	sprintf(name,"%s%d%s%d", "S",500000+iflav*100,"ct",lcat);  sprintf(title, "%s%s", thevar, " MC events  after all cuts");  h = new TH1D(name, title, thebins, themin, themax  );  h->Sumw2();
	sprintf(name,"%s%d%s%d", "S",1400000+iflav*100,"ct",lcat);  sprintf(title, "%s%s", thevar, " data events after lepton cuts");  h = new TH1D(name, title, thebins, themin, themax );  h->Sumw2();
	sprintf(name,"%s%d%s%d", "S",1500000+iflav*100,"ct",lcat);  sprintf(title, "%s%s", thevar, " MC events after lepton cuts");  h = new TH1D(name, title, thebins, themin, themax  );  h->Sumw2();
      }
      //Plots for Mes for all bins divided per D decay
      sprintf(name,"%s%d","S93206ct",lcat);  sprintf(title, "mes data after all cuts: depleted");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
      sprintf(name,"%s%d","S94206ct",lcat);  sprintf(title, "mes data after all cuts: depleted");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
      sprintf(name,"%s%d","S95206ct",lcat);  sprintf(title, "mes data after all cuts: depleted");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
      sprintf(name,"%s%d","S194206ct",lcat);  sprintf(title, "mes data after all cuts: depleted");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
      sprintf(name,"%s%d","S195206ct",lcat);  sprintf(title, "mes data after all cuts: depleted");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
      sprintf(name,"%s%d","S193206ct",lcat);  sprintf(title, "mes data after all cuts: depleted");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
      //Plots for Mes for each bin divided per D decay
      for (int lo=1;lo<thebins+1;lo++) {
	
	sprintf(number, "%d",  200+lo);   //Bin flag
	sprintf(name,"%s%s%s%d" , "S4",number,"ct",lcat);  sprintf(title, "mes data");  h = new TH1D(name, title, 40, 5.2, 5.3);   h->Sumw2();
	sprintf(name,"%s%s%s%d" , "S44",number,"ct",lcat);  sprintf(title, "mes data");  h = new TH1D(name, title, 40, 5.2, 5.3);   h->Sumw2();
	sprintf(name,"%s%s%s%d" , "S54",number,"ct",lcat);  sprintf(title, "mes data");  h = new TH1D(name, title, 40, 5.2, 5.3);   h->Sumw2();
	sprintf(name,"%s%s%s%d" , "S34",number,"ct",lcat);  sprintf(title, "mes data");  h = new TH1D(name, title, 40, 5.2, 5.3);   h->Sumw2();
	sprintf(name,"%s%s%s%d" , "S14",number,"ct",lcat);  sprintf(title, "mes data");  h = new TH1D(name, title, 40, 5.2, 5.3);   h->Sumw2();
	sprintf(name,"%s%s%s%d" , "S154",number,"ct",lcat);  sprintf(title, "mes data");  h = new TH1D(name, title, 40, 5.2, 5.3);   h->Sumw2();
	sprintf(name,"%s%s%s%d" , "S144",number,"ct",lcat);  sprintf(title, "mes data");  h = new TH1D(name, title, 40, 5.2, 5.3);   h->Sumw2();
	sprintf(name,"%s%s%s%d" , "S134",number,"ct",lcat);  sprintf(title, "mes data");  h = new TH1D(name, title, 40, 5.2, 5.3);   h->Sumw2();
	//      cout<< name <<" "<<endl;
      }
    }
  }
}

void BkgbrekTool::Fitmes(int cat, int cut, TString flag, int semilep){

  // fit to mes distribution and fill of the Mx plots...
  gROOT->cd();
  char name[100], namedi[100], le[100];
  double sigs[1000], errsigs[1000], tempbin, tempbinchb, tempbinb0os;
  double tempbinb0ss, temperr, temperrchb, temperrb0os;
  double temperrb0ss, chid = 0.174; 
  double foomean, foosigma, fooalpha, foon, usedmean, usedsigma;
  double usedalpha, usedn;  int i; 
  double usedmeand, usedsigmad;
  double usedalphad, usednd; 
  char nameou[100], namebch[100], nameb0os[100], nameb0ss[100];
  int addcut = 0; char tmpNo[100];
  sprintf(tmpNo,"NoComp");
  if(cut) addcut = 1;
  int const theb = thebins + 1;
  int const theb2 = thebins2 + 1;
  int NCat =34;
  sprintf(le, "B");
  if(cut == 0) {
    sprintf(nameou,"%s%s%s","Entries_Nocut_",flag.Data(),".txt");
  } else if (cut == 1) {
    sprintf(nameou,"%s%s%s","Entries_Cut_",flag.Data(),".txt");
  }
  
  if(semilep ==1) {
    NCat =10;     sprintf(le, "S");
    if(cut == 0) {
      sprintf(nameou,"%s%s%s","Entries_SNocut_",flag.Data(),".txt");
    } else if (cut == 1) {
      sprintf(nameou,"%s%s%s","Entries_SCut_",flag.Data(),".txt");
    }
  } else if (semilep ==2) {
    sprintf(nameou,"%s%s%s","Entries_KsStu_",flag.Data(),".txt");
    NCat =14; 
    //    NCat =10; 
  }
  
  ofstream outFile(nameou,ios::app);
  char MyVar2[100]; int il; char psname[200];
  sprintf(MyVar2,"%s",thevar2);
  int group = 200 + cat * 1000;

  //Fix parameters to Overall distributions
  sprintf(name, "%s%d%s","B",group+90006+addcut*100000,"ct0");  
  double dummyd1,dummyd2;
  sighisto(dummyd1,dummyd2,(TH1D*)gDirectory->Get(name),usedmeand,usedsigmad,usedalphad,usednd,1,-11111111.,-1111111.,-1111111.,5.,-1111111.);

  for (int y=0; y<NCat; y++) {

    if (semilep ==2 && (y ==1 || y ==2 || y==10)) {
      cout<<"Skip dummy cat"<<endl;
    } else {
    //extracting the signal fit parameters 
    sprintf(name, "%s%d%s%d",le,group+90006+addcut*100000,"ct",y);  
    cout << ((TH1D*)gDirectory->Get(name))->Integral() << endl;
    double dummy1,dummy2;
    //    sighisto(dummy1,dummy2,(TH1D*)gDirectory->Get(name),usedmean,usedsigma,usedalpha,usedn,1,-11111111.,-1111111.,-1111111.,5.,-1111111.);
    sighisto(dummy1,dummy2,(TH1D*)gDirectory->Get(name),usedmean,usedsigma,usedalpha,usedn,1,usedmeand,usedsigmad,usedalphad,usednd,-1111111);

    cout << "mes result for: " << name << " is, MEAN " << usedmean << " SIGMA " << usedsigma << " ALPHA " << usedalpha << " N " << usedn << endl;
    outFile<<dummy1<<" "<<dummy2<<endl;    
    /*
      c6 = new TCanvas("c6", "c6",1200,900);   
      c6->Clear(); c6->cd();
      ((TH1D*)gDirectory->Get(name))->Draw();
      sprintf(psname,"%s%s",name,".eps");
      c6->Print(psname);
    */
    if (semilep ==2) {
      fNall[y][cat-4] = dummy1;
    }
    for(int iflav =3; iflav<6; iflav++) {
      int title = cat * 100000 + addcut * 1000000 + iflav*100;
      sprintf(name, "%s%d%s%d",le,title,"ct",y);  
      for (i = 1;  i < theb; i++) { 
	sprintf(namedi ,  "%s%d%s%d",le,group+i+iflav*10000+addcut*100000,"ct",y);  
	//	sighisto(sigs[i-1],errsigs[i-1],(TH1D*)gDirectory->Get(namedi),foomean,foosigma,fooalpha,foon,1,-11111111.,-1111111.,-1111111.,-1111111.,-1111111.);
	sighisto(sigs[i-1],errsigs[i-1],(TH1D*)gDirectory->Get(namedi),foomean,foosigma,fooalpha,foon,1,usedmean,usedsigma,usedalpha,usedn,-1111111);
	((TH1D*)gDirectory->Get(name))->SetBinContent(i, sigs[i-1]);
	((TH1D*)gDirectory->Get(name))->SetBinError(i, errsigs[i-1]);
      }
    }
    
    //Old stuff
    sprintf(name,     "%s%d%s%d",le,cat*100000+600+addcut*1000000,"ct",y);  
    sprintf(namebch , "%s%d%s%d",le,cat*100000+300+addcut*1000000,"ct",y);  
    sprintf(nameb0os, "%s%d%s%d",le,cat*100000+400+addcut*1000000,"ct",y);  
    sprintf(nameb0ss, "%s%d%s%d",le,cat*100000+500+addcut*1000000,"ct",y);  
    for(int k=1;k<theb;k++){
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
  }
  }
}

void BkgbrekTool::overlap(int cut, TString dir, TString flag){
  
  char name[200];  int addcut = 0; 
  if(cut==1) addcut = 1;
  tl.SetNDC(kTRUE); char line[200];
  
  gROOT->cd();
  //Old stuff
  if (cut ==2) {
    cout<<"The problem"<<endl;
    int const theb = thebins + 1;
    int const theb2 = thebins2 + 1;
    double tempbin, tempbinchb, tempbinb0os;
    double tempbinb0ss, temperr, temperrchb, temperrb0os;
    double temperrb0ss, chid = 0.174; 
    char namebch[100], nameb0os[100], nameb0ss[100];
    char name[100], le[100];
    sprintf(le,"B");
    int NCat = 34;
    for (int y=0; y<NCat; y++) {
      sprintf(namebch , "%s%d%s%d",le,400000+300+addcut*1000000,"ct",y);  
      sprintf(nameb0os, "%s%d%s%d",le,400000+400+addcut*1000000,"ct",y);  
      sprintf(nameb0ss, "%s%d%s%d",le,400000+500+addcut*1000000,"ct",y);  
      sprintf(name,     "%s%d%s%d",le,400000+600+addcut*1000000,"ct",y);  
      for(int k=1;k<theb;k++){
	for(int kk=1;kk<theb2;kk++){
	  tempbinchb = ((TH2D*)gDirectory->Get(namebch))->GetBinContent(k,kk);
	  tempbinb0os = ((TH2D*)gDirectory->Get(nameb0os))->GetBinContent(k,kk);
	  tempbinb0ss = ((TH2D*)gDirectory->Get(nameb0ss))->GetBinContent(k,kk);
	  //	  temperrchb = ((TH2D*)gDirectory->Get(namebch))->GetBinError(k,kk);
	  //	  temperrb0os = ((TH2D*)gDirectory->Get(nameb0os))->GetBinError(k,kk);
	  //	  temperrb0ss = ((TH2D*)gDirectory->Get(nameb0ss))->GetBinError(k,kk);
	  tempbin = tempbinchb + ((1-chid)/(1-2*chid)) * tempbinb0os - (chid/(1-2*chid)) * tempbinb0ss;
	  //	  temperr = sqrt(temperrchb*temperrchb + ((1-chid)/(1-2*chid)) * ((1-chid)/(1-2*chid)) * tempbinb0os + (chid/(1-2*chid)) * (chid/(1-2*chid))* tempbinb0ss);
	  ((TH2D*)gDirectory->Get(name))->SetBinContent(k,kk, tempbin);
	}
      }
    }
  }
  c6 = new TCanvas("c6", "c6",1200,900);   
  c6->Clear(); c6->Divide(3,3); 
  for(int iDd=0; iDd<8; iDd++) {
    strstream Lpi, Lpiz, LK, LKs;
    sprintf(name,"%s%d%s%d" ,"B", 400600 + addcut * 1000000,"ct",iDd);
    c6->cd(iDd+1);
    if(MA[iDd][0] != 0) Lpi  << MA[iDd][0] << "#pi" <<ends;
    if(MA[iDd][1] != 0) Lpiz  << MA[iDd][1] << "#pi^{0}" <<ends;
    if(MA[iDd][2] != 0) LK  << MA[iDd][2] << "K" <<ends;
    if(MA[iDd][3] != 0) LKs  << MA[iDd][3] << "K^{0}" <<ends;
    sprintf(line,"%s%s%s%s",LK.str(),LKs.str(),Lpi.str(),Lpiz.str());
    tl.SetTextSizePixels(50.1);
    if(cut<2)    ((TH1D*)gDirectory->Get(name))->Draw();   
    if(cut == 2)  ((TH2D*)gDirectory->Get(name))->Draw("CONT1");   
    tl.DrawLatex(0.22, 0.75, line);    gPad->Update();
  }
  sprintf(name, "%s%s%s%s", dir.Data(),"/",flag.Data(),"_fi.ps");
  c6->Print(name);  
  
  c7 = new TCanvas("c7", "c7",1200,900);    
  c7->Clear(); c7->Divide(3,3);
  for(int iDd=8; iDd<16; iDd++) {
    strstream Lpi, Lpiz, LK, LKs;
    //    Lpi   <<""<<ends;    Lpiz  <<""<<ends;    LK    <<""<<ends;
    //    LKs   <<""<<ends;
    sprintf(name,"%s%d%s%d" ,"B", 400600 + addcut * 1000000,"ct",iDd);
    c7->cd(iDd-7);
    if(MA[iDd][0] != 0) Lpi  << MA[iDd][0] << "#pi" <<ends;
    if(MA[iDd][1] != 0) Lpiz  << MA[iDd][1] << "#pi^{0}" <<ends;
    if(MA[iDd][2] != 0) LK  << MA[iDd][2] << "K" <<ends;
    if(MA[iDd][3] != 0) LKs  << MA[iDd][3] << "K^{0}" <<ends;
    sprintf(line,"%s%s%s%s",LK.str(),LKs.str(),Lpi.str(),Lpiz.str());
    tl.SetTextSizePixels(50.1);
    if(cut<2)    ((TH1D*)gDirectory->Get(name))->Draw(); 
    if(cut==2)    ((TH2D*)gDirectory->Get(name))->Draw("CONT1"); 
    tl.DrawLatex(0.22, 0.75, line);    gPad->Update();
  }
  sprintf(name, "%s%s%s%s", dir.Data(),"/",flag.Data(),"_sc.ps");
  c7->Print(name);  
  c2 = new TCanvas("c2", "c2",1200,900);   
  c2->Clear(); c2->Divide(3,3); 
  for(int iDd=16; iDd<24; iDd++) {
    strstream Lpi, Lpiz, LK, LKs;
    //    Lpi   <<""<<ends;    Lpiz  <<""<<ends;    LK    <<""<<ends;
    //    LKs   <<""<<ends;
    sprintf(name,"%s%d%s%d" ,"B", 400600 + addcut * 1000000,"ct",iDd);
    c2->cd(iDd-15);
    if(MA[iDd][0] != 0) Lpi  << MA[iDd][0] << "#pi" <<ends;
    if(MA[iDd][1] != 0) Lpiz  << MA[iDd][1] << "#pi^{0}" <<ends;
    if(MA[iDd][2] != 0) LK  << MA[iDd][2] << "K" <<ends;
    if(MA[iDd][3] != 0) LKs  << MA[iDd][3] << "K^{0}" <<ends;
    sprintf(line,"%s%s%s%s",LK.str(),LKs.str(),Lpi.str(),Lpiz.str());
    tl.SetTextSizePixels(50.1);
    if(cut<2)    ((TH1D*)gDirectory->Get(name))->Draw(); 
    if(cut==2)    ((TH2D*)gDirectory->Get(name))->Draw("CONT1"); 
    tl.DrawLatex(0.22, 0.75, line);    gPad->Update();
  }
  sprintf(name, "%s%s%s%s", dir.Data(),"/",flag.Data(),"_th.ps");
  c2->Print(name);  
  
  c3 = new TCanvas("c3", "c3",1200,900);    
  c3->Clear(); c3->Divide(3,3);
  for(int iDd=24; iDd<31; iDd++) {
    strstream Lpi, Lpiz, LK, LKs;
    //    Lpi   <<""<<ends;    Lpiz  <<""<<ends;    LK    <<""<<ends;
    //    LKs   <<""<<ends;
    sprintf(name,"%s%d%s%d" ,"B", 400600 + addcut * 1000000,"ct",iDd);
    c3->cd(iDd-23);
    if(MA[iDd][0] != 0) Lpi  << MA[iDd][0] << "#pi" <<ends;
    if(MA[iDd][1] != 0) Lpiz  << MA[iDd][1] << "#pi^{0}" <<ends;
    if(MA[iDd][2] != 0) LK  << MA[iDd][2] << "K" <<ends;
    if(MA[iDd][3] != 0) LKs  << MA[iDd][3] << "K^{0}" <<ends;
    sprintf(line,"%s%s%s%s",LK.str(),LKs.str(),Lpi.str(),Lpiz.str());
    tl.SetTextSizePixels(50.1);
    if(cut<2)    ((TH1D*)gDirectory->Get(name))->Draw(); 
    if(cut==2)    ((TH2D*)gDirectory->Get(name))->Draw("CONT1"); 
    tl.DrawLatex(0.22, 0.75, line);    gPad->Update();
  }
  sprintf(name, "%s%s%s%s", dir.Data(),"/",flag.Data(),"_fo.ps");
  c3->Print(name);  
  
  c3 = new TCanvas("c3", "c3",1200,900);    
  c3->Clear(); c3->Divide(3,3);
  for(int iDd=31; iDd<34; iDd++) {
    strstream Lpi, Lpiz, LK, LKs;
    //    Lpi   <<""<<ends;    Lpiz  <<""<<ends;    LK    <<""<<ends;
    //    LKs   <<""<<ends;
    sprintf(name,"%s%d%s%d" ,"B", 400600 + addcut * 1000000,"ct",iDd);
    c3->cd(iDd-30);
    if(iDd == 31) Lpi  << "D0 events" <<ends;
    if(iDd == 32) Lpiz  << "D+ events" <<ends;
    if(iDd == 33) LK  << "Other events" <<ends;
    sprintf(line,"%s%s%s",LK.str(),Lpi.str(),Lpiz.str());
    tl.SetTextSizePixels(50.1);
    if(cut<2)    ((TH1D*)gDirectory->Get(name))->Draw(); 
    if(cut==2)    ((TH2D*)gDirectory->Get(name))->Draw("CONT1"); 
    tl.DrawLatex(0.22, 0.75, line);    gPad->Update();
  }
  sprintf(name, "%s%s%s%s", dir.Data(),"/",flag.Data(),"_ft.ps");
  c3->Print(name);  
}

void BkgbrekTool::overlapS(int cut, TString dir, TString flag){
  
  char name[200];  int addcut = 0; 
  if(cut) addcut = 1;
  tl.SetNDC(kTRUE);
  //  char Sltg[200]; 

  c6 = new TCanvas("c6", "c6",1200,900);   
  c6->Clear(); c6->Divide(3,3); 
  for(int iDd=0; iDd<9; iDd++) {
    sprintf(name,"%s%d%s%d" ,"S", 400600 + addcut * 1000000,"ct",iDd+1);
    c6->cd(iDd+1);
//    sprintf(Sltg,"\0");
    strstream Sltg;
    if(iDd == 0) Sltg << "D l #nu\0" <<ends;
    if(iDd == 1) Sltg << "D^{*} l #nu\0"<< ends;
    if(iDd == 3) Sltg <<"D_{2}^{*} l #nu\0"<< ends;
    if(iDd == 4) Sltg << "D_{1} l #nu\0"<< ends;
    if(iDd == 5) Sltg << "Other D (#gamma) decays\0"<< ends;
    if(iDd == 6) Sltg << "Xu resonant\0"<< ends;
    if(iDd == 7) Sltg << "Xu nonreso\0"<< ends;
    if(iDd == 8) Sltg << "Other not vcb\0"<< ends;
    tl.SetTextSizePixels(50.1);
    ((TH1D*)gDirectory->Get(name))->Draw();tl.DrawLatex(0.22, 0.75,Sltg.str());
    gPad->Update();
  }
  sprintf(name, "%s%s%s%s", dir.Data(),"/",flag.Data(),"_fi_Sem.ps");
  c6->Print(name);  
}

void BkgbrekTool::sameP(int cut, TString dir, TString flag, int Norm){
  
  int addcut = 0; char mm[100]; char myvar[100];
  if(cut) addcut = 1;
  tl.SetNDC(kTRUE);tl.SetTextSizePixels(50.1);

  double allb0=((TH1D*)gDirectory->Get("B400600ct0"))->Integral();
  double allbc=((TH1D*)gDirectory->Get("B500600ct0"))->Integral();
  double allb0d0=((TH1D*)gDirectory->Get("B400600ct1"))->Integral();
  double allbcd0=((TH1D*)gDirectory->Get("B500600ct1"))->Integral();
  double allb0dc=((TH1D*)gDirectory->Get("B400600ct2"))->Integral();
  double allbcdc=((TH1D*)gDirectory->Get("B500600ct2"))->Integral();
  
  double scal1,scal2,scal3;
  //Need to add another swith in case of 'cut'

  if(Norm == 0) {
    scal1 = 1;    scal2 = 1;    scal3 = 1;
  } else if(Norm == 1) { //Normalize to == number of eve
    scal1 = allbc/allb0; scal2 = allbc/allb0; scal3 = allb0/allbc;
  } else { //Normalize to area-events dc-bc
    scal1 = fNall[0][1]/fNall[0][0]; scal2 = fNall[1][1]/fNall[1][0];
    scal3 = fNall[2][0]/fNall[2][1];
  }
  sprintf(myvar,"%s",thevar);
  sprintf(mm,"mm2");
  int mis= strcmp(myvar,mm);
  
  c5 = new TCanvas("c5", "c5",600,1200);   
  c5->Clear(); c5->Divide(1,3);
  char Sltg[100], Sltg2[100], name[100], name2[100], name3[100];
  for(int adu = 0; adu<3; adu++) {
    c5->cd(adu+1);
    if (adu ==1) {
      sprintf(Sltg,"B+ (D0)");  sprintf(Sltg2,"B0 (D0)");  
      sprintf(name,"%s%d%s%d","B",500600 + addcut * 1000000,"ct",adu);
      sprintf(name2,"%s%d%s%d","B",400600 + addcut * 1000000,"ct",adu);
      setHist(((TH1D*)gDirectory->Get(name)),1,20,1,1,1);
      ((TH1D*)gDirectory->Get(name2))->Scale(scal2);
      setHist(((TH1D*)gDirectory->Get(name2)),2,24,1,1,1);
      ((TH1D*)gDirectory->Get(name))->Draw();
      ((TH1D*)gDirectory->Get(name2))->Draw("same");gPad->Update();
      //      if(!mis) {
      //	TPaveStats *st = (TPaveStats*)gPad->GetPrimitive("stats");
      //	st->SetX1NDC(0.1);  st->SetX2NDC(0.3);
      //	((TH1D*)gDirectory->Get(name2))->Draw("sames");  
      //      }
      tl.SetTextColor(1);
      tl.DrawLatex(0.22, 0.75,Sltg);
      tl.SetTextColor(2);
      tl.DrawLatex(0.22, 0.85,Sltg2);  gPad->Update();
    } else if (adu == 2) {
      sprintf(Sltg,"B+ (D+)");  sprintf(Sltg2,"B0 (D+)");  
      sprintf(name,"%s%d%s%d","B",400600 + addcut * 1000000,"ct",adu);
      sprintf(name2,"%s%d%s%d","B",500600 + addcut * 1000000,"ct",adu);
      setHist(((TH1D*)gDirectory->Get(name)),2,20,1,1,1);
      setHist(((TH1D*)gDirectory->Get(name2)),1,24,1,1,scal3);
      if(mis) {
	((TH1D*)gDirectory->Get(name))->Draw();
	((TH1D*)gDirectory->Get(name2))->Draw("same");
      } else {
	((TH1D*)gDirectory->Get(name2))->Draw();
	((TH1D*)gDirectory->Get(name))->Draw("same");gPad->Update();
	//	TPaveStats *st = (TPaveStats*)gPad->GetPrimitive("stats");
	//	st->SetX1NDC(0.1);  st->SetX2NDC(0.3);
	//	((TH1D*)gDirectory->Get(name2))->Draw("sames");  
      }
      tl.SetTextColor(2);
      tl.DrawLatex(0.22, 0.75,Sltg2);
      tl.SetTextColor(1);
      tl.DrawLatex(0.22, 0.85,Sltg);  gPad->Update();
    } else {
      sprintf(Sltg,"B+");  sprintf(Sltg2,"B0");  
      sprintf(name,"%s%d%s%d","B",500600 + addcut * 1000000,"ct",adu);
      sprintf(name2,"%s%d%s%d","B",400600 + addcut * 1000000,"ct",adu);
      setHist(((TH1D*)gDirectory->Get(name)),1,20,1,1,1);
      setHist(((TH1D*)gDirectory->Get(name2)),2,24,1,1,scal1);
      ((TH1D*)gDirectory->Get(name))->Draw();
      ((TH1D*)gDirectory->Get(name2))->Draw("same");gPad->Update();
      //      if(!mis) {
      //	TPaveStats *st = (TPaveStats*)gPad->GetPrimitive("stats");
      //	st->SetX1NDC(0.1);  st->SetX2NDC(0.3);
      //	((TH1D*)gDirectory->Get(name2))->Draw("sames");  
      //      }
      tl.SetTextColor(1);
      tl.DrawLatex(0.22, 0.75,Sltg);
      tl.SetTextColor(2);
      tl.DrawLatex(0.22, 0.85,Sltg2);  gPad->Update();
    }
  }
  sprintf(name3,"%s%s%s%s",dir.Data(),"/",flag.Data(),"_All_B.ps");
  c5->Print(name3);  

  if(Norm == 0) {
    scal1 = 1;
    scal2 = 1;
  } else if(Norm == 1) {
    scal1 = allbcd0/allbcdc;
    scal2 = allb0d0/allb0dc;
  } else {
    scal1 = fNall[1][1]/fNall[2][1];
    scal2 = fNall[1][0]/fNall[2][0];
  }
  c5->Clear(); c5->Divide(1,2);
  c5->cd(1); 
  sprintf(Sltg,"D0 (B+)");  sprintf(Sltg2,"D+ (B+)");  
  sprintf(name,"%s%d%s","B",500600 + addcut * 1000000,"ct1");
  sprintf(name2,"%s%d%s","B",500600 + addcut * 1000000,"ct2");
  setHist(((TH1D*)gDirectory->Get(name)),1,20,1,1,1);
  setHist(((TH1D*)gDirectory->Get(name2)),2,24,1,1,scal1);
  ((TH1D*)gDirectory->Get(name2))->Draw();
  ((TH1D*)gDirectory->Get(name))->Draw("same");gPad->Update();
  //  if(!mis) {
  //    TPaveStats *st = (TPaveStats*)gPad->GetPrimitive("stats");
  //    st->SetX1NDC(0.1);  st->SetX2NDC(0.3);
  //    ((TH1D*)gDirectory->Get(name))->Draw("sames");  
  //  }
  tl.SetTextColor(1);
  tl.DrawLatex(0.22, 0.75,Sltg);
  tl.SetTextColor(2);
  tl.DrawLatex(0.22, 0.85,Sltg2);  gPad->Update();

  c5->cd(2); 
  sprintf(Sltg,"D0 (B0)");  sprintf(Sltg2,"D+ (B0)");  
  sprintf(name,"%s%d%s","B",400600 + addcut * 1000000,"ct1");
  sprintf(name2,"%s%d%s","B",400600 + addcut * 1000000,"ct2");
  setHist(((TH1D*)gDirectory->Get(name)),1,20,1,1,1);
  setHist(((TH1D*)gDirectory->Get(name2)),2,24,1,1,scal2);
  ((TH1D*)gDirectory->Get(name))->Draw();
  ((TH1D*)gDirectory->Get(name2))->Draw("same");gPad->Update();
  //  if(!mis) {
  //    TPaveStats *st = (TPaveStats*)gPad->GetPrimitive("stats");
  //    st->SetX1NDC(0.1);  st->SetX2NDC(0.3);
  //    ((TH1D*)gDirectory->Get(name))->Draw("sames");  
  //  }
  tl.SetTextColor(1);  tl.DrawLatex(0.22, 0.75, Sltg);
  tl.SetTextColor(2);  tl.DrawLatex(0.22, 0.85, Sltg2);  gPad->Update();

  sprintf(name3,"%s%s%s%s",dir.Data(),"/",flag.Data(),"_All_D0Dc.ps");
  c5->Print(name3);  
}

void BkgbrekTool::sameK(int cut, TString dir, TString flag, int Norm){

  int addcut = 0; char mm[100]; char myvar[100];
  char name3[100];
  if(cut) addcut = 1;
  tl.SetNDC(kTRUE); tl.SetTextSizePixels(25.);
  tl.SetTextAngle(0);
  double allb0=((TH1D*)gDirectory->Get("B400600ct0"))->Integral();
  double allbc=((TH1D*)gDirectory->Get("B500600ct0"))->Integral();

  double allb0kpi=((TH1D*)gDirectory->Get("B400600ct3"))->Integral();
  double allbckpi=((TH1D*)gDirectory->Get("B500600ct3"))->Integral();
  double allb0kpio=((TH1D*)gDirectory->Get("B400600ct4"))->Integral();
  double allbckpio=((TH1D*)gDirectory->Get("B500600ct4"))->Integral();
  double allb0kl=((TH1D*)gDirectory->Get("B400600ct5"))->Integral();
  double allbckl=((TH1D*)gDirectory->Get("B500600ct5"))->Integral();
  double allb0k=((TH1D*)gDirectory->Get("B400600ct6"))->Integral();
  double allbck=((TH1D*)gDirectory->Get("B500600ct6"))->Integral();
  
  double scal1,scal2,scal3,scal4,scal5,scal6;
  //Need to add another swith in case of 'cut'

  if(Norm == 0) {
    scal1 = 1;    scal2 = 1;
    scal3 = 1;    scal4 = 1;
    scal5 = 1;    scal6 = 1;
  } else { //Normalize to area-events dc-bc
    scal1 = 1/allb0kpi;  scal2 = 1/allbckpi;
    scal3 = 1/allb0kpio; scal4 = 1/allbckpio;
    scal5 = 1/allb0kl;   scal6 = 1/allbckl;
  }
  sprintf(myvar,"%s",thevar);
  sprintf(mm,"mm2");
  int mis= strcmp(myvar,mm);

  c5 = new TCanvas("c5", "c5",600,1200);   
  c5->Clear(); c5->Divide(1,2);
  c5->cd(1); 
  char KSltg2[100], KSltg[100], KSltg3[100], Kname[100], Kname3[100], Kname4[100], Kname2[100];

  sprintf(KSltg,"B+; Ks->#pi#pi");      sprintf(KSltg2,"B+; Ks->#pi^{0}#pi^{0}");  
  sprintf(KSltg3,"B+; Kl");  
  sprintf(Kname,"%s%d%s" ,"B",500600 + addcut * 1000000,"ct3");
  sprintf(Kname2,"%s%d%s","B",500600 + addcut * 1000000,"ct4");
  sprintf(Kname3,"%s%d%s","B",500600 + addcut * 1000000,"ct5");
  sprintf(Kname4,"%s%d%s","B",500600 + addcut * 1000000,"ct6");

  setHist(((TH1D*)gDirectory->Get(Kname)) ,1,24,1,1,scal2);
  setHist(((TH1D*)gDirectory->Get(Kname2)),2,20,1,1,scal4);
  setHist(((TH1D*)gDirectory->Get(Kname3)),3,21,1,1,scal6);
  setHist(((TH1D*)gDirectory->Get(Kname4)),4,25,1,1,1);

  printHist(((TH1D*)gDirectory->Get(Kname)),((TH1D*)gDirectory->Get(Kname2)),((TH1D*)gDirectory->Get(Kname3)));
  gPad->Update();  tl.SetTextColor(1);
  tl.DrawLatex(0.22, 0.75, KSltg);
  tl.SetTextColor(2);
  tl.DrawLatex(0.22, 0.85, KSltg2);
  tl.SetTextColor(3);
  tl.DrawLatex(0.22, 0.65, KSltg3);gPad->Update();


  c5->cd(2); 
  sprintf(KSltg,"B^{0}: Ks->#pi#pi"); sprintf(KSltg2,"B^{0}: Ks->#pi^{0}#pi^{0}");  
  sprintf(KSltg3,"B^{0}: Kl");  
  sprintf(Kname,"%s%d%s" ,"B",400600 + addcut * 1000000,"ct3");
  sprintf(Kname2,"%s%d%s","B",400600 + addcut * 1000000,"ct4");
  sprintf(Kname3,"%s%d%s","B",400600 + addcut * 1000000,"ct5");
  sprintf(Kname4,"%s%d%s","B",400600 + addcut * 1000000,"ct6");

  setHist(((TH1D*)gDirectory->Get(Kname)) ,1,24,1,1,scal2);
  setHist(((TH1D*)gDirectory->Get(Kname2)),2,20,1,1,scal4);
  setHist(((TH1D*)gDirectory->Get(Kname3)),3,21,1,1,scal6);
  setHist(((TH1D*)gDirectory->Get(Kname4)),4,25,1,1,1);

  printHist(((TH1D*)gDirectory->Get(Kname)),((TH1D*)gDirectory->Get(Kname2)),((TH1D*)gDirectory->Get(Kname3)));
  gPad->Update();  tl.SetTextColor(1);
  tl.DrawLatex(0.22, 0.75, KSltg);
  tl.SetTextColor(2);
  tl.DrawLatex(0.22, 0.85, KSltg2);
  tl.SetTextColor(3);
  tl.DrawLatex(0.22, 0.65, KSltg3);gPad->Update();


  sprintf(name3,"%s%s%s%s",dir.Data(),"/",flag.Data(),"_All_Ksstu.ps");
  c5->Print(name3);  

}

void BkgbrekTool::Brek(int cut, TString dir, TString flag, int Norm){

  double scal1, scal2, scal3;
  int addcut = 0; char mm[100]; char myvar[100];
  char name3[100];
  if(cut) addcut = 1;
  tl.SetNDC(kTRUE); tl.SetTextSizePixels(25.);
  tl.SetTextAngle(0);
  
  c5 = new TCanvas("c5", "c5",600,800);   
  c5->Clear();  c5->cd(); 
  char  Kname[100], Kname3[100], Kname4[100], Kname2[100], Kname5[100], Kname6[100], Kname7[100], Kname8[100], Kname9[100], KSltg[100], Kname11[100], Kname12[100], Kname13[100];
  char  aKname[100], aKname3[100], aKname4[100], aKname2[100], aKname5[100], aKname6[100], aKname7[100], aKname8[100], aKname9[100], aKSltg[100], aKname11[100], aKname12[100], aKname13[100];

  sprintf(Kname8,"%s%d%s" ,"B",500600 + addcut * 1000000,"ct0");
  //  sprintf(Kname9,"%s%d%s" ,"B",500600 + addcut * 1000000,"ct2");
  sprintf(Kname,"%s%d%s" ,"B",500600 + addcut * 1000000,"ct3");
  sprintf(Kname2,"%s%d%s","B",500600 + addcut * 1000000,"ct4");
  sprintf(Kname3,"%s%d%s","B",500600 + addcut * 1000000,"ct5");
  sprintf(Kname4,"%s%d%s","B",500600 + addcut * 1000000,"ct6");
  sprintf(Kname5,"%s%d%s","B",500600 + addcut * 1000000,"ct7");
  sprintf(Kname6,"%s%d%s","B",500600 + addcut * 1000000,"ct8");
  sprintf(Kname7,"%s%d%s","B",500600 + addcut * 1000000,"ct9");
  sprintf(Kname11,"%s%d%s","B",500600 + addcut * 1000000,"ct11");
  sprintf(Kname12,"%s%d%s","B",500600 + addcut * 1000000,"ct12");
  sprintf(Kname13,"%s%d%s","B",500600 + addcut * 1000000,"ct13");
  sprintf(aKname8,"%s%d%s","B",400600 + addcut * 1000000,"ct0");
  //  sprintf(aKname9,"%s%d%s","B",400600 + addcut * 1000000,"ct2");
  sprintf(aKname,"%s%d%s" ,"B",400600 + addcut * 1000000,"ct3");
  sprintf(aKname2,"%s%d%s","B",400600 + addcut * 1000000,"ct4");
  sprintf(aKname3,"%s%d%s","B",400600 + addcut * 1000000,"ct5");
  sprintf(aKname4,"%s%d%s","B",400600 + addcut * 1000000,"ct6");
  sprintf(aKname5,"%s%d%s","B",400600 + addcut * 1000000,"ct7");
  sprintf(aKname6,"%s%d%s","B",400600 + addcut * 1000000,"ct8");
  sprintf(aKname7,"%s%d%s","B",400600 + addcut * 1000000,"ct9");
  sprintf(aKname11,"%s%d%s","B",400600 + addcut * 1000000,"ct11");
  sprintf(aKname12,"%s%d%s","B",400600 + addcut * 1000000,"ct12");
  sprintf(aKname13,"%s%d%s","B",400600 + addcut * 1000000,"ct13");


  double allDsemib0  =((TH1D*)gDirectory->Get(aKname11))->Integral();
  double allDsemibc  =((TH1D*)gDirectory->Get(Kname11))->Integral();
  double allDssemib0 =((TH1D*)gDirectory->Get(aKname12))->Integral();
  double allDssemibc =((TH1D*)gDirectory->Get(Kname12))->Integral();
  double allDsssemib0=((TH1D*)gDirectory->Get(aKname13))->Integral();
  double allDsssemibc=((TH1D*)gDirectory->Get(Kname13))->Integral();
  

  TH1D *hbg8 =  (TH1D*)gDirectory->Get(Kname8);
  TH1D *hsig8 = (TH1D*)gDirectory->Get(aKname8);
  TH1D *hsub8 = new TH1D("All B0", "", hsig8->GetNbinsX(), hsig8->GetBinLowEdge(1), hsig8->GetBinLowEdge(hsig8->GetNbinsX()+1)); 
  hsub8->Add(hsig8, hbg8, 1., 1.);

  TH1D *hbg =  (TH1D*)gDirectory->Get(Kname);
  TH1D *hsig = (TH1D*)gDirectory->Get(aKname);
  TH1D *hsub = new TH1D("All B1", "", hsig->GetNbinsX(), hsig->GetBinLowEdge(1), hsig->GetBinLowEdge(hsig->GetNbinsX()+1)); 
  hsub->Add(hsig, hbg, 1., 1.);

  TH1D *hbg2 =  (TH1D*)gDirectory->Get(Kname2);
  TH1D *hsig2 = (TH1D*)gDirectory->Get(aKname2);
  TH1D *hsub2 = new TH1D("All B2", "", hsig2->GetNbinsX(), hsig2->GetBinLowEdge(1), hsig2->GetBinLowEdge(hsig2->GetNbinsX()+1)); 
  hsub2->Add(hsig2, hbg2, 1., 1.);

  TH1D *hbg3 =  (TH1D*)gDirectory->Get(Kname3);
  TH1D *hsig3 = (TH1D*)gDirectory->Get(aKname3);
  TH1D *hsub3 = new TH1D("All B3", "", hsig3->GetNbinsX(), hsig3->GetBinLowEdge(1), hsig3->GetBinLowEdge(hsig3->GetNbinsX()+1)); 
  hsub3->Add(hsig3, hbg3, 1., 1.);

  TH1D *hbg4 =  (TH1D*)gDirectory->Get(Kname4);
  TH1D *hsig4 = (TH1D*)gDirectory->Get(aKname4);
  TH1D *hsub4 = new TH1D("All B4", "", hsig4->GetNbinsX(), hsig4->GetBinLowEdge(1), hsig4->GetBinLowEdge(hsig4->GetNbinsX()+1)); 
  hsub4->Add(hsig4, hbg4, 1., 1.);

  TH1D *hbg5 =  (TH1D*)gDirectory->Get(Kname5);
  TH1D *hsig5 = (TH1D*)gDirectory->Get(aKname5);
  TH1D *hsub5 = new TH1D("All B5", "", hsig5->GetNbinsX(), hsig5->GetBinLowEdge(1), hsig5->GetBinLowEdge(hsig5->GetNbinsX()+1)); 
  hsub5->Add(hsig5, hbg5, 1., 1.);

  TH1D *hbg6 =  (TH1D*)gDirectory->Get(Kname6);
  TH1D *hsig6 = (TH1D*)gDirectory->Get(aKname6);
  TH1D *hsub6 = new TH1D("All B6", "", hsig6->GetNbinsX(), hsig6->GetBinLowEdge(1), hsig6->GetBinLowEdge(hsig6->GetNbinsX()+1)); 
  hsub6->Add(hsig6, hbg6, 1., 1.);

  TH1D *hbg7 =  (TH1D*)gDirectory->Get(Kname7);
  TH1D *hsig7 = (TH1D*)gDirectory->Get(aKname7);
  TH1D *hsub7 = new TH1D("All B7", "", hsig7->GetNbinsX(), hsig7->GetBinLowEdge(1), hsig7->GetBinLowEdge(hsig7->GetNbinsX()+1)); 
  hsub7->Add(hsig7, hbg7, 1., 1.);

  //  TH1D *hbg9 =  (TH1D*)gDirectory->Get(Kname9);
  //  TH1D *hsig9 = (TH1D*)gDirectory->Get(aKname9);
  //  TH1D *hsub9 = new TH1D("All B8", "", hsig9->GetNbinsX(), hsig9->GetBinLowEdge(1), hsig9->GetBinLowEdge(hsig9->GetNbinsX()+1)); 
  //  hsub9->Add(hsig9, hbg9, 1., 1.);

  TH1D *hbg11 =  (TH1D*)gDirectory->Get(Kname11);
  TH1D *hsig11 = (TH1D*)gDirectory->Get(aKname11);
  TH1D *hsub11 = new TH1D("All B11", "", hsig11->GetNbinsX(), hsig11->GetBinLowEdge(1), hsig11->GetBinLowEdge(hsig11->GetNbinsX()+1)); 
  hsub11->Add(hsig11, hbg11, 1., 1.);

  TH1D *hbg12 =  (TH1D*)gDirectory->Get(Kname12);
  TH1D *hsig12 = (TH1D*)gDirectory->Get(aKname12);
  TH1D *hsub12 = new TH1D("All B12", "", hsig12->GetNbinsX(), hsig12->GetBinLowEdge(1), hsig12->GetBinLowEdge(hsig12->GetNbinsX()+1)); 
  hsub12->Add(hsig12, hbg12, 1., 1.);

  TH1D *hbg13 =  (TH1D*)gDirectory->Get(Kname13);
  TH1D *hsig13 = (TH1D*)gDirectory->Get(aKname13);
  TH1D *hsub13 = new TH1D("All B13", "", hsig13->GetNbinsX(), hsig13->GetBinLowEdge(1), hsig13->GetBinLowEdge(hsig13->GetNbinsX()+1)); 
  hsub13->Add(hsig13, hbg13, 1., 1.);

  setHist(((TH1D*)gDirectory->Get("All B0")),1,24,1,1,1);
  setHist(((TH1D*)gDirectory->Get("All B1")),2,25,1,1,1);
  setHist(((TH1D*)gDirectory->Get("All B2")),3,26,1,1,1);
  setHist(((TH1D*)gDirectory->Get("All B3")),4,20,1,1,1);
  setHist(((TH1D*)gDirectory->Get("All B4")),5,21,1,1,1);
  setHist(((TH1D*)gDirectory->Get("All B5")),6,22,1,1,1);
  setHist(((TH1D*)gDirectory->Get("All B6")),7,29,1,1,1);
  setHist(((TH1D*)gDirectory->Get("All B7")),8,30,1,1,1);
  //  setHist(((TH1D*)gDirectory->Get("All B8")),9,31,1,1,1);

  if(Norm == 0) {
    scal1 = 1;    scal2 = 1;
    scal3 = 1;  
  } else { //Normalize to area-events dc-bc
    double allD =    allDsemibc+allDsemib0;
    double allDs =    allDssemibc+allDssemib0;
    double allDss =    allDsssemibc+allDsssemib0;
    scal1 = 1/allD;  scal2 = 1/allDs;
    scal3 = 1/allDss;
  }


  setHist(((TH1D*)gDirectory->Get("All B11")),1,21,1,1,scal1);
  setHist(((TH1D*)gDirectory->Get("All B12")),2,25,1,1,scal2);
  setHist(((TH1D*)gDirectory->Get("All B13")),3,26,1,1,scal3);

  gPad->Update(); gStyle->SetOptStat(0); char nappa[100];
  double mb[8];   double ib[8];
  for(int fii=0; fii<8; fii++) { 
    sprintf(nappa,"%s%d","All B",fii);
    if(fii == 0) {
      ((TH1D*)gDirectory->Get(nappa))->Draw();  
    } else {
      ((TH1D*)gDirectory->Get(nappa))->Draw("same");
    }
    mb[fii] = ((TH1D*)gDirectory->Get(nappa))->GetMean();
    ib[fii] = ((TH1D*)gDirectory->Get(nappa))->Integral();
  }
  gPad->Update();

  TLegend *leg;
  TLegendEntry *legge; 
  leg = new TLegend(0.01,0.7,0.88,0.89);
  leg->SetFillStyle(0); leg->SetBorderSize(0.); leg->SetTextSize(0.05); 
  leg->SetFillColor(0); 
  legge = leg->AddEntry((TH1D*)gDirectory->Get("All B0"), "All", "p"); 
  legge = leg->AddEntry((TH1D*)gDirectory->Get("All B1"), "Ks #pi#pi", "p"); 
  legge->SetTextColor(2);
  legge = leg->AddEntry((TH1D*)gDirectory->Get("All B2"), "Ks #pi0#pi0", "p"); 
  legge->SetTextColor(3);
  legge = leg->AddEntry((TH1D*)gDirectory->Get("All B3"), "Kl", "p"); 
  legge->SetTextColor(4);
  legge = leg->AddEntry((TH1D*)gDirectory->Get("All B4"), "Additional #nu", "p"); 
  legge->SetTextColor(9);
  legge = leg->AddEntry((TH1D*)gDirectory->Get("All B5"), "K missed", "p"); 
  legge->SetTextColor(6);
  legge = leg->AddEntry((TH1D*)gDirectory->Get("All B6"), "K misid", "p"); 
  legge->SetTextColor(7);
  legge = leg->AddEntry((TH1D*)gDirectory->Get("All B7"), "other", "p"); 
  legge->SetTextColor(8);
  //  legge = leg->AddEntry((TH1D*)gDirectory->Get("All B8"), "Additional #nu", "p"); 
  //  legge->SetTextColor(9);
  leg->Draw();gPad->Update();
  
  sprintf(name3,"%s%s%s%s",dir.Data(),"/",flag.Data(),"_All_KsstuTotalm.ps");
  c5->Print(name3); 

  gPad->Update(); gStyle->SetOptStat(0);
  double mbs[3];   double ibs[3];
  for(int sii=0; sii<3; sii++) { 
    sprintf(nappa,"%s%d","All B",sii+11);
    mbs[sii] = ((TH1D*)gDirectory->Get(nappa))->GetMean();
    ibs[sii] = ((TH1D*)gDirectory->Get(nappa))->Integral();
  }
  ((TH1D*)gDirectory->Get("All B12"))->Draw();  
  ((TH1D*)gDirectory->Get("All B11"))->Draw("same");
  ((TH1D*)gDirectory->Get("All B13"))->Draw("same");
  gPad->Update();

  leg = new TLegend(0.01,0.7,0.88,0.89);
  leg->SetFillStyle(0); leg->SetBorderSize(0.); leg->SetTextSize(0.05); 
  leg->SetFillColor(0); 
  legge = leg->AddEntry((TH1D*)gDirectory->Get("All B11"), "Dl#nu", "p"); 
  legge = leg->AddEntry((TH1D*)gDirectory->Get("All B12"), "D^{*}l#nu", "p"); 
  legge->SetTextColor(2);
  legge = leg->AddEntry((TH1D*)gDirectory->Get("All B13"), "D^{**}l#nu", "p"); 
  legge->SetTextColor(3);
  leg->Draw();gPad->Update();
  
  sprintf(name3,"%s%s%s%s",dir.Data(),"/",flag.Data(),"_All_DsemilepTotalm.ps");
  c5->Print(name3); 
  char line[200];
  TLatex tl;  tl.SetNDC(kTRUE); tl.SetTextSize(0.08); tl.SetTextColor(1); 

  char tmpcan[100];
  c5 = new TCanvas("c5", "c5",600,800);   
  c5->Clear();  c5->Divide(3,3); 
  for(int ican=0; ican<8; ican++) {
    c5->cd(ican+1);
    //    gStyle->SetOptStat(001111);
    gStyle->SetOptStat(0);
    sprintf(tmpcan,"%s%d","All B",ican);
    ((TH1D*)gDirectory->Get(tmpcan))->Draw();
    leg = new TLegend(0.01,0.7,0.5,0.89);
    leg->SetFillStyle(0); leg->SetBorderSize(0.); leg->SetTextSize(0.05); 
    leg->SetFillColor(0); 

    if(ican == 0) {
      sprintf(KSltg,"I");     
    } else if(ican == 1) {
      sprintf(KSltg,"II");  
    } else if(ican == 2) {
      sprintf(KSltg,"III");  
    } else if(ican == 3) {
      sprintf(KSltg,"IV"); 
    } else if(ican == 4) {
      sprintf(KSltg,"V"); 
    } else if(ican == 5) {
      sprintf(KSltg,"VI"); 
    } else if(ican == 6) {
      sprintf(KSltg,"VII"); 
    } else if(ican == 7) {
      sprintf(KSltg,"VIII"); 
    }
    sprintf(line,"%s%5.3f%s%5.1f","M ",mb[ican],"; I ",ib[ican]);
    tl.DrawLatex(0.1, 0.91, line);
    legge = leg->AddEntry((TH1D*)gDirectory->Get(tmpcan),KSltg, "p"); 
    leg->SetTextSize(0.1); leg->Draw();gPad->Update();
  }
  sprintf(name3,"%s%s%s%s",dir.Data(),"/",flag.Data(),"_All_KsstuTotalmEach.ps");

  c5->Print(name3); 
  c5 = new TCanvas("c5", "c5",600,800);   
  c5->Clear();  c5->Divide(1,3); 
  for(int ican=0; ican<3; ican++) {
    c5->cd(ican+1);
    //    gStyle->SetOptStat(001111);
    gStyle->SetOptStat(0);
    sprintf(tmpcan,"%s%d","All B",11+ican);
    ((TH1D*)gDirectory->Get(tmpcan))->Draw();
    leg = new TLegend(0.01,0.7,0.5,0.89);
    leg->SetFillStyle(0); leg->SetBorderSize(0.); leg->SetTextSize(0.05); 
    leg->SetFillColor(0); 
    if(ican == 0) {
      sprintf(KSltg,"I");     
    } else if(ican == 1) {
      sprintf(KSltg,"II");  
    } else if(ican == 2) {
      sprintf(KSltg,"III");  
    }
    sprintf(line,"%s%5.3f%s%5.1f","M ",mbs[ican],"; I ",ibs[ican]);
    tl.DrawLatex(0.1, 0.91, line);
    legge = leg->AddEntry((TH1D*)gDirectory->Get(tmpcan),KSltg, "p"); 
    leg->SetTextSize(0.1); leg->Draw();gPad->Update();
  }
  sprintf(name3,"%s%s%s%s",dir.Data(),"/",flag.Data(),"_All_DsemilepTotalmEach.ps");
  c5->Print(name3); 
}
void BkgbrekTool::BrekB(int cut, TString dir, TString flag, int Norm, int fl){

  double scal1, scal2, scal3;
  int addcut = 0; char mm[100]; char myvar[100];
  char name3[100];
  if(cut) addcut = 1;
  tl.SetNDC(kTRUE); tl.SetTextSizePixels(25.);
  tl.SetTextAngle(0);
  
  c5 = new TCanvas("c5", "c5",600,800);   
  c5->Clear();  c5->cd(); 
  char KSltg2[100], KSltg[100], KSltg3[100], KSltg4[100], KSltg5[100], KSltg6[100], Kname[100], Kname3[100], Kname4[100], Kname2[100], Kname5[100], Kname6[100], Kname7[100], Kname8[100], Kname9[100], KSltg7[100], KSltg8[100], KSltg9[100], Kname11[100], Kname12[100], KSltg11[100], KSltg12[100], KSltg13[100], Kname13[100];

  sprintf(KSltg,"B All");      sprintf(KSltg3,"B; Ks->#pi^{0}#pi^{0}");  
  sprintf(KSltg2,"B; Ks->#pi#pi");   sprintf(KSltg4,"B; Kl"); 
  sprintf(KSltg5,"B; additional #nu");  sprintf(KSltg6,"B; K missed"); 
  sprintf(KSltg7,"B; misidentified");  sprintf(KSltg8,"B; ??"); 
  sprintf(KSltg11,"Dl#nu");  sprintf(KSltg12,"D^{*}l#nu"); sprintf(KSltg9,"B; additional #nu"); 
  sprintf(KSltg13,"D^{**}l#nu"); 
  int pino;
  if(fl) {
    pino =300000;
  } else {
    pino =400000;
  }
  sprintf(Kname8,"%s%d%s","B",100600 + pino + addcut * 1000000,"ct0");
  sprintf(Kname, "%s%d%s","B",100600 + pino + addcut * 1000000,"ct3");
  sprintf(Kname2,"%s%d%s","B",100600 + pino + addcut * 1000000,"ct4");
  sprintf(Kname3,"%s%d%s","B",100600 + pino + addcut * 1000000,"ct5");
  sprintf(Kname4,"%s%d%s","B",100600 + pino + addcut * 1000000,"ct6");
  sprintf(Kname5,"%s%d%s","B",100600 + pino + addcut * 1000000,"ct7");
  sprintf(Kname6,"%s%d%s","B",100600 + pino + addcut * 1000000,"ct8");
  sprintf(Kname7,"%s%d%s","B",100600 + pino + addcut * 1000000,"ct9");
  //  sprintf(Kname9,"%s%d%s","B",100600 + pino + addcut * 1000000,"ct2");

  sprintf(Kname11,"%s%d%s","B",100600 + pino + addcut * 1000000,"ct11");
  sprintf(Kname12,"%s%d%s","B",100600 + pino + addcut * 1000000,"ct12");
  sprintf(Kname13,"%s%d%s","B",100600 + pino + addcut * 1000000,"ct13");

  double allDsemi  =((TH1D*)gDirectory->Get(Kname11))->Integral();
  double allDssemi =((TH1D*)gDirectory->Get(Kname12))->Integral();
  double allDsssemi=((TH1D*)gDirectory->Get(Kname13))->Integral();

  setHist(((TH1D*)gDirectory->Get(Kname8)),1,24,1,1,1);
  setHist(((TH1D*)gDirectory->Get(Kname)),2,25,1,1,1);
  setHist(((TH1D*)gDirectory->Get(Kname2)),3,26,1,1,1);
  setHist(((TH1D*)gDirectory->Get(Kname3)),4,20,1,1,1);
  setHist(((TH1D*)gDirectory->Get(Kname4)),5,21,1,1,1);
  setHist(((TH1D*)gDirectory->Get(Kname5)),6,22,1,1,1);
  setHist(((TH1D*)gDirectory->Get(Kname6)),7,29,1,1,1);
  setHist(((TH1D*)gDirectory->Get(Kname7)),8,30,1,1,1);

  if(Norm == 0) {
    scal1 = 1;    scal2 = 1;
    scal3 = 1;  
  } else { //Normalize to area-events dc-bc
    double allD =    allDsemi;
    double allDs =    allDssemi;
    double allDss =    allDsssemi;
    scal1 = 1/allD;  scal2 = 1/allDs;
    scal3 = 1/allDss;
  }

  setHist(((TH1D*)gDirectory->Get(Kname11)),1,21,1,1,scal1);
  setHist(((TH1D*)gDirectory->Get(Kname12)),2,25,1,1,scal2);
  setHist(((TH1D*)gDirectory->Get(Kname13)),3,26,1,1,scal3);

  gPad->Update(); gStyle->SetOptStat(0);
  ((TH1D*)gDirectory->Get(Kname8))->Draw();  
  ((TH1D*)gDirectory->Get(Kname))->Draw("same");
  ((TH1D*)gDirectory->Get(Kname2))->Draw("same");
  ((TH1D*)gDirectory->Get(Kname3))->Draw("same"); 
  ((TH1D*)gDirectory->Get(Kname4))->Draw("same");
  ((TH1D*)gDirectory->Get(Kname5))->Draw("same"); 
  ((TH1D*)gDirectory->Get(Kname6))->Draw("same");
  ((TH1D*)gDirectory->Get(Kname7))->Draw("same");
  //  ((TH1D*)gDirectory->Get(Kname9))->Draw("same");
  gPad->Update();

  TLegend *leg;
  TLegendEntry *legge; 
  leg = new TLegend(0.01,0.7,0.88,0.89);
  leg->SetFillStyle(0); leg->SetBorderSize(0.); leg->SetTextSize(0.05); 
  leg->SetFillColor(0); 
  legge = leg->AddEntry((TH1D*)gDirectory->Get(Kname8), "All", "p"); 
  legge = leg->AddEntry((TH1D*)gDirectory->Get(Kname), "Ks #pi#pi", "p"); 
  legge->SetTextColor(2);		       
  legge = leg->AddEntry((TH1D*)gDirectory->Get(Kname2), "Ks #pi0#pi0", "p"); 
  legge->SetTextColor(3);		      
  legge = leg->AddEntry((TH1D*)gDirectory->Get(Kname3), "Kl", "p"); 
  legge->SetTextColor(4);		       
  legge = leg->AddEntry((TH1D*)gDirectory->Get(Kname4), "Additional #nu", "p"); 
  legge->SetTextColor(9);		       
  legge = leg->AddEntry((TH1D*)gDirectory->Get(Kname5), "K missed", "p"); 
  legge->SetTextColor(6);		       
  legge = leg->AddEntry((TH1D*)gDirectory->Get(Kname6), "K misid", "p"); 
  legge->SetTextColor(7);		       
  legge = leg->AddEntry((TH1D*)gDirectory->Get(Kname7), "other", "p"); 
  legge->SetTextColor(8);
  leg->Draw();gPad->Update();

  char * pinof;
  if(fl) {
    pinof ="b0";
  } else {
    pinof ="chb";
  }
  sprintf(name3,"%s%s%s%s%s%s",dir.Data(),"/",flag.Data(),"_All_KsstuTotalm_",pinof,".ps");
  c5->Print(name3);  

  gPad->Update(); gStyle->SetOptStat(0);
  ((TH1D*)gDirectory->Get(Kname12))->Draw();  
  ((TH1D*)gDirectory->Get(Kname11))->Draw("same");
  ((TH1D*)gDirectory->Get(Kname13))->Draw("same");
  gPad->Update();

  leg = new TLegend(0.01,0.7,0.88,0.89);
  leg->SetFillStyle(0); leg->SetBorderSize(0.); leg->SetTextSize(0.05); 
  leg->SetFillColor(0); 
  legge = leg->AddEntry((TH1D*)gDirectory->Get(Kname11), "Dl#nu", "p"); 
  legge = leg->AddEntry((TH1D*)gDirectory->Get(Kname12), "D^{*}l#nu", "p"); 
  legge->SetTextColor(2);		       
  legge = leg->AddEntry((TH1D*)gDirectory->Get(Kname13), "D^{**}l#nu", "p"); 
  legge->SetTextColor(3);		      
  leg->Draw();gPad->Update();

  sprintf(name3,"%s%s%s%s%s%s",dir.Data(),"/",flag.Data(),"_All_DsemilepTotalm_",pinof,".ps");
  c5->Print(name3);  

 char tmpcan[100];
  c5 = new TCanvas("c5", "c5",600,800);   
  c5->Clear();  c5->Divide(3,3); 
  for(int ican=0; ican<8; ican++) {
    c5->cd(ican+1);
    if(ican == 0) {
      sprintf(tmpcan,"%s%d%s","B",100600 + pino + addcut * 1000000,"ct0");
    } else if (ican == 1) {
      sprintf(tmpcan,"%s%d%s","B",100600 + pino + addcut * 1000000,"ct3");
    } else if (ican == 2) {
      sprintf(tmpcan,"%s%d%s","B",100600 + pino + addcut * 1000000,"ct4");
    } else if (ican == 3) {
      sprintf(tmpcan,"%s%d%s","B",100600 + pino + addcut * 1000000,"ct5");
    } else if (ican == 4) {
      sprintf(tmpcan,"%s%d%s","B",100600 + pino + addcut * 1000000,"ct6");
    } else if (ican == 5) {
      sprintf(tmpcan,"%s%d%s","B",100600 + pino + addcut * 1000000,"ct7");
    } else if (ican == 6) {
      sprintf(tmpcan,"%s%d%s","B",100600 + pino + addcut * 1000000,"ct8");
    }
    if(ican == 0) {
      sprintf(KSltg,"I");     
    } else if(ican == 1) {
      sprintf(KSltg,"II");  
    } else if(ican == 2) {
      sprintf(KSltg,"III");  
    } else if(ican == 3) {
      sprintf(KSltg,"IV"); 
    } else if(ican == 4) {
      sprintf(KSltg,"V"); 
    } else if(ican == 5) {
      sprintf(KSltg,"VI"); 
    } else if(ican == 6) {
      sprintf(KSltg,"VII"); 
    }
    //    gStyle->SetOptStat(001111);
    gStyle->SetOptStat(0);
    ((TH1D*)gDirectory->Get(tmpcan))->Draw();
    leg = new TLegend(0.01,0.7,0.4,0.89);
    leg->SetFillStyle(0); leg->SetBorderSize(0.); leg->SetTextSize(0.05); 
    leg->SetFillColor(0); 
    legge = leg->AddEntry((TH1D*)gDirectory->Get(tmpcan),KSltg, "p"); 
    leg->SetTextSize(0.1);leg->Draw();gPad->Update();
  }
  sprintf(name3,"%s%s%s%s%s%s",dir.Data(),"/",flag.Data(),"_All_KsstuTotalmEach_",pinof,".ps");
  c5->Print(name3); 

  c5 = new TCanvas("c5", "c5",600,800);   
  c5->Clear();  c5->Divide(1,3); 
  for(int ican=0; ican<3; ican++) {
    c5->cd(ican+1);
    if(ican == 0) {
      sprintf(tmpcan,"%s%d%s","B",100600 + pino + addcut * 1000000,"ct11");
    } else if (ican == 1) {
      sprintf(tmpcan,"%s%d%s","B",100600 + pino + addcut * 1000000,"ct12");
    } else if (ican == 2) {
      sprintf(tmpcan,"%s%d%s","B",100600 + pino + addcut * 1000000,"ct13");
    } 
    if(ican == 0) {
      sprintf(KSltg,"I");     
    } else if(ican == 1) {
      sprintf(KSltg,"II");  
    } else if(ican == 2) {
      sprintf(KSltg,"III");  
    }
    //    gStyle->SetOptStat(001111);
    gStyle->SetOptStat(0);
    ((TH1D*)gDirectory->Get(tmpcan))->Draw();
    leg = new TLegend(0.01,0.7,0.4,0.89);
    leg->SetFillStyle(0); leg->SetBorderSize(0.); leg->SetTextSize(0.05); 
    leg->SetFillColor(0); 
    legge = leg->AddEntry((TH1D*)gDirectory->Get(tmpcan),KSltg, "p"); 
    leg->SetTextSize(0.1);leg->Draw();gPad->Update();
  }
  sprintf(name3,"%s%s%s%s%s%s",dir.Data(),"/",flag.Data(),"_All_DsemilepTotalmEach_",pinof,".ps");
  c5->Print(name3); 

}

int BkgbrekTool::retDcat(int fDk,int fDks,int fDpi,int fDpiz) {
  
  int ip, ipiz, ik, iks, num;
  num = -99;
  for(int icep =0; icep <31; icep++) {
    ip = MA[icep][0]; ipiz = MA[icep][1]; ik = MA[icep][2]; iks = MA[icep][3];
    if((GfDpi == ip) && (GfDpiz == ipiz) && (GfDk == ik) && (GfDks == iks)) {
      num = icep;   
    }
  }
  return num;
}

int BkgbrekTool::retDScat(int vxbtyp) {
  
  int numS = -99;
  for(int icepS =0; icepS <9; icepS++) {
    if(TMath::Abs(vxbtyp) == icepS)  numS = icepS;   
  }
  return numS;
}


int BkgbrekTool::hist(double mx){

  // categories
  int bin = int((mx-themin)/((themax-themin)/thebins)+1);
  if (mx>themax || mx==themax) bin = thebins;
if (mx<themin || mx==themin) bin = 1;
  return bin;
}

int BkgbrekTool::hist2(double mnu){

  // categories
  int bin = int((mnu-themin2)/((themax2-themin2)/thebins2)+1);
  if (mnu>themax2 || mnu==themax2) bin = thebins2;
if (mnu<themin2 || mnu==themin2) bin = 1;
  return bin;
}

void BkgbrekTool::sighisto(double&signal, double& signalErr, TH1D *histo, double &resmean, double &ressigma, double &resalpha, double &resn, int fixpar, double mean, double sigma, double alpha,  double n, double argus){

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
void BkgbrekTool::dumpCuts() {
  cout << "====================================" << endl;
  cout << "Cut file " << fCutFile << endl; 
  cout << "------------------------------------" << endl;
  cout << "deSignal:          " << DESIGNALLO << " ... " << DESIGNALHI << endl;
  cout << "pcms:              " << CUTPCMS   << endl;
  cout << "intPurity:         " << INTPURITY << endl;
  cout << "Mx max cut:        " << MXMAXCUT  << endl;
  cout << "Mm2 cut:           " << MM2CUT    << endl;
  cout << "D from D:          " << DFROMD    << endl;
  cout << "====================================" << endl;
}

void BkgbrekTool::initRest() {
  // -- CUTS
  DESIGNALLO = -0.1; 
  DESIGNALHI =  0.1; 
  CUTPCMS    = 1.;
  INTPURITY  = 0.;
  MXMAXCUT  = 5.;
  MM2CUT  = .5;
  DFROMD  = 1;
}

// ----------------------------------------------------------------------
void BkgbrekTool::readCuts(TString filename, int dump) {
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



TH1D* BkgbrekTool::BGsub(char *hist, TCut tmpcutpa, double min, double max, char *dvar) {

  //Trying to get those bkg sub histos
  recoilAnalysis a;
  char name[200];   char var[200];
  //Plot of Mes  
  TH1 *pino = new TH1D("mes", "Mes distribution", 100, 5.21, 5.29); pino->Sumw2();
  ((TTree*)gDirectory->Get("events"))->Draw("mes >> pino",tmpcutpa);
  TH1D *hmes = (TH1D*)gDirectory->Get("pino");
  double bgScale = a.getBgScale(hmes);

  TCut bgcut = "(5.21 <= mes) && (mes <= 5.26)";
  TCut sgcut = "((5.27 <= mes) && (mes <= 5.29)) && ((-0.1 <= de) && (de <= 0.1))";
  TCut bkgcut = bgcut && tmpcutpa;
  TCut sigcut = sgcut && tmpcutpa;

  sig = new TH1D("sig", "X distribution", 100, min, max); sig->Sumw2();
  bkg = new TH1D("bkg", "X distribution", 100, min, max); bkg->Sumw2();

  sprintf(var,"%s%s",dvar," >> bkg"); ((TTree*)gDirectory->Get("events"))->Draw(var,bkgcut,"l");
  TH1D *hbg =  (TH1D*)gDirectory->Get("bkg");

  sprintf(var,"%s%s",dvar," >> sig"); ((TTree*)gDirectory->Get("events"))->Draw(var,sigcut,"l");
  TH1D *hsig = (TH1D*)gDirectory->Get("sig");

  sprintf(name, "sub_%s", hist);
  TH1D *hsub = new TH1D(name, "", hsig->GetNbinsX(), hsig->GetBinLowEdge(1), hsig->GetBinLowEdge(hsig->GetNbinsX()+1)); 
  hsub->Add(hsig, hbg, 1., -1.*bgScale);
  return hsub;
}


void BkgbrekTool::Init(TTree *tree)
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
   fChain->SetBranchAddress("GfDnu",&GfDnu); 
   fChain->SetBranchAddress("GfDpiz",&GfDpiz);
   fChain->SetBranchAddress("GfDkmiss",&GfDkmiss);  
   fChain->SetBranchAddress("GfDk",&GfDk);  
   fChain->SetBranchAddress("GfDks",&GfDks); 
   fChain->SetBranchAddress("GfDkl",&GfDkl); 
   fChain->SetBranchAddress("GfDkspipi",&GfDkspipi); 
   fChain->SetBranchAddress("GfDkspiopio",&GfDkspiopio); 
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
   fChain->SetBranchAddress("wdeltam",&wdeltam);
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

Bool_t BkgbrekTool::Notify()
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
   b_GfDnu=      fChain->GetBranch("GfDnu");  
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
   b_isDupli = fChain->GetBranch("isDupli");
   b_ValMap = fChain->GetBranch("ValMap");
   b_vub = fChain->GetBranch("vub");
   b_vcb = fChain->GetBranch("vcb");
   b_vxbtyp = fChain->GetBranch("vxbtyp");
   b_other = fChain->GetBranch("other");
   b_bgcat = fChain->GetBranch("bgcat");
   b_wdeltam = fChain->GetBranch("wdeltam");
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
void BkgbrekTool::setHist(TH1 *h, Int_t color, Int_t symbol, Double_t size, Double_t width, Double_t scale) {
  h->SetLineColor(color);   h->SetLineWidth(width);
  h->SetMarkerColor(color); h->SetMarkerStyle(symbol);  h->SetMarkerSize(size); 
  //  h->SetStats(kFALSE);  
  if(scale != 0) { 
    h->Scale(scale); 
  }
  //  h->SetFillStyle(0); h->SetFillColor(color);
}


// ----------------------------------------------------------------------
void BkgbrekTool::printHist(TH1 *h1, TH1 *h2, TH1 *h3) {

  double m1 = h1->GetMaximum();
  double m2 = h2->GetMaximum();
  double m3 = h3->GetMaximum();
  cout<<m1<<" "<<m2<<" "<<m3<<endl;
  double conf;

  if(m1 > m2 && m1 >m3 && m2 > m3) conf = 0;
  if(m1 > m2 && m1 >m3 && m3 > m2) conf = 1;
  if(m2 > m1 && m2 >m3 && m1 > m3) conf = 2;
  if(m2 > m1 && m2 >m3 && m3 > m1) conf = 3;
  if(m3 > m1 && m3 >m2 && m2 > m1) conf = 4;
  if(m3 > m1 && m3 >m2 && m1 > m2) conf = 5;
  gPad->Update();
  if(conf == 0) {
    h1->Draw();
    h2->Draw("same");
    h3->Draw("same");
  } else if (conf == 1) {
    h1->Draw();
    h3->Draw("same");
    h2->Draw("same");
  } else if (conf == 2) {
    h2->Draw();
    h1->Draw("same");
    h3->Draw("same");
  } else if (conf == 3) {
    h2->Draw();
    h3->Draw("same");
    h1->Draw("same");
  } else if (conf == 4) {
    h3->Draw();
    h2->Draw("same");
    h1->Draw("same");
  } else if (conf == 5) {
    h3->Draw();
    h1->Draw("same");
    h2->Draw("same");
  }
  gPad->Update();
}
