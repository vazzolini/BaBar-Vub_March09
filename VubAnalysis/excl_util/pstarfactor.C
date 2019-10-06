#define pstarfactor_cxx
#include <fstream.h>
#include "../RecoilAnalysis/mesData.hh"
#include "../RecoilAnalysis/recoilAnalysis.hh"
#include "pstarfactor.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"

void pstarfactor::Loop(int leptype, int imode, int isvub)
{
//   In a ROOT session, you can do:
//      Root > .L pstarfactor.C
//      Root > pstarfactor t
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


  char name[100],flavcat[100];
  bool islept;
  double cutflav = 0;    
  double pscut;

  lepton = leptype;

  fHistFile->cd();

   if (fChain == 0) return;

   Int_t nentries = Int_t(fChain->GetEntries());
   //   nentries = 10000;

   cout << "N events: " << nentries << endl;

   int countvub(0);

   Int_t nbytes = 0, nb = 0;
   for (Int_t jentry=0; jentry<nentries;jentry++) {
      Int_t ientry = LoadTree(jentry); //in case of a TChain, ientry is the entry number in the current file
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (jentry%1000 == 0) cout << "Event " << jentry << endl;
      // if (Cut(ientry) < 0) continue;

      islept = ((nle > 0) && (tlab<2.37) && (tlab>0.36) && (plab>0.5));  
      if(leptype == 0) islept = ((nel > 0) && (tlab<2.37) && (tlab>0.36) && (plab>0.5));
      if(leptype == 1) islept = ((nmu > 0) && (tlab<2.37) && (tlab>0.36) && (plab>0.5));

      bool Bmode = (brecocharge == 0);
      if(imode<0) Bmode = !(brecocharge == 0);
      if(Bmode==0) continue;
      
      if(pur<.06) continue;

      if((vcb||vub) && !isvub){
	((TH1D*)gDirectory->Get("mesvcb"))->Fill(mes);
      }
      if(Gvxbtyp==imode){
	((TH1D*)gDirectory->Get("mesvub"))->Fill(mes);
      }

      if(islept && !(TMath::Abs(brecocharge)!=0 && (lcharge + brecoflav)!=0)){
	 int flav =  lcharge + brecoflav; // charge correlation
	 // flavor category (3 = charged B, 4 = neutral B OS, 5 = neutral B SS)  
	 cutflav = 5;
	 if(TMath::Abs(brecocharge)) cutflav = 3;
	 if(TMath::Abs(brecocharge)==0 && flav==0) cutflav = 4;	  
	 
	 if(cutflav==3) sprintf(flavcat, "%s", "bch");
	 if(cutflav==4) sprintf(flavcat, "%s", "bos");
	 if(cutflav==5) sprintf(flavcat, "%s", "bss");

	 for (int j=10;j<26;j++){
	   
	   pscut = j/10.;

	   if(pcms > pscut){
	     if((vcb||vub) && !isvub){
	       sprintf(name, "%s%s%d","mesvcblep",flavcat,j);
	       ((TH1D*)gDirectory->Get(name))->Fill(mes);	       	   
	     }
	     if(Gvxbtyp==imode){
	       sprintf(name, "%s%s%d","mesvublep",flavcat,j);
	       ((TH1D*)gDirectory->Get(name))->Fill(mes);	       	   
	     }	   
	   } 
	 }
      }

   }

}

void pstarfactor::Bookhist(int imode)
{
  
  char name[100],title[100];
  sprintf(name,"%s%d%s","output_",imode,".root");
  fHistFile = new TFile(name, "RECREATE");
  fHistFile->cd();
  
  TH1 *h;
  sprintf(name,"mesvub");  sprintf(title, "mes vub breco");  h = new TH1D(name, title, 40, 5.2, 5.3 );  h->Sumw2();
  sprintf(name,"mesvcb");  sprintf(title, "mes vcb breco");  h = new TH1D(name, title, 40, 5.2, 5.3 );  h->Sumw2();

  for(int i=10; i<26; i++){
     sprintf(name,"%s%d","mesvublepbch",i);  sprintf(title, "mes vub lep bch");  h = new TH1D(name, title, 40, 5.2, 5.3 );  h->Sumw2();
     sprintf(name,"%s%d","mesvcblepbch",i);  sprintf(title, "mes vcb lep bch");  h = new TH1D(name, title, 40, 5.2, 5.3 );  h->Sumw2();
     sprintf(name,"%s%d","mesvublepbos",i);  sprintf(title, "mes vub lep bos");  h = new TH1D(name, title, 40, 5.2, 5.3 );  h->Sumw2();
     sprintf(name,"%s%d","mesvcblepbos",i);  sprintf(title, "mes vcb lep bos");  h = new TH1D(name, title, 40, 5.2, 5.3 );  h->Sumw2();
     sprintf(name,"%s%d","mesvublepbss",i);  sprintf(title, "mes vub lep bss");  h = new TH1D(name, title, 40, 5.2, 5.3 );  h->Sumw2();
     sprintf(name,"%s%d","mesvcblepbss",i);  sprintf(title, "mes vcb lep bss");  h = new TH1D(name, title, 40, 5.2, 5.3 );  h->Sumw2();
  }

  sprintf(name,"pstarfact");  sprintf(title, "pstarfactor vs pstar cut");  h = new TH1D(name, title, 16, 0.95, 2.55 );  h->Sumw2();  

}

void pstarfactor::FitMes()
{
  fHistFile->cd();
  char name[100], namebch[100], nameb0os[100], nameb0ss[100];
 
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
  mesData themes;
  double thesigma = -1111111;
  double themean =  -1111111;
  double theargus = -1111111;
  double thealpha = -1111111;
  double then = 5.;
  double resmean, ressigma, resalpha, resn;
  recoilAnalysis k;

  mesData tempthemes(*(k.vubMes(((TH1D*)gDirectory->Get("mesvcb")), resmean, ressigma, resalpha, resn, 1,1, themean, thesigma, thealpha, then, theargus)));
  nvcb=tempthemes.theSig(); nvcberr=tempthemes.theErrSig(); 
  cout << "nvcb = " << nvcb << " +- " <<  nvcberr;

  mesData tempthemes2(*(k.vubMes(((TH1D*)gDirectory->Get("mesvub")), resmean, ressigma, resalpha, resn, 1,1, themean, thesigma, thealpha, then, theargus)));
  nvub=tempthemes2.theSig(); nvuberr=tempthemes2.theErrSig(); 
  cout << "nvub = " << nvub << " +- " <<  nvuberr;
  

  for (int i = 10;  i < 26; i++) { 
    
    sprintf(namebch, "%s%d","mesvcblepbch",i);  
    sprintf(nameb0os, "%s%d","mesvcblepbos",i);  
    sprintf(nameb0ss, "%s%d","mesvcblepbss",i);  
    
    // mixing correction      
    
    mesData tempthemes3(*(k.vubMes(((TH1D*)gDirectory->Get(namebch)), resmean, ressigma, resalpha, resn, 1,1, themean, thesigma, thealpha, then, theargus)));
    tempbinchb=tempthemes3.theSig(); temperrchb=tempthemes3.theErrSig(); 
    cout << "nvcbbch" << i-10 << " = " << tempbinchb << " +- " << temperrchb;

    mesData tempthemes4(*(k.vubMes(((TH1D*)gDirectory->Get(nameb0os)), resmean, ressigma, resalpha, resn, 1,1, themean, thesigma, thealpha, then, theargus)));
    tempbinb0os=tempthemes4.theSig(); temperrb0os=tempthemes4.theErrSig(); 
    cout << "nvcbb0os" << i-10 << " = " << tempbinb0os << " +- " << temperrb0os;

    mesData tempthemes5(*(k.vubMes(((TH1D*)gDirectory->Get(nameb0ss)), resmean, ressigma, resalpha, resn, 1,1, themean, thesigma, thealpha, then, theargus)));
    tempbinb0ss=tempthemes5.theSig(); temperrb0ss=tempthemes5.theErrSig(); 
    cout << "nvcb0ss" << i-10 << " = " << tempbinb0ss << " +- " << temperrb0ss;

    tempbin = tempbinchb + ((1-chid)/(1-2*chid)) * tempbinb0os - (chid/(1-2*chid)) * tempbinb0ss;
    temperr = sqrt(temperrchb*temperrchb + ((1-chid)/(1-2*chid)) * ((1-chid)/(1-2*chid)) * tempbinb0os + (chid/(1-2*chid)) * (chid/(1-2*chid))* tempbinb0ss);
    
    nvcblep[i-10] = tempbin;
    nvcbleperr[i-10] = temperr;
    cout << "nvcb" << i-10 << " = " << nvcblep[i-10] << " +- " << nvcbleperr[i-10];
  }

  for (int i = 10;  i < 26; i++) { 
    
    sprintf(namebch, "%s%d","mesvublepbch",i);  
    sprintf(nameb0os, "%s%d","mesvublepbos",i);  
    sprintf(nameb0ss, "%s%d","mesvublepbss",i);  
    
    // mixing correction      
    
    mesData tempthemes(*(k.vubMes(((TH1D*)gDirectory->Get(namebch)), resmean, ressigma, resalpha, resn, 1,1, themean, thesigma, thealpha, then, theargus)));
    tempbinchb=tempthemes.theSig(); temperrchb=tempthemes.theErrSig(); 
    cout << "nvubbch" << i-10 << " = " << tempbinchb << " +- " << temperrchb;

    mesData tempthemes2(*(k.vubMes(((TH1D*)gDirectory->Get(nameb0os)), resmean, ressigma, resalpha, resn, 1,1, themean, thesigma, thealpha, then, theargus)));
    tempbinb0os=tempthemes2.theSig(); temperrb0os=tempthemes2.theErrSig(); 
    cout << "nvubb0os" << i-10 << " = " << tempbinb0os << " +- " << temperrb0os;

    mesData tempthemes3(*(k.vubMes(((TH1D*)gDirectory->Get(nameb0ss)), resmean, ressigma, resalpha, resn, 1,1, themean, thesigma, thealpha, then, theargus)));
    tempbinb0ss=tempthemes3.theSig(); temperrb0ss=tempthemes3.theErrSig(); 
    cout << "nvub0ss" << i-10 << " = " << tempbinb0ss << " +- " << temperrb0ss;

    tempbin = tempbinchb + ((1-chid)/(1-2*chid)) * tempbinb0os - (chid/(1-2*chid)) * tempbinb0ss;
    temperr = sqrt(temperrchb*temperrchb + ((1-chid)/(1-2*chid)) * ((1-chid)/(1-2*chid)) * tempbinb0os + (chid/(1-2*chid)) * (chid/(1-2*chid))* tempbinb0ss);
    
    nvublep[i-10] = tempbin;
    nvubleperr[i-10] = temperr;
    cout << "nvub" << i-10 << " = " << nvublep[i-10] << " +- " << nvubleperr[i-10];
  }

}

void pstarfactor::Finalize(int imode)
{

  TCanvas c0("c0","--c0--",472,0,800,900);
  double temppstar, temppstarerr, tempvubeff, tempcut, totvub, totvcb;
  char name[100];
  cout << endl;
  cout << endl;
  cout << "breco number of Vub events = " << nvub << " +- " << nvuberr << endl;
  cout << "breco number of Vcb events = " << nvcb << " +- " << nvcberr << endl;

  cout << endl;
  cout << endl;
  
  sprintf(name, "%s%d%s%d","pstar.dat_",lepton,"_",imode);        
  ofstream outfile(name);

  for (int i = 25;  i > 9; i--) { 
    temppstar = (nvublep[i-10]/nvub) / (nvcblep[i-10]/nvcb); 
    tempvubeff = nvublep[i-10]/nvub;
    temppstarerr = (sqrt(tempvubeff * (1-tempvubeff)/nvub) * (nvuberr/sqrt(nvub))) / tempvubeff;
    temppstarerr = temppstarerr * temppstar;
    tempcut = i/10.;
    
    cout << " pcms > " << tempcut << "    PSTAR FACTOR = " << temppstar << " +- " << temppstarerr << endl;
    outfile << tempcut << "  " << temppstar << endl;
    ((TH1D*)gDirectory->Get("pstarfact"))->SetBinContent(i-9,temppstar);
    ((TH1D*)gDirectory->Get("pstarfact"))->SetBinError(i-9,temppstarerr);
     
  }
  
  outfile.close();
  ((TH1D*)gDirectory->Get("pstarfact"))->SetStats(0);
  ((TH1D*)gDirectory->Get("pstarfact"))->Draw("pe");

  sprintf(name,"%s%d%s%d","genepsf_",lepton,"_",imode);   
  c0.SaveAs(name);

}
