#define MesPlots_cxx
#include "MesPlots.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void MesPlots::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L meshisto.C
//      Root > meshisto t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Int_t nbytes = 0, nb = 0;

   char name[100], title[100];
   sprintf(name,"mesplot.root");

   fHistFile = new TFile(name, "RECREATE");
   fHistFile->cd();

   sprintf(name,"allcutseta"); sprintf(title, "m_{ES} for B^{+} #rightarrow #eta l #nu"); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
   sprintf(name,"allcutsetap"); sprintf(title, "m_{ES} for B^{+} #rightarrow #eta' l #nu"); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();
   sprintf(name,"nsldata"); sprintf(title, "m_{ES} for B^{+} after lepton cut"); h = new TH1D(name, title, 40, 5.2, 5.3); h->Sumw2();

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (jentry%10000 == 0) cout << "Event " << jentry << "\n";

      Bool_t lPYes = (TMath::Abs(brecocharge) == 1)&&!(TMath::Abs(brecocharge)!=0&&(lcharge + brecoflav)!=0)&&(nlept500>0)&&(nle<=1)&&(tlab<2.37)&&(tlab>0.36)&&(plab>0.5)&&((isele==1&&pcms>0.5)||(isele==0&&pcms>0.8));

      Bool_t AllCutnomm2Eta = (xcharge+brecocharge)==0&&(nrecoEta>0)&&((modeEta[indexbestEta]==1&&nchg==1&&barembestEta>0.505&&barembestEta<0.585)||(modeEta[indexbestEta]==2&&nchg==3&&barembestEta>0.53&&barembestEta<0.56)||(modeEta[indexbestEta]==3&&nchg==1&&barembestEta>0.51&&barembestEta<0.58))&&(mm2bestPi0>1.5||mm2bestPi0<-10.);

      Bool_t AllCutnomm2Etap = (xcharge+brecocharge)==0&&(nrecoEtap>0)&&((modeEtap[indexbestEtap]==1&&nchg==3&&barembestEtap>0.93&&barembestEtap<0.98&&TMath::Abs(Rho0massdauEtap[indexbestEtap]-0.775)<0.18)||(modeEtap[indexbestEtap]==2&&nchg==3&&barembestEtap>0.94&&barembestEtap<0.97&&TMath::Abs(EtamassdauEtap[indexbestEtap]-0.54775)<0.04)||(modeEtap[indexbestEtap]==3&&nchg==5&&barembestEtap>0.935&&barembestEtap<0.975&&TMath::Abs(EtamassdauEtap[indexbestEtap]-0.54775)<0.015)||(modeEtap[indexbestEtap]==4&&nchg==3&&barembestEtap>0.91&&barembestEtap<1.00&&TMath::Abs(EtamassdauEtap[indexbestEtap]-0.54775)<0.035))&&(GammamomdauEtap[indexbestEtap]>0.35||GammamomdauEtap[indexbestEtap]<-10.);

      Bool_t mm2Etacut = (mm2bestEta>-0.5&&mm2bestEta<0.5);

      Bool_t mm2Etapcut = (mm2bestEtap>-0.5&&mm2bestEtap<0.5);

      if(lPYes){   
	sprintf(name,"nsldata");
        ((TH1D*)gDirectory->Get(name))->Fill(mes);
	if(AllCutnomm2Eta&&mm2Etacut){
	  sprintf(name,"allcutseta");
	  ((TH1D*)gDirectory->Get(name))->Fill(mes);
	}
	if(AllCutnomm2Etap&&mm2Etapcut){
	  sprintf(name,"allcutsetap");
	  ((TH1D*)gDirectory->Get(name))->Fill(mes);
	}
      }

   }

  fHistFile->Write();
  fHistFile->Close();
  delete fHistFile;

}
