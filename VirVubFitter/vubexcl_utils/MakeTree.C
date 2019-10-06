#define MakeTree_cxx
#include "MakeTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void MakeTree::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L MakeTree.C
//      Root > MakeTree t
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

   char name[100];
   Int_t modebestEta, modebestEtap;
   Double_t Rho0massdaubestEtap, EtamassdaubestEtap;

   TFile *f = new TFile("RooRootFiles/mydata.root","RECREATE");

   TTree *ntp2 = new TTree("ntp2","test tree");
   // Breco 
   ntp2->Branch("mes",&mes,"mes/D");
   ntp2->Branch("ass_deltapB",&ass_deltapB,"ass_deltapB/D");
   ntp2->Branch("brecocharge",&brecocharge,"brecocharge/I");
   ntp2->Branch("brecoflav",&brecoflav,"brecoflav/I");
   // Best lepton
   ntp2->Branch("nle",&nle,"nle/I");
   ntp2->Branch("nlept500",&nle,"nlept500/I"); 
   ntp2->Branch("lcharge",&lcharge,"lcharge/I");
   ntp2->Branch("isele",&isele,"isele/I");
   ntp2->Branch("pcms",&pcms,"pcms/D");
   ntp2->Branch("plab",&plab,"plab/D");
   ntp2->Branch("tlab",&tlab,"tlab/D");
   // MC truth
   ntp2->Branch("Gvxbtyp",&Gvxbtyp,"Gvxbtyp/I");
   // ETA L NU
   ntp2->Branch("nrecoEta",&nrecoEta,"nrecoEta/I");
   ntp2->Branch("indexbestEta",&indexbestEta,"indexbestEta/I");
   ntp2->Branch("modebestEta",&modebestEta,"modebestEta/I");
   ntp2->Branch("barembestEta",&barembestEta,"barembestEta/D");  
   // ETA PRIME L NU
   ntp2->Branch("nrecoEtap",&nrecoEtap,"nrecoEtap/I");
   ntp2->Branch("indexbestEtap",&indexbestEtap,"indexbestEtap/I");
   ntp2->Branch("modebestEtap",&modebestEtap,"modebestEtap/I");
   ntp2->Branch("barembestEtap",&barembestEtap,"barembestEtap/D");
   ntp2->Branch("Rho0massdaubestEtap",&Rho0massdaubestEtap,"Rho0massdaubestEtap/D");
   ntp2->Branch("EtamassdaubestEtap",&EtamassdaubestEtap,"EtamassdaubestEtap/D");

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (jentry%10000 == 0) cout << "Event " << jentry << "\n";

      modebestEta = -999;
      modebestEtap = -999;
      Rho0massdaubestEtap = -999;
      EtamassdaubestEtap = -999;

      if(nrecoEta>0&&indexbestEta>-999) modebestEta = modeEta[indexbestEta];
      if(nrecoEtap>0&&indexbestEtap>-999){ 
	modebestEtap = modeEtap[indexbestEtap];
	Rho0massdaubestEtap = Rho0massdauEtap[indexbestEtap];
        EtamassdaubestEtap = EtamassdauEtap[indexbestEtap];
      }

      ntp2->Fill();

   }

   f->cd();
   ntp2->Write();
   f->Close();

}
