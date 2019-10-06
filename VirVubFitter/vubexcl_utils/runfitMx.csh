#!/bin/csh -f
bbrroot -b<<EOF

gSystem->Load("libPhysics.so");
gROOT->SetStyle("Plain");

TChain *newtree = new TChain("ntp1");
newtree->Add("/nfs/farm/babar/AWG24/francesco/test_prod/root/signal/SP-3617-BSemiExcl-R14-*.root");
newtree->Add("/nfs/farm/babar/AWG24/francesco/test_prod/root/signal/SP-3618-BSemiExcl-R14-*.root");
newtree->Add("/nfs/farm/babar/AWG24/francesco/test_prod/root/eta/SP-4759-new-Run1-Run4-1.root");
newtree->Add("/nfs/farm/babar/AWG24/francesco/test_prod/root/etap/SP-4760-new-Run1-Run4-1.root");
cout << "LOADING DATA... DONE\n";

TFile *newfile = new TFile("RooRootFiles/mydata.root","recreate");
newtree->Write();
newfile->Close();
cout << "MERGING DATA... DONE\n";

gROOT->LoadMacro("MakeTree.C");
TFile g("RooRootFiles/mydata.root");
MakeTree pluto(ntp1);
pluto.Loop();
cout << "SELECTING DATA FOR FIT... DONE\n";

cout << "STARTING FIT...\n";
gROOT->Macro("fitMx.C");
cout << "FIT DONE\n";

.q
