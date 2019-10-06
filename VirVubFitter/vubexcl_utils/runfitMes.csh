#!/bin/csh -f

rm -rf RooTxtFiles/fitsettings.txt
echo $1 > RooTxtFiles/fitsettings.txt

bbrroot -b<<EOF

gSystem->Load("libPhysics.so");
gROOT->SetStyle("Plain");

ifstream samplefile;
samplefile.open("RooTxtFiles/fitsettings.txt");
char sample[100];
samplefile >> sample;
if(!samplefile.good()) break;
samplefile.close();

TChain *newtree = new TChain("ntp1");

// SIGNAL SP8
if(!strcmp(sample,"signal")) newtree->Add("/nfs/farm/babar/AWGsemilep01/dorazioa/signalSP8/SP-633*.root");

// COCKTAIL SP8
if(!strcmp(sample,"cocktail")) newtree->Add("/nfs/farm/babar/AWGsemilep01/dorazioa/cocktailSP8/SP-2223-*.root");

// GENERIC SP5-SP6
if(!strcmp(sample,"generic")) newtree->Add("/nfs/farm/babar/AWG24/francesco/test_prod/root/genb0/*.root");
if(!strcmp(sample,"generic")) newtree->Add("/nfs/farm/babar/AWG24/francesco/test_prod/root/genbch/*.root"); 

cout << "LOADING DATA... DONE\n";

TFile *newfile = new TFile("RooRootFiles/mergeddata.root","recreate");
newtree->Write();
newfile->Close();
cout << "MERGING DATA... DONE\n";

gROOT->LoadMacro("MakeTree.C");
TFile g("RooRootFiles/mergeddata.root");
MakeTree pluto(ntp1);
pluto.Loop();
cout << "SELECTING DATA FOR FIT... DONE\n";

cout << "STARTING FIT...\n";
gROOT->Macro("fitMes.C");
cout << "FIT DONE\n";

.q
