#!/bin/csh -f

bbrroot -b<<EOF

gROOT->SetStyle("Plain");

gSystem->Load("/u/ec/gallof/vol2/testfit-ana31/shlib/$BFARCH/libRecoilAnalysis.so");

TFile g1("/nfs/farm/babar/AWGsemilep01/gallof/fit/fitetap/data/testdatafitresult.root");

//TFile g1("mesplot.root");

gROOT->Macro("mesplot.C");

.q
