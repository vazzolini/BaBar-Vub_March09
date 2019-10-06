#!/bin/csh -f

bbrroot -b<<EOF

gROOT->SetStyle("Plain");

TFile g1("/nfs/farm/babar/AWGsemilep01/gallof/fit/fitetap/data/testdatafitresult.root");

//TFile g1("mesplot.root");

gROOT->Macro("makefinalplots.C");

.q
