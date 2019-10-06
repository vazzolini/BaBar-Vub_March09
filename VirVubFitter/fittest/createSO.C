
void createSO(){
  // LOAD ROOFIT LIBRARIES 
  gSystem->Load("../../RooFitCore/tmp/libRooFitCore.so");
  gSystem->Load("../../RooFitModels/tmp/libRooFitModels.so");
  
  //SETTING ROOFIT INCLUDE PATH
  gSystem->AddIncludePath("-I../../RooFitCore/tmp/");
  gSystem->AddIncludePath("-I../../RooFitModels/tmp/");

  // COMPILING PDF's
  //  gROOT->ProcessLine(".L gauz.cxx++");
  //  gROOT->ProcessLine(".L ccb.cxx++");
  gROOT->ProcessLine(".L thosig.cxx++");
}
