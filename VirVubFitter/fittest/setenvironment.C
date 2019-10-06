
{
  cout<<"!=======================>  Loading RooFit Libraries <================ "<<endl;
  gSystem->Load("../../RooFitCore/tmp/libRooFitCore.so");
  gSystem->Load("../../RooFitModels/tmp/libRooFitModels.so");

  //  cout<<"!=======================>  Loading Gauz Stuff...    <=============== "<<endl;
  //  gSystem->Load("gauz_cxx.so");
  //  cout<<"!=======================>  Loading CCB Stuff...    <=============== "<<endl;
  //  gSystem->Load("ccb_cxx.so");
  cout<<"!=======================>  Loading Thosig Stuff...    <=============== "<<endl;
  gSystem->Load("thosig_cxx.so");
  
  
  using namespace RooFit;
}
