#include <fstream.h>
#include <iomanip.h>
#include <math.h>
#include <stdlib.h>
#include "TROOT.h"


#include "TLatex.h"
#include "thecomp.hh"
#include <iostream.h>



int main(int argc, char *argv[]) {

  TString fileName, inputFileDat(""), inputChainDat(""), inputFileMC(""),inputChainMC(""), dirName(""),cutFile(""), varName;
  Int_t bins(0);
  Double_t min(0), max(10); 
  char * end;
  // -- command line arguments
  for (int i = 0; i < argc; i++){
    if(strcmp(argv[i],"-FD")  == 0) {inputFileDat    = TString(argv[++i]); }    
    if(strcmp(argv[i],"-cD")  == 0) {inputChainDat   = TString(argv[++i]); }    
    if(strcmp(argv[i],"-FMC")  == 0) {inputFileMC    = TString(argv[++i]); }    
    if(strcmp(argv[i],"-cMC")  == 0) {inputChainMC   = TString(argv[++i]); }    
    if(strcmp(argv[i],"-C")  == 0) {cutFile      = TString(argv[++i]); }    
    if(strcmp(argv[i],"-D")  == 0) {dirName      = TString(argv[++i]);  }
    if(strcmp(argv[i],"-V")  == 0) {varName      = TString(argv[++i]);  }
    if(strcmp(argv[i],"-b")  == 0) bins          = atoi(argv[++i]);                // printout
    if(strcmp(argv[i],"-min")  == 0) min         = strtod(argv[++i],&end);                // printout
    if(strcmp(argv[i],"-max")  == 0) max         = strtod(argv[++i],&end);                // printout
  }
  
  cout << inputFileDat.Data() << endl;
  cout << inputFileMC.Data() << endl;
  cout << cutFile.Data() << endl;
  cout << dirName.Data() << endl;
  cout << varName.Data() << endl;
  cout << bins << endl;
  cout << min << endl;
  cout << max << endl;

  TChain *chainDat = new TChain("events");
  if(strcmp(inputFileDat,"")){
    chainDat->Add(inputFileDat);
  }
  TChain *chainMC = new TChain("events");
  if(strcmp(inputFileMC,"")){
    chainMC->Add(inputFileMC);
  }
  char pName[2000]; 
  char buffer[200];
  if(strcmp(inputChainDat,"")){
     cout << "Chaining ... " << inputChainDat << endl;
     sprintf(buffer, "%s", inputChainDat.Data());
     ifstream is(buffer);  
     cout << buffer << endl;
     while(is.getline(buffer, 200, '\n')){
       if (buffer[0] == '#') continue; 
       sscanf(buffer, "%s", pName); 
       cout << "   Add: " << buffer << endl; 
       chainDat->Add(pName); 
     }
     is.close();     
  }
  if(strcmp(inputChainMC,"")){
     cout << "Chaining ... " << inputChainMC << endl;
     sprintf(buffer, "%s", inputChainMC.Data());
     ifstream is2(buffer);  
     cout << buffer << endl;
     while(is2.getline(buffer, 200, '\n')){
       if (buffer[0] == '#') continue; 
       sscanf(buffer, "%s", pName); 
       cout << "   Add: " << buffer << endl; 
       chainMC->Add(pName); 
     }
     is2.close();     
  }

  char var[100],dir[100];

  sprintf(var, "%s",varName.Data());
  sprintf(dir, "%s",dirName.Data());
  thecomp h(var, dir, min, max, bins);              
  h.Init(chainDat);
  h.readCuts(cutFile);
  h.Bookhist();             
  h.Loop(1000000000,1);
  h.Init(chainMC);
  h.Loop(1000000000,2);                
  cout<<"h101201 "<<((TH1D*)gDirectory->Get("h101201"))->GetEntries()<<endl;
  cout<<"h101202 "<<((TH1D*)gDirectory->Get("h101202"))->GetEntries()<<endl;
  cout<<"h102201 "<<((TH1D*)gDirectory->Get("h102201"))->GetEntries()<<endl;
  cout<<"h102202 "<<((TH1D*)gDirectory->Get("h102202"))->GetEntries()<<endl;
  //TCanvas c3;
  //c3.Divide(2,2);
  //c3.cd(1);
  //((TH1D*)gDirectory->Get("h101201"))->Draw();
  //c3.cd(2);
  //((TH1D*)gDirectory->Get("h101202"))->Draw();
  //c3.cd(3);
  //((TH1D*)gDirectory->Get("h102201"))->Draw();
  //c3.cd(4);
  //((TH1D*)gDirectory->Get("h102202"))->Draw();
  //c3.SaveAs("mesplots.ps");
  h.Fitmes(1,1);
  h.Fitmes(2,1);
  
  h.overlap(1,dirName);


}
