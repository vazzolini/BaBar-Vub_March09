#include <fstream.h>
#include <iomanip.h>
#include <math.h>
#include <stdlib.h>

#include "TLatex.h"
#include "theanal.hh"
#include <iostream.h>

int main(int argc, char *argv[]) {

  TString fileName, inputFile(""), inputChain(""), dirName(""),cutFile(""), varName;
  Int_t bins(0);
  Double_t min(0), max(10); 
  char * end;
  // -- command line arguments
  for (int i = 0; i < argc; i++){
    if(strcmp(argv[i],"-F")  == 0) {inputFile    = TString(argv[++i]); }    
    if(strcmp(argv[i],"-c")  == 0) {inputChain   = TString(argv[++i]); }    
    if(strcmp(argv[i],"-C")  == 0) {cutFile      = TString(argv[++i]); }    
    if(strcmp(argv[i],"-D")  == 0) {dirName      = TString(argv[++i]);  }
    if(strcmp(argv[i],"-V")  == 0) {varName      = TString(argv[++i]);  }
    if(strcmp(argv[i],"-b")  == 0) bins          = atoi(argv[++i]);                // printout
    if(strcmp(argv[i],"-min")  == 0) min         = strtod(argv[++i],&end);                // printout
    if(strcmp(argv[i],"-max")  == 0) max         = strtod(argv[++i],&end);                // printout
  }
  
  cout << inputFile.Data() << endl;
  cout << cutFile.Data() << endl;
  cout << dirName.Data() << endl;
  cout << varName.Data() << endl;
  cout << bins << endl;
  cout << min << endl;
  cout << max << endl;

  TChain *chain = new TChain("ntp1");
  if(strcmp(inputFile,"")){
    chain->Add(inputFile);
  }
  if(strcmp(inputChain,"")){
     cout << "Chaining ... " << inputChain << endl;
     char pName[2000]; 
     char buffer[200];
     sprintf(buffer, "%s", inputChain.Data());
     ifstream is(buffer);  
     cout << buffer << endl;
     while(is.getline(buffer, 200, '\n')){
       if (buffer[0] == '#') continue; 
       sscanf(buffer, "%s", pName); 
       cout << "   Add: " << buffer << endl; 
       chain->Add(pName); 
     }
     is.close();     
  }

  char var[100],dir[100];
  sprintf(var, "%s",varName.Data());
  sprintf(dir, "%s",dirName.Data());
  theanal h(var, dir, min, max, bins);              
  h.Init(chain);
  h.readCuts(cutFile);
  h.Bookhist();             

  h.Loop(100000000);                
  h.Fitmes(1,0);
  h.Fitmes(2,0);
  h.Fitmes(3,0);
  h.Fitmes(4,0);
  h.Fitmes(1,1);
  h.Fitmes(2,1);
  h.Fitmes(3,1);
  h.Fitmes(4,1);
  
  h.overlap(0,1.,dirName);
  h.overlap(1,1.,dirName);

}
