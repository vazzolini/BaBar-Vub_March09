#include <iostream>
#include "toyclass.h"
#include "TString.h"
using namespace std;

class toyclass;

int main( int argc, char* argv[]) {
  if(argc < 10) {
    cout << "------> Usage: " << endl;
    cout << "./toygen -start 10 -stop 20 -sfd 2.5 -sfg 2.5 -sfs 1.5" << endl;
    cout << endl;
    return 0;
  }

  Int_t start = 0;
  Int_t stop = 0;  
  Double_t sfd(0), sfg(0), sfs(0);
  
  for(int i = 1; i < argc; i++) {
    if( strcmp(argv[i],"-start") == 0) start = atoi(argv[++i]);
    if( strcmp(argv[i],"-stop") == 0) stop = atoi(argv[++i]);
    if( strcmp(argv[i],"-sfd") == 0) sfd = atof(argv[++i]);
    if( strcmp(argv[i],"-sfg") == 0) sfg = atof(argv[++i]);
    if( strcmp(argv[i],"-sfs") == 0) sfs = atof(argv[++i]);
  }

  string input;
  
  //   cout << "Provide full path for root file with mes Fits saved " << endl;
  //   cout << " (default ./PDFFile.root )" << endl;
  //   cin >> input;
  //   if(input.size() == 0)
  //     input = "PDFFile.root";
  
  input = "PDFfile.root";  
  cout << "from main input is " << input << endl;

  toyclass *a = new toyclass(TString(input));

  //  Double_t NEW_LUMI_DATA(0), NEW_LUMI_GENERIC(0), NEW_LUMI_SIGNAL(0);


  cout << " Lumi DATA in PDFfile.root is " << a->LUMI_DATA_OLD << " fb^-1" << endl;
  //  cout << " --->   Please input new LUMI for DATA (in fb^-1) ";
  sfd = 10;
  cout << " Scale Factor for data is " << sfd << endl;
  cout << endl;

  cout << " Lumi GENERIC in PDFfile.root (Valid for PStar, Vcb and VcbOther samples) is " << a->LUMI_GENERIC_OLD*1e-6 << " fb^-1" << endl;
  //  cout << " --->   Please input new LUMI for PStar, Vcb and VcbOth (in fb^-1) ";
  sfg = 4;
  cout << " Scale Factor for MC Generic is " << sfg << endl;
  cout << endl;


  cout << " Lumi GENERIC in PDFfile.root (Valid for VubIN and VubOUT) is " << a->LUMI_SIGNAL_OLD*1e-6 << " fb^-1 " << endl;
  //  cout << " --->   Please input new LUMI for VubIN and VubOUT (in fb^-1) ";
  sfs = 10;
  cout << " Scale Factor is for Signal MC is  " << sfs << endl;
  cout << endl;

  a->Loop(sfd, sfg, sfs, start, stop);
  
  return 0;
}
