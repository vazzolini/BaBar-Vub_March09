#define lumi_cxx
#include "lumi.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"

void lumi::Loop()
{
//   In a Root session, you can do:
//      Root > .L lumi.C
//      Root > lumi t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(i);  // read all branches
//by  b_branchname->GetEntry(i); //read only this branch
   if (fChain == 0) return;
   int runs[30000];
   int goodruns[30000];
   double mx[1000];
   for (int i=0; i< 30000; i++) runs[i]=0;
   Int_t nentries = Int_t(fChain->GetEntries());
   ofstream outfile("runlumi.dat");	          
   ofstream outfile2("runlumi_missing.dat");	          
   ofstream outfile3("runlumi_bad.dat");	          
   int temprun = 0;
   Int_t nbytes = 0, nb = 0;

   int irun = 0;
   float therun;
   char times[100];
   char tableName[1000], buffer[200], fname[100];

   TH1D goodrun("goodrun","goodrun",30000,0.5,30000.5);
   TH1D myrun("myrun","myrun",30000,.5,30000.5);
   TH1D missingrun("missingrun","missingrun",30000,.5,30000.5);
   
   sprintf(buffer, "%s", "/afs/slac.stanford.edu/g/babar/www/Physics/BaBarData/GoodRuns/good_all.txt");

   ifstream is(buffer);

   while (is.getline(buffer, 200, '\n')) {
     if (buffer[0] == '#') {continue;}
     sscanf(buffer, "%f %s", &therun, &times);
     goodruns[irun] = therun;     
     goodrun.Fill(therun);
     irun++;
   }   
   goodrun.Draw();
   cout << "total runs" << goodrun.Integral()<< endl;
   c0.SaveAs("goodrun.ps");

   for (Int_t jentry=0; jentry<nentries;jentry++) {
      Int_t ientry = LoadTree(jentry); //in case of a TChain, ientry is the entry number in the current file
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;      
      if(run!=temprun) {
	outfile << run  << endl;
	temprun = run;
	runs[run]++;
      }

   } 

   for (int i=1; i< 30000; i++){
     if(runs[i]>1) cout << "run " << i << " is duplicated" << endl;
     if(runs[i]>0) myrun.Fill(i);
     
   }
   myrun.Draw();
   cout << "my runs" << myrun.Integral()<< endl;
   c0.SaveAs("myrun.ps");

   outfile.close();

   missingrun.Add(&goodrun,1);
   missingrun.Add(&myrun,-1); 

   missingrun.Draw();
   cout << "missing runs" << missingrun.Integral()<< endl;
   c0.SaveAs("missingrun.ps");

   for (int i=1; i< 30000; i++){
     if(missingrun.GetBinContent(i)>.5)   outfile2 << i  << endl;
     if(missingrun.GetBinContent(i)<-.5)   outfile3 << i  << endl;
   }
   outfile2.close();
   outfile3.close();
   
}
