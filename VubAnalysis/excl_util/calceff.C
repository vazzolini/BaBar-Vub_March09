#define calceff_cxx
#include <fstream.h>
#include "calceff.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"

void calceff::Loop(int ID, int submode)
{
//   In a ROOT session, you can do:
//      Root > .L calceff.C
//      Root > calceff t
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
  
  
  char name[100],namevar[100],title[100],axistitle[100],submodename[100],submodestring[100],rootname[100];
  double rangemomtight = 0.2;
  double rangephitight = 0.025;
  double rangereso = 0.03;
  double lowbinmass = 0.11;
  double highbinmass = 0.16;
  int bins = 50;
  int binsmodes = 1;

  sprintf(submodename,"");
  sprintf(submodestring,"");

  if(ID == 11 || ID == -13) {
    rangemomtight = 0.1; 
  }
  if(ID == 13 || ID == -12 || ID == -14 || ID == -15 || ID == 19 || ID == -19) {
    rangephitight = 0.05;
  }
  if(ID == 13 || ID == -14 || ID == -12 || ID == -15 || ID == 19) rangereso = 0.08;
  if(ID == -19)  rangereso = 0.11;
  if(ID == 11){
    sprintf(name,"pi");      
    sprintf(namevar,"#pi");      
  }
  if(ID == -11){
    sprintf(name,"pi0");      
    sprintf(namevar,"#pi^{0}");      
  }
  if(ID == -13){
    sprintf(name,"rho0");      
    sprintf(namevar,"#rho^{0}");      
    lowbinmass = 0.3;
    highbinmass = 1.5;      
  }
  if(ID == 13){
    sprintf(name,"rho");      
    sprintf(namevar,"#rho");      
    lowbinmass = 0.3;
    highbinmass = 1.5;      
  }
  if(ID == -14){
    sprintf(name,"omega");      
    sprintf(namevar,"#omega");      
    lowbinmass = 0.68;
    highbinmass = 0.88;      
  }
  if(ID == -12){
    sprintf(name,"eta");      
    sprintf(namevar,"#eta");      
    if(submode==1) sprintf(submodename,"#eta #rightarrow #gamma#gamma");
    if(submode==2) sprintf(submodename,"#eta #rightarrow #pi#pi#pi^{0}");
    if(submode==3) sprintf(submodename,"#eta #rightarrow #pi^{0}#pi^{0}#pi^{0}");
    lowbinmass = 0.45;
    highbinmass = 0.65;      
    binsmodes = 3;
  }
  if(ID == -15){
    sprintf(name,"etap");      
    sprintf(namevar,"#eta'");      
    if(submode==1) sprintf(submodename,"#eta' #rightarrow #rho^{0}#gamma");
    if(submode==2) sprintf(submodename,"#eta' #rightarrow #eta#pi#pi,#eta #rightarrow #gamma#gamma");
    if(submode==3) sprintf(submodename,"#eta' #rightarrow #eta#pi#pi,#eta #rightarrow #pi#pi#pi^{0}");
    if(submode==4) sprintf(submodename,"#eta' #rightarrow #eta#pi#pi,#eta #rightarrow #pi^{0}#pi^{0}#pi^{0}");
    lowbinmass = 0.86;
    highbinmass = 1.06;      
    binsmodes = 4;
  }
  if(ID == -19){
    sprintf(name,"a0");      
    sprintf(namevar,"a^{0}_{0}");      
    if(submode==1) sprintf(submodename,"a^{0}_{0} #rightarrow #eta#pi^{0},#eta #rightarrow #gamma#gamma");
    if(submode==2) sprintf(submodename,"a^{0}_{0} #rightarrow #eta#pi^{0},#eta #rightarrow #pi#pi#pi^{0}");
    if(submode==3) sprintf(submodename,"a^{0}_{0} #rightarrow #eta#pi^{0},#eta #rightarrow #pi^{0}#pi^{0}#pi^{0}");
    lowbinmass = 0.75;
    highbinmass = 1.2;      
    binsmodes = 3;
    bins = 25;
  }
  if(ID == 19){
    sprintf(name,"a0p");      
    sprintf(namevar,"a^{+}_{0}");      
    if(submode==1) sprintf(submodename,"a^{0}_{+} #rightarrow #eta#pi,#eta #rightarrow #gamma#gamma");
    if(submode==2) sprintf(submodename,"a^{0}_{+} #rightarrow #eta#pi,#eta #rightarrow #pi#pi#pi^{0}");
    if(submode==3) sprintf(submodename,"a^{0}_{+} #rightarrow #eta#pi,#eta #rightarrow #pi^{0}#pi^{0}#pi^{0}");
    lowbinmass = 0.75;
    highbinmass = 1.2;      
    binsmodes = 3;
    bins = 25;
  }

  double rangemom = rangemomtight * 3;
  double rangephi = rangephitight * 3;


  sprintf(rootname,"%s%s%i%s","output_",name,submode,".root");
  fHistFile = new TFile(rootname, "RECREATE");
  fHistFile->cd();

  if(submode!=0) sprintf(submodestring,"%s%s%s"," (",submodename,")");
  sprintf(title,"%s%s%s","reconstructed mass for ",namevar,submodestring);
  TH1D *massplot = new TH1D("mass",title,bins,lowbinmass,highbinmass);
  sprintf(axistitle,"%s%s%s","mass(",namevar,")(GeV)");  
  massplot->SetXTitle(axistitle);

  sprintf(title,"%s%s%s","mass resolution for ",namevar,submodestring);
  TH1D *resomassplot = new TH1D("resomass",title,bins,-rangereso,rangereso);
  sprintf(axistitle,"%s%s%s","#sigma_{mass(",namevar,")}(GeV)");  
  resomassplot->SetXTitle(axistitle);

  sprintf(title,"%s%s%s","#theta resolution for ",namevar,submodestring);
  TH1D *resothetaplot = new TH1D("resotheta",title,bins,-rangephi,rangephi);
  sprintf(axistitle,"%s%s%s","#sigma_{#theta(",namevar,")}(GeV)");  
  resothetaplot->SetXTitle(axistitle);

  sprintf(title,"%s%s%s","#phi resolution for ",namevar,submodestring);
  TH1D *resophiplot = new TH1D("resophi",title,bins,-rangephi,rangephi);
  sprintf(axistitle,"%s%s%s","#sigma_{#phi(",namevar,")}(GeV)");  
  resophiplot->SetXTitle(axistitle);

  sprintf(title,"%s%s%s","mom resolution for ",namevar,submodestring);
  TH1D *resomomplot = new TH1D("resomom",title,bins,-rangemom,rangemom);
  sprintf(axistitle,"%s%s%s","#sigma_{mom(",namevar,")}(GeV)");  
  resomomplot->SetXTitle(axistitle);

  sprintf(title,"%s%s","submodes composition for ",namevar);
  TH1D *submodesplot = new TH1D("submodesplot",title,binsmodes,0.5,0.5+binsmodes);
  sprintf(axistitle,"%s%s%s","recon. modes(",namevar,")");  
  submodesplot->SetXTitle(axistitle);

  TH1D *entr = new TH1D("entr","",1,0.5,1.5);  
  

  if (fChain == 0) return;
  
  Int_t nentries = Int_t(fChain->GetEntries());
  
  Int_t nbytes = 0, nb = 0;

  int nevents(0);

  for (Int_t jentry=0; jentry<nentries;jentry++) {
    Int_t ientry = LoadTree(jentry); //in case of a TChain, ientry is the entry number in the current file
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%10000 == 0) cout << "Event " << jentry << endl;
    // if (Cut(ientry) < 0) continue;

    bool cuts = pcms>1 && mes>5.27 && brecocharge == 0;
    if(ID < 0) cuts = pcms>1 && mes>5.27 && brecocharge != 0;
    double themass = mpi0;
    double thetheta = thpi0;
    double thephi = phipi0;
    double themom = mompi0;
    double themode(submode); 
    if(ID == 11){
      themass = .139;
      thetheta = thpi;
      thephi = phipi;
      themom = mompi;
    }
    if(ID == -11){
      lowbinmass = 0.11;
      highbinmass = 0.16;      
    }
    if(ID == -13){
      themass = mrho0;
      thetheta = thrho0;
      thephi = phirho0;
      themom = momrho0;
    }
    if(ID == 13){
      themass = mrho;
      thetheta = thrho;
      thephi = phirho;
      themom = momrho;
    }
    if(ID == -14){
      themass = momega;
      thetheta = thomega;
      thephi = phiomega;
      themom = momomega;
    }
    if(ID == -12){
      themass = meta;
      thetheta = theta;
      thephi = phieta;
      themom = mometa;
      if(submode!=0) themode = modeeta;
    }
    if(ID == -15){
      themass = metap;
      thetheta = thetap;
      thephi = phietap;
      themom = mometap;
      if(submode!=0) themode = modeetap;
    }
    if(ID == 19){
      themass = ma0p;
      thetheta = tha0p;
      thephi = phia0p;
      themom = moma0p;
      if(submode!=0) themode = modea0p;
    }
    if(ID == -19){
      themass = ma0;
      thetheta = tha0;
      thephi = phia0;
      themom = moma0;
      if(submode!=0) themode = modea0;
    }
    
    if(cuts && Gvxbtyp == ID) {
      nevents++;
      entr->Fill(1.);
      if(themass>0 && TMath::Abs(thetheta-txhadgen)<rangephi && TMath::Abs(thephi-fxhadgen)<rangephi && TMath::Abs(themom-pxhadgen)<rangemom){

	if(TMath::Abs(thetheta-txhadgen)<rangephitight && TMath::Abs(thephi-fxhadgen)<rangephitight && TMath::Abs(themom-pxhadgen)<rangemomtight) 
	  submodesplot->Fill(themode);

	if(themode == submode){
	  if(TMath::Abs(thetheta-txhadgen)<rangephitight && TMath::Abs(thephi-fxhadgen)<rangephitight && TMath::Abs(themom-pxhadgen)<rangemomtight) {
	    massplot->Fill(themass);                  
	    resomassplot->Fill(themass-mxhadgen);
	  }
	  if(TMath::Abs(thephi-fxhadgen)<rangephitight && TMath::Abs(themom-pxhadgen)<rangemomtight) resothetaplot->Fill(thetheta-txhadgen);
	  if(TMath::Abs(thetheta-txhadgen)<rangephitight && TMath::Abs(themom-pxhadgen)<rangemomtight) resophiplot->Fill(thephi-fxhadgen);
	  if(TMath::Abs(thetheta-txhadgen)<rangephitight && TMath::Abs(thephi-fxhadgen)<rangephitight) resomomplot->Fill(themom-pxhadgen);
	}

      }
    }

  }
  
  cout << endl;
  cout << "Total number of events is: " << nevents << endl;
  cout << endl;

  gROOT->SetStyle("Plain");
  TCanvas *c1 = new TCanvas("c1"," ",200,10,800,800);   
  c1->Clear();

  sprintf(title,"%s%s%i%s","mass_",name,submode,".eps");
  massplot->Draw();
  c1->SaveAs(title);

  sprintf(title,"%s%s%i%s","resomass_",name,submode,".eps");
  resomassplot->Draw();
  c1->SaveAs(title);

  sprintf(title,"%s%s%i%s","resotheta_",name,submode,".eps");
  resothetaplot->Draw();
  c1->SaveAs(title);

  sprintf(title,"%s%s%i%s","resophi_",name,submode,".eps");
  resophiplot->Draw();
  c1->SaveAs(title);

  sprintf(title,"%s%s%i%s","resomom_",name,submode,".eps");
  resomomplot->Draw();
  c1->SaveAs(title);

  sprintf(title,"%s%s%s","submodes_",name,".eps");
  submodesplot->Draw();
  c1->SaveAs(title);

  fHistFile->Write();
  fHistFile->Close();


}

