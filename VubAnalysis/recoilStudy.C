#define recoilStudy_cxx
#include "recoilStudy.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLorentzVector.h"

recoilStudy::recoilStudy(){
}

void recoilStudy::Loop()
{
//   In a Root session, you can do:
//      Root > .L recoilStudy.C
//      Root > recoilStudy t
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
   char name[100];

   Int_t nentries = Int_t(fChain->GetEntries());

   Int_t nbytes = 0, nb = 0;

   sprintf(name, "%s%s%s", thedir.Data(),"/", histfile.Data()); 
   fHistFile = new TFile(name, "RECREATE");
   fHistFile->cd();
   cout << "Opened " << fHistFile->GetName() << endl;

   bookHist();

   char name[100];
   TLorentzVector p4t(0., 0., 0., 0.);
   TLorentzVector p4Recoil(0., 0., 0., 0.);

   for (Int_t jentry=0; jentry<nentries;jentry++) {
      Int_t ientry = LoadTree(jentry); //in case of a TChain, ientry is the entry number in the current file
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (jentry%1000 == 0) cout << "Event " << jentry << endl;
      if (ientry == 0) cout << "File " << fChain->GetCurrentFile()->GetName() << endl;

     // fill indexes vector
      
      PIPMASS = 0.13957; 
      PIZMASS = 0.13498;
      KAPMASS = 0.49368;
      KAZMASS = 0.49767;
      DSTARMASS  = 2.010;
      BPMASS  = 5.2790;
      BZMASS  = 5.2794;
      BMASS   = 5.2792;
      BQMASS  = 4.800;
      A0      = 1.29;

      findbestB();

      int modeB = modeB0[indexbestB];
      if(bestB0==0) {
	modeB = modeChB[indexbestB]; 
      }  

      if(bestB0) { 
	mes = mseB0[indexbestB]; 
	indexdau[0] = d1B0Index[indexbestB] - 1; 
	indexdau[1] = d2B0Index[indexbestB] - 1; 
	indexdau[2] = d3B0Index[indexbestB] - 1; 
	indexdau[3] = d4B0Index[indexbestB] - 1; 
	indexdau[4] = d5B0Index[indexbestB] - 1; 
	indexdau[5] = d6B0Index[indexbestB] - 1; 
	indexdau[6] = d7B0Index[indexbestB] - 1; 
	lunddau[0] = d1B0Lund[indexbestB]; 
	lunddau[1] = d2B0Lund[indexbestB]; 
	lunddau[2] = d3B0Lund[indexbestB]; 
	lunddau[3] = d4B0Lund[indexbestB]; 
	lunddau[4] = d5B0Lund[indexbestB]; 
	lunddau[5] = d6B0Lund[indexbestB]; 
	lunddau[6] = d7B0Lund[indexbestB]; 
      }else{
	mes = mseChB[indexbestB]; 
	indexdau[0] = d1ChBIndex[indexbestB] - 1; 
	indexdau[1] = d2ChBIndex[indexbestB] - 1; 
	indexdau[2] = d3ChBIndex[indexbestB] - 1; 
	indexdau[3] = d4ChBIndex[indexbestB] - 1; 
	indexdau[4] = d5ChBIndex[indexbestB] - 1; 
	indexdau[5] = d6ChBIndex[indexbestB] - 1; 
	indexdau[6] = d7ChBIndex[indexbestB] - 1; 
	lunddau[0] = d1ChBLund[indexbestB]; 
	lunddau[1] = d2ChBLund[indexbestB]; 
	lunddau[2] = d3ChBLund[indexbestB]; 
	lunddau[3] = d4ChBLund[indexbestB]; 
	lunddau[4] = d5ChBLund[indexbestB]; 
	lunddau[5] = d6ChBLund[indexbestB]; 
	lunddau[6] = d7ChBLund[indexbestB]; 
      }	

      for (int i=1; i<7; i++){

	mk4Vector(p4t, 0.,0.,0.,0.);
	if(TMath::Abs(lunddau[i]) == 211) mk4Vector(p4t, momentumTrk[indexdau[i]], thetaTrk[indexdau[i]], phiTrk[indexdau[i]], PIPMASS);
	if(TMath::Abs(lunddau[i]) == 321) mk4Vector(p4t, momentumTrk[indexdau[i]], thetaTrk[indexdau[i]], phiTrk[indexdau[i]], KAPMASS);
	if(TMath::Abs(lunddau[i]) == 111) mk4Vector(p4t, pPi0[indexdau[i]], thPi0[indexdau[i]], phiPi0[indexdau[i]], PIZMASS);
	if(TMath::Abs(lunddau[i]) == 310) mk4Vector(p4t, pKs[indexdau[i]], thKs[indexdau[i]], phiKs[indexdau[i]], KAZMASS);
      	p4Recoil += p4t;
	
      }
      
      if (bestB0) { 
	mk4Vector(p4t, pB0[indexbestB], thB0[indexbestB], phiB0[indexbestB], BZMASS);
      }else{
	mk4Vector(p4t, pChB[indexbestB], thChB[indexbestB], phiChB[indexbestB], BPMASS);	
      }
      //      p4Recoil += p4t;
      mk4Vector(p4t, pDstar[indexdau[0]], thDstar[indexdau[0]], phiDstar[indexdau[0]], DSTARMASS);
      //p4Recoil += p4t;
      
      double massrecoil = p4Recoil.Mag();
      
      fHistFile->cd();
       
      fHistFile->cd();  sprintf(name, "h%d", modeB); ((TH1D*)gDirectory->Get(name))->Fill(mes);
      if(mes > 5.27) {
	fHistFile->cd("sg"); 
	sprintf(name, "mx%d", modeB); ((TH1D*)gDirectory->Get(name))->Fill(massrecoil);
	sprintf(name, "mxfit%d", modeB); ((TH1D*)gDirectory->Get(name))->Fill(massrecoil);
      }
      if(mes > 5.21 && mes < 5.26){
	fHistFile->cd("bg"); 
	sprintf(name, "mx%d", modeB); ((TH1D*)gDirectory->Get(name))->Fill(massrecoil);
	sprintf(name, "mxfit%d", modeB); ((TH1D*)gDirectory->Get(name))->Fill(massrecoil);
      }
      // if (Cut(ientry) < 0) continue;
      p4Recoil.SetPx(0.);
      p4Recoil.SetPy(0.);
      p4Recoil.SetPz(0.);
      p4Recoil.SetE(0.);
   }
   fHistFile->cd();
   fHistFile->Write();
   fHistFile->Close();
   delete fHistFile;
}

void recoilStudy::readintpur(){    

   char tableName[1000], buffer[200], fname[100];
   float bmode, dmode, sig, bkg, pur, sb;
   ifstream is("tables/tablepurity.dat");
   int mode;
   while (is.getline(buffer, 200, '\n')) {
     if (buffer[0] == '#') {continue;}
     sscanf(buffer, "%f %f %f %f %f %f", &bmode, &dmode, &sig, &bkg, &pur, &sb);
     mode = (dmode+100) * 100 + bmode-10000;
     brecosig[mode] = sig;	
     brecobkg[mode] = bkg;
     brecointpur[mode] = pur; 
   }
	
}

void recoilStudy::findbestB()
{
    double tmpintpur=-999;
    indexbestB=-999;
    bestB0=1;
    
    for (int iB0=0; iB0<nB0; iB0++){
      int mode = modeB0[iB0]-10000;
      if(brecointpur[mode]>tmpintpur) {
	indexbestB = iB0;    
        tmpintpur = brecointpur[mode];
//	tmpintpur = intpurB0[iB0];   // old 
      }
    }
    
    for (int iChB=0; iChB<nChB; iChB++){
      int mode = modeChB[iChB]-10000;
      if(brecointpur[mode]>tmpintpur) {
	indexbestB = iChB;    
        tmpintpur = brecointpur[mode];
//	tmpintpur = intpurChB[iChB];  // old
	bestB0=0;
      }
    }

}


void recoilStudy::bookHist(){
  char name[100],title[100];
  fHistFile->mkdir("sg","Signal");
  fHistFile->mkdir("bg","Sideband"); 
  for (int i = 0; i < 600; ++i) {
    int thebmode = i%100;
    if(thebmode<54) {
      fHistFile->cd();  
      sprintf(name, "h%d", 11000+i);  sprintf(title, "mes, mode = %d", 11000+i);  h = new TH1D(name, title, 40, 5.2, 5.3); 
      sprintf(name, "h%d", 12000+i);  sprintf(title, "mes, mode = %d", 12000+i);  h = new TH1D(name, title, 40, 5.2, 5.3); 
      sprintf(name, "h%d", 13000+i);  sprintf(title, "mes, mode = %d", 13000+i);  h = new TH1D(name, title, 40, 5.2, 5.3); 
      sprintf(name, "h%d", 14000+i);  sprintf(title, "mes, mode = %d", 14000+i);  h = new TH1D(name, title, 40, 5.2, 5.3); 
      sprintf(name, "h%d", 15000+i);  sprintf(title, "mes, mode = %d", 15000+i);  h = new TH1D(name, title, 40, 5.2, 5.3); 
      fHistFile->cd("sg");
      sprintf(name, "mx%d", 11000+i);  sprintf(title, "mx recoil, mode = %d", 11000+i);  h = new TH1D(name, title, 40, 0., 4.); 
      sprintf(name, "mx%d", 12000+i);  sprintf(title, "mx recoil, mode = %d", 12000+i);  h = new TH1D(name, title, 40, 0., 4.); 
      sprintf(name, "mx%d", 13000+i);  sprintf(title, "mx recoil, mode = %d", 13000+i);  h = new TH1D(name, title, 40, 0., 4.); 
      sprintf(name, "mx%d", 14000+i);  sprintf(title, "mx recoil, mode = %d", 14000+i);  h = new TH1D(name, title, 40, 0., 4.); 
      sprintf(name, "mx%d", 15000+i);  sprintf(title, "mx recoil, mode = %d", 15000+i);  h = new TH1D(name, title, 40, 0., 4.); 
      fHistFile->cd("bg");
      sprintf(name, "mx%d", 11000+i);  sprintf(title, "mx recoil, mode = %d", 11000+i);  h = new TH1D(name, title, 40, 0., 4.); 
      sprintf(name, "mx%d", 12000+i);  sprintf(title, "mx recoil, mode = %d", 12000+i);  h = new TH1D(name, title, 40, 0., 4.); 
      sprintf(name, "mx%d", 13000+i);  sprintf(title, "mx recoil, mode = %d", 13000+i);  h = new TH1D(name, title, 40, 0., 4.); 
      sprintf(name, "mx%d", 14000+i);  sprintf(title, "mx recoil, mode = %d", 14000+i);  h = new TH1D(name, title, 40, 0., 4.); 
      sprintf(name, "mx%d", 15000+i);  sprintf(title, "mx recoil, mode = %d", 15000+i);  h = new TH1D(name, title, 40, 0., 4.); 
       
      fHistFile->cd("sg");
      sprintf(name, "mxfit%d", 11000+i);  sprintf(title, "mx recoil, mode = %d", 11000+i);  h = new TH1D(name, title, 100, 0., 4.); 
      sprintf(name, "mxfit%d", 12000+i);  sprintf(title, "mx recoil, mode = %d", 12000+i);  h = new TH1D(name, title, 100, 0., 4.); 
      sprintf(name, "mxfit%d", 13000+i);  sprintf(title, "mx recoil, mode = %d", 13000+i);  h = new TH1D(name, title, 100, 0., 4.); 
      sprintf(name, "mxfit%d", 14000+i);  sprintf(title, "mx recoil, mode = %d", 14000+i);  h = new TH1D(name, title, 100, 0., 4.); 
      sprintf(name, "mxfit%d", 15000+i);  sprintf(title, "mx recoil, mode = %d", 15000+i);  h = new TH1D(name, title, 100, 0., 4.); 
      fHistFile->cd("bg");
      sprintf(name, "mxfit%d", 11000+i);  sprintf(title, "mx recoil, mode = %d", 11000+i);  h = new TH1D(name, title, 100, 0., 4.); 
      sprintf(name, "mxfit%d", 12000+i);  sprintf(title, "mx recoil, mode = %d", 12000+i);  h = new TH1D(name, title, 100, 0., 4.); 
      sprintf(name, "mxfit%d", 13000+i);  sprintf(title, "mx recoil, mode = %d", 13000+i);  h = new TH1D(name, title, 100, 0., 4.); 
      sprintf(name, "mxfit%d", 14000+i);  sprintf(title, "mx recoil, mode = %d", 14000+i);  h = new TH1D(name, title, 100, 0., 4.); 
      sprintf(name, "mxfit%d", 15000+i);  sprintf(title, "mx recoil, mode = %d", 15000+i);  h = new TH1D(name, title, 100, 0., 4.); 
    }

  }
  

}

void allbgsub(const char *dir, int opt, int file = 2) {
  char name[100],name2[100],name3[100],bmode[100],dmode[100], dmode2[100], dmode3[100];
  int min,  max;
  // dstar:  file 2, modes 13000 ... 14000
  if (file == 2) {
    min = 13001.; 
    max = 13400.;
    sprintf(dmode2, "DstarX"); 
    sprintf(dmode3, "%s%s","B0->",dmode2); 
  }

  // dc:  file 4, modes 12000 ... 13000
  if (file == 4) {
    min = 12001.; 
    max = 12500.;
    sprintf(dmode2, "DcX"); 
    sprintf(dmode3, "%s%s","B0->",dmode2); 
  }

  // dstar0:  file 3, modes 14000 ... 15000
  if (file == 3) {
    min = 14001.; 
    max = 15400.;
    sprintf(dmode2, "Dstar0X"); 
    sprintf(dmode3, "%s%s","B+->",dmode2); 
   }

  // d0:  file 5, modes 11000 ... 12000
  if (file == 5) {
    min = 11001.; 
    max = 11400.;
    sprintf(dmode2, "D0X"); 
    sprintf(dmode3, "%s%s","B+->",dmode2); 
  }

  
  sprintf(name, "%s%s%d%s", dir,"/recoil-",file,".html");
  ofstream OUT(name);
  OUT << "<!DOCTYPE HTML PUBLIC \"-//IETF//DTD HTML//EN\">" << endl;
  OUT << "<html>" << endl;
  OUT << "<head>" << endl;
  OUT << "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\">" << endl;
  OUT << "<title>All B modes: </title>" << endl;
  OUT << "</head>" << endl;
  OUT << "<body bgcolor=white>" << endl;
  OUT << "<br>" << endl;
  OUT << "<br>" << endl;
  OUT << "<center><b><font size=+4> SemiExclusive Reco, X system study: " << dmode3 << " </font></b></center>" << endl;
  OUT << "<br>" << endl;
  OUT << "<br>" << endl;
  OUT << "<br>" << endl;
  OUT << "<br>" << endl;
  OUT << "<caption></caption>" << endl;
  OUT << "<table border=1  align=\"center\">" << endl;
  OUT << "<tr>" << endl;
  OUT << "<td><b><font color=red> <font size=+2>B mode</font></font></b></td>" 
      << "<td><b><font color=red> <font size=+2>D mode </font></font></b></td>"
      << "<td><b><font color=blue> <font size=+2> All plots  </font></font></b></td>"
      << "<td><b><font color=blue> <font size=+2> Signal & Sidebands  </font></font></b></td>"
      << "<td><b><font color=blue> <font size=+2> Sideband subtracted</font></font></b></td>"
      << endl;
  OUT << "</tr>" << endl;

  for (int mode = min; mode < max; ++mode) {
    
    if(file == 4){
      if(mode%100>53) continue;
    }else{
      if(mode%100>53) continue;
    }
    if(file == 3) {
      if(mode>14353&&mode<15001) continue;
    }

    bgsub(mode, opt, dir);
    Dmode(mode,dmode);
    Bmode(mode,bmode);
    
    sprintf(name,"%d%s",mode,".eps.gz");
    sprintf(name2,"%s%d%s","sigside",mode,".eps.gz");
    sprintf(name3,"%s%d%s","sidesub",mode,".eps.gz");

    OUT << "<tr><td>" << bmode << "</td>" 
	<< "<td>" << dmode << "</td>"
	<< "<td><a href=\"" << name << "\">" << "allplots" << "</td>"
	<< "<td><a href=\"" << name2 << "\">"  << "plot" << "</td>"
	<< "<td><a href=\"" << name3 << "\">"  << "plot" << "</td>"
	<< "</tr>" << endl;
    
    }

  OUT << "</table>" << endl;
  OUT << "</body>" << endl;
  OUT << "</html>" << endl;
  

  //allD modes added up
    
  sprintf(name, "%s%s%d%s", dir,"/recoilallmodes-",file,".html");
  ofstream OUT(name);
  OUT << "<!DOCTYPE HTML PUBLIC \"-//IETF//DTD HTML//EN\">" << endl;
  OUT << "<html>" << endl;
  OUT << "<head>" << endl;
  OUT << "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\">" << endl;
  OUT << "<title>All B modes: </title>" << endl;
  OUT << "</head>" << endl;
  OUT << "<body bgcolor=white>" << endl;
  OUT << "<br>" << endl;
  OUT << "<br>" << endl;
  OUT << "<center><b><font size=+4> SemiExclusive Reco, X system study, D modes added up: " << dmode3 << " </font></b></center>" << endl;
  OUT << "<br>" << endl;
  OUT << "<br>" << endl;
  OUT << "<br>" << endl;
  OUT << "<br>" << endl;
  OUT << "<caption></caption>" << endl;
  OUT << "<table border=1  align=\"center\">" << endl;
  OUT << "<tr>" << endl;
  OUT << "<td><b><font color=red> <font size=+2>B mode</font></font></b></td>" 
      << "<td><b><font color=blue> <font size=+2> All plots  </font></font></b></td>"
      << "<td><b><font color=blue> <font size=+2> Signal & Sidebands  </font></font></b></td>"
      << "<td><b><font color=blue> <font size=+2> Sideband subtracted</font></font></b></td>"
      << endl;
  OUT << "</tr>" << endl;

  for (int mode = min; mode < min+53; ++mode) {
    
    cout << mode << endl;
    if(file == 4){
      if(mode%100>53) continue;
    }else{
      if(mode%100>53) continue;
    }
 
    bgsub(mode, opt, dir, 1);
    Bmode(mode,bmode);
    
    sprintf(name,"%s%d%s","add",mode,".eps.gz");
    sprintf(name2,"%s%d%s","sigsideadd",mode,".eps.gz");
    sprintf(name3,"%s%d%s","sidesubadd",mode,".eps.gz");

    OUT << "<tr><td>" << bmode << "</td>" 
	<< "<td><a href=\"" << name << "\">" << "allplots" << "</td>"
	<< "<td><a href=\"" << name2 << "\">"  << "plot" << "</td>"
	<< "<td><a href=\"" << name3 << "\">"  << "plot" << "</td>"
	<< "</tr>" << endl;
    
  }
  
  OUT << "</table>" << endl;
  OUT << "</body>" << endl;
  OUT << "</html>" << endl;
  
  if (file == 2) {

    //full spectrum  added up for Dstar
    
    sprintf(name, "%s%s%d%s", dir,"/recoilallspec-",file,".html");
    ofstream OUT(name);
    OUT << "<!DOCTYPE HTML PUBLIC \"-//IETF//DTD HTML//EN\">" << endl;
    OUT << "<html>" << endl;
    OUT << "<head>" << endl;
    OUT << "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\">" << endl;
    OUT << "<title>All B modes: </title>" << endl;
    OUT << "</head>" << endl;
    OUT << "<body bgcolor=white>" << endl;
    OUT << "<br>" << endl;
    OUT << "<br>" << endl;
    OUT << "<center><b><font size=+4> SemiExclusive Reco, X system study, D modes added up full spectrum: "  << dmode3 << " </font></b></center>" << endl;
    OUT << "<br>" << endl;
    OUT << "<br>" << endl;
    OUT << "<br>" << endl;
    OUT << "<br>" << endl;
    OUT << "<caption></caption>" << endl;
    OUT << "<table border=1  align=\"center\">" << endl;
    OUT << "<tr>" << endl;
    OUT << "<td><b><font color=red> <font size=+2>B mode</font></font></b></td>" 
	<< "<td><b><font color=blue> <font size=+2> All plots  </font></font></b></td>"
	<< "<td><b><font color=blue> <font size=+2> Signal & Sidebands  </font></font></b></td>"
	<< "<td><b><font color=blue> <font size=+2> Sideband subtracted</font></font></b></td>"
	<< endl;
    OUT << "</tr>" << endl;
    
    for (int mode = min; mode < min+31; ++mode) {
      

      cout << mode << endl;
      
      bgsub(mode, opt, dir, 2);
      Bmode(mode,bmode);
      
      sprintf(name,"%s%d%s","multi",mode,".eps.gz");
      sprintf(name2,"%s%d%s","sigsidemulti",mode,".eps.gz");
      sprintf(name3,"%s%d%s","sidesubmulti",mode,".eps.gz");
      
      OUT << "<tr><td>" << bmode << "</td>" 
	  << "<td><a href=\"" << name << "\">" << "allplots" << "</td>"
	  << "<td><a href=\"" << name2 << "\">"  << "plot" << "</td>"
	  << "<td><a href=\"" << name3 << "\">"  << "plot" << "</td>"
	  << "</tr>" << endl;
      
    }
    
    OUT << "</table>" << endl;
    OUT << "</body>" << endl;
    OUT << "</html>" << endl;
    
  }
  

}
  
  

void bgsub(int histo, int opt, const char *dir, bool added = 0) {
  
  char name[100],name2[100],name3[100],name4[100], nameps[100], hist[100],command[1000];
  recoilAnalysis k;
  sprintf(hist,"%d",histo);
  sprintf(name,"%s%s","mx",hist);
  if (opt) sprintf(name,"%s%s","mxfit",hist);
  sprintf(name2,"%s%s","h",hist);
  if (added==1) {
    sprintf(name,"%s%s","mxtot",hist);
    if (opt) sprintf(name,"%s%s","mxfittot",hist);
    sprintf(name2,"%s%s","htot",hist);
  }
  if (added==2) {
    sprintf(name,"%s%s","mxall",hist);
    if (opt) sprintf(name,"%s%s","mxfitall",hist);
    sprintf(name2,"%s%s","hall",hist);
  }
  c0.Clear();
  c0.Divide(2,2);
  c0.cd(1);
  ((TH1D*)gDirectory->Get(name2))->SetStats(0);
  cout << name << " " << name2 << endl;
  TH1D *sub = k.bgSubtracted(name,"",name2,1);
  ((TH1D*)gDirectory->Get(name2))->Draw();
  c0.cd(2);
  double bgscale = k.getBgScale(((TH1D*)gDirectory->Get(name2)));
  sprintf(name3,"%s%s","bg/",name);
  sprintf(name4,"%s%s","sg/",name);
  ((TH1D*)gDirectory->Get(name3))->SetStats(0);
  ((TH1D*)gDirectory->Get(name3))->SetXTitle("MX (GeV)");
  TH1D *sig((TH1D*)gDirectory->Get(name4));  
  TH1D *side((TH1D*)gDirectory->Get(name3));  
  side->Add(side,side,bgscale,0.);
  sig->SetXTitle("MX (GeV)");
  sig->Draw();
  side->SetFillStyle(3005);
  side->SetLineColor(kBlue);
  side->SetFillColor(kBlue);  
  side->SetXTitle("MX (GeV)");
  side->DrawCopy("same");
  c0.cd(3);
  sub->SetStats(0);
  sub->SetXTitle("MX (GeV)");
  sub->Draw("pe");
  sprintf(nameps,"%s%s%d%s",dir,"/",histo,".eps");
  if (added==1) sprintf(nameps,"%s%s%d%s",dir,"/add",histo,".eps");
  if (added==2) sprintf(nameps,"%s%s%d%s",dir,"/multi",histo,".eps");
  c0.SaveAs(nameps);
  sprintf(command, "gzip -f %s", nameps);
  cout << command << endl;
  gSystem->Exec(command);
  c0.Clear();
  c0.Divide(1);
  sig->SetXTitle("MX (GeV)");
  sig->Draw();
  side->SetXTitle("MX (GeV)");
  side->DrawCopy("same");
  sprintf(nameps,"%s%s%d%s",dir,"/sigside",histo,".eps");
  if (added==1) sprintf(nameps,"%s%s%d%s",dir,"/sigsideadd",histo,".eps");
  if (added==2) sprintf(nameps,"%s%s%d%s",dir,"/sigsidemulti",histo,".eps");
  c0.SaveAs(nameps);
  sprintf(command, "gzip -f %s", nameps);
  cout << command << endl;
  gSystem->Exec(command);
  c0.Clear();
  c0.Divide(1);
  sub->Draw("pe");
  sprintf(nameps,"%s%s%d%s",dir,"/sidesub",histo,".eps");
  if (added==1) sprintf(nameps,"%s%s%d%s",dir,"/sidesubadd",histo,".eps");
  if (added==2) sprintf(nameps,"%s%s%d%s",dir,"/sidesubmulti",histo,".eps");
  c0.SaveAs(nameps);  
  sprintf(command, "gzip -f %s", nameps);
  cout << command << endl;
  gSystem->Exec(command);

}
  
void  Bmode(int theid, char filename[100]){
  
  char filename[100];

  int bmode = theid%100;  
  
  if(bmode == 1)  sprintf(filename,"%s%s","B->D","pi");
  if(bmode == 2)  sprintf(filename,"%s%s","B->D","k");
  if(bmode == 3)  sprintf(filename,"%s%s","B->D","pipi0_<1.5GeV");
  if(bmode == 4)  sprintf(filename,"%s%s","B->D","kpi0_<1.5GeV");
  if(bmode == 5)  sprintf(filename,"%s%s","B->D","piks");
  if(bmode == 6)  sprintf(filename,"%s%s","B->D","kks");
  if(bmode == 7)  sprintf(filename,"%s%s","B->D","pi2pi0_<1.5GeV");
  if(bmode == 8)  sprintf(filename,"%s%s","B->D","k2pi0_<1.5GeV");
  if(bmode == 9)  sprintf(filename,"%s%s","B->D","3pi_<1.5GeV");
  if(bmode == 10)  sprintf(filename,"%s%s","B->D","k2pi_<1.5GeV");
  if(bmode == 11)  sprintf(filename,"%s%s","B->D","2kpi_Ds");
  if(bmode == 12)  sprintf(filename,"%s%s","B->D","omegah");
  if(bmode == 13)  sprintf(filename,"%s%s","B->D","k2pipi0_<2.2GeV");
  if(bmode == 14)  sprintf(filename,"%s%s","B->D","2kpipi0_Ds*");
  if(bmode == 15)  sprintf(filename,"%s%s","B->D","pipi0ks");
  if(bmode == 16)  sprintf(filename,"%s%s","B->D","kpi0ks_<1.8GeV");
  if(bmode == 17)  sprintf(filename,"%s%s","B->D","k2pi0ks_1.8-2.2GeV");
  if(bmode == 18)  sprintf(filename,"%s%s","B->D","2ksX");
  if(bmode == 19)  sprintf(filename,"%s%s","B->D","3pi2pi0_<2.2GeV");
  if(bmode == 20)  sprintf(filename,"%s%s","B->D","k2pi2pi0_<2.2GeV");
  if(bmode == 21)  sprintf(filename,"%s%s","B->D","2kpi2pi0_Ds*");
  if(bmode == 22)  sprintf(filename,"%s%s","B->D","5pi_<2.3GeV");
  if(bmode == 23)  sprintf(filename,"%s%s","B->D","k4p_<2.7GeV");
  if(bmode == 24)  sprintf(filename,"%s%s","B->D","2K3pi_<2.7GeV");
  if(bmode == 25)  sprintf(filename,"%s%s","B->D","5pipi0_<2.2GeV");
  if(bmode == 26)  sprintf(filename,"%s%s","B->D","k4pipi0_<2.2GeV");
  if(bmode == 27)  sprintf(filename,"%s%s","B->D","2k3pipi0_<2.5GeV");
  if(bmode == 28)  sprintf(filename,"%s%s","B->D","3piks_D*");
  if(bmode == 29)  sprintf(filename,"%s%s","B->D","3pikspi0_D*");
  if(bmode == 30)  sprintf(filename,"%s%s","B->D","k2piks_D*");
  if(bmode == 31)  sprintf(filename,"%s%s","B->D","D*_Dpi0");
  if(bmode == 32)  sprintf(filename,"%s%s","B->D","pipi0_>1.5GeV");
  if(bmode == 33)  sprintf(filename,"%s%s","B->D","kpi0_>1.5GeV");
  if(bmode == 34)  sprintf(filename,"%s%s","B->D","pi2pi0_1.5-2GeV");
  if(bmode == 35)  sprintf(filename,"%s%s","B->D","k2pi0_>1.5GeV");
  if(bmode == 36)  sprintf(filename,"%s%s","B->D","3pi_1.5-2GeV");
  if(bmode == 37)  sprintf(filename,"%s%s","B->D","k2pi_>1.5GeV");
  if(bmode == 38)  sprintf(filename,"%s%s","B->D","2kpi_K*");
  if(bmode == 39)  sprintf(filename,"%s%s","B->D","2kpi_other");
  if(bmode == 40)  sprintf(filename,"%s%s","B->D","3pipi0_<1.6GeV");
  if(bmode == 41)  sprintf(filename,"%s%s","B->D","3pipi0_1.6-2.2GeV");
  if(bmode == 42)  sprintf(filename,"%s%s","B->D","k2pipi0_>2.2GeV");
  if(bmode == 43)  sprintf(filename,"%s%s","B->D","2kpipi0_other");
  if(bmode == 44)  sprintf(filename,"%s%s","B->D","kpi0ks_>1.8GeV");
  if(bmode == 45)  sprintf(filename,"%s%s","B->D","3pi2pi0_>2.2GeV");
  if(bmode == 46)  sprintf(filename,"%s%s","B->D","k2pi2pi0_>2.2GeV");
  if(bmode == 47)  sprintf(filename,"%s%s","B->D","2kpi2pi0_other");
  if(bmode == 48)  sprintf(filename,"%s%s","B->D","5pi_>2.3GeV");
  if(bmode == 49)  sprintf(filename,"%s%s","B->D","k4p_>2.7GeV");
  if(bmode == 50)  sprintf(filename,"%s%s","B->D","2K3pi_>2.7GeV");
  if(bmode == 51)  sprintf(filename,"%s%s","B->D","5pipi0_>2.2GeV");
  if(bmode == 52)  sprintf(filename,"%s%s","B->D","3piks_noD*");
  if(bmode == 53)  sprintf(filename,"%s%s","B->D","3pikspi0_noD*");

  return;

}

void Dmode(int theid, char filename[100]){
  
  char filename[100];
  char seed[100];

  int dmode = theid/100;
  
  if(dmode == 110)   sprintf(filename,"D0->kpi");
  if(dmode == 111)   sprintf(filename,"D0->kpipi0");
  if(dmode == 112)   sprintf(filename,"D0->k3pi");
  if(dmode == 113)   sprintf(filename,"D0->kspipi");
  if(dmode == 130)   sprintf(filename,"D*,D0->kpi");
  if(dmode == 131)   sprintf(filename,"D*,D0->kpipi0");
  if(dmode == 132)   sprintf(filename,"D*,D0->k3pi");
  if(dmode == 133)   sprintf(filename,"D*,D0->kspipi");
  if(dmode == 120)   sprintf(filename,"Dc->kspi");
  if(dmode == 121)   sprintf(filename,"Dc->kpipi");
  if(dmode == 122)   sprintf(filename,"Dc->kspipi0");
  if(dmode == 123)   sprintf(filename,"Dc->kpipipi0");
  if(dmode == 124)   sprintf(filename,"Dc->kspipipi");
  if(dmode == 140)   sprintf(filename,"D*0->D0pi0,D0->kpi");
  if(dmode == 141)   sprintf(filename,"D*0->D0pi0,D0->kpipi0");
  if(dmode == 142)   sprintf(filename,"D*0->D0pi0,D0->k3pi");
  if(dmode == 143)   sprintf(filename,"D*0->D0pi0,D0->kspipi");
  if(dmode == 150)   sprintf(filename,"D*0->D0gamma,D0->kpi");
  if(dmode == 151)   sprintf(filename,"D*0->D0gamma,D0->kpipi0");
  if(dmode == 152)   sprintf(filename,"D*0->D0gamma,D0->k3pi");
  if(dmode == 153)   sprintf(filename,"D*0->D0gamma,D0->kspipi");

  return;

}

void sum(const char *file){
  char name[100], name2[100], title[100];

  TFile *fHistFile2;
  fHistFile2 = new TFile(file, "UPDATE");
  fHistFile2->cd();
  for (int i = 0; i < 100; ++i) {
    int thebmode = i;
    if(thebmode<54) {
      fHistFile2->cd();  
      sprintf(name, "htot%d", 11000+i);  sprintf(title, "mes, mode = %d", 11000+i);  h = new TH1D(name, title, 40, 5.2, 5.3); 
      sprintf(name, "htot%d", 12000+i);  sprintf(title, "mes, mode = %d", 12000+i);  h = new TH1D(name, title, 40, 5.2, 5.3); 
      sprintf(name, "htot%d", 13000+i);  sprintf(title, "mes, mode = %d", 13000+i);  h = new TH1D(name, title, 40, 5.2, 5.3); 
      sprintf(name, "htot%d", 14000+i);  sprintf(title, "mes, mode = %d", 14000+i);  h = new TH1D(name, title, 40, 5.2, 5.3); 
      fHistFile2->cd("sg");
      sprintf(name, "mxtot%d", 11000+i);  sprintf(title, "mx recoil, mode = %d", 11000+i);  h = new TH1D(name, title, 40, 0., 4.); 
      sprintf(name, "mxtot%d", 12000+i);  sprintf(title, "mx recoil, mode = %d", 12000+i);  h = new TH1D(name, title, 40, 0., 4.); 
      sprintf(name, "mxtot%d", 13000+i);  sprintf(title, "mx recoil, mode = %d", 13000+i);  h = new TH1D(name, title, 40, 0., 4.); 
      sprintf(name, "mxtot%d", 14000+i);  sprintf(title, "mx recoil, mode = %d", 14000+i);  h = new TH1D(name, title, 40, 0., 4.); 
      fHistFile2->cd("bg");
      sprintf(name, "mxtot%d", 11000+i);  sprintf(title, "mx recoil, mode = %d", 11000+i);  h = new TH1D(name, title, 40, 0., 4.); 
      sprintf(name, "mxtot%d", 12000+i);  sprintf(title, "mx recoil, mode = %d", 12000+i);  h = new TH1D(name, title, 40, 0., 4.); 
      sprintf(name, "mxtot%d", 13000+i);  sprintf(title, "mx recoil, mode = %d", 13000+i);  h = new TH1D(name, title, 40, 0., 4.); 
      sprintf(name, "mxtot%d", 14000+i);  sprintf(title, "mx recoil, mode = %d", 14000+i);  h = new TH1D(name, title, 40, 0., 4.); 
      fHistFile2->cd("sg");
      sprintf(name, "mxfittot%d", 11000+i);  sprintf(title, "mx recoil, mode = %d", 11000+i);  h = new TH1D(name, title, 100, 0., 4.); 
      sprintf(name, "mxfittot%d", 12000+i);  sprintf(title, "mx recoil, mode = %d", 12000+i);  h = new TH1D(name, title, 100, 0., 4.); 
      sprintf(name, "mxfittot%d", 13000+i);  sprintf(title, "mx recoil, mode = %d", 13000+i);  h = new TH1D(name, title, 100, 0., 4.); 
      sprintf(name, "mxfittot%d", 14000+i);  sprintf(title, "mx recoil, mode = %d", 14000+i);  h = new TH1D(name, title, 100, 0., 4.); 
      fHistFile2->cd("bg");
      sprintf(name, "mxfittot%d", 11000+i);  sprintf(title, "mx recoil, mode = %d", 11000+i);  h = new TH1D(name, title, 100, 0., 4.); 
      sprintf(name, "mxfittot%d", 12000+i);  sprintf(title, "mx recoil, mode = %d", 12000+i);  h = new TH1D(name, title, 100, 0., 4.); 
      sprintf(name, "mxfittot%d", 13000+i);  sprintf(title, "mx recoil, mode = %d", 13000+i);  h = new TH1D(name, title, 100, 0., 4.); 
      sprintf(name, "mxfittot%d", 14000+i);  sprintf(title, "mx recoil, mode = %d", 14000+i);  h = new TH1D(name, title, 100, 0., 4.); 

      fHistFile2->cd();
      for (int j=0; j<4; j++){	
	sprintf(name, "htot%d", 11000+i); sprintf(name2, "h%d", 11000+i+j*100); ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get(name2)));
	sprintf(name, "htot%d", 13000+i); sprintf(name2, "h%d", 13000+i+j*100); ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get(name2)));
	sprintf(name, "htot%d", 14000+i); sprintf(name2, "h%d", 14000+i+j*100); ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get(name2)));
	sprintf(name, "htot%d", 14000+i); sprintf(name2, "h%d", 15000+i+j*100); ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get(name2)));
      }
      for (int j=0; j<5; j++){
	sprintf(name, "htot%d", 12000+i); sprintf(name2, "h%d", 12000+i+j*100); ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get(name2)));
      }													      

      fHistFile2->cd("sg");
      for (int j=0; j<4; j++){	
	sprintf(name, "mxtot%d", 11000+i); sprintf(name2, "mx%d", 11000+i+j*100); ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get(name2)));
	sprintf(name, "mxtot%d", 13000+i); sprintf(name2, "mx%d", 13000+i+j*100); ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get(name2)));
	sprintf(name, "mxtot%d", 14000+i); sprintf(name2, "mx%d", 14000+i+j*100); ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get(name2)));
	sprintf(name, "mxtot%d", 14000+i); sprintf(name2, "mx%d", 15000+i+j*100); ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get(name2)));
	sprintf(name, "mxfittot%d", 11000+i); sprintf(name2, "mxfit%d", 11000+i+j*100); ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get(name2)));
	sprintf(name, "mxfittot%d", 13000+i); sprintf(name2, "mxfit%d", 13000+i+j*100); ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get(name2)));
	sprintf(name, "mxfittot%d", 14000+i); sprintf(name2, "mxfit%d", 14000+i+j*100); ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get(name2)));
	sprintf(name, "mxfittot%d", 14000+i); sprintf(name2, "mxfit%d", 15000+i+j*100); ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get(name2)));
      }
      for (int j=0; j<5; j++){
	sprintf(name, "mxtot%d", 12000+i); sprintf(name2, "mx%d", 12000+i+j*100); ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get(name2)));
 	sprintf(name, "mxfittot%d", 12000+i); sprintf(name2, "mxfit%d", 12000+i+j*100); ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get(name2)));
      }													      

      fHistFile2->cd("bg");
      for (int j=0; j<4; j++){	
	sprintf(name, "mxtot%d", 11000+i); sprintf(name2, "mx%d", 11000+i+j*100); ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get(name2)));
	sprintf(name, "mxtot%d", 13000+i); sprintf(name2, "mx%d", 13000+i+j*100); ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get(name2)));
	sprintf(name, "mxtot%d", 14000+i); sprintf(name2, "mx%d", 14000+i+j*100); ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get(name2)));
	sprintf(name, "mxtot%d", 14000+i); sprintf(name2, "mx%d", 15000+i+j*100); ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get(name2)));
	sprintf(name, "mxfittot%d", 11000+i); sprintf(name2, "mxfit%d", 11000+i+j*100); ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get(name2)));
	sprintf(name, "mxfittot%d", 13000+i); sprintf(name2, "mxfit%d", 13000+i+j*100); ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get(name2)));
	sprintf(name, "mxfittot%d", 14000+i); sprintf(name2, "mxfit%d", 14000+i+j*100); ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get(name2)));
	sprintf(name, "mxfittot%d", 14000+i); sprintf(name2, "mxfit%d", 15000+i+j*100); ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get(name2)));
      }
      for (int j=0; j<5; j++){
	sprintf(name, "mxtot%d", 12000+i); sprintf(name2, "mx%d", 12000+i+j*100); ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get(name2)));
 	sprintf(name, "mxfittot%d", 12000+i); sprintf(name2, "mxfit%d", 12000+i+j*100); ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get(name2)));
      }													      
      

    }
    
  }
  
  fHistFile2->Write();
  fHistFile2->Close();

    
}

void sumregions(const char *file){
  char name[100], name2[100], title[100];
  double coeff[32];
  int    hist[32];

  coeff[0] = 0.;  //B->Dpi
  coeff[1] = 0.;  //B->Dk
  coeff[2] = 3.;  //B->Dpipi0 (Add mode 32)
  coeff[3] = 1.8;  //B->Dkpi0 (add mode 33)
  coeff[4] = 0.;  //B->Dpiks
  coeff[5] = 0.;  //B->Dkks
  coeff[6] = 4.;  //B->Dpi2pi0 (Add mode 34)
  coeff[7] = 3.;  //B->Dk2pi0 (add mode 35)
  coeff[8] = 1.5;  //B->D3pi (add mode 36)
  coeff[9] = 2.2;  //B->Dk2pi (add mode 37)
  coeff[10] = 2.5;  //B->D2kpi_Ds (only added by mode 39, no K*)
  coeff[11] = 2.5;  //B->D3pipi0 (no omega pi included) mode 40-41
  coeff[12] = 4.5;  //B->Dk2pipi0
  coeff[13] = 2.5;  //B->D2kpipi0_Ds* (??? don't understand mass window)
  coeff[14] = 0.;  //B->Dpipi0ks
  coeff[15] = 1.;  //B->Dkpi0ks_ (one is fine)
  coeff[16] = 0.;  //B->Dk2pi0ks_1.8-2.2GeV (no breakdown)
  coeff[17] = 0.;  //B->D2ksX
  coeff[18] = 1.5;  //B->D3pi2pi0 ???
  coeff[19] = 30.;  //B->Dk2pi2pi0
  coeff[20] = 6.;  //B->D2kpi2pi0_Ds*  mass window??????
  coeff[21] = 0.3;//B->D5pi_
  coeff[22] = 2.5; //B->Dk4p
  coeff[23] = 3.; //B->D2K3pi
  coeff[24] = 8.;  //B->D5pipi0
  coeff[25] = 0.;  //B->Dk4pipi0_<2.2GeV
  coeff[26] = 0.;  //B->D2k3pipi0_<2.5GeV
  coeff[27] = 2.5;  //B->D3piks_D*
  coeff[28] = 0.;  //B->D3pikspi0_D* high mx is not there
  coeff[29] = 0.;  //B->Dk2piks_D*
  coeff[30] = 0.;  //B->DD*_Dpi0
  coeff[31] = 0.;  //B->Dpipi0_>1.5GeV

  hist[0] = 1;  //B->Dpi
  hist[1] = 2;  //B->Dk
  hist[2] = 32;   //B->Dpipi0 (Add mode 32)
  hist[3] = 33;  //B->Dkpi0 (add mode 33)
  hist[4] = 5;  //B->Dpiks
  hist[5] = 6;  //B->Dkks
  hist[6] = 34;  //B->Dpi2pi0 (Add mode 34)
  hist[7] = 35;  //B->Dk2pi0 (add mode 35)
  hist[8] = 36;  //B->D3pi (add mode 36)
  hist[9] = 37;  //B->Dk2pi (add mode 37)
  hist[10] = 39;  //B->D2kpi_Ds (only added by mode 39, no K*)
  hist[11] = 12;  //B->D3pipi0 (no omega pi included) mode 40-41
  hist[12] = 42;  //B->Dk2pipi0
  hist[13] = 43;  //B->D2kpipi0_Ds* (??? don't understand mass window)
  hist[14] = 14;  //B->Dpipi0ks
  hist[15] = 44;  //B->Dkpi0ks_ (one is fine)
  hist[16] = 17;  //B->Dk2pi0ks_1.8-2.2GeV (no breakdown)
  hist[17] = 18;  //B->D2ksX
  hist[18] = 45;  //B->D3pi2pi0 ???
  hist[19] = 46;  //B->Dk2pi2pi0
  hist[20] = 47;  //B->D2kpi2pi0_Ds*  mass window??????
  hist[21] = 48;//B->D5pi_
  hist[22] = 49; //B->Dk4p
  hist[23] = 50; //B->D2K3pi
  hist[24] = 51;  //B->D5pipi0
  hist[25] = 26;  //B->Dk4pipi0_<2.2GeV
  hist[26] = 27;  //B->D2k3pipi0_<2.5GeV
  hist[27] = 52;  //B->D3piks_D*
  hist[28] = 29;  //B->D3pikspi0_D* high mx is not there
  hist[29] = 30;  //B->Dk2piks_D*
  hist[30] = 31;  //B->DD*_Dpi0
  hist[31] = 32;  //B->Dpipi0_>1.5GeV

  TFile *fHistFile2;  
  fHistFile2 = new TFile(file, "UPDATE");
  fHistFile2->cd();
  for (int i = 1; i < 32; ++i) {
    int thebmode = i;
    if(thebmode<54) {
      c0.Clear();
      c0.Divide(2,3);      
      fHistFile2->cd();  
      sprintf(name, "hall%d", 13000+i);  sprintf(title, "mes, mode = %d", 13000+i);  h = new TH1D(name, title, 40, 5.2, 5.3); 
      sprintf(name, "hall%d", 13000+i); sprintf(name2, "htot%d", 13000+i); ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get(name2)));

      fHistFile2->cd("sg");
      sprintf(name, "mxall%d", 13000+i);  sprintf(title, "mx recoil, mode = %d", 13000+i);  h = new TH1D(name, title, 40, 0., 4.); 
      sprintf(name, "mxall%d", 13000+i); sprintf(name2, "mxtot%d", 13000+i); ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get(name2)));
      sprintf(name, "mxfitall%d", 13000+i);  sprintf(title, "mx recoil, mode = %d", 13000+i);  h = new TH1D(name, title, 100, 0., 4.); 
      sprintf(name, "mxfitall%d", 13000+i); sprintf(name2, "mxfittot%d", 13000+i); ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get(name2)));

      fHistFile2->cd("bg");
      sprintf(name, "mxall%d", 13000+i);  sprintf(title, "mx recoil, mode = %d", 13000+i);  h = new TH1D(name, title, 40, 0., 4.); 
      sprintf(name, "mxall%d", 13000+i); sprintf(name2, "mxtot%d", 13000+i); ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get(name2)));
      sprintf(name, "mxfitall%d", 13000+i);  sprintf(title, "mx recoil, mode = %d", 13000+i);  h = new TH1D(name, title, 100, 0., 4.); 
      sprintf(name, "mxfitall%d", 13000+i); sprintf(name2, "mxfittot%d", 13000+i); ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get(name2)));

      fHistFile2->cd();
      sprintf(name, "hall%d", 13000+i); sprintf(name2, "h%d", 13000+hist[i-1]); ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get(name2)),coeff[i-1]);
      c0.cd(1);
      ((TH1D*)gDirectory->Get(name))->Draw();

      fHistFile2->cd("sg");
      sprintf(name, "mxall%d", 13000+i); sprintf(name2, "mxtot%d", 13000+hist[i-1]); ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get(name2)),coeff[i-1]);

      c0.cd(2);
      ((TH1D*)gDirectory->Get(name))->Draw();

      sprintf(name, "mxfitall%d", 13000+i); sprintf(name2, "mxfittot%d", 13000+hist[i-1],coeff[i-1]); ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get(name2)),coeff[i-1]);
  
      c0.cd(3);
      ((TH1D*)gDirectory->Get(name))->Draw();

      fHistFile2->cd("bg");
      sprintf(name, "mxall%d", 13000+i); sprintf(name2, "mxtot%d", 13000+hist[i-1]); ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get(name2)),coeff[i-1]);
      c0.cd(4);
      ((TH1D*)gDirectory->Get(name))->Draw();

      sprintf(name, "mxfitall%d", 13000+i); sprintf(name2, "mxfittot%d", 13000+hist[i-1],coeff[i-1]); ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get(name2)),coeff[i-1]);
  
      
      c0.cd(5);
      ((TH1D*)gDirectory->Get(name))->Draw();
      
      sprintf(name, "plots%d.eps", i);
      c0.SaveAs(name);

      //B->D3pipi0 (no omega pi included) mode 40-41
      //B->D2kpi (other is missing!!!!)
    }

  }
  fHistFile2->Write();
  fHistFile2->Close();

}
