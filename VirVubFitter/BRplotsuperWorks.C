void BRplotsuperWorks(TString thepath, TString thefile, TString theoutpath, Int_t err=-1, Int_t staterr=0) {
  //modify to rescale Bauer for the bug in acceptance
  //err=1 plot BRBR with error bars (default)
  //err=2 plot absolute value or errors
  //err=3 plot relative errors
  //staterr=0 only data statistics (default)
  //staterr=1 data+MC statistics
  char filetsC[100];
  char name[100], name2[100], nameeps[100], title[100];
  char buffertC[100];
  char canvtitl[100], leghead[100], nameY[100];
  double pbf, error;
  int bin;
  
  //  gROOT->SetStyle("Plain");
  gROOT->SetStyle("BABAR");
  TCanvas  *c1 = new TCanvas("canv","canv",800,600);
  c1->cd(1);
  
  for(int ij=staterr; ij<=staterr; ij++){
    int lines = 0;
    if(err<0){
      //      sprintf(name,"%s%i","BRBR",ij);  sprintf(title, "BRBR"); TH1D* h = new TH1D(name, title, 8, -1+0.1*ij,15+0.1*ij);
      sprintf(name,"%s%i","BRBR",ij);  sprintf(title, "BRBR"); TH1D* h = new TH1D(name, title, 150, -0.1+0.01*ij,15.1+0.01*ij);
    }else{
      sprintf(name,"%s%i","BRBR",ij);  sprintf(title, "BRBR"); TH1D* h = new TH1D(name, title, 15, 0.5+0.1*ij,15.5+0.1*ij);
    }
    if(ij==1){
      sprintf(filetsC,"%s%s%i",thepath.Data(),thefile.Data(),15);
    }else if(ij==2){
      sprintf(filetsC,"%s%s%i",thepath.Data(),thefile.Data(),17);
    }else if(ij==3){
      sprintf(filetsC,"%s%s%i",thepath.Data(),thefile.Data(),186);
    }
    cout << "Open file  " << filetsC << ij << endl;
    ifstream tsC(filetsC);
    sprintf(name,"%s%i","BRBR",ij);
    while (tsC.getline(buffertC, 200, '\n')) {
      sscanf(buffertC, "%lf %lf",&pbf,&error);
      int binn=lines*20;
      if(lines==0) binn=2;
      cout << binn << endl;
      lines++;
      cout << pbf << " " << error << endl;
      if(TMath::Abs(err)==2){
	sprintf(nameY,"Stat error (absolute)");
	((TH1D*)gDirectory->Get(name))->SetBinContent(binn,error);
	cout << "Absoluteerror " << error << endl;
      }else if(TMath::Abs(err)==3){
	sprintf(nameY,"Stat error (relative)");
	if(brbr !=0){
	  ((TH1D*)gDirectory->Get(name))->SetBinContent(binn,error/pbf);
	  cout << "Relativeerror " << error/pbf <<  endl;
	}
      }else if(TMath::Abs(err)==1){
	sprintf(nameY,"#Delta Br(#bar{B} #rightarrow X_{u} l #bar{#nu}) (10^{-3})");
	((TH1D*)gDirectory->Get(name))->SetBinContent(binn,pbf);
	((TH1D*)gDirectory->Get(name))->SetBinError(binn,error);
	cout << "Delta BF " << pbf << endl;
      }
    }
  }
  ((TH1D*)gDirectory->Get(name))->SetMinimum(0.);
  ((TH1D*)gDirectory->Get(name))->SetXTitle("q^{2}_{cut} (GeV^{2}/c^{4})");
  ((TH1D*)gDirectory->Get(name))->SetYTitle(nameY);
  //  ((TH1D*)gDirectory->Get(name))->SetTitleOffset(1.3,"Y");
  //  ((TH1D*)gDirectory->Get(name))->SetTitleOffset(1.3,"X");
  //  ((TH1D*)gDirectory->Get(name))->SetMarkerSize(1);  ((TH1D*)gDirectory->Get(name))->SetLineWidth(2);
  if(TMath::Abs(err)==2){
    ((TH1D*)gDirectory->Get(name))->SetMaximum(0.8);
  }else if(TMath::Abs(err)==1){
    ((TH1D*)gDirectory->Get(name))->SetMaximum(2.5);
    ((TH1D*)gDirectory->Get(name))->SetMinimum(0.000001);
  }else if(TMath::Abs(err)==3){
    ((TH1D*)gDirectory->Get(name))->SetMaximum(2.6);
  }
  
  //plot filled histograms
  gStyle->SetOptStat(0);  gStyle->SetOptTitle(0);
  
  TLegendEntry *legge; 
  TLegend *leg;
  leg = new TLegend(0.8,0.65,0.9,0.75);
  leg->SetBorderSize(0.); leg->SetTextSize(0.07);
  //make text size 0.08 to improve reading 
  leg->SetFillStyle(4000);  leg->SetFillColor(0); 
  
  //  ((TH1D*)gDirectory->Get("BRBR3"))->SetMarkerStyle(22); ((TH1D*)gDirectory->Get("BRBR3"))->SetMarkerColor(3);
  //  ((TH1D*)gDirectory->Get("BRBR3"))->SetLineColor(3);
  
  ((TH1D*)gDirectory->Get("BRBR2"))->SetMarkerStyle(21); ((TH1D*)gDirectory->Get("BRBR2"))->SetMarkerColor(2);
  ((TH1D*)gDirectory->Get("BRBR2"))->SetLineColor(2);
  
  //  ((TH1D*)gDirectory->Get("BRBR1"))->SetMarkerStyle(20); ((TH1D*)gDirectory->Get("BRBR1"))->SetMarkerColor(1);
  //  ((TH1D*)gDirectory->Get("BRBR1"))->SetLineColor(1);
  
  //  ((TH1D*)gDirectory->Get("BRBR3"))->Draw("p");
  //  ((TH1D*)gDirectory->Get("BRBR2"))->Draw("psame");
  //  ((TH1D*)gDirectory->Get("BRBR1"))->Draw("psame");

  //  ((TH1D*)gDirectory->Get("BRBR2"))->Draw("p");
  ((TH1D*)gDirectory->Get(name))->Draw("p");

  //  legge = leg->AddEntry(((TH1D*)gDirectory->Get("BRBR3")),  "Mx < 1.86", "p"); 
  //  legge = leg->AddEntry(((TH1D*)gDirectory->Get("BRBR2")),  "Mx < 1.7", "p"); 
  //  legge = leg->AddEntry(((TH1D*)gDirectory->Get("BRBR1")),  "Mx < 1.5", "p"); 
  legge = leg->SetHeader("(a)");
  leg->Draw();
  //  BABARSmartLabel(0.9,0.9,1,"~1");
  //  BABARSmartLabel(0.9,0.7,1,"a)");

  gPad->Update();
  sprintf(nameeps,"%s%s%s%i%i%s",theoutpath.Data(),thefile.Data(),"_s_",TMath::Abs(err),staterr,".eps");
  c1->Print(nameeps);
  sprintf(nameeps,"%s%s%s%i%i%s",theoutpath.Data(),thefile.Data(),"_s_",TMath::Abs(err),staterr,".gif");
  c1->Print(nameeps);
}
