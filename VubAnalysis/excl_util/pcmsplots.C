void pcmsplots(){

  TFile g("output_0412/csx-allsignal.root");
  
  TH1D* pcmspi = new TH1D("pcmspi","     generated p^{*} for #pi, #pi^{0}",50,0.,2.7);
  TH1D* pcmsgenpi = new TH1D("pcmsgenpi","     generated p^{*} for #pi, #pi^{0}",50,0.,2.7);
  TH1D* pcmsmm2pi = new TH1D("pcmsmm2pi","     generated p^{*} for #pi, #pi^{0}",50,0.,2.7); 

  TH1D* pcmsrho = new TH1D("pcmsrho","     generated p^{*} for #rho^{+}, #rho^{0}",50,0.,2.7);
  TH1D* pcmsgenrho = new TH1D("pcmsgenrho","     generated p^{*} for #rho^{+}, #rho^{0}",50,0.,2.7);
  TH1D* pcmsmm2rho = new TH1D("pcmsmm2rho","     generated p^{*} for #rho^{+}, #rho^{0}",50,0.,2.7); 

  TH1D* pcmsomega = new TH1D("pcmsomega","     generated p^{*} for #omega",50,0.,2.7);
  TH1D* pcmsgenomega = new TH1D("pcmsgenomega","     generated p^{*} for #omega",50,0.,2.7);
  TH1D* pcmsmm2omega = new TH1D("pcmsmm2omega","     generated p^{*} for #omega",50,0.,2.7); 

  TH1D* pcmseta = new TH1D("pcmseta","     generated p^{*} for #eta",50,0.,2.7);
  TH1D* pcmsgeneta = new TH1D("pcmsgeneta","     generated p^{*} for #eta",50,0.,2.7);
  TH1D* pcmsmm2eta = new TH1D("pcmsmm2eta","     generated p^{*} for #eta",50,0.,2.7); 

  TH1D* pcmsetap = new TH1D("pcmsetap","     generated p^{*} for #eta'",50,0.,2.7);
  TH1D* pcmsgenetap = new TH1D("pcmsgenetap","     generated p^{*} for #eta'",50,0.,2.7);
  TH1D* pcmsmm2etap = new TH1D("pcmsmm2etap","     generated p^{*} for #eta'",50,0.,2.7); 

  TH1D* pcmsa0 = new TH1D("pcmsa0","     generated p^{*} for a^{+}_{0}, a^{0}_{0}",20,0.,2.7);
  TH1D* pcmsgena0 = new TH1D("pcmsgena0","     generated p^{*} for a^{+}_{0}, a^{0}_{0}",20,0.,2.7);
  TH1D* pcmsmm2a0 = new TH1D("pcmsmm2a0","     generated p^{*} for a^{+}_{0}, a^{0}_{0}",20,0.,2.7); 

  char cuts[1000];

  events.Project("pcmsgenpi","pcmsgen","abs(Gvxbtyp)==11&&mes>5.27");
  events.Project("pcmspi","pcmsgen","abs(Gvxbtyp)==11&&mes>5.27&&pcms>0&&((tcmsgen-tcms)*(tcmsgen-tcms)+(fcmsgen-fcms)*(fcmsgen-fcms)<.01)&&mes>5.27");
  sprintf(cuts,"abs(Gvxbtyp)==11&&mes>5.27&&pcms>0&&((tcmsgen-tcms)*(tcmsgen-tcms)+(fcmsgen-fcms)*(fcmsgen-fcms)<.01)&&mes>5.27&&((abs(phipi0-fxhadgen)<0.05&&abs(thpi0-txhadgen)<0.05&&abs(mompi0-pxhadgen)<0.15&&mm2pi0<.5)||(abs(phipi-fxhadgen)<0.05&&abs(thpi-txhadgen)<0.05&&abs(mompi-pxhadgen)<0.2&&mm2pi<.5))");
  events.Project("pcmsmm2pi","pcmsgen",cuts);

  drawplots(pcmspi,pcmsgenpi,pcmsmm2pi,"pcmspi0.eps");

  events.Project("pcmsgenrho","pcmsgen","abs(Gvxbtyp)==13&&mes>5.27");
  events.Project("pcmsrho","pcmsgen","abs(Gvxbtyp)==13&&mes>5.27&&pcms>0&&((tcmsgen-tcms)*(tcmsgen-tcms)+(fcmsgen-fcms)*(fcmsgen-fcms)<.01)&&mes>5.27");
  sprintf(cuts,"abs(Gvxbtyp)==13&&mes>5.27&&pcms>0&&((tcmsgen-tcms)*(tcmsgen-tcms)+(fcmsgen-fcms)*(fcmsgen-fcms)<.01)&&mes>5.27&&((abs(phirho-fxhadgen)<0.05&&abs(thrho-txhadgen)<0.05&&abs(momrho-pxhadgen)<0.15&&mm2rho<.5)||(abs(phirho0-fxhadgen)<0.05&&abs(thrho0-txhadgen)<0.05&&abs(momrho0-pxhadgen)<0.15&&mm2rho0<.5))");
  events.Project("pcmsmm2rho","pcmsgen",cuts); 
 
  drawplots(pcmsrho,pcmsgenrho,pcmsmm2rho,"pcmsrho.eps");

  events.Project("pcmsgenomega","pcmsgen","abs(Gvxbtyp)==14&&mes>5.27"); 
  events.Project("pcmsomega","pcmsgen","abs(Gvxbtyp)==14&&mes>5.27&&pcms>0&&((tcmsgen-tcms)*(tcmsgen-tcms)+(fcmsgen-fcms)*(fcmsgen-fcms)<.01)&&mes>5.27"); 
  sprintf(cuts,"abs(Gvxbtyp)==14&&mes>5.27&&pcms>0&&((tcmsgen-tcms)*(tcmsgen-tcms)+(fcmsgen-fcms)*(fcmsgen-fcms)<.01)&&mes>5.27&&(abs(phiomega-fxhadgen)<0.05&&abs(thomega-txhadgen)<0.05&&abs(momomega-pxhadgen)<0.2&&mm2omega<.5)");  
  events.Project("pcmsmm2omega","pcmsgen",cuts); 
 
  drawplots(pcmsomega,pcmsgenomega,pcmsmm2omega,"pcmsomega.eps"); 
 
  events.Project("pcmsgeneta","pcmsgen","abs(Gvxbtyp)==12&&mes>5.27"); 
  events.Project("pcmseta","pcmsgen","abs(Gvxbtyp)==12&&mes>5.27&&pcms>0&&((tcmsgen-tcms)*(tcmsgen-tcms)+(fcmsgen-fcms)*(fcmsgen-fcms)<.01)&&mes>5.27"); 
  sprintf(cuts,"abs(Gvxbtyp)==12&&mes>5.27&&pcms>0&&((tcmsgen-tcms)*(tcmsgen-tcms)+(fcmsgen-fcms)*(fcmsgen-fcms)<.01)&&mes>5.27&&(abs(phieta-fxhadgen)<0.05&&abs(theta-txhadgen)<0.05&&abs(mometa-pxhadgen)<0.2&&mm2eta<.5)");
  events.Project("pcmsmm2eta","pcmsgen",cuts); 
 
  drawplots(pcmseta,pcmsgeneta,pcmsmm2eta,"pcmseta.eps"); 

  events.Project("pcmsgenetap","pcmsgen","abs(Gvxbtyp)==15&&mes>5.27");  
  events.Project("pcmsetap","pcmsgen","abs(Gvxbtyp)==15&&mes>5.27&&pcms>0&&((tcmsgen-tcms)*(tcmsgen-tcms)+(fcmsgen-fcms)*(fcmsgen-fcms)<.01)&&mes>5.27");  
  sprintf(cuts,"abs(Gvxbtyp)==15&&mes>5.27&&pcms>0&&((tcmsgen-tcms)*(tcmsgen-tcms)+(fcmsgen-fcms)*(fcmsgen-fcms)<.01)&&mes>5.27&&(abs(phietap-fxhadgen)<0.05&&abs(thetap-txhadgen)<0.05&&abs(mometap-pxhadgen)<0.2&&mm2etap<.5)"); 
  events.Project("pcmsmm2etap","pcmsgen",cuts); 
 
  drawplots(pcmsetap,pcmsgenetap,pcmsmm2etap,"pcmsetap.eps");   

  events.Project("pcmsgena0","pcmsgen","abs(Gvxbtyp)==19&&mes>5.27"); 
  events.Project("pcmsa0","pcmsgen","abs(Gvxbtyp)==19&&mes>5.27&&pcms>0&&((tcmsgen-tcms)*(tcmsgen-tcms)+(fcmsgen-fcms)*(fcmsgen-fcms)<.01)&&mes>5.27"); 
  sprintf(cuts,"abs(Gvxbtyp)==19&&mes>5.27&&pcms>0&&((tcmsgen-tcms)*(tcmsgen-tcms)+(fcmsgen-fcms)*(fcmsgen-fcms)<.01)&&mes>5.27&&((abs(phia0-fxhadgen)<0.05&&abs(tha0-txhadgen)<0.05&&abs(moma0-pxhadgen)<0.2&&mm2a0<.5)||(abs(phia0p-fxhadgen)<0.05&&abs(tha0p-txhadgen)<0.05&&abs(moma0p-pxhadgen)<0.2&&mm2a0p<.5))"); 
  events.Project("pcmsmm2a0","pcmsgen",cuts); 
 
  drawplots(pcmsa0,pcmsgena0,pcmsmm2a0,"pcmsa0.eps"); 
 


}

void drawplots(TH1D* historec, TH1D* histogen, TH1D* histomm2, char* outputname){

  c0.Clear(); c0.cd();
  TPad *fPads[2];
  
  fPads[0]= new TPad("pad1", "", 0.00, 0.35, 0.99, 0.99);   fPads[0]->Draw(); 
  fPads[1]= new TPad("pad2", "", 0.00, 0.00, 0.99, 0.33);   fPads[1]->Draw();  
  
  fPads[0]->cd(); 
  shrinkPad(0.001, 0.1); 

  
  histogen->SetLineColor(kBlack);
  histogen->SetFillColor(kYellow);
  histogen->SetLineWidth(1.5);
  historec->SetLineColor(kBlack);
  historec->SetFillColor(kBlue);
  historec->SetLineWidth(1.5);
  histomm2->SetLineColor(kBlack); 
  histomm2->SetFillColor(kRed); 
  histomm2->SetLineWidth(1.5); 

  histogen->SetStats(0);
  TLegendEntry *legge;  
  TLegend *leg; 
  leg = new TLegend(0.17,0.6,0.40,0.79); 
  leg->SetFillStyle(0); leg->SetBorderSize(0.); leg->SetTextSize(0.05);  
  leg->SetFillColor(0);  
  legge = leg->AddEntry(histogen, "no cut", "f");  
  legge = leg->AddEntry(historec, "rec. lepton", "f");  
  legge = leg->AddEntry(histomm2, "allcuts", "f");   

  histogen->Draw();
  historec->Draw("same");
  histomm2->Draw("same"); 
  leg->Draw(); 

  TH1D effi("effi","#epsilon(lepton reco)",20,0.,2.7);
  TH1D effimm2("effimm2","#epsilon(lept, reso, mm2)",20,0.,2.7); 

  for(int i=1;i<21;i++){

    double eff(0),erreff(0),effmm2(0),erreffmm2(0);
    if(histogen->GetBinContent(i)>0){

      eff = historec->GetBinContent(i)/histogen->GetBinContent(i);
      erreff = sqrt(eff*(1-eff)/histogen->GetBinContent(i));             
      effi->SetBinContent(i,eff);
      effi->SetBinError(i,erreff);

      effmm2 = histomm2->GetBinContent(i)/histogen->GetBinContent(i); 
      erreffmm2 = sqrt(effmm2*(1-effmm2)/histogen->GetBinContent(i));              
      effimm2->SetBinContent(i,effmm2); 
      effimm2->SetBinError(i,erreffmm2); 

    }
  }

  

  fPads[1]->cd(); 
  //shrinkPad(0.001, 0.2); 
  shrinkPad(0.4, 0.1, 0.1, 0.001);

  gPad->SetGridx(1);  gPad->SetGridy(1);
  
  effi.SetStats(0);
  effi.SetMinimum(0.); 
  effi.SetMaximum(0.75);  
  effi.SetXTitle("true p^{*}(GeV)");
  effi.SetYTitle("efficiency"); 
  effi.SetTitle("");  
  effi.SetLabelSize(.08,"X");
  effi.SetLabelSize(.08,"Y"); 
  effi.SetTitleOffset(.9,"X");
  effi.SetTitleOffset(.6,"Y"); 
  effi.SetTitleSize(.1,"X"); 
  effi.SetTitleSize(.07,"Y");  
  effi.SetMarkerColor(kBlue);

  effimm2.SetMarkerColor(kRed); 
 
  effi.Draw("pe");
  effimm2.Draw("samepe"); 

  
  c0.SaveAs(outputname);
  
}
