gSystem.Load("libPhysics.so");
gSystem.Load("../RecoilAnalysis/libRecoilAnalysis.so");

void history(int bins, TString dir, TString dirin, int islumi = 1){
  float xvalues[20], errxvalues[20];
  float lumi[20], errlumi[20];
  float meslepyvalues[20], mesleperryvalues[20];
  float mesallcutsyvalues[20], mesallcutserryvalues[20];
  float mesDpiyvalues[20], mesDpierryvalues[20];
  
  TFile *file[20];
  char name[100];
  char namein[100];
  for (int i=1; i<bins+1; i++){
    sprintf(name,"%s%d%s",dir.Data(),i,".root");
    file[i-1] = new TFile(name);
  }
  file[1].cd();
  recoilAnalysis r;
  for (int i=1; i<bins+1; i++){
    cout << "FILE " << i << endl;
    if (islumi) {sprintf(namein,"%s%d%s",dirin.Data(),i,"-lumi");}
    else {sprintf(namein,"%s%d%s",dirin.Data(),i,"-totevents"); }
    ifstream infile(namein);
    infile >> lumi[i-1];
    errlumi[i-1] = lumi[i-1]/75.;
    file[i-1].cd();
    xvalues[i-1] = i;
    errxvalues[i-1] = 0;
    mesData tempthemes(*(r.newMes(mesrecoil, 1, 1)));
    meslepyvalues[i-1] = tempthemes.theSig()/lumi[i-1];
    mesleperryvalues[i-1] = tempthemes.theErrSig()/lumi[i-1];
    mesData tempthemes(*(r.newMes(mesallcuts, 1, 1)));
    mesallcutsyvalues[i-1] = tempthemes.theSig()/lumi[i-1];
    mesallcutserryvalues[i-1] = tempthemes.theErrSig()/lumi[i-1];
    file[i-1].cd("breco");    
    mesData tempthemes(*(r.newMes(h13001, 1, 0)));
    mesDpiyvalues[i-1] = tempthemes.theSig()/lumi[i-1];
    mesDpierryvalues[i-1] = tempthemes.theErrSig()/lumi[i-1];
    
  }

  TGraphErrors *meslpeyields;
  meslepyields = new TGraphErrors(bins,xvalues,meslepyvalues,errxvalues,mesleperryvalues);
  TGraphErrors *mesallyields;
  mesallyields = new TGraphErrors(bins,xvalues,mesallcutsyvalues,errxvalues,mesallcutserryvalues);
  TGraphErrors *mesDpiyields;
  mesDpiyields = new TGraphErrors(bins,xvalues,mesDpiyvalues,errxvalues,mesDpierryvalues);
  TGraphErrors *thelumi;
  thelumi = new TGraphErrors(bins,xvalues,lumi,errxvalues,errlumi);
  c0.Clear();
  meslepyields->Draw("AP");
  c0.SaveAs("hist_lep.ps");
  c0.Clear();
  mesallyields->Draw("AP");
  c0.SaveAs("hist_all.ps");
  c0.Clear();
  mesDpiyields->Draw("AP");
  c0.SaveAs("hist_Dpi.ps");
  c0.Clear();
  thelumi->Draw("AP");
  c0.SaveAs("hist_lumi.ps");
}
